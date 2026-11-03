/**
 * @file model_cislunar_observation.cpp
 * @brief 用固定采样窗口插值轨道与钟差，迭代求解单程光时和伪距。
 * @details 内部几何长度以 6378136.3 m 归一化，时刻采用相对天数；公开结果使用 SI 单位。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include "utils_constants.hpp"
#include <sofam.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace observation
{

/**
 * @brief 为当前接收钟时选取采样窗口，并计算切比雪夫系数。
 * @param[in,out] series 读取 days、values、is_clock，向后更新 window_number 和 window_start；观测须按钟时递增调用。
 * @param[in] clock_tag_day 接收端钟时距参考历元的天数，用于固定本条观测的插值窗口。
 * @return 当前窗口的起始时刻、跨度及 17 行系数；分量和单位与 series 相同。
 * @details 轨道窗口每次前进 9 点；钟差从第二段起按 7 点步长、提前一行取 10 点。
 * 先以 Neville 方法求 17 个 Chebyshev–Gauss–Lobatto 节点值，再计算 16 阶展开。
 */
CislunarInterpolant PrepareCislunarInterpolant(CislunarSeries& series, double clock_tag_day)
{
    const int count = static_cast<int>(series.days.size());
    const int stride = series.is_clock ? 7 : 9;
    const int last_window = count / stride;
    double end_day = series.days[series.window_start + 9];

    if (clock_tag_day < series.days.front() || clock_tag_day > series.days.back())
        throw std::runtime_error("Clock tag outside satellite table: " + std::to_string(series.satellite_id));

    // 时刻递增时推进窗口；末段取文件最后十点。
    while (end_day < clock_tag_day && series.window_number <= last_window)
    {
        ++series.window_number;
        int start = stride * (series.window_number - 1);

        if (series.is_clock) --start;

        if (series.window_number == last_window + 1 || start + 10 > count)
        {
            start = count - 10;
        }
        else
        {
            if (start + 10 == count) ++series.window_number;
        }

        series.window_start = start;
        end_day = series.days[start + (series.is_clock ? 8 : 9)];
    }

    CislunarInterpolant interpolant;
    interpolant.start_day = series.days[series.window_start];
    interpolant.span_days = series.days[series.window_start + 9] - interpolant.start_day;
    interpolant.coefficients = Eigen::MatrixXd::Zero(17, series.values.cols());
    double nodes[10], cgl[17], polynomial[17][17];

    for (int i = 0; i < 10; ++i) nodes[i] = -1.0 + 2.0 * i / 9.0;

    for (int j = 0; j <= 16; ++j)
    {
        cgl[j] = -std::cos(pi * j / 16.0);
        polynomial[j][0] = 1.0;
        polynomial[j][1] = cgl[j];

        for (int k = 2; k <= 16; ++k)
            polynomial[j][k] = 2.0 * cgl[j] * polynomial[j][k - 1] - polynomial[j][k - 2];
    }

    // 每个分量独立做 Neville 插值，再由 CGL 节点值求切比雪夫系数。
    for (int component = 0; component < series.values.cols(); ++component)
    {
        double values[17];

        for (int k = 0; k <= 16; ++k)
        {
            double work[10];

            for (int j = 0; j < 10; ++j) work[j] = series.values(series.window_start + j, component);

            for (int i = 0; i < 10; ++i)
            {
                for (int j = i - 1; j >= 0; --j)
                    work[j] = work[j + 1] + (work[j + 1] - work[j]) * (cgl[k] - nodes[i]) / (nodes[i] - nodes[j]);
            }

            values[k] = work[0];
        }

        for (int i = 0; i <= 16; ++i)
        {
            double coefficient = (values[0] * polynomial[0][i] + values[16] * polynomial[16][i]) / 16.0;

            for (int j = 1; j < 16; ++j) coefficient += values[j] * polynomial[j][i] * 2.0 / 16.0;

            interpolant.coefficients(i, component) = coefficient;
        }
    }

    return interpolant;
}

/**
 * @brief 在固定的切比雪夫窗口内计算轨道或钟差。
 * @param[in] interpolant 起始时刻、跨度和 0—16 阶系数，单位与原始采样值相同。
 * @param[in] day 待求时刻距参考 TDB 历元的天数。
 * @return 按原始列顺序排列的六维轨道或二维钟差；单位为 m、m/s 或 s。
 * @details 使用 T0=1、T1=x 的递推，零阶系数取一半；收发事件共用接收钟时选定的窗口。
 */
Eigen::VectorXd EvaluateCislunarInterpolant(const CislunarInterpolant& interpolant, double day)
{
    const double x = ((day - interpolant.start_day) - interpolant.span_days / 2.0) * 2.0 / interpolant.span_days;
    double polynomial[17] = {1.0, x};

    for (int j = 2; j <= 16; ++j) polynomial[j] = 2.0 * x * polynomial[j - 1] - polynomial[j - 2];

    Eigen::VectorXd value = interpolant.coefficients.row(0).transpose() / 2.0;

    for (int j = 1; j <= 16; ++j) value += interpolant.coefficients.row(j).transpose() * polynomial[j];

    return value;
}

/**
 * @brief 取得指定天体的 BCRS 位置，用于端点平移和引力时延。
 * @param[in] reference_mjd 参考 TDB MJD 的整数日、小数日。
 * @param[in] day 距参考历元的 TDB 天数。
 * @param[in] body JPL 天体编号，1—11，3 为地球，10 为月球，11 为太阳。
 * @return BCRS 位置，单位 m；历表由主程序预先载入。
 */
static Eigen::Vector3d CislunarBodyPosition(const std::array<double, 2>& reference_mjd, double day, int body)
{
    const double jd1 = START_OF_MJD + reference_mjd[0];
    const double jd2 = reference_mjd[1] + day;

    if (jd1 + jd2 < JPLEPH_R1.r1.ss[0] || jd1 + jd2 > JPLEPH_R1.r1.ss[1])
        throw std::runtime_error("Observation event outside JPL ephemeris.");

    double state_km[6];
    pleph(jd1, jd2, body, 12, state_km);
    return Eigen::Vector3d(state_km[0], state_km[1], state_km[2]) * 1000.0;
}

/**
 * @brief 求解一条星间单程观测的收发 TDB 时刻和理论伪距。
 * @param[in] reference_mjd 参考 TDB MJD 的整数日和小数日，用于查询已载入的 JPL 二进制历表。
 * @param[in] clock_tag_day 接收端钟时距参考历元的天数。
 * @param[in] transmit_orbit 发射星 E-BCRS 轨道插值，位置 m、速度 m/s。
 * @param[in] receive_orbit 接收星 E-BCRS 轨道插值，位置 m、速度 m/s。
 * @param[in] transmit_clock 发射星钟差插值，两列为相对论钟差和物理钟偏，单位 s。
 * @param[in] receive_clock 接收星钟差插值，定义及单位与 transmit_clock 相同。
 * @return 收发时刻、几何距离、引力时延距离、两端钟差和伪距，以及迭代次数与末次变化量。
 * @details 先反求接收 TDB，再迭代发射 TDB；单次时刻变化不超过 1e-12 s 时收敛。
 * 两端 E-BCRS 位置分别加各自事件的地球 BCRS 位置。引力天体取接收 TDB，跳过地球，
 * 太阳项的对数分子、分母均加 2GM/c²。伪距等于几何距离、引力时延距离和接收减发射钟差距离之和。
 * 设备时延、天线偏移和测量噪声取零。十轮内未收敛时报告错误。
 */
OneWayRange ComputeCislunarOneWayRange(const std::array<double, 2>& reference_mjd, double clock_tag_day,
    const CislunarInterpolant& transmit_orbit, const CislunarInterpolant& receive_orbit,
    const CislunarInterpolant& transmit_clock, const CislunarInterpolant& receive_clock)
{
    const double c = SPEED_OF_LIGHT_MPS;
    const double length_scale = 6378136.3;
    const double time_scale = 806.810991306733;
    const double speed_scale = 7905.36614984607;
    const double c_norm = c / speed_scale;
    OneWayRange result;
    result.receive_day = clock_tag_day;
    Eigen::VectorXd receive_bias;

    // =================================== 01 由接收钟时反求接收 TDB ===================================
    for (int iteration = 1; iteration <= 10; ++iteration)
    {
        receive_bias = EvaluateCislunarInterpolant(receive_clock, result.receive_day);
        const double next = clock_tag_day - (receive_bias(0) + receive_bias(1)) / DAYSEC;
        result.receive_update_s = std::abs(next - result.receive_day) * DAYSEC;
        result.receive_day = next;
        result.receive_iterations = iteration;

        if (result.receive_update_s <= 1e-12) break;
    }

    if (result.receive_update_s > 1e-12) throw std::runtime_error("Receive TDB did not converge in 10 iterations.");

    const Eigen::VectorXd receive_state = EvaluateCislunarInterpolant(receive_orbit, result.receive_day);
    const Eigen::Vector3d receive_bcrs = (receive_state.head<3>()
        + CislunarBodyPosition(reference_mjd, result.receive_day, 3)) / length_scale;

    // =================================== 02 准备接收时刻的引力天体 ===================================
    const double gm_earth = 3.986004415e14;
    std::array<Eigen::Vector3d, 11> body_positions;

    for (int body = 1; body <= 11; ++body)
    {
        if (body == 3) continue;

        body_positions[body - 1] = CislunarBodyPosition(reference_mjd, result.receive_day, body) / length_scale;
    }

    // =================================== 03 迭代发射 TDB 与传播光时 ===================================
    result.transmit_day = result.receive_day;
    double delay_s = 0.0;

    for (int iteration = 1; iteration <= 10; ++iteration)
    {
        const Eigen::VectorXd transmit_state = EvaluateCislunarInterpolant(transmit_orbit, result.transmit_day);
        const Eigen::Vector3d transmit_bcrs = (transmit_state.head<3>()
            + CislunarBodyPosition(reference_mjd, result.transmit_day, 3)) / length_scale;
        const double range_norm = (receive_bcrs - transmit_bcrs).norm();
        result.geometric_m = range_norm * length_scale;
        delay_s = 0.0;

        // 端点到天体的距离与传播距离共同给出一阶引力时延。
        for (int body = 1; body <= 11; ++body)
        {
            if (body == 3) continue;

            const double r1 = (transmit_bcrs - body_positions[body - 1]).norm() * length_scale;
            const double r2 = (receive_bcrs - body_positions[body - 1]).norm() * length_scale;
            const double gm = gm_earth * MAJOR_BODY_EARTH_GM_RATIO[body - 1];
            const double solar_term = body == 11 ? 2.0 * gm / (c * c) : 0.0;
            const double denominator = r1 + r2 - result.geometric_m + solar_term;

            if (denominator <= 0.0) throw std::runtime_error("Invalid gravitational time-delay geometry.");

            delay_s += 2.0 * gm / (c * c * c)
                * std::log((r1 + r2 + result.geometric_m + solar_term) / denominator);
        }

        const double light_time_s = range_norm / c_norm * time_scale + delay_s;
        const double next = result.receive_day - light_time_s / DAYSEC;
        result.transmit_update_s = std::abs(next - result.transmit_day) * DAYSEC;
        result.transmit_day = next;
        result.transmit_iterations = iteration;

        if (result.transmit_update_s <= 1e-12) break;
    }

    if (result.transmit_update_s > 1e-12) throw std::runtime_error("Transmit TDB did not converge in 10 iterations.");

    // =================================== 04 组合传播距离与两端钟差 ===================================
    const Eigen::VectorXd transmit_bias = EvaluateCislunarInterpolant(transmit_clock, result.transmit_day);
    result.gravity_m = c * delay_s;
    result.propagation_m = result.geometric_m + result.gravity_m;
    result.transmit_physical_s = transmit_bias(1);
    result.transmit_relativistic_s = transmit_bias(0);
    result.receive_physical_s = receive_bias(1);
    result.receive_relativistic_s = receive_bias(0);
    result.clock_m = c * ((receive_bias(0) - transmit_bias(0)) + (receive_bias(1) - transmit_bias(1)));
    result.pseudorange_m = result.propagation_m + result.clock_m;
    return result;
}

}
