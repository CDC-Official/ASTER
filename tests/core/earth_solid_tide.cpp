/**
 * @file earth_solid_tide.cpp
 * @brief 对照 NAO 的固定点加速度，并用独立势函数梯度检查固体潮及动力学接入。
 */
#include "accel_complete.hpp"
#include "model_propagation.hpp"
#include "orbit_geometry.hpp"
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

/**
 * @brief 计算扣除永久部分后的二阶潮汐势，供数值梯度独立检查加速度。
 * @param[in] position_m 地固系卫星位置 [x, y, z]，单位 m。
 * @param[in] bodies_m 同一地固系的太阳和月球位置，顺序为 [太阳 x,y,z, 月球 x,y,z]，单位 m。
 * @return 在原始 GGM03C 上额外叠加的引力势，单位 m²/s²。
 * @details 按二阶勒让德势函数计算日月贡献，并减去 CSR 给出的永久 C20 势。
 * 使用固定的 SI 参数，数值求导过程与生产代码的向量加速度表达式分别计算。
 */
static long double TidePotential(const std::array<long double, 3>& position_m, const double bodies_m[6])
{
    // 固定点参考采用 GGM03C 半径、k2=0.3 和 DE440 日月 GM 的 TT 相容值。
    const long double reference_radius_m = 6378136.3L;
    const long double tt_scale = (1.0L - 6.969290134e-10L) / (1.0L - 1.550519768e-8L);
    const long double gm[2] = {1.32712440041279419e20L * tt_scale, 4.902800118e12L * tt_scale};
    const long double radius = std::sqrt(position_m[0] * position_m[0]
        + position_m[1] * position_m[1] + position_m[2] * position_m[2]);
    long double potential = 0.0L;

    // 势函数仅使用夹角的二阶勒让德多项式，随后由有限差分求三个位置偏导。
    for (int body = 0; body < 2; ++body)
    {
        long double distance_squared = 0.0L, dot_product = 0.0L;

        for (int axis = 0; axis < 3; ++axis)
        {
            const long double value = bodies_m[3 * body + axis];
            distance_squared += value * value;
            dot_product += position_m[axis] * value;
        }

        const long double distance = std::sqrt(distance_squared);
        const long double cos_angle = dot_product / (radius * distance);
        potential += 0.3L * gm[body] * std::pow(reference_radius_m, 5)
            * (3.0L * cos_angle * cos_angle - 1.0L) / (2.0L * std::pow(distance, 3) * std::pow(radius, 3));
    }

    // 原始 GGM03C 含有 -4.173e-9 的永久 C20；附加势中加入其相反数。
    const long double sin_latitude = position_m[2] / radius;
    potential += 3.986004415e14L * std::pow(reference_radius_m, 2) * 4.173e-9L
        * std::sqrt(5.0L) * (3.0L * sin_latitude * sin_latitude - 1.0L) / (2.0L * std::pow(radius, 3));
    return potential;
}

int main(int argc, char** argv)
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    using namespace constants::Scale;
    const double acceleration_unit = UNIT_GM_M3_S2 / (UNIT_LENGTH_M * UNIT_LENGTH_M);

    // 固定点参考与实际动力学检查共用本轮的 EOP、JPL 路径。
    if (argc != 4)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected EOP, JPL and NAO reference paths.", bastion::EventType::ERROR);
        return 1;
    }

    std::ifstream reference(argv[3]);

    if (!reference)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open NAO solid-tide reference.", bastion::EventType::ERROR);
        return 1;
    }

    // 比较普通低轨、高纬和高轨位置；每行保存三个位置及 NAO 完整潮汐加速度。
    int records = 0;
    double maximum_nao_difference = 0.0, maximum_gradient_difference = 0.0;
    double record[12];

    while (reference >> record[0])
    {
        for (int column = 1; column < 12; ++column) reference >> record[column];

        if (!reference)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Incomplete NAO solid-tide reference row.", bastion::EventType::ERROR);
            return 1;
        }

        // 生产接口使用积分长度单位，参考文件中的位置和加速度均为 SI 单位。
        double satellite[3], sun[3], moon[3], acceleration[3];

        for (int axis = 0; axis < 3; ++axis)
        {
            satellite[axis] = record[axis] / UNIT_LENGTH_M;
            sun[axis] = record[axis + 3] / UNIT_LENGTH_M;
            moon[axis] = record[axis + 6] / UNIT_LENGTH_M;
        }

        Earth_Solid_Tide_Acceleration(tracer, acceleration, satellite, sun, moon);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        // 先加回解析永久项，再与 NAO 完整二阶潮比较，将模型约定与算法差异分开。
        const double radius = std::sqrt(record[0] * record[0] + record[1] * record[1] + record[2] * record[2]);
        const double z_fraction = record[2] * record[2] / (radius * radius);
        const double factor = 1.5 * 3.986004415e14 * std::sqrt(5.0) * 4.173e-9
            * std::pow(6378136.3, 2) / std::pow(radius, 5);

        for (int axis = 0; axis < 3; ++axis)
        {
            const double actual = acceleration[axis] * acceleration_unit;
            const double permanent = factor * record[axis] * (5.0 * z_fraction - (axis == 2 ? 3.0 : 1.0));
            maximum_nao_difference = std::max(maximum_nao_difference, std::fabs(actual + permanent - record[axis + 9]));

            // 五点中心差分对独立势函数求导，100 m 间隔用于平衡截断和舍入。
            const std::array<long double, 3> position = {record[0], record[1], record[2]};
            long double potentials[4];
            const int offsets[4] = {-2, -1, 1, 2};

            for (int sample = 0; sample < 4; ++sample)
            {
                auto shifted = position;
                shifted[axis] += offsets[sample] * 100.0L;
                potentials[sample] = TidePotential(shifted, record + 3);
            }

            const double expected = static_cast<double>((potentials[0] - 8.0L * potentials[1]
                + 8.0L * potentials[2] - potentials[3]) / 1200.0L);
            maximum_gradient_difference = std::max(maximum_gradient_difference, std::fabs(actual - expected));
        }

        ++records;
    }

    if (records != 5 || maximum_nao_difference > 5e-20 || maximum_gradient_difference > 5e-17)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Fixed-point NAO or potential-gradient solid-tide mismatch.", bastion::EventType::ERROR);
        return 1;
    }

    // 在地轴两极直接计算，检查向量公式的输出及其与独立势函数梯度的关系。
    const double sun[3] = {1.4e11 / UNIT_LENGTH_M, 4e10 / UNIT_LENGTH_M, 2e10 / UNIT_LENGTH_M};
    const double moon[3] = {3e8 / UNIT_LENGTH_M, -2e8 / UNIT_LENGTH_M, 1e8 / UNIT_LENGTH_M};
    const double bodies_m[6] = {1.4e11, 4e10, 2e10, 3e8, -2e8, 1e8};

    for (double z_m : {-7000000.0, 7000000.0})
    {
        const double pole[3] = {0.0, 0.0, z_m / UNIT_LENGTH_M};
        double acceleration[3];
        Earth_Solid_Tide_Acceleration(tracer, acceleration, pole, sun, moon);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        for (int axis = 0; axis < 3; ++axis)
        {
            std::array<long double, 3> plus = {0.0L, 0.0L, z_m}, minus = plus;
            plus[axis] += 10.0L;
            minus[axis] -= 10.0L;
            const double expected = static_cast<double>((TidePotential(plus, bodies_m) - TidePotential(minus, bodies_m)) / 20.0L);

            if (!std::isfinite(acceleration[axis]) || std::fabs(acceleration[axis] * acceleration_unit - expected) > 5e-17)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Solid-tide acceleration at the pole differs from potential gradient.", bastion::EventType::ERROR);
                return 1;
            }
        }
    }

    // 准备同一物理状态，分别检查两套右函数实际加入的固体潮贡献。
    const int jpl_status = Read_JPLEPH(argv[2]);

    if (jpl_status <= 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot load solid-tide test JPL ephemeris.", bastion::EventType::ERROR);
        return 1;
    }

    OrbitPropagationJob job;
    std::string eop_path = argv[1];
    job.opn_EOP.Read_EOP(tracer, eop_path);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    job.opn_model.use_earth_NSG = true;
    job.opn_model.earth_zonal_order = 2;
    job.opn_model.earth_tesseral_degree = 2;
    const double state[6] = {7000000.0 / UNIT_LENGTH_M, 2000000.0 / UNIT_LENGTH_M, 1000000.0 / UNIT_LENGTH_M,
        -2000.0 * UNIT_TIME_S / UNIT_LENGTH_M, 6500.0 * UNIT_TIME_S / UNIT_LENGTH_M, 1500.0 * UNIT_TIME_S / UNIT_LENGTH_M};
    orbit::Cartesian cart;
    cart.Set(state[0], state[1], state[2], state[3], state[4], state[5], true);
    orbit::NonSingular2 nsg2;
    cart.ToNonSingular2(tracer, nsg2);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    double nsg_state[6];
    nsg2.ToArray(nsg_state);
    aEOP eop;
    const double tt = 55256.0, utc = TT2UTC(tt);
    job.opn_EOP.Get_EOP(tracer, utc, eop, 1);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    Matrix3d gcrs_to_itrs;
    GCRS_to_ITRS_Matrix(utc, eop, gcrs_to_itrs);
    double sun_gcrs[6], moon_gcrs[6];
    const int sun_status = Get_JPLEPH(tt, 11, 3, sun_gcrs);
    if (sun_status <= 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot query test Sun position.", bastion::EventType::ERROR);
        return 1;
    }

    const int moon_status = Get_JPLEPH(tt, 10, 3, moon_gcrs);
    if (moon_status <= 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot query test Moon position.", bastion::EventType::ERROR);
        return 1;
    }

    // 直接调用核心函数，得到在所有历表选项下都应加入的同一份 JPL 潮汐修正。
    const Vector3d satellite_itrs = gcrs_to_itrs * Vector3d(state[0], state[1], state[2]);
    const Vector3d sun_itrs = gcrs_to_itrs * Vector3d(sun_gcrs[0], sun_gcrs[1], sun_gcrs[2]);
    const Vector3d moon_itrs = gcrs_to_itrs * Vector3d(moon_gcrs[0], moon_gcrs[1], moon_gcrs[2]);
    double tide_itrs[3];
    Earth_Solid_Tide_Acceleration(tracer, tide_itrs, satellite_itrs.data(), sun_itrs.data(), moon_itrs.data());
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    const Vector3d expected_gcrs = gcrs_to_itrs.transpose() * Vector3d(tide_itrs[0], tide_itrs[1], tide_itrs[2]);
    Matrix3d inertial_to_rtn;
    RTN_Matrix_from_Cart(state[0], state[1], state[2], state[3], state[4], state[5], inertial_to_rtn);
    const Vector3d expected_rtn = inertial_to_rtn * expected_gcrs;
    double zero[3] = {0.0}, expected_rates[6], kepler_rates[6];
    GaussPlanetaryEquation_RTNtoNSG2(tracer, nsg2, expected_rtn.data(), expected_rates);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    GaussPlanetaryEquation_RTNtoNSG2(tracer, nsg2, zero, kepler_rates);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    double maximum_rhs_difference = 0.0, maximum_nsg_difference = 0.0;
    const int modes[][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 0}, {0, 1}};

    for (const auto& mode : modes)
    {
        job.opn_model.solar_gravity_config = mode[0];
        job.opn_model.lunar_gravity_config = mode[1];
        double off[6], on[6], nsg_off[6], nsg_on[6];
        job.opn_model.UseSolidTide = false;
        Acceleration_Cart(tracer, job, tt, 6, state, off);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        Acceleration_NSG2(tracer, job, tt, 6, nsg_state, nsg_off);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        // 只切换固体潮，日月直接引力及其历表选项均保留。
        job.opn_model.UseSolidTide = true;
        Acceleration_Cart(tracer, job, tt, 6, state, on);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        Acceleration_NSG2(tracer, job, tt, 6, nsg_state, nsg_on);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        for (int axis = 0; axis < 3; ++axis)
        {
            maximum_rhs_difference = std::max(maximum_rhs_difference,
                std::fabs(on[axis + 3] - off[axis + 3] - expected_gcrs[axis]) * acceleration_unit);
        }

        for (int component = 0; component < 6; ++component)
        {
            maximum_nsg_difference = std::max(maximum_nsg_difference,
                std::fabs(nsg_on[component] - nsg_off[component] - (expected_rates[component] - kepler_rates[component])));
        }
    }

    if (maximum_rhs_difference > 2e-14 || maximum_nsg_difference > 2e-14)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Solid-tide contribution differs between core calculation and dynamics.", bastion::EventType::ERROR);
        return 1;
    }

    // 无效位置应保留固体潮函数的追踪现场，输出数组保持调用前的值。
    bastion::Overwatch invalid_tracer(__FILE__, __func__, __LINE__);
    const double invalid_position[3] = {0.0, 0.0, 0.0};
    double unchanged[3] = {1.0, 2.0, 3.0};
    Earth_Solid_Tide_Acceleration(invalid_tracer, unchanged, invalid_position, sun, moon);
    std::ostringstream trace;
    invalid_tracer.ReportTrace(trace);

    if (invalid_tracer.Code() == bastion::ERR_NOERROR || trace.str().find("Earth_Solid_Tide_Acceleration") == std::string::npos
        || unchanged[0] != 1.0 || unchanged[1] != 2.0 || unchanged[2] != 3.0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid solid-tide position lost its error trace or changed output.", bastion::EventType::ERROR);
        return 1;
    }

    std::cout << std::scientific << std::setprecision(17)
        << "nao_acceleration_max_mps2=" << maximum_nao_difference << '\n'
        << "potential_gradient_max_mps2=" << maximum_gradient_difference << '\n'
        << "cartesian_rhs_max_mps2=" << maximum_rhs_difference << '\n'
        << "nsg2_rhs_max_normalized=" << maximum_nsg_difference << '\n';
    tracer.CheckOut(__func__);
    return 0;
}
