/**
 * @file main_prediction.cpp
 * @brief 直接调用ASTER的F、G级数，输出给定初始状态的完整预报序列。
 * @details 输入采用GCRS位置m、速度m/s和相对历元时间s；库内计算采用ASTER归一化尺度。
 * 本程序只做轨道预报，不构造观测，也不进行初轨迭代。位置直接取自现有系数，
 * 速度由同一位置表达式的五点中心差分得到，差分步长另存以便核查。
 */
#include "core_iod.hpp"
#include "data_gravityfield.hpp"
#include <array>

using bastion::Overwatch;
using constants::Scale::UNIT_LENGTH_M;
using constants::Scale::UNIT_TIME_S;

/**
 * @brief 用ASTER现有F、G表达式计算相对参考历元的位置。
 * @param[in] model 二体模型1或二体加J2模型2。
 * @param[in] epoch_mjd_utc 参考历元，单位UTC MJD；这两种自治模型不依赖该日期。
 * @param[in] tau_norm 相对参考历元的有符号时间，除以UNIT_TIME_S。
 * @param[in] state_norm 六维GCRS参考状态，位置除以UNIT_LENGTH_M，速度乘UNIT_TIME_S/UNIT_LENGTH_M。
 * @return 三维GCRS位置，除以UNIT_LENGTH_M。
 */
Eigen::Vector3d PredictPosition(int model, double epoch_mjd_utc, double tau_norm,
                                const std::array<double, 6>& state_norm)
{
    // 库接口接收非常量指针；复制参考状态，保证所有时刻从同一初始状态展开。
    std::array<double, 6> state = state_norm;
    double coefficients[4];
    Lagrange_Coefficients(model, epoch_mjd_utc, tau_norm, state.data(), coefficients);

    // 赤道平面的两个分量使用f/g，极轴分量使用fz/gz。
    return Eigen::Vector3d(
        coefficients[0] * state[0] + coefficients[1] * state[3],
        coefficients[0] * state[1] + coefficients[1] * state[4],
        coefficients[2] * state[2] + coefficients[3] * state[5]);
}

int main(int argc, char** argv)
{
    Overwatch tracer(__FILE__, __func__, __LINE__);
    std::cout << std::setprecision(17);

    // ============================== 01 提供与编译库一致的常数 ==============================
    if (argc == 2 && std::string(argv[1]) == "--constants")
    {
        std::cout << "{\"gm_m3_s2\":" << constants::Scale::UNIT_GM_M3_S2
                  << ",\"unit_length_m\":" << UNIT_LENGTH_M
                  << ",\"unit_time_s\":" << UNIT_TIME_S
                  << ",\"j2_normalized\":" << constants::GGM03C::NORMALIZED_J2
                  << ",\"reference_radius_m\":" << constants::GGM03C::REFERENCE_RADIUS << "}\n";
        tracer.CheckOut(__func__);
        return 0;
    }

    if (argc != 3)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Usage: aster_test_iod_lagrange_prediction input.txt output.csv, or --constants", bastion::ERROR);
        tracer.EmergencyExit();
    }

    // ============================== 02 读取物理初态与输出时刻 ==============================
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    int count = 0;
    input >> count;

    if (!input || !output || count <= 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open prediction input/output or invalid case count", bastion::ERROR);
        tracer.EmergencyExit();
    }

    output << std::setprecision(17)
           << "case_id,time_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,velocity_step_s\n";

    for (int index = 0; index < count; ++index)
    {
        std::string identifier;
        int model, samples;
        double epoch_mjd_utc, velocity_step_s;
        std::array<double, 6> state;
        input >> identifier >> model >> epoch_mjd_utc >> samples >> velocity_step_s;

        for (double& component : state) input >> component;

        if (!input || (model != IOD_DYNAMICS_TWOBODY && model != IOD_DYNAMICS_J2)
            || samples <= 0 || !(velocity_step_s > 0.0))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid case header: " + identifier, bastion::ERROR);
            tracer.EmergencyExit();
        }

        // 归一化只在接口边界执行；所有输出再还原为m和m/s。
        for (int j = 0; j < 3; ++j) state[j] /= UNIT_LENGTH_M;

        for (int j = 3; j < 6; ++j) state[j] *= UNIT_TIME_S / UNIT_LENGTH_M;

        const double derivative_step = velocity_step_s / UNIT_TIME_S;

        // ============================== 03 对同一初态逐时展开 ==============================
        for (int sample = 0; sample < samples; ++sample)
        {
            double time_s;
            input >> time_s;

            if (!input || !std::isfinite(time_s))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid prediction time: " + identifier, bastion::ERROR);
                tracer.EmergencyExit();
            }

            // 主比较使用库直接输出的位置，包含正向和反向的完整时间序列。
            const double tau = time_s / UNIT_TIME_S;
            const Eigen::Vector3d position = PredictPosition(model, epoch_mjd_utc, tau, state);

            // 同一表达式的时间导数作为辅助速度结果，后续用半步长重复核查。
            const Eigen::Vector3d velocity = (
                -PredictPosition(model, epoch_mjd_utc, tau + 2.0 * derivative_step, state)
                + 8.0 * PredictPosition(model, epoch_mjd_utc, tau + derivative_step, state)
                - 8.0 * PredictPosition(model, epoch_mjd_utc, tau - derivative_step, state)
                + PredictPosition(model, epoch_mjd_utc, tau - 2.0 * derivative_step, state)
            ) / (12.0 * derivative_step);

            if (!position.allFinite() || !velocity.allFinite())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite F/G prediction: " + identifier, bastion::ERROR);
                tracer.EmergencyExit();
            }

            output << identifier << ',' << time_s;

            for (int j = 0; j < 3; ++j) output << ',' << position[j] * UNIT_LENGTH_M;

            for (int j = 0; j < 3; ++j) output << ',' << velocity[j] * UNIT_LENGTH_M / UNIT_TIME_S;

            output << ',' << velocity_step_s << '\n';
        }
    }

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed to save F/G prediction sequence", bastion::ERROR);
        tracer.EmergencyExit();
    }

    tracer.CheckOut(__func__);
    return 0;
}
