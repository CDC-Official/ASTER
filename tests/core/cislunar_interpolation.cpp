/**
 * @file cislunar_interpolation.cpp
 * @brief 用解析多项式检查状态插值的窗口和端点，并核对 JPL 拆分日期的亚微秒响应。
 * @details 状态使用绝对 TDB 日历历元，位置、速度和钟差分别与解析值比较。
 * 采样序列跨越午夜；唯一命令行参数为测试所用的 JPL 历表路径。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include <cmath>
#include <iostream>

int main(int argc, char** argv)
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);

    // 从 2029 年 5 月 19 日 23:40 起，每分钟采样一次，覆盖午夜前后的插值窗口。
    std::vector<Cislunar::StateRecord> states;

    for (int i = 0; i < 50; ++i)
    {
        Cislunar::StateRecord state;
        const int minute = 23 * 60 + 40 + i;
        state.year = 2029; state.month = 5; state.day = 19 + minute / 1440;
        state.hour = minute % 1440 / 60;
        state.minute = minute % 60;

        // 位置为三次多项式，速度为其时间导数，钟差为线性函数。
        const double elapsed_s = i * 60.0, t = elapsed_s / 1000.0;
        state.position_ebcrs_m = {3e8 + 1e6 * t - 2e5 * t * t + 1e4 * t * t * t,
            -2e8 + 2e5 * t + 3e3 * t * t, 5e7 - 7e4 * t};
        state.velocity_ebcrs_mps = {1000.0 - 400.0 * t + 30.0 * t * t, 200.0 + 6.0 * t, -70.0};
        state.clock_offset_s = 1e-4 + elapsed_s * 1e-10;
        states.push_back(state);
    }

    // 同时检查首末端点、采样点附近以及中部和跨日窗口中的非节点时刻。
    for (double reference_s : {0.0, 5.0, 540.0, 545.0, 1045.0, 1205.0, 2930.0, 2940.0})
    {
        const double elapsed_s = reference_s == 0.0 || reference_s == 2940.0 ? reference_s : reference_s - 1.3;
        const double seconds = 23 * 3600.0 + 40 * 60.0 + elapsed_s;
        Cislunar::StateRecord query;
        query.year = 2029; query.month = 5; query.day = 19 + static_cast<int>(seconds / 86400.0);
        query.hour = static_cast<int>(std::fmod(seconds, 86400.0) / 3600.0);
        query.minute = static_cast<int>(std::fmod(seconds, 3600.0) / 60.0);
        query.second = std::fmod(seconds, 60.0);

        Cislunar::InterpolateState(tracer, states, query);
        if (tracer.Code() != bastion::ERR_NOERROR) return 1;

        // 直接计算查询时刻的解析状态，分别核对三个位置和速度分量。
        const double t = elapsed_s / 1000.0;
        const std::array<double, 3> expected_position = {3e8 + 1e6 * t - 2e5 * t * t + 1e4 * t * t * t,
            -2e8 + 2e5 * t + 3e3 * t * t, 5e7 - 7e4 * t};
        const std::array<double, 3> expected_velocity = {1000.0 - 400.0 * t + 30.0 * t * t, 200.0 + 6.0 * t, -70.0};

        for (int i = 0; i < 3; ++i)
        {
            if (std::abs(query.position_ebcrs_m[i] - expected_position[i]) > 2e-5 ||
                std::abs(query.velocity_ebcrs_mps[i] - expected_velocity[i]) > 1e-9)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "State interpolation mismatch at elapsed seconds " +
                    std::to_string(elapsed_s) + ", component " + std::to_string(i), bastion::EventType::ERROR);
                return 1;
            }
        }

        if (std::abs(query.clock_offset_s - (1e-4 + elapsed_s * 1e-10)) > 1e-17)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Clock interpolation mismatch at elapsed seconds " +
                std::to_string(elapsed_s), bastion::EventType::ERROR);
            return 1;
        }
    }

    // 1 微秒内地球约移动 3 厘米；拆分日期应保留这一变化。
    if (argc != 2 || Read_JPLEPH(argv[1]) < 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected one readable JPL ephemeris path.", bastion::EventType::ERROR);
        return 1;
    }

    double first[6], next[6];
    pleph(2462276.5, 0.25, 3, 12, first);
    pleph(2462276.5, 0.25 + 1e-6 / 86400.0, 3, 12, next);

    for (int i = 0; i < 3; ++i)
    {
        const double difference_m = (next[i] - first[i]) * 1000.0;
        const double expected_m = first[i + 3] * 1e-3;

        if (!std::isfinite(difference_m) || !std::isfinite(expected_m) || std::abs(difference_m - expected_m) > 1e-4)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "JPL split date lost the one-microsecond displacement.", bastion::EventType::ERROR);
            return 1;
        }
    }

    tracer.CheckOut(__func__);
    std::cout << "Passed: state interpolation across windows and midnight, endpoints, and split-date ephemeris.\n";
    return 0;
}
