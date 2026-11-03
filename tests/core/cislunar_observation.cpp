/**
 * @file cislunar_observation.cpp
 * @brief 检查状态插值、引力时延、跨日光行时及钟差符号。
 * @details 使用解析多项式和固定 BCRS 端点构造可独立核对的结果。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include "utils_constants.hpp"
#include <sofa.h>
#include <sofam.h>
#include <cmath>
#include <iostream>
#include <sstream>

int main(int argc, char* argv[])
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);

    if (argc != 2 || Read_JPLEPH(argv[1]) < 0) return 1;

    // =================================== 01 多项式插值与亚微秒历元 ===================================
    std::vector<Cislunar::StateRecord> states;

    for (int i = 0; i < 10; ++i)
    {
        Cislunar::StateRecord state;
        state.year = 2029; state.month = 5; state.day = 20;
        state.second = i * 2.0;
        state.position_ebcrs_m = {1000.0 + state.second * state.second, 0.0, 0.0};
        state.velocity_ebcrs_mps = {2.0 * state.second, 0.0, 0.0};
        state.clock_offset_s = 1e-6 + state.second * 1e-9;
        states.push_back(state);
    }

    Cislunar::StateRecord query = states.front();
    query.second = 7.123456789123;
    Cislunar::InterpolateState(tracer, states, query);
    if (tracer.Code() != bastion::ERR_NOERROR) tracer.EmergencyExit();

    if (std::abs(query.position_ebcrs_m[0] - (1000.0 + query.second * query.second)) > 1e-10
        || std::abs(query.velocity_ebcrs_mps[0] - 2.0 * query.second) > 1e-12
        || std::abs(query.clock_offset_s - (1e-6 + query.second * 1e-9)) > 1e-20) return 2;

    Cislunar::StateRecord shifted = query;
    shifted.second += 1e-10;
    Cislunar::InterpolateState(tracer, states, shifted);
    if (tracer.Code() != bastion::ERR_NOERROR) tracer.EmergencyExit();

    if (std::abs((shifted.position_ebcrs_m[0] - query.position_ebcrs_m[0]) - 2.0 * query.second * 1e-10) > 2e-12) return 3;

    // =================================== 02 单天体解析引力时延 ===================================
    std::array<std::array<double, 3>, 11> bodies{};
    double delay_s = 0.0;
    Cislunar::ComputeGravityDelay(tracer, {2e8, 0.0, 0.0}, {3e8, 0.0, 0.0}, bodies, {3}, delay_s);
    if (tracer.Code() != bastion::ERR_NOERROR) tracer.EmergencyExit();

    const double c = SPEED_OF_LIGHT_MPS;
    const double expected_delay = 2.0 * MAJOR_BODY_GM[2] / (c * c * c) * std::log(1.5);

    if (std::abs(delay_s - expected_delay) > 1e-24) return 4;

    // =================================== 03 跨日发射与接收钟读数 ===================================
    Cislunar::StateRecord receive;
    receive.year = 2029; receive.month = 5; receive.day = 20;
    receive.second = 0.2;
    receive.clock_offset_s = -0.3;
    double origin, mjd, earth[6];
    iauCal2jd(receive.year, receive.month, receive.day, &origin, &mjd);
    pleph(origin + mjd, receive.second / DAYSEC, 3, 12, earth);

    // 两个端点在 BCRS 中静止且相距一个光秒；E-BCRS 坐标随地球平移。
    const std::array<double, 3> transmit_bcrs{earth[0] * 1000.0, earth[1] * 1000.0, earth[2] * 1000.0};
    receive.position_ebcrs_m = {c, 0.0, 0.0};
    const auto query_transmit = [&](bastion::Overwatch&, Cislunar::StateRecord& state)
    {
        double jd0, day, body[6];
        iauCal2jd(state.year, state.month, state.day, &jd0, &day);
        pleph(jd0 + day, (state.hour * 3600.0 + state.minute * 60.0 + state.second) / DAYSEC, 3, 12, body);

        for (int i = 0; i < 3; ++i) state.position_ebcrs_m[i] = transmit_bcrs[i] - body[i] * 1000.0;

        state.clock_offset_s = 2e-6;
    };
    Cislunar::OneWayRangeRecord observation;
    Cislunar::ComputeOneWayRange(tracer, receive, query_transmit, {}, observation);
    if (tracer.Code() != bastion::ERR_NOERROR) tracer.EmergencyExit();

    if (std::abs(observation.geometric_m - c) > 5e-5 || std::abs(observation.light_time_s - 1.0) > 2e-13
        || observation.transmit_state.day != 19 || observation.transmit_state.hour != 23
        || observation.transmit_state.minute != 59 || std::abs(observation.transmit_state.second - 59.2) > 2e-13
        || observation.receive_clock_calendar[2] != 19 || std::abs(observation.receive_clock_second - 59.9) > 1e-13
        || std::abs(observation.pseudorange_m - c * (1.0 - 0.3 - 2e-6)) > 5e-5) return 5;

    // =================================== 04 错误调用保留追踪位置 ===================================
    bastion::Overwatch error_tracer(__FILE__, __func__, __LINE__);
    query.second = 30.0;
    Cislunar::InterpolateState(error_tracer, states, query);
    std::ostringstream trace;
    error_tracer.ReportTrace(trace);

    if (error_tracer.Code() == bastion::ERR_NOERROR || trace.str().find("InterpolateState") == std::string::npos) return 6;

    tracer.CheckOut(__func__);
    std::cout << "Passed: interpolation, epoch precision, gravity, midnight events, clock sign and error trace.\n";
    return 0;
}
