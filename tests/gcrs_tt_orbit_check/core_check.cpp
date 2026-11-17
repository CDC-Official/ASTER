/**
 * @file core_check.cpp
 * @brief 核对零 EOP 的启用边界、DAT 保留，以及额外行星的 TT 第三体加速度。
 */
#include "ASTER.h"
#include <cmath>
#include <iostream>
#include <limits>

/**
 * @brief 使用两日 EOP 样本和 DE440，检查缺测策略及牛顿差分引力。
 * @param[in] argc 参数个数，应为 3。
 * @param[in] argv 依次为程序名、两日 EOP 样本路径和 DE440 路径。
 * @return 0 表示全部检查通过，非零表示文件、右函数或数值检查失败。
 */
int main(int argc, char** argv)
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);

    if (argc != 3) return 2;

    // 严格策略下，未读资料的查询须保留输出并报错。
    EOP strict;
    aEOP value;
    value.xp = 99;
    bastion::Overwatch strict_tracer(__FILE__, __func__, __LINE__);
    strict.Get_EOP(strict_tracer, 61892.0, value);
    if (strict_tracer.Code() == ERR_NOERROR || value.xp != 99) return 1;

    // 显式零 EOP 策略保留 SOFA 的 DAT；非法时刻仍报错。
    OrbitPropagationJob job;
    job.opn_EOP.SetZeroFallback(true);
    job.opn_EOP.Get_EOP(tracer, 61892.25, value);
    if (tracer.Code() != ERR_NOERROR) return 1;

    if (value.xp != 0 || value.yp != 0 || value.dx00 != 0 || value.dy00 != 0
        || value.dut1 != 0 || value.dat != 37 || job.opn_EOP.ZeroFallbackCount() != 1) return 1;

    bastion::Overwatch nan_tracer(__FILE__, __func__, __LINE__);
    job.opn_EOP.Get_EOP(nan_tracer, std::numeric_limits<double>::quiet_NaN(), value);
    if (nan_tracer.Code() == ERR_NOERROR) return 1;

    // 有资料时仍插值；最后一天零时可用，其后缺少插值端点才回退。
    std::string eop_file = argv[1];
    job.opn_EOP.Read_EOP(tracer, eop_file);
    if (tracer.Code() != ERR_NOERROR) return 1;

    job.opn_EOP.Get_EOP(tracer, 61892.5, value);
    if (tracer.Code() != ERR_NOERROR) return 1;

    if (std::abs(value.xp - 2 * DAS2R) > 1e-20 || std::abs(value.dut1 - 0.2) > 1e-15
        || job.opn_EOP.ZeroFallbackCount() != 1) return 1;

    job.opn_EOP.Get_EOP(tracer, 61893.0, value);
    if (tracer.Code() != ERR_NOERROR) return 1;

    if (std::abs(value.xp - 3 * DAS2R) > 1e-20) return 1;

    job.opn_EOP.Get_EOP(tracer, 61893.25, value);
    if (tracer.Code() != ERR_NOERROR) return 1;

    if (value.xp != 0 || value.dut1 != 0 || value.dat != 37 || job.opn_EOP.ZeroFallbackCount() != 2) return 1;

    // 用同一 Cartesian 位置比较无行星与逐颗行星的加速度，独立写出长双精度差分公式。
    const int loaded = Read_JPLEPH(argv[2]);
    if (loaded <= 0) return 1;

    const double tt = UTC2TT(61892.0);
    const double state[6] = {10, 20, 30, 0, 0, 0};
    double base[6];
    Acceleration_Cart(tracer, job, tt, 6, state, base);
    if (tracer.Code() != ERR_NOERROR) return 1;

    for (const int body : {1, 2, 4, 5, 6, 7, 8, 9})
    {
        job.opn_model.cartesian_planetary_body_id = {body};
        double actual[6], planet[6];
        Acceleration_Cart(tracer, job, tt, 6, state, actual);
        if (tracer.Code() != ERR_NOERROR) return 1;

        const int queried = Get_JPLEPH(tt, body, 3, planet);
        if (queried <= 0) return 1;

        // 从行星到目标与地心的距离分别计算直接项和间接项。
        long double distance2 = 0, planet2 = 0;

        for (int axis = 0; axis < 3; ++axis)
        {
            const long double delta = static_cast<long double>(planet[axis]) - state[axis];
            distance2 += delta * delta;
            planet2 += static_cast<long double>(planet[axis]) * planet[axis];
        }

        const long double gm = constants::DE440::GM_TDB[body - 1] * constants::TDB_TO_TT_SCALE / constants::Scale::UNIT_GM_M3_S2;

        for (int axis = 0; axis < 3; ++axis)
        {
            const long double expected = gm * ((planet[axis] - static_cast<long double>(state[axis])) / std::pow(distance2, 1.5L)
                - planet[axis] / std::pow(planet2, 1.5L));

            if (std::abs((actual[axis + 3] - base[axis + 3]) - expected) > 2e-19L) return 1;
        }
    }

    tracer.CheckOut(__func__);
    std::cout << "EOP coverage/DAT and eight planetary accelerations passed.\n";
    return 0;
}
