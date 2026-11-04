/**
 * @file gravity_conventions.cpp
 * @brief 用有量纲二体公式、J2 公式和尺度变换核对 TT/TDB 引力参数。
 */
#include "accel_complete.hpp"
#include "model_propagation.hpp"
#include <iostream>

int main(int argc, char** argv)
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);

    // 两套右函数共用同一份 EOP 和 JPL 资料，仅打开地球中心引力。
    if (argc != 3 || Read_JPLEPH(argv[2]) < 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected EOP and readable JPL paths.", bastion::EventType::ERROR);
        return 1;
    }

    OrbitPropagationJob earth_job;
    std::string eop_path = argv[1];
    earth_job.opn_EOP.Read_EOP(tracer, eop_path);
    if (tracer.Code() != ERR_NOERROR) return 1;

    OrbitPropagationJob_EBCRS cislunar_job;
    cislunar_job.opn_EOP = earth_job.opn_EOP;

    // 直接给定 m、m/s 初值，检查还原后的加速度，避免归一化抵消错误。
    const double position_m[3] = {7000000.0, 2000000.0, 1000000.0};
    const double state[6] = {position_m[0] / ae, position_m[1] / ae, position_m[2] / ae,
        -2000.0 * UNITTIME / ae, 6500.0 * UNITTIME / ae, 1500.0 * UNITTIME / ae};
    const double radius = sqrt(position_m[0] * position_m[0] + position_m[1] * position_m[1] + position_m[2] * position_m[2]);
    double derivative_tt[6], derivative_tdb[6];

    Acceleration_Cart(tracer, earth_job, 55256.0, 6, state, derivative_tt);
    if (tracer.Code() != ERR_NOERROR) return 1;

    Acceleration_Cart_EBCRS(tracer, cislunar_job, 55256.0, 6, state, derivative_tdb);
    if (tracer.Code() != ERR_NOERROR) return 1;

    // IERS TT 值与 DE440 TDB 值分别进入对应右函数。
    for (int i = 0; i < 3; ++i)
    {
        const double expected_tt = -3.986004415e14 * position_m[i] / pow(radius, 3);
        const double expected_tdb = -3.98600435507e14 * position_m[i] / pow(radius, 3);

        if (fabs(derivative_tt[i + 3] * GE / (ae * ae) - expected_tt) > 2e-14 ||
            fabs(derivative_tdb[i + 3] * GE / (ae * ae) - expected_tdb) > 2e-14)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "TT/TDB point-mass acceleration mismatch.", bastion::EventType::ERROR);
            return 1;
        }
    }

    // J2 解析式独立使用 GGM03C 文件中的 C20、GM 和参考半径。
    const double j2 = sqrt(5.0) * 0.484169325971e-3;
    const double factor = 1.5 * 3.986004415e14 * j2 * pow(6378136.3, 2) / pow(radius, 5);
    const double z_fraction = position_m[2] * position_m[2] / (radius * radius);
    const double expected_j2[3] = {
        factor * position_m[0] * (5.0 * z_fraction - 1.0),
        factor * position_m[1] * (5.0 * z_fraction - 1.0),
        factor * position_m[2] * (5.0 * z_fraction - 3.0)};

    // 长度、时间和 GM 同乘尺度因子时，加速度应除以该因子。
    for (double scale : {1.0, 1.0 / TDB_TO_TT_SCALE})
    {
        double position_norm[3], acceleration_norm[3];

        for (int i = 0; i < 3; ++i) position_norm[i] = position_m[i] * scale / ae;

        const int status = Earth_Nonspherical_Gravitational_Acceleration(acceleration_norm, position_norm, 1, 2, 0,
            GGM03C_GM * scale, GGM03C_REFERENCE_RADIUS * scale);
        if (status < 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "GGM03C evaluation failed.", bastion::EventType::ERROR);
            return 1;
        }

        // 解析式检查参考半径的次数和输出的加速度单位。
        for (int i = 0; i < 3; ++i)
        {
            if (fabs(acceleration_norm[i] * GE / (ae * ae) - expected_j2[i] / scale) > 2e-16)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "GGM03C J2 acceleration or scaling mismatch.", bastion::EventType::ERROR);
                return 1;
            }
        }
    }

    std::cout << "TT/TDB central gravity and GGM03C J2 checks passed.\n";
    tracer.CheckOut(__func__);
    return 0;
}
