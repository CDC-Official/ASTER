/**
 * @file gravity_conventions.cpp
 * @brief 用有量纲二体、J2 公式及尺度变换核对引力参数和月球地心相对摄动。
 */
#include "accel_complete.hpp"
#include "accel_drag.hpp"
#include "model_propagation.hpp"
#include "orbit_geometry.hpp"
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
    const double state[6] = {position_m[0] / constants::Scale::UNIT_LENGTH_M, position_m[1] / constants::Scale::UNIT_LENGTH_M, position_m[2] / constants::Scale::UNIT_LENGTH_M,
        -2000.0 * constants::Scale::UNIT_TIME_S / constants::Scale::UNIT_LENGTH_M, 6500.0 * constants::Scale::UNIT_TIME_S / constants::Scale::UNIT_LENGTH_M, 1500.0 * constants::Scale::UNIT_TIME_S / constants::Scale::UNIT_LENGTH_M};
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

        if (fabs(derivative_tt[i + 3] * constants::Scale::UNIT_GM_M3_S2 / (constants::Scale::UNIT_LENGTH_M * constants::Scale::UNIT_LENGTH_M) - expected_tt) > 2e-14 ||
            fabs(derivative_tdb[i + 3] * constants::Scale::UNIT_GM_M3_S2 / (constants::Scale::UNIT_LENGTH_M * constants::Scale::UNIT_LENGTH_M) - expected_tdb) > 2e-14)
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
    for (double scale : {1.0, 1.0 / constants::TDB_TO_TT_SCALE})
    {
        double position_norm[3], acceleration_norm[3];

        for (int i = 0; i < 3; ++i) position_norm[i] = position_m[i] * scale / constants::Scale::UNIT_LENGTH_M;

        const int status = Earth_Nonspherical_Gravitational_Acceleration(acceleration_norm, position_norm, 1, 2, 0,
            constants::GGM03C::GM * scale, constants::GGM03C::REFERENCE_RADIUS * scale);
        if (status < 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "GGM03C evaluation failed.", bastion::EventType::ERROR);
            return 1;
        }

        // 解析式检查参考半径的次数和输出的加速度单位。
        for (int i = 0; i < 3; ++i)
        {
            if (fabs(acceleration_norm[i] * constants::Scale::UNIT_GM_M3_S2 / (constants::Scale::UNIT_LENGTH_M * constants::Scale::UNIT_LENGTH_M) - expected_j2[i] / scale) > 2e-16)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "GGM03C J2 acceleration or scaling mismatch.", bastion::EventType::ERROR);
                return 1;
            }
        }
    }

    // 在参考距离及两倍参考距离处，直接用 SI 光压公式核对加速度和平方反比关系。
    for (double distance_factor : {1.0, 2.0})
    {
        const double satellite[3] = {10000000.0 / constants::Scale::UNIT_LENGTH_M, 0.0, 0.0};
        const double sun[3] = {(10000000.0 + distance_factor * 149597870000.0) / constants::Scale::UNIT_LENGTH_M, 0.0, 0.0};
        double acceleration[3];
        Solar_Radition_Acceleration(acceleration, satellite, sun, 0.02 * constants::Conversion::IAM_TO_NORMALIZED_AM, 0.3);

        const double expected = -1.3 * 0.02 * 4.56e-6 / (distance_factor * distance_factor);
        const double actual = acceleration[0] * constants::Scale::UNIT_GM_M3_S2 / (constants::Scale::UNIT_LENGTH_M * constants::Scale::UNIT_LENGTH_M);

        if (fabs(actual - expected) > 1e-20 || acceleration[1] != 0.0 || acceleration[2] != 0.0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Solar radiation pressure or normalization mismatch.", bastion::EventType::ERROR);
            return 1;
        }
    }

    // 背日侧分别取圆柱阴影内、外两点，检查遮挡判定和光压加速度。
    for (bool shadowed : {true, false})
    {
        const double satellite[3] = {-10000000.0 / constants::Scale::UNIT_LENGTH_M,
            (shadowed ? 0.0 : 8000000.0) / constants::Scale::UNIT_LENGTH_M, 0.0};
        const double sun[3] = {149597870000.0 / constants::Scale::UNIT_LENGTH_M, 0.0, 0.0};
        double acceleration[3];
        Solar_Radition_Acceleration(acceleration, satellite, sun, 0.02 * constants::Conversion::IAM_TO_NORMALIZED_AM, 0.3);

        const bool zero_acceleration = acceleration[0] == 0.0 && acceleration[1] == 0.0 && acceleration[2] == 0.0;

        if (in_shadow(satellite, sun) != shadowed || zero_acceleration != shadowed ||
            !std::isfinite(acceleration[0]) || !std::isfinite(acceleration[1]) || !std::isfinite(acceleration[2]))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cylindrical Earth shadow or solar radiation acceleration mismatch.", bastion::EventType::ERROR);
            return 1;
        }
    }

    // 固定同一个物理状态、Cd、面质比和密度，以 SI 公式检查共转大气阻力。
    const double cd = 2.2, area_mass = 0.02, density = 1e-12;
    double drag[3];
    AtmosphereDrag_Acceleration(drag, state, cd, area_mass * constants::Conversion::IAM_TO_NORMALIZED_AM,
        density * pow(constants::Scale::UNIT_LENGTH_M, 3) / constants::Scale::UNIT_MASS_KG);

    const double relative_velocity[3] = {-2000.0 + 7.292115e-5 * position_m[1],
        6500.0 - 7.292115e-5 * position_m[0], 1500.0};
    const double relative_speed = sqrt(relative_velocity[0] * relative_velocity[0] +
        relative_velocity[1] * relative_velocity[1] + relative_velocity[2] * relative_velocity[2]);

    for (int i = 0; i < 3; ++i)
    {
        const double expected = -0.5 * cd * area_mass * density * relative_speed * relative_velocity[i];
        const double actual = drag[i] * constants::Scale::UNIT_GM_M3_S2 / pow(constants::Scale::UNIT_LENGTH_M, 2);

        if (fabs(actual - expected) > 1e-20)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Atmospheric rotation or drag normalization mismatch.", bastion::EventType::ERROR);
            return 1;
        }
    }

    // 用只有 J2 的月球场核对地心间接项，解析公式独立于球谐递推实现。
    const double lunar_j2 = 2.0e-4;
    moon_gravity_field.GM = 4.9028e12;
    moon_gravity_field.R_ref = 1738000.0;
    moon_gravity_field.max_degree = 2;
    moon_gravity_field.is_ready = true;
    moon_gravity_field.C.assign(3, std::vector<double>(3, 0.0));
    moon_gravity_field.S.assign(3, std::vector<double>(3, 0.0));
    moon_gravity_field.C[2][0] = -lunar_j2 / sqrt(5.0);
    earth_job.opn_model.lunar_gravity_config = LUNAR_GRAVITY_JPL;
    earth_job.opn_model.lunar_zonal_order = 2;
    earth_job.opn_model.lunar_tesseral_degree = 0;

    // JPL 查询沿用近地右函数的 TT 近似，给解析参考式和两个状态方程相同的月球位置、姿态。
    aEOP lunar_eop;
    const double lunar_utc = TT2UTC(55256.0);
    earth_job.opn_EOP.Get_EOP(tracer, lunar_utc, lunar_eop, 1);
    if (tracer.Code() != ERR_NOERROR) return 1;

    const TimesMJD lunar_times = Utc2Else(lunar_eop.dut1, lunar_utc);
    double moon_state[6], omega, inclination, lambda;
    const int moon_status = Get_JPLEPH(lunar_times.TT, 10, 3, moon_state);
    const int libration_status = Get_Libration(lunar_times.TT, omega, inclination, lambda);

    if (moon_status < 0 || libration_status < 0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Lunar ephemeris or libration evaluation failed.", bastion::EventType::ERROR);
        return 1;
    }

    const Eigen::Matrix3d lunar_rotation = Rz(lambda) * Rx(inclination) * Rz(omega);
    const Eigen::Vector3d moon_position = Eigen::Vector3d(moon_state[0], moon_state[1], moon_state[2]) * constants::Scale::UNIT_LENGTH_M;
    const double acceleration_scale = constants::Scale::UNIT_GM_M3_S2 / pow(constants::Scale::UNIT_LENGTH_M, 2);
    double max_cart_error = 0.0, max_nsg2_error = 0.0;

    // 分别检查近地与高轨状态；高轨处中心引力较弱，能更清楚地检出微小间接项的遗漏。
    for (double semimajor_axis : {7000000.0, 200000000.0})
    {
        orbit::Keplerian lunar_kep;
        lunar_kep.Set(semimajor_axis / constants::Scale::UNIT_LENGTH_M, 0.1, 0.6, 0.4, 0.2, 0.8, true);
        orbit::NonSingular2 lunar_nsg2;
        lunar_kep.ToNonSingular2(tracer, lunar_nsg2);
        if (tracer.Code() != ERR_NOERROR) return 1;

        // 两条右函数使用同一 NonSingular2 对象还原出的状态，消除输入转换差。
        double lunar_cart[6], lunar_elements[6];
        lunar_nsg2.ToCartesian(tracer, lunar_cart);
        if (tracer.Code() != ERR_NOERROR) return 1;

        lunar_nsg2.ToArray(lunar_elements);
        const Eigen::Vector3d satellite_position = Eigen::Vector3d(lunar_cart[0], lunar_cart[1], lunar_cart[2]) * constants::Scale::UNIT_LENGTH_M;
        const Eigen::Vector3d fixed_positions[2] = {
            lunar_rotation * (satellite_position - moon_position), lunar_rotation * (-moon_position)};
        Eigen::Vector3d j2_acceleration[2];

        // 按有量纲 J2 闭式公式分别计算卫星处和地心处加速度，单位 m/s²。
        for (int point = 0; point < 2; ++point)
        {
            const Eigen::Vector3d& position = fixed_positions[point];
            const double r = position.norm();
            const double z2_r2 = position[2] * position[2] / (r * r);
            const double coefficient = 1.5 * moon_gravity_field.GM * lunar_j2 * pow(moon_gravity_field.R_ref, 2) / pow(r, 5);
            j2_acceleration[point] = coefficient * Eigen::Vector3d(position[0] * (5.0 * z2_r2 - 1.0),
                position[1] * (5.0 * z2_r2 - 1.0), position[2] * (5.0 * z2_r2 - 3.0));
        }

        const Eigen::Vector3d expected_relative = lunar_rotation.transpose() * (j2_acceleration[0] - j2_acceleration[1]);
        double without_lunar_nsg[6], with_lunar_nsg[6];
        earth_job.opn_model.use_lunar_NSG = false;
        Acceleration_Cart(tracer, earth_job, 55256.0, 6, lunar_cart, without_lunar_nsg);
        if (tracer.Code() != ERR_NOERROR) return 1;

        earth_job.opn_model.use_lunar_NSG = true;
        Acceleration_Cart(tracer, earth_job, 55256.0, 6, lunar_cart, with_lunar_nsg);
        if (tracer.Code() != ERR_NOERROR) return 1;

        for (int i = 0; i < 3; ++i)
        {
            const double actual = (with_lunar_nsg[i + 3] - without_lunar_nsg[i + 3]) * acceleration_scale;
            max_cart_error = std::max(max_cart_error, fabs(actual - expected_relative[i]));
        }

        // 将解析相对加速度加入点质量摄动，经现有高斯方程核对 NonSingular2 路径。
        double lunar_point_mass[3], expected_rates[6], actual_rates[6];
        Thirdbody_Gravity_Acceleration(lunar_point_mass, lunar_cart, moon_state, constants::Scale::NORMALIZED_MOON_GM_TT);
        const Eigen::Vector3d expected_perturbation = Eigen::Vector3d(lunar_point_mass[0], lunar_point_mass[1], lunar_point_mass[2])
            + expected_relative / acceleration_scale;
        Eigen::Matrix3d inertial_to_rtn;
        RTN_Matrix_from_Cart(lunar_cart[0], lunar_cart[1], lunar_cart[2], lunar_cart[3], lunar_cart[4], lunar_cart[5], inertial_to_rtn);
        const Eigen::Vector3d expected_rtn = inertial_to_rtn * expected_perturbation;
        GaussPlanetaryEquation_RTNtoNSG2(tracer, lunar_nsg2, expected_rtn.data(), expected_rates);
        if (tracer.Code() != ERR_NOERROR) return 1;

        Acceleration_NSG2(tracer, earth_job, 55256.0, 6, lunar_elements, actual_rates);
        if (tracer.Code() != ERR_NOERROR) return 1;

        for (int i = 0; i < 6; ++i) max_nsg2_error = std::max(max_nsg2_error, fabs(actual_rates[i] - expected_rates[i]));
    }

    std::cout << "Lunar relative J2: Cartesian error [m/s^2] = " << max_cart_error
              << ", NonSingular2 derivative error [normalized] = " << max_nsg2_error << '\n';

    if (max_cart_error > 3e-15 || max_nsg2_error > 5e-16)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Geocentric lunar J2 direct-minus-indirect mismatch.", bastion::EventType::ERROR);
        return 1;
    }


    std::cout << "TT/TDB gravity, GGM03C J2, lunar relative J2, solar radiation and drag checks passed.\n";
    tracer.CheckOut(__func__);
    return 0;
}
