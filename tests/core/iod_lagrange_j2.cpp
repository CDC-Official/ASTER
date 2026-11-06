/**
 * @file iod_lagrange_j2.cpp
 * @brief 用独立J2加速度检查GL位置时间展开的三次项。
 * @details 从标准Cartesian加速度作五点差分得到加速度的时间导数，
 * 再与现有Fz、Gz共同形成的位置三次项比较。检查覆盖二体及四个含J2分支，
 * 采用五组不同位置、速度，并包含赤道面内的零z分量对照。
 */
#include "core_iod.hpp"
#include "data_gravityfield.hpp"
#include <array>

namespace
{
using State = std::array<double, 6>;

/**
 * @brief 计算标准地球质点引力和可选J2引力的z分量。
 * @param[in] position 三维GCRS位置，除以Scale::UNIT_LENGTH_M。
 * @param[in] use_j2 是否加入GGM03C配套J2项。
 * @return z方向加速度，以Scale::UNIT_LENGTH_M/Scale::UNIT_TIME_S²归一化。
 */
double AccelerationZ(const Eigen::Vector3d& position, bool use_j2)
{
    const double radius = position.norm();
    double acceleration = -position.z() / std::pow(radius, 3);

    if (use_j2)
    {
        // 直接由标准轴对称J2加速度计算，不复用F、G的展开系数。
        const double j = 1.5 * constants::GGM03C::NORMALIZED_J2;
        acceleration += j * position.z() / std::pow(radius, 5)
                        * (5.0 * position.z() * position.z() / (radius * radius) - 3.0);
    }

    return acceleration;
}

/**
 * @brief 根据Fz、Gz取得某一相对时刻的z位置。
 * @param[in] state 六维GCRS参考状态，位置以UNIT_LENGTH_M、速度以UNIT_LENGTH_M/UNIT_TIME_S归一化。
 * @param[in] model IOD_DYNAMICS_*模型编号。
 * @param[in] tau 以UNIT_TIME_S归一化的相对参考历元时间。
 * @return 以UNIT_LENGTH_M归一化的z位置。
 */
double PositionZ(const State& state, int model, double tau)
{
    // 原接口接收非常量指针，因此用副本传入，所有时刻仍从同一初态展开。
    State copy = state;
    double coefficients[4];
    Lagrange_Coefficients(model, 60400.0, tau, copy.data(), coefficients);
    return coefficients[2] * state[2] + coefficients[3] * state[5];
}

/**
 * @brief 提取Fz、Gz共同形成的位置三次项系数。
 * @param[in] state 六维归一化GCRS参考状态，分量顺序为x、y、z、vx、vy、vz。
 * @param[in] model IOD_DYNAMICS_*模型编号。
 * @param[in] h 正的归一化时间步长。
 * @return z位置对归一化时间的三次项系数，等于第三阶时间导数的1/6。
 */
double CubicCoefficientZ(const State& state, int model, double h)
{
    // 正反时刻相减消去偶次项，再用两种步长消去五次项。
    const double small = (0.5 * (PositionZ(state, model, h) - PositionZ(state, model, -h))
                          - state[5] * h) / std::pow(h, 3);
    const double large = (0.5 * (PositionZ(state, model, 2.0 * h) - PositionZ(state, model, -2.0 * h))
                          - state[5] * 2.0 * h) / std::pow(2.0 * h, 3);
    return (4.0 * small - large) / 3.0;
}
}

int main()
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    const std::array<State, 5> states = {{
        {1.1, 0.0, 0.0, 0.03, 0.7, 0.4},
        {0.9, 0.4, 0.8, -0.25, 0.6, 0.2},
        {3.0, -1.0, 0.8, 0.1, 0.35, -0.2},
        {6.4, 0.4, -0.2, -0.03, 0.25, 0.14},
        {1.2, 0.1, 0.0, 0.05, 0.8, 0.0}
    }};
    double maximum_difference = 0.0;

    for (std::size_t index = 0; index < states.size(); ++index)
    {
        const State& state = states[index];
        const Eigen::Vector3d position(state[0], state[1], state[2]);
        const Eigen::Vector3d velocity(state[3], state[4], state[5]);

        for (int model = IOD_DYNAMICS_TWOBODY; model <= IOD_DYNAMICS_J2_SOLAR_LUNAR; ++model)
        {
            // 沿初始速度方向对独立加速度求导，得到本次需要的三阶位置导数。
            const bool use_j2 = model != IOD_DYNAMICS_TWOBODY;
            const double step = 0.001;
            const double jerk = (
                AccelerationZ(position - 2.0 * step * velocity, use_j2)
                - 8.0 * AccelerationZ(position - step * velocity, use_j2)
                + 8.0 * AccelerationZ(position + step * velocity, use_j2)
                - AccelerationZ(position + 2.0 * step * velocity, use_j2)) / (12.0 * step);

            // 四个摄动分支共用同一J2的Fz/Gz三次展开，按其合成位置检查。
            const double coefficient = CubicCoefficientZ(state, model, 0.02);
            const double difference = std::abs(6.0 * coefficient - jerk);
            maximum_difference = std::max(maximum_difference, difference);

            if (!std::isfinite(difference) || difference > 2e-8)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Wrong cubic z coefficient: state=" + std::to_string(index)
                    + ", model=" + std::to_string(model) + ", difference=" + std::to_string(difference), bastion::ERROR);
                tracer.ReportEvents(std::cerr);
                return 1;
            }
        }
    }

    std::cout << std::setprecision(17) << "25 cubic-coefficient checks passed; maximum normalized jerk difference = "
              << maximum_difference << '\n';
    tracer.CheckOut(__func__);
    return 0;
}
