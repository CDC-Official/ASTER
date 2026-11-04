/**
 * @file orbit_stm.hpp
 * @brief 六状态 NSG2 解析状态转移矩阵及其近似设置。
 */
#pragma once
#include "orbit_basic.hpp"
#include "data_gravityfield.hpp"
#include "model_ephemeris_analytic.hpp"

namespace determination
{
struct AnalyticStmSetting
{
    // 二体、J2 长期项和日月平均项分别由设置启用。
    bool use_J2 = false, use_lunisolar = false;

    // J2 势项在程序尺度下的强度，包含配套 GM 比和参考半径平方比。
    double j2 = constants::GGM03C::NORMALIZED_J2;

    // 日月 GM 与距离都换算到积分尺度，组成平均项中的 GM/r³。
    double lunar_mass = constants::Scale::NORMALIZED_MOON_GM_TT, solar_mass = constants::Scale::NORMALIZED_SUN_GM_TT;
    double lunar_distance = constants::MeanOrbit::LUNAR_SEMIMAJOR_AXIS / constants::Scale::UNIT_LENGTH_M;
    double solar_distance = constants::MeanOrbit::SOLAR_SEMIMAJOR_AXIS / constants::Scale::UNIT_LENGTH_M;

    // 此平均近似将日月轨道面取为平黄道面，倾角采用同一个 J2000 平黄赤交角。
    double lunar_sin_i = sin(constants::MeanOrbit::OBLIQUITY), solar_sin_i = sin(constants::MeanOrbit::OBLIQUITY);
};

StateMatrix AnalyticStmNSG2(const StateVector &state, double elapsed_norm, const AnalyticStmSetting &setting);
} // namespace determination
