/**
 * @file orbit_stm.hpp
 * @brief 六状态 NSG2 解析状态转移矩阵及其近似设置。
 */
#pragma once
#include "orbit_basic.hpp"

namespace determination
{
struct AnalyticStmSetting
{
    // 导师解析近似：二体 + J2 长期项 + 日月平均项；不等于完整动力学变分方程。
    // 默认仅二体；与导师默认对比时应明确打开两项，不能暗中跟随力模型。
    bool use_J2 = false, use_lunisolar = false;
    // j2 是无量纲的地球 J2 系数，仅在启用 J2 长期项时使用。
    double j2 = 0.00108263;
    // 名称沿用来源公式；实际保存日月 GM 除以地球 GE 的比值，不是以 kg 表示的质量。
    double lunar_mass = 4902.7934e9 / GE, solar_mass = 1.32712438e20 / GE;
    // 平均距离除以地球半径 ae，成为归一化长度；与上面的 GM 比值组成 GM/r³ 强度。
    double lunar_distance = 384747981.0 / ae, solar_distance = 1.495980225898274e11 / ae;
    // 这里保存日月轨道倾角的正弦，不是角度；作为平均长期项的固定近似，不随历元更新。
    double lunar_sin_i = 0.39777697, solar_sin_i = 0.39777697;
};

StateMatrix AnalyticStmNSG2(const StateVector &state, double elapsed_norm, const AnalyticStmSetting &setting);
} // namespace determination
