/**
 * @file analyt_solar_radiation.cpp
 * @brief 太阳光压无地影分支的解析短周期公式。
 * @details 输入包含 Cr、归一化面积质量比和太阳距离的幅值；先计算组合角改正，再转换为 NSG2 六槽。
 * 长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad；输出为短周期改正，不是状态时间导数。
 */

// 用于瞬转平的解析短周期公式的定向移植与本库适配，原公式作者 Jingshi Tang。
// Copyright (c) 2020 国家航天局空间碎片监测与应用中心。All rights reserved.
// 仅保留 O2M 实际使用的代数项；轨道、时间、重力及迭代流程由 CDC_lib 提供。
#include "analyt_terms.hpp"
#include <cmath>
#include <array>

namespace orbit::averaging_detail
{
/**
 * @brief 计算无地影光压的组合短周期中间项。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] beta1 包含 Cr、归一化面积质量比及太阳距离因子的光压幅值，按归一化动力学尺度准备；物理 AMR 不直接传入。
 * @param[in] apow 当前卫星归一化半长轴幂表 apow[j]=a^j，含下标 0..10，零槽为一。
 * @param[in] epow 当前卫星偏心率幂表 epow[j]=e^j，含下标 0..5，零槽为一。
 * @param[in] xe2 sqrt(1-e^2)，无量纲。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] se 偏近点角的倍角表 sin(j*E)，含下标 0..4；零槽为零。
 * @param[in] ce 偏近点角的倍角表 cos(j*E)，含下标 0..4；零槽为一。
 * @param[in] A 卫星轨道面与太阳方向的余弦投影组合，无量纲，由本次变换首轮地影几何确定。
 * @param[in] B 卫星轨道面与太阳方向的正弦投影组合，无量纲，与 A 使用同一首轮几何。
 * @param[in] solarPar 当前太阳方向的 si/cip1/cim1 及与卫星角度组合的正余弦；对应当前均根数估计。
 * @param[out] keps2_srp 覆写六槽；前 3 槽为 delta_a、delta_e、delta_i，后 3 槽为 sin(i)*delta_Omega、e*(delta_omega+delta_Omega)、delta_lambda，半长轴尺度 ae，角度 rad。
 * @param[in] short_flag RadiationFull=11 或 RadiationCombined=12；瞬转平采用 Combined。
 * @details 来自 analyt_radiationpressure.cpp 的 AnaSRPShortKep。Full 分支增加 e/2 并记录尚未修订的平近点角项警告，随后贯穿到 Combined 的角组合；两者都不能按普通 Kepler 六槽解释。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
static void RadiationKeplerShort(bastion::Overwatch& tracer, double beta1, const std::vector<double>& apow, const std::vector<double>& epow, double xe2,
    double si, double ci, const double* se, const double* ce, double A, double B,
    const ThirdBodyArguments& solarPar, double keps2_srp[], int short_flag)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 光压分支支持 Full 和 Combined；Full 先调整偏近点角项，再进入相同的组合输出。
    if (short_flag != RadiationFull && short_flag != RadiationCombined)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown flag for SRP short-period terms!", bastion::ERROR);
        return;
    }

    // H1/H2 将当前卫星倾角和太阳方向组合，A/B 由外层地影几何步骤传入。
    double H1, H2;

    // 两方向分别采用太阳组合相角的余弦和正弦。
    H1 = si * (solarPar.cip1 * (solarPar.c_mp - solarPar.c_pm) + solarPar.cim1 * (solarPar.c_mm - solarPar.c_pp)) +
         2.0 * ci * solarPar.si * (solarPar.c_m - solarPar.c_p);
    H2 = si * (solarPar.cip1 * (solarPar.s_mp - solarPar.s_pm) + solarPar.cim1 * (solarPar.s_mm - solarPar.s_pp)) +
         2.0 * ci * solarPar.si * (solarPar.s_m - solarPar.s_p);

    // se/ce[j] 是 jE 的正余弦；cee 先取 cos(E)。
    double cee;
    cee = ce[1];

    // Full 分支在 cee 中增加 e/2。
    if (short_flag == RadiationFull)
    {
        cee += epow[1] / 2.0;
    }

    // beta1 已包含 Cr、归一化面积质量比和太阳距离因子，再乘 a² 得 ba2。
    double ba2, s[6];
    ba2 = beta1 * apow[2];

    // 先计算半长轴和偏心率相关项，使用一、二倍 E 的三角表。
    s[0] = 2.0 * ba2 * apow[1] * (A * cee + B * xe2 * se[1]);
    s[1] = ba2 * ((1.0 - epow[2]) * A * ce[2] / 4.0 + xe2 * B * (-epow[1] * se[1] / 2.0 + se[2] / 4.0));

    // 倾角相关项由 H1/H2 与偏近点角谐波组合，并含原 xe2 缩放。
    s[2] = ba2 / 4.0 * ((cee - epow[1] * ce[2] / 4.0) * H1 -
                        ((1.0 - epow[2] / 2.0) * se[1] - epow[1] * se[2] / 4.0) * H2 / xe2);

    // 节点相关中间槽仍保留 sin(i) 缩放，后面的 NSG2 包装直接组合它。
    s[3] = ba2 / 4.0 * (((1.0 - epow[2] / 2.0) * se[1] - epow[1] * se[2] / 4.0) * H1 / xe2 +
                        (cee - epow[1] * ce[2] / 4.0) * H2);

    // 这一槽与偏心率和节点槽组合为 e*(delta_omega+delta_Omega)。
    s[4] = ba2 * (xe2 * A * (-epow[1] * se[1] / 2.0 + se[2] / 4.0) + B * (epow[1] * cee - ce[2] / 4.0));

    // 平经度相关周期项使用相同 A/B 方向，但幅值为 -5*ba2。
    s[5] = -5.0 * ba2 * (A * ((1.0 - epow[2] / 2.0) * se[1] - epow[1] * se[2] / 4.0) +
                         B * xe2 * (-ce[1] + epow[1] * ce[2] / 4.0));

    // tih=sin(i)/(1+cos(i))=tan(i/2)，用于末三槽组合。
    double tih;
    tih = si / (1.0 + ci);

    // 前三槽先传递半长轴、偏心率及倾角相关结果。
    keps2_srp[0] = s[0];
    keps2_srp[1] = s[1];
    keps2_srp[2] = s[2];

    // Full 报告平近点角项尚未修订，再贯穿到 Combined；输出为组合六槽。
    switch (short_flag)
    {
        // 原实现指出该分支无地影平近点角项尚未修订，警告后继续下面组合。
        case RadiationFull:
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Short-period term for mean anomaly WITHOUT eclipse is NOT revised yet.", bastion::WARNING);

            // Combined 输出用于 NSG2；这里没有新增地影边界计算。
        case RadiationCombined:

            // 下标 3 保留 sin(i)*delta_Omega，外层按半倾角公式将它组合为 h/k 改正。
            keps2_srp[3] = s[3];

            // 此槽已是 e*(delta_omega+delta_Omega)，不再作为普通 delta_omega 解释。
            keps2_srp[4] = epow[1] * tih * s[3] + s[4];

            // 末槽合成 delta_lambda，避免最终转换时重新除以 e。
            keps2_srp[5] = tih * s[3] + epow[1] / (1.0 + xe2) * s[4] + s[5];
            break;

            // 其他分支标志通过 tracer 报错。
        default:
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown flag for SRP short-period terms!", bastion::ERROR);
            return;
            break;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 将无地影光压组合项转换为 NSG2 短周期改正。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] beta1 包含 Cr、归一化面积质量比及太阳距离因子的光压幅值，按归一化动力学尺度准备；物理 AMR 不直接传入。
 * @param[in] apow 当前卫星归一化半长轴幂表 apow[j]=a^j，含下标 0..10，零槽为一。
 * @param[in] epow 当前卫星偏心率幂表 epow[j]=e^j，含下标 0..5，零槽为一。
 * @param[in] xe2 sqrt(1-e^2)，无量纲。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] cih 当前卫星轨道半倾角的 cos(i/2)，无量纲；要求 i<pi。
 * @param[in] so 当前卫星升交点赤经的 sin(Omega)，无量纲。
 * @param[in] co 当前卫星升交点赤经的 cos(Omega)，无量纲。
 * @param[in] swo 当前卫星近地点经度的 sin(omega+Omega)，无量纲。
 * @param[in] cwo 当前卫星近地点经度的 cos(omega+Omega)，无量纲。
 * @param[in] se 偏近点角的倍角表 sin(j*E)，含下标 0..4；零槽为零。
 * @param[in] ce 偏近点角的倍角表 cos(j*E)，含下标 0..4；零槽为一。
 * @param[in] A 卫星轨道面与太阳方向的余弦投影组合，无量纲，由本次变换首轮地影几何确定。
 * @param[in] B 卫星轨道面与太阳方向的正弦投影组合，无量纲，与 A 使用同一首轮几何。
 * @param[in] solarPar 当前太阳方向的 si/cip1/cim1 及与卫星角度组合的正余弦；对应当前均根数估计。
 * @param[out] nsg2s2_srp 覆写六个 NSG2 短周期改正 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @param[in] short_flag RadiationFull=11 或 RadiationCombined=12；瞬转平采用 Combined，有地影时由外层跳过本函数。
 * @details 来自 analyt_radiationpressure.cpp 的 AnaSRPDShortNSG2。几何投影 A/B 在本次瞬转平首轮固定，其余太阳组合角按当前均根数估计更新。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void RadiationShort(bastion::Overwatch& tracer, double beta1, const std::vector<double>& apow, const std::vector<double>& epow, double xe2,
    double si, double ci, double cih, double so, double co, double swo, double cwo,
    const double* se, const double* ce, double A, double B, const ThirdBodyArguments& solarPar,
    double nsg2s2_srp[], int short_flag)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 先计算光压组合中间槽；输出为短周期改正，六槽用于组合 NSG2。
    double keps2_srp[6];
    RadiationKeplerShort(tracer, beta1, apow, epow, xe2, si, ci, se, ce, A, B, solarPar, keps2_srp, short_flag);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 前三槽为 delta_a、delta_xi、delta_eta，偏心率向量以 omega+Omega 定向。
    nsg2s2_srp[0] = keps2_srp[0];
    nsg2s2_srp[1] = cwo * keps2_srp[1] - swo * keps2_srp[4];
    nsg2s2_srp[2] = swo * keps2_srp[1] + cwo * keps2_srp[4];

    // h/k 将倾角改正与缩放节点改正组合，cos(i/2) 分母仍使用对应公式。
    nsg2s2_srp[3] = 0.5 * cih * co * keps2_srp[2] - 0.5 / cih * so * keps2_srp[3];
    nsg2s2_srp[4] = 0.5 * cih * so * keps2_srp[2] + 0.5 / cih * co * keps2_srp[3];

    // 末槽已经组合为 delta_lambda，直接传递。
    nsg2s2_srp[5] = keps2_srp[5];

    tracer.CheckOut(__func__);
}
} // namespace orbit::averaging_detail
