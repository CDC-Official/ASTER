/**
 * @file analyt_terms.hpp
 * @brief 瞬转平所用的 J2、J22、坐标架、日月及光压解析短周期公式接口。
 * @details 各公式接收当前均根数估计的轨道几何，覆写短周期改正，由 OsculatingToMean 组织迭代和累计。
 * 长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad；输入相角均已由调用方按历元准备。
 * NSG2 六槽顺序为 [a,xi,eta,h,k,lambda]，xi=e*cos(omega+Omega)、eta=e*sin(omega+Omega)，
 * h=sin(i/2)*cos(Omega)、k=sin(i/2)*sin(Omega)、lambda=M+omega+Omega。
 */
#ifndef CDC_ANALYT_TERMS_HPP
#define CDC_ANALYT_TERMS_HPP

#include "utils_overwatch.hpp"
#include <vector>

namespace orbit::averaging_detail
{
constexpr int RadiationFull = 11, RadiationCombined = 12;

/** @brief 简化岁差及章动主项，均为 TT 历元下的 rad 角度。 */
struct FrameArguments
{
    double thetaA;   // 线性岁差角。
    double Dtheta;   // 简化章动角主项。
    double Depsilon; // 黄赤交角章动主项。
};

/**
 * @brief 当前第三体的位置尺度及相对于卫星轨道的方向三角量。
 * @details 每个 UTC 历元以 UT1=UTC 准备日月自身量；每轮均根数估计更新含卫星角度的组合量。
 */
struct ThirdBodyArguments
{
    double r = 0; // 第三体地心距离，以 ae 归一化。
    double beta2 = 0, beta3 = 0, beta4 = 0; // 分别为 (GM'/GE)/r^3、/r^4、/r^5。
    double Omega = 0, u = 0; // 第三体地心 J2000 赤道轨道的升交点赤经、纬度幅角，rad。
    double n = 0; // 第三体平运动，单位为 1/UNITTIME；双快变量公式使用月球值。
    double si = 0, cip1 = 0, cim1 = 0; // sin(i')、1+cos(i')、1-cos(i')。

    // theta=Omega_sat-Omega；pp/pm/mp/mm 依次表示 omega+u+theta、omega-u+theta、omega+u-theta、omega-u-theta。
    double s_pp = 0, c_pp = 0, s_pm = 0, c_pm = 0;
    double s_mp = 0, c_mp = 0, s_mm = 0, c_mm = 0;
    double s_p = 0, c_p = 0, s_m = 0, c_m = 0; // sin/cos(omega+u)、sin/cos(omega-u)。
};

/**
 * @brief 卫星与第三体方向投影的三个无量纲乘积 A^2、B^2、A*B。
 * @details ThirdBodyShort 按当前卫星估计和第三体方向写入；ThirdBodyDoubleShort 只读它们，计算四次半长轴项。
 * 两次调用须使用同一轮卫星几何和同一个第三体；双快变量结果覆写单快变量六槽，调用方最终只累计一次。
 */
struct ThirdBodyProducts
{
    double AA = 0, BB = 0, AB = 0;
};

/**
 * @brief 计算 J2 一阶 NSG2 短周期项及独立的二阶半长轴项。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] h 当前 NSG2 分量 sin(i/2)*cos(Omega)，无量纲。
 * @param[in] k 当前 NSG2 分量 sin(i/2)*sin(Omega)，无量纲。
 * @param[in] jp 1.5*J2/p^2，p=a*(1-e^2) 为 ae 归一化半通径。
 * @param[in] ja 1.5*J2/a，a 为 ae 归一化半长轴。
 * @param[in] a 当前均根数估计的半长轴，以 ae 归一化。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] f1e 近圆组合因子 1/(1+sqrt(1-e^2))。
 * @param[in] f2e 近圆组合因子 (1+2*sqrt(1-e^2))/(1+sqrt(1-e^2))^2。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] w 当前卫星近地点幅角 omega，rad。
 * @param[in] wo 当前卫星近地点经度 omega+Omega，rad。
 * @param[in] f 当前卫星真近点角，rad。
 * @param[in] ar 当前卫星半长轴与地心距离之比 a/r，无量纲。
 * @param[in] fmesf 方程中心 f-M 回绕到 [-pi,pi) 后加 e*sin(f)，rad。
 * @param[out] nsg2s1_J2 覆写六个 NSG2 短周期改正 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @param[out] as2_J2 覆写 J2 二阶半长轴短周期项，以 ae 归一化，由调用方另行扣除一次。
 * @details 来自 analyt_j2.cpp 的 AnaJ2ShortNSG2。六槽和二阶半长轴分别写出，调用方各累计或扣除一次。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void J2Short(bastion::Overwatch& tracer, double h, double k, double jp, double ja, double a,
    double e, double f1e, double f2e, double si, double ci, double w, double wo, double f, double ar, double fmesf,
    double nsg2s1_J2[], double& as2_J2);

/**
 * @brief 计算 J2 二阶平经度短周期项。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] j2p4 J2^2/p^4，p 为 ae 归一化半通径，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] si2 sin(i)^2，无量纲。
 * @param[in] si4 sin(i)^4，无量纲。
 * @param[in] s2f2w sin(2*f+2*omega)，无量纲。
 * @param[in] s4f4w sin(4*f+4*omega)，无量纲。
 * @param[out] ls2_J2 覆写 J2 二阶平经度短周期项，rad，由调用方另行扣除一次。
 * @details 独立输出 delta_lambda，不写入 J2 一阶六槽。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void J2SecondOrderLongitude(bastion::Overwatch& tracer, double j2p4, double ci, double si2, double si4, double s2f2w, double s4f4w, double& ls2_J2);

/**
 * @brief 计算用于替换 J22 半长轴槽的高阶短周期补充项。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] j22a 未归一化 J22 幅值除以 a，a 以 ae 归一化。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] o22 Omega 减格林尼治角及 J22 系数半相位之和，rad。
 * @param[in] w 当前卫星近地点幅角 omega，rad。
 * @param[in] m 当前卫星平近点角 M，rad。
 * @param[in] alpha 地球自转角速度与卫星平运动之比 omega_E/n，无量纲。
 * @param[out] as2_J22 覆写 J22 高阶半长轴短周期项，以 ae 归一化，用于替换 J22Short 输出的半长轴槽。
 * @param[in] resonance_id 由密切输入固定的共振类别；1/2/3 为 alpha 约等于 1/2、1、3/2，其余类别按调用方的频率表解释，-1 为未命中。
 * @details 来自 analyt_j22.cpp 的 AnaJ22ShortSMAExpansion。本项替换 J22Short 的第 0 槽，不与该槽重复相加。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void J22SemimajorShort(bastion::Overwatch& tracer, double j22a, double e, double ci, double o22, double w, double m, double alpha, double& as2_J22, int resonance_id);

/**
 * @brief 将 J22 组合中间项转换为 NSG2 短周期改正。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] h 当前 NSG2 分量 sin(i/2)*cos(Omega)，无量纲。
 * @param[in] k 当前 NSG2 分量 sin(i/2)*sin(Omega)，无量纲。
 * @param[in] j22a2 未归一化 J22 幅值除以 a^2，a 以 ae 归一化。
 * @param[in] a 当前均根数估计的半长轴，以 ae 归一化。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] cih 当前卫星轨道半倾角的 cos(i/2)，无量纲；要求 i<pi。
 * @param[in] so 当前卫星升交点赤经的 sin(Omega)，无量纲。
 * @param[in] co 当前卫星升交点赤经的 cos(Omega)，无量纲。
 * @param[in] o22 Omega 减格林尼治角及 J22 系数半相位之和，rad。
 * @param[in] swo 当前卫星近地点经度的 sin(omega+Omega)，无量纲。
 * @param[in] cwo 当前卫星近地点经度的 cos(omega+Omega)，无量纲。
 * @param[in] w 当前卫星近地点幅角 omega，rad。
 * @param[in] m 当前卫星平近点角 M，rad。
 * @param[in] alpha 地球自转角速度与卫星平运动之比 omega_E/n，无量纲。
 * @param[out] nsg2s2_J22 覆写六个 NSG2 短周期改正 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @param[in] resonance_id 由密切输入固定的共振类别；1/2/3 为 alpha 约等于 1/2、1、3/2，其余类别按调用方的频率表解释，-1 为未命中。
 * @details 来自 analyt_j22.cpp 的 AnaJ22ShortNSG2。内部使用 Combined 组合，偏心率方向取 omega+Omega；半长轴可随后由 J22SemimajorShort 替换。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void J22Short(bastion::Overwatch& tracer, double h, double k, double j22a2, double a, double e,
    double si, double ci, double cih, double so, double co, double o22, double swo, double cwo, double w, double m, double alpha,
    double nsg2s2_J22[], int resonance_id);

/**
 * @brief 将岁差章动与 J2 耦合的组合短周期项转换为 NSG2。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] frameArg 当前 TT 历元的简化岁差角 thetaA 和章动主项 Dtheta、Depsilon，均为 rad。
 * @param[in] ja 3*J2/a，a 为 ae 归一化半长轴。
 * @param[in] jp2 1.5*J2/p^2，p 为 ae 归一化半通径。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] xe2 sqrt(1-e^2)，无量纲。
 * @param[in] xe3 (1-e^2)^(-3/2)，无量纲。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] sih 当前卫星轨道半倾角的 sin(i/2)，无量纲。
 * @param[in] cih 当前卫星轨道半倾角的 cos(i/2)，无量纲；要求 i<pi。
 * @param[in] so 当前卫星升交点赤经的 sin(Omega)，无量纲。
 * @param[in] co 当前卫星升交点赤经的 cos(Omega)，无量纲。
 * @param[in] s2w sin(2*omega)，无量纲。
 * @param[in] c2w cos(2*omega)，无量纲。
 * @param[in] swo 当前卫星近地点经度的 sin(omega+Omega)，无量纲。
 * @param[in] cwo 当前卫星近地点经度的 cos(omega+Omega)，无量纲。
 * @param[in] sf sin(f)，无量纲。
 * @param[in] cf cos(f)，无量纲。
 * @param[in] fm 回绕到 [-pi,pi) 的方程中心 f-M，rad。
 * @param[in] fmef fm+e*sin(f)，rad。
 * @param[in] ar3 当前卫星距离比 (a/r)^3，无量纲。
 * @param[out] nsg2s2_frameJ2 覆写六个 NSG2 短周期改正 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @details 来自 analyt_frameinduced.cpp 的 AnaFrameJ2ShortNSG2。内部使用 Combined 的缩放节点和近地点组合，按半倾角关系构造 h/k。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void FrameShort(bastion::Overwatch& tracer, const FrameArguments& frameArg, double ja, double jp2,
    double e, double xe2, double xe3, double si, double ci, double sih, double cih,
    double so, double co, double s2w, double c2w, double swo, double cwo, double sf, double cf, double fm, double fmef, double ar3,
    double nsg2s2_frameJ2[]);

/**
 * @brief 计算第三体单快变量 NSG2 短周期项及当前方向乘积。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] apow 当前卫星归一化半长轴幂表 apow[j]=a^j，含下标 0..10，零槽为一。
 * @param[in] ra 当前卫星地心距离与半长轴之比 r/a，无量纲。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] je2 1-e^2，无量纲。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] sih 当前卫星轨道半倾角的 sin(i/2)，无量纲。
 * @param[in] cih 当前卫星轨道半倾角的 cos(i/2)，无量纲；要求 i<pi。
 * @param[in] so 当前卫星升交点赤经的 sin(Omega)，无量纲。
 * @param[in] co 当前卫星升交点赤经的 cos(Omega)，无量纲。
 * @param[in] swo 当前卫星近地点经度的 sin(omega+Omega)，无量纲。
 * @param[in] cwo 当前卫星近地点经度的 cos(omega+Omega)，无量纲。
 * @param[in] sm 平近点角的倍角表 sin(j*M)，含下标 0..4；零槽为零。
 * @param[in] cm 平近点角的倍角表 cos(j*M)，含下标 0..4；零槽为一。
 * @param[in] sf 真近点角的倍角表 sin(j*f)，含下标 0..4；零槽为零。
 * @param[in] cf 真近点角的倍角表 cos(j*f)，含下标 0..4；零槽为一。
 * @param[in] se 偏近点角的倍角表 sin(j*E)，含下标 0..4；零槽为零。
 * @param[in] ce 偏近点角的倍角表 cos(j*E)，含下标 0..4；零槽为一。
 * @param[in] lspar 当前第三体的 beta2/beta3/beta4、轨道方向及卫星组合相角；地心距离以 ae 归一化，角度 rad，平运动 n 以 1/UNITTIME 表示。
 * @param[out] products 覆写当前方向投影的 AA=A^2、BB=B^2、AB=A*B，无量纲，供同一配置的双快变量四次半长轴项使用。
 * @param[out] nsg2s2_lunisolar 覆写六个 NSG2 短周期改正 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @details 来自 analyt_lunisolar.cpp 的 AnaLunisolarShortNSG2。包含 Legendre 二、三次六槽及四次半长轴项；双快变量分支仍需本函数先准备相同配置的三个乘积。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void ThirdBodyShort(bastion::Overwatch& tracer, const std::vector<double>& apow, double ra,
    double e, double je2, double si, double ci, double sih, double cih, double so, double co, double swo, double cwo,
    const double sm[], const double cm[], const double sf[], const double cf[],
    const double se[], const double ce[], const ThirdBodyArguments& lspar, ThirdBodyProducts& products,
    double nsg2s2_lunisolar[]);

/**
 * @brief 计算卫星和第三体双快变量 NSG2 短周期项。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] apow 当前卫星归一化半长轴幂表 apow[j]=a^j，含下标 0..10，零槽为一。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] je2 1-e^2，无量纲。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] sih 当前卫星轨道半倾角的 sin(i/2)，无量纲。
 * @param[in] cih 当前卫星轨道半倾角的 cos(i/2)，无量纲；要求 i<pi。
 * @param[in] so 当前卫星升交点赤经的 sin(Omega)，无量纲。
 * @param[in] co 当前卫星升交点赤经的 cos(Omega)，无量纲。
 * @param[in] swo 当前卫星近地点经度的 sin(omega+Omega)，无量纲。
 * @param[in] cwo 当前卫星近地点经度的 cos(omega+Omega)，无量纲。
 * @param[in] sm 平近点角的倍角表 sin(j*M)，含下标 0..4；零槽为零。
 * @param[in] cm 平近点角的倍角表 cos(j*M)，含下标 0..4；零槽为一。
 * @param[in] lspar 当前第三体的 beta2/beta3/beta4、轨道方向及卫星组合相角；地心距离以 ae 归一化，角度 rad，平运动 n 以 1/UNITTIME 表示。
 * @param[in] products 由同一轮卫星估计和同一第三体的 ThirdBodyShort 准备的 AA、BB、AB，仅用于四次半长轴项；本函数保持其值不变。
 * @param[out] nsg2s2_lunisolar 覆写六个 NSG2 短周期改正 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @details 来自 analyt_lunisolar.cpp 的 AnaLunisolarShortNSG2_DoubleFastArg。读取月球平运动 n 和当前几何乘积，覆写六槽以替代单快变量结果；调用方只累计本结果一次。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
void ThirdBodyDoubleShort(bastion::Overwatch& tracer, const std::vector<double>& apow, double e, double je2,
    double si, double ci, double sih, double cih, double so, double co, double swo, double cwo,
    const double sm[], const double cm[], const ThirdBodyArguments& lspar, const ThirdBodyProducts& products,
    double nsg2s2_lunisolar[]);

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
    double nsg2s2_srp[], int short_flag);
} // namespace orbit::averaging_detail
#endif
