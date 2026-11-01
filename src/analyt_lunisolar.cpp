/**
 * @file analyt_lunisolar.cpp
 * @brief 第三体引力的单快变量和双快变量解析短周期公式。
 * @details 单快变量输出当前几何乘积；双快变量读取同一配置的三个乘积并覆写 NSG2 六槽，调用方只累计选定结果一次。
 * 长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad；输出为短周期改正，不是状态时间导数。
 */

// 用于瞬转平的解析短周期公式的定向移植与本库适配，原公式作者 Jingshi Tang。
// Copyright (c) 2020 国家航天局空间碎片监测与应用中心。All rights reserved.
// 仅保留 O2M 实际使用的代数项；轨道、时间、重力及迭代流程由 CDC_lib 提供。
#include "analyt_terms.hpp"
#include <cmath>
#include <array>

namespace orbit::averaging_detail {
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
    double nsg2s2_lunisolar[])
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // s2/s3/s4 分别保存 Legendre 二、三、四次展开贡献，六槽暂为组合中间量。
    double s2[6] = {0.0}, s3[6] = {0.0}, s4[6] = {0.0};

    // apow[j]=a^j，beta2/beta3/beta4 已含第三体归一化 GM 与距离幂。
    double a    = apow[1];
    double ba3  = 1.5*lspar.beta2*apow[3];

    // 高一次 Legendre 项多一个 a/第三体距离因子，这里保留各阶幅值。
    double ba4  = 1.5*lspar.beta3*apow[4];
    double ba5  = 1.5*lspar.beta4*apow[5];

    // tih=tan(i/2)，cip1=1+cos(i)，都用于不直接除以 sin(i) 的方向组合。
    double tih  = sih/cih;
    double cip1 = 1.0+ci;

    // je2=1-e²，先从它恢复 e²、e³、e⁴，再求 sqrt(1-e²)。
    double e2   = 1.0-je2;
    double e3   = e2*e;
    double e4   = e2*e2;

    // xe2 供角变量的无奇异组合使用，不改变输入的归一化尺度。
    double xe2  = sqrt(je2);

    // 第三体 1+cos(i') 与组合相角余弦相乘，i' 为第三体轨道倾角。
    double cp_cmp = lspar.cip1 * lspar.c_mp;
    double cp_cpm = lspar.cip1 * lspar.c_pm;

    // 同一 1+cos(i') 权重的正弦组合随后用于 B 方向量。
    double cp_smp = lspar.cip1 * lspar.s_mp;
    double cp_spm = lspar.cip1 * lspar.s_pm;

    // 另一组采用第三体 1-cos(i')，保留原来的和差角槽位。
    double cm_cmm = lspar.cim1 * lspar.c_mm;
    double cm_cpp = lspar.cim1 * lspar.c_pp;

    // 这一对正弦量与上面余弦量共同描述相反节点方向。
    double cm_smm = lspar.cim1 * lspar.s_mm;
    double cm_spp = lspar.cim1 * lspar.s_pp;

    // 不含节点差的余弦组合乘第三体 sin(i')，补充倾斜轨道方向分量。
    double s_cm = lspar.si * lspar.c_m;
    double s_cp = lspar.si * lspar.c_p;

    // 对应的正弦组合使用相同 sin(i') 权重。
    double s_sm = lspar.si * lspar.s_m;
    double s_sp = lspar.si * lspar.s_p;

    // A/B 及带数字后缀的量是卫星-第三体方向投影的分组中间量。
    double A0, B0, dA, dB, A, B, A20, dA2, A2, B20, dB2, B2, A1, B1;

    // 后缀 _si 表示该部分尚未乘入卫星 sin(i)，用于避免后面再除以 sin(i)。
    double A0_si, B0_si, A20_si, B20_si;

    // 先保留方向平方、乘积及其方向变化组合，后面按阶展开共用。
    double AA,BB,AB,AB1,A1B,AA1,BB1;

    // 用半倾角与第三体方向先构造 A0、B0 中未乘 sin(i) 的部分。
    A0_si =  0.25*(tih*(cp_cmp+cm_cmm) + 2.0*(s_cm-s_cp));
    B0_si = -0.25*(tih*(cp_smp+cm_smm) + 2.0*(s_sm-s_sp));

    // 补上卫星 sin(i)，得到这一组实际方向投影。
    A0 =  si*A0_si;
    B0 =  si*B0_si;

    // 另一组投影采用 1+cos(i) 权重，随后与 A0/B0 相加。
    dA =  0.25*cip1*(cp_cpm + cm_cpp);
    dB = -0.25*cip1*(cp_spm + cm_spp);

    // A、B 是两部分投影之和；后面的 AA、BB、AB 都由这一组当前几何得到。
    A  =  A0+dA;
    B  =  B0+dB;

    // A20/B20 同样先保留未乘 sin(i) 的形式。
    A20_si = 0.25*tih*(cp_smp + cm_smm);
    B20_si = 0.25*tih*(cp_cmp + cm_cmm);

    // 补上 sin(i) 后再与 dA2/dB2 合成 A2/B2。
    A20 = si*A20_si;
    B20 = si*B20_si;

    // 两支投影满足原来的交叉关系 dA2=dB、dB2=-dA。
    dA2 =  dB;
    dB2 = -dA;

    // 合成用于节点相关中间式的 A2、B2。
    A2 = A20+dA2;
    B2 = B20+dB2;

    // A1/B1 为当前卫星倾角方向变化的组合，后面与 A/B 构成交叉乘积。
    A1 =  0.25*(si*(cp_cmp + cm_cmm - cp_cpm - cm_cpp) + 2.0*ci*(s_cm - s_cp));
    B1 = -0.25*(si*(cp_smp + cm_smm - cp_spm - cm_spp) + 2.0*ci*(s_sm - s_sp));

    // AA=A²、BB=B²，是 Legendre 二次部分的平方方向项。
    AA  = A*A;
    BB  = B*B;

    // AB 保存两方向交叉项，AB1 保存 B1 与 A 的乘积。
    AB  = A*B;
    AB1 = A*B1;

    // 其余乘积继续保留 A1/B1 对方向变化的影响。
    A1B = A1*B;
    AA1 = A*A1;

    // BB1 完成这组一阶方向乘积，与 AA1 成对使用。
    BB1 = B*B1;

    // Is1..Is3 是倾角改正所需的组合，保留重排后无显式 sin(i) 分母的原式。
    double Is1, Is2, Is3;

    Is1 = -0.5*(A0_si*B2 + A20_si*B  + B0_si*dB + B20_si*dA + ci*(A0_si*(A+dA)-B0_si*(B+dB)))+0.5*tih*(dA*dA-dB*dB);

    // 其余两个组合仍读取未乘 sin(i) 的中间量，避免重新引入小倾角除法。
    Is2 =  0.5*(A0_si*A2 + A20_si*dA - B0_si*B2 - B20_si*dB) - ci*(A0_si*B+B0_si*dA)+tih*dA*dB;
    Is3 = A0_si*A2 + B0_si*B2 + A20_si*dA + B20_si*dB;

    // S1..S3 将 A²、B² 和 AB 组合成 Legendre 二次方向系数。
    double S1, S2, S3;
    S1 = -1.0/3+0.5*(AA+BB);

    // 平方差与交叉项分别对应两类角谐波权重。
    S2 = 0.5*(AA-BB);
    S3 = AB;

    // G 组把方向权重与偏近点角 E 的谐波相接，se/ce[j] 为 sin/cos(jE)。
    double G1, G2, G5, G6, G7, G8, G9, G10;
    double ece = e*ce[1];

    // 先准备 e*sin(E)，与上一行 e*cos(E) 配对。
    double ese = e*se[1];

    // G1/G2 分别进入半长轴、偏心率相关项，使用同一 S1..S3 几何。
    G1  = S1*(-e2-2.0*ece+0.5*e2*ce[2]) + S2*(-e2-2.0*ece+(1.0-e2/2.0)*ce[2]) + S3*xe2*(-2.0*ese+se[2]);
    G2  = S1*je2*(-e-2.0*ce[1]+0.5*e*ce[2]) + S2*je2*(-e+3.0*ce[1]-1.5*e*ce[2]+ce[3]/3.0) + S3*xe2*((3.0-e2/2)*se[1]-1.5*e*se[2]+(1.0-e2/2.0)*se[3]/3);

    // G8/G9 汇集一至三倍 E 的谐波，后面被倾角及节点组合共用。
    G8  = xe2*(5.0*ece-(1.0+e2)*ce[2]+e*ce[3]/3.0);
    G9  = (5.0-2.5*e2)*ese-(1.0+e2/2)*se[2]+(1.0-e2/2)*e*se[3]/3.0;

    // G10 补充偏心率谐波，G5 随后结合方向变化乘积形成节点相关量。
    G10 = (2.0-0.75*e2)*ese-0.75*e2*se[2]+e3*se[3]/12.0;
    G5  = -0.5*(AB1+A1B)*G8 + 0.5*(AA1-BB1)*G9 + (AA1+BB1)*G10;

    // G6/G7 用于近地点和平经度的组合输出，仍保留原 xe2 因子。
    G6  = S1*((-2.0+e2)*se[1]+e*se[2]/2.0) + S2*((-3.0+3.0*e2)*se[1]-e*se[2]/2.0+se[3]/3.0) + S3/xe2*((3.0-5.5*e2)*ce[1]+(1.0+e2)*e*ce[2]/2-(1.0-e2/2)*ce[3]/3);
    G7  = S1*G10+0.5*S2*G9-0.5*S3*G8;

    // 二次展开前两槽保存半长轴和偏心率相关改正，幅值均带 ba3。
    s2[0] = 2.0*a*ba3*G1;
    s2[1] = ba3*G2;

    // 倾角与缩放节点槽分别由 Is 组合及 G5 给出，节点槽相当于 sin(i)*delta_Omega。
    s2[2] = ba3*(Is1*G8+Is2*G9+Is3*G10)/xe2;
    s2[3] = ba3*(-G5)/xe2;

    // 后两槽已组合 e*(delta_omega+delta_Omega) 与 delta_lambda，不能按普通 Kepler 六槽读取。
    s2[4] = ba3*(e*tih*(-G5)/xe2+xe2*G6);
    s2[5] = ba3*(-tih*G5/xe2+e*xe2/(1.0+xe2)*G6+7.0*G7);

    // 进入 Legendre 三次部分，先准备三次方向多项式所需系数。
    double S4, S5, S6, S7;
    double AB54 = 1.25*(AA+BB)-1.0;

    // A²-3B²、3A²-B² 与 A/B 相乘形成三次谐波权重。
    double AB35 = (AA-BB*3.0)*5.0/12;
    double A3B5 = (AA*3.0-BB)*5.0/12;

    // AB15 保存倾角方向乘积之和，供下方 G13 等使用。
    double AB15 = 2.5*(AA1+BB1);

    // 先组合 A 方向的三次项，再组合 B 方向的三次项。
    S4  = AB54*A;
    S5  = AB35*A;

    // B 方向与上组 A 方向使用不同平方组合，不交换系数。
    S6  = AB54*B;
    S7  = A3B5*B;

    // G13/G14/G17/G18 将三次方向项与 A1/B1 的变化连接。
    double G13, G14, G17, G18;
    G13 = AB15*A + AB54*A1;

    // 继续准备 A 向和 B 向不同三次谐波的方向变化权重。
    G14 = (AA1-3.0*BB1)*5.0/6*A + AB35*A1;
    G17 = AB15*B + AB54*B1;

    // G18 补齐 B 向另一支三次谐波的权重。
    G18 = (3.0*AA1-BB1)*5.0/6*B + A3B5*B1;

    // A*A2、B*B2 供三次倾角改正的组合使用。
    double AA2, BB2;
    AA2 = A*A2;

    // BB2 与 AA2 同属当前几何，随后也写入供外层复用的 products。
    BB2 = B*B2;

    // Is6/Is7 给出三次部分的倾角方向组合，保持原先较长的代数表达式。
    double Is6, Is7;
    Is6 = 2.5*(A0_si*A2*dB+B0_si*(AA2+B2*(B+dB))+A20_si*dA*dB + B20_si*dB*dB) + (B20_si+ci*A0_si-tih*dA)*AB54;

    // 第二个组合含平方差与交叉项，仍按原括号及顺序求值。
    Is7 = 5.0/12*(A0_si*(6.0*dB*dB-3.0*(A+dA)*dA) + B0_si*(6.0*AA2-2.0*BB2+(B+3.0*dB)*dA) + A20_si*6.0*A*dB + B20_si*(3.0*AA-BB-2.0*B*dB) + 3.0*ci*(A0_si*(AA-3.0*BB+(A+dA)*dA)-3.0*B0_si*(B+dB)*dA)) + 1.25*tih*(3.0*dB*dB-dA*dA)*dA;

    // ra=r/a；三次径向因子用于三次展开的半长轴槽。
    double ra3 = ra*ra*ra;

    // 本移植沿用 NAO 将 Is4/Is5 置零的实现，不将其解释为已补全的通用公式。
    double Is4 = 0.0, Is5 = 0.0;

    // 三次展开先给半长轴和偏心率相关项，使用 f/M 的不同谐波表。
    s3[0] = 2.0*a*ba4*(S4*(ra3*cf[1] + 2.5*(e+0.75*e3)) + S5*(ra3*cf[3]+4.375*e3) + S6*ra3*sf[1] + S7*ra3*sf[3]);
    s3[1] = ba4*(S4*(-cm[2])/4.0 + S5*(2.25*cm[2]+0.375*cm[4]) + S6*(-sm[2])/4.0 + S7*(2.25*sm[2]+0.375*sm[4]));

    // 倾角与缩放节点项各自结合 Is 或 G 方向系数，使用当前 M 的一、三倍角。
    s3[2] = ba4*(-Is4*sm[1]-Is5*sm[3]/3.0+Is6*cm[1]+Is7*cm[3]/3.0);
    s3[3] = ba4*(G13*sm[1]+G14*sm[3]/3.0-G17*cm[1]-G18*cm[3]/3.0);

    // 末两槽保持三次近地点及平经度组合，最后与二次结果一起转 NSG2。
    s3[4] = ba4*(S4*(-sm[2])/4.0 + S5*(-2.25*sm[2]+0.375*sm[4]) + S6*cm[2]/4.0 + S7*(2.25*cm[2]-0.375*cm[4]));
    s3[5] = 3.0*ba4*(S4*(-3.0*sm[1]) + S5*(-sm[3]) + S6*3.0*cm[1] + S7*cm[3]);

    // Legendre 四次部分在本实现仅写半长轴槽，先准备四次方向系数。
    double S8, S9, S10, S11, S12;
    double AAAA, BBBB, AABB;

    // ra4 随后由 ra3 再乘 ra，不改变径向定义。
    double ra4;

    // A⁴、B⁴、A²B² 是四次方向多项式的基础量。
    AAAA = AA*AA;
    BBBB = BB*BB;
    AABB = AA*BB;

    // S8/S9 收集和、差两种平方组合，保留参考展开系数。
    S8   = 0.2-AA-BB+0.875*(AAAA+BBBB+2*AABB);
    S9   = 7.0/6*(AAAA-BBBB)-AA+BB;

    // 其余权重分别与四倍角余弦及二、四倍角正弦配合。
    S10  = 7.0/24*(AAAA+BBBB)-7.0/4*AABB;
    S11  = 7.0/3*AB*(AA+BB)-2.0*AB;

    // S12 补齐交叉方向项，ra4 给出四次径向因子。
    S12  = 7.0/6*AB*(AA-BB);
    ra4  = ra3*ra;

    // 四次贡献只进入 delta_a，式中扣除对应偏心率平均部分。
    s4[0] = 2.5*ba5*a*(S8*(ra4-1.0-5.0*e2-1.875*e4) + S9*(ra4*cf[2]-5.25*e2-2.625*e4) + S10*(ra4*cf[4]-7.875*e4) + S11*ra4*sf[2] + S12*ra4*sf[4]);

    // 只保存双快变量四次半长轴项实际使用的三个几何乘积。
    // 这些值属于本轮卫星估计及当前第三体方向，调用方在本轮内传给双快变量公式。
    products.AA = AA;
    products.BB = BB;
    products.AB = AB;

    // 最终前三槽为 delta_a、delta_xi、delta_eta；只有半长轴还加四次贡献。
    nsg2s2_lunisolar[0] = s2[0] + s3[0] + s4[0];
    nsg2s2_lunisolar[1] = cwo*(s2[1]+s3[1]) - swo*(s2[4]+s3[4]);
    nsg2s2_lunisolar[2] = swo*(s2[1]+s3[1]) + cwo*(s2[4]+s3[4]);

    // h/k 共同组合倾角与缩放节点项，偏心率方向则已在前三槽按 omega+Omega 投影。
    nsg2s2_lunisolar[3] = 0.5*cih*co*(s2[2]+s3[2]) - 0.5*so/cih*(s2[3]+s3[3]);
    nsg2s2_lunisolar[4] = 0.5*cih*so*(s2[2]+s3[2]) + 0.5*co/cih*(s2[3]+s3[3]);

    // 平经度直接加二、三次组合槽，不重复补节点或近地点项。
    nsg2s2_lunisolar[5] = s2[5] + s3[5];

    tracer.CheckOut(__func__);
}

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
    double nsg2s2_lunisolar[])
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 双快变量使用卫星 M 与第三体 u'；这里规定各谐角最大倍数而不是数组长度。
    const size_t maxMultiW  = 4;
    const size_t maxMultiM  = 4;

    // 第三体纬度幅角最高三倍，后面的数组都包含零倍槽。
    const size_t maxMultiUp = 3;

    // 二、三次展开幅值与单快变量接口各自使用对应公式，不能直接交换 ba4 系数。
    double ba3  = 1.5*lspar.beta2*apow[3];
    double ba4  = lspar.beta3*apow[4];

    // 先准备 e² 及 1-2.5e²，后面各倍频共享此截断组合。
    double e2   = 1.0-je2;
    double e25  = 1.0-2.5*e2;

    // sip/cip 是第三体 sin(i')、cos(i')，从已准备的参数恢复。
    double sip  = lspar.si;
    double cip  = lspar.cip1-1.0;

    // 第三体节点方向用于构造卫星与第三体节点差。
    double sop  = sin(lspar.Omega);
    double cop  = cos(lspar.Omega);

    // sdo/cdo 分别为 sin(Omega-Omega')、cos(Omega-Omega')。
    double sdo  = so*cop - co*sop;
    double cdo  = co*cop + so*sop;

    // A1/A2 是卫星与第三体方向的两支和差投影，与单快变量同名量含义不同。
    double A1   =  0.5*(sip*si+(1.0+cip*ci)*cdo);
    double A2   = -0.5*(sip*si-(1.0-cip*ci)*cdo);

    // A3/A4 保留节点差正弦与两轨道倾角和差的组合。
    double A3   = -0.5*(cip+ci)*sdo;
    double A4   =  0.5*(cip-ci)*sdo;

    // 后缀 p 表示这里对卫星倾角 i 的一阶方向导数。
    double A1p  =  0.5*(sip*ci - cip*si*cdo);
    double A2p  = -A1p;

    // 另一对倾角导数具有相同表达式，因此沿用 A4p=A3p。
    double A3p  =  0.5*si*sdo;
    double A4p  =  A3p;

    // 后缀 pp 在这里表示对节点差的一阶方向导数，并不是 p 的二阶导数。
    double A1pp = -0.5*(1.0+cip*ci)*sdo;
    double A2pp = -0.5*(1.0-cip*ci)*sdo;

    // 节点方向导数用 cdo，与上一对含 sdo 的导数共同进入角改正。
    double A3pp = -0.5*(cip+ci)*cdo;
    double A4pp =  0.5*(cip-ci)*cdo;

    // alpha[0]=n'/n，后续槽存各组合频率的无量纲分母。
    double alpha[21] = {0.0};
    alpha[0] = lspar.n*sqrt(apow[3]);

    // 先为 1..4 倍卫星 M 准备 1 减/加 2alpha/k 两组分母。
    for (size_t k = 1; k <= 4; k++)
    {
        double xa = 2.0*alpha[0]/(double)k;
        alpha[k]    = 1.0 - xa;

        // 加号分母位于 k+4 槽，后续与加第三体纬度相角配对。
        alpha[k+4]  = 1.0 + xa;
    }

    // 其余槽准备 u' 的一、三倍与卫星 M 倍频组合所需分母。
    alpha[11] = 1.0-3.0*alpha[0];
    alpha[12] = 1.0-alpha[0]/3.0;
    alpha[13] = 1.0-alpha[0]/4.0;

    // 保留 1-1.5alpha、1-0.75alpha 等各频率组合及下标。
    alpha[14] = 1.0-1.5*alpha[0];
    alpha[15] = 1.0-0.75*alpha[0];
    alpha[16] = 1.0+3.0*alpha[0];

    // 对应的加号频率组合单独保存，不复用减号分母。
    alpha[17] = 1.0+alpha[0]/3.0;
    alpha[18] = 1.0+alpha[0]/4.0;
    alpha[19] = 1.0+1.5*alpha[0];

    // 最后一槽补全 1+0.75alpha；本函数没有额外增加共振筛选。
    alpha[20] = 1.0+0.75*alpha[0];

    // 由 omega+Omega 和 Omega 的三角量恢复 sin/cos(omega)。
    double sw = swo*co - cwo*so;
    double cw = cwo*co + swo*so;

    // 先准备 2omega 的三角量，再准备第三体 2u'。
    double s2w = 2.0*sw*cw;
    double c2w = 2.0*cw*cw-1.0;

    // 两种快角分别进入后面 jM+2omega 加减 2u' 的相角表。
    double s2u = sin(lspar.u*2.0);
    double c2u = cos(lspar.u*2.0);

    // smpu/cmpu 从 +2u' 起步，数组下标继续增加 M。
    double smpu[5]  = {s2u};
    double cmpu[5]  = {c2u};

    // smmu/cmmu 从 -2u' 起步，保持另一快角方向。
    double smmu[5]  = {-s2u};
    double cmmu[5]  = {c2u};

    // smw/cmw 从 2omega 起步，与不含第三体角的谐项对应。
    double smw[5]   = {s2w};
    double cmw[5]   = {c2w};

    // 这组零槽为 2omega+2u'，随后递推增加 M。
    double smwpu[5] = {s2w*c2u+c2w*s2u};
    double cmwpu[5] = {c2w*c2u-s2w*s2u};

    // 这组零槽为 2omega-2u'，节点方向权重仍由 A 组给出。
    double smwmu[5] = {s2w*c2u-c2w*s2u};
    double cmwmu[5] = {c2w*c2u+s2w*s2u};

    for (size_t k = 1; k <= 4; k++)
    {
        // 先递推 kM+2u' 的正余弦，均读取上一槽。
        smpu[k]     =  smpu[k-1]*cm[1] +  cmpu[k-1]*sm[1];
        cmpu[k]     =  cmpu[k-1]*cm[1] -  smpu[k-1]*sm[1];

        // 再递推 kM-2u'，与加角表分开保存。
        smmu[k]     =  smmu[k-1]*cm[1] +  cmmu[k-1]*sm[1];
        cmmu[k]     =  cmmu[k-1]*cm[1] -  smmu[k-1]*sm[1];

        // 不含第三体快角的一组为 kM+2omega。
        smw[k]      =   smw[k-1]*cm[1] +   cmw[k-1]*sm[1];
        cmw[k]      =   cmw[k-1]*cm[1] -   smw[k-1]*sm[1];

        // 随后生成 kM+2omega+2u'。
        smwpu[k]    = smwpu[k-1]*cm[1] + cmwpu[k-1]*sm[1];
        cmwpu[k]    = cmwpu[k-1]*cm[1] - smwpu[k-1]*sm[1];

        // 最后生成 kM+2omega-2u'，保留原数组更新次序。
        smwmu[k]    = smwmu[k-1]*cm[1] + cmwmu[k-1]*sm[1];
        cmwmu[k]    = cmwmu[k-1]*cm[1] - smwmu[k-1]*sm[1];
    }

    // 用 atan2 恢复 omega、M，角度仅用于下方整数倍谐角，不改变输入根数。
    double w    = atan2(sw, cw);
    double M    = atan2(sm[1], cm[1]);

    // up 是第三体纬度幅角 u'，单位 rad。
    double up   = lspar.u;

    // 谐角倍率为 M:0..4、w:0..4、u':0..3；固定尺寸避免为每条轨道分配三维容器。
    std::array<std::array<std::array<double,4>,5>,5> sxp{},cxp{},sxm{},cxm{};

    // 三维槽依次对应 M、omega、u' 倍数；零倍也参与准备。
    for (size_t j = 0; j <= maxMultiM; j++)
    {
        for (size_t k = 0; k <= maxMultiW; k++)
        {
            // 先组合 jM+komega，再在内层加入或减去 q*u'。
            double jMkw = (double)j*M+(double)k*w;

            for (size_t q = 0; q <= maxMultiUp; q++)
            {
                // px/mx 是第三体快角正负两支，同一 j、k、q 保持一一对应。
                double px = jMkw + (double)q*up;
                double mx = jMkw - (double)q*up;

                // 先保存加 q*u' 的正余弦，再保存减 q*u' 的正余弦。
                sxp[j][k][q] = sin(px);
                cxp[j][k][q] = cos(px);

                // 减号表用于后面 alpha 减号频率组，符号关系保持原实现。
                sxm[j][k][q] = sin(mx);
                cxm[j][k][q] = cos(mx);
            }
        }
    }

    // 只有本地除数 si_div 参与小倾角缓解，不改输入 sin(i) 及其他几何量。
    double si_div = si;

    // abs(sin i)<1e-6 时原实现把除数置为正 1e-6 并警告，阈值和符号策略照旧。
    if (std::abs(si) < 1e-6)
    {
        si_div = 1e-6;
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Zero inclination mitigation is activated!", bastion::WARNING_AVERAGING_ZEROINC_MITIGATION_ACTIVATED);
    }

    // C/S 组由 A1..A4 的平方与交叉项组成，分别作为余弦、正弦谐项权重。
    double C0   = 0.5*(A1*A1+A2*A2+A3*A3+A4*A4)-1.0/3;
    double C1   = A1*A2 + A3*A4;

    // 继续构造平方差权重，区分两支第三体快角方向。
    double C2   = A1*A2 - A3*A4;
    double C3   = 0.5*(A1*A1-A3*A3);

    // C4 完成余弦权重，S1 开始保存交叉方向的正弦权重。
    double C4   = 0.5*(A2*A2-A4*A4);
    double S1   = A1*A4 - A2*A3;

    // 各 S 项的加减号区分 A1/A2 与 A3/A4 的相位组合。
    double S2   = A1*A4 + A2*A3;
    double S3   = A1*A3;

    // S4 补齐二次方向系数，随后准备倾角导数和节点差导数。
    double S4   = A2*A4;

    // C*p、S*p 由乘积求导组成，表示各方向权重对卫星倾角的变化。
    double C0p  = A1*A1p+A2*A2p+A3*A3p+A4*A4p;
    double C1p  = A1*A2p+A1p*A2+A3*A4p+A3p*A4;

    // 倾角变化继续进入余弦组，保持原乘积展开次序。
    double C2p  = A1*A2p+A1p*A2-A3*A4p-A3p*A4;
    double C3p  = A1*A1p-A3*A3p;

    // 最后一个余弦导数后，转到同一倾角导数的正弦组。
    double C4p  = A2*A2p-A4*A4p;
    double S1p  = A1*A4p+A1p*A4-A2*A3p-A2p*A3;

    // 正弦方向权重同样逐项使用 A*p 组合。
    double S2p  = A1*A4p+A1p*A4+A2*A3p+A2p*A3;
    double S3p  = A1*A3p+A1p*A3;

    // S4p 补齐本组倾角导数，没有对已准备的 A 项重新求值。
    double S4p  = A2*A4p+A2p*A4;

    // C*pp、S*pp 对应节点差的一阶变化，不能当作二阶倾角导数。
    double C0pp = A1*A1pp+A2*A2pp+A3*A3pp+A4*A4pp;
    double C1pp = A1pp*A2+A1*A2pp+A3pp*A4+A3*A4pp;

    // 余弦组的节点变化由 A*pp 与 A 本身交叉相乘。
    double C2pp = A1pp*A2+A1*A2pp-A3pp*A4-A3*A4pp;
    double C3pp = A1*A1pp-A3*A3pp;

    // 继续完成余弦方向，再准备正弦方向的节点差导数。
    double C4pp = A2*A2pp-A4*A4pp;
    double S1pp = A1pp*A4+A1*A4pp-A2pp*A3-A2*A3pp;

    // S 组保留各乘积的节点方向符号，与 p 组的导数方向不同。
    double S2pp = A1pp*A4+A1*A4pp+A2pp*A3+A2*A3pp;
    double S3pp = A1pp*A3+A1*A3pp;

    // S4pp 补齐二次方向的节点变化，供倾角/节点耦合槽使用。
    double S4pp = A2pp*A4+A2*A4pp;

    // 以下为 Legendre 二次的双快变量改正，后面仍须组合成最终 NSG2。
    double s2[6] = {0.0};

    // 半长轴槽按 C0..C4、S1..S4 权重叠加 jM+2omega 加减 2u' 的谐项。
    s2[0] = 2.0*ba3*apow[1]*(C0*(-2.0*e*cm[1]-e2/2.0*cm[2])
                            +C1*(-e*(cmmu[1]/alpha[1]+cmpu[1]/alpha[5])-e2/4.0*(cmmu[2]/alpha[2]+cmpu[2]/alpha[6]))

                            // 继续累计 C2 组余弦方向贡献，仍属于同一 s2[0]。
                            +C2*(-3.0*e*cmw[1]+e25*cmw[2]+e*cmw[3]+e2*cmw[4])
                            +C3*(-3.0*e*cmwmu[1]/alpha[1]+e25*cmwmu[2]/alpha[2]+e*cmwmu[3]/alpha[3]+e2*cmwmu[4]/alpha[4])

                            // 继续累计 C4 组余弦方向贡献，仍属于同一 s2[0]。
                            +C4*(-3.0*e*cmwpu[1]/alpha[5]+e25*cmwpu[2]/alpha[6]+e*cmwpu[3]/alpha[7]+e2*cmwpu[4]/alpha[8])
                            +S1*(e*(smmu[1]/alpha[1]-smpu[1]/alpha[5])+e2/4.0*(smmu[2]/alpha[2]-smpu[2]/alpha[6]))

                            // 继续累计 S2 组正弦方向贡献，仍属于同一 s2[0]。
                            +S2*(-3.0*e*smw[1]+e25*smw[2]+e*smw[3]+e2*smw[4])
                            +S3*(-3.0*e*smwmu[1]/alpha[1]+e25*smwmu[2]/alpha[2]+e*smwmu[3]/alpha[3]+e2*smwmu[4]/alpha[4])

                            // 继续累计 S4 组正弦方向贡献，仍属于同一 s2[0]。
                            +S4*(-3.0*e*smwpu[1]/alpha[5]+e25*smwpu[2]/alpha[6]+e*smwpu[3]/alpha[7]+e2*smwpu[4]/alpha[8]));

    // 第一个偏心率相关槽使用三维余弦/正弦表，仍按原谐角倍数和频率分母求和。
    s2[1] = ba3*(C0*(-2.0*cxm[1][1][0]-0.5*e*cxm[2][1][0])
                 +C1*(-(cxm[1][1][2]/alpha[1]+cxp[1][1][2]/alpha[5])-0.25*e*(cxm[2][1][2]/alpha[2]+cxp[2][1][2]/alpha[6]))

                 // 继续累计 C2 组余弦方向贡献，仍属于同一 s2[1]。
                 +C2*(3.0*cxm[1][1][0]+e*cxm[2][1][0]-1.5*e*cxm[2][3][0]+cxm[3][3][0]/3.0+0.5*e*cxm[4][3][0])
                 +C3*(3.0*cxm[1][1][2]/alpha[1]+e*cxm[2][1][2]/alpha[2]-1.5*e*cxm[2][3][2]/alpha[2]+cxm[3][3][2]/3.0/alpha[3]+0.5*e*cxm[4][3][2]/alpha[4])

                 // 继续累计 C4 组余弦方向贡献，仍属于同一 s2[1]。
                 +C4*(3.0*cxp[1][1][2]/alpha[5]+e*cxp[2][1][2]/alpha[6]-1.5*e*cxp[2][3][2]/alpha[6]+cxp[3][3][2]/3.0/alpha[7]+0.5*e*cxp[4][3][2]/alpha[8])
                 +S1*((sxm[1][1][2]/alpha[1]-sxp[1][1][2]/alpha[5])+0.25*e*(sxm[2][1][2]/alpha[2]-sxp[2][1][2]/alpha[6]))

                 // 继续累计 S2 组正弦方向贡献，仍属于同一 s2[1]。
                 +S2*(3.0*sxm[1][1][0]+e*sxm[2][1][0]-1.5*e*sxm[2][3][0]+sxm[3][3][0]/3.0+0.5*e*sxm[4][3][0])
                 +S3*(3.0*sxm[1][1][2]/alpha[1]+e*sxm[2][1][2]/alpha[2]-1.5*e*sxm[2][3][2]/alpha[2]+sxm[3][3][2]/3.0/alpha[3]+0.5*e*sxm[4][3][2]/alpha[4])

                 // 继续累计 S4 组正弦方向贡献，仍属于同一 s2[1]。
                 +S4*(3.0*sxp[1][1][2]/alpha[5]+e*sxp[2][1][2]/alpha[6]-1.5*e*sxp[2][3][2]/alpha[6]+sxp[3][3][2]/3.0/alpha[7]+0.5*e*sxp[4][3][2]/alpha[8]));

    // 第二个偏心率相关槽与上一槽配对，最终一起按节点方向投影到 xi/eta。
    s2[2] = ba3*(C0*(-2.0*sxm[1][1][0]-0.5*e*sxm[2][1][0])
                 +C1*(-(sxm[1][1][2]/alpha[1]+sxp[1][1][2]/alpha[5])-0.25*e*(sxm[2][1][2]/alpha[2]+sxp[2][1][2]/alpha[6]))

                 // 继续累计 C2 组余弦方向贡献，仍属于同一 s2[2]。
                 +C2*(-3.0*sxm[1][1][0]-e*sxm[2][1][0]-1.5*e*sxm[2][3][0]+sxm[3][3][0]/3.0+0.5*e*sxm[4][3][0])
                 +C3*(-3.0*sxm[1][1][2]/alpha[1]-e*sxm[2][1][2]/alpha[2]-1.5*e*sxm[2][3][2]/alpha[2]+sxm[3][3][2]/3.0/alpha[3]+0.5*e*sxm[4][3][2]/alpha[4])

                 // 继续累计 C4 组余弦方向贡献，仍属于同一 s2[2]。
                 +C4*(-3.0*sxp[1][1][2]/alpha[5]-e*sxp[2][1][2]/alpha[6]-1.5*e*sxp[2][3][2]/alpha[6]+sxp[3][3][2]/3.0/alpha[7]+0.5*e*sxp[4][3][2]/alpha[8])
                 +S1*(-(cxm[1][1][2]/alpha[1]-cxp[1][1][2]/alpha[5])-0.25*e*(cxm[2][1][2]/alpha[2]-cxp[2][1][2]/alpha[6]))

                 // 继续累计 S2 组正弦方向贡献，仍属于同一 s2[2]。
                 +S2*(3.0*cxm[1][1][0]+e*cxm[2][1][0]+1.5*e*cxm[2][3][0]-cxm[3][3][0]/3.0-0.5*e*cxm[4][3][0])
                 +S3*(3.0*cxm[1][1][2]/alpha[1]+e*cxm[2][1][2]/alpha[2]+1.5*e*cxm[2][3][2]/alpha[2]-cxm[3][3][2]/3.0/alpha[3]-0.5*e*cxm[4][3][2]/alpha[4])

                 // 继续累计 S4 组正弦方向贡献，仍属于同一 s2[2]。
                 +S4*(3.0*cxp[1][1][2]/alpha[5]+e*cxp[2][1][2]/alpha[6]+1.5*e*cxp[2][3][2]/alpha[6]-cxp[3][3][2]/3.0/alpha[7]-0.5*e*cxp[4][3][2]/alpha[8]));

    // 倾角槽结合节点差导数 pp 和 ci 权重，除数使用已检查的 si_div。
    s2[3] = ba3/si_div*(ci*(C2*(-6.0*e*cxm[1][2][0]+cxm[2][2][0]+2.0*e/3*cxm[3][2][0])
                        +C3*(-6.0*e*cxm[1][2][2]/alpha[1]+cxm[2][2][2]/alpha[2]+2.0*e/3*cxm[3][2][2]/alpha[3])

                        // 继续累计 C4 组余弦方向贡献，仍属于同一 s2[3]。
                        +C4*(-6.0*e*cxp[1][2][2]/alpha[5]+cxp[2][2][2]/alpha[6]+2.0*e/3*cxp[3][2][2]/alpha[7])
                        +S2*(-6.0*e*sxm[1][2][0]+sxm[2][2][0]+2.0*e/3*sxm[3][2][0])

                        // 继续累计 S3 组正弦方向贡献，仍属于同一 s2[3]。
                        +S3*(-6.0*e*sxm[1][2][2]/alpha[1]+sxm[2][2][2]/alpha[2]+2.0*e/3*sxm[3][2][2]/alpha[3])
                        +S4*(-6.0*e*sxp[1][2][2]/alpha[5]+sxp[2][2][2]/alpha[6]+2.0*e/3*sxp[3][2][2]/alpha[7]))

                    // 继续同一 s2[3] 的方向耦合组，保持各频率分母与原括号顺序。
                    -(C0pp*(-2.0*e*sm[1])
                      +C1pp*(-1.0*e*(sxm[1][0][2]/alpha[1]+sxp[1][0][2]/alpha[5]))

                      // 继续累计 C2pp 组余弦方向对节点差的一阶变化贡献，仍属于同一 s2[3]。
                      +C2pp*(-3.0*e*sxm[1][2][0]+0.5*sxm[2][2][0]+e*sxm[3][2][0]/3.0)
                      +C3pp*(-3.0*e*sxm[1][2][2]/alpha[1]+0.5*sxm[2][2][2]/alpha[2]+e*sxm[3][2][2]/3.0/alpha[3])

                      // 继续累计 C4pp 组余弦方向对节点差的一阶变化贡献，仍属于同一 s2[3]。
                      +C4pp*(-3.0*e*sxp[1][2][2]/alpha[5]+0.5*sxp[2][2][2]/alpha[6]+e*sxp[3][2][2]/3.0/alpha[7])
                      +S1pp*(-e*(cxm[1][0][2]/alpha[1]-cxp[1][0][2]/alpha[5]))

                      // 继续累计 S2pp 组正弦方向对节点差的一阶变化贡献，仍属于同一 s2[3]。
                      +S2pp*(3.0*e*cxm[1][2][0]-0.5*cxm[2][2][0]-e*cxm[3][2][0]/3.0)
                      +S3pp*(3.0*e*cxm[1][2][2]/alpha[1]-0.5*cxm[2][2][2]/alpha[2]-e*cxm[3][2][2]/3.0/alpha[3])

                      // 继续累计 S4pp 组正弦方向对节点差的一阶变化贡献，仍属于同一 s2[3]。
                      +S4pp*(3.0*e*cxp[1][2][2]/alpha[5]-0.5*cxp[2][2][2]/alpha[6]-e*cxp[3][2][2]/3.0/alpha[7])));

    // 缩放节点相关槽使用倾角导数 p；后面与 e 和半倾角共同组合 NSG2。
    s2[4] = ba3*(C0p*(-2.0*e*sm[1])
                 +C1p*(-e*(sxm[1][0][2]/alpha[1]+sxp[1][0][2]/alpha[5]))

                 // 继续累计 C2p 组余弦方向对倾角的一阶变化贡献，仍属于同一 s2[4]。
                 +C2p*(-3.0*e*sxm[1][2][0]+0.5*sxm[2][2][0]+e*sxm[3][2][0]/3.0)
                 +C3p*(-3.0*e*sxm[1][2][2]/alpha[1]+0.5*sxm[2][2][2]/alpha[2]+e*sxm[3][2][2]/3.0/alpha[3])

                 // 继续累计 C4p 组余弦方向对倾角的一阶变化贡献，仍属于同一 s2[4]。
                 +C4p*(-3.0*e*sxp[1][2][2]/alpha[5]+0.5*sxp[2][2][2]/alpha[6]+e*sxp[3][2][2]/3.0/alpha[7])
                 +S1p*(-e*(cxm[1][0][2]/alpha[1]-cxp[1][0][2]/alpha[5]))

                 // 继续累计 S2p 组正弦方向对倾角的一阶变化贡献，仍属于同一 s2[4]。
                 +S2p*(3.0*e*cxm[1][2][0]-0.5*cxm[2][2][0]-e*cxm[3][2][0]/3.0)
                 +S3p*(3.0*e*cxm[1][2][2]/alpha[1]-0.5*cxm[2][2][2]/alpha[2]-e*cxm[3][2][2]/3.0/alpha[3])

                 // 继续累计 S4p 组正弦方向对倾角的一阶变化贡献，仍属于同一 s2[4]。
                 +S4p*(3.0*e*cxp[1][2][2]/alpha[5]-0.5*cxp[2][2][2]/alpha[6]-e*cxp[3][2][2]/3.0/alpha[7]));

    // 平经度中额外出现频率倒数的组合，先准备减号频率组。
    double xa1 = 3.5+3.0/alpha[1];
    double xa2 = 2.0+1.5/alpha[2];
    double xa3 = 7.0/6.0+1.0/alpha[3];

    // 加号频率组对应同一类组合，保持各槽独立。
    double xa5 = 3.5+3.0/alpha[5];
    double xa6 = 2.0+1.5/alpha[6];
    double xa7 = 7.0/6.0+1.0/alpha[7];

    // 二次平经度相关槽先累计各方向谐项，最后还要补半倾角乘节点槽。
    s2[5]= -ba3*(C0*(-13.0*e*sm[1])
                +C1*(-xa1*e*smmu[1]/alpha[1]-xa5*e*smpu[1]/alpha[5])

                // 继续累计 C2 组余弦方向贡献，仍属于同一 s2[5]。
                +C2*(-19.5*e*smw[1]+3.5*smw[2]+13.0/6*e*smw[3])
                +C3*(-3.0*xa3*e*smwmu[1]/alpha[1]+xa2*smwmu[2]/alpha[2]+xa3*e*smwmu[3]/alpha[3])

                // 继续累计 C4 组余弦方向贡献，仍属于同一 s2[5]。
                +C4*(-3.0*xa5*e*smwpu[1]/alpha[5]+xa6*smwpu[2]/alpha[6]+xa7*e*smwpu[3]/alpha[7])
                +S1*(-xa1*e*cmmu[1]/alpha[1]+xa5*e*cmpu[1]/alpha[5])

                // 继续累计 S2 组正弦方向贡献，仍属于同一 s2[5]。
                +S2*(19.5*e*cmw[1]-3.5*cmw[2]-13.0/6*e*cmw[3])
                +S3*(3.0*xa1*e*cmwmu[1]/alpha[1]-xa2*cmwmu[2]/alpha[2]-xa3*e*cmwmu[3]/alpha[3])

                // 继续累计 S4 组正弦方向贡献，仍属于同一 s2[5]。
                +S4*(3.0*xa5*e*cmwpu[1]/alpha[5]-xa6*cmwpu[2]/alpha[6]-xa7*e*cmwpu[3]/alpha[7]));

    // 三次展开从二次 C/S 与 A1..A4 相乘构造，保留原谐波截断。
    double C5   = A1*(C0-4.0/15)+(C1*A2+C2*A2+C3*A1+S1*A4+S2*A4+S3*A3)/2.0;
    double C6   = (C1*A1+C3*A2+S1*A3+S3*A4)/2.0;

    // 这一对构造同一第三体方向的另一组奇次谐波余弦权重。
    double C7   = (C2*A1+C3*A2-S2*A3-S3*A4)/2.0;
    double C8   = (C3*A1-S3*A3)/2.0;

    // C9..C12 对应另一第三体快角方向，不能与 C5..C8 交换。
    double C9   = A2*(C0-4.0/15)+(C1*A1+C2*A1+C4*A2-S1*A3+S2*A3+S4*A4)/2.0;
    double C10  = (C1*A2+C4*A1-S1*A4+S4*A3)/2.0;

    // 继续补齐第二方向余弦权重。
    double C11  = (C2*A2+C4*A1-S2*A4-S4*A3)/2.0;
    double C12  = (C4*A2-S4*A4)/2.0;

    // 正弦权重同样由 A 组与二次 C/S 交叉形成。
    double S5   = A3*(C0-4.0/15)+(C1*A4-C2*A4-C3*A3-S1*A2+S2*A2+S3*A1)/2.0;
    double S6   = (C1*A3-C3*A4-S1*A1+S3*A2)/2.0;

    // 继续准备第一方向奇次谐波的正弦权重。
    double S7   = (C2*A3+C3*A4+S2*A1+S3*A2)/2.0;
    double S8   = (C3*A3+S3*A1)/2.0;

    // 第二方向正弦权重使用对应 A4 和交叉系数。
    double S9   = A4*(C0-4.0/15)+(C1*A3-C2*A3-C4*A4+S1*A1+S2*A1+S4*A2)/2.0;
    double S10  = (C1*A4-C4*A3+S1*A2+S4*A1)/2.0;

    // 补齐三次正弦权重，下面分别准备倾角和节点差导数。
    double S11  = (C2*A4+C4*A3+S2*A2+S4*A1)/2.0;
    double S12  = (C4*A4+S4*A2)/2.0;

    // 三次 C*p 通过原乘积展开连接二次倾角导数与 A*p。
    double C5p  = C0p*A1+(C0-4.0/15)*A1p+0.5*(C1p*A2+C2p*A2+C3p*A1+S1p*A4+S2p*A4+S3p*A3+C1*A2p+C2*A2p+C3*A1p+S1*A4p+S2*A4p+S3*A3p);
    double C6p  = 0.5*(C1p*A1+C3p*A2+S1p*A3+S3p*A4+C1*A1p+C3*A2p+S1*A3p+S3*A4p);

    // 继续第一方向的倾角导数，保持两类乘积都存在。
    double C7p  = 0.5*(C2p*A1+C3p*A2-S2p*A3-S3p*A4+C2*A1p+C3*A2p-S2*A3p-S3*A4p);
    double C8p  = 0.5*(C3p*A1-S3p*A3+C3*A1p-S3*A3p);

    // 第二方向余弦权重也同时包含二次权重导数和 A 的导数。
    double C9p  = C0p*A2+(C0-4.0/15)*A2p+0.5*(C1p*A1+C2p*A1+C4p*A2-S1p*A3+S2p*A3+S4p*A4+C1*A1p+C2*A1p+C4*A2p-S1*A3p+S2*A3p+S4*A4p);
    double C10p = 0.5*(C1p*A2+C4p*A1-S1p*A4+S4p*A3+C1*A2p+C4*A1p-S1*A4p+S4*A3p);

    // 补齐第二方向余弦权重的倾角导数。
    double C11p = 0.5*(C2p*A2+C4p*A1-S2p*A4-S4p*A3+C2*A2p+C4*A1p-S2*A4p-S4*A3p);
    double C12p = 0.5*(C4p*A2-S4p*A4+C4*A2p-S4*A4p);

    // 对应正弦权重的倾角变化同样保留原乘积求导展开。
    double S5p  = C0p*A3+(C0-4.0/15)*A3p+0.5*(C1p*A4-C2p*A4-C3p*A3-S1p*A2+S2p*A2+S3p*A1+C1*A4p-C2*A4p-C3*A3p-S1*A2p+S2*A2p+S3*A1p);
    double S6p  = 0.5*(C1p*A3-C3p*A4-S1p*A1+S3p*A2+C1*A3p-C3*A4p-S1*A1p+S3*A2p);

    // 继续第一方向正弦权重的倾角导数。
    double S7p  = 0.5*(C2p*A3+C3p*A4+S2p*A1+S3p*A2+C2*A3p+C3*A4p+S2*A1p+S3*A2p);
    double S8p  = 0.5*(C3p*A3+S3p*A1+C3*A3p+S3*A1p);

    // 第二方向倾角变化使用 S9 所对应的 A4 与二次系数。
    double S9p  = C0p*A4+(C0-4.0/15)*A4p+0.5*(C1p*A3-C2p*A3-C4p*A4+S1p*A1+S2p*A1+S4p*A2+C1*A3p-C2*A3p-C4*A4p+S1*A1p+S2*A1p+S4*A2p);
    double S10p = 0.5*(C1p*A4-C4p*A3+S1p*A2+S4p*A1+C1*A4p-C4*A3p+S1*A2p+S4*A1p);

    // 补齐三次正弦权重的倾角导数，供后面的缩放节点项使用。
    double S11p = 0.5*(C2p*A4+C4p*A3+S2p*A2+S4p*A1+C2*A4p+C4*A3p+S2*A2p+S4*A1p);
    double S12p = 0.5*(C4p*A4+S4p*A2+C4*A4p+S4*A2p);

    // 以下 pp 仍是节点差的一阶方向变化，保持源实现的全部组合项。
    double C5pp  = (C0-4.0/15)*A1pp+0.5*(C1pp*A2+C2pp*A2+C3pp*A1+S1pp*A4+S2pp*A4+S3pp*A3+C1*A2pp+C2*A2pp+C3*A1pp+S1*A4pp+S2*A4pp+S3*A3pp);
    double C6pp  = 0.5*(C1pp*A1+C3pp*A2+S1pp*A3+S3pp*A4+C1*A1pp+C3*A2pp+S1*A3pp+S3*A4pp);

    // 第一方向余弦导数继续使用 A*pp 与 C/S 的节点变化。
    double C7pp  = 0.5*(C2pp*A1+C3pp*A2-S2pp*A3-S3pp*A4+C2*A1pp+C3*A2pp-S2*A3pp-S3*A4pp);
    double C8pp  = 0.5*(C3pp*A1-S3pp*A3+C3*A1pp-S3*A3pp);

    // 第二方向的节点差导数与第一方向独立保存。
    double C9pp  = (C0-4.0/15)*A2pp+0.5*(C1pp*A1+C2pp*A1+C4pp*A2-S1pp*A3+S2pp*A3+S4pp*A4+C1*A1pp+C2*A1pp+C4*A2pp-S1*A3pp+S2*A3pp+S4*A4pp);
    double C10pp = 0.5*(C1pp*A2+C4pp*A1-S1pp*A4+S4pp*A3+C1*A2pp+C4*A1pp-S1*A4pp+S4*A3pp);

    // 补齐第二方向余弦权重的节点变化。
    double C11pp = 0.5*(C2pp*A2+C4pp*A1-S2pp*A4-S4pp*A3+C2*A2pp+C4*A1pp-S2*A4pp-S4*A3pp);
    double C12pp = 0.5*(C4pp*A2-S4pp*A4+C4*A2pp-S4*A4pp);

    // 正弦权重的节点变化同样保留原始乘积与符号。
    double S5pp  = (C0-4.0/15)*A3pp+0.5*(C1pp*A4-C2pp*A4-C3pp*A3-S1pp*A2+S2pp*A2+S3pp*A1+C1*A4pp-C2*A4pp-C3*A3pp-S1*A2pp+S2*A2pp+S3*A1pp);
    double S6pp  = 0.5*(C1pp*A3-C3pp*A4-S1pp*A1+S3pp*A2+C1*A3pp-C3*A4pp-S1*A1pp+S3*A2pp);

    // 继续第一方向正弦权重的节点差导数。
    double S7pp  = 0.5*(C2pp*A3+C3pp*A4+S2pp*A1+S3pp*A2+C2*A3pp+C3*A4pp+S2*A1pp+S3*A2pp);
    double S8pp  = 0.5*(C3pp*A3+S3pp*A1+C3*A3pp+S3*A1pp);

    // 第二方向正弦导数仍连接各自的 A*pp 和二次导数组合。
    double S9pp  = (C0-4.0/15)*A4pp+0.5*(C1pp*A3-C2pp*A3-C4pp*A4+S1pp*A1+S2pp*A1+S4pp*A2+C1*A3pp-C2*A3pp-C4*A4pp+S1*A1pp+S2*A1pp+S4*A2pp);
    double S10pp = 0.5*(C1pp*A4-C4pp*A3+S1pp*A2+S4pp*A1+C1*A4pp-C4*A3pp+S1*A2pp+S4*A1pp);

    // 补齐这组节点方向变化，后面不再重新准备几何系数。
    double S11pp = 0.5*(C2pp*A4+C4pp*A3+S2pp*A2+S4pp*A1+C2*A4pp+C4*A3pp+S2*A2pp+S4*A1pp);
    double S12pp = 0.5*(C4pp*A4+S4pp*A2+C4*A4pp+S4*A2pp);

    // 开始 Legendre 三次的双快变量项，使用 u' 的一倍、三倍谐角。
    double s3[6] = {0.0};

    // 半长轴槽叠加三次 C/S 权重，频率分母按各 M、omega、u' 倍数配对。
    s3[0] = 5.0*ba4*apow[1]*(C5*(cxm[1][1][1]/alpha[2]-0.5*e*cxm[2][1][1]/alpha[4])
                            +C6*(cxm[1][1][3]/alpha[11]-0.5*e*cxm[2][1][3]/alpha[14])

                            // 继续累计 C7 组余弦方向贡献，仍属于同一 s3[0]。
                            +C7*(-4.5*e*cxm[2][3][1]/alpha[4]+cxm[3][3][1]/alpha[12]+1.5*e*cxm[4][3][1]/alpha[13])
                            +C8*(-4.5*e*cxm[2][3][3]/alpha[14]+cxm[3][3][3]/alpha[2]+1.5*e*cxm[4][3][3]/alpha[15])

                            // 继续累计 C9 组余弦方向贡献，仍属于同一 s3[0]。
                            +C9*(cxp[1][1][1]/alpha[6]-0.5*e*cxp[2][1][1]/alpha[8])
                            +C10*(cxp[1][1][3]/alpha[16]-0.5*e*cxp[2][1][3]/alpha[19])

                            // 继续累计 C11 组余弦方向贡献，仍属于同一 s3[0]。
                            +C11*(-4.5*e*cxp[2][3][1]/alpha[8]+cxp[3][3][1]/alpha[17]+1.5*e*cxp[4][3][1]/alpha[18])
                            +C12*(-4.5*e*cxp[2][3][3]/alpha[19]+cxp[3][3][3]/alpha[6]+1.5*e*cxp[4][3][3]/alpha[20])

                            // 继续累计 S5 组正弦方向贡献，仍属于同一 s3[0]。
                            +S5*(sxm[1][1][1]/alpha[2]-0.5*e*sxm[2][1][1]/alpha[4])
                            +S6*(sxm[1][1][3]/alpha[11]-0.5*e*sxm[2][1][3]/alpha[14])

                            // 继续累计 S7 组正弦方向贡献，仍属于同一 s3[0]。
                            +S7*(-4.5*e*sxm[2][3][1]/alpha[4]+sxm[3][3][1]/alpha[12]+1.5*e*sxm[4][3][1]/alpha[13])
                            +S8*(-4.5*e*sxm[2][3][3]/alpha[14]+sxm[3][3][3]/alpha[2]+1.5*e*sxm[4][3][3]/alpha[15])

                            // 继续累计 S9 组正弦方向贡献，仍属于同一 s3[0]。
                            +S9*(sxp[1][1][1]/alpha[6]-0.5*e*sxp[2][1][1]/alpha[8])
                            +S10*(sxp[1][1][3]/alpha[16]-0.5*e*sxp[2][1][3]/alpha[19])

                            // 继续累计 S11 组正弦方向贡献，仍属于同一 s3[0]。
                            +S11*(-4.5*e*sxp[2][3][1]/alpha[8]+sxp[3][3][1]/alpha[17]+1.5*e*sxp[4][3][1]/alpha[18])
                            +S12*(-4.5*e*sxp[2][3][3]/alpha[19]+sxp[3][3][3]/alpha[6]+1.5*e*sxp[4][3][3]/alpha[20]));

    // 三次倾角槽使用节点差导数 pp 与 ci 权重，仍除以同一个 si_div。
    s3[3] = 5.0*ba4/3/si_div*(ci*(C5*cxm[1][1][1]/alpha[2]+C6*cxm[1][1][3]/alpha[11]+C7*cxm[3][3][1]/alpha[12]+C8*cxm[3][3][3]/alpha[2]
                              +C9*cxp[1][1][1]/alpha[6]+C10*cxp[1][1][3]/alpha[16]+C11*cxp[3][3][1]/alpha[17]+C12*cxp[3][3][3]/alpha[6]

                              // 继续累计 S5 组正弦方向贡献，仍属于同一 s3[3]。
                              +S5*sxm[1][1][1]/alpha[2]+S6*sxm[1][1][3]/alpha[11]+S7*sxm[3][3][1]/alpha[12]+S8*sxm[3][3][3]/alpha[2]
                              +S9*sxp[1][1][1]/alpha[6]+S10*sxp[1][1][3]/alpha[16]+S11*sxp[3][3][1]/alpha[17]+S12*sxp[3][3][3]/alpha[6])

                          // 继续同一 s3[3] 的方向耦合组，保持各频率分母与原括号顺序。
                          -(C5pp*sxm[1][1][1]/alpha[2]+C6pp*sxm[1][1][3]/alpha[11]+C7pp*sxm[3][3][1]/3.0/alpha[12]+C8pp*sxm[3][3][3]/3.0/alpha[2]
                            +C9pp*sxp[1][1][1]/alpha[6]+C10pp*sxp[1][1][3]/alpha[16]+C11pp*sxp[3][3][1]/3.0/alpha[17]+C12pp*sxp[3][3][3]/3.0/alpha[6]

                            // 继续累计 S5pp 组正弦方向对节点差的一阶变化贡献，仍属于同一 s3[3]。
                            -S5pp*cxm[1][1][1]/alpha[2]-S6pp*cxm[1][1][3]/alpha[11]-S7pp*cxm[3][3][1]/3.0/alpha[12]-S8pp*cxm[3][3][3]/3.0/alpha[2]
                            -S9pp*cxp[1][1][1]/alpha[6]-S10pp*cxp[1][1][3]/alpha[16]-S11pp*cxp[3][3][1]/3.0/alpha[17]-S12pp*cxp[3][3][3]/3.0/alpha[6]));

    // 三次缩放节点槽改用倾角导数 p，保持原源式各三角函数的选择。
    s3[4] = 2.5*ba4*(C5p*sxm[1][1][1]/alpha[2]+C6p*sxm[1][1][3]/alpha[11]+C7p*sxm[3][3][1]/3.0/alpha[12]+C8p*sxm[3][3][3]/3.0/alpha[2]
                     +C9p*sxp[1][1][1]/alpha[6]+C10p*sxp[1][1][3]/alpha[16]+C11p*sxp[3][3][1]/3.0/alpha[17]+C12p*sxp[3][3][3]/3.0/alpha[6]

                     // 继续累计 S5p 组正弦方向对倾角的一阶变化贡献，仍属于同一 s3[4]。
                     -S5p*cxm[1][1][1]/alpha[2]-S6p*cxm[1][1][3]/alpha[11]-S7p*cxm[3][3][1]/3.0/alpha[12]-S8p*cxm[3][3][3]/3.0/alpha[2]
                     -S9p*cxp[1][1][1]/alpha[6]-S10p*cxp[1][1][3]/alpha[16]-S11p*cxp[3][3][1]/3.0/alpha[17]-S12p*sxp[3][3][3]/3.0/alpha[6]);

    // 第一偏心率相关槽使用两倍、四倍 M 的三维谐角表。
    s3[1] = 2.5*ba4*(C5*(-0.25*cxm[2][2][1]/alpha[4])+C6*(-0.25*cxm[2][2][3]/alpha[14])
                     +C7*(2.25*cxm[2][2][1]/alpha[4]+0.375*cxm[4][4][1]/alpha[13])

                     // 继续累计 C8 组余弦方向贡献，仍属于同一 s3[1]。
                     +C8*(2.25*cxm[2][2][3]/alpha[14]+0.375*cxm[4][4][3]/alpha[15])
                     +C9*(-0.25*cxp[2][2][1]/alpha[8])+C10*(-0.25*cxp[2][2][3]/alpha[19])

                     // 继续累计 C11 组余弦方向贡献，仍属于同一 s3[1]。
                     +C11*(2.25*cxp[2][2][1]/alpha[8]+0.375*cxp[4][4][1]/alpha[18])
                     +C12*(2.25*cxp[2][2][3]/alpha[19]+0.375*cxp[4][4][3]/alpha[20])

                     // 继续累计 S5 组正弦方向贡献，仍属于同一 s3[1]。
                     +S5*(-0.25*sxm[2][2][1]/alpha[4])+S6*(-0.25*sxm[2][2][3]/alpha[14])
                     +S7*(2.25*sxm[2][2][1]/alpha[4]+0.375*sxm[4][4][1]/alpha[13])

                     // 继续累计 S8 组正弦方向贡献，仍属于同一 s3[1]。
                     +S8*(2.25*sxm[2][2][3]/alpha[14]+0.375*sxm[4][4][3]/alpha[15])
                     +S9*(-0.25*sxp[2][2][1]/alpha[8])+S10*(-0.25*sxp[2][2][3]/alpha[19])

                     // 继续累计 S11 组正弦方向贡献，仍属于同一 s3[1]。
                     +S11*(2.25*sxp[2][2][1]/alpha[8]+0.375*sxp[4][4][1]/alpha[18])
                     +S12*(2.25*sxp[2][2][3]/alpha[19]+0.375*sxp[4][4][3]/alpha[20]));

    // 第二偏心率相关槽与上一槽配对，正弦/余弦及系数均照原式。
    s3[2] = 2.5*ba4*(C5*(-0.25*sxm[2][2][1]/alpha[4])+C6*(-0.25*sxm[2][2][3]/alpha[14])
                     +C7*(-2.25*sxm[2][2][1]/alpha[4]+0.375*sxm[4][4][1]/alpha[13])

                     // 继续累计 C8 组余弦方向贡献，仍属于同一 s3[2]。
                     +C8*(-2.25*sxm[2][2][3]/alpha[14]+0.375*sxm[4][4][3]/alpha[15])
                     +C9*(-0.25*sxp[2][2][1]/alpha[8])+C10*(-0.25*sxp[2][2][3]/alpha[19])

                     // 继续累计 C11 组余弦方向贡献，仍属于同一 s3[2]。
                     +C11*(-2.25*sxp[2][2][1]/alpha[8]+0.375*sxp[4][4][1]/alpha[18])
                     +C12*(-2.25*sxp[2][2][3]/alpha[19]+0.375*sxp[4][4][3]/alpha[20])

                     // 继续累计 S5 组正弦方向贡献，仍属于同一 s3[2]。
                     +S5*(0.25*cxm[2][2][1]/alpha[4])+S6*(0.25*cxm[2][2][3]/alpha[14])
                     +S7*(2.25*cxm[2][2][1]/alpha[4]-0.375*cxm[4][4][1]/alpha[13])

                     // 继续累计 S8 组正弦方向贡献，仍属于同一 s3[2]。
                     +S8*(2.25*cxm[2][2][3]/alpha[14]-0.375*cxm[4][4][3]/alpha[15])
                     +S9*(0.25*cxp[2][2][1]/alpha[8])+S10*(0.25*cxp[2][2][3]/alpha[19])

                     // 继续累计 S11 组正弦方向贡献，仍属于同一 s3[2]。
                     +S11*(2.25*cxp[2][2][1]/alpha[8]-0.375*cxp[4][4][1]/alpha[18])
                     +S12*(2.25*cxp[2][2][3]/alpha[19]-0.375*cxp[4][4][3]/alpha[20]));

    // 三次平经度相关组合保留原来的 22.5 幅值及全部谐角选择。
    s3[5] = -22.5*ba4*(C5*sxm[1][1][1]/alpha[2]+C6*sxm[1][1][3]/alpha[11]+C7*sxm[3][3][1]/3.0/alpha[12]+C8*sxm[3][3][3]/3.0/alpha[2]
                       +C9*sxp[1][1][1]/alpha[6]+C10*sxp[1][1][3]/alpha[16]+C11*sxp[3][3][1]/3.0/alpha[17]+C12*sxp[3][3][3]/3.0/alpha[6]

                       // 继续累计 S5 组正弦方向贡献，仍属于同一 s3[5]。
                       -S5*cxm[1][1][1]/alpha[2]-S6*cxm[1][1][3]/alpha[11]-S7*cxm[3][3][1]/3.0/alpha[12]-S8*cxm[3][3][3]/3.0/alpha[2]
                       -S9*cxp[1][1][1]/alpha[6]-S10*cxp[1][1][3]/alpha[16]-S11*sxp[3][3][1]/3.0/alpha[17]-S12*sxp[3][3][3]/3.0/alpha[6]);

    // products 提供同一配置的 A²、B²、AB；四次半长轴项使用它们补齐几何多项式。
    double AA   = products.AA;
    double BB   = products.BB;

    // 用当前乘积再构造 A⁴、B⁴，供四次方向多项式使用。
    double AB   = products.AB;
    double AAAA = AA*AA;

    // B⁴ 与 A⁴ 采用相同当前几何，随后只进入四次半长轴贡献。
    double BBBB = BB*BB;

    // 四次部分只写半长轴槽，a6 为 a⁶。
    double s4[6] = {0.0};
    double a6 = apow[3]*apow[3];

    // 四次半长轴按二倍、四倍 M 的谐波组合，不额外改变其他槽。
    s4[0] = lspar.beta4*a6*(5.0/8*((7.0*(AAAA-BBBB)-6.0*(AA-BB))*cm[2]+2.0*AB*(7.0*(AA+BB)-6.0)*sm[2])
                           +35.0/32*((AAAA+BBBB-6.0*AB*AB)*cm[4]+4.0*AB*(AA-BB)*sm[4]));

    // 输出前三槽为 delta_a、delta_xi、delta_eta；偏心率组合还要加入缩放节点影响。
    nsg2s2_lunisolar[0] = s2[0]+s3[0]+s4[0];
    nsg2s2_lunisolar[1] = (s2[1]+s3[1])*co-(s2[2]+s3[2])*so - e*swo*(s2[4]+s3[4])/si_div;
    nsg2s2_lunisolar[2] = (s2[1]+s3[1])*so+(s2[2]+s3[2])*co + e*cwo*(s2[4]+s3[4])/si_div;

    // h/k 组合倾角槽与缩放节点槽，半倾角余弦保留原分母。
    nsg2s2_lunisolar[3] = 0.5*cih*co*(s2[3]+s3[3]) - 0.5*so/cih*(s2[4]+s3[4]);
    nsg2s2_lunisolar[4] = 0.5*cih*so*(s2[3]+s3[3]) + 0.5*co/cih*(s2[4]+s3[4]);

    // 末槽把半倾角乘节点项与平经度相关槽相加，得到 delta_lambda。
    nsg2s2_lunisolar[5] = sih/cih*(s2[4]+s3[4])+(s2[5]+s3[5]);

    tracer.CheckOut(__func__);
}

} // namespace orbit::averaging_detail
