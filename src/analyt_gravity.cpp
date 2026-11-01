/**
 * @file analyt_gravity.cpp
 * @brief J2、J22 及岁差章动与 J2 耦合的解析短周期公式。
 * @details 输入当前均根数估计的几何；NSG2 六槽覆写，J2 二阶半长轴和平经度另行输出，J22 高阶半长轴用于替换对应槽。
 * 长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad；输出为短周期改正，不是状态时间导数。
 */

// 解析短周期公式的定向移植与本库适配，原公式作者 Jingshi Tang。
// Copyright (c) 2020 国家航天局空间碎片监测与应用中心。All rights reserved.
// 仅保留 O2M 实际使用的代数项；轨道、时间、重力及迭代流程由 CDC_lib 提供。
#include "analyt_terms.hpp"
#include <cmath>

namespace orbit::averaging_detail {
namespace {
constexpr int KeplerFull = 0, KeplerCombined = 2;
constexpr int ResonanceHalf = 1, ResonanceSynchronous = 2, ResonanceThreeHalves = 3;
}

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
    double nsg2s1_J2[], double& as2_J2)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 准备 e²、e³ 和 1-e²，后面的径向项与近圆组合共用这些无量纲因子。
    double e2   = e*e;
    double e3   = e2*e;
    double je2  = 1.0-e2;

    // sin²(i)、sin⁴(i) 与 sin(2i) 控制 J2 对轨道平面方向的依赖。
    double si2  = si*si;
    double si4  = si2*si2;
    double s2i  = 2.0*si*ci;

    // ar=a/r；径向三次、四次幂均对应本轮当前卫星位置。
    double ar3  = ar*ar*ar;
    double ar4  = ar3*ar;

    // xe2=sqrt(1-e²)，xe3=(1-e²)^(-3/2)，与径向周期部分共同进入展开。
    double xe2  = sqrt(je2);
    double xe3  = 1.0/xe2/je2;

    // 这两个多项式分别收集 1-1.5sin²(i) 和 2-2.5sin²(i) 的重复权重。
    double si23 = 1.0-1.5*si2;
    double si25 = 2.0-2.5*si2;

    // w 是近地点幅角；保存 2w 的三角量，供含真近点角的谐波组合。
    double wx2  = w*2.0;
    double c2w  = cos(wx2);
    double s2w  = sin(wx2);

    // 先准备 cos(f) 及 cos(kf+2w)，其中 k 是真近点角倍数。
    double cf    = cos(f);
    double c2f2w = cos(2.0*f+wx2);
    double c1f2w = cos(1.0*f+wx2);

    // 三倍 f 的加角余弦后，转到不含 2w 的正弦谐波。
    double c3f2w = cos(3.0*f+wx2);
    double sf    = sin(f);
    double s2f   = sin(2.0*f);

    // sin(3f) 与含 2w 的正弦项在角变量公式中采用不同权重。
    double s3f   = sin(3.0*f);
    double s2f2w = sin(2.0*f+wx2);
    double s1f2w = sin(1.0*f+wx2);

    // 继续准备三至五倍 f 加 2w 的正弦，保留原展开的谐波次序。
    double s3f2w = sin(3.0*f+wx2);
    double s4f2w = sin(4.0*f+wx2);
    double s5f2w = sin(5.0*f+wx2);

    // 补充 sin(f-2w)，与上面的加角谐波共同构成近地点项。
    double sfm2w = sin(f-wx2);

    // as1 是一阶半长轴短周期项，ar³ 与 xe3 的差包含相应平均部分的扣除。
    double as1 = ja*(si23*2.0/3*(ar3-xe3)+si2*ar3*c2f2w);

    // 先组合 e*cos(f) 与重复径向多项式，后面仍按原乘加顺序计算。
    double ecf = e*cf;
    double xecf = cf*(3.0+3.0*ecf+ecf*ecf);

    // 该偏心率组合与 xecf 一起构成一阶偏心率项 es1。
    double e222 = (2.0+xe2-e2)/(1.0+xe2);
    double es1= jp*(si23/3.0*(e*e222+xecf)+si2/2.0*((e+xecf)*c2f2w-je2*(c1f2w+c3f2w/3.0)));

    // is1、is2 保留本展开的一阶、二阶倾角表达式。
    double is1 = jp*s2i*(e*c1f2w+c2f2w+e*c3f2w/3.0)/4.0;
    double is2 = jp*jp*s2i*si2/(-8.0)*c2f2w;

    // 将一至三倍 f 的加角正弦收集起来，再构成节点改正。
    double s123f = e*s1f2w+s2f2w+e*s3f2w/3.0;
    double os1 = -jp*ci*(fmesf-s123f/2.0);

    // os2 保留原二阶节点式；后面 NSG2 一阶输出实际使用 os1。
    double os2 = jp*jp*ci*si2/2.0*s2f2w;

    // ews1 已乘偏心率 e，最终组合偏心率向量时无需再除以 e。
    double ews1 = jp*(si25*fmesf*e+si23*((1.0-e2/4.0)*sf+e/2.0*s2f+e2*s3f/12.0)-e/2.0*s123f+si2/48*(3.0*e2*sfm2w-(12.0-45.0*e2)*s1f2w+60.0*e*s2f2w+(28.0+19.0*e2)*s3f2w+18.0*e*s4f2w+3.0*e2*s5f2w));

    // wms1 合并近地点幅角和平近点角的改正，后面与节点项合成平经度。
    double wms1 = -ci*os1 + jp*(si23*fmesf +f1e*si23*((1.0-e2/4.0)*e*sf+e2/2.0*s2f+e3/12.0*s3f) +si2*(0.75*e*s1f2w+0.75*s2f2w+0.25*e*s3f2w)

                                // 此续行仍属于同一 wms1 的高次谐波部分，不是独立输出。
                                +f1e*si2*(-(0.25+e2*5.0/16.0)*e*s1f2w+(7.0/12.0-e2/48.0)*e*s3f2w+0.375*e2*s4f2w+e3/16.0*(s5f2w+sfm2w)));

    // 带 b 的量是从对应一阶表达式中扣除的平均部分。
    double es1b  = jp*si2*f2e*je2*e*c2w/6.0;
    double is1b  = -jp*s2i*f2e*e2*c2w/12.0;

    // 节点、近地点及角变量组合的平均项保留各自的倾角系数。
    double os1b  = -jp*ci*f2e*e2*s2w/6.0;
    double ws1b  = jp*(si2*(0.75+je2*f2e)+ci*ci*e2*f2e)*s2w/6.0;

    // wms1b 对应上面的 wms1，最终在平经度槽中相减。
    double wms1b = jp*(si2*f1e*(0.75+f2e+e2*f2e/2.0)+(si23-si2)*f2e)*e2*s2w/6.0;

    // 用已乘 e 的近地点组合构造平近点角相关量，供二阶半长轴交叉项使用。
    double ems1  = e*wms1 - ews1;
    double ems1b = e*(wms1b - ws1b);

    // e⁴ 与 adcdl 收集二阶半长轴中的重复部分，不另设新的截断。
    double e4 = e2*e2;
    double adcdl = jp*jp*a*xe2*(si23*si23*((16.0+19.0*e2)/9+35.0*e4/18/je2+2.0*xe2/9)+si2*(1.0+2.0*e2/3)

                                // 这里继续同一 adcdl 的 sin⁴(i) 与 e⁴ 项，括号及求值次序使用对应公式。
                                +si4*(-5.0/6+25.0*e2/24+35.0*e4/16/je2)+si2*e2*c2w*(-(4.0-si2)/3*f2e+(5.0/6-7.0*si2/4)+e2/je2*7.0/3*si23)+si4*e4*(c2w*c2w-s2w*s2w)/32.0/je2);

    // 二阶半长轴由一阶改正交叉项和 adcdl 组成，调用方会单独扣除它。
    as2_J2 = -2.0/a*(as1+ja*xe3*si23)*as1
            + ja*(si23*(2.0*ar4*cf-e*xe3/je2)-4.0/je2*si2*ar3*(1.0+ecf/2.0)*sf*s2f2w+3.0*si2*ar4*cf*c2f2w)*(es1-es1b)

            // 依次加入倾角及角变量的交叉改正，仍属于同一个二阶半长轴输出。
            + ja*(s2i*ar3*(c2f2w-1.0))*(is1-is1b)
            + ja*(-2.0*si2*ar3*s2f2w)*(wms1-wms1b)

            // 最后的交叉项采用原 xe2、je2 分母，再扣去 adcdl。
            + ja*(2.0*si2*ar3*s2f2w/je2/je2*(-e*f2e/f1e+e3-xe2*(2.0*cf+ecf*cf))-ar4*sf*(2.0*si23+3.0*si2*c2f2w)/xe2)*(ems1-ems1b)
            - adcdl;

    // wo=omega+Omega 是偏心率向量的方向，先求其正余弦。
    double cwo = cos(wo);
    double swo = sin(wo);

    // 由输入 h/k 恢复 sin(i/2)、cos(i/2)，用于倾角分量组合。
    double sih = sqrt(h*h+k*k);
    double cih = sqrt(1.0-sih*sih);

    // 先从一阶偏心率、倾角改正中扣除各自平均部分。
    double es1t     = es1 - es1b;
    double is1t     = is1 - is1b;

    // 节点短周期并入 e*(delta_omega+delta_Omega)，随后与偏心率改正一起投影。
    double os1t     = os1 - os1b;
    double ewos1t   = ews1-e*ws1b+e*os1t;

    // 提前收集倾角因子，避免最终 h/k 公式直接除以 sin(i/2)。
    double is1t_sih = (jp*(e*c1f2w+c2f2w+e*c3f2w/3.0)/4.0 - (-jp*f2e*e2*c2w/12.0))*(4.0*cih*ci);

    // 前三槽为 delta_a、delta_xi、delta_eta，偏心率分量沿 wo 投影。
    nsg2s1_J2[0] = as1;
    nsg2s1_J2[1] = cwo*es1t - swo*ewos1t;
    nsg2s1_J2[2] = swo*es1t + cwo*ewos1t;

    // h/k 两槽组合倾角与节点改正，方向由当前 NSG2 的 h/k 决定。
    nsg2s1_J2[3] = cih/2.0*h*is1t_sih - k*os1t;
    nsg2s1_J2[4] = cih/2.0*k*is1t_sih + h*os1t;

    // 末槽是 delta_lambda=delta_M+delta_omega+delta_Omega，角度单位 rad。
    nsg2s1_J2[5] = os1t+wms1-wms1b;

    tracer.CheckOut(__func__);
}

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
void J2SecondOrderLongitude(bastion::Overwatch& tracer, double j2p4, double ci, double si2, double si4, double s2f2w, double s4f4w, double& ls2_J2)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // D 系数用于二阶节点项，零项与非零项均沿用参考展开。
    const double D42 = 0.0;
    const double D52 = 9.0/8.0;

    // 继续准备二倍与四倍纬度幅角谐波所需的 D 系数。
    const double D62 = 0.0;
    const double D74 = -9.0/32.0;

    // D84 完成节点项系数，G52 起准备近地点加平近点角的组合系数。
    const double D84 = -3.0/32.0;
    const double G52 = 0.0;

    // 这些 G 系数分别与 sin²(i)、sin⁴(i) 配合进入二倍角项。
    const double G62 = -15.0/4.0;
    const double G72 = 99.0/32.0;

    // 四倍角项保留原系数和后续固定分母，不重新归一化。
    const double G84 = 9.0/2.0;
    const double G94 = -3.0;

    // 最后一个 G 系数补全 sin⁴(i) 对四倍角谐波的权重。
    const double G104 = -78.0;

    // 分别计算二阶节点改正及近地点加平近点角的组合改正。
    double os2  = j2p4*ci*((D42+si2*D52+(1.0-1.5*si2)*D62)*s2f2w + (D74+si2*D84)*s4f4w);
    double wms2 = j2p4*((G52+si2*G62+si4*G72)*s2f2w + (G84/16.0+si2*G94/64.0+si4*G104/108.0)*s4f4w);

    // 两者相加得到独立二阶平经度项，外层只扣除一次。
    ls2_J2 = os2+wms2;

    tracer.CheckOut(__func__);
}

/**
 * @brief 计算 J22 的 Kepler 或 NSG2 组合中间短周期项。
 * @param[in,out] tracer 记录公式计算中的错误、适用范围警告及调用链。
 * @param[in] j22a2 未归一化 J22 幅值除以 a^2，a 以 ae 归一化。
 * @param[in] a 当前均根数估计的半长轴，以 ae 归一化。
 * @param[in] e 当前均根数估计的偏心率，无量纲，0<=e<1。
 * @param[in] si 当前卫星轨道倾角的 sin(i)，无量纲。
 * @param[in] ci 当前卫星轨道倾角的 cos(i)，无量纲。
 * @param[in] o22 Omega 减格林尼治角及 J22 系数半相位之和，rad。
 * @param[in] w 当前卫星近地点幅角 omega，rad。
 * @param[in] m 当前卫星平近点角 M，rad。
 * @param[in] alpha 地球自转角速度与卫星平运动之比 omega_E/n，无量纲。
 * @param[out] keps2_J22 覆写六槽；前 3 槽为 delta_a、delta_e、delta_i，后 3 槽含义由 short_flag 决定，半长轴尺度 ae，角度 rad。
 * @param[in] short_flag KeplerFull=0 输出 Kepler 改正，KeplerCombined=2 输出供 NSG2 变换使用的组合槽。
 * @param[in] resonance_id 由密切输入固定的共振类别；1/2/3 为 alpha 约等于 1/2、1、3/2，其余类别按调用方的频率表解释，-1 为未命中。
 * @details 来自 analyt_j22.cpp 的 AnaJ22ShortKep。Full 输出普通 Kepler 改正；Combined 的第 3/4/5 槽分别为 delta_Omega、e*delta_omega、delta_lambda。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
static void J22KeplerShort(bastion::Overwatch& tracer, double j22a2, double a, double e, double si, double ci, double o22, double w, double m, double alpha, double keps2_J22[], int short_flag,
    int resonance_id)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 只接受普通 Kepler 或供 NSG2 组合的输出标志，其余按原逻辑报错。
    if ((short_flag != KeplerFull) && (short_flag != KeplerCombined))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown flag for J22 short-period terms!", bastion::ERROR);
        return;
    }

    // 先准备 1+cos(i) 及平方，供正节点相位一组权重使用。
    double cip1     = 1.0 + ci;
    double cip12    = cip1*cip1;

    // 1-cos(i) 及平方供另一组节点相位使用，不能与上一组交换。
    double cim1     = 1.0 - ci;
    double cim12    = cim1*cim1;

    // sin²(i) 给出幅值权重，2w 和 2o22 给出 J22 二阶相角。
    double si2      = si*si;
    double ww       = w*2.0;

    // o22 是 Omega 减 J22 参考经度，阶数为 2，故此处取两倍。
    double oo       = o22*2.0;

    // s2 先累计组合中间槽，最后按 short_flag 确定六槽输出含义。
    double s2[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    // A/B 相角对应 2w 加减 2o22，数组下标继续增加 M 倍数。
    double mwoA[5], mwoB[5];
    double smwoA[5], cmwoA[5], smwoB[5], cmwoB[5];

    // 这对相角不含 w，只保存 jM 加减 2o22 及其正余弦。
    double moA[3], moB[3];
    double smoA[3], cmoA[3], smoB[3], cmoB[3];

    // 零倍 M 的初值固定近地点及节点相位，循环每轮再增加一个 M。
    mwoA[0] = ww + oo;
    mwoB[0] = ww - oo;

    // 不含 w 的两组从正负 2o22 开始，与 mwoA/B 的方向对应。
    moA[0]  = oo;
    moB[0]  = -oo;

    // 依次生成 1..4 倍 M，相角递推和三角量赋值顺序保持原样。
    for (int j = 1; j < 5; j++)
    {
        // 两组相角分别在上一槽基础上加 M，节点相位正负号不变。
        mwoA[j] = m + mwoA[j-1];
        mwoB[j] = m + mwoB[j-1];

        // 先保存 A 相角的正余弦，再处理 B 相角。
        smwoA[j] = sin(mwoA[j]);
        cmwoA[j] = cos(mwoA[j]);

        // B 组的节点符号相反，后面与 1-cos(i) 一组系数配合。
        smwoB[j] = sin(mwoB[j]);
        cmwoB[j] = cos(mwoB[j]);

        // 不含 w 的数组仅有 0..2 槽，因此只在 j<3 时更新。
        if (j < 3)
        {
            // 不含 w 的两组同样逐次增加 M。
            moA[j] = m + moA[j-1];
            moB[j] = m + moB[j-1];

            // 保存 jM+2o22 的正余弦，供对应频率组使用。
            smoA[j] = sin(moA[j]);
            cmoA[j] = cos(moA[j]);

            // 保存 jM-2o22 的正余弦，供另一方向的节点谐项使用。
            smoB[j] = sin(moB[j]);
            cmoB[j] = cos(moB[j]);
        }
    }

    // eh=e/2；下一组 af21=1-2alpha，alpha=omega_E/n。
    double eh = e/2.0;

    double af21     = 1.0-2.0*alpha;

    // 将倾角权重除以同一频率分母，供不同根数槽共用。
    double i1af21   = cip1/af21;
    double i2af21   = cip12/af21;

    // 继续准备 sin²(i) 权重及角变量的额外频率因子。
    double si2af21  = 3.0*si2/af21;
    double yjaf21   = 1.0-0.5/af21;

    // 半同步共振为 alpha=1/2，该类轨道保留这组慢项，不作为短周期扣除。
    if (resonance_id != ResonanceHalf)
    {
        // 本频率先累计半长轴、偏心率与倾角中间项。
        s2[0] += -i2af21*eh*cmwoA[1]  +si2af21*e*cmoA[1];
        s2[1] +=  i2af21/2.0*cmwoA[1] +si2af21*cmoA[1];
        s2[2] +=  i1af21*eh*cmwoA[1]  -3.0*e/af21*cmoA[1];

        // 节点与已乘 e 的近地点项使用正弦相角，符号与前三槽不同。
        s2[3] +=  i1af21*eh*smwoA[1]  +ci*3.0*e/af21*smoA[1];
        s2[4] += -i2af21/2.0*smwoA[1] +si2af21*smoA[1];

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
        if (short_flag == KeplerFull)
        {
            s2[5] += -i2af21*e*yjaf21*smwoA[1]  +2.0*e*si2af21*yjaf21*smoA[1];
        }
        else if (short_flag == KeplerCombined)
        {
            // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

            s2[5] += 0.75*j22a2*eh*(-i2af21/2.0*smwoA[1] +si2af21*smoA[1]) + 2.25*j22a2*(-i2af21*e*yjaf21*smwoA[1] +2.0*e*si2af21*yjaf21*smoA[1]);
        }
    }

    // 1-alpha 对应同步共振附近的分母。
    double af11     = 1.0 - alpha;
    double i1af11   = cip1/af11;

    // 为同一频率准备倾角平方权重和 sin²(i) 权重。
    double i2af11   = cip12/af11;
    double si2af11  = 3.0*si2/af11;

    // 该频率因子进入第六槽的角变量组合。
    double yjaf11   = 1.0 - 0.5/af11;

    // 同步共振轨道跳过这一组慢项，其余轨道照原顺序累计。
    if (resonance_id != ResonanceSynchronous)
    {
        // 本组先累计半长轴、偏心率和倾角中间改正。
        s2[0] +=  i2af11*cmwoA[2];
        s2[1] += -i2af11*eh*cmwoA[2] +1.5*e*si2af11*cmoA[2];
        s2[2] += -i1af11/2.0*cmwoA[2];

        // 节点与近地点组合先存入 s2，最后统一恢复输出尺度。
        s2[3] += -i1af11/2.0*smwoA[2];
        s2[4] += -2.5*e*i2af11*smwoA[2] +1.5*e*si2af11*smoA[2];

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
        if (short_flag == KeplerFull)
        {
            s2[5] += i2af11*yjaf11*smwoA[2];
        }
        else if (short_flag == KeplerCombined)
        {
            // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

            s2[5] += 2.25*j22a2*i2af11*yjaf11*smwoA[2];
        }
    }

    // 1-2alpha/3 对应 alpha=3/2 共振，使用三倍 M 相角。
    double af23     = 1.0 - 2.0*alpha/3;
    double i1af23   = cip1*7.0/af23*e/6.0;

    // 该组倾角平方系数及角变量频率因子继续独立保存。
    double i2af23   = cip12*7.0/af23/6.0;
    double yjaf23   = 1.0 - 0.5/af23;

    // alpha=3/2 类轨道保留这一组共振慢项。
    if (resonance_id != ResonanceThreeHalves)
    {
        // 三倍 M 组依次累计半长轴、偏心率与倾角中间槽。
        s2[0] +=  i2af23*cmwoA[3]*3.0*e;
        s2[1] +=  i2af23*cmwoA[3];
        s2[2] += -i1af23*cmwoA[3];

        // 节点及近地点组合采用同一三倍 M 的正弦相角。
        s2[3] += -i1af23*smwoA[3];
        s2[4] +=  i2af23*smwoA[3];

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
        if (short_flag == KeplerFull)
        {
            s2[5] += i2af23*e*2.0*yjaf23*smwoA[3];
        }
        else if (short_flag == KeplerCombined)
        {
            // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

            s2[5] += 0.75*j22a2*eh*i2af23*smwoA[3] + 2.25*j22a2*i2af23*e*2.0*yjaf23*smwoA[3];
        }
    }

    // 纯节点项采用 2o22 的三角量，sio 再将 sin²(i) 除以 alpha。
    double s2o = sin(oo);
    double c2o = cos(oo);
    double sio = si2/alpha;

    // 不含 M 的节点频率项进入倾角、节点及近地点中间槽。
    s2[2] +=  c2o/alpha;
    s2[3] += -ci/alpha*s2o;
    s2[4] += -sio*3.0*e*s2o;

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
    if (short_flag == KeplerFull)
    {
        s2[5] += -2.0*sio*s2o;
    }
    else if (short_flag == KeplerCombined)
    {
        // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

        s2[5] += 2.25*j22a2*(-2.0*sio*s2o);
    }

    // 四倍 M 补项使用 1-alpha/2 分母，只写偏心率和近地点组合。
    double cipeaf2 = cip12*4.25*e/(1.0-0.5*alpha);
    s2[1] += cipeaf2*cmwoA[4];

    // 这行正弦补项与上一行余弦补项配对。
    s2[4] += cipeaf2*smwoA[4];

    // B 向使用 1+2alpha、1+alpha，频率符号与共振组不同。
    double af1B     = 1.0+2.0*alpha;
    double af2B     = 1.0+alpha;

    // 将两个正号分母分别与 sin²(i) 权重组合。
    double si2af1B  = si2/af1B;
    double si2af2B  = si2/af2B;

    // 3e 和额外频率因子进入这组 B 向幅值与角改正。
    double ex3      = 3.0*e;
    double yjaf1B   = 1.0-0.5/af1B;

    // 先累计 B 向的半长轴、偏心率及倾角改正。
    s2[0] += si2af1B*ex3*cmoB[1];
    s2[1] += si2af1B*3.0*cmoB[1] +si2af2B*4.5*e*cmoB[2];
    s2[2] += ex3/af1B*cmoB[1];

    // 节点和近地点中间槽对应同一 B 相位的正弦。
    s2[3] += ci*ex3/af1B*smoB[1];
    s2[4] += si2af1B*3.0*smoB[1] +si2af2B*4.5*e*smoB[2];

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
    if (short_flag == KeplerFull)
    {
        s2[5] += si2af1B*ex3*2.0*yjaf1B*smoB[1];
    }
    else if (short_flag == KeplerCombined)
    {
        // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

        s2[5] += 0.75*j22a2*eh*si2af1B*3.0*smoB[1] +2.25*j22a2*si2af1B*ex3*2.0*yjaf1B*smoB[1];
    }

    // 按 1..4 倍 M 分别保存 B 向分母与倾角权重。
    double afB[5];
    double ci1afB[5];

    // afcoef 保留每个倍频的原展开系数及下标。
    double ci2afB[5];
    double afcoef[5] = {0.0, 1.0, 1.0, 7.0, 17.0};

    for (int j = 1; j < 5; j++)
    {
        // 先得到 1+2alpha/j，再分别乘 1-cos(i) 及其平方。
        afB[j] = 1.0+2.0*alpha/(double)j;
        ci1afB[j] = afcoef[j]*cim1/afB[j];

        // 平方倾角权重与一次权重服务不同根数槽。
        ci2afB[j] = afcoef[j]*cim12/afB[j];
    }

    // 将 B 向各余弦谐项依公式系数并入半长轴和偏心率中间槽。
    s2[0] += -ci2afB[1]*eh*cmwoB[1]     +ci2afB[2]*cmwoB[2]     +ci2afB[3]*eh*cmwoB[3];
    s2[1] +=  ci2afB[1]/2.0*cmwoB[1]    -ci2afB[2]*eh*cmwoB[2]  +ci2afB[3]/6.0*cmwoB[3]     +ci2afB[4]*e/4.0*cmwoB[4];

    // 倾角、节点中间项使用一次倾角权重。
    s2[2] += -ci1afB[1]*eh*cmwoB[1]     +ci1afB[2]/2.0*cmwoB[2] +ci1afB[3]*e/6.0*cmwoB[3];
    s2[3] += -ci1afB[1]*eh*smwoB[1]     +ci1afB[2]/2.0*smwoB[2] +ci1afB[3]*e/6.0*smwoB[3];

    // 已乘 e 的近地点组合使用各倍频的正弦谐项。
    s2[4] += -ci2afB[1]/2.0*smwoB[1] -ci2afB[2]*eh*5.0*smwoB[2] +ci2afB[3]/6.0*smwoB[3]     +ci2afB[4]*e/4.0*smwoB[4];

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
    if (short_flag == KeplerFull)
    {
        s2[5] += -ci2afB[1]*e*(1.0-0.5/afB[1])*smwoB[1] +ci2afB[2]*(1.0-0.5/afB[2])*smwoB[2] +ci2afB[3]*e/3.0*(1.0-0.5/afB[3])*smwoB[3];
    }
    else if (short_flag == KeplerCombined)
    {
        // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

        s2[5] += 0.75*j22a2*eh*(-ci2afB[1]/2.0*smwoB[1] + ci2afB[3]/6.0*smwoB[3])
                +2.25*j22a2*(-ci2afB[1]*e*(1.0-0.5/afB[1])*smwoB[1]+ci2afB[2]*(1.0-0.5/afB[2])*smwoB[2]+ci2afB[3]*e/3.0*(1.0-0.5/afB[3])*smwoB[3]);
    }

    // 求和后统一乘 1.5*J22/a²，半长轴槽另乘 a。
    double ja32 = 1.5*j22a2;

    s2[0] = ja32*s2[0]*a;

    // 偏心率、倾角槽分别恢复原式规定的尺度。
    s2[1] = ja32*s2[1]/2.0;
    s2[2] = ja32*s2[2]*si;

    // 节点槽此时仍是 delta_Omega，随后再与近地点及平近点角组合。
    s2[3] = ja32*s2[3];

    // Full 在本频率组保留普通 Kepler 角改正形式，第六槽仍按原公式处理。
    if (short_flag == KeplerFull)
    {
        s2[4] = -ci*s2[3] + ja32*s2[4]/(2.0*e);
        s2[5] = -(1.0-e*e/2.0)*(s2[4]+ci*s2[3]) + 1.5*ja32*s2[5];
    }
    else if (short_flag == KeplerCombined)
    {
        // Combined 先组合平经度相关项，后续不再由近地点改正额外除以 e。

        // 第四下标保存 e*delta_omega，第五下标已组合为 delta_lambda。
        s2[4] = -ci*s2[3]*e + ja32*s2[4]/2.0;
        s2[5] = cim1*s2[3] + s2[5];
    }

    // 按顺序写出选定语义的六槽；Combined 末槽已是 delta_lambda。
    for (int j = 0; j < 6; j++)
    {
        keps2_J22[j] = s2[j];
    }

    tracer.CheckOut(__func__);
}

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
void J22SemimajorShort(bastion::Overwatch& tracer, double j22a, double e, double ci, double o22, double w, double m, double alpha, double& as2_J22, int resonance_id)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // cos²(i)、sin²(i) 在末尾用于组合三个方向的半长轴累计量。
    double ci2 = ci*ci;
    double si2 = 1.0-ci2;

    // 近地点和相对节点参考角均取两倍，与 J22 阶数一致。
    double ww = w+w;
    double oo = o22+o22;

    // 四组余弦表区分是否含 2w，以及节点相位的加减号。
    double cxadd[7], cxsub[7], cyadd[7], cysub[7];
    double afdiv[7];

    // 逐个准备 0..6 倍 M 的相角、余弦与频率比。
    for (size_t k = 0; k < 7; k++)
    {
        // y 相角只包含 kM 加减 2o22。
        double yadd = m*(double)k+oo;
        double ysub = m*(double)k-oo;

        // x 相角在对应 y 相角上再加 2w。
        double xadd = yadd+ww;
        double xsub = ysub+ww;

        // 先保存含 2w 的两组余弦，再保存不含 2w 的两组。
        cxadd[k] = cos(xadd);
        cxsub[k] = cos(xsub);

        // 不含近地点角的余弦随后进入 as20。
        cyadd[k] = cos(yadd);
        cysub[k] = cos(ysub);

        // k=0 不除以 k，其余槽保存 2alpha/k；保留原三目分支。
        afdiv[k] = (k == 0) ? 0.0 : 2.0*alpha/(double)k;

        // 共振类只控制下列槽位置零，保留原分母表达式与置零时机。
        switch (resonance_id)
        {
            // 半同步对应第一倍频，正向含 w 与不含 w 两槽同时置零。
            case ResonanceHalf:

                cxadd[1] = 0.0;
                cyadd[1] = 0.0;
                break;

            // 同步对应第二倍频，沿用正向两槽的置零处理。
            case ResonanceSynchronous:

                cxadd[2] = 0.0;
                cyadd[2] = 0.0;
                break;

            // 其他类别不增加额外置零，继续使用参考实现的处理。
            default:

                break;
        }
    }

    // 分别累计正节点相位、负节点相位及不含近地点角的三组贡献。
    double as2add = 0.0;
    double as2sub = 0.0;

    // 三组保持独立累计，最后再分别乘各自倾角权重。
    double as20   = 0.0;

    // 先累计零阶、一阶偏心率的正向谐项。
    as2add += cxadd[2]/(1.0-alpha);
    as2add += e*(-0.5*cxadd[1]/(1.0-afdiv[1]) +3.5*cxadd[3]/(1.0-afdiv[3]));

    // 负向谐项使用 1+alpha 等正号分母。
    as2sub += cxsub[2]/(1.0+alpha);
    as2sub += e*(-0.5*cxsub[1]/(1.0+afdiv[1]) +3.5*cxsub[3]/(1.0+afdiv[3]));

    // 不含近地点角的一阶项同时组合正、负节点方向。
    as20   += 1.5*e*(cyadd[1]/(1.0-afdiv[1])+cysub[1]/(1.0+afdiv[1]));

    // 二阶偏心率项依原顺序分别累加到三组中。
    double e2 = e*e;
    as2add += e2*(-2.5*cxadd[2]/(1.0-afdiv[2]) +8.5*cxadd[4]/(1.0-afdiv[4]));
    as2sub += e2*(-2.5*cxsub[2]/(1.0+afdiv[2]) +8.5*cxadd[4]/(1.0+afdiv[4]));

    // 这一行补齐不含近地点角的 e² 项。
    as20   += 2.25*e2*(cyadd[2]/(1.0-afdiv[2]) +cysub[2]/(1.0+afdiv[2]));

    // 三阶项涉及 1、3、5 倍 M，保留原展开系数和加减顺序。
    double e3 = e2*e;
    as2add += e3*(cxadd[1]/(1.0-afdiv[1])/12.0 +cxadd[5]/(1.0-afdiv[5])*4225.0/240.0 -cxadd[3]/(1.0-afdiv[3])*369.0/48.0);
    as2sub += e3*(cxsub[1]/(1.0+afdiv[1])/12.0 +cxsub[5]/(1.0+afdiv[5])*4225.0/240.0 -cxsub[3]/(1.0+afdiv[3])*369.0/48.0);

    // 不含 w 的三阶项仍按两种节点方向组合。
    as20   += e3*((cyadd[1]/(1.0-afdiv[1])+cysub[1]/(1.0+afdiv[1]))*27.0/16.0 + (cyadd[3]/(1.0-afdiv[3])+cysub[3]/(1.0+afdiv[3]))*159.0/48.0);

    // 四阶项使用 2、4、6 倍 M 的既有相角表。
    double e4 = e3*e;
    as2add += e4*(cxadd[2]/(1.0-afdiv[2])*41.0/48.0 +cxadd[6]/(1.0-afdiv[6])*1599.0/48.0 -cxadd[4]/(1.0-afdiv[4])*115.0/6.0);
    as2sub += e4*(cxsub[2]/(1.0+afdiv[2])*41.0/48.0 +cxsub[6]/(1.0+afdiv[6])*1599.0/48.0 -cxsub[4]/(1.0+afdiv[4])*115.0/6.0);

    // 补齐不含 w 的四阶项，三个累计量至此准备完成。
    as20   += e4*((cyadd[2]/(1.0-afdiv[2])+cysub[2]/(1.0+afdiv[2]))*1.75 + (cyadd[4]/(1.0-afdiv[4])+cysub[4]/(1.0+afdiv[4]))*4.8125);

    // 以 (1+cos i)²、(1-cos i)²、2sin²(i) 合成半长轴输出，外层用它替换第零槽。
    as2_J22 = 1.5*j22a*((1.0+2.0*ci+ci2)*as2add + (1.0-2.0*ci+ci2)*as2sub + 2.0*si2*as20);

    tracer.CheckOut(__func__);
}

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
    double nsg2s2_J22[], int resonance_id)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 先取 Combined 中间槽，避免与普通 Kepler 角改正混淆。
    double keps2_j22[6] = {0.0};
    J22KeplerShort(tracer, j22a2, a, e, si, ci, o22, w, m, alpha, keps2_j22, KeplerCombined, resonance_id);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 前三槽为 delta_a、delta_xi、delta_eta，偏心率向量沿 omega+Omega 投影。
    nsg2s2_J22[0] = keps2_j22[0];
    nsg2s2_J22[1] = cwo*keps2_j22[1] - swo*(keps2_j22[4]+e*keps2_j22[3]);
    nsg2s2_J22[2] = swo*keps2_j22[1] + cwo*(keps2_j22[4]+e*keps2_j22[3]);

    // h/k 组合倾角改正与节点转动；输入 h、k 已为当前 NSG2 分量。
    nsg2s2_J22[3] = 0.5*cih*co*keps2_j22[2] - k*keps2_j22[3];
    nsg2s2_J22[4] = 0.5*cih*so*keps2_j22[2] + h*keps2_j22[3];

    // Combined 末槽已是 delta_lambda，直接传递，不重复加节点或近地点改正。
    nsg2s2_J22[5] = keps2_j22[5];

    tracer.CheckOut(__func__);
}

/**
 * @brief 计算简化岁差章动与 J2 耦合的组合短周期项。
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
 * @param[in] sf sin(f)，无量纲。
 * @param[in] cf cos(f)，无量纲。
 * @param[in] fm 回绕到 [-pi,pi) 的方程中心 f-M，rad。
 * @param[in] fmef fm+e*sin(f)，rad。
 * @param[in] ar3 当前卫星距离比 (a/r)^3，无量纲。
 * @param[out] keps2_frameJ2 覆写六槽；前 3 槽为 delta_a、delta_e、delta_i，后 3 槽含义由 short_flag 决定，半长轴尺度 ae，角度 rad。
 * @param[in] short_flag KeplerFull=0 输出 Kepler 改正，KeplerCombined=2 输出供 NSG2 变换使用的组合槽。
 * @details 来自 analyt_frameinduced.cpp 的 AnaFrameJ2ShortKep。Full 输出 Kepler 改正；Combined 的第 3/4/5 槽分别为 sin(i/2)*delta_Omega、e*(delta_omega+delta_Omega)、delta_lambda。
 * 输入采用当前 GCRS 均根数估计的几何，长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
static void FrameKeplerShort(bastion::Overwatch& tracer, const FrameArguments& frameArg, double ja, double jp2,
    double e, double xe2, double xe3, double si, double ci, double sih, double cih,
    double so, double co, double s2w, double c2w, double sf, double cf, double fm, double fmef, double ar3, double keps2_frameJ2[], int short_flag)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 零槽分别初始化 sin/cos(kf) 与 sin/cos(kf+2w)，k 为真近点角倍数。
    double skf[4] = {0.0}, ckf[4] = {1.0}, sf2w[6] = {s2w}, cf2w[6] = {c2w};

    // 按原顺序用角度递推准备含 2w 的倍角项。
    for (int k = 1; k < 6; k++)
    {
        // 正余弦均读上一槽，因此不会混入本次刚更新的值。
        sf2w[k] = sf2w[k-1]*cf + cf2w[k-1]*sf;
        cf2w[k] = cf2w[k-1]*cf - sf2w[k-1]*sf;

        // 不含 2w 的倍角只需 k<=3，数组范围与公式读取范围一致。
        if (k <= 3)
        {
            // 单独递推 sin(kf)、cos(kf)，供不含近地点相角的项使用。
            skf[k] = skf[k-1]*cf + ckf[k-1]*sf;
            ckf[k] = ckf[k-1]*cf - skf[k-1]*sf;
        }
    }

    // qa 为简化岁差角，dq、de 为章动主项，角度单位 rad。
    double qa = frameArg.thetaA;
    double dq = frameArg.Dtheta;

    // 将岁差加章动按当前节点方向投影成 Aq、Bq。
    double de = frameArg.Depsilon;
    double Aq = (qa+dq)*so - de*co;

    // Bq 与 Aq 配对，后面分别进入倾角和节点相关项。
    double Bq = (qa+dq)*co + de*so;

    // 先准备偏心率平方和 cos(i) 加权投影，供各谐波共用。
    double e2   = e*e;
    double ciA  = ci*Aq;

    // 继续准备 cos(i)*Bq 与 cos(2i)*Aq，两者方向作用不同。
    double ciB  = ci*Bq;
    double c2iA = (ci*ci-si*si)*Aq;

    // 补充 sin(f-2w)、cos(f-2w)，与前面加角倍频项共同使用。
    double sfj2w = sf*c2w-cf*s2w;
    double cfj2w = cf*c2w+sf*s2w;

    // 六个槽先保存几何展开，末尾再乘 J2 幅值并按输出标志组合。
    double s2[6] = {0.0};

    // 半长轴和偏心率相关展开先计算，径向平均部分包含 ar³-xe3。
    s2[0] = ciA*(ar3-xe3-ar3*cf2w[2]) - Bq*ar3*sf2w[2];
    s2[1] = ciA*((3+0.75*e2)*ckf[1] +1.5*e*ckf[2] +e2/4.0*ckf[3] -e2/8.0*cfj2w -(0.5+11.0*e2/8.0)*cf2w[1] -2.5*e*cf2w[2] -(28.0+17.0*e2)/24*cf2w[3] -0.75*e*cf2w[4]-e2/8.0*cf2w[5]) -Bq*(-e2/8.0*sfj2w +(0.5+11.0*e2/8.0)*sf2w[1] +2.5*e*sf2w[2] +(28.0+17.0*e2)/24*sf2w[3] +0.75*e*sf2w[4] +e2/8.0*sf2w[5]);

    // 倾角和节点相关展开共同使用 f-M+e*sin(f) 与一至三倍纬度谐角。
    s2[2] = -c2iA*(0.5*e*cf2w[1] +cf2w[2]/2.0 +e*cf2w[3]/6.0) -ciB*(fmef +e*sf2w[1]/2.0 +sf2w[2]/2.0 +e*sf2w[3]/6.0);
    s2[3] = c2iA*(fmef -e*sf2w[1]/2.0 -sf2w[2]/2.0 -e*sf2w[3]/6.0) +ciB*(e*cf2w[1]/2.0 +cf2w[2]/2.0 +e*cf2w[3]/6.0);

    // 后两槽准备近地点和平近点角组合，保留谐波系数及括号。
    s2[4] = ciA*(3.0*e*fm +(3.0+2.25*e2)*skf[1] +1.5*e*skf[2] +e2/4.0*skf[3] -e2/8.0*sfj2w +(0.5-0.875*e2)*sf2w[1] -1.5*e*sf2w[2] -(28.0+11.0*e2)/24*sf2w[3] -0.75*e*sf2w[4] -e2/8.0*sf2w[5]) -Bq*(e2/8.0*cfj2w +(0.5-0.875*e2)*cf2w[1] -1.5*e*cf2w[2] -(28.0+11.0*e2)/24*cf2w[3] -0.75*e*cf2w[4] -e2/8.0*cf2w[5]);
    s2[5] = ciA*(fmef -e/2.0*sf2w[1] -sf2w[2]/2.0 -e*sf2w[3]/6.0) +Bq*(e*cf2w[1]/2.0 +cf2w[2]/2.0 +e*cf2w[3]/6.0);

    // 前三槽恢复半长轴、偏心率、倾角改正的 ja/jp2 尺度。
    keps2_frameJ2[0] = ja*si*s2[0];
    keps2_frameJ2[1] = jp2*si*s2[1];
    keps2_frameJ2[2] = jp2*s2[2];

    // Full 分支给普通角根数改正，保留原式 e 和 sin(i) 分母。
    if (short_flag == KeplerFull)
    {
        // 先构造近地点相关量，再与节点项组合出普通角改正。
        double jpies4 = jp2*si/e*s2[4];
        keps2_frameJ2[3] = jp2/si*s2[3];

        // 这两槽仍是 Full 含义，不能直接按 NSG2 的组合槽读取。
        keps2_frameJ2[4] = -ci*keps2_frameJ2[3]+jpies4;
        keps2_frameJ2[5] = -xe2*jpies4 + 3.0*jp2*si*xe2*s2[5];
    }
    else if (short_flag == KeplerCombined)
    {
        // NSG2 包装使用 Combined 分支，以半倾角和已乘 e 的组合表达角改正。

        // tih=tan(i/2)，jps3/jpis4 保存未直接除以 sin(i) 或 e 的中间量。
        double tih = sih/cih;
        double jps3 = jp2*s2[3];
        double jpis4 = jp2*si*s2[4];

        // 末三槽为 sin(i/2)*delta_Omega、e*(delta_omega+delta_Omega)、delta_lambda。
        keps2_frameJ2[3] = jp2/(2.0*cih)*s2[3];
        keps2_frameJ2[4] = e*tih*jps3 + jpis4;

        // 节点与近地点组合准备后，再合成不含 e 分母的平经度槽。
        keps2_frameJ2[5] = tih*jps3 + e/(1.0+xe2)*jpis4 + 3.0*jp2*si*xe2*s2[5];
    }
    else
    {
        // 其他输出标志不受支持，通过 tracer 报错。
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown short_flag for FrameJ2 short-period terms!", bastion::ERROR);
        return;
    }

    tracer.CheckOut(__func__);
}

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
    double nsg2s2_frameJ2[])
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 先调用 Combined 形式，底层失败时不继续组合 NSG2 输出。
    double keps2_framej2[6] = {0.0};
    FrameKeplerShort(tracer, frameArg, ja, jp2, e, xe2, xe3, si, ci, sih, cih, so, co, s2w, c2w, sf, cf, fm, fmef, ar3, keps2_framej2, KeplerCombined);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 前三槽为 delta_a，以及沿 omega+Omega 投影的 delta_xi、delta_eta。
    nsg2s2_frameJ2[0] = keps2_framej2[0];
    nsg2s2_frameJ2[1] = cwo*keps2_framej2[1] - swo*keps2_framej2[4];
    nsg2s2_frameJ2[2] = swo*keps2_framej2[1] + cwo*keps2_framej2[4];

    // h/k 结合倾角及已缩放节点项，无需再给中间节点槽补乘 sin(i/2)。
    nsg2s2_frameJ2[3] = 0.5*cih*co*keps2_framej2[2] - so*keps2_framej2[3];
    nsg2s2_frameJ2[4] = 0.5*cih*so*keps2_framej2[2] + co*keps2_framej2[3];

    // 末槽已是 delta_lambda，直接传递。
    nsg2s2_frameJ2[5] = keps2_framej2[5];

    tracer.CheckOut(__func__);
}
} // namespace orbit::averaging_detail
