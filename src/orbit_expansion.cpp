/**
 * @file orbit_expansion.cpp
 * @brief Kaula 倾角函数、导数及 Hansen 偏心率展开。
 */
#include "orbit_expansion.hpp"
#include "math_basic.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

double Flmp(int l, int m, int p, double i)
{
    // 标量入口与导数入口共用同一半角公式；这里仅把函数值 F 返回给调用方。
    double F;
    double dF;
    double ddF;
    // 导数由共用函数写入局部变量，此标量接口不额外计算另一套公式。
    FlmpDerivatives(l, m, p, i, F, dF, ddF);
    return F;
}

void FlmpDerivatives(int l, int m, int p, double i, double& F, double& dF, double& ddF)
{
    // i 为 rad；整数索引须处于当前 Kaula 倾角函数定义域，非法值保持抛 invalid_argument。
    if (l < 0 || m < 0 || m > l || p < 0 || p > l || !std::isfinite(i))
        throw std::invalid_argument("Invalid Kaula inclination-function arguments");
    // 将原有 sin(i)^A(1+cos(i))^B 化为半角的非负整数幂，消除零倾角的 0/0。
    const double sh = sin(i/2);
    const double ch = cos(i/2);
    // scale 为每个求和项共有的阶乘系数；保留分母和除法的原有结合次序。
    const double scale = Factorial(l+m)/(pow(2., l)*Factorial(p)*Factorial(l-p));
    F = dF = ddF = 0;
    // k 的上下限保证原始半角幂 A、B 非负，不对零半角计算负幂。
    for (int k = std::max(0, l-m-2*p); k <= std::min(l-m, 2*l-2*p); ++k)
    {
        const int A = -l+m+2*p+2*k;
        const int B = 2*l-A;
        // parity 决定本项符号，c 包含符号、公共尺度和两个二项式系数。
        const int parity = k+(l-m+(l+m)%2)/2;
        const double c = (parity%2 ? -1. : 1.)*scale*
            Combinatorial(2*l-2*p,k)*Combinatorial(2*p,l-m-k);
        // 先求半角幂乘积，再乘 c，保持原来 c*term(A,B) 的运算树。
        F += c*(pow(sh,A)*pow(ch,B));
        // 一阶导数由 sin(i/2)、cos(i/2) 的乘积法则给出；指数为零时跳过相应项。
        if (A) dF += .5*c*A*(pow(sh,A-1)*pow(ch,B+1));
        if (B) dF -= .5*c*B*(pow(sh,A+1)*pow(ch,B-1));
        // 二阶中间项始终使用原非负幂；移低两次幂的两项分别检查 A、B，保护零倾角端点。
        ddF -= .25*c*(2*A*B+A+B)*(pow(sh,A)*pow(ch,B));
        if (A>=2) ddF += .25*c*A*(A-1)*(pow(sh,A-2)*pow(ch,B+2));
        if (B>=2) ddF += .25*c*B*(B-1)*(pow(sh,A+2)*pow(ch,B-2));
    }
}

// Hansen 标量展开公式改编自 Jingshi Tang, orbit_fundamental.cpp (2020)。
// Copyright (c) 2020 国家航天局空间碎片监测与应用中心。All rights reserved.
// 为保持短周期模型一致，导数采用参考程序的低阶截断，见公共接口说明。
void HansenSeries(int l, int p, int q, double ecc, int expan_order, double& X, double& dX)
{
    // 偏心率限定为椭圆轨道 0<=e<1；此截断表只覆盖 0..4 阶以及 |q-p|<=4。
    if (!std::isfinite(ecc) || ecc<0 || ecc>=1 || expan_order<0 || expan_order>4 || std::abs(q-p)>4)
        throw std::invalid_argument("Invalid Hansen-series arguments");
    // 未启用的幂保持零；各阶仍按原先的逐次相乘顺序形成，不改为通用 pow。
    double e1 = 0.0;
    double e2 = 0.0;
    // e3、e4 只在展开阶次达到对应值时启用，以此控制 X 中保留的偏心率项。
    double e3 = 0.0;
    double e4 = 0.0;
    switch (expan_order)
    {
        // 一阶仅启用 e，±1 谐波可出现最低阶贡献。
        case 1:
            e1 = ecc;
            break;
        // 二阶加入 e²，供中心谐波及 ±2 谐波使用。
        case 2:
            e1 = ecc;
            e2 = e1*e1;
            break;
        // 三阶进一步启用 e³，保留原递推乘积，供 ±1 与 ±3 谐波使用。
        case 3:
            e1 = ecc;
            e2 = e1*e1;
            // 三次幂仍由平方再乘一次 e 得到，供本阶奇次侧带使用。
            e3 = e2*e1;
            break;
        // 四阶在前三项之外启用 e⁴，中心、±2 与 ±4 谐波由此获得相应项。
        case 4:
            e1 = ecc;
            e2 = e1*e1;
            // 最高一项继续乘 e 得到，避免改变原多项式所用幂值的舍入顺序。
            e3 = e2*e1;
            e4 = e3*e1;
            break;
        // 零阶保留所有 e 次幂为零，但 dX 仍按下面来源程序的低阶截断公式返回。
        case 0:
            break;
        default: throw std::invalid_argument("Invalid Hansen-series order or index");
    }

    // 将整数索引转为 double，再准备多项式系数反复使用的 l、p 幂和交叉项。
    const double l1 = static_cast<double>(l);
    const double p1 = static_cast<double>(p);
    const double l2 = l1*l1;
    const double p2 = p1*p1;
    // 三、四次幂按低阶幂继续相乘，保留原系数计算的顺序。
    const double l3 = l2*l1;
    const double l4 = l2*l2;
    // p 的高次幂与 l 分开保存，便于核对正负谐波的符号差异。
    const double p3 = p2*p1;
    const double p4 = p2*p2;
    const double lp = p1*l1;
    // 以下交叉项仍对应 l*p²、l²*p、l²*p²，不在展开式内部重新组合。
    const double lp2  = l1*p2;
    const double l2p  = l2*p1;
    const double l2p2 = l2*p2;
    // 最后两项用于四阶系数中的 l³p 与 lp³。
    const double l3p  = l3*p1;
    const double lp3  = l1*p3;

    // q-p 选择相对中心谐波的偏移；X 保留所设阶次，而 dX 按参考程序使用单独低阶截断。
    // 不把 dX 改为 X 的完整多项式导数，以保持短周期模型与解析定轨的既有约定。
    switch (q-p)
    {
        // 中心谐波包含 1、e²、e⁴；dX 仅保留对 e² 项求导后的一次项。
        case 0:
            X  = 1.0 + e2*(l2+l1-4.0*p2)/4.0 + e4*(l4-2.0*l3-l2-8.0*l2p2+2.0*l1-9.0*p2+16.0*p4)/64.0;
            dX = e1*(l2+l1-4.0*p2)/2.0;
            break;
        // 正一阶侧带含 e、e³，dX 保留常数与 e²。
        case 1:
            X  = (2.0*p1-l1)*e1/2.0 - e3*(l3-l2-2.0*l2p-3.0*l1-5.0*lp-4.0*lp2+2.0*p1+10.0*p2+8.0*p3)/16.0;
            dX = p1 - l1/2.0 -3.0*e2/16*(l2*(l1-1.0-2.0*p1)-(3.0+5.0*p1+4.0*p2)*l1+p1*(2.0+10.0*p1+8.0*p2));
            break;
        // 负一阶侧带与正侧带有不同奇次索引符号，dX 同样保留常数与 e²。
        case -1:
            X  = (-l1-2.0*p1)*e1/2.0 - e3*(l3-l2+2.0*l2p-3.0*l1+5.0*lp-4.0*lp2-2.0*p1+10.0*p2-8.0*p3)/16.0;
            dX = -p1 - l1/2.0 -3.0*e2/16*(l2*(l1-1.0+2.0*p1)-(3.0-5.0*p1+4.0*p2)*l1-p1*(2.0-10.0*p1+8.0*p2));
            break;
        // 正二阶侧带含 e²、e⁴，dX 仅保留由 e² 得到的一次项。
        case 2:
            X  = (l2-4.0*lp-3.0*l1+4.0*p2+5.0*p1)*e2/8.0 + (l4-6.0*l3-4.0*l3p-l2+3.0*l2p+22.0*l1+47.0*lp+48.0*lp2+16.0*lp3-22.0*p1-64.0*p2-60.0*p3-16.0*p4)*e4/96.0;

            dX = (l2-4.0*lp-3.0*l1+4.0*p2+5.0*p1)*e1/4.0;
            break;
        // 负二阶侧带沿用对应负谐波系数，dX 仍截断到一次项。
        case -2:
            X  = (l2+4.0*lp-3.0*l1+4.0*p2-5.0*p1)*e2/8.0 + (l4-6.0*l3+4.0*l3p-l2-3.0*l2p+22.0*l1-47.0*lp+48.0*lp2-16.0*lp3+22.0*p1-64.0*p2+60.0*p3-16.0*p4)*e4/96.0;

            dX = (l2+4.0*lp-3.0*l1+4.0*p2-5.0*p1)*e1/4.0;
            break;
        // 负三阶侧带从 e³ 开始；来源近似将其 dX 截为零。
        case -3:
            X  = -e3*(l3-9.0*l2+6.0*l2p+17.0*l1-33.0*lp+12.0*lp2+26.0*p1-30.0*p2+8.0*p3)/48.0;
            dX = 0.0;
            break;
        // 正三阶侧带也从 e³ 开始；保持 dX=0 的既有低阶近似。
        case 3:
            X  = -e3*(l3-9.0*l2-6.0*l2p+17.0*l1+33.0*lp+12.0*lp2-26.0*p1-30.0*p2-8.0*p3)/48.0;
            dX = 0.0;
            break;
        // 负四阶侧带只含 e⁴；此处仍不加入其三次导数项。
        case -4:
            X  = e4*(l4-18.0*l3+8.0*l3p+95.0*l2-102.0*l2p+24.0*l2p2-142.0*l1+330.0*lp-192.0*lp2+32.0*lp3-206.0*p1+283.0*p2-120.0*p3+16.0*p4)/384.0;
            dX = 0.0;
            break;
        // 正四阶侧带只含 e⁴；与负四阶相同，dX 按来源约定为零。
        case 4:
            X  = e4*(l4-18.0*l3-8.0*l3p+95.0*l2+102.0*l2p+24.0*l2p2-142.0*l1-330.0*lp-192.0*lp2-32.0*lp3+206.0*p1+283.0*p2+120.0*p3+16.0*p4)/384.0;
            dX = 0.0;
            break;
        default: throw std::invalid_argument("Invalid Hansen-series order or index");
    }

}
