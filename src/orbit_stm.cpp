/**
 * @file orbit_stm.cpp
 * 2026-10-21：迁入导师 NSG2 几何偏导和六状态解析 STM，改用 Eigen。
 * 数学来源：NAO/gaussLib/orbit_geometry.cpp、orbit_statetransitionmatrix.cpp，Jingshi Tang。
 * 原文件 Copyright (c) 2020 国家航天局空间碎片监测与应用中心，All rights reserved.
 * 保留来源说明；本地科研复刻，不改变原始来源的权利归属。
 */
#include "orbit_stm.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

namespace determination
{
StateMatrix AnalyticStmNSG2(const StateVector &state, double dt, const AnalyticStmSetting &setting)
{
    // NSG2 使用归一化状态和时间；先检查椭圆轨道与非退化倾角定义域，非法输入保持原异常行为。
    if (!state.allFinite() || state[0] <= 0 || state[1] * state[1] + state[2] * state[2] >= 1 ||
        state[3] * state[3] + state[4] * state[4] >= 1 || !std::isfinite(dt))
        throw std::runtime_error("Invalid elliptic NSG2 state for analytic STM.");
    // Phi_ij = dz_i(t)/dz_j(t0)。t=t0 时为单位阵。
    // 二体运动中仅 lambda 随时间增加 n*dt；n=a^(-3/2)，故 d(lambda)/da=-3*n*dt/(2*a)。
    // 下面的 J2、日月项对平均根数长期变化率求偏导，再乘 dt；它们是一次近似，
    // 不是逐步积分的全动力学变分方程。不能把打开这些项等同于支持所有摄动力偏导。
    // 每个数组是一行对初始六根数 [a,xi,eta,h,k,lambda] 的偏导，先清零非二体贡献。
    std::array<double, 6> dad{};
    std::array<double, 6> dxd{};
    // dxd/ded 分别对应 xi/eta；名称沿用来源公式，以便逐项核对。
    std::array<double, 6> ded{};
    std::array<double, 6> dhd{};
    // dhd/dkd 对应轨道平面根数，dld 对应平均经度 lambda。
    std::array<double, 6> dkd{};
    std::array<double, 6> dld{};
    // a 已除以地球半径，n 采用归一化时间的倒数，以下 n*dt 无量纲。
    const double a = state[0];
    const double n = 1.0 / a / sqrt(a);

    // 二体情况下 a、xi、eta 不随时间变化，故相应的同分量偏导为一。
    dad[0] = 1.0;
    dxd[1] = 1.0;
    ded[2] = 1.0;
    // h、k 也保持不变；lambda 还包含平均角速度对 a 的变化。
    dhd[3] = 1.0;
    dkd[4] = 1.0;
    dld[0] = -1.5 * n / a * dt;
    dld[5] = 1.0;

    // NSG2 顺序为 a、xi、eta、h、k、lambda；先准备偏心率和轨道倾角的组合量。
    const double xi = state[1];
    const double eta = state[2];
    // h、k 编码轨道平面方向；平方和等于 sin²(i/2)。
    const double h = state[3];
    const double k = state[4];

    const double xi2 = xi * xi;
    const double et2 = eta * eta;
    // 后面按原公式重复使用交叉项与平面根数的平方，不重新展开乘法。
    const double xieta = xi * eta;
    const double h2 = h * h;
    const double k2 = k * k;
    // ci、si2 分别为 cos(i)、sin²(i)，无需先反解倾角再求三角函数。
    const double hk = h * k;
    const double ci = 1.0 - 2.0 * (h2 + k2); // cos(i) = 1-2(h^2+k^2)
    const double si2 = 1.0 - ci * ci;        // sin^2(i)
    // e2 与 je2 为偏心率平方及其余量，后续半通径和长期变化率共用。
    const double e2 = xi2 + et2;             // e2 = xi^2+eta^2
    const double je2 = 1.0 - e2;             // 1-e^2 = 1-xi^2-eta^2

    // J2 长期项：先计算升交点、近地点与平近点角的变化率，再对初始根数求偏导。
    if (setting.use_J2)
    {
        const double p = a * je2;                     // p = a(1-e^2)
        const double j2p2 = 1.5 * setting.j2 / p / p; // 3J2/2p^2
        const double jpndt = j2p2 * n * dt;           // (3J2/2p^2)n*dt

        // o1、w1、m1 分别对应 Omega、omega、M 的摄动变化率；wt1、l1 为其组合。
        const double o1 = -j2p2 * n * ci;
        const double w1 = j2p2 * n * (2.0 - 2.5 * si2);
        const double m1 = j2p2 * n * (1.0 - 1.5 * si2) * sqrt(je2);
        // wt1 合并节点与近地点进动，l1 再加平近点角项，分别进入偏心率矢量和平均经度。
        const double wt1 = o1 + w1;
        const double l1 = wt1 + m1;
        // 收集对 a、偏心率与倾角求导时反复出现的系数，并乘传播时长 dt。
        const double a35 = 3.5 / a;
        const double je4 = 4.0 / je2;
        const double hk4 = 4.0 * hk;
        // 累计进动角用于一阶状态转移，j5ci 则来自倾角方向的求导。
        const double o1dt = o1 * dt;
        const double wt1dt = wt1 * dt;
        // j5ci 与 mwtdt 分别收集倾角导数和偏心率导数中的公共系数。
        const double j5ci = 1.0 - 5.0 * ci;               // 1-5*cos(i) = 1-5[1-2(h^2+k^2)]
        const double mwtdt = (3.0 * m1 + 4.0 * wt1) * dt; // (3M1+4wt1)*dt
        // 以下组合用于 xi/eta、h/k 和 lambda 各行，避免重复展开同一偏导表达式。
        const double cidt4 = j5ci * jpndt * 4.0;                        // 4*{1-5[1-2(h^2+k^2)]}*(3J2/2p2)n*dt
        const double ei1xx = -4.0 * (3.0 * sqrt(je2) + 5.0) * ci + 4.0; // -4[(3sqrt(1-e2)+5]cos(i)+4
        // je4dt 和 exj4t 对应近地点进动的偏心率导数，eo1dt 则用于轨道平面两行。
        const double je4dt = je4 * wt1dt;                               // 4/(1-e2)*wt1*dt
        const double exj4t = -xieta * je4dt;                            // -eta*xi*4/(1-e2)*wt1*dt
        const double eo1dt = je4 * o1dt;

        // xi、eta 两行：轨道偏心率矢量随近地点和升交点进动而变化。
        dxd[0] += -eta * (-a35) * wt1dt;
        dxd[1] += -exj4t;
        dxd[2] += -wt1dt - et2 * je4dt;
        // xi 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dxd[3] += -eta * h * cidt4;
        dxd[4] += -eta * k * cidt4;

        ded[0] += xi * (-a35) * wt1dt;
        ded[1] += wt1dt + xi2 * je4dt;
        ded[2] += exj4t;
        // eta 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        ded[3] += xi * h * cidt4;
        ded[4] += xi * k * cidt4;

        // h、k 两行：升交点进动引起的轨道平面方向变化。
        dhd[0] += -k * (-a35) * o1dt;
        dhd[1] += -k * xi * eo1dt;
        dhd[2] += -k * eta * eo1dt;
        // h 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dhd[3] += -hk4 * jpndt;
        dhd[4] += -o1dt - 4.0 * k2 * jpndt;

        dkd[0] += h * (-a35) * o1dt;
        dkd[1] += h * xi * eo1dt;
        dkd[2] += h * eta * eo1dt;
        // k 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dkd[3] += o1dt + 4.0 * h2 * jpndt;
        dkd[4] += hk4 * jpndt;

        // lambda 行：Omega、omega、M 的长期变化共同改变平均经度。
        dld[0] += -a35 * l1 * dt;
        dld[1] += xi / je2 * mwtdt;
        dld[2] += eta / je2 * mwtdt;
        // lambda 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dld[3] += h * ei1xx * jpndt;
        dld[4] += k * ei1xx * jpndt;
    }

    // 日月引力采用平均长期项；这里的固定质量比、距离和轨道倾角来自解析 STM 设置。
    if (setting.use_lunisolar)
    {
        // 质量比除距离三次方给出平均引力强度；距离均按地球半径归一化。
        const double lunarBeta =
            setting.lunar_mass / setting.lunar_distance / setting.lunar_distance / setting.lunar_distance;
        const double solarBeta =
            setting.solar_mass / setting.solar_distance / setting.solar_distance / setting.solar_distance;
        // 固定日月轨道倾角只进入平均长期项，不随当前卫星状态更新。
        const double lunarSinI = setting.lunar_sin_i;
        const double solarSinI = setting.solar_sin_i;
        double W =
            lunarBeta * (1.0 - 1.5 * lunarSinI * lunarSinI) + solarBeta * (1.0 - 1.5 * solarSinI * solarSinI);
        // 保留先求日月加权和、再逐次乘 a 的原次序，避免重排解析近似的舍入路径。
        W = W * 0.75 * a * a * a;

        // 由日月平均引力强度 W 得到各长期变化率，再整理对偏心率、倾角的导数。
        const double w2p = W / sqrt(je2) * n;
        const double o2p = -w2p * (1.0 + 1.5 * e2);
        const double m2p = -W * (7.0 / 3.0 + e2) * n;
        // w2、o2、m2 分别对应近地点、节点与平近点角的平均变化率。
        const double w2 = w2p * (2.0 - 2.5 * si2 + e2 / 2.0);
        const double o2 = o2p * ci;
        const double m2 = m2p * (1.0 - 1.5 * si2);

        // 以下导数系数来自平均变化率对偏心率及平面根数求导，不是高阶时间积分项。
        const double o2pp = (1.0 - 0.375 * e2) / (1.0 + 1.5 * e2) * o2;
        const double wt2 = o2 + w2;
        const double wt2p = (8.0 - 3.0 * e2) / (2.0 + 3.0 * e2) * o2 + (3.0 - 2.5 * si2 - e2 / 2.0) * w2p;
        // 合并节点与近地点的偏导，再加入 M 项，形成 lambda 行所需系数。
        const double wt2pp = o2p + 5.0 * ci * w2p;
        const double l2p = wt2p + 6.0 * je2 / (7.0 + 3.0 * e2) * m2;
        const double l2pp = wt2pp + 3.0 * ci * m2p;

        // 长期变化率乘以 dt，得到本次传播区间内的一阶累计贡献。
        const double o2dt = o2 * dt;
        const double l2dt = (o2 + w2 + m2) * dt;
        const double o2pdt = o2p * dt;
        // 将各导数系数分别乘 dt；这些量只在本次 STM 近似中使用。
        const double o2ppdt = o2pp * dt;
        const double wt2dt = wt2 * dt;
        const double wt2pdt = wt2p * dt;
        // lambda 的导数系数也乘同一传播时长，保持所有行参考同一拟合历元。
        const double wt2ppdt = wt2pp * dt;
        const double l2pdt = l2p * dt;
        const double l2ppdt = l2pp * dt;

        // h、k 的四倍与偏心率余量比值是下面各列的公共因子。
        const double hx4 = h * 4.0;
        const double kx4 = k * 4.0;
        const double a23 = 1.5 / a;
        // xi/je2、eta/je2 对应偏心率方向求导，交叉项负责 xi 与 eta 的耦合。
        const double xije2 = xi / je2;
        const double etaje2 = eta / je2;
        const double xietaje2 = xieta / je2;

        // 按 xi、eta、h、k、lambda 行依次累加日月贡献；列始终对应六个初始根数。
        dxd[0] += -eta * a23 * wt2dt;
        dxd[1] += -xietaje2 * wt2pdt;
        dxd[2] += -wt2dt - et2 / je2 * wt2pdt;
        // xi 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dxd[3] += hx4 * eta * wt2ppdt;
        dxd[4] += kx4 * eta * wt2ppdt;

        // eta 行与 xi 行保留来源公式的正负号对应关系，分别累加两方向的进动响应。
        ded[0] += xi * a23 * wt2dt;
        ded[1] += wt2dt + xi2 / je2 * wt2pdt;
        ded[2] += xietaje2 * wt2pdt;
        // eta 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        ded[3] += -hx4 * xi * wt2ppdt;
        ded[4] += -kx4 * xi * wt2ppdt;

        // h、k 两行描述轨道平面变化；前三列是对 a、xi、eta 的敏感度。
        dhd[0] += -k * a23 * o2dt;
        dhd[1] += -kx4 * xije2 * o2ppdt;
        dhd[2] += -kx4 * etaje2 * o2ppdt;
        // h 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dhd[3] += 4.0 * hk * o2pdt;
        dhd[4] += -o2dt + 4.0 * k2 * o2pdt;

        // k 行继续使用同一节点变化率，不能与 h 行混用偏心率矢量的 wt2 系数。
        dkd[0] += h * a23 * o2dt;
        dkd[1] += hx4 * xije2 * o2ppdt;
        dkd[2] += hx4 * etaje2 * o2ppdt;
        // k 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dkd[3] += o2dt - 4.0 * h2 * o2pdt;
        dkd[4] += -4.0 * hk * o2pdt;

        // 平均经度行合并节点、近地点和平近点角贡献；最后两列仍为 h/k 耦合。
        dld[0] += a23 * l2dt;
        dld[1] += xije2 * l2pdt;
        dld[2] += etaje2 * l2pdt;
        // lambda 行的第 3、4 列对应初始 h、k；下标从零开始，以下为平面方向耦合。
        dld[3] += -hx4 * l2ppdt;
        dld[4] += -kx4 * l2ppdt;
    }

    // 将各根数的导数行装入 Phi，供观测偏导链 H*G*Phi 使用。
    StateMatrix phi;
    for (int j = 0; j < 6; ++j)
    {
        phi(0, j) = dad[j];
        phi(1, j) = dxd[j];
        phi(2, j) = ded[j];
        // 后三行对应 h、k、lambda；第 j 列始终是同一个初始根数。
        phi(3, j) = dhd[j];
        phi(4, j) = dkd[j];
        phi(5, j) = dld[j];
    }
    return phi;
}

} // namespace determination
