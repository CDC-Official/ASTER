/**
 * @file analyt_averaging.cpp
 * @brief 在给定摄动模型下，从 GCRS 密切 NSG2 根数迭代扣除解析短周期项。
 * @details 输入历元为 UTC MJD，岁差章动使用 TT，日月平位置及格林尼治角取 UT1=UTC。
 * 长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad；先准备模型，再按固定次数更新均根数估计。
 * 各摄动公式在 analyt_gravity、analyt_lunisolar 和 analyt_solar_radiation 中实现，错误通过 tracer 返回。
 */
#include "orbit_expansion.hpp"
#include "analyt_averaging.hpp"
#include "analyt_terms.hpp"
#include "model_ephemeris_analytic.hpp"
#include "utils_constants.hpp"
#include "data_gravityfield.hpp"
#include "math_basic.hpp"
#include "orbit_basic.hpp"
#include "orbit_statedef.hpp"
#include "time_transform.hpp"

#include <array>
#include <cmath>
#include <exception>
#include <iomanip>
#include <sstream>

using namespace constants::Scale;
using namespace constants::Conversion;
using constants::PI;

namespace orbit {
namespace {
using namespace averaging_detail;
using SixElements = std::array<double, 6>;

// 田谐解析模型使用 IAU 1982 平恒星时角速度，按下方 GMST 多项式在 J2000 的斜率计算。
// 系数可对照 SOFA gmst82；每儒略世纪的转角包含日周项和平均恒星时进动项。
// 该角速度包含平春分点运动，单位为每归一化时间的 rad。
const double GmstRatePerCentury = 230121.675315368986346;
const double RotationRate = GmstRatePerCentury * UNIT_TIME_S / (constants::DAYS_PER_CENTURY * constants::SECONDS_PER_DAY);

// 两表同一槽位给出 alpha=omega_E/n 的有理数分子和分母，槽位 0 不参与搜索。
// 槽位 2 对应同步共振 alpha=1，其余槽位给出各组地球自转/卫星轨道频率比。
constexpr std::array<int, 15> ResonanceNumerator{0, 1, 1, 3, 1, 2, 1, 3, 1, 1, 1, 1, 1, 1, 2};
constexpr std::array<int, 15> ResonanceDenominator{1, 2, 1, 2, 3, 3, 4, 4, 7, 9, 12, 13, 14, 15, 29};

/**
 * @brief 将角度回绕到 [0,2pi)。
 * @param[in] angle 待回绕的 rad 角度。
 * @return 回绕后的 rad 角度。
 * @details 采用 fmod 后为负余数补一个整周，用于相角及平经度。
 */
double Wrap_Angle(double angle)
{
    // fmod 允许负余数，因此负值还需补一个整周。
    angle = fmod(angle, 2 * PI);

    if (angle < 0)
    {
        return angle + 2 * PI;
    }

    // 余数非负时已经落入目标角度区间。
    return angle;
}

/**
 * @brief 按输入密切半长轴选择地球自转与卫星平运动的共振类别。
 * @param[in] a 密切半长轴，以 ae 归一化。
 * @return 频率表中首次命中的类别 1..14；未命中返回 -1。
 * @details alpha=omega_E/n；同步类别 2 的频率比容差为 0.005，其余为 0.002，一次变换内固定该类别。
 */
int Find_Resonance(double a)
{
    // n=a^(-3/2)，故 alpha=omega_E/n；这里的 a 已按 ae 归一化。
    const double alpha = RotationRate * a * sqrt(a);

    for (int j = 1; j < 15; ++j)
    {
        // 容差施加于无量纲频率比，同步共振采用 0.005，其余采用 0.002。
        double threshold = .002;

        if (j == 2)
        {
            threshold = .005;
        }

        // 按频率表顺序返回首次命中的类别；不命中时由下方返回 -1。
        const double reference_ratio = double(ResonanceNumerator[j]) / ResonanceDenominator[j];

        if (std::abs(alpha - reference_ratio) <= threshold)
        {
            return j;
        }
    }

    // 所有参考频率比都未命中时，后续按非预定义共振处理。
    return -1;
}

/**
 * @brief 从 GGM03C 完全归一化系数取得未归一化重力幅值及相位。
 * @param[in] l 重力次数，本变换使用 2..10。
 * @param[in] m 重力阶数，0..l。
 * @param[out] phase m>0 时为 atan2(S_lm,C_lm)，rad；m=0 时为零。
 * @return m=0 时返回有符号 C_l0；m>0 时返回 sqrt(C_lm^2+S_lm^2)，均已去归一化。
 * @details 使用 GGM03C_arr 及阶乘归一化因子，为 J2/J22 和 Kaula 求和准备同一组重力系数。
 */
double GravityAmplitudeAndPhase(int l, int m, double& phase)
{
    // GGM03C 系数从二次项开始按三角形下标存放，每槽依次为 C_lm、S_lm。
    const auto& cs = constants::GGM03C::COEFFICIENTS[l * (l + 1) / 2 - 3 + m];
    double order_factor = 2.;

    // 带谐没有独立正弦系数，归一化因子中不采用田谐的倍数 2。
    if (m == 0)
    {
        order_factor = 1.;
    }

    // 归一化因子由次数、阶数及阶乘比构成。
    const double normalization = sqrt(order_factor * (2 * l + 1) * Factorial(l - m) / Factorial(l + m));

    // 球谐势项包含配套 GM 和参考半径的 l 次幂，统一换成程序积分尺度。
    const double field_scale = (constants::GGM03C::GM / UNIT_GM_M3_S2) *
        pow(constants::GGM03C::REFERENCE_RADIUS / UNIT_LENGTH_M, l);

    if (m == 0)
    {
        // m=0 直接保留有符号 C_l0，不取幅值绝对值。
        phase = 0;
        return cs[0] * normalization * field_scale;
    }

    // 田谐系数将正弦、余弦两部分合成幅值与相位，后面统一参与 Kaula 相角。
    phase = atan2(cs[1], cs[0]);
    return hypot(cs[0], cs[1]) * normalization * field_scale;
}

/**
 * @brief 准备解析短周期公式使用的日月距离、方向及频率。
 * @param[in,out] tracer 记录辅助状态恢复轨道根数时的错误。
 * @param[in] mjd_ut1 日月平位置的 UT1 MJD，本瞬转平模型取 UT1=UTC。
 * @param[in] moon true 采用月球平位置、GM 和平运动，false 采用太阳。
 * @return 当前第三体的 ae 归一化距离、rad 轨道方向、质量距离幂及 1/UNITTIME 平运动；出错时返回零值并设置 tracer。
 * @details 复用 Cart_MeanSun/Moon；辅助速度仅用于在 mu=1 的归一化转换中恢复方向，并非真实日月速度。含卫星角度的组合量由每轮迭代另行更新。
 */
ThirdBodyArguments Prepare_ThirdBody(bastion::Overwatch& tracer, double mjd_ut1, bool moon)
{
    // 同时取原六根数转换已算出的辅助速度；单独一个位置不能确定轨道平面与近地点方向。
    double r[3], v[3];

    if (moon)
    {
        Cart_MeanMoon(mjd_ut1, r, v);
    }
    else
    {
        // 太阳使用同一 UT1 历元的平位置及辅助速度。
        Cart_MeanSun(mjd_ut1, r, v);
    }

    // v 不是解析历表对历元的真实导数或物理日月速度；仅在原转换的 mu=1 约定下恢复根数。
    // 两者均已归一化，因此直接构造 Cartesian，再由现有转换求出方向，不能把 v 换成真实速度。
    Cartesian cartesian(r[0], r[1], r[2], v[0], v[1], v[2], true);
    Keplerian keplerian;
    cartesian.ToKeplerian(tracer, keplerian);
    if (tracer.Code() != bastion::ERR_NOERROR)
    {
        return {};
    }

    // 第三体 Kepler 六槽为 [a,e,i,Omega,omega,M]；E 由平近点角和偏心率求得。
    std::array<double, 6> elements;
    keplerian.ToArray(elements.data());
    const double E = MeanAnomaly2EccentricAnomaly(elements[5], elements[1]);

    // 由第三体偏近点角恢复 ae 归一化地心距离。
    ThirdBodyArguments body{};
    body.r = elements[0] * (1 - elements[1] * cos(E));

    // 日月 GM 与数值传播共用 DE440 参数，先除以地球 GM，再按距离构成摄动强度。
    const double normalized_gm = moon ? NORMALIZED_MOON_GM_TT : NORMALIZED_SUN_GM_TT;

    // beta2 是归一化 GM 除以当前距离三次幂。
    body.beta2 = normalized_gm / pow(body.r, 3);

    // 更高次项依次增加一个距离分母；不在公式调用时重新改变归一化尺度。
    body.beta3 = body.beta2 / body.r;
    body.beta4 = body.beta3 / body.r;

    // i、Omega 与纬度幅角 u 描述第三体轨道方向，角度统一为 rad。
    const double i = elements[2];
    body.Omega = elements[3];
    body.u = elements[4] + EccentricAnomaly2TrueAnomaly(E, elements[1]);

    // 第三体平运动以归一化时间的倒数表示；双快变量月球公式会读取它。
    if (moon)
    {
        body.n = 13.065 * DEGREE_TO_RAD * UNIT_TIME_S / constants::SECONDS_PER_DAY;
    }
    else
    {
        // 太阳平运动取 365.25 日周期，再转为归一化时间单位。
        body.n = 2 * PI * UNIT_TIME_S / (365.25 * constants::SECONDS_PER_DAY);
    }

    // 后面的方向组合直接使用 sin(i)、1+cos(i)、1-cos(i)。
    body.si = sin(i);
    body.cip1 = 1 + cos(i);
    body.cim1 = 1 - cos(i);

    // 几何三角量在卫星根数准备好后才组合，这里只返回第三体自身参数。
    return body;
}

/**
 * @brief 按当前卫星节点和近地点更新第三体组合相角。
 * @param[in,out] body 读取第三体 Omega、u，覆写 s_pp/c_pp 等和差角表，自身距离及轨道方向保持不变。
 * @param[in] Omega 当前卫星 GCRS 升交点赤经，rad。
 * @param[in] omega 当前卫星近地点幅角，rad。
 * @details theta=Omega-Omega_body；分别准备 omega+/-u+/-theta 和 omega+/-u 的正余弦，供日月引力、光压和地影使用。
 */
void Update_ThirdBodyAngles(ThirdBodyArguments& body, double Omega, double omega)
{
    // theta 为卫星与第三体升交点差，plus/minus 分别为 omega 加减第三体纬度幅角。
    const double theta = Omega - body.Omega;
    const double plus = omega + body.u;
    const double minus = omega - body.u;

    // pp/pm/mp/mm 依次对应 plus+theta、minus+theta、plus-theta、minus-theta。
    body.s_pp = sin(plus + theta);
    body.c_pp = cos(plus + theta);

    // 减去第三体纬度幅角的这支仍加 theta，区分近地点和第三体相位。
    body.s_pm = sin(minus + theta);
    body.c_pm = cos(minus + theta);

    // 减去节点差的一组与上组共同描述第三体方向相对当前卫星轨道的投影。
    body.s_mp = sin(plus - theta);
    body.c_mp = cos(plus - theta);

    // 第二组同时减去 u 和 theta，与另外三组和差角共同描述第三体方向。
    body.s_mm = sin(minus - theta);
    body.c_mm = cos(minus - theta);

    // 不含节点差的组合供第三体和圆柱地影公式共同使用。
    body.s_p = sin(plus);
    body.c_p = cos(plus);

    // 最后一对只保留 omega-u，供无节点差的方向投影使用。
    body.s_m = sin(minus);
    body.c_m = cos(minus);
}

/**
 * @brief 计算三次及以上带谐和田谐引力的 NSG2 短周期改正。
 * @param[in,out] tracer 记录未处理田谐小分母的警告。
 * @param[in] options 读取 gravity_degree、hansen_order 及 use_jl/use_jlm，控制次数、偏心率截断和谐项选择。
 * @param[in] keplerian 当前 GCRS 均根数估计 [a,e,i,Omega,omega,M]，半长轴尺度 ae，四个角 rad。
 * @param[in] greenwich_angle 同历元的格林尼治角，UT1=UTC，rad。
 * @param[in] resonance 由密切输入固定的频率表类别，1..14 或未命中的 -1。
 * @return 覆写意义的六分量值 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @details 按 l/m/p/q 顺序求和，每轮重算倾角函数及 Hansen 系数；零频及选定共振慢项不扣除，其余小分母项按阈值跳过并警告。
 */
SixElements Compute_KaulaShortPeriod(bastion::Overwatch& tracer, const AveragingOptions& options,
                              const SixElements& keplerian, double greenwich_angle, int resonance)
{
    // Kepler 六槽为 [a,e,i,Omega,omega,M]；a 已归一化，其余角度均为 rad。
    const double a = keplerian[0];
    const double e = keplerian[1];
    const double i = keplerian[2];

    // 三个轨道角依次控制节点、近地点及卫星沿轨道的位置。
    const double Omega = keplerian[3];
    const double omega = keplerian[4];
    const double M = keplerian[5];

    // 倾角和半倾角因子用于各谐项，以及末尾从组合改正到 NSG2 的变换。
    const double si = sin(i);
    const double ci = cos(i);

    // 半倾角因子用于将节点改正组合成 h/k，输入定义域排除了 i=pi。
    const double cih = cos(i / 2);
    const double tih = tan(i / 2);

    // 低阶 Hansen 截断时取 xe=1；否则 xe=sqrt(1-e²)，je=xe²。
    double xe = 1;

    if (options.hansen_order >= 2)
    {
        xe = sqrt(1 - e * e);
    }

    // je 跟随当前截断分支的 xe 定义，不在低阶截断时另用 1-e²。
    const double je = xe * xe;

    // 当前均根数估计给出 n；frequency=(k*n-m*omega_E)/n 为无量纲小分母。
    const double alpha = RotationRate * a * sqrt(a);
    SixElements summed_terms{};
    bool skipped_small_divisor = false;

    // l 是重力展开次数，m 是阶数；m=0 为带谐，m>0 为田谐。
    // 严格保留 l、m、p、q 四重循环的求和次序，避免改变浮点累加过程。
    for (int l = 3; l <= options.gravity_degree; ++l)
    {
        for (int m = 0; m <= l; ++m)
        {
            // 对应开关关闭时略过整组系数，J2 和 J22 由主流程单独计算。
            if ((m == 0 && !options.use_jl) || (m > 0 && !options.use_jlm))
            {
                continue;
            }

            // J 是未归一化幅值；l-m 的奇偶性决定后面采用哪组正弦、余弦。
            double phase;
            const double J = GravityAmplitudeAndPhase(l, m, phase);
            const int parity = (l - m) % 2;

            // p 是倾角展开索引，取 0..l；F、dF 分别为 F_lmp(i) 及其一阶倾角导数。
            for (int p = 0; p <= l; ++p)
            {
                // 每个 p 都单独求倾角函数及导数，避免复用上一项的幅值。
                double F;
                double dF;
                double ddF;

                // 接口按当前倾角一次给出 F 及两阶导数，本处后续只用前两项。
                FlmpDerivatives(l, m, p, i, F, dF, ddF);

                // 接口同时提供二阶导数，但本处只将 F 和 dF 乘入当前重力幅值。
                F *= J;
                dF *= J;
                const int r = l - 2 * p;

                // q 是偏心率展开的谐波索引；k=l-2p+q 是平近点角 M 的整数倍频。
                for (int q = -options.hansen_order; q <= options.hansen_order; ++q)
                {
                    const int k = r + q;

                    // m=0、k=0 同时成立时是非短周期的零频项。
                    if (m == 0 && k == 0)
                    {
                        continue;
                    }

                    // 命中已选共振的项应保留在慢变量中，不能作为短周期扣除。
                    if (m > 0 && resonance >= 1 &&
                        m * ResonanceNumerator[resonance] == k * ResonanceDenominator[resonance])
                    {
                        continue;
                    }

                    // 每项使用当前估计的 alpha；同步共振容差与其他共振仍分别为 .005、.002。
                    const double frequency = k - m * alpha;
                    double threshold = .002;

                    // 同步共振的频率比容差大于其他类别，仍保留原取值。
                    if (resonance == 2)
                    {
                        threshold = .005;
                    }

                    // 零分母保护始终执行；已识别共振附近再按次数 l 放宽同一阈值。
                    if (m > 0 && (std::abs(frequency) < 1e-14 ||
                        (resonance >= 0 && std::abs(frequency) <= l * threshold)))
                    {
                        // 记下这一类跳项，求和结束后仅上报一次原警告。
                        skipped_small_divisor = true;
                        continue;
                    }

                    // den 是无量纲频率倒数，scale 再乘 a^(-l)。
                    const double den = 1 / frequency;
                    const double scale = den / pow(a, l);

                    // G 与 dG 由下方 Hansen 接口按这一组 l、p、q 索引求得。
                    double G;
                    double dG;

                    // G 为指定 Hansen 系数，dG 是对偏心率的导数；每轮按当前 e 重算。
                    HansenSeries(-l - 1, r, k, e, options.hansen_order, G, dG);
                    SixElements harmonic{};
                    harmonic[0] = 2 * k * F * G * scale * a;

                    // 近圆轨道改用 dG，避免在这一支中直接计算 G/e。
                    double eccentricity_factor;

                    if (e <= 1e-6)
                    {
                        eccentricity_factor = dG;
                    }
                    else
                    {
                        // e>1e-6 时直接使用 G/e。
                        eccentricity_factor = G / e;
                    }

                    // 先计算偏心率组合中的相消项，再乘倾角、频率及半长轴幅值。
                    harmonic[1] = (q * je * eccentricity_factor - e * xe * r * G / (1 + xe)) * F * scale;

                    // 小倾角使用两支极限式；一般倾角才直接除以 sin(i)。
                    double inclination_factor;

                    if (std::abs(si) <= 1e-6)
                    {
                        if (r == m)
                        {
                            // r=m 时使用这一支小倾角极限。
                            inclination_factor = -si / (1 + ci) * m * F;
                        }
                        else
                        {
                            // 其他 r/m 组合使用倾角导数给出的极限式。
                            inclination_factor = (r * ci - m) * dF;
                        }
                    }
                    else
                    {
                        // 一般倾角才除以 sin(i)，不把这条式子用于小倾角分支。
                        inclination_factor = (r * ci - m) * F / si;
                    }

                    // 倾角几何确定后，再乘偏心率函数及同一谐项的幅值。
                    harmonic[2] = inclination_factor / xe * G * scale;

                    // 后三槽暂存节点、近地点和平经度的组合中间项，不是普通 Kepler 角改正。
                    harmonic[3] = dF * G / xe * scale;
                    harmonic[4] = xe * F * dG * scale;
                    harmonic[5] = (2 + 2 * l - 3 * k * den) * F * G * scale;

                    // Kaula 相角包含 kM、(l-2p)omega、m(Omega-thetaG) 及重力系数相位。
                    const double angle = k * M + r * omega + m * (Omega - greenwich_angle) - phase;
                    double S;
                    double SS;

                    // l-m 为奇数时用 sin/-cos，为偶数时用 cos/sin。
                    if (parity != 0)
                    {
                        S = sin(angle);
                        SS = -cos(angle);
                    }
                    else
                    {
                        // 偶数奇偶类使用余弦和正弦，顺序与对应公式一致。
                        S = cos(angle);
                        SS = sin(angle);
                    }

                    // 前三槽和后三槽分别乘对应相角因子；仍按 0/3、1/4、2/5 顺序累计。
                    for (int slot = 0; slot < 3; ++slot)
                    {
                        summed_terms[slot] += harmonic[slot] * S;
                        summed_terms[slot + 3] += harmonic[slot + 3] * SS;
                    }
                }
            }
        }
    }

    // 未被预定义共振精确覆盖的小分母项只记录警告，不改变上述跳过策略。
    if (skipped_small_divisor)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__,
            "Unhandled tesseral small divisors were excluded from short-period terms.",
            bastion::WARNING_UNHANDLED_TESSERAL_ORBITAL_RESONANCE);
    }

    // NSG2 偏心率分量沿 omega+Omega 定向；sw 先组合 e*(delta_omega+delta_Omega)。
    const double swo = sin(omega + Omega);
    const double cwo = cos(omega + Omega);
    const double sw = e * tih * summed_terms[3] + summed_terms[4];

    // 输出六槽依次为 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]。
    // xi=e*cos(omega+Omega)、eta=e*sin(omega+Omega)，这里直接组合以避免再除以 e。
    SixElements nsg2_short_period{};
    nsg2_short_period[0] = summed_terms[0];

    // 偏心率两槽由幅值变化和方向转动一起构成。
    nsg2_short_period[1] = cwo * summed_terms[1] - swo * sw;
    nsg2_short_period[2] = swo * summed_terms[1] + cwo * sw;

    // h=sin(i/2)*cos(Omega)、k=sin(i/2)*sin(Omega)，共同吸收倾角和节点改正。
    nsg2_short_period[3] = .5 * cih * cos(Omega) * summed_terms[2] - sin(Omega) * summed_terms[3] / (2 * cih);
    nsg2_short_period[4] = .5 * cih * sin(Omega) * summed_terms[2] + cos(Omega) * summed_terms[3] / (2 * cih);

    // lambda=M+omega+Omega；末槽将三个中间项按无奇异组合相加。
    nsg2_short_period[5] = tih * summed_terms[3] + e / (1 + xe) * summed_terms[4] + summed_terms[5];
    return nsg2_short_period;
}

/**
 * @brief 判定当前轨道是否存在圆柱地影，并准备太阳投影 A/B。
 * @param[in,out] tracer 记录轨道交地警告及地影边界不收敛错误。
 * @param[in] sun 已按本轮卫星 Omega、omega 组合的太阳方向三角量。
 * @param[in] a 当前卫星半长轴，以 ae 归一化。
 * @param[in] e 当前卫星偏心率，无量纲。
 * @param[in] si sin(i)，无量纲。
 * @param[in] ci cos(i)，无量纲。
 * @param[out] A 写入卫星轨道面与太阳方向的余弦投影组合，无量纲。
 * @param[out] B 写入卫星轨道面与太阳方向的正弦投影组合，无量纲。
 * @return 判为有影时返回 true，否则 false；出错时返回 false 并设置 tracer，调用方须检查错误码。
 * @details 最多进行 100 次边界更新，用 1 rad 的迭代判据判断有无地影；投影退化或轨道交地时判无影。外层仅在首轮调用，并固定本次变换的地影结论及 A/B。
 */
bool Check_OrbitEclipse(bastion::Overwatch& tracer, const ThirdBodyArguments& sun,
                      double a, double e, double si, double ci, double& A, double& B)
{
    // 太阳组合三角量已按当前卫星 Omega、omega 更新，下面分别计算 A/B 投影。
    A = .25 * ((1 - ci) * (sun.cip1 * sun.c_mp + sun.cim1 * sun.c_mm) +
        (1 + ci) * (sun.cip1 * sun.c_pm + sun.cim1 * sun.c_pp) + 2 * si * sun.si * (sun.c_m - sun.c_p));

    // B 使用同组相角的正弦分量，整体负号来自参考地影方向定义。
    B = -.25 * ((1 - ci) * (sun.cip1 * sun.s_mp + sun.cim1 * sun.s_mm) +
        (1 + ci) * (sun.cip1 * sun.s_pm + sun.cim1 * sun.s_pp) + 2 * si * sun.si * (sun.s_m - sun.s_p));

    // K 和方向相角确定太阳在轨道面的投影，rp2 是地球半径与半通径之比的平方。
    const double K = hypot(A, B);
    const double angle = atan2(A, B);
    const double rp2 = pow(constants::Wgs84::EQUATORIAL_RADIUS / UNIT_LENGTH_M / (a * (1 - e * e)), 2);

    // 投影退化时判无影；半通径落入地球时警告并判无影。
    if (K < 1e-15)
    {
        return false;
    }

    if (rp2 > 1)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__,
            "Mean orbit intersects Earth; eclipse is deemed absent.", bastion::WARNING);
        return false;
    }

    // 先用半通径给出两端真近点角初值；sq<=-1 时两端从同一初值开始。
    const double sq = -sqrt(1 - rp2) / K;
    double f1;
    double f2;

    if (sq <= -1)
    {
        f1 = 1.5 * PI - angle;
        f2 = f1;
    }
    else
    {
        // 两端分别取反正弦对应的进入和离开方向初值。
        f1 = PI - asin(sq) - angle;
        f2 = 2 * PI + asin(sq) - angle;
    }

    // 最多迭代 100 次，依当前边界处半径更新圆柱地影方程。
    for (int iteration = 0; iteration < 100; ++iteration)
    {
        const double x1 = rp2 * pow(1 + e * cos(f1), 2);
        const double x2 = rp2 * pow(1 + e * cos(f2), 2);

        // 任一边界处落入地球时警告并判无影。
        if (x1 > 1 || x2 > 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "Mean orbit intersects Earth; eclipse is deemed absent.", bastion::WARNING);
            return false;
        }

        // 两个边界各自修正；超出反正弦定义域时返回无影。
        const double s1 = -sqrt(1 - x1) / K;
        const double s2 = -sqrt(1 - x2) / K;

        if (s1 <= -1 || s2 <= -1)
        {
            return false;
        }

        // 两端变化均不超过 1 rad 时判有影；本函数只判断有无地影，边界角不作为输出。
        const double next1 = PI - asin(s1) - angle;
        const double next2 = 2 * PI + asin(s2) - angle;

        if (std::abs(next1 - f1) <= 1 && std::abs(next2 - f2) <= 1)
        {
            return true;
        }

        // 尚未满足该判据时，以本轮两端位置继续下一次边界迭代。
        f1 = next1;
        f2 = next2;
    }

    // 100 次仍未满足地影判据时，记录错误供外层停止本次变换。
    tracer.LogEvent(__FILE__, __func__, __LINE__, "Eclipse-boundary iteration failed to converge", bastion::ERROR);
    return false;
}
} // namespace

/**
 * @brief 在给定摄动模型下，将 GCRS 密切 NSG2 根数转换为均根数。
 * @param[in,out] tracer 记录输入定义域、摄动模型警告及公式计算错误，失败时保留调用链。
 * @param[in] mjd_utc 轨道历元，UTC 的 MJD；岁差章动转 TT，日月平位置和格林尼治角取 UT1=UTC。
 * @param[in] osculating GCRS 密切 NSG2，按 [a,xi,eta,h,k,lambda] 排列；半长轴尺度 ae，平经度 rad，其余无量纲。
 * @param[in] options 选择 J2/J22、坐标架、JL/JLM、日月及光压项，指定固定迭代次数、展开阶数和物理 AMR(m^2/kg)、Cr。
 * @param[out] mean 成功时覆写同历元、同坐标系和同归一化尺度的均 NSG2 根数；失败时保持传入值。
 * @details 固定次数求解 z_mean=z_osc-short(z_mean)，每轮更新轨道几何及倾角/Hansen 函数。
 * xi=e*cos(omega+Omega)、eta=e*sin(omega+Omega)、h=sin(i/2)*cos(Omega)、k=sin(i/2)*sin(Omega)、lambda=M+omega+Omega。
 * 时间尺度 UNITTIME=sqrt(ae^3/GE)。输入要求椭圆、i<pi、a<=90000 km；启用 JLM 且输入 e>=0.5 时，本次调用关闭 JLM 并将最高次数设为 4。
 * 日月位置复用 Cart_MeanSun/Moon，每条记录独立准备；地影结果及其 A/B 投影在本次首轮固定，有影时省略光压短周期项。
 * J2 二阶半长轴和平经度独立扣除，J22 高阶半长轴替换对应槽，月球双快变量六槽替代单快变量结果。
 */
void OsculatingToMean(bastion::Overwatch& tracer, double mjd_utc,
    const NonSingular2& osculating, const AveragingOptions& options, NonSingular2& mean)
{
    // 调用链已有错误时立即返回；后续失败仍通过同一 tracer 上报。
    if (tracer.Code() != bastion::ERR_NOERROR)
    {
        return;
    }

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    try
    {
        // ==================== 01 检查输入与模型定义域 ====================
        // 原始密切六槽固定为 [a,xi,eta,h,k,lambda]，迭代中始终从这组输入扣除短周期项。
        std::array<double, 6> original;
        osculating.ToArray(original.data());

        for (double value : original)
        {
            if (!std::isfinite(value))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Nonfinite osculating state", bastion::ERROR);
                return;
            }
        }

        // 历元须有限，重力次数、固定迭代次数及 Hansen 展开阶数在各自支持范围内检查。
        // AMR 是物理 m²/kg，Cr 是无量纲光压系数；它们的有限性和范围在此统一检查。
        if (!std::isfinite(mjd_utc) || options.gravity_degree < 2 || options.gravity_degree > 10 ||
            options.iteration_count < 1 || options.iteration_count > 100 ||
            options.hansen_order < 0 || options.hansen_order > 4 ||
            !std::isfinite(options.srp_area_mass_m2_kg) || options.srp_area_mass_m2_kg < 0 ||
            !std::isfinite(options.radiation_coefficient) ||
            options.radiation_coefficient < 0 || options.radiation_coefficient > 2)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid averaging options or epoch", bastion::ERROR);
            return;
        }

        // hypot(xi,eta)=e、hypot(h,k)=sin(i/2)，因此这里只接受椭圆及 i<pi。
        if (original[0] <= 0 || hypot(original[1], original[2]) >= 1 ||
            hypot(original[3], original[4]) >= 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "O2M requires an elliptic NSG2 state with inclination below pi", bastion::ERROR);
            return;
        }

        // 半长轴由 ae 归一化量恢复为 km 后检查；超过硬上限时不进入解析展开。
        if (original[0] * UNIT_LENGTH_M / 1000 > 90000)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Semimajor axis exceeds 90000 km O2M limit", bastion::ERROR);
            return;
        }

        if (original[0] * UNIT_LENGTH_M / 1000 > 43000)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "Semimajor axis exceeds 43000 km; short-period accuracy may deteriorate.", bastion::WARNING);
        }

        // 只在 JLM 已开启时应用参考算法的大偏心率临时处理；复制配置以免修改调用者选项。
        AveragingOptions active_options = options;

        if (active_options.use_jlm && hypot(original[1], original[2]) >= .5)
        {
            active_options.use_jlm = false;
            active_options.gravity_degree = 4;
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "CRITICAL eccentricity >= 0.5: JLM disabled; gravity degree set to 4.",
                bastion::WARNING_AVERAGING_OVER_CRITICAL_ECCENTRICITY);
        }
        else if (active_options.use_jlm && hypot(original[1], original[2]) >= .25)
        {
            // 中等偏心率继续计算 JLM，并报告 Kaula 展开精度警告。
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "Eccentricity >= 0.25; Kaula expansion accuracy deteriorates.",
                bastion::WARNING_AVERAGING_OVER_MILD_ECCENTRICITY);
        }

        // 岁差章动短周期项与 J2 耦合；仅开启 frame 而关闭 J2 时忽略该项并警告。
        if (active_options.use_frame && !active_options.use_j2)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "Frame-induced perturbation requires J2 and is omitted.", bastion::WARNING);
        }

        // ==================== 02 准备时标、日月与重力参数 ====================
        // 本模型取 DUT1=0，即 UT1=UTC；日月平位置函数接收 UT1，不将 TT 代入其参数。
        // TT 用于下面的简化岁差章动。
        const double mjd_tt = UTC2TT(mjd_utc);
        const double ut1_centuries = (mjd_utc - constants::J2000_MJD) / constants::DAYS_PER_CENTURY;
        const double tt_centuries = (mjd_tt - constants::J2000_MJD) / constants::DAYS_PER_CENTURY;

        // 格林尼治角由 UT1 世纪数给出，回绕后为 [0,2pi) rad。
        const double greenwich_angle = Wrap_Angle(
            4.894961223731367 + GmstRatePerCentury * ut1_centuries +
            .000006770713945 * ut1_centuries * ut1_centuries -
            .000450876723432e-6 * ut1_centuries * ut1_centuries * ut1_centuries);

        // 简化章动仅保留月球升交点主项，月球节点经度先由度转换为 rad。
        const double lunar_node_longitude = (125.0446 - 1934.14 * tt_centuries) * DEGREE_TO_RAD;
        const FrameArguments frame{
            // thetaA 保留 TT 世纪的一次岁差项，Dtheta、Depsilon 是简化章动主项，均为 rad。
            2004.3109 * tt_centuries * ARCSEC_TO_RAD,
            -6.84 * sin(lunar_node_longitude) * ARCSEC_TO_RAD,
            9.20 * cos(lunar_node_longitude) * ARCSEC_TO_RAD
        };

        // 月球只在月球摄动开启时准备；太阳位置同时服务太阳摄动和光压模型。
        ThirdBodyArguments moon{};
        ThirdBodyArguments sun{};

        if (active_options.use_lunar)
        {
            moon = Prepare_ThirdBody(tracer, mjd_utc, true);
            if (tracer.Code() != bastion::ERR_NOERROR) return;
        }

        // 即使关闭太阳引力，只要启用光压仍需太阳轨道。
        if (active_options.use_solar || active_options.use_srp)
        {
            sun = Prepare_ThirdBody(tracer, mjd_utc, false);
            if (tracer.Code() != bastion::ERR_NOERROR) return;
        }

        // 解析轨道失败时停止本次计算，调用者的 mean 保持不变。
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            return;
        }

        double phase;
        double j22_coefficient_phase;

        // J2 取未归一化 C20 的负值；J22 将 C22/S22 合为幅值和相位。
        const double j2 = -GravityAmplitudeAndPhase(2, 0, phase);
        const double j22 = GravityAmplitudeAndPhase(2, 2, j22_coefficient_phase);

        // J22 相位折半后与格林尼治角合成参考经度，公式中使用 Omega-此经度。
        const double j22_reference_longitude = greenwich_angle + j22_coefficient_phase / 2;

        // 共振类别由输入密切半长轴固定，current 则从完整输入状态开始迭代。
        const int resonance = Find_Resonance(original[0]);
        NonSingular2 current = osculating;

        // 地影几何在本次瞬转平的第一轮固定；这些变量不跨输入记录复用。
        double eclipse_A = 0;
        double eclipse_B = 0;
        bool has_eclipse = false;

        // ==================== 03 按固定次数更新均根数 ====================
        // 每轮计算 short(current)，再从同一 original 扣除；不根据残差提前终止。
        for (int iteration = 0; iteration < active_options.iteration_count; ++iteration)
        {
            // -------------------- 当前估计的轨道几何 --------------------
            Keplerian keplerian;

            // 从当前均根数估计恢复 Kepler 几何；失败不使用本轮部分数据。
            current.ToKeplerian(tracer, keplerian);
            if (tracer.Code() != bastion::ERR_NOERROR)
            {
                return;
            }

            // 同时保留 Kepler 六槽与当前 NSG2 六槽，供不同解析接口读取。
            std::array<double, 6> elements;
            keplerian.ToArray(elements.data());
            std::array<double, 6> nsg2;
            current.ToArray(nsg2.data());

            // a 为归一化半长轴，e 为偏心率，i 为轨道倾角。
            const double a = elements[0];
            const double e = elements[1];
            const double i = elements[2];

            // Omega、omega、M 分别为升交点赤经、近地点幅角和平近点角，单位 rad。
            const double Omega = elements[3];
            const double omega = elements[4];
            const double M = elements[5];

            // 由当前平近点角解偏近点角 E，再转为真近点角 f，供各短周期公式使用。
            const double E = MeanAnomaly2EccentricAnomaly(M, e);
            const double f = EccentricAnomaly2TrueAnomaly(E, e);

            // je=1-e²、xe=sqrt(je)；p=a(1-e²) 是归一化半通径。
            const double je = 1 - e * e;
            const double xe = sqrt(je);
            const double p = a * je;

            // ar=a/r；n=a^(-3/2) 为归一化平运动，后续与地球自转频率组成比值。
            const double ar = 1 / (1 - e * cos(E));
            const double n = 1 / (a * sqrt(a));

            // 倾角及半倾角三角量进入 J2、J22、frame 和第三体的 NSG2 组合式。
            const double si = sin(i);
            const double ci = cos(i);

            // NSG2 的 h/k 使用半倾角，因此同时保存 sin(i/2) 和 cos(i/2)。
            const double sih = sin(i / 2);
            const double cih = cos(i / 2);

            // Omega 的方向余弦用于节点投影，omega+Omega 用于偏心率向量投影。
            const double so = sin(Omega);
            const double co = cos(Omega);

            // 偏心率向量以近地点经度 omega+Omega 为方向。
            const double swo = sin(omega + Omega);
            const double cwo = cos(omega + Omega);

            // fm 是回绕到 [-pi,pi) 的 f-M；fmesf 再加 e*sin(f)，保持原公式定义。
            const double fm = Wrap_Angle(f - M + PI) - PI;
            const double fmesf = fm + e * sin(f);

            // apow[j]=a^j、epow[j]=e^j，各自从零次幂 1 递推，供第三体及光压公式读取。
            std::vector<double> apow(11, 1);
            std::vector<double> epow(6, 1);

            // 逐次相乘生成 a¹ 到 a¹⁰。
            for (std::size_t j = 1; j < apow.size(); ++j)
            {
                apow[j] = apow[j - 1] * a;
            }

            // 偏心率幂单独递推，保持先准备 a 幂、再准备 e 幂的顺序。
            for (std::size_t j = 1; j < epow.size(); ++j)
            {
                epow[j] = epow[j - 1] * e;
            }

            // 三组倍角表按下标 j 存 j*M、j*f、j*E；零倍正弦为 0，零倍余弦为 1。
            double sm[5]{};
            double cm[5]{1};

            // f 的倍角按相同下标单独保存，避免与 M 或 E 的表混用。
            double sf[5]{};
            double cf[5]{1};

            // E 的倍角同样从零次项开始；公式接口直接读取数组。
            double se[5]{};
            double ce[5]{1};

            // 每轮由当前根数重算倍角表，不跨迭代或跨输入记录保存几何量。
            for (int j = 1; j < 5; ++j)
            {
                sm[j] = sin(j * M);
                cm[j] = cos(j * M);

                // 真近点角倍角供引力短周期项使用。
                sf[j] = sin(j * f);
                cf[j] = cos(j * f);

                // 偏近点角倍角供第三体及光压短周期项使用。
                se[j] = sin(j * E);
                ce[j] = cos(j * E);
            }

            // -------------------- 按模型顺序累计短周期项 --------------------
            // total 六槽为 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda]。
            SixElements total{};

            // J2 二阶半长轴和平经度单独保存，不并入 total，最后各扣除一次。
            double second_order_a = 0;
            double second_order_lambda = 0;

            // J2 一阶六槽与二阶半长轴共用当前轨道几何，nsg2[3]/[4] 分别为 h/k。
            if (active_options.use_j2)
            {
                SixElements term{};

                // a/e 与近圆因子随后连接倾角、纬度幅角和径向几何，均为归一化量。
                J2Short(tracer, nsg2[3], nsg2[4], 1.5 * j2 / (p * p), 1.5 * j2 / a,
                    a, e, 1 / (1 + xe), (1 + 2 * xe) / pow(1 + xe, 2),
                    si, ci, omega, omega + Omega, f, ar, fmesf, term.data(), second_order_a);
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 二阶平经度另算 delta_lambda；此调用不覆写上面的一阶六槽。
                J2SecondOrderLongitude(tracer, j2 * j2 / pow(p, 4), ci, si * si, pow(si, 4),
                    sin(2 * f + 2 * omega), sin(4 * f + 4 * omega), second_order_lambda);
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 按 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda] 累计本项。
                for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
            }

            // J22 首先计算 NSG2 六槽，alpha=RotationRate/n；共振类别仍由输入密切轨道确定。
            if (active_options.use_j22)
            {
                SixElements term{};

                // 节点相对 J22 参考经度定向，swo/cwo 则用于偏心率分量的方向组合。
                J22Short(tracer, nsg2[3], nsg2[4], j22 / (a * a), a, e, si, ci, cih, so, co,
                    Omega - j22_reference_longitude, swo, cwo, omega, M,
                    RotationRate / n, term.data(), resonance);
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 较高阶半长轴公式只覆写 term[0]，其余五槽继续保留 J22Short 的输出。
                J22SemimajorShort(tracer, j22 / a, e, ci, Omega - j22_reference_longitude,
                    omega, M, RotationRate / n, term[0], resonance);
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 按 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda] 累计本项。
                for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
            }

            // 简化岁差章动与 J2 的耦合项仅在两个开关同时开启时加入。
            if (active_options.use_j2 && active_options.use_frame)
            {
                SixElements term{};

                // xe 与 1/(xe*je) 分别是 sqrt(1-e²)、(1-e²)^(-3/2)。
                FrameShort(tracer, frame, 3 * j2 / a, 1.5 * j2 / (p * p),
                    e, xe, 1 / (xe * je), si, ci, sih, cih, so, co, sin(2 * omega), cos(2 * omega),
                    swo, cwo, sf[1], cf[1], fm, fmesf, pow(ar, 3), term.data());
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 按 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda] 累计本项。
                for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
            }

            // l>=3 的带谐/田谐由 Kaula 四重求和给出，在当前估计上重算倾角与 Hansen 函数。
            if (active_options.use_jl || active_options.use_jlm)
            {
                const SixElements term = Compute_KaulaShortPeriod(tracer, active_options, elements, greenwich_angle, resonance);
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 按 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda] 累计本项。
                for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
            }

            // 月球方向先与当前卫星节点及近地点组合，再计算月球引力短周期项。
            if (active_options.use_lunar)
            {
                Update_ThirdBodyAngles(moon, Omega, omega);
                SixElements term{};
                ThirdBodyProducts products{};

                // 单快变量项同时准备 products；apow 为 a 的幂，1/ar=r/a。
                ThirdBodyShort(tracer, apow, 1 / ar, e, je, si, ci, sih, cih, so, co, swo, cwo,
                    sm, cm, sf, cf, se, ce, moon, products, term.data());
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 双快变量公式需要上面已准备的 AA/BB/AB，并覆写六槽来替代单快变量输出。
                if (active_options.lunar_double_average)
                {
                    ThirdBodyDoubleShort(tracer, apow, e, je, si, ci, sih, cih, so, co, swo, cwo,
                        sm, cm, moon, products, term.data());
                    if (tracer.Code() != bastion::ERR_NOERROR) return;
                }

                // 无论采用哪一支，月球 NSG2 短周期只在此累计一次。
                for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
            }

            // 太阳方向同时供太阳引力与光压使用，任一模型开启时都在本轮更新。
            if (active_options.use_solar || active_options.use_srp)
            {
                Update_ThirdBodyAngles(sun, Omega, omega);
            }

            // 太阳引力与太阳光压在此仍是独立开关，不能把后者当作前者的一部分。
            if (active_options.use_solar)
            {
                SixElements term{};
                ThirdBodyProducts products{};

                // 太阳引力使用单快变量公式，输入几何和六槽输出约定与月球相同。
                ThirdBodyShort(tracer, apow, 1 / ar, e, je, si, ci, sih, cih, so, co, swo, cwo,
                    sm, cm, sf, cf, se, ce, sun, products, term.data());
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                // 按 [delta_a,delta_xi,delta_eta,delta_h,delta_k,delta_lambda] 累计本项。
                for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
            }

            // 光压地影只在第一轮判定并固定；存在地影时不加入这里的解析短周期项。
            if (active_options.use_srp)
            {
                if (iteration == 0)
                {
                    has_eclipse = Check_OrbitEclipse(tracer, sun, a, e, si, ci, eclipse_A, eclipse_B);
                    if (tracer.Code() != bastion::ERR_NOERROR) return;
                }

                // 地影结果在本次调用内固定，有影分支不累加此套解析光压项。
                if (!has_eclipse)
                {
                    SixElements term{};

                    // 压力乘物理面质比得到 m/s²，再除以程序加速度尺度。
                    // beta 的负号表示光压方向与地心指向太阳的方向相反。
                    const double beta = -active_options.radiation_coefficient *
                        active_options.srp_area_mass_m2_kg * constants::SolarRadiation::REFERENCE_PRESSURE /
                        (UNIT_GM_M3_S2 / UNIT_LENGTH_M / UNIT_LENGTH_M) *
                        pow(constants::SolarRadiation::REFERENCE_DISTANCE / UNIT_LENGTH_M / sun.r, 2);

                    // Combined 分支直接输出 NSG2 组合改正，使用第一轮得到的 eclipse_A/B。
                    RadiationShort(tracer, beta, apow, epow, xe, si, ci, cih, so, co, swo, cwo,
                        se, ce, eclipse_A, eclipse_B, sun, term.data(), RadiationCombined);
                    if (tracer.Code() != bastion::ERR_NOERROR) return;

                    // 按 NSG2 分量累计本轮光压短周期项。
                    for (int slot = 0; slot < 6; ++slot) total[slot] += term[slot];
                }
            }

            // -------------------- 扣除短周期并验收本轮状态 --------------------
            // 更新式为 z_mean=z_osc-short(z_mean)，被减数始终是固定的密切输入 original。
            // total 则由当前均根数估计得到，各摄动项已按 NSG2 六槽累计。
            SixElements next{};

            for (int slot = 0; slot < 6; ++slot)
            {
                next[slot] = original[slot] - total[slot];
            }

            // J2 两个二阶项在这里各扣一次；最终平经度回绕至 [0,2pi)。
            next[0] -= second_order_a;
            next[5] = Wrap_Angle(next[5] - second_order_lambda);

            // 先检查有限性，再检查椭圆与倾角定义域，避免把失效状态交给下一轮。
            for (double value : next)
            {
                if (!std::isfinite(value))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Nonfinite short-period correction", bastion::ERROR);
                    return;
                }
            }

            // 继续要求 e<1、sin(i/2)<1，与输入采用相同 NSG2 定义域。
            if (next[0] <= 0 || hypot(next[1], next[2]) >= 1 || hypot(next[3], next[4]) >= 1)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "O2M iteration left the elliptic NSG2 domain", bastion::ERROR);
                return;
            }

            // Set 的最后一项 true 表示仍采用归一化状态；先写内部估计，不提前覆盖 mean。
            current.Set(tracer, next[0], next[1], next[2], next[3], next[4], next[5], true);
            if (tracer.Code() != bastion::ERR_NOERROR) return;
        }

        // 各次公式调用及状态更新均已通过检查，固定次数完成后写回最终均根数。
        mean = current;
    }
    catch (const std::exception& error)
    {
        // 将标准库异常写入 tracer，失败时调用者的 mean 保持原值。
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), bastion::ERROR);
        return;
    }

    // 成功完成变换后退出当前调用层。
    tracer.CheckOut(__func__);
}

/**
 * @brief 显示瞬转平请求的摄动模型、展开设置和面积质量比。
 * @param[in,out] tracer 记录摘要写入失败；本方法只显示设置，不进行轨道计算。
 * @param[in] srp_amr_norm 与本对象 srp_area_mass_m2_kg 对应的原始归一化光压面积质量比，用于同时显示原始值和 m^2/kg 值。
 * @param[in] drag_amr_norm 配置中的归一化阻力面积质量比，只记录在日志中，不参与瞬转平计算。
 * @param[in,out] output 接收模型摘要的输出流，追加文本并保留调用前的流格式。
 * @param[in] precision 面积质量比及反射系数的有效数字数，默认 16。
 * @details 按重力次数、迭代次数、Hansen 阶数、八个模型开关、日月模型、月球模式、AMR 和 Cr 的顺序显示。
 * 路径、记录数和历元/坐标系说明由主程序写出；面积质量比的归一化因子为 MASS_OF_EARTH/ae^2。
 */
void AveragingOptions::Display(bastion::Overwatch& tracer, double srp_amr_norm, double drag_amr_norm, std::ostream& output, int precision) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 用局部文本流组织摘要，数值格式只作用于本段输出。
    std::ostringstream text;

    // 先显示展开次数与固定迭代设置，再按八个开关的固定顺序列出摄动项。
    text << "Requested model: GGM03C degree " << gravity_degree << ", iterations " << iteration_count
        << ", Hansen order " << hansen_order << "\n"
        << "Flags J2/J22/frame/JL/JLM/lunar/solar/SRP: ";

    text << use_j2 << '/' << use_j22 << '/' << use_frame << '/' << use_jl << '/'
        << use_jlm << '/' << use_lunar << '/' << use_solar << '/' << use_srp << '\n'
        << "Ephemerides lunar/solar: analytic_mean/analytic_mean\n";

    // 月球模式及光压 AMR 同时保留归一化值和物理 m²/kg 值。
    text << "Lunar short-period mode: " << (lunar_double_average ? "double_fastarg" : "single_fastarg")
        << "\nSRP_AMR_normalized: " << std::setprecision(precision) << srp_amr_norm
        << "\nSRP_AMR_m2_kg: " << srp_amr_norm * NORMALIZED_AM_TO_IAM;

    // 阻力 AMR 仅作配置留存，Cr 为参与光压计算的无量纲系数。
    text << "\nDrag_AMR_normalized (not used by O2M): " << drag_amr_norm
        << "\nReflectivity (Cr = 1 + reflection coefficient): " << radiation_coefficient << "\n\n";

    output << text.str();

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed to write O2M model settings.", bastion::ERROR);
        return;
    }

    tracer.CheckOut(__func__);
}
} // namespace orbit
