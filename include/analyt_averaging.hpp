/**
 * @file analyt_averaging.hpp
 * @brief 地球卫星瞬转平的模型选项及 GCRS NSG2 变换接口。
 * @details 通过解析短周期扣除求均根数，结果依赖所选摄动模型；长度尺度 ae、时间尺度 UNITTIME=sqrt(ae^3/GE)、角度 rad。
 */
#ifndef CDC_ANALYT_AVERAGING_HPP
#define CDC_ANALYT_AVERAGING_HPP

#include "orbit_statedef.hpp"
#include <ostream>

namespace orbit {

// 瞬转平仅删除解析短周期项；所得均根数与下面选定的摄动模型相关。
// 长度采用地球赤道半径 ae 归一化，时间采用 sqrt(ae^3/GE)，角度为 rad。
struct AveragingOptions
{
    int gravity_degree = 10; // JL/JLM 的最高次数；J2/J22 由各自开关单独控制。
    int iteration_count = 3; // 固定点更新次数，不按残差自适应终止。
    int hansen_order = 4; // Kaula 求和的 q 截断及 Hansen 偏心率展开阶数，范围 0..4。
    bool use_j2 = true;
    bool use_j22 = true;
    bool use_frame = true; // 岁差章动与 J2 的耦合短周期项；仅 use_j2 同时开启时计算。
    bool use_jl = true; // l=3..gravity_degree、m=0 的带谐短周期项。
    bool use_jlm = true; // l=3..gravity_degree、m>0 的田谐短周期项，不包含单列的 J22。
    bool use_lunar = true; // 月球第三体项；开启后才会使用 lunar_double_average。
    bool use_solar = true; // 太阳第三体项；关闭后若 use_srp 开启仍需准备太阳位置。
    bool use_srp = true;
    bool lunar_double_average = false; // 月球采用卫星/月球双快变量公式，覆写单快变量结果。
    // AMR 是物理面积质量比（m²/kg），不是归一化状态量；光压公式内部再转换单位。
    double srp_area_mass_m2_kg = 0.010000001460108767;
    double radiation_coefficient = 1.44; // Cr，吸收及反射的合成系数。

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
    void Display(bastion::Overwatch& tracer, double srp_amr_norm, double drag_amr_norm, std::ostream& output, int precision = 16) const;
};

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
    const NonSingular2& osculating, const AveragingOptions& options, NonSingular2& mean);

} // namespace orbit
#endif
