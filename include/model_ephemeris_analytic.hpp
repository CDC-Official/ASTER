//
//  ephemeris.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/28.
//

#ifndef CDC_MODEL_EPHEMERIS_ANALYTIC_HPP
#define CDC_MODEL_EPHEMERIS_ANALYTIC_HPP

#include "utils_constants.hpp"
#include "orbit_basic.hpp"
#include <iostream>
#include <cmath>

// 用于解析日月位置计算，以及解析 STM 中的日月平均摄动项。
// 参数取自汤靖师《人造卫星精密定轨(II)》2025 第 104 页的长期进动椭圆模型。
namespace constants::MeanOrbit
{
    // 解析月球椭圆的半长轴，也作为解析 STM 月球平均项中的距离，单位 m。
    constexpr double LUNAR_SEMIMAJOR_AXIS = 384747981.0;

    // 解析太阳椭圆的半长轴，也作为解析 STM 太阳平均项中的距离，由 AU 换成 m。
    constexpr double SOLAR_SEMIMAJOR_AXIS = 1.00000102 * constants::ASTRONOMICAL_UNIT;

    // J2000 平黄赤交角，用于将解析日月位置转到赤道系，并给出解析 STM 的平均轨道面倾角，单位 rad。
    const double OBLIQUITY = 23.4393 * constants::Conversion::DEGREE_TO_RAD;
}

//用于积分器的老版本，使用了static类型的MJD0，定义在RKF78.hpp里
void MeanMoonCart(double MJD0, double RKF_t, double* R_Moon);
void MeanSunCart(double MJD0, double RKF_t, double* R_Sun);
// 新版本：UT1 的 MJD 输入，返回地心 J2000 赤道系的归一化平位置，R 为三个分量。
// V 非空时另写三个辅助速度分量；仍由原六根数在 mu=1 下转换，不是历表对历元的真实导数，
// 也不是物理日月速度，仅供在同一 mu=1 约定下恢复轨道方向；位置以 ae、速度以 ae/UNITTIME 归一化。
void Cart_MeanMoon(double MJD_UT1, double* R_Moon, double* V_Moon = nullptr);
void Cart_MeanSun(double MJD_UT1, double* R_Sun, double* V_Sun = nullptr);
//不变椭圆
void Cart_OsculatingMoon(double MJD_UT1, double* R_Moon);
void Cart_OsculatingSun(double MJD_UT1, double* R_Sun);

#endif /* CDC_MODEL_EPHEMERIS_ANALYTIC_HPP */
