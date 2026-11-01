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
