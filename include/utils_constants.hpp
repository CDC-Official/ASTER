//
//  utils_constants.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#ifndef CDC_UTILS_CONSTANTS_HPP
#define CDC_UTILS_CONSTANTS_HPP

#include <cmath>
#include <sofam.h>

// TDB 相容的天体 GM 采用 Park 等（2021，AJ 161:105）DE440/DE441 论文表 2，单位 m³/s²。
// 顺序按 JPL 编号 1—11：水星、金星、地球、火星、木星、土星、天王星、海王星、冥王星、月球、太阳。
// 火星至冥王星采用表中的行星系统 GM；这些参数也是李昊翰程序注释所引用的来源。
constexpr double MAJOR_BODY_GM[11] = {
    2.2031868551e13,
    3.24858592000e14,
    3.98600435507e14,
    4.2828375816e13,
    1.26712764100000e17,
    3.7940584841800e16,
    5.794556400000e15,
    6.836527100580e15,
    9.75500000e11,
    4.902800118e12,
    1.32712440041279419e20
};

// IERS Conventions 2010 第 6.1 节的 TT 相容地球 GM，单位 m³/s²，供 GCRS/TT 动力学使用。
// GE 也定义程序的数值归一化尺度；E-BCRS/TDB 动力学另用 MAJOR_BODY_GM 中的物理参数。
constexpr double GEOCENTRIC_GRAVITATIONAL_CONSTANT = 3.986004415e14;
constexpr double GE = GEOCENTRIC_GRAVITATIONAL_CONSTANT;

// TT 与 TDB 相容参数的尺度比；IAU 定义的 L_G、L_B 取自 SOFA。
// 长度和 GM 均按此比例缩放，日月 GM 转为 TT 相容值后用于 GCRS 摄动。
constexpr double TDB_TO_TT_SCALE = (1.0 - ELG) / (1.0 - ELB);
constexpr double LUNAR_GRAVITATIONAL_CONSTANT = MAJOR_BODY_GM[9] * TDB_TO_TT_SCALE;
constexpr double HELIOCENTRIC_GRAVITATIONAL_CONSTANT = MAJOR_BODY_GM[10] * TDB_TO_TT_SCALE;
const double CONSTANT_OF_GRAVITATION = 6.672E-11;                      // Constant of gravitaion G m3kg-1s-2 跟老师一致
const double MASS_OF_EARTH = GEOCENTRIC_GRAVITATIONAL_CONSTANT / CONSTANT_OF_GRAVITATION;//地球质量 kg
//const double ae = 6378136.6;                                          // Equatorial radius of the Earth 之前的应该不自洽，纯二体的误差都比较大，达不到double的精度
const double
	ae = 6378137.0;                                          // Equatorial radius of the Earth 跟老师一致，二体预报精度跟double符合

const double EARTH_FLATTENING = 1.0 / 298.25642;                              //Flattening factor of the Earth
const double E_SQUARED = EARTH_FLATTENING * (2 - EARTH_FLATTENING);                // e^2=2f-f^2
const double J2 = 0.0010826359;                                          // Dynamic form-factor

const double pi = 3.141592653589793238462643;                          // 15～16 accuracy
const double UNITTIME = sqrt(pow(ae, 3) / GE);                          // unit time 跟老师一致
const double START_OF_MJD = 2400000.5;

// 真空光速，单位 m/s；归一化时除以各模型的速度尺度。
constexpr double SPEED_OF_LIGHT_MPS = 299792458.0;

// TDB 天体 GM 除以数值归一化常数 GE，供 E-BCRS 的牛顿项和后牛顿项共同使用。
// 地球项也保留这一系数，因为 TDB 地球 GM 与 GE 的 TT 数值不同。
constexpr double MAJOR_BODY_GM_NORMALIZED[11] = {
    MAJOR_BODY_GM[0] / GE,
    MAJOR_BODY_GM[1] / GE,
    MAJOR_BODY_GM[2] / GE,
    MAJOR_BODY_GM[3] / GE,
    MAJOR_BODY_GM[4] / GE,
    MAJOR_BODY_GM[5] / GE,
    MAJOR_BODY_GM[6] / GE,
    MAJOR_BODY_GM[7] / GE,
    MAJOR_BODY_GM[8] / GE,
    MAJOR_BODY_GM[9] / GE,
    MAJOR_BODY_GM[10] / GE
};

// GCRS/TT 动力学中，日月与地球 GM 使用同一时标约定。
constexpr double MOON_EARTH_MASS_RATIO = LUNAR_GRAVITATIONAL_CONSTANT / GE;
constexpr double SUN_EARTH_MASS_RATIO = HELIOCENTRIC_GRAVITATIONAL_CONSTANT / GE;
//速度归一化
const double MPS_TO_NORMALIZED_VELOCITY = UNITTIME / ae;
const double KMPS_TO_NORMALIZED_VELOCITY = 1000.0 * UNITTIME / ae;
const double NORMALIZED_VELOCITY_TO_MPS = 1.0 / MPS_TO_NORMALIZED_VELOCITY;
const double NORMALIZED_VELOCITY_TO_KMPS = 1.0 / KMPS_TO_NORMALIZED_VELOCITY;
//长度归一化
const double M_TO_NORMALIZED_LENGTH = 1.0 / ae;
const double KM_TO_NORMALIZED_LENGTH = 1000.0 * 1.0 / ae;
const double NORMALIZED_LENGTH_TO_M = 1.0 / M_TO_NORMALIZED_LENGTH;
const double NORMALIZED_LENGTH_TO_KM = 1.0 / KM_TO_NORMALIZED_LENGTH;
//时间归一化
const double SEC_TO_NORMALIZED_TIME = 1.0 / UNITTIME;
const double MIN_TO_NORMALIZED_TIME = 60.0 / UNITTIME;
const double HOUR_TO_NORMALIZED_TIME = 3600.0 / UNITTIME;
const double DAY_TO_NORMALIZED_TIME = 86400.0 / UNITTIME;
const double NORMALIZED_TIME_TO_SEC = 1.0 / SEC_TO_NORMALIZED_TIME;
const double NORMALIZED_TIME_TO_MIN = 1.0 / MIN_TO_NORMALIZED_TIME;
const double NORMALIZED_TIME_TO_HOUR = 1.0 / HOUR_TO_NORMALIZED_TIME;
const double NORMALIZED_TIME_TO_DAY = 1.0 / DAY_TO_NORMALIZED_TIME;
//角度转换
const double DEGREE_TO_RAD = pi / 180.0;
const double RAD_TO_DEGREE = 180.0 / pi;
const double ARCSEC_TO_RAD = pi / (180.0 * 3600.0);
const double RAD_TO_ARCSEC = (180.0 * 3600.0) / pi;

//面质比归一化 I 国际单位 m^2/kg
const double IAM_TO_NORMALIZED_AM = MASS_OF_EARTH / ae / ae;
const double NORMALIZED_AM_TO_IAM = 1.0 / IAM_TO_NORMALIZED_AM;

namespace Identifier
{
	// to define variables
	enum varID
	{
		VAR_UNDEFINED = 0,
		// 1-/2-digit numbers are for general type id
		VAR_CARTESIAN = 1,
		VAR_KEPLERIAN = 2,
		VAR_NONSINGULAR1 = 3,
		VAR_NONSINGULAR2 = 4,
		VAR_DRAG_COEF = 10,
		VAR_SRP_COEF = 11,
		VAR_PROCNOISE_GM1_RTN = 12,
		VAR_SATELLITE_MASS_THRUST = 13,
		// 3-/4-digit numbers are for specific type id
		VAR_CARTX = 101,
		VAR_CARTY,
		VAR_CARTZ,
		VAR_CARTVX,
		VAR_CARTVY,
		VAR_CARTVZ,
		VAR_CARTPOS,
		VAR_CARTVEL,
		VAR_SEMIMAJORAXIS = 201,
		VAR_ECCENTRICITY,
		VAR_INCLINATION,
		VAR_RAAN,
		VAR_ARGPERIGEE,
		VAR_MEANANOMALY,
		VAR_TRUEANOMALY,
		VAR_ECCANOMALY,
		VAR_R,
		VAR_U,
		VAR_A_NSG2 = 401,
		VAR_XI_NSG2,
		VAR_ETA_NSG2,
		VAR_H_NSG2,
		VAR_K_NSG2,
		VAR_LAMBDA_NSG2,
		VAR_U_NSG2,
		VAR_UTILDE_NSG2,
		VAR_H_GEODETIC = 501,
		VAR_LAMBDA_GEODETIC,
		VAR_PHI_GEODETIC,
		VAR_DRAG_BSTAR = 1001,
		VAR_DRAG_BSTAR_TLE = 1002,
		VAR_DRAG_CD = 1003,
		VAR_DRAG_AMR = 1004,
		VAR_DRAG_RHO = 1005,
		VAR_SRP_AMR = 1101,
		VAR_PROCNOISE_GM1_RADIAL = 1201,
		VAR_PROCNOISE_GM1_TRANSVERSE = 1202,
		VAR_PROCNOISE_GM1_NORMAL = 1203
	};
}
#endif /* CDC_UTILS_CONSTANTS_HPP */
