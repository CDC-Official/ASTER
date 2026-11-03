//
//  utils_constants.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#ifndef CDC_UTILS_CONSTANTS_HPP
#define CDC_UTILS_CONSTANTS_HPP

#include <cmath>

const double GE = 3.986004418E14;                                      // Geocentric gravitational constant(EGM96 value)
const double
	GEOCENTRIC_GRAVITATIONAL_CONSTANT = 3.986004418E14;      // Geocentric gravitational constant(EGM96 value) 跟老师一致
const double
	HELIOCENTRIC_GRAVITATIONAL_CONSTANT = 1.32712442076E20;  // HELIOCENTRIC gravitational constant(EGM96 value)
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

// 主要天体的 GM / GM_Earth，按 JPL 编号 1—11 排列：
// 水星、金星、地球、火星、木星、土星、天王星、海王星、冥王星、月球、太阳。
// 地月预报和测距引力时延共用这组比值，各模型分别提供所用的地球 GM。
constexpr double MAJOR_BODY_EARTH_GM_RATIO[11] = {
    0.0552730669322440,
    0.8149980859574224,
    1.0000000000000000,
    0.1074468866586270,
    317.8941938154875402,
    95.1845042355296442,
    14.5372555667923731,
    17.1513287281893625,
    0.0024473129307027,
    0.0123000370327340,
    332946.0487730671302415
};

const double MOON_EARTH_MASS_RATIO = 0.0123000383;  // 月球/地球质量比
const double SUN_EARTH_MASS_RATIO = HELIOCENTRIC_GRAVITATIONAL_CONSTANT / GEOCENTRIC_GRAVITATIONAL_CONSTANT;  //太阳/地球质量比
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
