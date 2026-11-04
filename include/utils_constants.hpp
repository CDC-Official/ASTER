/**
 * @file utils_constants.hpp
 * @brief 数学常数、天文参数及程序归一化尺度。
 * @details 为轨道预报、初轨计算、观测仿真和坐标转换提供常数。
 * 通用常数直接放在 constants 中；重力场、参考椭球和光压参数按计算用途分组。
 * Scale 定义程序的数值单位，Conversion 保存单位换算系数。来源和调用位置见 doc/11_constants.md。
 */
#ifndef CDC_UTILS_CONSTANTS_HPP
#define CDC_UTILS_CONSTANTS_HPP

#include <cmath>
#include <sofam.h>

namespace constants
{

// 圆周率，用于角度换算、轨道相位和角度回绕。
const double PI = 3.141592653589793238462643;

// 真空光速，用于光行时、伪距和后牛顿加速度计算；单位 m/s，采用 SI 定义值。
constexpr double SPEED_OF_LIGHT_MPS = 299792458.0;

// 万有引力常数，用于由地球 GM 求质量尺度，进而换算密度和面质比。
// 单位 m³/(kg·s²)，采用 CODATA 2022 推荐值：physics.nist.gov/constants。
constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;

// 天文单位，用于将解析太阳轨道的半长轴从 AU 换成 m。
// 采用 IAU 2012 Resolution B2 的固定长度定义，单位 m。
constexpr double ASTRONOMICAL_UNIT = 149597870700.0;

// JD 与 MJD 互换时加减的日数，也用于向 SOFA 和 JPL 历表接口传入历元；单位 d。
constexpr double MJD_OFFSET = DJM0;

// 将以日表示的时间间隔换成秒，用于积分步长、钟差和历元计算；取 SOFA 的日长定义。
constexpr double SECONDS_PER_DAY = DAYSEC;

// 将距 J2000.0 的日数换成儒略世纪，供解析日月历表、岁差和恒星时公式使用。
constexpr double DAYS_PER_CENTURY = DJC;

// 计算距 J2000.0 的时间间隔时使用的起点，分别以 JD 和 MJD 表示；取自 SOFA。
constexpr double J2000_JD = DJ00;
constexpr double J2000_MJD = DJM00;

// 将 TDB 相容的日月 GM 换成 TT 相容值时乘此比例。
// 将 TT 重力场的 GM 和参考半径用于 TDB 方程时，两者都除以此比例；L_G、L_B 取自 SOFA。
constexpr double TDB_TO_TT_SCALE = (1.0 - ELG) / (1.0 - ELB);

// 用于 E-BCRS/TDB 轨道预报中的天体引力，以及单程测距的引力传播时延。
// GM 采用 Park 等（2021，AJ 161:105）DE440/DE441 论文表 2。
namespace DE440
{
    // TDB 相容 GM，单位 m³/s²。顺序为水星、金星、地球、火星、木星、土星、天王星、海王星、冥王星、月球、太阳。
    // 下标为 JPL 天体编号减 1；火星至冥王星使用包含卫星的行星系统 GM。
    constexpr double GM_TDB[11] = {
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
}

// 用于近地轨道的地球质点引力和共转大气速度计算。
// 参数取自 IERS Conventions 2010 第 1 章及第 6.1 节。
namespace IERS2010
{
    // GCRS/TT 地球质点引力使用的 GM，同时定义轨道积分和根数转换的 GM 尺度；单位 m³/s²。
    constexpr double EARTH_GM_TT = 3.986004415e14;

    // 计算共转大气速度 omega × r，再与卫星速度相减得到阻力模型的相对速度。
    // 取表 1.1 的名义平均地球自转角速度，单位 rad/s。
    constexpr double EARTH_ROTATION_RATE = 7.292115e-5;
}

// 用于测站经纬度、高度与地固坐标之间的转换，也为大气高度和地影计算提供地球形状。
// 半径和扁率配套采用 WGS84（NGA.STND.0036），与 SOFA eform 的 WGS84 分支一致。
namespace Wgs84
{
    // 测站和大气参考椭球的赤道半径；圆柱地影模型也以此作为地球半径，单位 m。
    constexpr double EQUATORIAL_RADIUS = 6378137.0;

    // 扁率用于确定椭球形状，并计算大气参考面随纬度的变化；无量纲。
    constexpr double FLATTENING = 1.0 / 298.257223563;

    // 大地坐标转换用它计算卯酉圈曲率半径和地固 Z 分量；无量纲。
    constexpr double ECCENTRICITY_SQUARED = FLATTENING * (2.0 - FLATTENING);
}

// 用于地球非球形引力计算，也供初轨、解析 STM 和瞬转平计算取球谐参数。
// GM、参考半径与 data_gravityfield.hpp 中的系数配套，取自 CSR GGM03C 文件。
// 来源：https://download.csr.utexas.edu/puace/GGM03/GGM03_Archive.zip
namespace GGM03C
{
    // 确定球谐引力势及加速度的整体大小，单位 m³/s²，为 TT 相容值。
    constexpr double GM = 3.986004415e14;

    // 球谐展开中 (REFERENCE_RADIUS / r)^n 的参考半径，单位 m。
    constexpr double REFERENCE_RADIUS = 6378136.3;
}

// 用于数值轨道预报的光压加速度和瞬转平计算的光压摄动项。
// 压力乘反射系数和面质比得到加速度；参考距离与压力采用 Orekit 13.1 SolarRadiationPressure 的配套默认值。
namespace SolarRadiation
{
    // 压力标定处到太阳的距离；卫星处压力按 (REFERENCE_DISTANCE / 日星距离)^2 缩放，单位 m。
    constexpr double REFERENCE_DISTANCE = 149597870000.0;

    // 上述距离处的太阳辐射压力，单位 N/m²。
    constexpr double REFERENCE_PRESSURE = 4.56e-6;
}

// 轨道积分、根数转换和摄动计算使用的数值单位。输入物理量除以相应尺度，输出时乘回。
// 长度、时间、GM 尺度满足 UNIT_GM_M3_S2 × UNIT_TIME_S² / UNIT_LENGTH_M³ = 1，使近地二体方程中的引力参数为 1。
namespace Scale
{
    // 各天体 GM 除以此值后进入归一化引力方程，单位 m³/s²。
    constexpr double UNIT_GM_M3_S2 = constants::IERS2010::EARTH_GM_TT;

    // 位置、轨道半长轴和模型参考半径都除以此值；取 GGM03C 参考半径，单位 m。
    const double UNIT_LENGTH_M = constants::GGM03C::REFERENCE_RADIUS;

    // 时间尺度，单位 s。积分时间间隔和步长除以此值，速度则按 UNIT_LENGTH_M / UNIT_TIME_S 换算。
    const double UNIT_TIME_S = sqrt(pow(constants::Scale::UNIT_LENGTH_M, 3) / constants::Scale::UNIT_GM_M3_S2);

    // 用于大气密度和卫星面质比的归一化，由地球 GM/G 求得，单位 kg。
    const double UNIT_MASS_KG = constants::IERS2010::EARTH_GM_TT / constants::GRAVITATIONAL_CONSTANT;

    // TDB 天体 GM 除以程序 GM 尺度，供 E-BCRS 牛顿项、后牛顿项使用。
    // 地球的 TDB GM 与归一化尺度不同，因此地球项也保留相应系数。
    constexpr double NORMALIZED_BODY_GM_TDB[11] = {
        constants::DE440::GM_TDB[0] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[1] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[2] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[3] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[4] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[5] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[6] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[7] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[8] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[9] / constants::Scale::UNIT_GM_M3_S2,
        constants::DE440::GM_TDB[10] / constants::Scale::UNIT_GM_M3_S2
    };

    // 月球 GM 先换成 TT 相容值，再归一化，用于近地月球质点摄动、初轨修正和解析平均项。
    constexpr double NORMALIZED_MOON_GM_TT = (constants::DE440::GM_TDB[9] * constants::TDB_TO_TT_SCALE) / constants::Scale::UNIT_GM_M3_S2;

    // 太阳 GM 先换成 TT 相容值，再归一化，用于近地太阳质点摄动、初轨修正和解析平均项。
    constexpr double NORMALIZED_SUN_GM_TT = (constants::DE440::GM_TDB[10] * constants::TDB_TO_TT_SCALE) / constants::Scale::UNIT_GM_M3_S2;
}

// 用于读取轨道、观测和卫星参数时换算单位，以及输出计算结果时恢复物理单位。
// 每个系数都与待换算数值相乘，名称中的 TO 表示换算方向。
namespace Conversion
{
    // 输入的 m/s、km/s 换成积分速度，或将积分速度恢复为输出单位。
    const double MPS_TO_NORMALIZED_VELOCITY = constants::Scale::UNIT_TIME_S / constants::Scale::UNIT_LENGTH_M;
    const double KMPS_TO_NORMALIZED_VELOCITY = 1000.0 * constants::Scale::UNIT_TIME_S / constants::Scale::UNIT_LENGTH_M;
    const double NORMALIZED_VELOCITY_TO_MPS = 1.0 / constants::Conversion::MPS_TO_NORMALIZED_VELOCITY;
    const double NORMALIZED_VELOCITY_TO_KMPS = 1.0 / constants::Conversion::KMPS_TO_NORMALIZED_VELOCITY;

    // 轨道位置、测站位置及距离在 m、km 和积分长度单位之间换算。
    const double M_TO_NORMALIZED_LENGTH = 1.0 / constants::Scale::UNIT_LENGTH_M;
    const double KM_TO_NORMALIZED_LENGTH = 1000.0 * 1.0 / constants::Scale::UNIT_LENGTH_M;
    const double NORMALIZED_LENGTH_TO_M = 1.0 / constants::Conversion::M_TO_NORMALIZED_LENGTH;
    const double NORMALIZED_LENGTH_TO_KM = 1.0 / constants::Conversion::KM_TO_NORMALIZED_LENGTH;

    // 将预报时长、步长换成积分时间，或将积分时间恢复为秒、分、时、日。
    const double SEC_TO_NORMALIZED_TIME = 1.0 / constants::Scale::UNIT_TIME_S;
    const double MIN_TO_NORMALIZED_TIME = 60.0 / constants::Scale::UNIT_TIME_S;
    const double HOUR_TO_NORMALIZED_TIME = 3600.0 / constants::Scale::UNIT_TIME_S;
    const double DAY_TO_NORMALIZED_TIME = 86400.0 / constants::Scale::UNIT_TIME_S;
    const double NORMALIZED_TIME_TO_SEC = 1.0 / constants::Conversion::SEC_TO_NORMALIZED_TIME;
    const double NORMALIZED_TIME_TO_MIN = 1.0 / constants::Conversion::MIN_TO_NORMALIZED_TIME;
    const double NORMALIZED_TIME_TO_HOUR = 1.0 / constants::Conversion::HOUR_TO_NORMALIZED_TIME;
    const double NORMALIZED_TIME_TO_DAY = 1.0 / constants::Conversion::DAY_TO_NORMALIZED_TIME;

    // 将轨道根数、测站经纬度和测角观测换成 rad 参与计算，输出时恢复为度或角秒。
    const double DEGREE_TO_RAD = constants::PI / 180.0;
    const double RAD_TO_DEGREE = 180.0 / constants::PI;
    const double ARCSEC_TO_RAD = constants::PI / (180.0 * 3600.0);
    const double RAD_TO_ARCSEC = (180.0 * 3600.0) / constants::PI;

    // 将卫星面质比从 m²/kg 换成阻力、光压函数使用的归一化量，或作反向换算。
    const double IAM_TO_NORMALIZED_AM = constants::Scale::UNIT_MASS_KG / constants::Scale::UNIT_LENGTH_M / constants::Scale::UNIT_LENGTH_M;
    const double NORMALIZED_AM_TO_IAM = 1.0 / constants::Conversion::IAM_TO_NORMALIZED_AM;
}
}

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
