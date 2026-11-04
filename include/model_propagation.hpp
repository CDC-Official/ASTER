//
// Created by 陈酊存 on 2024/1/23.
//

#ifndef CDC_LIB_INCLUDE_MODEL_PROPAGATION_HPP_
#define CDC_LIB_INCLUDE_MODEL_PROPAGATION_HPP_

#include <iostream>
#include <algorithm>
#include <array>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
#include "utils_constants.hpp"
#include "orbit_statedef.hpp"
#include "frame_earth.hpp"
#include "model_atmosphere.hpp"
#include "accel_solar_radiation.hpp"
#include "accel_thirdbody.hpp"
#include "model_ephemeris_jpl.hpp"
#include "model_ephemeris_analytic.hpp"
#include "time_transform.hpp"
#include "accel_gravity.hpp"
#include "utils_overwatch.hpp"
#include "model_satellite.hpp"
#include "orbit_typedef.hpp"

using namespace Eigen;
using namespace Identifier;
using orbit::Keplerian;
using orbit::Cartesian;
using orbit::NonSingular2;

using std::cout;
using std::vector;
using std::endl;

// 日月第三体摄动的相关设置
#define SOLAR_GRAVITY_OFF 0
#define SOLAR_GRAVITY_JPL 1
#define SOLAR_GRAVITY_ANALYTICAL 2

#define LUNAR_GRAVITY_OFF 0
#define LUNAR_GRAVITY_JPL 1
#define LUNAR_GRAVITY_ANALYTICAL 2

// 输出星历的单位
#define STATE_UNIT_NORMALIZE 1
#define STATE_UNIT_KM_SEC_DEGREE 2
#define STATE_UNIT_M_SEC_DEGREE 3

// =============================================================================
// 用于 GCRS / TT 下的结构和函数
// =============================================================================

// 积分器专用参数集中放在这里；共同的起点、时长和步长见下方 OrbitPropagationSetting。
struct RkfSetting
{
	// RKF78 用同一组阶段导数构造两种阶次的结果，取两者之差估计本步截断误差。
	// 这里的 TE 是六分量归一化误差绝对值之和，不是以米计的整段轨道误差。
	double maxTE_norm = 1E-12;
	int iterCount_max = 10; // 自适应单步最多试算次数（含首次，须 >= 1），不是整段预报步数。
	// false：误差超限后减半步长重算；true：保留步长并警告，误差仍正常计算。
	bool fixed_step = false;
};

struct AdamsSetting
{
	// 多步法复用等间距历史节点的导数，减少每一步重新计算动力学的次数。
	int steps = 9; // 历史导数节点数，支持 2..9；不是 RKF78 的阶段数。
	// false 为 PE：预测(P)后计算末端导数(E)，留给下一步使用。
	// true 为 PECE：预测后计算导数，再改正(C)一次，并在改正状态重新计算导数。
	bool use_correction = false; // 只改正一次，不迭代求解隐式方程。
};

// GCRS/TT 预报的共同设置，以及各积分器的专用设置。
// 调用示例：job.opn_setting.rkf.maxTE_norm、job.opn_setting.adams.steps。
struct OrbitPropagationSetting
{
	double t0_MJD_UTC = 0; // 积分起点 UTC。
	double t0_MJD_TT = 0;  // 对应的积分起点 TT。
	double duration_norm = 0; // 归一化总时长；正值前向，负值反向，零值仅记录初态。
	double h0_norm = 0; // 归一化步长；自适应 RKF 作试探值，定步长 RKF／Adams 作固定值。
	int display_progress = 2; // 0 静默，1 每 10%，2 每步／试算；不抑制警告和错误。

	RkfSetting rkf; // RKF78 的误差控制；Adams 起步及不足整步的尾段也会用到。
	AdamsSetting adams; // Adams 的节点数及 PE／PECE 模式。

	/**
	 * @brief 显示当前积分器的步长、误差控制与起步设置。
	 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
	 * @param[in] integrator_id 积分器编号：101 为自适应 RKF78，102 为定步长 RKF78，201/202 为 Adams PE/PECE。
	 * @param[in,out] output 接收设置文本的输出流；函数保持该流原有的格式状态。
	 * @param[in] precision 普通设置值的有效数字位数；步长和 TOL 固定保留两位小数与两位科学记数法尾数。
	 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
	 * @details 读取 h0_norm、rkf 和 adams，不修改设置或执行积分。步长乘 UNITTIME 恢复为 s。
	 * TOL 比较内部归一化状态的局部截断误差估计；定步长 RKF78 仅将其作为告警阈值。
	 * Adams 显示 RKF78 起步、尾段所用的容差及试算次数上限。
	 */
	void Display(bastion::Overwatch& tracer, int integrator_id, std::ostream& output = std::cerr, int precision = 14, int label_width = 31) const;
};

struct OrbitPropagationNModel
{
	// 250529 数值法轨道预报的动力学模型Model，默认地球为中心天体，仅支持考虑日月第三体摄动
	/* ----------------------- controlling parameters ---------------------- */
	bool use_earth_NSG = false;                    // 是否使用中心天体非球形引力摄动，不使用即为二体质点
	bool use_SRP = false;                    // 是否使用光压摄动
	bool use_drag = false;                   // 是否使用大气阻力摄动
	int earth_zonal_order = 0;                     // 地球带谐项 C_n0 的最高 n；0 关闭该部分
	int earth_tesseral_degree = 0;                 // 地球田谐项的最高 n（含扇谐项，m=1..n）；0 关闭该部分
	int solar_gravity_config = SOLAR_GRAVITY_OFF;  // 太阳第三体质点摄动
	int lunar_gravity_config = LUNAR_GRAVITY_OFF;  // 月亮第三体的质点摄动

	// 需要注意 开启月球非球形必须要读JPL历表
	bool use_lunar_NSG = false;              // 是否使用月球非球形引力摄动（不含质点引力），进阶选项，仅在考虑JPL月球第三体后才能额外开启此选项（因为月固坐标转换必须要读JPL文件）
	int lunar_zonal_order = 0;                     // 带谐项 order
	int lunar_tesseral_degree = 0;                 // 田谐项 degree 理论上要小于等于order

	// 密度表注明的固定空间天气条件，仅用于说明资料；缺少表头时保留 NaN。
	double atmosphere_Ap = std::numeric_limits<double>::quiet_NaN();
	double atmosphere_F107 = std::numeric_limits<double>::quiet_NaN();

	// 卫星变量，赋值质量、面质比等参数，主要用于大气和光压
	model::common::SphereSatellite sat;

	/* ----------------------- 未启用预留参数 ---------------------- */
	bool UseSolidTide = false;               // 是否使用中心天体固体潮摄动
	bool UseOceanTide = false;               // 是否使用中心天体海潮摄动
	bool UsePoleTide = false;                // 是否使用中心天体极潮摄动
	bool UseNoise = false;                   // 是否使用过程噪声
	bool UseManeuver = false;                // 是否使用轨道机动
	bool UseEmpirical = false;               // 是否使用经验力模型
	bool SolveVarEq = false;                 // 是否积分变分方程

	/**
	 * @brief 显示实际启用的引力、阻力和光压模型及其物理参数。
	 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
	 * @param[in,out] output 接收模型文本的输出流；函数保持该流原有的格式状态。
	 * @param[in] precision 面积、AMR、系数和大气表头条件等普通数值的有效数字位数。
	 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
	 * @details 读取地月非球形开关及截断、日月历表模式、阻力和光压开关、sat 和密度表条件。
	 * 面积乘 ae² 恢复为 m²，AMR 乘 ae²/MASS_OF_EARTH 恢复为 m²/kg。
	 * 月球引力截断为 0 / 0 时注明仅考虑质点引力；缺少密度表头条件时显示 Not provided。
	 * 光压使用 Cannonball 模型，采用地球圆柱地影，进入阴影后光压为零。函数不修改动力学参数。
	 */
	void Display(bastion::Overwatch& tracer, std::ostream& output = std::cerr, int precision = 14, int label_width = 31) const;
};    //预报动力学模型设定，用于GCRS/TT
// 默认轨道预报动力学模型配置命名空间
namespace DefaultOrbitPropagationNModel
{
	// 参考导师的注释和设计思路：
	// 1. 目前预定义的模型为中等精度（moderate precision） 其实偏高
	// 2. 适用于四种轨道类型：低轨圆轨道（LEO）、中轨圆轨道（MEO）、高轨圆轨道（GEO）和大偏心率轨道（HEO）
	// 3. 预定义的参数包括天体引力、历表模型、重力场模型、固体潮、光压、大气等开关

	// 参考导师的思路，但为了避免大量栈空间和运行时内存消耗，采用工厂函数而非预定义const对象的方式
	// 这些都不是const，可以直接传入轨道预报函数使用！

	// moderate-precision OPN model for LEO
	// 参考导师配置：使用固体潮、光压、大气阻力，高精度重力场
	inline OrbitPropagationNModel getLEOModerate()
	{
		OrbitPropagationNModel model;

		// LEO特点：大气阻力重要，需要高精度重力场
		model.use_earth_NSG = true;                    // 使用非球形引力
		model.use_SRP = true;                    // 使用光压（参考导师配置）
		model.use_drag = true;                   // 使用大气阻力（关键）
		model.earth_zonal_order = 30;                  // 高精度带谐项（参考导师GravDeg=30）
		model.earth_tesseral_degree = 30;              // 高精度田谐项
		model.solar_gravity_config = SOLAR_GRAVITY_JPL;     // 使用JPL星历（参考导师）
		model.lunar_gravity_config = LUNAR_GRAVITY_JPL;     // 使用JPL星历

		// 固体潮，暂未支持，开了也没用（参考导师LEO配置使用固体潮）
		model.UseSolidTide = true;

		return model;
	}

	// moderate-precision OPN model for MEO
	// 参考导师配置：不使用固体潮和大气阻力，使用光压，中等精度重力场
	inline OrbitPropagationNModel getMEOModerate()
	{
		OrbitPropagationNModel model;

		// MEO特点：大气阻力可忽略，光压和第三体摄动重要
		model.use_earth_NSG = true;
		model.use_SRP = true;                    // 使用光压（参考导师配置）
		model.use_drag = false;                  // 不使用大气阻力（参考导师配置）
		model.earth_zonal_order = 10;                  // 中等精度（参考导师GravDeg=10）
		model.earth_tesseral_degree = 10;
		model.solar_gravity_config = SOLAR_GRAVITY_JPL;     // 使用JPL星历
		model.lunar_gravity_config = LUNAR_GRAVITY_JPL;

		// 不使用固体潮（参考导师MEO配置）
		model.UseSolidTide = false;

		return model;
	}

	// moderate-precision OPN model for GEO
	// 参考导师配置：与MEO类似，不使用固体潮和大气阻力，使用光压
	inline OrbitPropagationNModel getGEOModerate()
	{
		OrbitPropagationNModel model;

		// GEO特点：第三体摄动和光压最重要
		model.use_earth_NSG = true;
		model.use_SRP = true;                    // 使用光压（参考导师配置）
		model.use_drag = false;                  // 不使用大气阻力（参考导师配置）
		model.earth_zonal_order = 10;                  // 中等精度（参考导师GravDeg=10）
		model.earth_tesseral_degree = 10;
		model.solar_gravity_config = SOLAR_GRAVITY_JPL;     // 使用JPL星历
		model.lunar_gravity_config = LUNAR_GRAVITY_JPL;

		// 不使用固体潮，暂未支持，开了也没用（参考导师GEO配置）
		model.UseSolidTide = false;

		return model;
	}

	// moderate-precision OPN model for HEO
	// 参考导师配置：类似LEO，使用所有摄动，高精度重力场
	inline OrbitPropagationNModel getHEOModerate()
	{
		OrbitPropagationNModel model;

		// HEO特点：覆盖不同高度，需要考虑所有摄动
		model.use_earth_NSG = true;
		model.use_SRP = true;                    // 使用光压（参考导师配置）
		model.use_drag = true;                   // 使用大气阻力（近地点可能很低）
		model.earth_zonal_order = 30;                  // 高精度（参考导师GravDeg=30）
		model.earth_tesseral_degree = 30;
		model.solar_gravity_config = SOLAR_GRAVITY_JPL;     // 使用JPL星历
		model.lunar_gravity_config = LUNAR_GRAVITY_JPL;

		// 使用固体潮，暂未支持，开了也没用（参考导师HEO配置）
		model.UseSolidTide = true;

		return model;
	}

	// 通用接口函数 - 根据轨道类型获取对应配置 先写着 不一定用得上 最好在main里switch
	// 因为这里面不太好配置tracer
	inline OrbitPropagationNModel getModerate(int orbit_type)
	{
		switch (orbit_type)
		{
			case ORBIT_GEN_TYPE_LOW_NEAR_CIRCULAR:
				return getLEOModerate();
			case ORBIT_GEN_TYPE_MEDIUM_NEAR_CIRCULAR:
				return getMEOModerate();
			case ORBIT_GEN_TYPE_HIGH_NEAR_CIRCULAR:
				return getGEOModerate();
			case ORBIT_GEN_TYPE_HIGHLY_ECCENTRIC_ORBIT:
				return getHEOModerate();
			default:
				return OrbitPropagationNModel{}; // 默认全关闭配置
		}
	}
}


struct OrbitPropagationJob
{
	OrbitPropagationNModel opn_model; // 物理力模型。
	OrbitPropagationSetting opn_setting; // 共同预报设置及各积分器参数。
	EOP opn_EOP;
	vector<AtmosphereData> atm_list;
	ReferenceDensity refDensity;
	Keplerian initKep_norm;
	Cartesian initCart_norm;
};

// =============================================================================
// 用于 EBCRS / TDB 下的预报设置、力模型和任务。
// 本项目的 EBCRS 将 BCRS 原点纯几何平移到地心，保留 BCRS 轴向，时间参数仍为 TDB。
// 同一 TDB 时刻：r_EBCRS = r_BCRS - r_Earth_BCRS；
// 对时间求导：v_EBCRS = v_BCRS - v_Earth_BCRS，加速度也相应作差。
// 这不是完整的相对论 BCRS→GCRS 坐标变换，不能把 EBCRS 直接称为 GCRS。
// =============================================================================
struct OrbitPropagationSetting_EBCRS
{
	int display_progress = 2;      // 0 静默；1 每 10%；2 每个接受步及接受/拒绝试算误差

	double t0_MJD_TDB = 0.0;        // 积分起点 TDB
	double duration_norm = 0.0;     // 归一化总时长；正值前向，负值反向，零值仅记录初态
	double h0_norm = 0.0;           // 归一化试探／固定步长；整段预报取大小，方向由总时长决定
	RkfSetting rkf;                // 与 GCRS 共用误差／步长控制，但积分历元始终为 TDB。
	AdamsSetting adams;            // E-BCRS Cartesian 的等间距历史节点及 PE／PECE 模式。

	/**
	 * @brief 显示 E-BCRS/TDB 积分器的步长、误差控制与起步设置。
	 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
	 * @param[in] integrator_id 积分器编号：101 为自适应 RKF78，102 为定步长 RKF78，201/202 为 Adams PE/PECE。
	 * @param[in,out] output 接收设置文本的输出流；函数保持该流原有的格式状态。
	 * @param[in] precision 普通设置值的有效数字位数；步长和 TOL 固定保留两位小数与两位科学记数法尾数。
	 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
	 * @details 读取 h0_norm、rkf 和 adams，不修改设置或执行积分。步长乘 UNITTIME 恢复为 s，积分历元采用 TDB。
	 * TOL 比较内部归一化状态的局部截断误差估计；定步长 RKF78 仅将其作为告警阈值。
	 * Adams 显示 RKF78 起步、尾段所用的容差及试算次数上限。
	 */
	void Display(bastion::Overwatch& tracer, int integrator_id, std::ostream& output = std::cerr, int precision = 14, int label_width = 31) const;
};

struct OrbitPropagationNModel_EBCRS
{
	// EBCRS/TDB 数值法轨道预报的动力学模型Model。
	// 坐标原点为平移到地心的BCRS，时间变量为TDB，动力学模型后续主要服务于地月空间目标。
	bool use_earth_NSG = false;                    // 是否使用地球非球形引力摄动
	bool use_lunar_NSG = false;                    // 是否使用月球非球形引力摄动
	bool use_SRP = false;                          // 是否使用光压摄动
	bool use_drag = false;                         // 是否使用大气阻力摄动，地月空间通常关闭
	bool use_post_newtonian = false;               // 是否使用后牛顿相对论修正
	bool use_planetary_gravity = false;            // 是否使用除地月日之外的行星质点摄动
	vector<int> planetary_body_id;                 // 需要加入经典行星第三体摄动的 JPL 天体编号，不含地球、月球、太阳
	int earth_zonal_order = 0;                     // 地球带谐项 C_n0 的最高 n；0 关闭该部分
	int earth_tesseral_degree = 0;                 // 地球田谐项的最高 n（含扇谐项，m=1..n）；0 关闭该部分
	int lunar_zonal_order = 0;                     // 月球非球形引力带谐项 order
	int lunar_tesseral_degree = 0;                 // 月球非球形引力田谐项 degree
	int solar_gravity_config = SOLAR_GRAVITY_OFF;  // 太阳质点摄动，目前只支持 JPL
	int lunar_gravity_config = LUNAR_GRAVITY_OFF;  // 月球质点摄动，目前只支持 JPL

	// 密度表头条件供模型显示；NaN 表示原始资料没有提供对应值。
	double atmosphere_Ap = std::numeric_limits<double>::quiet_NaN();
	double atmosphere_F107 = std::numeric_limits<double>::quiet_NaN();

	model::common::SphereSatellite sat;            // 卫星质量、面积、面质比等参数

	/**
	 * @brief 显示 E-BCRS 实际启用的引力、阻力、光压及后牛顿修正。
	 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
	 * @param[in,out] output 接收模型文本的输出流；函数保持该流原有的格式状态。
	 * @param[in] precision 面积、AMR、系数和大气表头条件等普通数值的有效数字位数。
	 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
	 * @details 读取地月非球形开关及截断、日月与行星引力开关、行星 JPL 编号、阻力、光压、后牛顿开关、sat 和密度表条件。
	 * 面积乘 ae² 恢复为 m²，AMR 乘 ae²/MASS_OF_EARTH 恢复为 m²/kg。
	 * 月球引力截断为 0 / 0 时注明仅考虑质点引力；缺少密度表头条件时显示 Not provided。
	 * 日月与行星使用 JPL；摄动天体按 JPL 编号及英文名称显示，非法编号显示 Unknown。
	 * 后牛顿开关对应当前 EIH 修正，使用 11 个 JPL 大天体，独立于经典行星摄动列表。
	 * 光压采用 Cannonball 模型，采用地球圆柱地影，进入阴影后光压为零；函数不修改动力学参数或转换轨道坐标。
	 */
	void Display(bastion::Overwatch& tracer, std::ostream& output = std::cerr, int precision = 14, int label_width = 31) const;
};

struct OrbitPropagationJob_EBCRS
{
	// EBCRS/TDB 轨道预报任务总接口。后续 EBCRS 右函数和积分器应以此结构体作为统一入口。
	OrbitPropagationNModel_EBCRS opn_model;
	OrbitPropagationSetting_EBCRS opn_setting;
	EOP opn_EOP;
	vector<AtmosphereData> atm_list;
	ReferenceDensity refDensity;
	Cartesian initCart_norm;
	Keplerian initKep_norm;
};

#endif // CDC_LIB_INCLUDE_MODEL_PROPAGATION_HPP_
