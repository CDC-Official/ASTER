//
//  core_iod.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#ifndef CDC_CORE_IOD_HPP
#define CDC_CORE_IOD_HPP

#include <cmath>
#include <array>
#include <algorithm>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem> // C++17 提供

#include "utils_constants.hpp"
#include "frame_earth.hpp"
#include "model_ephemeris_analytic.hpp"
#include "io_common.hpp"
#include "math_basic.hpp"
#include "orbit_basic.hpp"
#include "orbit_statedef.hpp"
#include "utils_overwatch.hpp"
#include <sofa.h>

using namespace Eigen;
using namespace orbit;
using namespace bastion;

// IOD约束类型
#define IOD_CONSTRAINT_NONE 0
#define IOD_CONSTRAINT_ORIENTATION 1
#define IOD_CONSTRAINT_CIRCULAR 2
#define IOD_CONSTRAINT_ORIENTATION_CIRCULAR 3
// 方向条件用于排除原始 GL 叉乘方程允许的反向视线解。
// 类型 4 同时使用方向条件与指定偏心率条件：每次观测有 3 个基本行、
// 1 个偏心率行和 1 个方向行，共 5 行；指定偏心率时始终包含方向条件。
#define IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY 4
// IOD动力学模型
#define IOD_DYNAMICS_TWOBODY 1
#define IOD_DYNAMICS_J2 2
#define IOD_DYNAMICS_J2_SOLAR 3
#define IOD_DYNAMICS_J2_LUNAR 4
#define IOD_DYNAMICS_J2_SOLAR_LUNAR 5
// IOD迭代方式
#define IOD_ITERATION_SIMPLE 1
#define IOD_ITERATION_NEWTON 2
// IOD的观测数据类型，暂不支持混合类型数据
#define OBS_TYPE_RADEC 1
#define OBS_TYPE_AZEL 2
#define OBS_TYPE_RADEC_RHO 3
#define OBS_TYPE_AZEL_RHO 4
#define OBSTYPE_XYZ 5
// IOD是否输出迭代信息
#define IOD_PRINT_ITERATION_INFO_OFF 0
#define IOD_PRINT_ITERATION_INFO_ON 1
// IOD是否输出归一化结果
#define IOD_SOLUTION_NORMALIZE 1
#define IOD_SOLUTION_KMDEG 2
#define IOD_SOLUTION_MDEG 3
// IOD初值方式
#define IOD_INITIAL_VALUE_DEFAULT 1     //默认初值，书上那种
#define IOD_INITIAL_VALUE_SPECIFIED 2     //指定初值，需要定轨前提供的
namespace iod
{
	struct IodObservation
	{
		// 观测资料应该全部归一化
		int type;        // 观测资料类型，可以参考老师的写法
		int number;        // 观测条数
		ArrayXd time_MJDUTC;    //时间，MJDUTC
		ArrayXd tau;    //时间，做差后的归一化时间序列
		ArrayXd az_rad, el_rad;        // 方位角、高度角，rad
		ArrayXd ra_rad, dec_rad;    // 赤经、赤纬，rad
		ArrayXd rho;        // 测距
		ArrayXd X, Y, Z;    // 点位资料

		//FIXME: hardcode
		double sigma_RaDec = 3.0 * ARCSEC_TO_RAD; // 光学测角的噪声标准差，弧度制，默认3角秒
	};

	struct IodObserver
	{
		// 测站资料
		int observer_type;        // 测站的类型
		MatrixXd positions_gcrs_norm_emat;                // 各观测历元的测站 GCRS 位置，归一化长度；3×N，每列对应一条观测。
		Vector3d position_gcrs_epoch_norm;    // 参考历元的测站 GCRS 位置，归一化长度。
	};

	struct IodInitialValue
	{
		// IOD初值
		int type;//由宏来决定，1为默认初值，2为指定初值，需要定义下面两个之一
		Keplerian initKep; //norm
		Cartesian initCart;// norm
	};

	/**
	 * @brief 一轮 GL 迭代结束时的参考状态、改正量和测角拟合质量。
	 *
	 * 状态均位于 IodJob::epochMJDUTC（UTC MJD），坐标轴为 GCRS。
	 * r0v0_norm 前三项为位置/ae，后三项为速度/(ae/UNITTIME)；
	 * dr0v0_norm 使用相同尺度。三个 RMS 均为 rad，采用本条状态重新计算。
	 * 指定初值记为第 0 轮，默认线性初值记为第 1 轮；这两种起始记录的
	 * 改正量及最大改正量为 NaN，表示尚未发生一次状态更新。
	 */
	struct IodIterationRecord
	{
		int iteration = 0;
		std::array<double, 6> r0v0_norm = {};
		std::array<double, 6> dr0v0_norm = {};
		double rms_ra_cosdec_rad = 0.0; // RMS((RA_obs-RA_calc)*cos(Dec_obs))。
		double rms_dec_rad = 0.0;       // RMS(Dec_obs-Dec_calc)。
		double rms_angle_rad = 0.0;     // sqrt((RMS_RA^2+RMS_Dec^2)/2)。
		double max_dr0v0_norm = 0.0;    // 六个归一化改正分量的最大绝对值。
	};

	struct IodSolution
	{
		// IOD的解以及相关信息
		int flag;    // 收敛标识：正数为收敛轮次，-1 为达到迭代上限，-2 为数值失败。
		Keplerian kep;
		Cartesian cart;

		// 归一化残差，各种残差先不声明长度，应当在预处理时按照需要resize来调整
		ArrayXd resAz, resEl;
		ArrayXd resRa, resDec;        // O-C，rad：ΔRA·cos(Dec_obs)、ΔDec
		ArrayXd resRho;                // 测距
		ArrayXd resX, resY, resZ;    // 点位资料
		double RMS_angle;            // 角度量纲的RMS
		double RMS_length;            // 长度量纲的RMS
		double RMS_sum;                // 总RMS
		double COD[6];    //		6D的决定系数
		double resRaCOD1d, resRaCOD2d, resDecCOD1d, resDecCOD2d;//残差序列的一次、二次决定系数
		double angle_norm;            // 定轨历元时，测站矢量和轨道面的夹角，归一化

		double resRa_Q, resDec_Q, resRa_p, resDec_p;// Ljung–Box 遗留字段；当前未计算时为 -1。
		double J_likelihood;

		// 尝试性重构核心每次调用时清空，并依次保存起始状态及每轮更新结果。
		// display_iterations 只控制屏幕显示，迭代历史始终保存。
		std::vector<IodIterationRecord> iterations;

		// 逐轮显示文字单独留存，供任务报告和 GLML 所选来源重放；每次求解清空。
		// 它属于计算报告，不进入 tracer 事件日志。
		std::string iteration_report;

		/**
		 * @brief 输出 IOD 最终轨道及测角残差的均方根，可附带来源迭代信息。
		 * @param[in,out] output 接收结果文字的输出流；名称宽度为 40，状态按两行三分量排列。
		 * @param[in] include_iteration_info 为 true 时显示 flag、迭代次数和最大改正量；
		 *     为 false 时仅显示轨道与 RMS，适用于 GLML 边界处理后的最终轨道。
		 * @details 读取 cart、kep、resRa、resDec 和 RMS_angle；东西、南北 RMS 由保存的
		 * 残差分别求均方根。附带迭代信息时，读取 flag 和 iterations 的最后一条记录。
		 * Cartesian 从 ae、ae/UNITTIME 换算为 km、km/s，
		 * Keplerian 从 ae、rad 换算为 km、deg，最大改正量保留归一化单位。
		 * 东西、南北 RMS 分别占一行，每行依次显示 rad 和 arcsec；合成 RMS 也显示两种单位。
		 * 弧度值采用科学计数法，尾数保留两位小数；角秒值按定点格式保留两位小数。
		 * 轨道均属于调用方的参考历元和 GCRS，不重新传播或计算残差。
		 * 本函数不输出分区标题、任务耗时和内存；不修改结果及输出流的格式设置。
		 */
		void PrintSummary(std::ostream& output, bool include_iteration_info = true) const;
	};
	struct IodJob
	{
		int dynamicsModel;        // 动力学模型
		int constraintType;        // 约束类型
		int iterMethod;        // 迭代模型，简单迭代or微分改正
		double iterTOL;        // 认为改正量收敛的tolerance
		int iterCountMax;        // 最大迭代次数
		double epochMJDUTC;            // 定轨历元

		IodObservation obs;        // 观测数据
		IodObserver observer;
		IodSolution solution;
		IodInitialValue initValue;

		std::string fn_obs;
		std::string fn_solution;
		std::string fn_residual;
		std::string fn_report;
		std::string fn_EOP;

		std::string objectId;
		bool printIterInfo;

		int solutionUnit;//仅改变输出时的单位，不改变结构体的归一化单位

		//int iterCountMaxEcc;    //约束的偏心率，241127
		double constraintEcc;

		/**
		 * @brief 输出 GL 初轨确定的动力学、迭代、约束和观测误差设置。
		 * @param[in,out] output 接收 IOD settings 分区的输出流；名称宽度为 40，参数按键值对排列。
		 * @details 读取 dynamicsModel、iterMethod、constraintType、constraintEcc、initValue.type、
		 * iterCountMax、iterTOL 和 obs.sigma_RaDec。动力学按五种 IOD 宏定义显示，
		 * iterTOL 为六个归一化状态改正分量的最大绝对值阈值；sigma_RaDec 为东西、南北
		 * 两个独立高斯误差共同的标准差，由 rad 换算为 arcsec 后显示。
		 * 本函数不修改设置及输出流的格式设置，也不读取文件或执行定轨计算。
		 */
		void Display(std::ostream& output) const;

	};

}

// 根据动力学模型计算拉格朗日传播系数。
void Lagrange_Coefficients(int dynamics_model, double epoch_of_IOD, double tau, double* X, double* LagCoeff);

// 计算传播系数对参考位置的偏导数。
void PartialFG(double tau, double x0, double y0, double z0, double* PFG);

#endif
