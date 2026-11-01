//
//  integrator_rkf78.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/28.
//

#ifndef CDC_INTEGRATOR_RKF78_HPP
#define CDC_INTEGRATOR_RKF78_HPP

#include <iostream>
#include <cmath>
#include <iomanip>
#include "utils_constants.hpp"
#include "model_ephemeris_analytic.hpp"
#include "model_propagation.hpp"

// 13 阶段 RKF7(8) 系数。活动预报接口使用下表的八阶权重和 41/840 嵌入误差系数。
// 核对来源：boostorg/odeint boost-1.87.0，runge_kutta_fehlberg78.hpp。
// 单步实现位于 integrator_rkf78.cpp；整段预报调度位于 core_propagation.cpp。
const double A_rkf78[13] = {
	0.0, 2.0 / 27.0, 1.0 / 9.0, 1.0 / 6.0, 5.0 / 12.0, 0.5,
	5.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0, 1.0 / 3.0, 1.0, 0.0, 1.0 };
const double B_rkf78[13][12] = {
	{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	{ 2.0 / 27.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	{ 1.0 / 36.0, 1.0 / 12.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	{ 1.0 / 24.0, 0.0, 1.0 / 8.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
	{
		5.0 / 12.0,
		0.0,
		-25.0 / 16.0,
		25.0 / 16.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0,
		0.0
	},
	{ 1.0 / 20.0, 0.0, 0.0, 1.0 / 4.0, 1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	  0.0 },
	{ -25.0 / 108.0, 0.0, 0.0, 125.0 / 108.0, -65.0 / 27.0, 125.0 / 54.0, 0.0,
	  0.0, 0.0, 0.0, 0.0, 0.0 },
	{ 31.0 / 300.0, 0.0, 0.0, 0.0, 61.0 / 225.0, -2.0 / 9.0, 13.0 / 900.0, 0.0,
	  0.0, 0.0, 0.0, 0.0 },
	{ 2.0, 0.0, 0.0, -53.0 / 6.0, 704.0 / 45.0, -107.0 / 9.0, 67.0 / 90.0, 3.0,
	  0.0, 0.0, 0.0, 0.0 },
	{ -91.0 / 108.0, 0.0, 0.0, 23.0 / 108.0, -976.0 / 135.0, 311.0 / 54.0, -19.0 / 60.0,
	  17.0 / 6.0, -1.0 / 12.0, 0.0, 0.0, 0.0 },
	  // 此项为 45/164；已与 Boost.Odeint RKF78 的第 11 阶段系数核对。
	{ 2383.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -301.0 / 82.0,
	  2133.0 / 4100.0, 45.0 / 82.0, 45.0 / 164.0, 18.0 / 41.0, 0.0, 0.0 },
	{ 3.0 / 205.0, 0.0, 0.0, 0.0, 0.0, -6.0 / 41.0, -3.0 / 205.0, -3.0 / 41.0,
	  3.0 / 41.0, 6.0 / 41.0, 0.0, 0.0 },
	{ -1777.0 / 4100.0, 0.0, 0.0, -341.0 / 164.0, 4496.0 / 1025.0, -289.0 / 82.0,
	  2193.0 / 4100.0, 51.0 / 82.0, 33.0 / 164.0, 12.0 / 41.0, 0.0, 1.0 }};

//C_rkf78 8阶公式 等于书上P18 表1.3中的 c\hat
const double C_rkf78[13] = { 0.0, 0.0, 0.0, 0.0, 0.0, 34.0 / 105.0, 9.0 / 35.0, 9.0 / 35.0,
							 9.0 / 280.0, 9.0 / 280.0, 0.0, 41.0 / 840.0, 41.0 / 840.0 };  // 8 order parameters

/**
 * @brief GCRS/TT 下的 RKF78 单步积分（NSG2 或 Cartesian）。
 * @param t0_TT 当前单步起始 MJD(TT)，不是整段预报的固定起点。
 * @param X0 输入/输出六分量归一化状态。NSG2 顺序为 [a,xi,eta,h,k,lambda]，
 *           Cartesian 顺序为 [x,y,z,vx,vy,vz]；不进行单位或参考系转换。
 * @param h0_norm 有限、非零的归一化试探步长；正值前向，负值反向。
 * @param h_return 成功时返回实际接受的有符号步长；失败时为 0，X0 保持不变。
 * @note 容差和步长模式取自 thisjob.opn_setting.rkf。默认误差超限时减半重算；
 *       fixed_step=true 时接受给定步长并报告超限警告，仍返回实际误差估计。
 *       自适应试算达到上限、非有限数或步长无法推进 MJD 时，通过 tracer 报错并返回。
 *       调用方须检查 tracer.Code()，应用层自行决定是否 EmergencyExit。
 *       两个单步接口均要求已加载覆盖当前时段的 EOP；物理有效性仍受力模型限制。
 */
// 返回当前接受步的归一化局部误差估计；失败为 NaN，原六参数接口仍可使用。
void RKF78_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob,
	double t0_TT, double* X0, double h0_norm, double& h_return, double& accepted_TE_norm);

void RKF78_NSG2(bastion::Overwatch& tracer,
	const OrbitPropagationJob& thisjob,
	double t0_TT,
	double* X0,
	double h0_norm,
	double& h_return);

// Cartesian 单步接口；输入输出、误差和失败约定同上。
void RKF78_Cart(bastion::Overwatch& tracer,
	const OrbitPropagationJob& thisjob,
	double t0_TT,
	double* X0,
	double h0_norm,
	double& h_return);

// 需要监测积分精度时使用此重载：成功后额外返回接受步的局部截断误差估计。
// 该值是六个归一化分量误差绝对值之和，不是整段轨道误差；失败时返回 NaN。
void RKF78_Cart(bastion::Overwatch& tracer,
	const OrbitPropagationJob& thisjob,
	double t0_TT,
	double* X0,
	double h0_norm,
	double& h_return,
	double& accepted_TE_norm);

/**
 * @brief Advance one EBCRS/TDB normalized Cartesian state with one RKF78 step.
 *
 * @param tracer Runtime event logger and error-state carrier.
 * @param thisjob EBCRS orbit propagation job containing the force model,
 *                integration settings, EOP data, and satellite model.
 * @param t0_TDB Start epoch of the current RKF78 step, expressed as MJD(TDB).
 * @param X0 In/out normalized EBCRS Cartesian state. On input it is the state
 *           at t0_TDB; on output it is advanced by the accepted step.
 * @param h0_norm Initial trial step size in normalized time units.
 * @param h_return 成功时为实际接受的有符号归一化步长；失败时为 0，X0 保持不变。
 *
 * @note 自适应模式误差超限时减半重算；定步长模式保留步长并警告。
 *       达到试算上限、非有限数或 MJD 无法前进时，
 *       通过 tracer 报错返回，调用方须检查错误码，不在此处退出整个进程。
 *       须先加载 JPL 星历与覆盖当前时段的 EOP；步长可正可负。
 * @note This function performs only one RKF78 step. The outer propagation loop
 *       is responsible for accumulating time, clipping the final step, and
 *       recording ephemerides.
 * @note Each RKF78 stage passes the actual intermediate MJD(TDB) and the
 *       corresponding intermediate Cartesian state to Acceleration_Cart_EBCRS.
 */
// 返回当前接受步的归一化局部误差估计；失败为 NaN，原六参数接口仍可使用。
void RKF78_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& thisjob,
	double t0_TDB, double* X0, double h0_norm, double& h_return, double& accepted_TE_norm);

void RKF78_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& thisjob,
	double t0_TDB, double* X0, double h0_norm, double& h_return);

#endif /* CDC_INTEGRATOR_RKF78_HPP */
