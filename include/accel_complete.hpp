//
// Created by 陈酊存 on 2024/10/16.
//

#ifndef CDC_LIB_INCLUDE_ACCEL_COMPLETE_HPP_
#define CDC_LIB_INCLUDE_ACCEL_COMPLETE_HPP_

#include <iostream>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
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
#include "integrator_rkf78.hpp"
#include "orbit_equation.hpp"
#include "utils_overwatch.hpp"

using namespace orbit;

//向前声明OrbitPropagationJob结构体，避免互相包含头文件
struct OrbitPropagationJob;
struct OrbitPropagationJob_EBCRS;

/**
 * @brief Compute the GCRS/TT right-hand side for a normalized nonsingular-2 state.
 *
 * @param tracer Runtime event logger and error-state carrier.
 * @param thisjob Orbit propagation job containing the GCRS force-model switches,
 *                EOP data, atmosphere data, reference density, and satellite model.
 * @param MJD_TT Current right-hand-side evaluation epoch, expressed as MJD(TT).
 * @param n State dimension. The current implementation expects n = 6.
 * @param x Input normalized nonsingular-2 state at MJD_TT.
 * @param f Output normalized time derivative of the nonsingular-2 state.
 *
 * @note The caller, usually an integrator, is responsible for passing the
 *       actual evaluation epoch for the current stage. This function no longer
 *       reconstructs the epoch from an integration-step start time and a stage
 *       offset.
 * @note The perturbing Cartesian acceleration is first computed in the inertial
 *       frame, then projected to RTN/STW and converted to nonsingular-2 element
 *       derivatives through the Gauss planetary equations.
 */
void Acceleration_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob, const double MJD_TT, int n, const double* x, double* f);

/**
 * @brief Compute the GCRS/TT right-hand side for a normalized Cartesian state.
 *
 * @param tracer Runtime event logger and error-state carrier.
 * @param thisjob Orbit propagation job containing the GCRS force-model switches,
 *                EOP data, atmosphere data, reference density, and satellite model.
 * @param MJD_TT Current right-hand-side evaluation epoch, expressed as MJD(TT).
 * @param n State dimension. The current implementation expects n = 6.
 * @param x Input normalized Cartesian state at MJD_TT.
 * @param f Output normalized Cartesian state derivative, [vx, vy, vz, ax, ay, az].
 *
 * @note The state is interpreted in the geocentric inertial frame used by the
 *       existing GCRS/TT propagator. Position and velocity use the project's
 *       normalized Earth units.
 * @note The caller must pass the actual stage epoch. This keeps the dynamics
 *       interface independent of the particular integrator implementation.
 */
void Acceleration_Cart(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob, const double MJD_TT,  int n, const double* x, double* f);

/**
 * @brief Compute the EBCRS post-Newtonian correction from the EIH N-body equations.
 *
 * @param tracer Runtime event logger and error-state carrier.
 * @param major_body_bcrs_state_norm BCRS states of the 11 major bodies in the JPL order
 *                                   used by this project:
 *                                   Mercury, Venus, Earth, Mars, Jupiter, Saturn,
 *                                   Uranus, Neptune, Pluto, Moon, and Sun.
 *                                   Each row is [x, y, z, vx, vy, vz] in the
 *                                   project's normalized Earth units.
 * @param spacecraft_ebcrs_state_norm Spacecraft EBCRS state
 *                                    [x, y, z, vx, vy, vz] in normalized Earth units.
 *                                    EBCRS translates the BCRS origin to the geocenter without rotating its axes; the time argument remains TDB.
 * @param a_post_newtonian_ebcrs_norm Output post-Newtonian acceleration in EBCRS,
 *                                    in normalized acceleration units.
 *
 * @note This function implements the same first post-Newtonian EIH correction
 *       used by the reference Fortran routines FPNEWTON1 and FPNEWTON2 in
 *       EOM_EBCRS.FOR. It returns F_PN = F_PN^S - F_PN^E, where F_PN^S is the
 *       spacecraft BCRS correction and F_PN^E is the Earth BCRS correction.
 * @note The input major-body states must be barycentric, not geocentric. The
 *       spacecraft state is intentionally EBCRS because the function adds the
 *       Earth BCRS state internally to obtain the spacecraft BCRS state.
 * @note The model includes all 11 major bodies listed above. It does not include
 *       the optional Lense-Thirring term or the relativistic J2 term that appear
 *       in the reference Fortran file but are not added in its ordinary orbit
 *       propagation branch.
 */
void PostNewtonian_EIH_EBCRS(bastion::Overwatch& tracer,
	const double major_body_bcrs_state_norm[11][6],
	const double spacecraft_ebcrs_state_norm[6],
	double a_post_newtonian_ebcrs_norm[3]);

/**
 * @brief 计算 E-BCRS/TDB 下直角坐标状态的动力学右函数。
 * @param[in,out] tracer 记录函数调用和模型计算中的错误。
 * @param[in] thisjob 提供力模型开关、EOP、卫星参数和大气密度资料的地月预报任务。
 * @param[in] MJD_TDB 当前积分评估点的 TDB 时刻，单位 MJD。
 * @param[in] n 状态维数，当前使用 6 维状态。
 * @param[in] x 六维 E-BCRS 状态 [x,y,z,vx,vy,vz]，位置尺度 ae，速度尺度 ae/UNITTIME。
 * @param[out] f 六维状态导数，前三项为归一化速度，后三项为归一化加速度，时间尺度 UNITTIME。
 * @details 牛顿质点引力和 EIH 后牛顿项均采用 DE440 的 TDB 相容 GM，以 GE 归一化。
 *          JPL 历表须预先加载。地球姿态由当前 TDB 时刻和位置换算的 TT、UTC 计算；
 *          地球非球形项将 E-BCRS 位置近似旋转到地固轴系，采用 TDB 尺度的 GGM03C 参数。
 */
void Acceleration_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& thisjob,
	double MJD_TDB, int n, const double* x, double* f);

#endif //CDC_LIB_INCLUDE_ACCEL_COMPLETE_HPP_
