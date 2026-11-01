#include "accel_drag.hpp"
#include "utils_constants.hpp"
#include <Eigen/Dense>

using namespace Eigen;

int AtmosphereDrag_Acceleration(double* a_celestial,
	const double* rv_celestial,
	double CD,
	double AMR_norm,
	double rho_norm)
{
	/*
	 * a_celestial J2000下的大气阻力摄动加速度 归一化
	 * rv_celestial J2000下的位置速度 x y z vx vy vz
	 * MJD_UT1 计算太阳平轨道用，严格来说（似乎）要用UT1，但是大气的不确定因素太大，UTC也完全没问题
	 * CD 阻力系数，推荐为2.2
	 * AMR_norm 归一化面质比
	 * rho_norm 归一化大气密度
	 *
	 */
	double x = rv_celestial[0];
	double y = rv_celestial[1];
	//double z = rv_celestial[2];
	double vx = rv_celestial[3];
	double vy = rv_celestial[4];
	double vz = rv_celestial[5];

	//地球自传速度，可以设置为static避免重复计算
	// FIXME:地球自转速度，要定义好
	// deg/day (TT), Liu's book is consistent with IAU recommended value (not to be confused with per UT1 value)
	const double earthMeanRotationRate = 360.985647365;
	static double omega_a = earthMeanRotationRate * DEGREE_TO_RAD / (1.0 * DAY_TO_NORMALIZED_TIME);

	Vector3d v_vec(vx, vy, vz);
	Vector3d va_vec(-y * omega_a, x * omega_a, 0);//(1.217)
	Vector3d V_vec = v_vec - va_vec;//(1.214)
	double V_module = V_vec.norm();

	//1.213
	Vector3d F_6 = -0.5 * CD * AMR_norm * rho_norm * V_module * V_vec;

	a_celestial[0] = F_6[0];
	a_celestial[1] = F_6[1];
	a_celestial[2] = F_6[2];

	return 1;
}
