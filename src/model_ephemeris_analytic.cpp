//
//  ephemeris.cpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/28.
//

#include "model_ephemeris_analytic.hpp"

//用于积分器的老版本，使用了static类型的MJD0，定义在RKF78.hpp里
void MeanMoonCart(double MJD0,
	double RKF_t,
	double* R_Moon)  // XYZ of the Mean Moon in Earth J2000, using normalized unit
{
	/*
	MJD0为积分的起点儒略日，通常设置为一个全局的static变量，RKFt为积分中的归一化时间，从0开始，因此积分过程中的MJD实际为MJD0+RKFt（单位转化为天后）
	*/
	double T, d;  // Julian centuries and days from J2000.0
	double JD = MJD0 + RKF_t * UNITTIME / 86400.0 + 2400000.5;
	d = JD - 2451545.0;
	T = d / 36525.0;

	double a, e, i, OMEGA, omega, M;  // mean elements in ecliptic coordinate, normalized & rad
	a = 384747981 / ae;
	e = 0.05488;
	i = 5.1298 / 180. * pi;
	OMEGA = (125.0446 - 1934.14 * T) / 180. * pi;
	omega = (318.3087 + 6003.15 * T) / 180. * pi;
	M = (134.9634 + 13.065 * d) / 180. * pi;

	double eps = 23.4393 / 180 * pi;
	//求出月球在平黄道坐标系下的位置，在旋转Rx（-eps）转换一下就行

	double state_ecliptic_norm_array[6];
	KeplerianToCartsian(a, e, i, OMEGA, omega, M, state_ecliptic_norm_array);
	double x = state_ecliptic_norm_array[0];
	double y = state_ecliptic_norm_array[1];
	double z = state_ecliptic_norm_array[2];  // 3d position in ecliptic coordinate

	R_Moon[0] = x;
	R_Moon[1] = cos(-eps) * y + sin(-eps) * z;
	R_Moon[2] = -sin(-eps) * y + cos(-eps) * z;
}

void MeanSunCart(double MJD0,
	double RKF_t,
	double* R_Sun)    // XYZ of the Mean Sun in Earth J2000, using normalized unit
{
	/*
	MJD0为积分的起点儒略日，通常设置为一个全局的const变量，RKFt为积分中的归一化时间，从0开始，因此积分过程中的MJD实际为MJD0+RKFt（单位转化为天后）
	*/
	double T, d;  // Julian centuries and days from J2000.0
	double JD = MJD0 + RKF_t * UNITTIME / 86400.0 + 2400000.5;
	d = JD - 2451545.0;
	T = d / 36525.0;

	double a, e, i, OMEGA, omega, M;  // mean elements in equator coordinate, normalized & rad
	a = 1.00000102 * 1.4959787E11 / ae;
	e = 0.016709;
	i = 23.4393 / 180 * pi;
	OMEGA = 0.;
	omega = (282.9373 + 0.32 * T) / 180. * pi;
	M = (357.5291 + 0.9856 * d) / 180. * pi;

	double vector[6];
	KeplerianToCartsian(a, e, i, OMEGA, omega, M, vector);

	R_Sun[0] = vector[0];
	R_Sun[1] = vector[1];
	R_Sun[2] = vector[2];
}

//新版本，直接返回MJD对应的平位置
void Cart_MeanMoon(double MJD_UT1, double* R_Moon, double* V_Moon)
{
	/// XYZ of the Mean Moon in Earth J2000, using normalized unit
	/// - Parameters:
	///   - MJD_UT1: MJD_UT1 of a certain moment
	///   - Cart: XYZ of the Mean Moon in Earth J2000, using normalized unit
	//	231227 跟精密定轨part2 P94一样 平均轨道1
	double T, d;  // Julian centuries and days from J2000.0
	double JD = MJD_UT1 + 2400000.5;// Julian day
	d = JD - 2451545.0;
	T = d / 36525.0;

	double a, e, i, OMEGA, omega, M;  // mean elements in ecliptic coordinate, normalized & rad
	a = 384747981. / ae;
	e = 0.05488;
	i = 5.1298 / 180. * pi;
	OMEGA = (125.0446 - 1934.14 * T) / 180. * pi;
	omega = (318.3087 + 6003.15 * T) / 180. * pi;
	M = (134.9634 + 13.065 * d) / 180. * pi;

	double eps = 23.4393 / 180 * pi;
	//求出月球在平黄道坐标系下的位置，在旋转Rx（-eps）转换一下就行

	double state_ecliptic_norm_array[6];
	KeplerianToCartsian(a, e, i, OMEGA, omega, M, state_ecliptic_norm_array);
	double x = state_ecliptic_norm_array[0];
	double y = state_ecliptic_norm_array[1];
	double z = state_ecliptic_norm_array[2];  // 3d position in ecliptic coordinate


	R_Moon[0] = x;
	R_Moon[1] = cos(-eps) * y + sin(-eps) * z;
	R_Moon[2] = -sin(-eps) * y + cos(-eps) * z;

	// 原转换已在 state_ecliptic_norm_array 后三槽给出 mu=1 的归一化辅助速度，旧调用仍只取位置。
	// 这不是随历元变化的解析历表导数，也不能作为物理月球速度使用。
	if (V_Moon != nullptr)
	{
		// 辅助速度与位置使用同一个 Rx(-eps)，使两者共同给出赤道系的轨道平面。
		double vx = state_ecliptic_norm_array[3];
		double vy = state_ecliptic_norm_array[4];
		double vz = state_ecliptic_norm_array[5];

		// 输出以 ae/UNITTIME 归一化，仅供同一 mu=1 约定下恢复原根数的方向。
		V_Moon[0] = vx;
		V_Moon[1] = cos(-eps) * vy + sin(-eps) * vz;
		V_Moon[2] = -sin(-eps) * vy + cos(-eps) * vz;
	}

}

void Cart_MeanSun(double MJD_UT1, double* R_Sun, double* V_Sun)
{
	/// XYZ of the Mean Sun in Earth J2000, using normalized unit
	/// - Parameters:
	///   - MJD_UT1: MJD_UT1 of a certain moment
	///   - Cart: XYZ of the Mean Moon in Earth J2000, using normalized unit
	//	231227 跟精密定轨part2 P94一样 平均轨道1
	double T, d;  // Julian centuries and days from J2000.0
	double JD = MJD_UT1 + 2400000.5;
	d = JD - 2451545.0;
	T = d / 36525.0;

	double a, e, i, OMEGA, omega, M;  // mean elements in equator coordinate, normalized & rad
	a = 1.00000102 * 1.4959787E11 / ae;
	e = 0.016709;
	i = 23.4393 / 180 * pi;
	OMEGA = 0.;
	omega = (282.9373 + 0.32 * T) / 180. * pi;
	M = (357.5291 + 0.9856 * d) / 180. * pi;

	double vector[6];
	KeplerianToCartsian(a, e, i, OMEGA, omega, M, vector);

	R_Sun[0] = vector[0];
	R_Sun[1] = vector[1];
	R_Sun[2] = vector[2];

	// vector 后三槽是原六根数在 mu=1 下转换的归一化辅助速度，不增加另一套日月模型。
	// 它不是解析历表随历元的真实导数或物理太阳速度；旧二参调用不写这些分量。
	if (V_Sun != nullptr)
	{
		// 原根数已在赤道系中；直接输出，以 ae/UNITTIME 归一化，仅供恢复轨道方向。
		V_Sun[0] = vector[3];
		V_Sun[1] = vector[4];
		V_Sun[2] = vector[5];
	}
}

void Cart_OsculatingMoon(double MJD_UT1, double* R_Moon)
{
	/// XYZ of the Osculating Moon in Earth J2000, using normalized unit
	/// - Parameters:
	///   - MJD_UT1: MJD_UT1 of a certain moment
	///   - Cart: XYZ of the Mean Moon in Earth J2000, using normalized unit
	//	231227 跟精密定轨part2 P94一样 删去了平均项
	double d;  // Julian centuries and days from J2000.0
	double JD = MJD_UT1 + 2400000.5;// Julian day
	d = JD - 2451545.0;

	double a, e, i, OMEGA, omega, M;  // mean elements in ecliptic coordinate, normalized & rad
	a = 384747981. / ae;
	e = 0.05488;
	i = 5.1298 / 180. * pi;
	OMEGA = (125.0446) / 180. * pi;
	omega = (318.3087) / 180. * pi;
	M = (134.9634 + 13.065 * d) / 180. * pi;

	double eps = 23.4393 / 180 * pi;
	//求出月球在平黄道坐标系下的位置，在旋转Rx（-eps）转换一下就行

	double state_ecliptic_norm_array[6];
	KeplerianToCartsian(a, e, i, OMEGA, omega, M, state_ecliptic_norm_array);
	double x = state_ecliptic_norm_array[0];
	double y = state_ecliptic_norm_array[1];
	double z = state_ecliptic_norm_array[2];  // 3d position in ecliptic coordinate


	R_Moon[0] = x;
	R_Moon[1] = cos(-eps) * y + sin(-eps) * z;
	R_Moon[2] = -sin(-eps) * y + cos(-eps) * z;

}

void Cart_OsculatingSun(double MJD_UT1, double* R_Sun)
{
	//不变椭圆，只有M有变化
	/// XYZ of the Osculating Sun in Earth J2000, using normalized unit
	/// - Parameters:
	///   - MJD_UT1: MJD_UT1 of a certain moment
	///   - Cart: XYZ of the Mean Moon in Earth J2000, using normalized unit
	//	231227 跟精密定轨part2 P94一样 忽略了omega的平均运动
	double d;  // Julian centuries and days from J2000.0
	double JD = MJD_UT1 + 2400000.5;
	d = JD - 2451545.0;

	double a, e, i, OMEGA, omega, M;  // mean elements in equator coordinate, normalized & rad
	a = 1.00000102 * 1.4959787E11 / ae;
	e = 0.016709;
	i = 23.4393 / 180 * pi;
	OMEGA = 0.;
	omega = (282.9373) / 180. * pi;//不变
	M = (357.5291 + 0.9856 * d) / 180. * pi;//不变

	double vector[6];
	KeplerianToCartsian(a, e, i, OMEGA, omega, M, vector);

	R_Sun[0] = vector[0];
	R_Sun[1] = vector[1];
	R_Sun[2] = vector[2];
}
