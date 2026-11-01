//
// Created by 陈酊存 on 2024/1/23.
//

#include "accel_thirdbody.hpp"

//计算日月第三体摄动加速度，使用JPL数值历表，N代表数值
// 下面的四个老了已经不用了
int Solar_Gravity_Acceleration_N(double* a_celestial, const double* r_v_Norm, double MJD_TT)
{
	/*
	 * r_v_Norm 卫星的归一化位置速度，此处只需要位置信息
	 */
	//int nSun = 11;    //目标天体序号 SUN
	//int nMoon = 10; //目标天体序号 Moon
	//int nEarth = 3;    //中心天体序号 EARTH

	//需要在外部事先读JPL
	double JD_TT = MJD_TT + START_OF_MJD;
	double r_v_Sun_Norm[6];
	int a = Get_JPLEPH(JD_TT, 11, 3, r_v_Sun_Norm);

	Vector3d r_prime_vec(r_v_Sun_Norm[0], r_v_Sun_Norm[1], r_v_Sun_Norm[2]);
	Vector3d r_vec(r_v_Norm[0], r_v_Norm[1], r_v_Norm[2]);
	Vector3d Delta_vec = r_vec - r_prime_vec;

	double Delta_module = Delta_vec.norm();

	Vector3d F2 =
		-SUN_EARTH_MASS_RATIO
			* (Delta_vec / Delta_module / Delta_module / Delta_module + r_prime_vec / pow(r_prime_vec.norm(), 3));

	a_celestial[0] = F2[0];
	a_celestial[1] = F2[1];
	a_celestial[2] = F2[2];

	return 1;
}

int Solar_Gravity_Acceleration_A(double* a_celestial, const double* r_v_Norm, double MJD_UT1)
{
	/*
	 * r_v_Norm 卫星的归一化位置速度，此处只需要位置信息
	 */
	//int nSun = 11;    //目标天体序号 SUN
	//int nMoon = 10; //目标天体序号 Moon
	//int nEarth = 3;    //中心天体序号 EARTH

	//需要在外部事先读JPL

	double r_Sun_Norm[3];
	Cart_MeanSun(MJD_UT1, r_Sun_Norm);

	Vector3d r_prime_vec(r_Sun_Norm[0], r_Sun_Norm[1], r_Sun_Norm[2]);
	Vector3d r_vec(r_v_Norm[0], r_v_Norm[1], r_v_Norm[2]);
	Vector3d Delta_vec = r_vec - r_prime_vec;
	double Delta_module = Delta_vec.norm();

	Vector3d F2 =
		-SUN_EARTH_MASS_RATIO
			* (Delta_vec / Delta_module / Delta_module / Delta_module + r_prime_vec / pow(r_prime_vec.norm(), 3));

	a_celestial[0] = F2[0];
	a_celestial[1] = F2[1];
	a_celestial[2] = F2[2];

	return 1;
}

int Lunar_Gravity_Acceleration_N(double* a_celestial, const double* r_v_Norm, double MJD_TT)
{
	/*
	 * r_v_Norm 卫星的归一化位置速度，此处只需要位置信息
	 */
	//int nSun = 11;    //目标天体序号 SUN
	//int nMoon = 10; //目标天体序号 Moon
	//int nEarth = 3;    //中心天体序号 EARTH

	//需要在外部事先读JPL
	double JD_TT = MJD_TT + START_OF_MJD;
	double r_v_Moon_Norm[6];
	int a = Get_JPLEPH(JD_TT, 10, 3, r_v_Moon_Norm);

	Vector3d r_prime_vec(r_v_Moon_Norm[0], r_v_Moon_Norm[1], r_v_Moon_Norm[2]);
	Vector3d r_vec(r_v_Norm[0], r_v_Norm[1], r_v_Norm[2]);
	Vector3d Delta_vec = r_vec - r_prime_vec;

	double Delta_module = Delta_vec.norm();

	Vector3d F2 = -MOON_EARTH_MASS_RATIO
		* (Delta_vec / Delta_module / Delta_module / Delta_module + r_prime_vec / pow(r_prime_vec.norm(), 3));

	a_celestial[0] = F2[0];
	a_celestial[1] = F2[1];
	a_celestial[2] = F2[2];

	return 1;
}

int Lunar_Gravity_Acceleration_A(double* a_celestial, const double* r_v_Norm, double MJD_UT1)
{
	/*
	 * r_v_Norm 卫星的归一化位置速度，此处只需要位置信息
	 */
	//int nSun = 11;    //目标天体序号 SUN
	//int nMoon = 10; //目标天体序号 Moon
	//int nEarth = 3;    //中心天体序号 EARTH

	//需要在外部事先读JPL

	double r_Moon_Norm[3];
	Cart_MeanMoon(MJD_UT1, r_Moon_Norm);

	Vector3d r_prime_vec(r_Moon_Norm[0], r_Moon_Norm[1], r_Moon_Norm[2]);
	Vector3d r_vec(r_v_Norm[0], r_v_Norm[1], r_v_Norm[2]);
	Vector3d Delta_vec = r_vec - r_prime_vec;

	double Delta_module = Delta_vec.norm();

	Vector3d F2 = -MOON_EARTH_MASS_RATIO
		* (Delta_vec / Delta_module / Delta_module / Delta_module + r_prime_vec / pow(r_prime_vec.norm(), 3));

	a_celestial[0] = F2[0];
	a_celestial[1] = F2[1];
	a_celestial[2] = F2[2];

	return 1;
}

//优化后的第三体摄动加速度
int Thirdbody_Gravity_Acceleration(double* a_celestial, const double* r_v_Norm, const double* r_prime_Norm, double m_prime_Norm)
{
	/*
	 * 输出卫星、第三体归一化位置、第三体归一化质量，计算第三体摄动加速度
	 */
	//第三体位置
	Vector3d r_prime_vec(r_prime_Norm[0], r_prime_Norm[1], r_prime_Norm[2]);
	//卫星位置
	Vector3d r_vec(r_v_Norm[0], r_v_Norm[1], r_v_Norm[2]);
	//相对的Delta
	Vector3d Delta_vec = r_vec - r_prime_vec;

	double Delta_module = Delta_vec.norm();

	Vector3d F2 = -m_prime_Norm
		* (Delta_vec / Delta_module / Delta_module / Delta_module + r_prime_vec / pow(r_prime_vec.norm(), 3));

	a_celestial[0] = F2[0];
	a_celestial[1] = F2[1];
	a_celestial[2] = F2[2];
	return 1;
}
