//
// Created by 陈酊存 on 2024/1/23.
//

#include "accel_solar_radiation.hpp"

using namespace constants::Scale;

bool in_shadow(const double* r_Norm, const double* r_Sun_Norm)
{
	// 1.212 不考虑间断问题
	Vector3d r_vec(r_Norm[0], r_Norm[1], r_Norm[2]);
	Vector3d rs_vec(r_Sun_Norm[0], r_Sun_Norm[1], r_Sun_Norm[2]);

	double r_module = r_vec.norm();

	r_vec.normalize();	 // Normalized不改变本身，返回新的，而Normalize改变！
	rs_vec.normalize();	 // 原本是错的后来修改为normalized，后来修改为normalize

	double cos_phi = r_vec.dot(rs_vec);
	double sin_phi = sqrt(1.0 - cos_phi * cos_phi);

	if (cos_phi < 0 && sin_phi < constants::Wgs84::EQUATORIAL_RADIUS / UNIT_LENGTH_M / r_module)  // 满足1.212
	{
		return true;
	}
	else
		return false;
}

int Solar_Radition_Acceleration(double* a_Norm,
								const double* r_Norm,
								const double* r_Sun_Norm,
								double AM_Norm,
								double eta)
{
	/*
	 * 计算光压
	 * eta 反射系数 0～1
	 */
	// 卫星进入地球圆柱阴影后，太阳光压加速度为零。
	if (in_shadow(r_Norm, r_Sun_Norm))
	{
		a_Norm[0] = 0;
		a_Norm[1] = 0;
		a_Norm[2] = 0;
		return 1;
	}
	else
	{
		// 面质比先还原为 m²/kg，辐射压力乘面质比后再换成归一化加速度。
		double k0 = (1.0 + eta) * AM_Norm * UNIT_LENGTH_M * UNIT_LENGTH_M / UNIT_MASS_KG
			* constants::SolarRadiation::REFERENCE_PRESSURE / (UNIT_GM_M3_S2 / UNIT_LENGTH_M / UNIT_LENGTH_M);

		// 参考距离与参考压力配套，按距离平方反比得到卫星处的压力。
		const double Delta_s_squared = pow(constants::SolarRadiation::REFERENCE_DISTANCE / UNIT_LENGTH_M, 2);
		Vector3d Delta_vec(r_Norm[0] - r_Sun_Norm[0], r_Norm[1] - r_Sun_Norm[1], r_Norm[2] - r_Sun_Norm[2]);
		double Delta_module = Delta_vec.norm();

		Vector3d F5 = k0 * Delta_s_squared / Delta_module / Delta_module * Delta_vec.normalized();
		a_Norm[0] = F5[0];
		a_Norm[1] = F5[1];
		a_Norm[2] = F5[2];
		return 1;
	}
}
