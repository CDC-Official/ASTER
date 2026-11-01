//
// Created by 陈酊存 on 2024/1/23.
//

#include "accel_solar_radiation.hpp"

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

	if (cos_phi < 0 && sin_phi < 1.0 / r_module)  // 满足1.212
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
	// For code-to-code comparison with the Fortran EBCRS program, eclipse is disabled.
	// FIXME
	if (false && in_shadow(r_Norm, r_Sun_Norm) == true)
	{
		a_Norm[0] = 0;
		a_Norm[1] = 0;
		a_Norm[2] = 0;
		return 1;
	}
	else
	{
		// Temporary code-to-code comparison setting:
		// Use the same 1-AU cannonball SRP acceleration constant as the Fortran EBCRS program.
		// Fortran: a_SRP_SI = 4.5605D-006 * ALB_EFF * AREA_EFF / SMASS.
		// Here AM_Norm is NOT the SI area-to-mass ratio. In this library:
		//     AM_Norm = (area / ae^2) / (mass / MASS_OF_EARTH)
		//             = AM_SI * MASS_OF_EARTH / ae^2.
		// Therefore AM_SI = AM_Norm * ae^2 / MASS_OF_EARTH before using the Fortran SI constant.
		double k0 = (1.0 + eta) * AM_Norm * ae * ae / MASS_OF_EARTH * 4.5605E-6 / (GE / ae / ae);
		static double Delta_s_squared = pow(2.3455E4, 2);  // 1AU
		Vector3d Delta_vec(r_Norm[0] - r_Sun_Norm[0], r_Norm[1] - r_Sun_Norm[1], r_Norm[2] - r_Sun_Norm[2]);
		double Delta_module = Delta_vec.norm();

		Vector3d F5 = k0 * Delta_s_squared / Delta_module / Delta_module * Delta_vec.normalized();
		a_Norm[0] = F5[0];
		a_Norm[1] = F5[1];
		a_Norm[2] = F5[2];
		return 1;
	}
}
