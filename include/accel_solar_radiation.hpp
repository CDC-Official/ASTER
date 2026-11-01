//
// Created by 陈酊存 on 2024/1/23.
//

#ifndef CDC_ACCEL_SOLAR_RADIATION_HPP
#define CDC_ACCEL_SOLAR_RADIATION_HPP

#include <cmath>
#include "model_ephemeris_jpl.hpp"
#include "model_ephemeris_analytic.hpp"
#include "utils_constants.hpp"
#include <Eigen/Dense>

bool in_shadow(const double* r_Norm, const double* r_Sun_Norm);

int Solar_Radition_Acceleration(double* a_Norm,
	const double* r_Norm,
	const double* r_Sun_Norm,
	double AM_Norm,
	double eta);
using namespace Eigen;
#endif //CDC_ACCEL_SOLAR_RADIATION_HPP
