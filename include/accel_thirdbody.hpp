//
// Created by 陈酊存 on 2024/1/23.
//

#ifndef CDC_ACCEL_THIRDBODY_HPP
#define CDC_ACCEL_THIRDBODY_HPP

#include <iostream>
#include "model_ephemeris_jpl.hpp"
#include "model_ephemeris_analytic.hpp"
#include "utils_constants.hpp"
#include <Eigen/Dense>

using namespace Eigen;

//N代表使用JPL数值历表，A代表使用平均的分析历表
int Solar_Gravity_Acceleration_N(double* a_celestial, const double* r_v_Norm, double MJD_TT);
int Solar_Gravity_Acceleration_A(double* a_celestial, const double* r_v_Norm, double MJD_UT1);
int Lunar_Gravity_Acceleration_N(double* a_celestial, const double* r_v_Norm, double MJD_TT);
int Lunar_Gravity_Acceleration_A(double* a_celestial, const double* r_v_Norm, double MJD_UT1);

int Thirdbody_Gravity_Acceleration(double* a_celestial, const double* r_v_Norm, const double* r_prime_Norm, double m_prime_Norm);

//int Solar_Gravity_Acceleration(double* a_celestial,double
//	int Lunar_Gravity_Acceleration(double* a_celestial,
#endif //CDC_ACCEL_THIRDBODY_HPP
