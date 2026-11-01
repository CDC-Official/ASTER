//
//  orbital_operation.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#ifndef CDC_ORBIT_BASIC_HPP
#define CDC_ORBIT_BASIC_HPP

#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <sofa.h>
#include "math_basic.hpp"
#include "utils_constants.hpp"
using namespace Eigen;

namespace determination
{
using StateVector = Eigen::Matrix<double, 6, 1>;
using StateMatrix = Eigen::Matrix<double, 6, 6>;
} // namespace determination

const double DEFAULT_KEPLEREQ_TOL = 1E-18;

double TrueAnomaly2MeanAnomaly(double f, double e);
double TrueAnomaly2EccentricAnomaly(double f, double e);

double EccentricAnomaly2MeanAnomaly(double E, double e);
double EccentricAnomaly2TrueAnomaly(double E, double e);

double MeanAnomaly2EccentricAnomaly(double M, double e, double keplerTOL = DEFAULT_KEPLEREQ_TOL);
double MeanAnomaly2TrueAnomaly(double M, double e);

double SolveKepler(double e, double M, double x0, double TOL = DEFAULT_KEPLEREQ_TOL, int m = 30);  //微分改正解开普勒方程

double SolveGeneralizedKepler(double xi, double eta, double lambda,
	double TOL = DEFAULT_KEPLEREQ_TOL, int m = 30);//简单迭代解广义开普勒方程

void CartesianToKeplerian(double x, double y, double z, double vx, double vy, double vz, double* ele);  //星历根数转换

void KeplerianToCartsian(double a, double e, double i,
	double OMEGA, double omega, double M0, double* vec, double kelperTOL = 1E-18);

#endif /* orbital_transforms_hpp */
