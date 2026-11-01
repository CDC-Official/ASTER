//
//  math_basic.hpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#ifndef CDC_MATH_BASIC_HPP
#define CDC_MATH_BASIC_HPP

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <algorithm>


using namespace Eigen;

// n<0: NaN; overflow: +inf.
double Factorial(int n);

double Combinatorial(int n, int m);

//三个旋转矩阵
Matrix3d Rx(double a);

Matrix3d Ry(double a);

Matrix3d Rz(double a);

double Module_3D(double X, double Y, double Z);
//两个向量的夹角

double Angle_Between_Vectors(const VectorXd& vector1, const VectorXd& vector2);

//老版本的多项式拟合，应该不如下面的Polyfit_coeff全面
VectorXd LS_ndPoly_coeff(const int n, ArrayXd x, ArrayXd y);

VectorXd Polyfit_coeff(const int n, ArrayXd& x, ArrayXd& y, double& R_squared);

double Poly_derivatives(VectorXd coeff, double x);

#endif /* CDC_MATH_BASIC_HPP */
