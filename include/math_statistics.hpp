/**
 * @file math_statistics.hpp
 * @brief 误差统计接口与时间序列相关性检验。
 */
#pragma once
#include <Eigen/Dense>

using namespace Eigen;

void Error_estimation(double* err_x0,
	double* err_ad,
	ArrayXd x0,
	ArrayXd time,
	ArrayXd alpha,
	ArrayXd delta,
	ArrayXd X,
	ArrayXd Y,
	ArrayXd Z,
	double epochofIOD);

void Error_estimation_Cov(double* err_x0,
	double* err_ad,
	double epochofIOD,
	const ArrayXd& x,
	const ArrayXd& time,
	const ArrayXd& alpha,
	const ArrayXd& delta,
	const MatrixXd& R_obs,
	MatrixXd& mat_cov);

void autoCorrelation(const Eigen::ArrayXd& series, int k, double& value);

void ljungBoxTest(const Eigen::ArrayXd& series, int h, double& Q, double& p);
