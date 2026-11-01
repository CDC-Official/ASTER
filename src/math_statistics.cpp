/**
 * @file math_statistics.cpp
 * @brief 时间序列自相关与 Ljung-Box 检验。
 */
#include "math_statistics.hpp"
#include <boost/math/distributions/chi_squared.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

void autoCorrelation(const Eigen::ArrayXd& series, int k, double& value)
{
	// 计算时间序列在滞后k的自相关系数
	int n = series.size();  // 序列长度
	double mean = series.mean();  // 计算序列均值
	double num = 0.0, denom = 0.0;

	if (std::isnan(mean))
	{
		std::cout << "mean nan" << std::endl;
		std::cout << series;
	}

	// 计算自相关的分子和分母
	for (int i = 0; i < n - k; ++i)
	{
		num += (series(i) - mean) * (series(i + k) - mean);  // 分子
		denom += std::pow(series(i) - mean, 2);  // 分母
	}

	value = num / denom;  // 返回自相关系数
}

void ljungBoxTest(const Eigen::ArrayXd& series, int h, double& Q, double& p)
{
	try {
		// 进行Ljung-Box检验，计算Q值和p值 滞后h
		int n = series.size();
		std::vector<double> acf_values;

		// 计算滞后1到m的自相关系数
		for (int k = 1; k <= h; ++k)
		{
			double tempvalue;
			autoCorrelation(series, k, tempvalue);
			if (std::isnan(tempvalue))
			{
				std::cout << "temp value nan!";
				throw std::runtime_error("NaN in autocorrelation");
			}
			acf_values.push_back(tempvalue);
		}

		// 计算Ljung-Box统计量Q
		double tempQ = 0;
		for (int k = 1; k <= h; ++k)
		{
			tempQ += (acf_values[k - 1] * acf_values[k - 1]) / (n - k);
		}
		tempQ *= n * (n + 2);

		//得到最终的Q和p
		Q = tempQ;

		// 基本检查
		if (h <= 0 || std::isnan(Q) || Q < 0) {
			throw std::runtime_error("Invalid parameters");
		}

		// 计算p值
		boost::math::chi_squared dist(h);  // 卡方分布，自由度为h
		p = 1 - boost::math::cdf(dist, Q);

		// 检查p值是否有效
		if (std::isnan(p) || p < 0 || p > 1) {
			throw std::runtime_error("Invalid p-value");
		}

	} catch (...) {
		// 任何错误都设置为-1
		Q = -1.0;
		p = -1.0;
	}
}
