#include <stdexcept>
//
//  math_basic.cpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#include "math_basic.hpp"
#include <limits>

// 非负整数阶乘：逐项相乘，负输入仍返回 NaN，超过 double 范围时仍返回正无穷。
double Factorial(int n)
{
    if (n < 0) return std::numeric_limits<double>::quiet_NaN();
    double result = 1.0;
    // 乘下一整数前检查范围，避免先溢出再把无穷结果当作有效阶乘继续使用。
    for (int i = 2; i <= n; ++i)
    {
        if (result > std::numeric_limits<double>::max() / i)
            return std::numeric_limits<double>::infinity();
        // 保留从 2 到 n 的原累乘次序，n=0、1 时循环不执行并返回一。
        result *= i;
    }
    return result;
}

// 广义整数二项式系数：负 n 使用 C(n,m)=(-1)^m C(-n+m-1,m)。
double Combinatorial(int n, int m)
{
    if (m < 0 || (n >= 0 && m > n)) return 0.0;
    if (m == 0) return 1.0;
    // 先用 long long 形成非负上标，避免直接对最小 int 取负；再取较短的对称乘积。
    const long long top = n < 0 ? -static_cast<long long>(n) + m - 1 : n;
    const long long count = std::min(static_cast<long long>(m), top - m);
    double result = 1.0;
    // 按原顺序逐项乘 (top-k+1)/k；出现无穷时保留原来的提前结束行为。
    for (long long k = 1; k <= count; ++k)
    {
        result *= static_cast<double>(top - k + 1) / static_cast<double>(k);
        if (!std::isfinite(result)) break;
    }
    // 仅负上标且 m 为奇数时改变符号，非负上标继续返回通常的组合数。
    return n < 0 && m % 2 != 0 ? -result : result;
}

Matrix3d Rx(double a)
{
	// a rad
	Matrix3d r;
	r << 1, 0, 0,
		0, cos(a), sin(a),
		0, -sin(a), cos(a);
	return r;
}

Matrix3d Ry(double a)
{
	Matrix3d r;
	r << cos(a), 0, -sin(a),
		0, 1, 0,
		sin(a), 0, cos(a);
	return r;
}

Matrix3d Rz(double a)
{
	Matrix3d r;
	r << cos(a), sin(a), 0,
		-sin(a), cos(a), 0,
		0, 0, 1;
	return r;
}

double Module_3D(double X, double Y, double Z)
{
	//返回三维向量的模
	return sqrt(pow(X, 2) + pow(Y, 2) + pow(Z, 2));
}

double Angle_Between_Vectors(const VectorXd& vector1, const VectorXd& vector2)
{
	double dot_product = vector1.normalized().dot(vector2.normalized());
	dot_product = std::max(-1.0, std::min(1.0, dot_product)); // 截断到[-1, 1]，避免向量接近重合时略大于1，导致acos返回nan
	double angle_rad = acos(dot_product);
	return angle_rad;
}

VectorXd LS_ndPoly_coeff(const int n, ArrayXd x, ArrayXd y)
{
	// m组数据的最小二乘多项式拟合，次数为n,n>=1
	// 返回多项式系数，从0次到n次
	//应该不如下面的Polyfit_coeff

	const int m = x.size();
	MatrixXd mat_A(m, n + 1);
	MatrixXd b(m, 1);
	int i, j;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n + 1; j++)
		{
			mat_A(i, j) = pow(x(i), j);
		}
		b(i, 0) = y(i);
	}
	VectorXd coeff = mat_A.colPivHouseholderQr().solve(b);
	return coeff;
}

double Poly_derivatives(VectorXd coeff, double x)
{
	//已知多项式f(x)次数大于等于1,求导数f(x)在x处的值
	const int n = coeff.size() - 1;     // n为多项式次数，为size-1
	double y;
	if (n == 1)
	{
		y = coeff(1);
	}
	else
	{
		y = 0;
		int i;
		for (i = 1; i < n + 1; i++)
		{
			y += coeff(i) * i * pow(x, i - 1);
		}
	}

	return y;
}

VectorXd Polyfit_coeff(const int n, ArrayXd& x, ArrayXd& y, double& R_squared)
{
	// m组数据的最小二乘多项式拟合，次数为n,n>=1
	// 返回多项式系数，从0次到n次
	// 计算决定系数R^2, R_squared, coefficient of determination
	// R^2 = 1 - FVU = 1 - RSS / TSS
	// 已与MATLAB核对，一致

	const int m = x.size();
	MatrixXd mat_A(m, n + 1);
	MatrixXd b(m, 1);
	int i, j;

	for (i = 0; i < m; i++)
	{
		for (j = 0; j < n + 1; j++)
		{
			mat_A(i, j) = pow(x(i), j);
		}
		b(i, 0) = y(i);
	}
	VectorXd coeff = mat_A.colPivHouseholderQr().solve(b);

	double RSS, TSS;
	//计算y-f的二范数（还要再取平方）,f为拟合值

	ArrayXd f(m), temp1(m);     // numerator
	f.setZero();
	for (i = 0; i < n + 1; i++)
	{
		f += coeff[i] * pow(x, i);
	}
	temp1 = y - f;
	// array转vector，array本身没有norm函数
	VectorXd Numerator = VectorXd::Map(&temp1[0], temp1.size());
	RSS = pow(Numerator.norm(), 2);

	//计算y-y_bar的二范数（还要再取平方）,y_bar为原平均值
	ArrayXd temp2(m);
	temp2 = y - y.mean();
	// array转vector，array本身没有norm函数
	VectorXd Denominator = VectorXd::Map(&temp2[0], temp2.size());
	TSS = pow(Denominator.norm(), 2);

	R_squared = 1 - RSS / TSS;

	return coeff;
}
