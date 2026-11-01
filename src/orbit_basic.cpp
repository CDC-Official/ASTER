//
//  orbital_operation.cpp
//  Test_dylib
//  一些基本的轨道力学变换
//  Created by 陈酊存 on 2023/8/30.
//

#include "orbit_basic.hpp"
#include <algorithm>
#include <limits>

double MeanAnomaly2EccentricAnomaly(double M, double e, double keplerTOL)
{
	//平近点角->偏近点角，求解开普勒方程，keplerTOL默认 1E-16 最大迭代次数20
	//返回的E 在0-2pi
	double E;
	E = SolveKepler(e, M, M, keplerTOL, 20);
	return E;
}

double EccentricAnomaly2MeanAnomaly(double E, double e)
{
	//由偏近点角->平近点角 rad 仅适用于椭圆
	//M在0-2pi
	double M = E - e * sin(E);
	M = fmod(M, 2 * pi);
	if (M >= 0)
		return M;  //换为(0, 2pi)
	else
		return M + 2 * pi;
}

double EccentricAnomaly2TrueAnomaly(double E, double e)
{
	//E->f
	double f = atan2(sqrt(1 - e * e) * sin(E), cos(E) - e);
	if (f < 0)
	{
		f += 2 * pi;
	}
	return f;
}

double TrueAnomaly2MeanAnomaly(double f, double e)
{
	//f->M
	double E = TrueAnomaly2EccentricAnomaly(f, e);
	double M = EccentricAnomaly2MeanAnomaly(E, e);
	return M;
}

double TrueAnomaly2EccentricAnomaly(double f, double e)
{
	//f->E
	double cos_f = cos(f);
	double sin_f = sin(f);
	double r_divide_a = (1.0 - e * e) / (1.0 + e * cos_f);
	double cos_E = r_divide_a * cos_f + e;
	double sin_E = r_divide_a / sqrt(1 - e * e) * sin_f;
	double E = atan2(sin_E, cos_E);
	if (E < 0)
	{
		E += 2 * pi;
	}
	return E;
}

double MeanAnomaly2TrueAnomaly(double M, double e)
{
	// M->f
	double E = MeanAnomaly2EccentricAnomaly(M, e);
	double f = EccentricAnomaly2TrueAnomaly(E, e);
	return f;
}

double SolveKepler(double e, double M, double x0, double TOL,
	int m)  //离心率e, 平近点角M, 初值x0, TOL, 最大迭代次数
{
	double p0 = x0;
	double p;
	int count1 = 0;
	while (true)
	{
		count1 += 1;
		p = p0 - (p0 - e * sin(p0) - M) / (1 - e * cos(p0));
		if (count1 >= m || abs(p - p0) <= TOL)
			break;
		p0 = p;
	}

	if (p >= 0)
		return p;  //换为(0, 2pi)

	else
		return p + 2 * pi;
}

double SolveGeneralizedKepler(double xi, double eta, double lambda, double TOL, int m)  //离心率e, 平近点角M, 初值x0, TOL, 最大迭代次数
{
	//解广义Kepler方程，《轨道力学基础》P88 式(3.101)
	int count = 0;
	double uTilde, uTilde0;
	if (abs(xi) >= 1.0 || abs(eta) >= 1.0 || xi * xi + eta * eta >= 1.0)
	{
		std::cerr << "Invalid non-singular element xi and/or eta!" << std::endl;
		return 0.0;
	}
	//初值为lambda
	uTilde = lambda;
	// Simple iteration
	while (true)
	{
		count++;
		uTilde0 = uTilde;
		uTilde = lambda + xi * sin(uTilde0) - eta * cos(uTilde0);
		if (count >= m || abs(uTilde - uTilde0) <= TOL)
			break;
	}

	// 2026-10-21：高偏心率时简单迭代可能用完 m 轮仍未收敛，直接返回会使
	// 根数转换及定轨偏导不一致。仅在方程残差仍超限时，补用带区间保护的 Newton 法。
	// F(u)=u-xi*sin(u)+eta*cos(u)-lambda，F'(u)>=1-e>0；根必在 [lambda-e,lambda+e]。
	// 保留原本已收敛的计算路径；阈值同时考虑 double 在当前角度量级的分辨率。
	const double e = std::hypot(xi, eta);
	const double tolerance = std::max(TOL, 8 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(lambda)));
	double lo = lambda - e, hi = lambda + e;
	for (int extra = 0; extra < 80; ++extra)
	{
		const double f = uTilde - xi * sin(uTilde) + eta * cos(uTilde) - lambda;
		if (std::abs(f) <= tolerance) break;
		if (f > 0) hi = uTilde; else lo = uTilde;
		const double next = uTilde - f / (1 - xi * cos(uTilde) - eta * sin(uTilde));
		uTilde = (next > lo && next < hi) ? next : (lo + hi) / 2;
	}

	if (uTilde >= 0)
		return uTilde;  //换为(0, 2pi)

	else
		return uTilde + 2 * pi;

}

void KeplerianToCartsian(double a, double e, double i, double OMEGA, double omega,
	double M0, double* vec, double kelperTOL)
{
	double E;  //偏近点角
	Array3d P, Q, r, v;
	double r_mod;
	E = SolveKepler(e, M0, M0, kelperTOL, 10);
	P << cos(OMEGA) * cos(omega) - sin(OMEGA) * sin(omega) * cos(i),
		sin(OMEGA) * cos(omega) + cos(OMEGA) * sin(omega) * cos(i),
		sin(omega) * sin(i);
	Q << -cos(OMEGA) * sin(omega) - sin(OMEGA) * cos(omega) * cos(i),
		-sin(OMEGA) * sin(omega) + cos(OMEGA) * cos(omega) * cos(i),
		cos(omega) * sin(i);
	r = a * (cos(E) - e) * P + a * sqrt(1 - e * e) * sin(E) * Q;

	r_mod = sqrt(r(0) * r(0) + r(1) * r(1) + r(2) * r(2));    // r模长

	v = sqrt(a) / r_mod * (-sin(E) * P + sqrt(1 - e * e) * cos(E) * Q);

	vec[0] = r(0);
	vec[1] = r(1);
	vec[2] = r(2);
	vec[3] = v(0);
	vec[4] = v(1);
	vec[5] = v(2);
}

void CartesianToKeplerian(double x, double y, double z, double vx, double vy,
	double vz, double* ele)  //由位置矢量、速度矢量计算轨道根数
{
	//根数a, e, i, OMEGA(Ω),omega(ω),M0;
	//不再采用static方法，因为只有一个地址，在此调用会改变值，之前的指针仍然指向此处，容易出问题
	// static double elements[6];
	// //声明变量不可同时赋值，因为程序会忽略static变量的初始化，会不变
	double r0, v0, r0v0, a, e, i, OMEGA, omega, E0, M0, Pz, Qz, k, Rx, Ry, Rz;
	r0 = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
	v0 = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
	r0v0 = x * vx + y * vy + z * vz;

	a = 1.0 / (2.0 / r0 - pow(v0, 2) / 1.0);    //归一化，μ = GM = GE = 1

	e = sqrt(pow(1.0 - r0 / a, 2) + pow(r0v0, 2) / (1 * a));

	E0 = atan2(r0v0 / sqrt(1.0 * a), 1 - r0 / a);     //(-pi, pi)换为(0, 2pi)
	if (E0 < 0)
		E0 += 2 * pi;

	M0 = E0 - e * sin(E0);
	if (M0 < 0)
		M0 += 2 * pi;

	Pz = cos(E0) / r0 * z - sqrt(a / 1) * sin(E0) * vz;
	Qz = (sin(E0) / r0 * z + sqrt(a / 1) * (cos(E0) - e) * vz) /
		sqrt(1 - pow(e, 2));

	k = sqrt(1.0 / (1 * a * (1 - pow(e, 2))));
	Rx = k * (y * vz - z * vy);
	Ry = k * (z * vx - x * vz);
	Rz = k * (x * vy - y * vx);

	omega = atan2(Pz, Qz);
	if (omega < 0)
		omega += 2 * pi;

	OMEGA = atan2(Rx, -Ry);
	if (OMEGA < 0)
		OMEGA += 2 * pi;
	i = acos(Rz);  //(0, pi)
	ele[0] = a;
	ele[1] = e;
	ele[2] = i;
	ele[3] = OMEGA;
	ele[4] = omega;
	ele[5] = M0;
}
