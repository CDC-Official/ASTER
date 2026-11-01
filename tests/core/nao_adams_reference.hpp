// 测试专用：导师 NAO 的 AdamsPECE 原函数，仅更改函数限定名，保留正文。
// 源文件：NAO/gaussLib/maths_integrator.cpp，582 行开始；不参与生产库构建。
// 用薄适配层将原函数接到相同 NSG2 右函数，以隔离积分公式与力模型差异。
// 原函数片段 SHA256: fda0c8aba2596f87576e427a0e9e92c67b6362bcc0d41145c265f617af95e01e
#pragma once
#include "core_propagation.hpp"
#include "accel_complete.hpp"
namespace nao_reference {
using Detective = bastion::Overwatch;
using d_matrix = std::vector<std::vector<double>>;
constexpr int OPN_INTG_ADAMS_MAXSTAGE = 9;
constexpr int UTILS_EVENT_ERROR = bastion::ERROR, UTILS_EVENT_WARNING = bastion::WARNING;
constexpr int MATHS_FUNDAMENTAL_PARITY_EVEN = 0;
inline int IsParity(int n) { return n % 2; }
inline double BiCoef(int n, int k) {
    double numerator=1,denominator=1;
    if(k>n/2) k=n-k;
    for(int i=n-k+1;i<=n;++i) numerator*=i;
    for(int i=1;i<=k;++i) denominator*=i;
    return numerator/denominator;
}
struct OrbitPropagationNModel {
    const OrbitPropagationJob& job;
    struct Coefficients { bool isFirst=true; int lastK=0; double bex[9]{},bim[9]{}; } cache_AdamsCoef;
};
inline void AccelCompleteNSG2(Detective& tracer, double t0, double t, OrbitPropagationNModel& model,
                            int n, double* state, double* result) {
    Acceleration_NSG2(tracer, model.job, t0+t*NORMALIZED_TIME_TO_DAY,n,state,result);
}
inline void Step(Detective& sherlock,
	int n,
	int k,
	double h,
	double t0,
	OrbitPropagationNModel& opnMod,
	double& t,
	double x[],
	d_matrix& yt,
	bool useCorrection)
{
	// ************************************************************************
	//    Function: AdamsPECE
	// Description: Using Adams method to numerically integrate ODE, with or without correction.
	//    Last Rev: Jingshi @ 2020/01/24
	//
	//       Calls: TODO: acceleration routine
	//   Called By: TODO: propagation routine
	//
	//       Input: n - int, variable dimension
	//              k - int, step number
	//              h - double, step size, in normalize unit
	//              t0 - double, initial epoch, MJD
	//              opnMod - an object of OrbitPropagationNModel, also output, by reference to avoid expensive data copy and is subject to possible modification
	//              t - double, by reference, also output, current epoch, presumably wrt t0 in normalized unit
	//                          it comes in as current epoch and comes out also as current epoch after a step is moved
	//              x - double array of size n, variables to be integrated, correponding to t and is subjected to change
	//      Output: opnMod - an object of OrbitPropagationNModel, also input, by reference to avoid expensive data copy and is subject to possible modification
	//              t - double, by reference, also output, current epoch, presumably wrt t0 in normalized unit
	//                          it comes in as current epoch and comes out also as current epoch after a step is moved
	//              x - double array of size n, variables to be integrated, correponding to t and is subjected to change
	//              yt - d_matrix, with size k*n, acceleration of previous k steps, shifted 1 step forward upon exit, by reference
	//      Return: None
	//
	//        Note: 1. The input/output variable t is passed by reference. It is NOT the current epoch but normalized time since t0.
	//                 The current epoch can be calculated in the form of t0+t1/86400. It is expected that such design helps to preserve more digits.
	//                 Upon exit, t is updated to t+h.
	//              2. Acceleration (YHC) interface:
	//          TODO:  accelCompleteNSG2(double t0, double t1, orbitPropagateNModel opnMod, int n, double x1[n], double y1[n]);
	//              3. This code is an *almost* direct translation from the Fortran library Delta/GIANT.
	//                 Note that the index starts from 1 in Fortran while from 0 in C/C++.
	//              4. The full Adams PECE method uses Bashforth's formula (explicit) as predictor and Moulton's fourmula (implicit) as corrector.
	//              5. In case of PECE, both acceleration calculations employ the complete force models, thus a little bit time-consuming.
	//                 TODO: future plan is to use truncated model (e.g two-body) in PE when useCorrection=true.
	//              6. For yt, the largest index is the latest while the lowest the farthest.
	// ************************************************************************

	sherlock.CheckIn(__FILE__, __func__, __LINE__);

	const double gex[OPN_INTG_ADAMS_MAXSTAGE] =
		{ 1.0, 0.5, 5.0 / 12, 0.375, 251.0 / 720, 95.0 / 288, 19087.0 / 60480, 5257.0 / 17280, 1070017.0 / 3628800 };
	const double gim[OPN_INTG_ADAMS_MAXSTAGE] =
		{ 1.0, -0.5, -1.0 / 12, -1.0 / 24, -19.0 / 720, -3.0 / 160, -863.0 / 60480, -275.0 / 24192,
		  -33953.0 / 3628800 };

	// coefficients of Adams method, calculated once and remains unchanged if k remains unchanged
	double bex[OPN_INTG_ADAMS_MAXSTAGE], bim[OPN_INTG_ADAMS_MAXSTAGE];
	// binomial coefficients, calculated to calculate bex and bim, not needed elsewhere
	double cbi[OPN_INTG_ADAMS_MAXSTAGE][OPN_INTG_ADAMS_MAXSTAGE];

//    double x1[n],y1[n];
	double* x1 = new double[n]{ 0.0 };
	double* y1 = new double[n]{ 0.0 };

	// verify input
	if (k > OPN_INTG_ADAMS_MAXSTAGE)
	{
		sherlock.LogEvent(__FILE__,
			__func__,
			__LINE__,
			"The step number of Adams method is not supported!",
			UTILS_EVENT_ERROR);
		delete[] x1;
		delete[] y1;
		return;
	}

	// copy AdamsCoef from cache
	for (int i = 0; i < OPN_INTG_ADAMS_MAXSTAGE; i++)
	{
		bex[i] = opnMod.cache_AdamsCoef.bex[i];
		bim[i] = opnMod.cache_AdamsCoef.bim[i];
	}

	if (!opnMod.cache_AdamsCoef.isFirst && (k != opnMod.cache_AdamsCoef.lastK))
	{
		sherlock.LogEvent(__FILE__,
			__func__,
			__LINE__,
			"Step number is changed since first run. Coefficients will be recalculated",
			UTILS_EVENT_WARNING);
		opnMod.cache_AdamsCoef.isFirst = true;
	}

	// calculate coefficients if isFirst=true
	if (opnMod.cache_AdamsCoef.isFirst)
	{
		for (int l = 0; l < k; l++)
		{
			bex[l] = 0.0;   // to be accumulated, so initilized here
			bim[l] = 0.0;   // to be accumulated, so initilized here
			for (int m = l; m < k; m++)
			{
				cbi[m][l] = BiCoef(m, l);
				bex[l] = bex[l] + cbi[m][l] * gex[m];
				bim[l] = bim[l] + cbi[m][l] * gim[m];
			}   // end-for of m
			double powl = (IsParity(l) == MATHS_FUNDAMENTAL_PARITY_EVEN ? 1.0 : -1.0);
			bex[l] = powl * bex[l];
			bim[l] = powl * bim[l];
		}       // end-for of l
		opnMod.cache_AdamsCoef.isFirst = false;
		opnMod.cache_AdamsCoef.lastK = k;
		// copy AdamsCoef to cache
		for (int i = 0; i < OPN_INTG_ADAMS_MAXSTAGE; i++)
		{
			opnMod.cache_AdamsCoef.bex[i] = bex[i];
			opnMod.cache_AdamsCoef.bim[i] = bim[i];
		}
	}
	t = t + h;

	// predictor (Bashforth's formula)
	for (int i = 0; i < n; i++)
	{
		x1[i] = x[i];
	}
	for (int l = 0; l <= k - 1; l++)
	{
		for (int i = 0; i < n; i++)
		{
			x1[i] = x1[i] + h * bex[l] * yt[k - l - 1][i];
		}
	}

	if (useCorrection)
	{
		// evaluate y1, if requested
		AccelCompleteNSG2(sherlock, t0, t, opnMod, n, x1, y1);
		if (sherlock.Code() != ERR_NOERROR)
		{
			delete[] x1;
			delete[] y1;
			return;
		}

		// corrector (Moulton's formula)
		// x1 = x
		for (int i = 0; i < n; i++)
		{
			x1[i] = x[i];
		}
		// x1 = x1+h*bim*F
		for (int l = 1; l <= k - 1; l++)
		{
			for (int i = 0; i < n; i++)
			{
//                x1[i] = x1[i]+h*bim[l]*yt[k-l-1][i];      // a possible bug corrected on 2021/09/24.
				x1[i] = x1[i] + h * bim[l] * yt[k - l][i];
			}
		}
		// x1 = x1+h*bim*F1
		for (int i = 0; i < n; i++)
		{
			x1[i] = x1[i] + h * bim[0] * y1[i];
		}
	}   // end-if of useCorrection

	// evaluate and update x and yt
	for (int i = 0; i <= k - 2; i++)
	{
		for (int j = 0; j < n; j++)
		{
			yt[i][j] = yt[i + 1][j];
		}
	}   // end-for of shifting yt

	AccelCompleteNSG2(sherlock, t0, t, opnMod, n, x1, y1);
	if (sherlock.Code() != ERR_NOERROR)
	{
		delete[] x1;
		delete[] y1;
		return;
	}
	for (int j = 0; j < n; j++)
	{
		yt[k - 1][j] = y1[j];
		x[j] = x1[j];
	}

	delete[] x1;
	delete[] y1;
	sherlock.CheckOut(__func__);
}
} // namespace nao_reference
