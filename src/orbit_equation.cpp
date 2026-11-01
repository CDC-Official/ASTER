//
// Created by 陈酊存 on 2024/10/16.
//

#include "orbit_equation.hpp"

void GaussPlanetaryEquation_RTNtoNSG2(bastion::Overwatch& tracer,
	const NonSingular2& nsg2,
	const double* RTN,
	double* f)
{
	// ************************************************************************
	//    Function: GaussPlanetaryEquation_RTNtoNSG2
	// Description: Gauss planetary equation - decompose RTN acceleration to NSG2 elements.
	//    Last Rev: Jingshi @ 2020/01/27
	//
	//       Calls: NonSingular2.Get, sin, cos, pow, sqrt
	//   Called By: When Gauss planetary equation is needed, e.g, in AccelComplete
	//
	//       Input: nsg2 - NonSingular2 object, const and by reference
	//              RTN - double[3], const, RTN acceleration
	//      Output: dxdt - double[6], NSG2 acceleration
	//      Return: None
	//
	//        Note: None yet.
	// ************************************************************************
	// 直接借用老师的代码

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double a = nsg2.GetValue(tracer, VAR_SEMIMAJORAXIS);
	double xi = nsg2.GetValue(tracer, VAR_XI_NSG2);
	double eta = nsg2.GetValue(tracer, VAR_ETA_NSG2);
	double h = nsg2.GetValue(tracer, VAR_H_NSG2);
	double k = nsg2.GetValue(tracer, VAR_K_NSG2);
	double u = nsg2.GetValue(tracer, VAR_U_NSG2);
	double uTilde = nsg2.GetValue(tracer, VAR_UTILDE_NSG2);

	double n = pow(a, -1.5);
	double e2 = xi * xi + eta * eta;
	double je2 = 1.0 - e2;
	double xe2 = sqrt(je2);
	double xje2 = xe2 + je2;              // xje2 = sqrt(1-e2)+(1-e2) = sqrt(1-e2)*(1+sqrt(1-e2))
	double cih = sqrt(1.0 - h * h - k * k);    // cos(i/2) = sqrt(1-h*h-k*k)
	double ra = nsg2.GetValue(tracer, VAR_R) / a;
	double rp = ra / je2;
	double su = sin(u);
	double cu = cos(u);
	double esf = xi * su - eta * cu;       // e*sin(f) = xi*sin(u)-eta*cos(u)
	double ecf = xi * cu + eta * su;       // e*cos(f) = xi*cos(u)+eta*sin(u)
	double ese = xi * sin(uTilde) - eta * cos(uTilde);   // e*sin(E) = xi*sin(u_tilde)-eta*cos(u_tilde)
	double tihswf = (h * su - k * cu) / cih;  // tan(i/2)*sin(w+f) = [h*sin(u)-k*cos(u)]/cos(i/2)
	double sihcwf = k * su + h * cu;        // sin(i/2)*cos(w+f) = k*sin(u)+h*cos(u)

	double RR = RTN[0];
	double TT = RTN[1];
	double NN = RTN[2];

	double xena = xe2 / n / a;                   // sqrt(1-e2)/(na)
	double Wnaeir = NN / (2.0 * n * a * xe2 * cih) * ra;   // N/[2na*sqrt(1-e2)*cos(i/2)]*(r/a)
	double esexje2 = ese / xje2;                  // [xi*sin(u_tilde)-eta*cos(u_tilde)]/sqrt(1-e2)/(1+sqrt(1-e2))
	double Nrpiwf = NN * rp * tihswf;              // N*(r/p)*(h*sin(u)-k*cos(u))/cos(i/2)

	f[0] = 2.0 / n / xe2 * (RR * esf + TT / rp);
	f[1] = xena * (RR * su + TT * (cos(uTilde) + cu - eta * esexje2) + Nrpiwf * (-eta));
	f[2] = xena * (-RR * cu + TT * (sin(uTilde) + su + xi * esexje2) + Nrpiwf * xi);
	f[3] = Wnaeir * (cu - h * sihcwf);
	f[4] = Wnaeir * (su - k * sihcwf);
	f[5] = n - xena * (RR * 2.0 * xe2 * rp + (RR * ecf - TT * (1.0 + rp) * esf) / (1.0 + xe2) - Nrpiwf);

	tracer.CheckOut(__func__);
}
