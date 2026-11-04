/**
 * @file orbit_geometry.cpp
 * 2026-10-21：迁入导师 NSG2 几何偏导和六状态解析 STM，改用 Eigen。
 * 数学来源：NAO/gaussLib/orbit_geometry.cpp、orbit_statetransitionmatrix.cpp，Jingshi Tang。
 * 原文件 Copyright (c) 2020 国家航天局空间碎片监测与应用中心，All rights reserved.
 * 保留来源说明；本地科研复刻，不改变原始来源的权利归属。
 */
#include "orbit_geometry.hpp"
#include "math_basic.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

using constants::PI;

EquatorialCoordinates LOS_to_Equatorial(const Eigen::Vector3d& line_of_sight)
{
    const double xy = std::hypot(line_of_sight[0], line_of_sight[1]);
    const double rho = std::hypot(xy, line_of_sight[2]);
    if (rho == 0)
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return {rho, nan, nan};
    }
    double alpha = xy == 0 ? 0 : std::atan2(line_of_sight[1], line_of_sight[0]);
    if (alpha < 0)
        alpha += 2 * PI;
    // 极小负角加一周后可能舍入成 2π，将其归回零以保持半开区间。
    if (alpha >= 2 * PI)
        alpha = 0;
    return {rho, alpha, std::atan2(line_of_sight[2], xy)};
}

void RelativePosition_to_Equatorial(
	double x, double y, double z, double X, double Y, double Z,
	double* equator)  //地心直角矢量转化为站心赤道坐标(rho,alpha,delta)
{
    const auto coordinates = LOS_to_Equatorial(Eigen::Vector3d(x - X, y - Y, z - Z));
    equator[0] = coordinates.rho;
    equator[1] = coordinates.right_ascension;
    equator[2] = coordinates.declination;
}

void RelativePosition_to_Horizontal(Matrix3d& r, double x, double y, double z, double X, double Y, double Z, double* rAh)
{
	// 对应教材式（1.60）的逆过程：先从目标位置减去测站位置。
	// r 为原轴系到站心北—西—天地平系的矩阵；调用方传入正向矩阵的转置。
	// 再按式（1.50）恢复 A、h，其中西向分量带负号，所以方位角用 atan2(-W,N)。
	// 输出距离沿用输入长度单位，A 为 [0,2π) 内的弧度，h 也为弧度。


	double rho, A, h;
	Vector3d L_c, H_c;    //计算得到的单位矢量
	L_c << x - X, y - Y, z - Z;
	rho = L_c.norm();  //先不单位化，求得rho
	L_c.normalize();   //单位化求角量
	H_c = r * L_c;       // 旋转得到北—西—天地平系中的单位视线。
	A = atan2(-H_c[1], H_c[0]);
	if (A < 0)
		A += 2 * PI;  //(0,2pi)
	h = asin(H_c[2]);

	rAh[0] = rho;
	rAh[1] = A;
	rAh[2] = h;
}

void Orbit_Normal(double i, double Omega, Vector3d& vector)
{
	//由i，Omega（rad）计算轨道面单位法向量,输出到一个长度为3的Eigen Vector/Array等等,使用一般数组应该会报错！
	vector[0] = sin(i) * sin(Omega);
	vector[1] = -sin(i) * cos(Omega);
	vector[2] = cos(i);
}

double Visibility_Para_Theta0(double R, double r, double h0)
{
	//返回P17,弧度单位
	//输入为归一化单位、弧度单位
	double theta0 = acos(R / r * cos(h0)) - h0;
	return theta0;
}

bool Visibility(double X, double Y, double Z, double x, double y, double z, double theta0)
{
	//输入的长度单位都为归一化单位，角度单位为弧度
	double R = Module_3D(X, Y, Z);
	double r = Module_3D(x, y, z);
	double rho = Module_3D(X - x, Y - y, Z - z);
	if (R <= 0.0 || r <= 0.0 || !std::isfinite(theta0)) return false;
	const double cosine = (R * R + r * r - rho * rho) / (2 * R * r);
	if (!std::isfinite(cosine)) return false;
	// 天顶附近的舍入可能使余弦略超出 [-1,1]；截回定义域，避免误删可见点。
	double theta = acos(std::clamp(cosine, -1.0, 1.0));
	if (theta < theta0)
		return true;
	else
		return false;
}


void RTN_Matrix_from_Kep(double i, double OMEGA, double omega, double f, Matrix3d& inertial_to_rtn_emat)
{

	// 使用kepler根数计算r t w三个单位矢量组成的矩阵rtnMat，用于将xyz方向的摄动加速度投影到 RTN
	// inertial_to_rtn_emat = [rhat^T that^T what^T]^T
	// 必须传入归一化rad
	// 参考POD ppt part2 P9
	double u = omega + f;
	double cos_OMEGA = cos(OMEGA);
	double sin_OMEGA = sin(OMEGA);
	double cos_i = cos(i);
	double sin_i = sin(i);
	double cos_u = cos(u);
	double sin_u = sin(u);

	Vector3d Phat, Qhat, rhat, that, what;
	Phat << cos_OMEGA, sin_OMEGA, 0;
	Qhat << -sin_OMEGA * cos_i, cos_OMEGA * cos_i, sin_i;

	rhat = cos_u * Phat + sin_u * Qhat;
	that = -sin_u * Phat + cos_u * Qhat;
	what << sin_OMEGA * sin_i, -cos_OMEGA * sin_i, cos_i;// what = rhat x that = Phat x Qhat

	inertial_to_rtn_emat << rhat[0], rhat[1], rhat[2],
		that[0], that[1], that[2],
		what[0], what[1], what[2];
}

void RTN_Matrix_from_Cart(double x, double y, double z, double vx, double vy, double vz, Matrix3d& inertial_to_rtn_emat)
{

	// 使用直角坐标计算r t w三个单位矢量组成的矩阵rtnMat，用于将xyz方向的摄动加速度投影到 RTN
	// inertial_to_rtn_emat = [rhat^T that^T what^T]^T
	// 必须传入归一化norm
	// 参考POD ppt part2 P9
	Vector3d rhat, v, that, what;

	rhat << x, y, z;
	v << vx, vy, vz;

	what = rhat.cross(v);
	what.normalize();
	rhat.normalize();
	that = what.cross(rhat);
	//t本身就是单位化，不过保险期间再归一化一下也无所谓
	that.normalize();

	inertial_to_rtn_emat << rhat[0], rhat[1], rhat[2],
		that[0], that[1], that[2],
		what[0], what[1], what[2];
}

namespace determination
{
StateMatrix GeometryPartialsNSG2(bastion::Overwatch &tracer, const orbit::NonSingular2 &state)
{
    // 同步近点角并求当前 Cartesian；几何偏导只描述同一历元的表示转换，不承担传播。
    orbit::NonSingular2 nsg2 = state;
    nsg2.SyncAnomaly(tracer);
    // 在本地副本中生成状态，不修改调用方传入的根数或近点角。
    orbit::Cartesian cart;
    nsg2.ToCartesian(tracer, cart);
    std::array<double, 6> z;
    nsg2.ToArray(z.data());
    // 根数顺序为 [a,xi,eta,h,k,lambda]，Cartesian 顺序为 [x,y,z,vx,vy,vz]。
    std::array<double, 6> pv;
    cart.ToArray(pv.data());
    const double a = z[0];
    const double xi = z[1];
    // xi、eta 是偏心率矢量分量；h、k 编码轨道平面的倾角与节点方向。
    const double eta = z[2];
    const double h = z[3];
    const double k = z[4];
    // cih=cos(i/2)；接近完全逆行时本组根数的几何偏导不再适用，保持原异常边界。
    const double e2 = xi * xi + eta * eta;
    const double cih = std::sqrt(1 - h * h - k * k);
    if (a <= 0 || e2 >= 1 || !std::isfinite(cih) || cih <= 1e-10)
        throw std::runtime_error("NSG2 geometry outside elliptic non-retrograde domain.");
    // 准备轨道半径和轨道面内角度 u；以下长度、速度和时间全部保持归一化。
    const double r = nsg2.GetValue(tracer, Identifier::VAR_R);
    const double u = nsg2.GetValue(tracer, Identifier::VAR_U_NSG2);
    // 正弦、余弦分别复用，避免在位置与速度各列中重复求值。
    const double su = sin(u);
    const double cu = cos(u);
    // p 为半通径，psr=sqrt(p)；以下三个比值是根数对位置/速度求导时的尺度组合。
    const double p = a * (1 - e2);
    const double psr = sqrt(p);
    const double rp = r / p;
    // ap、ar 分别是 a/p、a/r，仍为无量纲。
    const double ap = a / p;
    const double ar = a / r;
    // xe2 汇集偏心率余量，n 是归一化二体平均角速度；保留原式的求值次序。
    const double xe2 = sqrt(1 - e2) / (1 + sqrt(1 - e2));
    const double n = pow(a, -1.5);
    // pos、vel 分别对应 G 的前三行与后三行，不把米制数据混入此几何矩阵。
    const Eigen::Vector3d pos(pv[0], pv[1], pv[2]);
    const Eigen::Vector3d vel(pv[3], pv[4], pv[5]);
    // 平面根数的平方和交叉项供平面基矢的导数使用。
    const double h2 = h * h;
    const double k2 = k * k;
    const double hk = h * k;
    // Ph/Pk、Qh/Qk 是轨道面基矢对 h/k 的导数，描述平面微转引起的位置与速度变化。
    const Eigen::Vector3d Ph(0, 2 * k, 2 * hk / cih);
    const Eigen::Vector3d Pk(-4 * k, 2 * h, -2 * (1 - h2 - 2 * k2) / cih);
    // Q 基矢的导数与 P 基矢配套，用同一 h/k 约定组装两条平面根数列。
    const Eigen::Vector3d Qh(2 * k, -4 * h, 2 * (1 - 2 * h2 - k2) / cih);
    const Eigen::Vector3d Qk(2 * h, 0, -2 * hk / cih);
    // suet、cuxi 为轨道面速度中的组合；AA..DD 是 xi/eta 对位置的系数。
    const double suet = su + eta;
    const double cuxi = cu + xi;
    // AA、BB 分别乘 pos、vel，形成 xi 列的位置三行。
    const double AA = ap * (-cuxi - rp * suet * (xi * su - eta * cu));
    const double BB = a * r / psr * (su + ar * xe2 * eta + rp * suet);
    // CC、DD 对应 eta 列的位置三行，保留与来源公式一致的正负号。
    const double CC = -ap * (suet - rp * cuxi * (xi * su - eta * cu));
    const double DD = -a * r / psr * (cu + ar * xe2 * xi + rp * cuxi);
    // A2、B2 形成 xi 列的速度三行，单位为归一化速度对无量纲 xi 的导数。
    const double A2 = (suet / p - (ap * suet + ar * su + ar * ar * xe2 * eta) / r) / psr;
    const double B2 = rp * cu + ap * xi;
    // C2、D2 对应 eta 列的速度三行，仍只描述当前历元的几何转换。
    const double C2 = -(cuxi / p - (ap * cuxi + ar * cu + ar * ar * xe2 * xi) / r) / psr;
    const double D2 = rp * su + ap * eta;
    // G 的列顺序为 [a,xi,eta,h,k,lambda]，行顺序为 [x,y,z,vx,vy,vz]，下标从零开始。
    // block<行数,列数>(起始行,起始列) 选取矩阵局部；例如 block<3,1>(3,0) 是 a 列的速度三行。
    StateMatrix G;
    G.block<3, 1>(0, 0) = pos / a;
    G.block<3, 1>(3, 0) = -0.5 * vel / a;
    // xi 列：位置与速度分别使用上面整理的两组系数。
    G.block<3, 1>(0, 1) = AA * pos + BB * vel;
    G.block<3, 1>(3, 1) = A2 * pos + B2 * vel;
    // eta 列同样按位置三行、速度三行成对填写，六行均在当前历元。
    G.block<3, 1>(0, 2) = CC * pos + DD * vel;
    G.block<3, 1>(3, 2) = C2 * pos + D2 * vel;
    // h 列来自轨道面基矢对 h 的变化，不再叠加从拟合历元传播产生的动力学敏感度。
    G.block<3, 1>(0, 3) = r * (cu * Ph + su * Qh);
    G.block<3, 1>(3, 3) = (-suet * Ph + cuxi * Qh) / psr;
    // k 列使用同一半通径和面内角度，仅替换相应的基矢导数。
    G.block<3, 1>(0, 4) = r * (cu * Pk + su * Qk);
    G.block<3, 1>(3, 4) = (-suet * Pk + cuxi * Qk) / psr;
    // lambda 列等价于二体轨道上的微小时移：dr/dlambda=v/n，dv/dlambda=-r/(n*r³)。
    G.block<3, 1>(0, 5) = vel / n;
    G.block<3, 1>(3, 5) = -pos / (n * r * r * r);
    return G;
}
} // namespace determination
