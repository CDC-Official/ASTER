/**
 * @file orbit_geometry.hpp
 * @brief 视线、观测几何、轨道方向及瞬时 NSG2 几何偏导。
 */
#pragma once
#include "orbit_statedef.hpp"

// 同一赤道直角坐标系中的相对视线（目标位置减观测者位置）；距离沿用输入长度单位。
struct EquatorialCoordinates
{
    double rho; // 距离标量，单位与输入相对位置向量一致。
    double right_ascension; // rad，范围 [0, 2π)。
    double declination;    // rad，范围 [-π/2, π/2]。
};
// 输入应为有限向量。零视线返回零距离及 NaN 角度；极点约定 RA=0，
// 但此时 RA 无几何意义，调用方仍须自行判断测角偏导等计算是否适用。
EquatorialCoordinates LOS_to_Equatorial(const Eigen::Vector3d& line_of_sight);

// 目标与观测者须处于同一赤道坐标系，输出依次为距离、RA(rad)、Dec(rad)。
void RelativePosition_to_Equatorial(double x, double y, double z, double X, double Y, double Z, double* equator);
// 目标与测站须同系、同历元、同长度单位；先作差得到站心视线。
// 输入矩阵 r 将原轴系分量转为教材 O-NWZ（北—西—天）分量。
// 输出依次为距离、方位角、高度角；距离沿用输入单位，角度为 rad。
void RelativePosition_to_Horizontal(Matrix3d& r, double x, double y, double z, double X, double Y, double Z, double* rAh);

void Orbit_Normal(double i, double Omega, Vector3d& vector);

double Visibility_Para_Theta0(double R, double r, double h0);//可见性参数\theta_0

bool Visibility(double X, double Y, double Z, double x, double y, double z, double theta0);//判断可见性

// RTN 矩阵将输入惯性轴系中的分量投影到径向 R、横向 T、法向 N。
// R 沿位置，N 沿轨道角动量，T=N×R；三行分别为三个单位向量的转置。
// T 一般不与速度方向重合。根数角度为 rad，Cartesian 位置速度须使用一致单位。
void RTN_Matrix_from_Kep(double i, double OMEGA, double omega, double f, Matrix3d& inertial_to_rtn_emat);

void RTN_Matrix_from_Cart(double x, double y, double z, double vx, double vy, double vz, Matrix3d& inertial_to_rtn_emat);

namespace determination
{
// 当前历元的几何雅可比 G=d[x,y,z,vx,vy,vz]/d[a,xi,eta,h,k,lambda]。
// 两边均归一化；只反映瞬时根数转换，不包含从拟合历元到观测历元的动力学。
StateMatrix GeometryPartialsNSG2(bastion::Overwatch &tracer, const orbit::NonSingular2 &state);
} // namespace determination
