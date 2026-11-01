//
//  frame_earth.hpp
//  Created by 陈酊存 on 2023/8/30.
//  240924 大规模进行了修改，主要是封装成class，删去了大量底层代码
//  跑了一下IOD没发现问题，之前的代码被备份在老的include文件夹里
//

#ifndef CDC_FRAME_EARTH_HPP
#define CDC_FRAME_EARTH_HPP

#include <Eigen/Dense>
#include "earth_eop.hpp"

using namespace Eigen;
namespace frame { class Geodetic; }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 仍然有大用的非类函数代码（广泛用于类中） @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

// 构造 GCRS 到 ITRS 的位置分量旋转矩阵，不改变长度单位。
// position_itrs = rotation_gcrs_to_itrs * position_gcrs。
// 输入为 MJD(UTC)；EOP 的 dut1、dat 为秒，xp、yp、dx00、dy00 为 rad。
// 速度转换还需考虑矩阵的时间导数，不能只把速度乘以该矩阵。
void GCRS_to_ITRS_Matrix(double MJD_UTC, const aEOP& EOP_needed, Matrix3d& gcrs_to_itrs_emat);

// 上述位置旋转的逆变换，使用矩阵转置。
void ITRS_to_GCRS_Matrix(double MJD_UTC, const aEOP& EOP_needed, Matrix3d& itrs_to_gcrs_emat);

// TOD 使用瞬时真赤道、真春分点；取 SOFA 的 NPB 转置得到 GCRS 分量。
void TOD_to_GCRS_Matrix(double MJD_UTC, double dat, Matrix3d& tod_to_gcrs_emat);

// 教材 O-NWZ 站心地平系（北—西—天）到 GCRS 轴向的旋转，不含原点平移。
// 方位角从北向东增加；经纬度单位 rad，dut1/dat 单位 s。
// 当前用大地经纬度代替天文经纬度，并且本接口未单独施加极移修正。
void Horizontal_to_GCRS_Matrix(double geodetic_lambda_rad, double geodetic_phi_rad,
	double MJD_UTC, double dut1, double dat, Matrix3d& horizontal_to_gcrs_emat);

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 很少用的非类函数代码，但是仍然保留 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void ICRF_difference(const frame::Geodetic& geodetic_class, double* ICRF_dot,
	double MJD_UTC, double dt_sec, const aEOP& EOP_needed);

// 原始二维数组接口，方向及单位约定与 Eigen 矩阵接口一致。
void GCRS_to_ITRS_Matrix(double utc, double dut1, double dat, double xp,
	double yp, double dx00, double dy00, double rc2it[3][3]);  // matrix

#endif // CDC_FRAME_EARTH_HPP
