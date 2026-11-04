//
//  frame_earth.cpp
//  Created by 陈酊存 on 2023/8/30.
//  240924 大规模进行了修改，主要是封装成class，删去了大量底层代码，之前的代码被备份在老的include文件夹里
//

#include "frame_earth.hpp"

#include "utils_constants.hpp"
#include "math_basic.hpp"
#include <sofa.h>

using constants::PI;

void GCRS_to_ITRS_Matrix(double MJD_UTC, const aEOP& EOP_needed, Matrix3d& gcrs_to_itrs_emat)
{
	//
	// utc    MJD     MJD     for a certain  moment
	// dut1   second  UT1-UTC
	// dat    second  bulletin C, TAI-UTC
	// xp,yp,dx00,dy00 rad
	// This function modify rc2ti[3][3],the matrix rorating GCRS to ITRS
	// Polar motion included


	double rc2it[3][3];
	double tai = MJD_UTC + EOP_needed.dat / constants::SECONDS_PER_DAY;     // MJD
	double tt = tai + 32.184 / constants::SECONDS_PER_DAY;     // MJD
	double ut1 = MJD_UTC + EOP_needed.dut1 / constants::SECONDS_PER_DAY;     // MJD

	// IAU2000B method
	double x, y, s;
	double rc2i[3][3], rc2ti[3][3], rpom[3][3];

	// CIP and CIO, *******************************IAU2000B**
	iauXys00b(constants::MJD_OFFSET, tt, &x, &y, &s);

	// add CIP corrections
	x += EOP_needed.dx00;
	y += EOP_needed.dy00;

	// GCRS to CIRS matrix
	iauC2ixys(x, y, s, rc2i);

	// ERA
	double era = iauEra00(constants::MJD_OFFSET, ut1);

	// Form celestial-terrestrial matrix (no polar motion yet)
	iauCr(rc2i, rc2ti);
	iauRz(era, rc2ti);

	// Polar motion matrix (TIRS->ITRS,IERS2003)
	iauPom00(EOP_needed.xp, EOP_needed.yp, iauSp00(constants::MJD_OFFSET, tt), rpom);

	// Form celestial-terrestrial matrix(including polar motion)
	iauRxr(rpom, rc2ti, rc2it);
	gcrs_to_itrs_emat << rc2it[0][0], rc2it[0][1], rc2it[0][2],
		rc2it[1][0], rc2it[1][1], rc2it[1][2],
		rc2it[2][0], rc2it[2][1], rc2it[2][2];
}


// 假设你已经包含了 SOFA 头文件和 Eigen 头文件
// #include <sofa.h>
// #include <Eigen/Dense>

/*
void ICRF_to_ITRF_Matrix(double MJD_UTC, const aEOP& EOP_needed, Matrix3d& C2T_mat)
{
	//
	 //utc    MJD
	// dut1   second  UT1-UTC
	 //dat    second  TAI-UTC
	 //xp,yp  rad     Polar motion
	 //dx00,dy00 rad  CIP correction (Celestial Pole Offsets)

	// Function: Calculates the Rotation Matrix from J2000 (EME2000) to ITRF


	// 1. Time Systems Processing
	double rc2it[3][3];
	double tai = MJD_UTC + EOP_needed.dat / 86400.0;     // MJD (TAI)
	double tt = tai + 32.184 / 86400.0;                  // MJD (TT)
	double ut1 = MJD_UTC + EOP_needed.dut1 / 86400.0;    // MJD (UT1)

	// 2. Celestial Intermediate Pole (CIP) & Origin (CIO) - IAU2000B
	double x, y, s;
	double rc2i[3][3], rc2ti[3][3], rpom[3][3];

	// Calculate X, Y, s using IAU2000B model
	iauXys00b(2400000.5, tt, &x, &y, &s);

	// Apply CIP corrections (EOP dx, dy)
	// 【注意】: 确保 EOP_needed.dx00/dy00 单位是【弧度】。
	// 如果文件给的是 mas (毫角秒)，请在此处除以 (1000*3600*180/PI)
	x += EOP_needed.dx00;
	y += EOP_needed.dy00;

	// Form GCRS to CIRS matrix (Celestial to Intermediate)
	iauC2ixys(x, y, s, rc2i);

	// 3. Earth Rotation Angle (ERA)
	double era = iauEra00(2400000.5, ut1);

	// Form CIRS to TIRS matrix (Intermediate to Terrestrial Intermediate)
	iauCr(rc2i, rc2ti);   // Copy rc2i to rc2ti
	iauRz(era, rc2ti);    // Rotate Z-axis by ERA

	// 4. Polar Motion (TIRS -> ITRS)
	// Calculate Polar Motion Matrix
	double sp = iauSp00(2400000.5, tt); // TIO locator s'
	iauPom00(EOP_needed.xp, EOP_needed.yp, sp, rpom);

	// Form GCRS to ITRS matrix (This is your original matrix)
	// rc2it_GCRS = rpom * rc2ti
	iauRxr(rpom, rc2ti, rc2it);

	// =================================================================
	// 🛑 核心修复区：从 GCRS 转回 J2000 (Frame Bias)
	// =================================================================
	// 你的导师用的是 J2000 (EME2000)，而上面算出的是基于 GCRS 的。
	// 两者差一个微小的 Frame Bias (约 23 mas, 0.8米)。
	// 我们需要计算: Matrix_Final = Matrix_GCRS_to_ITRF * Matrix_J2000_to_GCRS

	double r_gcrs2j2k[3][3];  // Matrix: GCRS -> J2000
	double r_j2k2gcrs[3][3];  // Matrix: J2000 -> GCRS
	double r_final[3][3];     // Matrix: J2000 -> ITRF

	// 使用 iauPmat00 计算 J2000.0 时刻的 "Precession-Bias" 矩阵
	// 因为起止时间都是 J2000 (51544.5)，岁差为0，这就纯粹是 Frame Bias 矩阵
	// 注意：iauPmat00 返回的是 GCRS -> Mean J2000
	iauPmat00(2400000.5, 51544.5, r_gcrs2j2k);

	// 转置得到 J2000 -> GCRS
	iauTr(r_gcrs2j2k, r_j2k2gcrs);

	// 级联矩阵：先从 J2000 转到 GCRS，再从 GCRS 转到 ITRF
	// R_total = R(GCRS->ITRF) * R(J2000->GCRS)
	iauRxr(rc2it, r_j2k2gcrs, r_final);

	// 5. Output
	C2T_mat << r_final[0][0], r_final[0][1], r_final[0][2],
		r_final[1][0], r_final[1][1], r_final[1][2],
		r_final[2][0], r_final[2][1], r_final[2][2];
}
*/



void ITRS_to_GCRS_Matrix(double MJD_UTC, const aEOP& EOP_needed, Matrix3d& itrs_to_gcrs_emat)
{
	/*
	 utc    MJD     MJD     for a certain  moment
	 dut1   second  UT1-UTC
	 dat    second  bulletin C, TAI-UTC
	 xp,yp,dx00,dy00 rad
	 This function modify rc2ti[3][3],the matrix rorating GCRS to ITRS
	 Polar motion included
	 */
	Matrix3d gcrs_to_itrs_emat;
	GCRS_to_ITRS_Matrix(MJD_UTC, EOP_needed, gcrs_to_itrs_emat);
	//T2C_Mat = C2T_Mat.inverse();
	//直接求逆==求转置
	itrs_to_gcrs_emat = gcrs_to_itrs_emat.transpose();
}

void TOD_to_GCRS_Matrix(double MJD_UTC, double dat, Matrix3d& tod_to_gcrs_emat)
{
	// TOD 是瞬时真赤道、真春分点轴系，不是以 CIO 为起点的 CIRS。
	// iauPnm00b 返回 NPB：GCRS 到 TOD，包含框架偏置、岁差和章动。
	// 取转置后得到 TOD 到 GCRS 的旋转，不是到严格的 J2000 平赤道系。
	// 输入为 MJD(UTC) 和 TAI-UTC 秒数；内部使用 TT 计算岁差章动。


	double MJD_TT = MJD_UTC + dat / DAYSEC + 32.184 / DAYSEC;  // TT = UTC + DAT + 32s.184

	double rbpn[3][3];    //虽然写的是bpn，但是确实是NPB：V(date) = rbpn * V(GCRS)
	// 对应教材岁差章动环节；SOFA 此处还包含框架偏置 B，实际矩阵为 NPB。
	iauPnm00b(DJM0, MJD_TT, rbpn);

	tod_to_gcrs_emat << rbpn[0][0], rbpn[0][1], rbpn[0][2],
		rbpn[1][0], rbpn[1][1], rbpn[1][2],
		rbpn[2][0], rbpn[2][1], rbpn[2][2];

	// 原地转置，避免赋值时读写同一矩阵导致元素相互覆盖。
	tod_to_gcrs_emat.transposeInPlace();
}

void Horizontal_to_GCRS_Matrix(double geodetic_lambda_rad, double geodetic_phi_rad,
	double MJD_UTC, double dut1, double dat, Matrix3d& horizontal_to_gcrs_emat)
{
	// 按《航天动力学引论》式（1.50）、（1.53）构造站心地平视线的旋转。
	// 地平系为教材 O-NWZ（北、西、天顶），方位角 A 从北向东增加，
	// 因而单位视线为 [cos(h)cos(A), -cos(h)sin(A), sin(h)]。
	// 这里只变换轴向，不包含站心到地心的原点平移。
	// geodetic_lambda_rad、geodetic_phi_rad 为 rad，MJD_UTC 为 UTC 历元，dut1、dat 为秒。
	// 当前以大地经纬度代替天文经纬度，不计垂线偏差，也未单独施加极移修正。

	Matrix3d horizontal_to_tod_emat, gcrs_to_tod_emat;
	double MJD_UT1 = MJD_UTC + dut1 / DAYSEC;
	double MJD_TT = MJD_UTC + dat / DAYSEC + 32.184 / DAYSEC;  // TT = UTC + DAT + 32s.184

	//计算格林尼治真恒星时GAST(rad),输入为UT1(MJD)
	double GAST = iauGst00b(DJM0, MJD_UT1);
	double S = GAST + geodetic_lambda_rad;  //地方恒星时(rad), GAST + geodetic_lambda_rad
	// 教材式（1.53）中的 ZR：站心地平系到瞬时真赤道系。
	horizontal_to_tod_emat = Rz(PI - S) * Ry(0.5 * PI - geodetic_phi_rad);

	double rbpn[3][3];    //虽然写的是bpn，但是确实是NPB：V(date) = rbpn * V(GCRS)
	// 对应教材岁差章动环节；SOFA 此处还包含框架偏置 B，实际矩阵为 NPB。
	iauPnm00b(DJM0, MJD_TT, rbpn);
	gcrs_to_tod_emat << rbpn[0][0], rbpn[0][1], rbpn[0][2],
		rbpn[1][0], rbpn[1][1], rbpn[1][2],
		rbpn[2][0], rbpn[2][1], rbpn[2][2];

	// 先用教材 ZR 转到 TOD，再用 (NPB)^T 转到 GCRS。
	horizontal_to_gcrs_emat = gcrs_to_tod_emat.transpose() * horizontal_to_tod_emat;
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 很少用的非类函数代码，先注释掉 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
/*
void ICRF_difference(const frame::Geodetic& geodetic_class, double* ICRF_dot,
	double MJD_UTC, double dt_sec, const EOP& EOP_needed)
{
	//差分得到大地坐标系测站在ICRF坐标系关于时间的导数即速度：[X(MJD+dt)-X(MJD)]/dt 单位归一化！！
	// dt_sec 单位s

	// step1 计算t0=MJD处、t1=MJD+dt处的坐标，归一化单位
	frame::ICRF ICRF_t0, ICRF_t1;
	geodetic_class.ToICRF(MJD_UTC, EOP_needed, ICRF_t0);
	geodetic_class.ToICRF(MJD_UTC + dt_sec / 86400.0, EOP_needed, ICRF_t1);

	// step2 转为普通数组,归一化单位
	double ICRFArr_t0[3], ICRFArr_t1[3];
	ICRF_t0.ToArray(ICRFArr_t0);
	ICRF_t1.ToArray(ICRFArr_t1);

	// step3 差分得到测站速度 归一化单位
	ICRF_dot[0] = (ICRFArr_t1[0] - ICRFArr_t0[0]) / (dt_sec / 86400.0);// X_dot norm
	ICRF_dot[1] = (ICRFArr_t1[1] - ICRFArr_t0[1]) / (dt_sec / 86400.0);// Y_dot norm
	ICRF_dot[2] = (ICRFArr_t1[2] - ICRFArr_t0[2]) / (dt_sec / 86400.0);// Z_dot norm
}
*/

void GCRS_to_ITRS_Matrix(double utc, double dut1, double dat, double xp,
	double yp, double dx00, double dy00, double rc2it[3][3])
{
	/*
	 utc    MJD     MJD     for a certain  moment
	 dut1   second  UT1-UTC
	 dat    second  bulletin C, TAI-UTC
	 xp,yp,dx00,dy00 rad
	 This function modify rc2ti[3][3],the matrix rorating GCRS to ITRS
	 Polar motion included

	 */
	double tai = utc + dat / constants::SECONDS_PER_DAY;     // MJD
	double tt = tai + 32.184 / constants::SECONDS_PER_DAY;     // MJD
	double ut1 = utc + dut1 / constants::SECONDS_PER_DAY;     // MJD

	// IAU2000B method
	double x, y, s;
	double rc2i[3][3], rc2ti[3][3], rpom[3][3];

	// CIP and CIO, *******************************IAU2000B**
	iauXys00b(constants::MJD_OFFSET, tt, &x, &y, &s);

	// add CIP corrections
	x += dx00;
	y += dy00;

	// GCRS to CIRS matrix
	iauC2ixys(x, y, s, rc2i);

	// ERA
	double era = iauEra00(constants::MJD_OFFSET, ut1);

	// Form celestial-terrestrial matrix (no polar motion yet)
	iauCr(rc2i, rc2ti);
	iauRz(era, rc2ti);

	// Polar motion matrix (TIRS->ITRS,IERS2003)
	iauPom00(xp, yp, iauSp00(constants::MJD_OFFSET, tt), rpom);

	// Form celestial-terrestrial matrix(including polar motion)
	iauRxr(rpom, rc2ti, rc2it);
}
