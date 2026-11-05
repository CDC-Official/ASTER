/**
 * @file frame_definition.cpp
 * @brief 地面站址的大地坐标、ITRS 和 GCRS 位置对象及六个方向的转换。
 * @details 大地坐标采用 WGS84。对象未归一化时高度、XYZ 为 m，角度为 deg；
 * 椭球计算使用按 constants::Scale::UNIT_LENGTH_M 归一化的高度、XYZ 和 rad 制角度。
 * 地固与天球位置转换接收 UTC MJD 和同历元 EOP，保持输入的长度单位。
 */
#include "frame_definition.hpp"
#include "frame_earth.hpp"
#include "utils_constants.hpp"
#include "orbit_basic.hpp"
#include "math_basic.hpp"
#include <cmath>
#include <sofa.h>

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 类定义代码 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
frame::Geodetic::Geodetic()
{
	//默认构造函数,默认没有归一化
	m_coord[0] = 0;
	m_coord[1] = 0;
	m_coord[2] = 0;
	m_is_normalized = false;
}

frame::Geodetic::Geodetic(double height, double lambda, double phi, bool is_normalized)
{
	m_coord[0] = height;
	m_coord[1] = lambda;
	m_coord[2] = phi;
	m_is_normalized = is_normalized;

}

frame::Geodetic::~Geodetic()
{
	//do nothing
}

void frame::Geodetic::Set(double height, double lambda, double phi, bool is_normalized)
{
	m_coord[0] = height;
	m_coord[1] = lambda;
	m_coord[2] = phi;
	m_is_normalized = is_normalized;
}

double frame::Geodetic::GetValue(Overwatch& tracer, int idx) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value;
	switch (idx)
	{
		case Identifier::VAR_H_GEODETIC:
			value = m_coord[0];//高度
			break;
		case Identifier::VAR_LAMBDA_GEODETIC:
			value = m_coord[1];//经度
			break;
		case Identifier::VAR_PHI_GEODETIC:
			value = m_coord[2];//纬度
			break;
		default:
		{
			value = 99999;
			std::string msg;
			msg = "Wrong index when get value from a geodetic class, return 99999 in default!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			break;
		}

	}

	tracer.CheckOut(__func__);
	return value;
}

void frame::Geodetic::Normalize(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);
	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length smaller than 0!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;// 出错直接return
	}

	if (m_is_normalized == true)
	{
		std::string msg = "Geodetic coordinates already normalized, this normalization is no longer performed!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else if (m_is_normalized == false)
	{
		m_coord[0] /= uLength;
		m_coord[1] *= DD2R;
		m_coord[2] *= DD2R;
		m_is_normalized = true;
	}

	tracer.CheckOut(__func__);
}

void frame::Geodetic::Restore(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);
	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length smaller than 0!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;// 出错直接return
	}

	if (m_is_normalized == false)
	{
		std::string msg = "Geodetic coordinates already restored, this restoration is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else if (m_is_normalized == true)
	{
		m_coord[0] *= uLength;
		m_coord[1] /= DD2R;
		m_coord[2] /= DD2R;
		m_is_normalized = false;
	}
	tracer.CheckOut(__func__);
}

/**
 * @brief 将 WGS84 大地经纬高转换为 ITRS 地固位置。
 * @param[in,out] tracer 记录本次转换的调用链和输入错误。
 * @param[out] itrs_position 地固位置 [X,Y,Z]，以 constants::Scale::UNIT_LENGTH_M 归一化。
 * @details 当前对象须已归一化：高度除以 UNIT_LENGTH_M，经度和大地纬度为 rad。
 * 使用 constants::Wgs84 的赤道半径与扁率；失败时保留传入的输出对象。
 */
void frame::Geodetic::ToITRS(Overwatch& tracer, frame::ITRS& itrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 椭球计算使用归一化高度和弧度制经纬度。
	if (!m_is_normalized || !std::isfinite(m_coord[0]) || !std::isfinite(m_coord[1])
		|| !std::isfinite(m_coord[2]) || std::fabs(m_coord[2]) > constants::PI / 2.0)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected normalized, finite geodetic coordinates with latitude in [-pi/2, pi/2].", EventType::ERROR);
		return;
	}

	// WGS84 的曲率半径先除以数值长度单位，与输入高度配套。
	const double height = m_coord[0], lambda = m_coord[1], phi = m_coord[2];
	const double N = constants::Wgs84::EQUATORIAL_RADIUS / constants::Scale::UNIT_LENGTH_M
		* pow(1 - constants::Wgs84::ECCENTRICITY_SQUARED * pow(sin(phi), 2), -0.5);

	// 由椭球法线方向及高度计算地固 XYZ，输出仍为归一化位置。
	const double X = (N + height) * cos(phi) * cos(lambda);
	const double Y = (N + height) * cos(phi) * sin(lambda);
	const double Z = (N * (1 - constants::Wgs84::ECCENTRICITY_SQUARED) + height) * sin(phi);
	itrs_position.Set(X, Y, Z, true);

	tracer.CheckOut(__func__);
}

void frame::Geodetic::ToArray(double* geodetic_array) const
{
	geodetic_array[0] = m_coord[0];
	geodetic_array[1] = m_coord[1];
	geodetic_array[2] = m_coord[2];
}

/**
 * @brief 将 WGS84 大地经纬高转换为指定时刻的 GCRS 地心天球位置。
 * @param[in,out] tracer 记录椭球转换和地球定向旋转的调用链及错误。
 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
 * @param[out] gcrs_position GCRS 位置 [X,Y,Z]，以 constants::Scale::UNIT_LENGTH_M 归一化。
 * @details 当前对象须已归一化：高度除以 UNIT_LENGTH_M，角度为 rad。
 * 依次调用 ToITRS 和 ITRS::ToGCRS；采用与初轨程序相同的 IAU 2000B 地球定向旋转。
 * 失败时保留传入的输出对象。
 */
void frame::Geodetic::ToGCRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::GCRS& gcrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 先把 WGS84 大地站址转换为地固 XYZ。
	frame::ITRS itrs_position;
	this->ToITRS(tracer, itrs_position);
	if (tracer.Code() != ERR_NOERROR) return;

	// 再按观测时刻和 EOP 旋转到 GCRS，供观测方程与外部初轨算法使用。
	itrs_position.ToGCRS(tracer, MJD_UTC, EOP_needed, gcrs_position);
	if (tracer.Code() != ERR_NOERROR) return;

	tracer.CheckOut(__func__);
}

frame::ITRS::ITRS()
{
	m_coord[0] = 0;
	m_coord[1] = 0;
	m_coord[2] = 0;
	m_is_normalized = false;
}

frame::ITRS::ITRS(double X, double Y, double Z, bool is_normalized)
{

	m_coord[0] = X;
	m_coord[1] = Y;
	m_coord[2] = Z;
	m_is_normalized = is_normalized;
}

frame::ITRS::~ITRS()
{
	//do nothing
}

void frame::ITRS::Set(double X, double Y, double Z, bool is_normalized)
{
	m_coord[0] = X;
	m_coord[1] = Y;
	m_coord[2] = Z;
	m_is_normalized = is_normalized;
}

double frame::ITRS::GetValue(Overwatch& tracer, int idx) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value;
	switch (idx)
	{
		case Identifier::VAR_CARTX:
			value = m_coord[0];//X
			break;
		case Identifier::VAR_CARTY:
			value = m_coord[1];//Y
			break;
		case Identifier::VAR_CARTZ:
			value = m_coord[2];//Z
			break;
		default:
		{
			value = 99999;
			std::string msg = "Wrong index, return 99999 in default!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			break;
		}

	}
	tracer.CheckOut(__func__);
	return value;
}

void frame::ITRS::Normalize(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length smaller than 0!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;// 出错直接return
	}

	if (m_is_normalized == true)
	{
		std::string msg = "ITRS coordinates already normalized, this normalization is no longer performed!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else if (m_is_normalized == false)
	{
		m_coord[0] /= uLength;
		m_coord[1] /= uLength;
		m_coord[2] /= uLength;
		m_is_normalized = true;
	}

	tracer.CheckOut(__func__);
}

void frame::ITRS::Restore(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length smaller than 0!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;// 出错直接return
	}
	if (m_is_normalized == false)
	{
		std::string msg = "ITRS coordinates already restored, this restoration is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else if (m_is_normalized == true)
	{
		m_coord[0] *= uLength;
		m_coord[1] *= uLength;
		m_coord[2] *= uLength;
		m_is_normalized = false;
	}

	tracer.CheckOut(__func__);
}

void frame::ITRS::ToArray(double* itrs_array) const
{
	itrs_array[0] = m_coord[0];
	itrs_array[1] = m_coord[1];
	itrs_array[2] = m_coord[2];
}

/**
 * @brief 将 ITRS 地固位置转换为 WGS84 大地经纬高。
 * @param[in,out] tracer 记录本次反算的调用链和输入错误。
 * @param[out] geodetic_position 大地坐标 [高度,经度,大地纬度]；高度以
 * constants::Scale::UNIT_LENGTH_M 归一化，角度为 rad。
 * @details 当前对象须为以 UNIT_LENGTH_M 归一化的有限、非零地固位置。
 * 使用 SOFA iauGc2gde 和 WGS84 椭球参数；自转轴上经度取 0，失败时保留输出。
 */
void frame::ITRS::ToGeodetic(Overwatch& tracer, frame::Geodetic& geodetic_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 地心零向量没有测站经纬高；其余输入须为有限的归一化地固坐标。
	if (!m_is_normalized || !std::isfinite(m_coord[0]) || !std::isfinite(m_coord[1])
		|| !std::isfinite(m_coord[2]) || (m_coord[0] == 0.0 && m_coord[1] == 0.0 && m_coord[2] == 0.0))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected a normalized, finite, nonzero ITRS position.", EventType::ERROR);
		return;
	}

	// SOFA 的椭球半径与 XYZ 使用相同尺度，反算高度也保持该尺度。
	double xyz[3] = {m_coord[0], m_coord[1], m_coord[2]};
	double height, lambda, phi;
	const int status = iauGc2gde(constants::Wgs84::EQUATORIAL_RADIUS / constants::Scale::UNIT_LENGTH_M,
		constants::Wgs84::FLATTENING, xyz, &lambda, &phi, &height);

	if (status != 0 || !std::isfinite(height) || !std::isfinite(lambda) || !std::isfinite(phi))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "SOFA WGS84 geodetic conversion failed.", EventType::ERROR);
		return;
	}

	// 类对象以高度、经度、纬度保存，恢复米和度由调用方执行。
	geodetic_position.Set(height, lambda, phi, true);
	tracer.CheckOut(__func__);
}

/**
 * @brief 将 ITRS 地固位置旋转为指定时刻的 GCRS 地心天球位置。
 * @param[in,out] tracer 记录本次旋转的调用链和输入错误。
 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
 * @param[out] gcrs_position GCRS 位置 [X,Y,Z]，长度单位及归一化状态与当前对象相同。
 * @details 复用 ITRS_to_GCRS_Matrix 的 IAU 2000B 位置旋转；失败时保留输出对象。
 */
void frame::ITRS::ToGCRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::GCRS& gcrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 旋转需要有限的位置、UTC 时刻及 s、rad 制的地球定向参数。
	if (!std::isfinite(m_coord[0]) || !std::isfinite(m_coord[1]) || !std::isfinite(m_coord[2])
		|| !std::isfinite(MJD_UTC) || !std::isfinite(EOP_needed.dut1) || !std::isfinite(EOP_needed.dat)
		|| !std::isfinite(EOP_needed.xp) || !std::isfinite(EOP_needed.yp)
		|| !std::isfinite(EOP_needed.dx00) || !std::isfinite(EOP_needed.dy00))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected finite position, UTC MJD and EOP values.", EventType::ERROR);
		return;
	}

	// 复用现有地球定向矩阵，位置旋转不改变长度单位。
	Matrix3d itrs_to_gcrs_emat;
	ITRS_to_GCRS_Matrix(MJD_UTC, EOP_needed, itrs_to_gcrs_emat);
	const Vector3d input(m_coord[0], m_coord[1], m_coord[2]);
	const Vector3d result = itrs_to_gcrs_emat * input;

	// 写入目标坐标，归一化标志随输入位置保留。
	gcrs_position.Set(result[0], result[1], result[2], m_is_normalized);
	tracer.CheckOut(__func__);
}

frame::GCRS::GCRS()
{
	m_coord[0] = 0;
	m_coord[1] = 0;
	m_coord[2] = 0;
	m_is_normalized = false;
}

frame::GCRS::GCRS(double X, double Y, double Z, bool is_normalized)
{
	m_coord[0] = X;
	m_coord[1] = Y;
	m_coord[2] = Z;
	m_is_normalized = is_normalized;
}

frame::GCRS::~GCRS()
{
	// do nothing
}

void frame::GCRS::Set(double X, double Y, double Z, bool is_normalized)
{
	m_coord[0] = X;
	m_coord[1] = Y;
	m_coord[2] = Z;
	m_is_normalized = is_normalized;
}

double frame::GCRS::GetValue(Overwatch& tracer, int idx) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value;
	switch (idx)
	{
		case Identifier::VAR_CARTX:
			value = m_coord[0];//X
			break;
		case Identifier::VAR_CARTY:
			value = m_coord[1];//Y
			break;
		case Identifier::VAR_CARTZ:
			value = m_coord[2];//Z
			break;
		default:
		{
			value = 99999;
			std::string msg = "Wrong index, return 99999 in default!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			break;
		}

	}
	tracer.CheckOut(__func__);
	return value;
}

void frame::GCRS::Normalize(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length smaller than 0!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;// 出错直接return
	}

	if (m_is_normalized == true)
	{
		std::string msg = "GCRS coordinates already normalized, this normalization is no longer performed!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else if (m_is_normalized == false)
	{
		m_coord[0] /= uLength;
		m_coord[1] /= uLength;
		m_coord[2] /= uLength;
		m_is_normalized = true;
	}

	tracer.CheckOut(__func__);
}

void frame::GCRS::Restore(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length smaller than 0!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;// 出错直接return
	}
	if (m_is_normalized == false)
	{
		std::string msg = "GCRS coordinates already restored, this restoration is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else if (m_is_normalized == true)
	{
		m_coord[0] *= uLength;
		m_coord[1] *= uLength;
		m_coord[2] *= uLength;
		m_is_normalized = false;
	}

	tracer.CheckOut(__func__);
}

void frame::GCRS::ToArray(double* gcrs_array) const
{
	gcrs_array[0] = m_coord[0];
	gcrs_array[1] = m_coord[1];
	gcrs_array[2] = m_coord[2];
}

/**
 * @brief 将指定时刻的 GCRS 地心天球位置转换为 WGS84 大地经纬高。
 * @param[in,out] tracer 记录地球定向旋转和椭球反算的调用链及错误。
 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
 * @param[out] geodetic_position 大地坐标 [高度,经度,大地纬度]；高度以
 * constants::Scale::UNIT_LENGTH_M 归一化，角度为 rad。
 * @details 当前对象须已按 UNIT_LENGTH_M 归一化。依次调用 ToITRS 和 ITRS::ToGeodetic，
 * 采用 IAU 2000B 位置旋转及 WGS84 椭球；失败时保留输出对象。
 */
void frame::GCRS::ToGeodetic(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::Geodetic& geodetic_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 椭球反算使用按 UNIT_LENGTH_M 归一化的位置。
	if (!m_is_normalized)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Normalize GCRS before transforming to geodetic.", EventType::ERROR);
		return;
	}

	// 先按同一时刻的地球姿态转回 ITRS，再在 WGS84 椭球上求经纬高。
	frame::ITRS itrs_position;
	this->ToITRS(tracer, MJD_UTC, EOP_needed, itrs_position);
	if (tracer.Code() != ERR_NOERROR) return;

	itrs_position.ToGeodetic(tracer, geodetic_position);
	if (tracer.Code() != ERR_NOERROR) return;

	tracer.CheckOut(__func__);
}

/**
 * @brief 将 GCRS 地心天球位置旋转为指定时刻的 ITRS 地固位置。
 * @param[in,out] tracer 记录本次旋转的调用链和输入错误。
 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
 * @param[out] itrs_position 地固位置 [X,Y,Z]，长度单位及归一化状态与当前对象相同。
 * @details 复用 GCRS_to_ITRS_Matrix 的 IAU 2000B 位置旋转；失败时保留输出对象。
 */
void frame::GCRS::ToITRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::ITRS& itrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 旋转需要有限的位置、UTC 时刻及 s、rad 制的地球定向参数。
	if (!std::isfinite(m_coord[0]) || !std::isfinite(m_coord[1]) || !std::isfinite(m_coord[2])
		|| !std::isfinite(MJD_UTC) || !std::isfinite(EOP_needed.dut1) || !std::isfinite(EOP_needed.dat)
		|| !std::isfinite(EOP_needed.xp) || !std::isfinite(EOP_needed.yp)
		|| !std::isfinite(EOP_needed.dx00) || !std::isfinite(EOP_needed.dy00))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected finite position, UTC MJD and EOP values.", EventType::ERROR);
		return;
	}

	// 复用现有地球定向矩阵，位置旋转不改变长度单位。
	Matrix3d gcrs_to_itrs_emat;
	GCRS_to_ITRS_Matrix(MJD_UTC, EOP_needed, gcrs_to_itrs_emat);
	const Vector3d input(m_coord[0], m_coord[1], m_coord[2]);
	const Vector3d result = gcrs_to_itrs_emat * input;

	// 写入目标坐标，归一化标志随输入位置保留。
	itrs_position.Set(result[0], result[1], result[2], m_is_normalized);
	tracer.CheckOut(__func__);
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 仍然有大用的中层代码（广泛用于类中） @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
