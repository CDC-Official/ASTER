#include "frame_definition.hpp"
#include "frame_earth.hpp"
#include "utils_constants.hpp"
#include "orbit_basic.hpp"
#include "math_basic.hpp"
#include <cmath>

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

void frame::Geodetic::ToITRS(Overwatch& tracer, frame::ITRS& itrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (m_is_normalized == false)
	{
		std::string msg = "Normalize the geodetic coordindates before transforming to ITRS!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
        return;
	}
	else
	{
		// 以 WGS84 参考椭球将测站的大地经纬高转换为 ITRS XYZ。
		// height 和输出 XYZ 均以 ae 归一化，lambda、phi 为 rad；phi 是大地纬度。
		double height = m_coord[0], lambda = m_coord[1], phi = m_coord[2];
		double N = constants::Wgs84::EQUATORIAL_RADIUS / constants::Scale::UNIT_LENGTH_M * pow(1 - constants::Wgs84::ECCENTRICITY_SQUARED * pow(sin(phi), 2), -0.5);    // N 为以 ae 归一化的卯酉圈曲率半径。
		double X = (N + height) * cos(phi) * cos(lambda);
		double Y = (N + height) * cos(phi) * sin(lambda);
		double Z = (N * (1 - constants::Wgs84::ECCENTRICITY_SQUARED) + height) * sin(phi);

		itrs_position.Set(X, Y, Z, true);
	}
	tracer.CheckOut(__func__);
}

void frame::Geodetic::ToArray(double* geodetic_array) const
{
	geodetic_array[0] = m_coord[0];
	geodetic_array[1] = m_coord[1];
	geodetic_array[2] = m_coord[2];
}

void frame::Geodetic::ToGCRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::GCRS& gcrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	//需要归一化，否则在ToITRS方法中会报错
	frame::ITRS itrs_position;
	this->ToITRS(tracer, itrs_position);
	itrs_position.ToGCRS(tracer, MJD_UTC, EOP_needed, gcrs_position);

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

void frame::ITRS::ToGeodetic(Overwatch& tracer, frame::Geodetic& geodetic_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 必须要归一化
	if (m_is_normalized == false)
	{
		std::string msg = "Normalize the ITRS coordindates before transforming to geodetic!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		// ITRS XYZ 与高度都以 ae 归一化，经纬度结果为 rad。
		// 使用与正转换相同的现有椭球常量；这是坐标表示转换，不是基准转换。
		double X = m_coord[0], Y = m_coord[1], Z = m_coord[2];
		// 计算经度 (lambda)
		double lambda = atan2(Y, X);
		// p 是位置到地球自转轴的距离，以 ae 归一化。
		double p = sqrt(X * X + Y * Y);
		// 用椭球偏心率构造大地纬度的初始估计，再迭代求解。
		double theta = atan2(Z, p * (1 - constants::Wgs84::ECCENTRICITY_SQUARED));
		// phi 为大地纬度；每次按当前 phi 更新曲率半径 N。
		double sin_phi, N;
		double phi = theta;

		while (true)
		{
			sin_phi = sin(phi);
			N = constants::Wgs84::EQUATORIAL_RADIUS / constants::Scale::UNIT_LENGTH_M / sqrt(1 - constants::Wgs84::ECCENTRICITY_SQUARED * sin_phi * sin_phi);  // 卯酉圈曲率半径，以 ae 归一化。

			double phi_next = atan2(Z + constants::Wgs84::ECCENTRICITY_SQUARED * N * sin_phi, p);
			if (fabs(phi_next - phi) < 1e-14)
			{  // 收敛条件
				break;
			}
			phi = phi_next;
		}
		// 计算大地高 (height)
		double height = p / cos(phi) - N;
		// atan2 和迭代结果已经是 rad，直接按归一化状态保存。
		// 需要米制高度和度制经纬度时，由 Restore(tracer, ae) 统一恢复。
		geodetic_position.Set(height, lambda, phi, true);
	}

	tracer.CheckOut(__func__);
}

void frame::ITRS::ToGCRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::GCRS& gcrs_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	//归不归一化其实无所谓，不过最好在归一化下进行
	Matrix3d itrs_to_gcrs_emat;
	ITRS_to_GCRS_Matrix(MJD_UTC, EOP_needed, itrs_to_gcrs_emat);
	Vector3d itrs_evec, gcrs_evec;
	itrs_evec << m_coord[0], m_coord[1], m_coord[2];
	gcrs_evec = itrs_to_gcrs_emat * itrs_evec;
	gcrs_position.Set(gcrs_evec[0], gcrs_evec[1], gcrs_evec[2], m_is_normalized);

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

void frame::GCRS::ToGeodetic(Overwatch& tracer,
	double MJD_UTC,
	const aEOP& EOP_needed,
	frame::Geodetic& geodetic_position) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (m_is_normalized == false)
	{
		std::string msg = "Normalize GCRS before transforming to geodetic!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		frame::ITRS itrs_position;
		this->ToITRS(tracer, MJD_UTC, EOP_needed, itrs_position);
		itrs_position.ToGeodetic(tracer, geodetic_position);
	}

	tracer.CheckOut(__func__);
}

void frame::GCRS::ToITRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, frame::ITRS& itrs_position) const
{
	//归不归一化其实无所谓，不过最好在归一化下进行
	Matrix3d gcrs_to_itrs_emat;
	GCRS_to_ITRS_Matrix(MJD_UTC, EOP_needed, gcrs_to_itrs_emat);
	Vector3d itrs_evec, gcrs_evec;
	gcrs_evec << m_coord[0], m_coord[1], m_coord[2];
	itrs_evec = gcrs_to_itrs_emat * gcrs_evec;
	itrs_position.Set(itrs_evec[0], itrs_evec[1], itrs_evec[2], m_is_normalized);
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 仍然有大用的中层代码（广泛用于类中） @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
