//
//  orbit_statedef.cpp
//  CDC_lib
//
//  Created by 陈酊存 on 2023/9/5.
//

#include "orbit_statedef.hpp"
using namespace constants::Conversion;
using constants::PI;

typedef bool b;
using namespace Identifier;

orbit::Cartesian::Cartesian()
{
	//默认构造函数,默认没有归一化
	//m_t_MJD_UTC = 0;
	m_x = 0;
	m_y = 0;
	m_z = 0;
	m_vx = 0;
	m_vy = 0;
	m_vz = 0;
	m_r = 0;
	m_v = 0;
	m_is_normalized = false;
}

orbit::Cartesian::Cartesian(
	double x,
	double y,
	double z,
	double vx,
	double vy,
	double vz,
	bool is_normalized)
{
	//m_t_MJD_UTC = t_MJD_UTC;从参数列表中删除了t_MJD_UTC
	m_x = x;
	m_y = y;
	m_z = z;
	m_vx = vx;
	m_vy = vy;
	m_vz = vz;
	m_r = sqrt(x * x + y * y + z * z);
	m_v = sqrt(vx * vx + vy * vy + vz * vz);
	m_is_normalized = is_normalized;
}

orbit::Cartesian::~Cartesian()
{
	//std::cout << "Cartesian destructed." << std::endl;
}

void orbit::Cartesian::Show(std::ostream& os, int precision) const
{
	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_x;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_y;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_z;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_vx;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_vy;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_vz;
	//std::cout << std::endl;
}

void orbit::Cartesian::Set(double x, double y, double z, double vx, double vy, double vz, bool is_normalized)
{
	//m_t_MJD_UTC = t_MJD_UTC;从参数列表中删除了t_MJD_UTC
	m_x = x;
	m_y = y;
	m_z = z;
	m_vx = vx;
	m_vy = vy;
	m_vz = vz;
	m_r = sqrt(x * x + y * y + z * z);
	m_v = sqrt(vx * vx + vy * vy + vz * vz);
	m_is_normalized = is_normalized;
}

double orbit::Cartesian::GetValue(Overwatch& tracer, int varId) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value;
	switch (varId)
	{
		case VAR_CARTX:
			value = m_x;
			break;
		case VAR_CARTY:
			value = m_y;
			break;
		case VAR_CARTZ:
			value = m_z;
			break;
		case VAR_CARTVX:
			value = m_vx;
			break;
		case VAR_CARTVY:
			value = m_vy;
			break;
		case VAR_CARTVZ:
			value = m_vz;
			break;
		case VAR_CARTPOS:
		case VAR_R:
			value = m_r;
			break;
		case VAR_CARTVEL:
			value = m_v;
			break;

		default:
			std::string msg = "Wrong index, return 99999 in default!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			value = 99999;
			break;
	}
	tracer.CheckOut(__func__);
	return value;
}

void orbit::Cartesian::Normalize(Overwatch& tracer, double uLength, double uTime)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0 || uTime <= 0)
	{
		std::string msg = "Invalid unit length and/or unit time!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	if (m_is_normalized)
	{
		std::string msg = "Cartesian already normalized, this normalization is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		double uVelocity = uLength / uTime;

		m_x /= uLength;
		m_y /= uLength;
		m_z /= uLength;
		m_vx /= uVelocity;
		m_vy /= uVelocity;
		m_vz /= uVelocity;
		m_r /= uLength;
		m_v /= uVelocity;
		m_is_normalized = true;
	}

	tracer.CheckOut(__func__);
}

void orbit::Cartesian::Normalized(Overwatch& tracer, double uLength, double uTime, Cartesian& new_cart) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0 || uTime <= 0)
	{
		std::string msg = "Invalid unit length and/or unit time!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	if (m_is_normalized)
	{
		std::string msg = "Cartesian already normalized, give value to return object directly.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);

		new_cart.Set(this->m_x,
			this->m_y,
			this->m_z,
			this->m_vx,
			this->m_vy,
			this->m_vz, true);
	}
	else
	{
		double uVelocity = uLength / uTime;
		new_cart.Set(this->m_x / uLength,
			this->m_y / uLength,
			this->m_z / uLength,
			this->m_vx / uVelocity,
			this->m_vy / uVelocity,
			this->m_vz / uVelocity, true);
	}
	tracer.CheckOut(__func__);
}

void orbit::Cartesian::Restore(Overwatch& tracer, double uLength, double uTime)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Cartesian already restored, this restoration is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		double uVelocity = uLength / uTime;

		m_x *= uLength;
		m_y *= uLength;
		m_z *= uLength;

		m_vx *= uVelocity;
		m_vy *= uVelocity;
		m_vz *= uVelocity;

		m_r *= uLength;
		m_v *= uVelocity;
		m_is_normalized = false;
	}

	tracer.CheckOut(__func__);
}

void orbit::Cartesian::ToArray(double* cart_array) const
{
	cart_array[0] = m_x;
	cart_array[1] = m_y;
	cart_array[2] = m_z;
	cart_array[3] = m_vx;
	cart_array[4] = m_vy;
	cart_array[5] = m_vz;
}

void orbit::Cartesian::ToKeplerian(Overwatch& tracer, double* kep_array) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Cartesian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		CartesianToKeplerian(m_x, m_y, m_z, m_vx, m_vy, m_vz, kep_array);
	}
	tracer.CheckOut(__func__);
}

void orbit::Cartesian::ToKeplerian(Overwatch& tracer, orbit::Keplerian& kep_class) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Cartesian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		double temp_kep[6];
		CartesianToKeplerian(m_x, m_y, m_z, m_vx, m_vy, m_vz, temp_kep);
		kep_class.Set(temp_kep[0], temp_kep[1], temp_kep[2], temp_kep[3], temp_kep[4], temp_kep[5], true);
	}
	tracer.CheckOut(__func__);
}

void orbit::Cartesian::ToNonSingular2(Overwatch& tracer, double* nsg2_array) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Cartesian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	Keplerian kep_class;
	this->ToKeplerian(tracer, kep_class);
	if (tracer.Code() != ERR_NOERROR) return;// 这里其实不用做异常追踪，因为进入到这里的逻辑会保证调用函数时不会报错，但是保险起见还是加一下

	kep_class.ToNonSingular2(tracer, nsg2_array);
	if (tracer.Code() != ERR_NOERROR) return;

	tracer.CheckOut(__func__);
}

void orbit::Cartesian::ToNonSingular2(Overwatch& tracer, orbit::NonSingular2& nsg2_class) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Cartesian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	Keplerian kep_class;
	this->ToKeplerian(tracer, kep_class);
	if (tracer.Code() != ERR_NOERROR) return;

	kep_class.ToNonSingular2(tracer, nsg2_class);
	if (tracer.Code() != ERR_NOERROR) return;

	tracer.CheckOut(__func__);
}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


orbit::Keplerian::Keplerian()
{
	//默认构造函数,默认没有归一化
	//m_t_MJD_UTC = 0;
	m_a = 0;
	m_e = 0;
	m_i = 0;
	m_O = 0;
	m_o = 0;
	m_M = 0;
	m_is_normalized = false;
}

orbit::Keplerian::Keplerian(double a, double e, double i,
	double OMEGA, double omega, double M, bool is_normalized)//自定义构造函数
{
	//m_t_MJD_UTC = 0;
	m_a = a;
	m_e = e;
	m_i = i;
	m_O = OMEGA;
	m_o = omega;
	m_M = M;

	if (is_normalized)
	{
		//转换为0-2pi
		m_O = fmod(OMEGA, 2 * PI);
		if (m_O < 0)
			m_O += 2 * PI;

		m_o = fmod(omega, 2 * PI);
		if (m_o < 0)
			m_o += 2 * PI;

		m_M = fmod(M, 2 * PI);
		if (m_M < 0)
			m_M += 2 * PI;
	}

	m_is_normalized = is_normalized;
}

orbit::Keplerian::~Keplerian()
{
	//std::cout << "Keplerian destructed." << std::endl;
}

void orbit::Keplerian::Set(double a, double e, double i, double OMEGA, double omega, double M, bool is_normalized)
{
	//m_t_MJD_UTC = 0;
	m_is_fullySet = false; // Invalidate cached anomalies and radius.
	m_a = a;
	m_e = e;
	m_i = i;
	m_O = OMEGA;
	m_o = omega;
	m_M = M;
	m_is_normalized = is_normalized;
	if (m_is_normalized)
	{
		//转换为0-2pi
		m_O = fmod(OMEGA, 2 * PI);
		if (m_O < 0)
			m_O += 2 * PI;

		m_o = fmod(omega, 2 * PI);
		if (m_o < 0)
			m_o += 2 * PI;

		m_M = fmod(M, 2 * PI);
		if (m_M < 0)
			m_M += 2 * PI;
	}

}

void orbit::Keplerian::Show(std::ostream& os, int precision) const
{
	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_a;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_e;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_i;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_O;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_o;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_M;

	//os << std::endl;
}

void orbit::Keplerian::Normalize(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	if (m_is_normalized)
	{
		std::string msg = "Keplerian already norlamlized, this normalization is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		m_a /= uLength;
		m_i *= DEGREE_TO_RAD;
		m_O *= DEGREE_TO_RAD;
		m_o *= DEGREE_TO_RAD;
		m_M *= DEGREE_TO_RAD;

		if (m_is_fullySet)
		{
			//首次归一化时不会执行，因为计算这些值必须要先归一化；不排除后续又归一化
			m_r /= uLength;
			m_E *= DEGREE_TO_RAD;
			m_f *= DEGREE_TO_RAD;
			m_u *= DEGREE_TO_RAD;
		}
		m_is_normalized = true;
	}
	tracer.CheckOut(__func__);
}

void orbit::Keplerian::Normalized(Overwatch& tracer, double uLength, Keplerian& new_kep) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	if (m_is_normalized)
	{
		new_kep.Set(this->m_a, this->m_e, this->m_i, this->m_O, this->m_o, this->m_M, true);
	}
	else
	{
		new_kep.Set(this->m_a / uLength, this->m_e, this->m_i * DEGREE_TO_RAD,
			this->m_O * DEGREE_TO_RAD, this->m_o * DEGREE_TO_RAD, this->m_M * DEGREE_TO_RAD, true);
	}

	tracer.CheckOut(__func__);

}

void orbit::Keplerian::Restore(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (uLength <= 0)
	{
		std::string msg = "Invalid unit length!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	if (!m_is_normalized)
	{
		std::string msg = "Keplerian already restored, this restoration is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		m_a *= uLength;
		m_i *= RAD_TO_DEGREE;
		m_O *= RAD_TO_DEGREE;
		m_o *= RAD_TO_DEGREE;
		m_M *= RAD_TO_DEGREE;

		if (m_is_fullySet)
		{
			m_r *= uLength;
			m_E *= RAD_TO_DEGREE;
			m_f *= RAD_TO_DEGREE;
			m_u *= RAD_TO_DEGREE;
		}
		m_is_normalized = false;
	}
	tracer.CheckOut(__func__);
}

double orbit::Keplerian::GetValue(Overwatch& tracer, int varId) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value;
	switch (varId)
	{
		case VAR_SEMIMAJORAXIS:
			value = m_a;
			break;
		case VAR_ECCENTRICITY:
			value = m_e;
			break;
		case VAR_INCLINATION:
			value = m_i;
			break;
		case VAR_RAAN:
			value = m_O;
			break;
		case VAR_ARGPERIGEE:
			value = m_o;
			break;
		case VAR_MEANANOMALY:
			value = m_M;
			break;
		case VAR_ECCANOMALY:
			if (m_is_fullySet)
			{
				value = m_E;
			}
			else
			{
				value = 99999;
				std::string msg = "Eccentric anomaly E is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return value;
			}
			break;
		case VAR_TRUEANOMALY:
			if (m_is_fullySet)
			{
				value = m_f;
			}
			else
			{
				value = 99999;
				std::string msg = "True anomaly f is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return value;
			}
			break;
		case VAR_R:
			if (m_is_fullySet)
			{
				value = m_r;
			}
			else
			{
				value = 99999;
				std::string msg = "Geocentric distance r is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return value;
			}
			break;
		case VAR_U:
			if (m_is_fullySet)
			{
				value = m_u;
			}
			else
			{
				value = 99999;
				std::string msg = "u = w + f is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return value;
			}
			break;
		default:
			value = 99999;
			std::string msg = "Wrong index, return 99999 in default!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return value;
	}
	tracer.CheckOut(__func__);
	return value;
}

void orbit::Keplerian::ToCartesian(Overwatch& tracer, double* cart_array) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Keplerian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		KeplerianToCartsian(m_a, m_e, m_i, m_O, m_o, m_M, cart_array);
	}

	tracer.CheckOut(__func__);
}

void orbit::Keplerian::ToCartesian(Overwatch& tracer, orbit::Cartesian& cart_class) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Keplerian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		double temp_cart[6];
		KeplerianToCartsian(m_a, m_e, m_i, m_O, m_o, m_M, temp_cart);
		cart_class.Set(temp_cart[0], temp_cart[1], temp_cart[2], temp_cart[3], temp_cart[4], temp_cart[5], true);
	}
	tracer.CheckOut(__func__);
}

void orbit::Keplerian::ToNonSingular2(Overwatch& tracer, double* nsg2_array, bool strict_check) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);
	if (!this->ToNonSingular2Array(tracer, nsg2_array, strict_check)) return;
	tracer.CheckOut(__func__);
}

bool orbit::Keplerian::ToNonSingular2Array(Overwatch& tracer, double* nsg2_array, bool strict_check) const
{

	if (!m_is_normalized)
	{
		std::string msg = "Keplerian is not normalized, you must normalize it first before transform!";
		tracer.LogEvent(__FILE__, "ToNonSingular2", __LINE__, msg, EventType::ERROR);
		return false;
	}

	if (strict_check && (m_a <= 0.0))
	{
		std::string msg = "Negative semi-major axis for ellipse!";
		tracer.LogEvent(__FILE__, "ToNonSingular2", __LINE__, msg, EventType::ERROR);
		return false;
	}
	if (strict_check && (m_e < 0.0))
	{
		std::string msg = "Invalid eccentricity < 0 !";
		tracer.LogEvent(__FILE__, "ToNonSingular2", __LINE__, msg, EventType::ERROR);
		return false;
	}
	if (strict_check && (m_e >= 1.0))
	{
		std::string
			msg = "Keplerian to NonSingular2 convertion finds a hyperbola of e = " + std::to_string(m_e) + ".";
		tracer.LogEvent(__FILE__, "ToNonSingular2", __LINE__, msg, EventType::WARNING);
	}

	double sini2 = sin(m_i / 2.0);
	double ow = m_O + m_o;
	nsg2_array[0] = m_a;
	nsg2_array[1] = m_e * cos(ow);
	nsg2_array[2] = m_e * sin(ow);
	// 260204 修复m_O和m_o写反的bug
	nsg2_array[3] = sini2 * cos(m_O);
	nsg2_array[4] = sini2 * sin(m_O);
	nsg2_array[5] = fmod(ow + m_M + 2 * PI, 2 * PI);

	return true;
}

void orbit::Keplerian::ToNonSingular2(Overwatch& tracer, orbit::NonSingular2& nsg2_class, bool strict_check) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double nsg2_array[6];
	if (!this->ToNonSingular2Array(tracer, nsg2_array, strict_check)) return;

	nsg2_class.Set(tracer, nsg2_array[0], nsg2_array[1], nsg2_array[2],
		nsg2_array[3], nsg2_array[4], nsg2_array[5], true);
	if (tracer.Code() != ERR_NOERROR) return;

	tracer.CheckOut(__func__);
}

void orbit::Keplerian::ToArray(double* kep_array) const
{
	kep_array[0] = m_a;
	kep_array[1] = m_e;
	kep_array[2] = m_i;
	kep_array[3] = m_O;
	kep_array[4] = m_o;
	kep_array[5] = m_M;
}

void orbit::Keplerian::SyncAnomaly(Overwatch& tracer)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "Keplerian is not normalized, you must normalize it first before synchronize anomaly!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{
		m_E = MeanAnomaly2EccentricAnomaly(m_M, m_e);
		m_f = EccentricAnomaly2TrueAnomaly(m_E, m_e);
		m_u = fmod(m_o + m_f, 2 * PI);
		if (m_u < 0)
		{
			m_u += 2 * PI;
		}

		m_r = m_a * (1.0 - m_e * cos(m_E));
		m_is_fullySet = true;
	}
	tracer.CheckOut(__func__);
}

bool orbit::Keplerian::IsNormalized() const
{
	return m_is_normalized;
}

orbit::NonSingular2::NonSingular2()
{
	//默认构造函数,默认没有归一化
	m_a = 0;
	m_xi = 0;
	m_eta = 0;
	m_h = 0;
	m_k = 0;
	m_lambda = 0;
	m_is_normalized = false;
	// 以下可有可无
	//	m_r = 0;
	//	m_uTilde = 0;
	//	m_u = 0;
	//	m_is_fullySet = false;

}

orbit::NonSingular2::~NonSingular2()
{
	// do nothing
}

orbit::NonSingular2::NonSingular2(Overwatch& tracer, double a, double xi, double eta,
	double h, double k, double lambda, bool is_normalized)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	m_a = a;
	m_xi = xi;
	m_eta = eta;
	m_h = h;
	m_k = k;
	m_lambda = lambda;
	m_is_normalized = is_normalized;

	if (m_is_normalized)//定义的时候自动更新anomaly
	{
		this->SyncAnomaly(tracer);
		if (tracer.Code() != ERR_NOERROR) return;
	}

	tracer.CheckOut(__func__);
}

void orbit::NonSingular2::Set(Overwatch& tracer, double a, double xi, double eta,
	double h, double k, double lambda, bool is_normalized)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);
	//定义的时候自动更新anomaly
	m_a = a;
	m_xi = xi;
	m_eta = eta;
	m_h = h;
	m_k = k;
	m_lambda = lambda;
	m_is_normalized = is_normalized;
	if (m_is_normalized)
	{
		this->SyncAnomaly(tracer);
		if (tracer.Code() != ERR_NOERROR) return;
	}
	tracer.CheckOut(__func__);
}

void orbit::NonSingular2::SyncAnomaly(Overwatch& tracer)
{
	// 参考《轨道设计基础》 P88
	// 同步生成uTidle和u，uTilde = m_Eow = E + omega + w	u = fow = f + omega + w
	// verify input elements
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)//没归一化做不了，得先归一化
	{
		std::string msg = "Normalize this nsg2 before sync anomaly!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else
	{

		if (m_strict_check && (m_a <= 0.0))
		{
			std::string msg = "Elliptic semi-major axis is non-positive!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return;
		}

		if (m_strict_check && (m_h > 1.0 || m_k > 1.0 || m_h < -1.0 || m_k < -1.0))
		{
			std::string msg = "Elliptic non-singular elements h or k is greater than 1!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return;
		}

		double e2 = m_xi * m_xi + m_eta * m_eta;//e^2
		if (e2 < 1.0)
		{
			m_uTilde = SolveGeneralizedKepler(m_xi, m_eta, m_lambda);
			// uTilde 可能已归入一周，lambda 则保留累计角度；须按 2π 而非 π 还原角度差。
			// 保留小角差原值，避免改变未跨圈算例的正常计算路径。
			double delta_angle = m_uTilde - m_lambda;
			while (delta_angle > PI) delta_angle -= 2.0 * PI;
			while (delta_angle < -PI) delta_angle += 2.0 * PI;

			double xces = m_xi * cos(m_uTilde) + m_eta * sin(m_uTilde);
			double ra = 1.0 - xces;
			m_r = ra * m_a;

			// 由偏近点经度与平经度之差求真经度修正，状态量仍使用归一化单位。
			double suuTilde = (delta_angle) * (1.0 - xces / (1.0 + sqrt(1.0 - e2))) / ra;

			// 现在 suuTilde 应该是一个很小的数 (sin值)，不会越界了
			// 增加一层保险是个好习惯
			if (suuTilde > 1.0) suuTilde = 1.0;
			if (suuTilde < -1.0) suuTilde = -1.0;

			m_u = m_uTilde + asin(suuTilde);
			m_is_fullySet = true;
		}
		else
		{
			std::string
				msg = "Anomaly (longitude) sync is skipped for non-ellipse conic e = " + std::to_string(sqrt(e2)) + ".";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
		}
	}
	tracer.CheckOut(__func__);
}

double orbit::NonSingular2::GetValue(Overwatch& tracer, int varId) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value;
	switch (varId)
	{
		case VAR_SEMIMAJORAXIS:
		case VAR_A_NSG2:
			value = m_a;
			break;
		case VAR_XI_NSG2:
			value = m_xi;
			break;
		case VAR_ETA_NSG2:
			value = m_eta;
			break;
		case VAR_H_NSG2:
			value = m_h;
			break;
		case VAR_K_NSG2:
			value = m_k;
			break;
		case VAR_LAMBDA_NSG2:
			value = m_lambda;
			break;
		case VAR_U_NSG2:
			if (m_is_fullySet)
			{
				value = m_u;
			}
			else
			{
				value = 99999;
				std::string msg = "u=f+w+Omega is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			}
			break;
		case VAR_UTILDE_NSG2:
			if (m_is_fullySet)
			{
				value = m_uTilde;
			}
			else
			{
				value = 99999;
				std::string msg = "u_tilde=E+w+Omega is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			}
			break;
		case VAR_R:
			if (m_is_fullySet)
			{
				value = m_r;
			}
			else
			{
				value = 99999;
				std::string msg = "Geocentric distance r is not set yet, 99999 is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			}
			break;
		default:
			value = 99999;
			std::string msg = "Wrong index, return 99999 in default!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
			break;
	}

	tracer.CheckOut(__func__);
	return value;
}

void orbit::NonSingular2::Normalize(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (m_is_normalized)
	{
		std::string msg = "NonSingular2 already norlamlized, this normalization is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		m_a /= uLength;
		m_lambda *= DEGREE_TO_RAD;
		m_is_normalized = true;

		// 顺便通过sync更新其余导出量
		SyncAnomaly(tracer);
		if (tracer.Code() != ERR_NOERROR) return;//层层return
	}
	tracer.CheckOut(__func__);
}

void orbit::NonSingular2::Restore(Overwatch& tracer, double uLength)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_is_normalized)
	{
		std::string msg = "NonSingular2 already restored, this restoration is no longer performed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		m_a *= uLength;
		m_lambda *= RAD_TO_DEGREE;
		m_is_normalized = false;

		if (m_is_fullySet)
		{
			m_r *= uLength;
			m_uTilde *= RAD_TO_DEGREE;
			m_u *= RAD_TO_DEGREE;
		}
	}

	tracer.CheckOut(__func__);
}

void orbit::NonSingular2::ToArray(double* nsg2_array) const
{
	nsg2_array[0] = m_a;
	nsg2_array[1] = m_xi;
	nsg2_array[2] = m_eta;
	nsg2_array[3] = m_h;
	nsg2_array[4] = m_k;
	nsg2_array[5] = m_lambda;
}

void orbit::NonSingular2::Show(std::ostream& os, int precision) const
{
	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_a;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_xi;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_eta;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_h;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_k;

	os << std::setiosflags(std::ios::left | std::ios::showpoint)
	   << std::setw(precision + 6) << std::setprecision(precision) << m_lambda;

	//os << std::endl;
}

void orbit::NonSingular2::ToCartesian(Overwatch& tracer, orbit::Cartesian& cart_class) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double cart_array[6];
	if (!this->ToCartesianArray(tracer, cart_array)) return;

	cart_class.Set(cart_array[0], cart_array[1], cart_array[2],
		cart_array[3], cart_array[4], cart_array[5], true);

	tracer.CheckOut(__func__);
}

void orbit::NonSingular2::ToCartesian(Overwatch& tracer, double* cart_array) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);
	if (!this->ToCartesianArray(tracer, cart_array)) return;
	tracer.CheckOut(__func__);
}

bool orbit::NonSingular2::ToCartesianArray(Overwatch& tracer, double* cart_array) const
{
// 参考《轨道设计基础》 P89

	if (!m_is_normalized)//没归一化做不了，得先归一化
	{
		std::string msg = "Normalize before converting to cartesian!";
		tracer.LogEvent(__FILE__, "ToCartesian", __LINE__, msg, EventType::ERROR);
		return false;
	}
	else
	{
		if (m_strict_check && (m_a <= 0.0))
		{
			std::string msg = "Elliptic semi-major axis is non-positive!";
			tracer.LogEvent(__FILE__, "ToCartesian", __LINE__, msg, EventType::ERROR);
			return false;
		}

		if (m_strict_check && (m_h > 1.0 || m_k > 1.0 || m_h < -1.0 || m_k < -1.0))
		{
			std::string msg = "Elliptic non-singular elements h or k is greater than 1!";
			tracer.LogEvent(__FILE__, "ToCartesian", __LINE__, msg, EventType::ERROR);
			return false;
		}

		if (m_strict_check && (m_xi >= 1.0 || m_eta >= 1.0))
		{
			std::string msg = "Elliptic non-singular element xi or eta is not less than 1!";
			tracer.LogEvent(__FILE__, "ToCartesian", __LINE__, msg, EventType::ERROR);
			return false;
		}

		if (m_is_fullySet)
		{
			Vector3d P, Q, r, v;//P*,Q*

			double cosu = cos(m_u);
			double sinu = sin(m_u);
			double hk = m_h * m_k;
			double h2 = m_h * m_h;
			double k2 = m_k * m_k;
			double cosi2 = sqrt(1.0 - h2 - k2);
			double p = m_a * (1.0 - m_xi * m_xi - m_eta * m_eta);

			P << 1.0 - 2.0 * k2, 2.0 * hk, -2.0 * m_k * cosi2;
			Q << 2.0 * hk, 1.0 - 2.0 * h2, 2.0 * m_h * cosi2;

			r = cosu * P + sinu * Q;
			r *= m_r;
			//默认GM=1
			v = -(sinu + m_eta) * P + (cosu + m_xi) * Q;
			v *= sqrt(1.0 / p);

			cart_array[0] = r[0];
			cart_array[1] = r[1];
			cart_array[2] = r[2];
			cart_array[3] = v[0];
			cart_array[4] = v[1];
			cart_array[5] = v[2];
		}
		else
		{
			std::string msg = "This nsg2 is not fully set yet, cannot covert to Cartestian!";
			tracer.LogEvent(__FILE__, "ToCartesian", __LINE__, msg, EventType::ERROR);
			return false;
		}
	}

	return true;
}

void orbit::NonSingular2::ToKeplerian(Overwatch& tracer, orbit::Keplerian& kep_class,
	double epsZeroEcc, double epsZeroInc) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double kep_array[6];
	if (!this->ToKeplerianArray(tracer, kep_array, epsZeroEcc, epsZeroInc)) return;

	kep_class.Set(kep_array[0], kep_array[1], kep_array[2],
		kep_array[3], kep_array[4], kep_array[5], true);

	tracer.CheckOut(__func__);
}

void orbit::NonSingular2::ToKeplerian(Overwatch& tracer, double* kep_array,
	double epsZeroEcc, double epsZeroInc) const
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);
	if (!this->ToKeplerianArray(tracer, kep_array, epsZeroEcc, epsZeroInc)) return;
	tracer.CheckOut(__func__);
}

bool orbit::NonSingular2::ToKeplerianArray(Overwatch& tracer, double* kep_array,
	double epsZeroEcc, double epsZeroInc) const
{
	// 参考《轨道设计基础》 P89

	if (!m_is_normalized)//没归一化做不了，得先归一化
	{
		std::string msg = "Normalize before converting to keplerian!";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::ERROR);
		return false;
	}

	if (m_strict_check && m_a <= 0.0)
	{
		std::string msg = "Elliptic semi-major axis is non-positive!";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::ERROR);
		return false;
	}

	if (m_strict_check && (m_h > 1.0 || m_k > 1.0 || m_h < -1.0 || m_k < -1.0))
	{
		std::string msg = "Elliptic non-singular elements h or k is greater than 1!";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::ERROR);
		return false;
	}

	//
	double e = sqrt(m_xi * m_xi + m_eta * m_eta);
	double sini2 = sqrt(m_h * m_h + m_k * m_k);//sin(i/2)
	double i = 2.0 * asin(sini2);

	if (e >= 1.0)
	{
		std::string
			msg = "NonSingular2 to Keplerian convertion finds a parabola/hyperbola of e = " + std::to_string(e) + ".";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::WARNING);
	}

	if (m_strict_check && (i < 0.0 || i >= PI))
	{

		/* in practice, only PI will be found */
		std::string msg = "Orbital inclination should be within [0,PI)!";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::ERROR);
		return false;

	}

	kep_array[0] = m_a;
	kep_array[1] = e;
	kep_array[2] = i;

	if ((e > epsZeroEcc) && (i > epsZeroInc))
	{
		// NEITHER small e NOR small i, process o, w and m normally
		double ow = fmod(atan2(m_eta / e, m_xi / e) + 2 * PI, 2 * PI);   // o+w 先转到(pi,3pi]再mod！
		kep_array[3] = fmod(atan2(m_k / sini2, m_h / sini2) + 2 * PI, 2 * PI);    // o
		kep_array[4] = fmod(ow - kep_array[3] + 2 * PI, 2 * PI);                  // w
		kep_array[5] = fmod(m_lambda - ow + 2 * PI, 2 * PI);                // lambda-(o+w)
	}
	else if ((e <= epsZeroEcc) && (i <= epsZeroInc))
	{
		// BOTH small e AND small i, set o=w=0
		kep_array[3] = 0.0;
		kep_array[4] = 0.0;
		kep_array[5] = m_lambda;
		std::string msg = "For zero eccentricity and zero inclination, countermeasures are taken.";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::WARNING);
	}
	else if ((e <= epsZeroEcc) && (i > epsZeroInc))
	{
		// small e but NOT small i, process o normally but w=0
		kep_array[3] = fmod(atan2(m_k / sini2, m_h / sini2) + 2 * PI, 2 * PI);
		kep_array[4] = 0.0;
		kep_array[5] = fmod(m_lambda - kep_array[3] + 2 * PI, 2 * PI);
		std::string msg = "For zero eccentricity, countermeasures are taken.";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::WARNING);
	}
	else if ((e > epsZeroEcc) && (i <= epsZeroInc))
	{
		// NOT small e but small i, process w normally but o=0
		kep_array[3] = 0.0;
		kep_array[4] = fmod(atan2(m_eta / e, m_xi / e) + 2 * PI, 2 * PI);
		kep_array[5] = fmod(m_lambda - kep_array[4] + 2 * PI, 2 * PI);
		std::string msg = "For zero inclination, countermeasures are taken.";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::WARNING);
	}
	else
	{
		std::string msg = "Unaccessible logic detected!";
		tracer.LogEvent(__FILE__, "ToKeplerian", __LINE__, msg, EventType::ERROR);
		return false;
	}

	return true;
}
