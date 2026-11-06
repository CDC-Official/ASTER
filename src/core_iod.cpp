//
//  core_iod.cpp
//  Test_dylib
//
//  Created by 陈酊存 on 2023/8/30.
//

#include "orbit_geometry.hpp"
#include "math_statistics.hpp"
#include "core_iod.hpp"
#include "data_gravityfield.hpp"
#include <sstream>

using namespace constants::Scale;
using namespace constants::Conversion;

using namespace iod;

/**
 * @brief 输出 GL 初轨确定的动力学、迭代、约束和观测误差设置。
 * @param[in,out] output 接收 IOD settings 分区的输出流；名称宽度为 40，参数按键值对排列。
 * @details 读取 dynamicsModel、iterMethod、constraintType、constraintEcc、initValue.type、
 * iterCountMax、iterTOL 和 obs.sigma_RaDec。动力学按五种 IOD 宏定义显示，
 * iterTOL 为六个归一化状态改正分量的最大绝对值阈值；sigma_RaDec 为东西、南北
 * 两个独立高斯误差共同的标准差，由 rad 换算为 arcsec 后显示。
 * 本函数不修改设置及输出流的格式设置，也不读取文件或执行定轨计算。
 */
void iod::IodJob::Display(std::ostream& output) const
{
	// 按实际宏定义取得显示名称，模型编号与输入关键词中的取值一致。
	const char* dynamics_name = "Unknown";

	switch (dynamicsModel)
	{
		case IOD_DYNAMICS_TWOBODY: dynamics_name = "Two-body"; break;
		case IOD_DYNAMICS_J2: dynamics_name = "J2"; break;
		case IOD_DYNAMICS_J2_SOLAR: dynamics_name = "J2 + Sun"; break;
		case IOD_DYNAMICS_J2_LUNAR: dynamics_name = "J2 + Moon"; break;
		case IOD_DYNAMICS_J2_SOLAR_LUNAR: dynamics_name = "J2 + Sun + Moon"; break;
	}

	const char* method_name = "Unknown";

	switch (iterMethod)
	{
		case IOD_ITERATION_SIMPLE: method_name = "Simple iteration"; break;
		case IOD_ITERATION_NEWTON: method_name = "Differential correction"; break;
	}

	const char* constraint_name = "Unknown";

	switch (constraintType)
	{
		case IOD_CONSTRAINT_NONE: constraint_name = "None"; break;
		case IOD_CONSTRAINT_ORIENTATION: constraint_name = "Orientation"; break;
		case IOD_CONSTRAINT_CIRCULAR: constraint_name = "Circular"; break;
		case IOD_CONSTRAINT_ORIENTATION_CIRCULAR: constraint_name = "Orientation + Circular"; break;
		case IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY: constraint_name = "Orientation + Eccentricity"; break;
	}

	const char* initial_state_name = "Unknown";

	switch (initValue.type)
	{
		case IOD_INITIAL_VALUE_DEFAULT: initial_state_name = "Default"; break;
		case IOD_INITIAL_VALUE_SPECIFIED: initial_state_name = "Specified"; break;
	}

	// 在独立文字流中设置格式，使调用方的输出流保持原有精度和对齐方式。
	const int label_width = 40;
	std::ostringstream settings;
	settings << "\n===================================== IOD settings =====================================\n\n"
			 << std::left << std::setw(label_width) << "- Dynamics model" << ": " << dynamicsModel << " — " << dynamics_name
			 << '\n' << std::setw(label_width) << "- Iteration method" << ": " << method_name
			 << '\n' << std::setw(label_width) << "- Constraints" << ": " << constraint_name;

	// 指定偏心率属于当前约束，在同一行给出实际参与方程的数值。
	if (constraintType == IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY)
	{
		settings << " (e = " << std::defaultfloat << std::setprecision(8) << constraintEcc << ')';
	}

	settings << '\n' << std::setw(label_width) << "- Initial state" << ": " << initial_state_name
			 << '\n' << std::setw(label_width) << "- Max iterations" << ": " << iterCountMax
			 << '\n' << std::setw(label_width) << "- State correction limit [norm]" << ": "
			 << std::scientific << std::setprecision(2) << iterTOL
			 << "\n\n" << std::setw(label_width) << "- Observation error standard deviation" << ": "
			 << std::defaultfloat << std::setprecision(8) << obs.sigma_RaDec / ARCSEC_TO_RAD << " arcsec (East / North)\n";
	output << settings.str();
}

/**
 * @brief 输出 IOD 最终轨道及测角残差的均方根，可附带来源迭代信息。
 * @param[in,out] output 接收结果文字的输出流；名称宽度为 40，状态按两行三分量排列。
 * @param[in] include_iteration_info 为 true 时显示 flag、迭代次数和最大改正量；
 *     为 false 时仅显示轨道与 RMS，适用于 GLML 边界处理后的最终轨道。
 * @details 读取 cart、kep、resRa、resDec 和 RMS_angle；东西、南北 RMS 由保存的
 * 残差分别求均方根。附带迭代信息时，读取 flag 和 iterations 的最后一条记录。
 * Cartesian 从 ae、ae/UNITTIME 换算为 km、km/s，
 * Keplerian 从 ae、rad 换算为 km、deg，最大改正量保留归一化单位。
 * 东西、南北 RMS 分别占一行，每行依次显示 rad 和 arcsec；合成 RMS 也显示两种单位。
 * 弧度值采用科学计数法，尾数保留两位小数；角秒值按定点格式保留两位小数。
 * 轨道均属于调用方的参考历元和 GCRS，不重新传播或计算残差。
 * 本函数不输出分区标题、任务耗时和内存；不修改结果及输出流的格式设置。
 */
void iod::IodSolution::PrintSummary(std::ostream& output, bool include_iteration_info) const
{
	// 将已保存的归一化状态复制成显示用数组，结果对象仍保留计算时的单位。
	double cart_values[6], kep_values[6];
	cart.ToArray(cart_values);
	kep.ToArray(kep_values);

	for (int j = 0; j < 6; ++j)
	{
		cart_values[j] *= j < 3 ? UNIT_LENGTH_M / 1000.0 : UNIT_LENGTH_M / (1000.0 * UNIT_TIME_S);
		kep_values[j] *= j == 0 ? UNIT_LENGTH_M / 1000.0 : (j >= 2 ? RAD_TO_DEGREE : 1.0);
	}

	const int label_width = 40;
	const int value_width = 18;
	const std::string continuation(label_width + 2, ' ');
	std::ostringstream summary;
	summary << std::left;

	// GLML 的边界处理不产生新的迭代记录，调用方可单独说明来源状态。
	if (include_iteration_info)
	{
		summary << std::setw(label_width) << "- Iteration Status" << ": ";

		if (flag > 0)
		{
			summary << "CONVERGED";
		}
		else if (flag == -1)
		{
			summary << "NOT CONVERGED — ITERATION LIMIT";
		}
		else
		{
			summary << "NUMERICAL FAILURE";
		}

		// 迭代次数取最后一条记录的轮次；起始记录尚无改正时显示 N/A。
		summary << '\n' << std::setw(label_width) << "- Iterations" << ": ";

		if (!iterations.empty())
		{
			summary << iterations.back().iteration;
		}
		else
		{
			summary << "N/A";
		}

		summary << '\n' << std::setw(label_width) << "- Max absolute correction [norm]" << ": ";

		if (!iterations.empty() && std::isfinite(iterations.back().max_dr0v0_norm))
		{
			summary << std::scientific << std::setprecision(6) << iterations.back().max_dr0v0_norm;
		}
		else
		{
			summary << "N/A";
		}

		summary << "\n\n";
	}

	// 两种状态各排成两行，续行的左方括号与第一行严格对齐。
	summary << std::setw(label_width) << "- Final Cartesian [km, km/s]" << ": [ "
			<< std::right << std::fixed << std::setprecision(6)
			<< std::setw(value_width) << cart_values[0] << ", " << std::setw(value_width) << cart_values[1] << ", "
			<< std::setw(value_width) << cart_values[2] << " ]\n"
			<< continuation << "[ " << std::setw(value_width) << cart_values[3] << ", "
			<< std::setw(value_width) << cart_values[4] << ", " << std::setw(value_width) << cart_values[5] << " ]\n"
			<< '\n' << std::left << std::setw(label_width) << "- Final Keplerian [km, deg]" << ": [ " << std::right
			<< std::setw(value_width) << kep_values[0] << ", " << std::setw(value_width) << kep_values[1] << ", "
			<< std::setw(value_width) << kep_values[2] << " ]\n"
			<< continuation << "[ " << std::setw(value_width) << kep_values[3] << ", "
			<< std::setw(value_width) << kep_values[4] << ", " << std::setw(value_width) << kep_values[5] << " ]\n";

	// 两个方向各占一行，依次显示弧度值及换算后的角秒值。
	summary << '\n' << std::left;

	if (resRa.size() > 0 && resRa.size() == resDec.size())
	{
		const double rms_east = std::sqrt(resRa.square().mean());
		const double rms_north = std::sqrt(resDec.square().mean());

		summary << std::setw(label_width) << "- RMS East [rad / arcsec]" << ": "
				<< std::scientific << std::setprecision(2) << rms_east << " / "
				<< std::fixed << std::setprecision(2) << rms_east / ARCSEC_TO_RAD
				<< '\n' << std::setw(label_width) << "- RMS North [rad / arcsec]" << ": "
				<< std::scientific << std::setprecision(2) << rms_north << " / "
				<< std::fixed << std::setprecision(2) << rms_north / ARCSEC_TO_RAD;
	}
	else
	{
		summary << std::setw(label_width) << "- RMS East [rad / arcsec]" << ": N/A / N/A\n"
				<< std::setw(label_width) << "- RMS North [rad / arcsec]" << ": N/A / N/A";
	}

	summary << '\n' << std::setw(label_width) << "- RMS [rad / arcsec]" << ": "
			<< std::scientific << std::setprecision(2) << RMS_angle << " / "
			<< std::fixed << std::setprecision(2) << RMS_angle / ARCSEC_TO_RAD << '\n';
	output << summary.str();
}

/**
 * @brief 由参考状态和相对时间计算GL的拉格朗日传播系数。
 * @param[in] dynamics_model 取IOD_DYNAMICS_*，选择二体、J2及解析日月摄动组合。
 * @param[in] epoch_of_IOD 轨道参考历元，单位UTC MJD；用于解析日月位置计算。
 * @param[in] tau 相对参考历元的时间差，除以Scale::UNIT_TIME_S。
 * @param[in] X 六维GCRS参考状态[x,y,z,vx,vy,vz]，位置除以Scale::UNIT_LENGTH_M，
 *     速度乘Scale::UNIT_TIME_S/Scale::UNIT_LENGTH_M；本函数不修改该数组。
 * @param[out] LagCoeff 四个归一化传播系数[f,g,fz,gz]；x/y使用f、g，z使用fz、gz。
 * @details 二体部分保留至时间六次项，J2摄动部分保留至四次项。
 * gz的J2附加三次项为(3J2/2)*tau³*(-2/r0⁵)/6，按位置三阶导数除以3!取得。
 * 日月距离在对应观测历元使用解析位置计算，所有系数与输入的归一化尺度配套。
 */
void Lagrange_Coefficients(int dynamics_model, double epoch_of_IOD, double tau, double* X, double* LagCoeff)
{
	// 从参考状态构造径向距离倒数、位置速度内积和时间幂，供所选模型展开使用。
	double J, u1, u2, u3, u4, u5, u6, u7, u8, u9, sigma, sigma2, v02;
	u1 = 1.0 / sqrt(X[0] * X[0] + X[1] * X[1] + X[2] * X[2]);
	u2 = u1 * u1;
	u3 = u2 * u1;
	u4 = u3 * u1;
	u5 = u4 * u1;
	u6 = u5 * u1;
	u7 = u6 * u1;
	u8 = u7 * u1;
	u9 = u8 * u1;
	sigma = X[0] * X[3] + X[1] * X[4] + X[2] * X[5];
	sigma2 = sigma * sigma;
	v02 = X[3] * X[3] + X[4] * X[4] + X[5] * X[5];    // v

	// J2 与 GGM03C 的 GM、参考半径及 C20 配套，换算到当前积分尺度。
	J = 1.5 * constants::GGM03C::NORMALIZED_J2;

	double tau2, tau3, tau4, tau5, tau6;
	tau2 = tau * tau;
	tau3 = tau2 * tau;
	tau4 = tau3 * tau;
	tau5 = tau4 * tau;
	tau6 = tau5 * tau;

	double z0 = X[2], dz0 = X[5];

	switch (dynamics_model)
	{
		case IOD_DYNAMICS_TWOBODY:
		{
			// F
			LagCoeff[0] =
				1 + tau2 / 2.0 * (-u3) +
				tau3 / 6.0 * (3 * u5 * sigma) +
				tau4 / 24.0 * (u5 * (3 * v02 - 2 * u1 - 15 * u2 * sigma2)) +
				tau5 / 120.0 * u7 * (15 * sigma * (-3 * v02 + 2 * u1 + 7 * u2 * sigma2)) +
				tau6 / 720.0 * u7 * (u2 * sigma2 * (630 * v02 - 420 * u1 - 945 * u2 * sigma2)
									 - (22 * u2 - 66 * u1 * v02 + 45 * pow(v02, 2)));
			// G
			LagCoeff[1] =
				tau +
				tau3 / 6.0 * (-u3) +
				tau4 / 24.0 * (6 * u5 * sigma) +
				tau5 / 120.0 * u5 * (9 * v02 - 8 * u1 - 45 * u2 * sigma2) +
				tau6 / 720.0 * u7 * (30 * sigma * (-6 * v02 + 5 * u1 + 14 * u2 * sigma2));
			// Fz 仅二体下是不用区分的
			LagCoeff[2] = LagCoeff[0];
			// Gz
			LagCoeff[3] = LagCoeff[1];

			break;
		}
		case IOD_DYNAMICS_J2:
		{
			//  F
			LagCoeff[0] =
				1 + tau2 / 2.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau3 / 6.0 *
				(3 * u5 * sigma +
				 J * (5 * (u7 - 7 * u9 * pow(z0, 2)) * sigma + 10 * u7 * z0 * dz0)) +
				tau4 / 24.0 *
				(u5 * (3 * v02 - 2 * u1 - 15 * u2 * sigma2) +
				 J * (6 * u8 * (4 * u2 * pow(z0, 2) - 1) -
					  5 * u7 * (7 * u2 * pow(z0, 2) - 1) * v02 +
					  10 * u7 * pow(dz0, 2) +
					  35 * u9 * (9 * u2 * pow(z0, 2) - 1) * sigma2 -
					  140 * u9 * sigma * z0 * dz0)) +
				tau5 / 120.0 * u7 *
				(15 * sigma * (-3 * v02 + 2 * u1 + 7 * u2 * sigma2)) +
				tau6 / 720.0 * u7 *
				(u2 * sigma2 *
				 (630 * v02 - 420 * u1 - 945 * u2 * sigma2) -
				 (22 * u2 - 66 * u1 * v02 + 45 * pow(v02, 2)));
			// G
			LagCoeff[1] =
				tau + tau3 / 6.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau4 / 24.0 *
				(6 * u5 * sigma + J * (20 * u7 * z0 * dz0 -
									   10 * u7 * (7 * u2 * pow(z0, 2) - 1) * sigma)) +
				tau5 / 120.0 * u5 * (9 * v02 - 8 * u1 - 45 * u2 * sigma2) +
				tau6 / 720.0 * u7 *
				(30 * sigma * (-6 * v02 + 5 * u1 + 14 * u2 * sigma2));
			// Fz
			LagCoeff[2] =
				LagCoeff[0] +
				J * (tau2 / 2.0 * (-2 * u5) + tau3 / 6.0 * (10 * u7 * sigma) +
					 tau4 / 24.0 * u7 * (10 * v02 - 6 * u1 - 70 * u2 * sigma2));
			// Gz
			LagCoeff[3] = LagCoeff[1] +
						  J * (tau3 / 6.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));

			break;
		}
		case IOD_DYNAMICS_J2_SOLAR:
		{
			// F
			LagCoeff[0] =
				1 + tau2 / 2.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau3 / 6.0 *
				(3 * u5 * sigma +
				 J * (5 * (u7 - 7 * u9 * pow(z0, 2)) * sigma + 10 * u7 * z0 * dz0)) +
				tau4 / 24.0 *
				(u5 * (3 * v02 - 2 * u1 - 15 * u2 * sigma2) +
				 J * (6 * u8 * (4 * u2 * pow(z0, 2) - 1) -
					  5 * u7 * (7 * u2 * pow(z0, 2) - 1) * v02 +
					  10 * u7 * pow(dz0, 2) +
					  35 * u9 * (9 * u2 * pow(z0, 2) - 1) * sigma2 -
					  140 * u9 * sigma * z0 * dz0)) +
				tau5 / 120.0 * u7 *
				(15 * sigma * (-3 * v02 + 2 * u1 + 7 * u2 * sigma2)) +
				tau6 / 720.0 * u7 *
				(u2 * sigma2 *
				 (630 * v02 - 420 * u1 - 945 * u2 * sigma2) -
				 (22 * u2 - 66 * u1 * v02 + 45 * pow(v02, 2)));
			// G
			LagCoeff[1] =
				tau + tau3 / 6.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau4 / 24.0 *
				(6 * u5 * sigma + J * (20 * u7 * z0 * dz0 -
									   10 * u7 * (7 * u2 * pow(z0, 2) - 1) * sigma)) +
				tau5 / 120.0 * u5 * (9 * v02 - 8 * u1 - 45 * u2 * sigma2) +
				tau6 / 720.0 * u7 *
				(30 * sigma * (-6 * v02 + 5 * u1 + 14 * u2 * sigma2));
			// Fz
			LagCoeff[2] =
				LagCoeff[0] +
				J * (tau2 / 2.0 * (-2 * u5) + tau3 / 6.0 * (10 * u7 * sigma) +
					 tau4 / 24.0 * u7 * (10 * v02 - 6 * u1 - 70 * u2 * sigma2));
			// Gz
			LagCoeff[3] = LagCoeff[1] +
						  J * (tau3 / 6.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));

			//这里开始考虑日摄动，要计算每一时刻的日位置
			double MJD = epoch_of_IOD + tau * NORMALIZED_TIME_TO_DAY;  //对应时刻的MJD
			double r_sun[3], u3_sun, r0_sun;
			Cart_MeanSun(MJD, r_sun);  //该函数返回归一化单位
			r0_sun = sqrt(r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]);
			u3_sun = 1.0 / r0_sun / r0_sun / r0_sun;
			LagCoeff[0] +=
				tau2 / 2.0 * (-NORMALIZED_SUN_GM_TT * u3_sun)
				+ tau4 / 24.0 * (u3 * NORMALIZED_SUN_GM_TT * u3_sun);    // F
			LagCoeff[1] +=
				tau3 / 6.0 * (-NORMALIZED_SUN_GM_TT * u3_sun);// G

			break;
		}
		case IOD_DYNAMICS_J2_LUNAR:
		{
			// F
			LagCoeff[0] =
				1 + tau2 / 2.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau3 / 6.0 *
				(3 * u5 * sigma +
				 J * (5 * (u7 - 7 * u9 * pow(z0, 2)) * sigma + 10 * u7 * z0 * dz0)) +
				tau4 / 24.0 *
				(u5 * (3 * v02 - 2 * u1 - 15 * u2 * sigma2) +
				 J * (6 * u8 * (4 * u2 * pow(z0, 2) - 1) -
					  5 * u7 * (7 * u2 * pow(z0, 2) - 1) * v02 +
					  10 * u7 * pow(dz0, 2) +
					  35 * u9 * (9 * u2 * pow(z0, 2) - 1) * sigma2 -
					  140 * u9 * sigma * z0 * dz0)) +
				tau5 / 120.0 * u7 *
				(15 * sigma * (-3 * v02 + 2 * u1 + 7 * u2 * sigma2)) +
				tau6 / 720.0 * u7 *
				(u2 * sigma2 *
				 (630 * v02 - 420 * u1 - 945 * u2 * sigma2) -
				 (22 * u2 - 66 * u1 * v02 + 45 * pow(v02, 2)));
			// G
			LagCoeff[1] =
				tau + tau3 / 6.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau4 / 24.0 *
				(6 * u5 * sigma + J * (20 * u7 * z0 * dz0 -
									   10 * u7 * (7 * u2 * pow(z0, 2) - 1) * sigma)) +
				tau5 / 120.0 * u5 * (9 * v02 - 8 * u1 - 45 * u2 * sigma2) +
				tau6 / 720.0 * u7 *
				(30 * sigma * (-6 * v02 + 5 * u1 + 14 * u2 * sigma2));
			// Fz
			LagCoeff[2] =
				LagCoeff[0] +
				J * (tau2 / 2.0 * (-2 * u5) + tau3 / 6.0 * (10 * u7 * sigma) +
					 tau4 / 24.0 * u7 * (10 * v02 - 6 * u1 - 70 * u2 * sigma2));
			// Gz
			LagCoeff[3] = LagCoeff[1] +
						  J * (tau3 / 6.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));

			//这里开始考虑日摄动，要计算每一时刻的日位置
			double MJD = epoch_of_IOD + tau * NORMALIZED_TIME_TO_DAY;  //对应时刻的MJD
			double r_moon[3], u3_moon, r0_moon;
			Cart_MeanMoon(MJD, r_moon);     //该函数返回归一化单位
			r0_moon = sqrt(r_moon[0] * r_moon[0] + r_moon[1] * r_moon[1] + r_moon[2] * r_moon[2]);
			u3_moon = 1.0 / r0_moon / r0_moon / r0_moon;

			//日月分开加，便于剥离；只有FG的低次幂有三体摄动
			LagCoeff[0] += tau2 / 2.0 * (-NORMALIZED_MOON_GM_TT * u3_moon)
						   + tau4 / 24.0 * (u3 * NORMALIZED_MOON_GM_TT * u3_moon);    // F
			LagCoeff[1] += tau3 / 6.0
						   * (-NORMALIZED_MOON_GM_TT * u3_moon);// G

			break;
		}
		case IOD_DYNAMICS_J2_SOLAR_LUNAR:
		{
			// F
			LagCoeff[0] =
				1 + tau2 / 2.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau3 / 6.0 *
				(3 * u5 * sigma +
				 J * (5 * (u7 - 7 * u9 * pow(z0, 2)) * sigma + 10 * u7 * z0 * dz0)) +
				tau4 / 24.0 *
				(u5 * (3 * v02 - 2 * u1 - 15 * u2 * sigma2) +
				 J * (6 * u8 * (4 * u2 * pow(z0, 2) - 1) -
					  5 * u7 * (7 * u2 * pow(z0, 2) - 1) * v02 +
					  10 * u7 * pow(dz0, 2) +
					  35 * u9 * (9 * u2 * pow(z0, 2) - 1) * sigma2 -
					  140 * u9 * sigma * z0 * dz0)) +
				tau5 / 120.0 * u7 *
				(15 * sigma * (-3 * v02 + 2 * u1 + 7 * u2 * sigma2)) +
				tau6 / 720.0 * u7 *
				(u2 * sigma2 *
				 (630 * v02 - 420 * u1 - 945 * u2 * sigma2) -
				 (22 * u2 - 66 * u1 * v02 + 45 * pow(v02, 2)));
			// G
			LagCoeff[1] =
				tau + tau3 / 6.0 * (-u3 + J * (5 * u7 * pow(z0, 2) - u5)) +
				tau4 / 24.0 *
				(6 * u5 * sigma + J * (20 * u7 * z0 * dz0 -
									   10 * u7 * (7 * u2 * pow(z0, 2) - 1) * sigma)) +
				tau5 / 120.0 * u5 * (9 * v02 - 8 * u1 - 45 * u2 * sigma2) +
				tau6 / 720.0 * u7 *
				(30 * sigma * (-6 * v02 + 5 * u1 + 14 * u2 * sigma2));
			// Fz
			LagCoeff[2] =
				LagCoeff[0] +
				J * (tau2 / 2.0 * (-2 * u5) + tau3 / 6.0 * (10 * u7 * sigma) +
					 tau4 / 24.0 * u7 * (10 * v02 - 6 * u1 - 70 * u2 * sigma2));
			// Gz
			LagCoeff[3] = LagCoeff[1] +
						  J * (tau3 / 6.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));
			//这里开始考虑日月摄动，要计算每一时刻的日月位置;日月分开
			//日
			double MJD = epoch_of_IOD + tau * NORMALIZED_TIME_TO_DAY;  //对应时刻的MJD
			double r_sun[3], u3_sun, r0_sun;
			Cart_MeanSun(MJD, r_sun);  //该函数返回归一化单位
			r0_sun = sqrt(r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]);
			u3_sun = 1.0 / r0_sun / r0_sun / r0_sun;
			LagCoeff[0] +=
				tau2 / 2.0 * (-NORMALIZED_SUN_GM_TT * u3_sun)
				+ tau4 / 24.0 * (u3 * NORMALIZED_SUN_GM_TT * u3_sun);    // F
			LagCoeff[1] +=
				tau3 / 6.0
				* (-NORMALIZED_SUN_GM_TT * u3_sun);                                                        // G
			//月
			double r_moon[3], u3_moon, r0_moon;
			Cart_MeanMoon(MJD, r_moon);     //该函数返回归一化单位
			r0_moon = sqrt(r_moon[0] * r_moon[0] + r_moon[1] * r_moon[1] + r_moon[2] * r_moon[2]);
			u3_moon = 1.0 / r0_moon / r0_moon / r0_moon;

			//日月分开加，便于剥离；只有FG的低次幂有三体摄动
			LagCoeff[0] += tau2 / 2.0 * (-NORMALIZED_MOON_GM_TT * u3_moon)
						   + tau4 / 24.0 * (u3 * NORMALIZED_MOON_GM_TT * u3_moon);    // F
			LagCoeff[1] += tau3 / 6.0
						   * (-NORMALIZED_MOON_GM_TT * u3_moon);  // G

			break;
		}
		default:
		{
			std::cerr << "Wrong dynamics model, two body model is used in default!" << std::endl;
			// F
			LagCoeff[0] =
				1 + tau2 / 2.0 * (-u3) +
				tau3 / 6.0 * (3 * u5 * sigma) +
				tau4 / 24.0 * (u5 * (3 * v02 - 2 * u1 - 15 * u2 * sigma2)) +
				tau5 / 120.0 * u7 * (15 * sigma * (-3 * v02 + 2 * u1 + 7 * u2 * sigma2)) +
				tau6 / 720.0 * u7 * (u2 * sigma2 * (630 * v02 - 420 * u1 - 945 * u2 * sigma2)
									 - (22 * u2 - 66 * u1 * v02 + 45 * pow(v02, 2)));
			// G
			LagCoeff[1] =
				tau +
				tau3 / 6.0 * (-u3) +
				tau4 / 24.0 * (6 * u5 * sigma) +
				tau5 / 120.0 * u5 * (9 * v02 - 8 * u1 - 45 * u2 * sigma2) +
				tau6 / 720.0 * u7 * (30 * sigma * (-6 * v02 + 5 * u1 + 14 * u2 * sigma2));
			// Fz 仅二体下是不用区分的
			LagCoeff[2] = LagCoeff[0];
			// Gz
			LagCoeff[3] = LagCoeff[1];

			break;
		}
	}
}


void PartialFG(double tau, double x0, double y0, double z0, double* PFG)
{
	// PFG=partialF&G F1 F2 F3 G1 G2 G3
	double F_common =
		1.5 * pow(tau, 2) * pow(pow(x0, 2) + pow(y0, 2) + pow(z0, 2), -2.5);
	double G_common =
		0.5 * pow(tau, 3) * pow(pow(x0, 2) + pow(y0, 2) + pow(z0, 2), -2.5);

	PFG[0] = x0 * F_common;     // F1
	PFG[1] = y0 * F_common;     // F2
	PFG[2] = z0 * F_common;     // F3
	PFG[3] = x0 * G_common;     // G1
	PFG[4] = y0 * G_common;     // G2
	PFG[5] = z0 * G_common;     // G3
}

