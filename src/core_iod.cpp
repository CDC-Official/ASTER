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
 * @brief 计算各观测历元所需的日月地心距离倒数三次方。
 * @param[in,out] tracer 记录模型、历元、相对时间或第三体距离中的错误。
 * @param[in] dynamics_model 取 IOD_DYNAMICS_*，仅计算该模型启用的太阳或月球。
 * @param[in] epoch_mjd_utc 轨道参考历元，单位 UTC MJD。
 * @param[in] tau N 个观测相对参考历元的时间差，以 UNITTIME 归一化；顺序与观测一致。
 * @param[out] distances N 条第三体距离记录，距离以 ae 归一化；未启用的第三体填零。
 * @details 在 epoch_mjd_utc+tau*NORMALIZED_TIME_TO_DAY 处调用平均日月位置模型，
 * 使用其地心赤道坐标的长度。UTC 历元直接代入解析模型的近似与 IOD 传播系数一致。
 * 各条记录不依赖待估轨道，可在相同历元、相对时间和动力学模型的候选间只读共享。
 */
void Prepare_ThirdBodyDistances(Overwatch &tracer, int dynamics_model, double epoch_mjd_utc,
                               const ArrayXd &tau, std::vector<IodThirdBodyDistances> &distances)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 每条观测使用与传播系数完全相同的历元表达式，输入应为有限值和有效模型。
	if (dynamics_model < IOD_DYNAMICS_TWOBODY || dynamics_model > IOD_DYNAMICS_J2_SOLAR_LUNAR ||
		!std::isfinite(epoch_mjd_utc) || !tau.allFinite())
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid IOD model, reference epoch or normalized time offsets for third-body distances.", EventType::ERROR);
		return;
	}

	const bool use_sun = dynamics_model == IOD_DYNAMICS_J2_SOLAR || dynamics_model == IOD_DYNAMICS_J2_SOLAR_LUNAR;
	const bool use_moon = dynamics_model == IOD_DYNAMICS_J2_LUNAR || dynamics_model == IOD_DYNAMICS_J2_SOLAR_LUNAR;
	distances.assign(tau.size(), IodThirdBodyDistances{});

	for (Eigen::Index k = 0; k < tau.size(); ++k)
	{
		const double mjd = epoch_mjd_utc + tau(k) * NORMALIZED_TIME_TO_DAY;

		if (!std::isfinite(mjd))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite third-body epoch in observation " + std::to_string(k + 1), EventType::ERROR);
			return;
		}

		// 太阳与月球分别计算；距离的平方根和连续除法按传播函数的顺序执行。
		if (use_sun)
		{
			double position[3];
			Cart_MeanSun(mjd, position);
			const double radius = sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
			distances[k].sun_inverse_radius_cubed = 1.0 / radius / radius / radius;
		}

		if (use_moon)
		{
			double position[3];
			Cart_MeanMoon(mjd, position);
			const double radius = sqrt(position[0] * position[0] + position[1] * position[1] + position[2] * position[2]);
			distances[k].moon_inverse_radius_cubed = 1.0 / radius / radius / radius;
		}

		// 第三体距离必须有限且为正；未启用的分量保留零值。
		if ((use_sun && (!std::isfinite(distances[k].sun_inverse_radius_cubed) || distances[k].sun_inverse_radius_cubed <= 0.0)) ||
			(use_moon && (!std::isfinite(distances[k].moon_inverse_radius_cubed) || distances[k].moon_inverse_radius_cubed <= 0.0)))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid third-body distance in observation " + std::to_string(k + 1), EventType::ERROR);
			return;
		}
	}

	tracer.CheckOut(__func__);
}

/**
 * @brief 由参考状态和相对时间计算 GL 的拉格朗日传播系数。
 * @param[in] dynamics_model 取 IOD_DYNAMICS_*，选择二体、J2 及平均日月摄动组合。
 * @param[in] epoch_of_IOD 轨道参考历元，单位 UTC MJD。
 * @param[in] tau 观测相对参考历元的时间差，以 UNITTIME 归一化。
 * @param[in] X 六维 GCRS 参考状态 [x,y,z,vx,vy,vz]，位置以 ae、速度以 ae/UNITTIME 归一化。
 * @param[out] LagCoeff 四个传播系数 [f,g,fz,gz]，与归一化状态和时间配套使用。
 * @param[in] third_body_distances 可选的本条观测日月距离记录，须由相同模型、参考历元和 tau 准备。
 *     传入 nullptr 时现场计算所需的日月位置；二体及 J2 模型不读取此记录。
 * @details 以参考状态展开 f/g，x/y 与 z 分量分别使用 [f,g] 和 [fz,gz]。
 * 复用第三体距离仅替代与待估轨道无关的解析位置及距离计算，后续浮点表达式保持不变。
 */
void Lagrange_Coefficients(int dynamics_model, double epoch_of_IOD, double tau, double* X, double* LagCoeff,
                           const IodThirdBodyDistances* third_body_distances)
{
	// 从当前参考状态构造径向倒数、位置速度内积和时间幂，供所选动力学展开使用。
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
						  J * (tau3 / 2.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));

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
						  J * (tau3 / 2.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));

			// 读取本次观测的太阳距离，用于 f/g 中的太阳摄动项。
			double MJD = epoch_of_IOD + tau * NORMALIZED_TIME_TO_DAY;  // 观测历元，UTC MJD。
			double u3_sun;

			if (third_body_distances != nullptr)
			{
				u3_sun = third_body_distances->sun_inverse_radius_cubed;
			}
			else
			{
				double r_sun[3];
				Cart_MeanSun(MJD, r_sun);
				const double r0_sun = sqrt(r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]);
				u3_sun = 1.0 / r0_sun / r0_sun / r0_sun;
			}

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
						  J * (tau3 / 2.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));

			// 读取本次观测的月球距离，用于 f/g 中的月球摄动项。
			double MJD = epoch_of_IOD + tau * NORMALIZED_TIME_TO_DAY;  // 观测历元，UTC MJD。
			double u3_moon;

			if (third_body_distances != nullptr)
			{
				u3_moon = third_body_distances->moon_inverse_radius_cubed;
			}
			else
			{
				double r_moon[3];
				Cart_MeanMoon(MJD, r_moon);
				const double r0_moon = sqrt(r_moon[0] * r_moon[0] + r_moon[1] * r_moon[1] + r_moon[2] * r_moon[2]);
				u3_moon = 1.0 / r0_moon / r0_moon / r0_moon;
			}


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
						  J * (tau3 / 2.0 * (-2 * u5) + tau4 / 24.0 * (20 * u7 * sigma));
			// 太阳、月球的距离按观测历元复用，两个摄动贡献分别加入 f/g。
			double MJD = epoch_of_IOD + tau * NORMALIZED_TIME_TO_DAY;  // 观测历元，UTC MJD。
			double u3_sun;

			if (third_body_distances != nullptr)
			{
				u3_sun = third_body_distances->sun_inverse_radius_cubed;
			}
			else
			{
				double r_sun[3];
				Cart_MeanSun(MJD, r_sun);
				const double r0_sun = sqrt(r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]);
				u3_sun = 1.0 / r0_sun / r0_sun / r0_sun;
			}

			LagCoeff[0] +=
				tau2 / 2.0 * (-NORMALIZED_SUN_GM_TT * u3_sun)
				+ tau4 / 24.0 * (u3 * NORMALIZED_SUN_GM_TT * u3_sun);    // F
			LagCoeff[1] +=
				tau3 / 6.0
				* (-NORMALIZED_SUN_GM_TT * u3_sun);                                                        // G

			// 月球摄动使用同一观测历元。
			double u3_moon;

			if (third_body_distances != nullptr)
			{
				u3_moon = third_body_distances->moon_inverse_radius_cubed;
			}
			else
			{
				double r_moon[3];
				Cart_MeanMoon(MJD, r_moon);
				const double r0_moon = sqrt(r_moon[0] * r_moon[0] + r_moon[1] * r_moon[1] + r_moon[2] * r_moon[2]);
				u3_moon = 1.0 / r0_moon / r0_moon / r0_moon;
			}


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

