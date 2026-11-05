/**
 * @file core_iod_solver.cpp
 * @brief 使用单站或多站 RA/Dec 观测构造 GL 方程，并迭代求解参考历元的轨道。
 *
 * @details
 * 计算流程包括首轮状态求解、单次观测基本方程与约束的构造、全部观测的方程组装、
 * QR 最小二乘求解，以及更新后轨道的残差评价。InitialOrbitDetermination 组织整个迭代过程。
 * 简单迭代直接求下一轮状态，微分改正求六维状态改正量；附加约束用于微分改正。
 *
 * @par 数据约定
 * 观测时刻与参考历元使用 UTC MJD；RA、Dec 和测角残差使用 rad。
 * job.obs.tau 由调用方预先计算，保存各次观测相对参考历元的时间差，以 UNITTIME 归一化。
 * 状态和测站位置使用 GCRS，长度单位为 ae，时间单位为 UNITTIME，速度单位为 ae/UNITTIME。
 * 六维状态按 [x,y,z,vx,vy,vz] 排列。观测、视线和测站位置按同一下标逐条对应。
 *
 * @par 结果与错误
 * 最终轨道、残差、质量指标及每轮迭代记录写入 job.solution。
 * 函数通过 tracer 记录错误，调用方在返回后检查 tracer.Code()。
 */
#include "core_iod_solver.hpp"
#include "orbit_geometry.hpp"
#include "math_statistics.hpp"
#include <limits>
#include <sstream>

using namespace constants::Scale;
using constants::PI;

using namespace Eigen;
using namespace iod;
using namespace bastion;

/**
 * @brief 根据站心相对位置与观测视线平行的条件，构造三个 GL 基本方程。
 * @param[in,out] tracer 记录输入状态、传播系数或方程中的数值错误。
 * @param[in] iteration_method 迭代方法，取 IOD_ITERATION_SIMPLE 或 IOD_ITERATION_NEWTON。
 * @param[in] tau 观测时刻相对参考历元的时间差，除以 UNITTIME 后的无量纲值。
 * @param[in] los 本次观测的 GCRS 单位视线，无量纲。
 * @param[in] Rk_gcrs_norm 本次观测时刻的测站 GCRS 位置，已除以 ae。
 * @param[in] r0v0_norm 参考历元的六维 GCRS 状态，顺序为 [x,y,z,vx,vy,vz]。
 *     位置以 ae 归一化，速度以 ae/UNITTIME 归一化；微分改正在该状态附近线性化。
 * @param[in] fg 四个传播系数，顺序为 [f,g,fz,gz]，与归一化状态和时间配套使用。
 *     f、g 用于 x/y 分量，fz、gz 用于 z 分量。
 * @param[out] matCoeff 3×6 的基本方程系数矩阵，六列按位置、速度分量排列。
 *     简单迭代时乘下一轮完整状态，微分改正时乘六维状态改正量。
 * @param[out] vecRHS 三维方程右端，采用归一化长度尺度，与 matCoeff 的三行对应。
 *
 * @details
 * 消去未知站心距离后，三个基本行表达相对位置与视线平行，通常只有两个独立几何条件。
 * 简单迭代固定本轮传播系数求状态；微分改正使用 PartialFG 给出的低阶偏导近似，
 * 并从右端扣除当前状态的方程左端，求出所需的状态改正量。
 */
void iod_solver::Build_BasicEquations(Overwatch &tracer, int iteration_method, double tau,
	const Vector3d &los, const Vector3d &Rk_gcrs_norm,
	const VectorXd &r0v0_norm, const double fg[4],
	MatrixXd &matCoeff, VectorXd &vecRHS)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 检查参考状态、视线、测站位置与时间差是否为有效数值。
	if (r0v0_norm.size() != 6 || !r0v0_norm.allFinite() || !los.allFinite() ||
		!Rk_gcrs_norm.allFinite() || !std::isfinite(tau))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Basic GL equations require a finite six-component state, LOS, station "
			"position and tau.",
			ERROR);
		return;
	}

	// 四个传播系数都需要是有限数值，才能构造方程。
	for (int j = 0; j < 4; ++j)
	{
		if (!std::isfinite(fg[j]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite Lagrange coefficient at index " + std::to_string(j), ERROR);
			return;
		}
	}

	// 提取 x/y 与 z 分量各自使用的位置和速度传播系数。
	const double f = fg[0], g = fg[1], fz = fg[2], gz = fg[3];

	// GL 消去了每次观测的站心距离，要求传播位置减测站位置与单位视线平行。
	// 下面 D 的三行依次给出 ν*x-λ*z、ν*y-μ*z、μ*x-λ*y；
	// 其中 (λ,μ,ν)=los，因此基本方程可以写成 D*(r(t)-R(t))=0。
	MatrixXd matD(3, 3), matT(3, 6);
	matD << los(2), 0, -los(0),
			0, los(2), -los(1),
			los(1), -los(0), 0;

	// T 将参考状态映射到观测时刻的位置：r(t)=T*[r0;v0]。
	// x、y 分量使用 f/g，z 分量使用 fz/gz；这些系数均采用归一化时间。
	matT << f, 0, 0, g, 0, 0,
			0, f, 0, 0, g, 0,
			0, 0, fz, 0, 0, gz;

	// matPosition 等于 D*T，将参考状态映射到三个几何方程的左端。
	MatrixXd matPosition(3, 6);
	matPosition << f * los(2), 0, -fz * los(0), g * los(2), 0, -gz * los(0),
				   0, f * los(2), -fz * los(1), 0, g * los(2), -gz * los(1),
				   f * los(1), -f * los(0), 0, g * los(1), -g * los(0), 0;

	// 测站位置形成右端 D*R；三行使用与 matPosition 相同的分量顺序。
	vecRHS.resize(3);
	vecRHS << los(2) * Rk_gcrs_norm(0) - los(0) * Rk_gcrs_norm(2),
		los(2) * Rk_gcrs_norm(1) - los(1) * Rk_gcrs_norm(2),
		los(1) * Rk_gcrs_norm(0) - los(0) * Rk_gcrs_norm(1);

	if (iteration_method == IOD_ITERATION_SIMPLE)
	{
		// 简单迭代先把本轮 f/g 当作已知系数，直接求下一轮状态。
		// 此时右端只有 D*R，待求解的是 matCoeff*r0v0_next=vecRHS。
		matCoeff = matPosition;
	}
	else if (iteration_method == IOD_ITERATION_NEWTON)
	{
		if (!(r0v0_norm.head<3>().norm() > 0.0))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "PartialFG is undefined for zero geocentric radius.", ERROR);
			return;
		}
		// 微分改正需要知道方程随六个参考状态分量变化多少。
		// T 是直接乘状态的部分；S 是 f/g 自身随参考位置变化引起的部分，
		// 因而基本偏导矩阵为 D*(S+T)。PartialFG 提供关于参考位置的低阶偏导近似：
		// PFG[0..2]=df/dr0，PFG[3..5]=dg/dr0，均相对于归一化状态求导。
		double PFG[6];
		PartialFG(tau, r0v0_norm(0), r0v0_norm(1), r0v0_norm(2), PFG);
		MatrixXd matS1(3, 6), matS2(3, 6);
		for (int j = 0; j < 3; ++j)
		{
			matS1.row(j) << PFG[0] * r0v0_norm(j), PFG[1] * r0v0_norm(j), PFG[2] * r0v0_norm(j), 0,
				0, 0;
			matS2.row(j) << PFG[3] * r0v0_norm(j + 3), PFG[4] * r0v0_norm(j + 3),
				PFG[5] * r0v0_norm(j + 3), 0, 0, 0;
		}
		matCoeff = matD * (matS1 + matS2 + matT);

		// 当前状态已产生 D*T*r0v0，右端扣除这一部分，得到使基本方程
		// 趋于零所需的变化量。下面求解的未知数是 dr0v0，而非完整状态。
		for (int j = 0; j < 3; ++j) vecRHS(j) -= (matPosition.row(j) * r0v0_norm)(0);
	}
	else
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown GL iteration method.", ERROR);
		return;
	}

	// 确认组装后的系数与右端均为有限数值。
	if (!matCoeff.allFinite() || !vecRHS.allFinite())
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite GL basic coefficients or RHS.", ERROR);
		return;
	}

	tracer.CheckOut(__func__);
}

/**
 * @brief 在三个基本行之后，追加所选约束的线性化方程。
 * @param[in,out] tracer 记录轨道状态、站心距离或约束方程中的数值错误。
 * @param[in] constraint_type 约束组合，取 IOD_CONSTRAINT_*。
 *     无约束、方向、圆轨道、方向＋圆轨道、方向＋偏心率对应的总行数 m 分别为 3、4、5、6、5。
 * @param[in] constraint_ecc 指定偏心率，无量纲；仅在方向＋偏心率组合中使用。
 * @param[in] los 本次观测的 GCRS 单位视线，无量纲。
 * @param[in] Rk_gcrs_norm 本次观测时刻的测站 GCRS 位置，已除以 ae。
 * @param[in] r0v0_norm 当前参考历元的六维 GCRS 状态，顺序为 [x,y,z,vx,vy,vz]。
 *     位置以 ae 归一化，速度以 ae/UNITTIME 归一化；约束在此状态附近线性化。
 * @param[in] fg 本次观测的传播系数 [f,g,fz,gz]，与基本方程使用的系数一致。
 * @param[in,out] matCoeff 输入为 3×6 的基本方程系数矩阵；输出为 m×6 矩阵，前三行保持不变。
 * @param[in,out] vecRHS 输入为三个基本行的右端；输出为 m 维右端，与 matCoeff 同步追加约束行。
 *
 * @details
 * 方向条件要求目标位于观测视线的正向。圆轨道条件为速度模平方等于地心距离的倒数，
 * 以及位置与速度点积为零，这里采用归一化引力参数 mu=1。
 * 方向＋偏心率组合按三个基本行、一个偏心率行、一个方向行的顺序排列。
 * 各约束行乘六维状态改正量，右端为当前状态需要补偿的约束差额。
 * 约束按自身的归一化表达式构造，并随每条观测重复加入总方程组。
 */
void iod_solver::Append_ConstraintEquations(Overwatch &tracer, int constraint_type,
	double constraint_ecc, const Vector3d &los,
	const Vector3d &Rk_gcrs_norm,
	const VectorXd &r0v0_norm, const double fg[4],
	MatrixXd &matCoeff, VectorXd &vecRHS)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 按约束组合确定单次观测需要的总行数。
	int rows = 3;

	switch (constraint_type)
	{
	case IOD_CONSTRAINT_NONE:
		tracer.CheckOut(__func__);
		return;
	case IOD_CONSTRAINT_ORIENTATION:
		rows = 4;
		break;
	case IOD_CONSTRAINT_CIRCULAR:
		rows = 5;
		break;
	case IOD_CONSTRAINT_ORIENTATION_CIRCULAR:
		rows = 6;
		break;
	case IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY:
		// 在 GL-O 的方向修正上增加指定偏心率条件，共保留五行方程。
		rows = 5;
		break;
	default:
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown GL constraint type.", ERROR);
		return;
	}

	if (matCoeff.rows() != 3 || matCoeff.cols() != 6 || vecRHS.size() != 3 ||
		r0v0_norm.size() != 6 || !r0v0_norm.allFinite() || !los.allFinite() ||
		!Rk_gcrs_norm.allFinite() || !matCoeff.allFinite() || !vecRHS.allFinite())
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Constraints require three finite basic rows and a finite six-component "
			"reference state.",
			ERROR);
		return;
	}

	// conservativeResize 保留已经建立的三个基本行，再为附加条件准备位置。
	matCoeff.conservativeResize(rows, 6);
	vecRHS.conservativeResize(rows);

	const VectorXd &x = r0v0_norm; // x 是当前六维参考状态的只读引用，供下面的标量公式使用。
	const double f = fg[0], g = fg[1], fz = fg[2], gz = fg[3];
	const double r0 = x.head<3>().norm();
	const double v0 = x.tail<3>().norm();

	if (!(r0 > 0.0) || !std::isfinite(r0) || !std::isfinite(v0))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid geocentric radius or speed in GL constraints.", ERROR);
		return;
	}

	if (constraint_type == IOD_CONSTRAINT_CIRCULAR ||
		constraint_type == IOD_CONSTRAINT_ORIENTATION_CIRCULAR)
	{
		// 归一化引力参数 μ=1 时，圆轨道满足 v0²=1/r0 且 r0·v0=0。
		// 第一条限制速度大小，第二条限制径向速度；分别线性化后占两行。
		// 若同时使用方向约束，则先留出第 3 行给方向条件。
		const int row = constraint_type == IOD_CONSTRAINT_CIRCULAR ? 3 : 4;
		for (int j = 0; j < 3; ++j)
		{
			matCoeff(row, j) = x(j) / std::pow(r0, 3);
			matCoeff(row, j + 3) = 2 * x(j + 3);
			matCoeff(row + 1, j) = x(j + 3);
			matCoeff(row + 1, j + 3) = x(j);
		}
		vecRHS(row) = 1 / r0 - std::pow(v0, 2);
		vecRHS(row + 1) = -x(0) * x(3) - x(1) * x(4) - x(2) * x(5);
	}

	if (constraint_type == IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY)
	{
		if (!std::isfinite(constraint_ecc))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "The eccentricity constraint requires finite target eccentricity.", ERROR);
			return;
		}
		// 该组合在 GL-O 上增加一个偏心率条件，方向条件由下方的公共段构造。
		// 按零起始下标，偏心率条件位于第 3 行，方向条件位于第 4 行。
		// 两个约束行在每次观测的方程组中各出现一次。
		// G=(r0*v0*v0-1)^2+(2/r0-v0*v0)*r0r0dot^2-constraint_ecc^2。
		// 六个偏导与右端 -G 使用同一个函数；归一化 mu=1，v0 始终表示速度模。
		// 下方注释保留已修正的错误项，便于对照速度幂次求导中的区别。
		// 正确速度偏导已无 1/v0，因此零速度状态也能计算这一行。
		const double x0 = x(0), y0 = x(1), z0 = x(2), vx0 = x(3), vy0 = x(4), vz0 = x(5);
		const double r0r0dot = x0 * vx0 + y0 * vy0 + z0 * vz0;
		matCoeff(3, 0) =
			-2.0 * vx0 * v0 * v0 * r0r0dot
			- 2 * x0 * pow(r0r0dot, 2) / pow(r0, 3)
			+ 4 * vx0 * r0r0dot / r0
			// 原错误项来自 (r0*v0-1)^2；正确项须对 (r0*v0*v0-1)^2 求导。
			// + 2 * x0 * v0 * (-1.0 + r0 * v0) / r0;
			+ 2 * x0 * v0 * v0 * (-1.0 + r0 * v0 * v0) / r0;

		matCoeff(3, 1) =
			-2.0 * vy0 * v0 * v0 * r0r0dot
			- 2 * y0 * pow(r0r0dot, 2) / pow(r0, 3)
			+ 4 * vy0 * r0r0dot / r0
			// 原错误项来自 (r0*v0-1)^2；正确项须对 (r0*v0*v0-1)^2 求导。
			// + 2 * y0 * v0 * (-1.0 + r0 * v0) / r0;
			+ 2 * y0 * v0 * v0 * (-1.0 + r0 * v0 * v0) / r0;

		matCoeff(3, 2) =
			-2.0 * vz0 * v0 * v0 * r0r0dot
			- 2 * z0 * pow(r0r0dot, 2) / pow(r0, 3)
			+ 4 * vz0 * r0r0dot / r0
			// 原错误项来自 (r0*v0-1)^2；正确项须对 (r0*v0*v0-1)^2 求导。
			// + 2 * z0 * v0 * (-1.0 + r0 * v0) / r0;
			+ 2 * z0 * v0 * v0 * (-1.0 + r0 * v0 * v0) / r0;

		matCoeff(3, 3) =
			-2.0 * x0 * v0 * v0 * r0r0dot
			- 2 * vx0 * pow(r0r0dot, 2)
			+ 4 * x0 * r0r0dot / r0
			// 原错误项对速度模求导；正确项使用 d(v0*v0)/dvx0=2*vx0。
			// + 2 * vx0 * r0 * (-1.0 + r0 * v0) / v0;
			+ 4 * vx0 * r0 * (-1.0 + r0 * v0 * v0);

		matCoeff(3, 4) =
			-2.0 * y0 * v0 * v0 * r0r0dot
			- 2 * vy0 * pow(r0r0dot, 2)
			+ 4 * y0 * r0r0dot / r0
			// 原错误项对速度模求导；正确项使用 d(v0*v0)/dvy0=2*vy0。
			// + 2 * vy0 * r0 * (-1.0 + r0 * v0) / v0;
			+ 4 * vy0 * r0 * (-1.0 + r0 * v0 * v0);

		matCoeff(3, 5) =
			-2.0 * z0 * v0 * v0 * r0r0dot
			- 2 * vz0 * pow(r0r0dot, 2)
			+ 4 * z0 * r0r0dot / r0
			// 原错误项对速度模求导；正确项使用 d(v0*v0)/dvz0=2*vz0。
			// + 2 * vz0 * r0 * (-1.0 + r0 * v0) / v0;
			+ 4 * vz0 * r0 * (-1.0 + r0 * v0 * v0);

		// 偏心率条件在当前状态的差额，等于目标偏心率平方减当前表达式。
		vecRHS(3) = constraint_ecc * constraint_ecc + (v0 * v0 - 2.0 / r0) * pow(r0r0dot, 2) -
					pow(r0 * v0 * v0 - 1, 2);
	}

	// 这三种组合都包含方向条件。指定偏心率只限制轨道形状，不能替代
	// 观测视线的正向条件，因此方向＋偏心率组合也在这里构造方向行。
	if (constraint_type == IOD_CONSTRAINT_ORIENTATION ||
		constraint_type == IOD_CONSTRAINT_ORIENTATION_CIRCULAR ||
		constraint_type == IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY)
	{
		// 叉乘为零允许目标落在视线的正向或反向。方向条件要求
		// los·(r(t)-R(t))=|r(t)-R(t)|，用站心距离的正值选取视线正向。
		// 本方向行采用固定本轮 f/g 的线性化近似。
		const double dx = f * x(0) + g * x(3) - Rk_gcrs_norm(0);
		const double dy = f * x(1) + g * x(4) - Rk_gcrs_norm(1);
		const double dz = fz * x(2) + gz * x(5) - Rk_gcrs_norm(2);
		const double rho0 = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));
		if (!(rho0 > 0.0) || !std::isfinite(rho0))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__,
				"Direction constraint has zero or non-finite observer-to-target distance.", ERROR);
			return;
		}
		const int row = constraint_type == IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY ? 4 : 3;
		matCoeff.row(row) << f * (los(0) - dx / rho0), f * (los(1) - dy / rho0),
			fz * (los(2) - dz / rho0), g * (los(0) - dx / rho0), g * (los(1) - dy / rho0),
			gz * (los(2) - dz / rho0);
		vecRHS(row) =
			rho0 +
			(Rk_gcrs_norm(0) * los(0) + Rk_gcrs_norm(1) * los(1) + Rk_gcrs_norm(2) * los(2)) -
			(f * los(0) * x(0) + f * los(1) * x(1) + fz * los(2) * x(2) + g * los(0) * x(3) +
			 g * los(1) * x(4) + gz * los(2) * x(5));
	}

	// 确认组装后的系数与右端均为有限数值。
	if (!matCoeff.allFinite() || !vecRHS.allFinite())
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite GL constraint coefficients or RHS.", ERROR);
		return;
	}

	tracer.CheckOut(__func__);
}

/**
 * @brief 构造第 k 条观测的基本方程及所选约束方程。
 * @param[in,out] tracer 记录方程构造中的错误，并补充出错观测的编号和时刻。
 * @param[in] job 提供参考历元、观测 UTC MJD 时刻、逐条测站位置和算法设置。
 *     job.obs.tau(k) 为第 k 条观测相对参考历元的时间差，以 UNITTIME 归一化。
 *     测站位置为 GCRS 坐标，已除以 ae。
 * @param[in] k 观测下标，从 0 开始，取值为 0 到 N-1。
 * @param[in] los 第 k 条观测的 GCRS 单位视线，无量纲。
 * @param[in] r0v0_norm 参考历元的六维 GCRS 状态，顺序为 [x,y,z,vx,vy,vz]。
 *     位置以 ae 归一化，速度以 ae/UNITTIME 归一化。
 * @param[out] matCoeff m×6 的系数矩阵，前三行为基本方程，其后为约束方程。
 * @param[out] vecRHS m 维方程右端，与 matCoeff 逐行对应；各行使用相应方程的归一化尺度。
 * @param[in] third_body_distances 可选的 N 条只读日月距离记录，按观测顺序对应；
 *     须由相同 dynamicsModel、epochMJDUTC 和 obs.tau 准备，传入 nullptr 时现场计算日月位置。
 *
 * @details
 * 先根据当前参考状态计算该观测时刻的传播系数，再构造三个基本方程。
 * 微分改正按约束类型追加相应行；简单迭代只构造三个基本行。
 * 每次使用第 k 条观测对应的测站位置，因此同一目标的不同测站观测可以共同组装求解。
 */
void iod_solver::Build_ObservationEquations(Overwatch &tracer, const IodJob &job, int k,
	const Vector3d &los, const VectorXd &r0v0_norm,
	MatrixXd &matCoeff, VectorXd &vecRHS,
	const std::vector<iod::IodThirdBodyDistances>* third_body_distances)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 复用数组必须与本次观测逐条对应；后续按观测下标直接读取。
	if (third_body_distances != nullptr && third_body_distances->size() != static_cast<std::size_t>(job.obs.number))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Third-body distance count does not match IOD observations.", ERROR);
		return;
	}

	// 每列给出一条观测对应的测站位置，单位为 ae。
	const MatrixXd &R_gcrs_norm = job.observer.positions_gcrs_norm_emat;

	if (k < 0 || k >= job.obs.number || k >= job.obs.time_MJDUTC.size() || k >= job.obs.tau.size() ||
		R_gcrs_norm.rows() != 3 || k >= R_gcrs_norm.cols() || r0v0_norm.size() != 6 ||
		!r0v0_norm.allFinite() || !(r0v0_norm.head<3>().norm() > 0.0))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Cannot prepare GL observation at zero-based index " + std::to_string(k) +
			": check time/station dimensions and nonzero reference radius.",
			ERROR);
		return;
	}

	// 同一条记录决定时间、视线及观测者；测站编号变化时，计算仍使用对应列。
	const Vector3d Rk_gcrs_norm = R_gcrs_norm.col(k);
	const double tau = job.obs.tau(k);

	if (!std::isfinite(tau))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite normalized time for observation " + std::to_string(k + 1), ERROR);
		return;
	}

	// 将六维状态复制到传播系数函数使用的数组中。
	double state[6];
	for (int j = 0; j < 6; ++j) state[j] = r0v0_norm(j);

	double fg[4]; // [f,g,fz,gz]：归一化状态与归一化时间下的拉格朗日系数。
	Lagrange_Coefficients(job.dynamicsModel, job.epochMJDUTC, tau, state, fg,
		third_body_distances == nullptr ? nullptr : &(*third_body_distances)[k]);

	// 先建立三条共同的几何条件，再按约束种类追加相应行。
	Build_BasicEquations(tracer, job.iterMethod, tau, los, Rk_gcrs_norm, r0v0_norm, fg, matCoeff, vecRHS);

	if (tracer.Code() != ERR_NOERROR)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Basic-equation failure in observation " + std::to_string(k + 1) +
			", MJD UTC=" + std::to_string(job.obs.time_MJDUTC(k)),
			ERROR);
		return;
	}

	// 简单迭代的入口只接受无约束配置；它直接求状态，不进入改正量约束分支。
	if (job.iterMethod == IOD_ITERATION_NEWTON)
	{
		const double constraint_ecc =
			job.constraintType == IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY ? job.constraintEcc : 0.0;
		Append_ConstraintEquations(tracer, job.constraintType, constraint_ecc, los, Rk_gcrs_norm, r0v0_norm, fg, matCoeff, vecRHS);
		if (tracer.Code() != ERR_NOERROR)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__,
				"Constraint failure in observation " + std::to_string(k + 1) +
				", MJD UTC=" + std::to_string(job.obs.time_MJDUTC(k)),
				ERROR);
			return;
		}
	}

	tracer.CheckOut(__func__);
}

/**
 * @brief 将全部观测的方程按顺序拼接，形成待求解的线性系统。
 * @param[in,out] tracer 记录方程构造中的错误及出错函数的调用链。
 * @param[in] job 提供观测条数 N、UTC MJD 时刻、逐条 GCRS 测站位置，以及动力学、约束和迭代设置。
 *     job.obs.tau 为 N 维归一化相对时间数组，各元素对应观测时刻减去参考历元，再除以 UNITTIME。
 *     测站位置已除以 ae。
 * @param[in] L 3×N 的 GCRS 单位视线矩阵，无量纲；各列与 job 中的观测顺序一致。
 * @param[in] r0v0_norm 参考历元的六维 GCRS 状态，顺序为 [x,y,z,vx,vy,vz]。
 *     位置以 ae 归一化，速度以 ae/UNITTIME 归一化。
 * @param[out] matCoeff mN×6 的系数矩阵，每条观测占连续 m 行，六列依次对应位置和速度分量。
 * @param[out] vecRHS mN 维方程右端，行顺序与 matCoeff 一致。
 *     简单迭代求解 matCoeff*r0v0_next=vecRHS；微分改正求解 matCoeff*dr0v0=vecRHS。
 * @param[in] third_body_distances 可选的 N 条只读日月距离记录，按观测顺序对应；
 *     须由相同 dynamicsModel、epochMJDUTC 和 obs.tau 准备，传入 nullptr 时现场计算日月位置。
 *
 * @details
 * 简单迭代每条观测使用三个基本行。微分改正按约束组合使用 3、4、5 或 6 行。
 * 各行由对应的归一化方程构造，轨道约束随每条观测重复加入。
 * 输出是供 QR 最小二乘求解的长矩阵，不是 6×6 的法矩阵。
 */
void iod_solver::Build_LinearSystem(Overwatch &tracer, const IodJob &job, const MatrixXd &L,
	const VectorXd &r0v0_norm, MatrixXd &matCoeff,
	VectorXd &vecRHS,
	const std::vector<iod::IodThirdBodyDistances>* third_body_distances)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 复用数组必须与本次观测逐条对应；后续按观测下标直接读取。
	if (third_body_distances != nullptr && third_body_distances->size() != static_cast<std::size_t>(job.obs.number))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Third-body distance count does not match IOD observations.", ERROR);
		return;
	}

	const int obs_count = job.obs.number;
	if (obs_count < 1 || L.rows() != 3 || L.cols() != obs_count || !L.allFinite())
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "GL requires one finite LOS column per observation.", ERROR);
		return;
	}

	int rows = 3;
	// 每个记录占用相同的行数；简单迭代只构造原始 3N 个基本方程。

	if (job.iterMethod == IOD_ITERATION_NEWTON)
	{
		switch (job.constraintType)
		{
		case IOD_CONSTRAINT_NONE:
			rows = 3;
			break;
		case IOD_CONSTRAINT_ORIENTATION:
			rows = 4;
			break;
		case IOD_CONSTRAINT_CIRCULAR:
			rows = 5;
			break;
		case IOD_CONSTRAINT_ORIENTATION_CIRCULAR:
			rows = 6;
			break;
		case IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY:
			rows = 5;
			break;
		default:
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown GL constraint type.", ERROR);
			return;
		}
	}

	// 为全部观测分配 mN 行；六列对应参考历元的位置和速度分量。
	matCoeff.resize(static_cast<Eigen::Index>(rows) * obs_count, 6);
	vecRHS.resize(static_cast<Eigen::Index>(rows) * obs_count);

	// 逐条构造小方程组，再按观测顺序填入总矩阵。
	MatrixXd matObservation;
	VectorXd vecObservationRHS;

	for (int k = 0; k < obs_count; ++k)
	{
		const Vector3d los = L.col(k);

		// 用该条观测的视线、时刻和测站位置构造完整的小方程组。
		Build_ObservationEquations(tracer, job, k, los, r0v0_norm, matObservation, vecObservationRHS, third_body_distances);
		if (tracer.Code() != ERR_NOERROR) return;

		// 第 k 条记录放在 rows*k 起始的连续行中，所有行共同求六个未知量。
		matCoeff.block(static_cast<Eigen::Index>(rows) * k, 0, rows, 6) = matObservation;
		vecRHS.segment(static_cast<Eigen::Index>(rows) * k, rows) = vecObservationRHS;
	}

	tracer.CheckOut(__func__);
}

/**
 * @brief 根据初值类型，取得参考历元的起始位置和速度。
 * @param[in,out] tracer 记录初值读取或线性求解中的错误。
 * @param[in] job 提供初值类型、参考历元、归一化相对时间和逐条测站位置。
 *     历元采用 UTC MJD，测站位置为 GCRS 坐标，已除以 ae。
 *     job.obs.tau 为 N 维数组，保存各次观测相对参考历元的时间差，以 UNITTIME 归一化。
 *     使用指定初值时，job.initValue.initCart 应保存参考历元的归一化状态。
 * @param[in] L 3×N 的 GCRS 单位视线矩阵，无量纲；第 k 列对应第 k 条观测，
 *     三个分量为 [cos(Dec)*cos(RA), cos(Dec)*sin(RA), sin(Dec)]。
 * @param[out] r0v0_norm 六维 GCRS 起始状态，顺序为 [x,y,z,vx,vy,vz]，归一化。
 * @param[out] iteration 起始状态对应的轮次：默认策略求解一次线性方程，取 1；直接读取指定初值时取 0。
 *
 * @details
 * 默认策略令 f=fz=1、g=gz=tau，以 r(t)≈r0+tau*v0 近似目标运动，
 * 拼接各次观测的三个基本方程，通过最小二乘求出起始状态。
 * 初值策略由 job.initValue.type 选择，可在此处增加其他初值方法。
 */
void iod_solver::Solve_FirstIteration(Overwatch &tracer, const IodJob &job, const MatrixXd &L,
	VectorXd &r0v0_norm, int &iteration)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 弧段的观测数目
	const int obs_count = job.obs.number;

	// 首轮迭代，resize为六维状态向量
	r0v0_norm.resize(6);

	// 根据初值类型，选择线性求解或直接读取指定状态。
	switch (job.initValue.type)
	{
	case IOD_INITIAL_VALUE_DEFAULT:
	{
		// 尚无参考轨道时，使用教材默认策略，先采用 f=fz=1、g=gz=tau，拟合近似直线运动。
		MatrixXd matCoeff(3 * obs_count, 6), matObservation;
		VectorXd vecRHS(3 * obs_count), vecObservationRHS;

		// 基本方程的简单迭代形式使用固定传播系数，这里提供六维占位状态。
		const VectorXd placeholder_state = VectorXd::Zero(6);

		// 逐条构造三个基本方程，并拼接为整段观测的方程组。
		for (int k = 0; k < obs_count; ++k)
		{
			const double tau = job.obs.tau(k);
			const double fg[4] = {1.0, tau, 1.0, tau};

			const Vector3d los = L.col(k);
			const Vector3d Rk_gcrs_norm = job.observer.positions_gcrs_norm_emat.col(k);

			// 由本次观测的视线与测站位置构造三个基本方程。
			Build_BasicEquations(tracer, IOD_ITERATION_SIMPLE, tau, los, Rk_gcrs_norm,
				placeholder_state, fg, matObservation, vecObservationRHS);
			if (tracer.Code() != ERR_NOERROR) return;

			// 将本次观测的三行填入总矩阵中对应的位置。
			matCoeff.block(3 * k, 0, 3, 6) = matObservation;
			vecRHS.segment(3 * k, 3) = vecObservationRHS;
		}

		// 使用列主元 QR 拟合整段观测，得到供后续迭代使用的起始状态。
		r0v0_norm = matCoeff.colPivHouseholderQr().solve(vecRHS);
		iteration = 1;
		break;
	}

	case IOD_INITIAL_VALUE_SPECIFIED:
	{
		// 指定初值已位于本次参考历元，直接读取归一化直角位置和速度。
		job.initValue.initCart.ToArray(r0v0_norm.data());
		iteration = 0;
		break;
	}

	default:
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown IOD initial value type.", ERROR);
		return;
	}

	// 检查起始状态是否有效，后续传播需要有限状态和非零地心距离。
	if (!r0v0_norm.allFinite() || !(r0v0_norm.head<3>().norm() > 0.0))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "First GL state has non-finite components or zero radius.", ERROR);
		return;
	}

	tracer.CheckOut(__func__);
}

/**
 * @brief 根据当前参考轨道重新计算观测量，得到 O-C 残差及 RMS。
 * @param[in,out] tracer 记录传播结果、站心距离或残差中的数值错误。
 * @param[in,out] job 输入观测 UTC MJD 时刻、RA/Dec(rad)、逐条 GCRS 测站位置及动力学设置。
 *     job.obs.tau 为 N 维归一化相对时间数组，用于将参考状态传播至各次观测时刻。
 *     测站位置已除以 ae。结果写入 job.solution.resRa、resDec、RMS_angle 和 RMS_sum，单位均为 rad。
 *     resRa 为赤经差的带符号主值乘 cos(Dec_obs)，resDec 为 Dec_obs-Dec_calc。
 * @param[in] r0v0_norm 参考历元的六维 GCRS 状态，顺序为 [x,y,z,vx,vy,vz]。
 *     位置以 ae 归一化，速度以 ae/UNITTIME 归一化。
 * @param[out] rms_ra_cosdec_rad N 个东西方向残差 resRa 的均方根，单位 rad。
 * @param[out] rms_dec_rad N 个南北方向残差 resDec 的均方根，单位 rad。
 * @param[in] third_body_distances 可选的 N 条只读日月距离记录，按观测顺序对应；
 *     须由相同 dynamicsModel、epochMJDUTC 和 obs.tau 准备，传入 nullptr 时现场计算日月位置。
 *
 * @details
 * 每次调用都根据输入状态重新计算传播系数与观测量，使状态、残差和 RMS 对应同一条轨道。
 * 两个方向分别对 N 个残差分量计算 RMS，合成 RMS 对全部 2N 个分量计算：
 * RMS_angle=sqrt((sum(resRa²)+sum(resDec²))/(2N))，RMS_sum 保存相同数值。
 */
void iod_solver::Compute_Residuals(Overwatch &tracer, IodJob &job, const VectorXd &r0v0_norm,
	double &rms_ra_cosdec_rad, double &rms_dec_rad,
	const std::vector<iod::IodThirdBodyDistances>* third_body_distances)
{
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 复用数组必须与本次观测逐条对应；后续按观测下标直接读取。
	if (third_body_distances != nullptr && third_body_distances->size() != static_cast<std::size_t>(job.obs.number))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Third-body distance count does not match IOD observations.", ERROR);
		return;
	}

	const int obs_count = job.obs.number;

	if (obs_count < 1 || r0v0_norm.size() != 6 || !r0v0_norm.allFinite() ||
		!(r0v0_norm.head<3>().norm() > 0.0) || job.obs.time_MJDUTC.size() != obs_count ||
		job.obs.tau.size() != obs_count || !job.obs.tau.allFinite() ||
		job.obs.ra_rad.size() != obs_count || job.obs.dec_rad.size() != obs_count ||
		job.observer.positions_gcrs_norm_emat.rows() != 3 ||
		job.observer.positions_gcrs_norm_emat.cols() != obs_count)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Residual calculation requires matching observations, stations and reference state.",
			ERROR);
		return;
	}

	// 每列给出一条观测对应的测站位置，单位为 ae。
	const MatrixXd &R_gcrs_norm = job.observer.positions_gcrs_norm_emat;
	job.solution.resRa.resize(obs_count);
	job.solution.resDec.resize(obs_count);

	// 将六维状态复制到传播系数函数使用的数组中。
	double state[6];

	for (int j = 0; j < 6; ++j) state[j] = r0v0_norm(j);

	// 分别累计两个方向及全部测角分量的残差平方和。
	double sum_ra = 0.0, sum_dec = 0.0, sum_angle = 0.0;

	for (int k = 0; k < obs_count; ++k)
	{
		// 用当前状态重新计算fg系数，使下面的计算角度对应本轮更新后的轨道。
		const double tau = job.obs.tau(k);
		double fg[4];
		Lagrange_Coefficients(job.dynamicsModel, job.epochMJDUTC, tau, state, fg,
			third_body_distances == nullptr ? nullptr : &(*third_body_distances)[k]);

		// 由当前状态和fg系数计算目标在观测时刻的位置
		const Vector3d r_gcrs_norm(fg[0] * state[0] + fg[1] * state[3],
								   fg[0] * state[1] + fg[1] * state[4],
								   fg[2] * state[2] + fg[3] * state[5]);

	   // 由目标与测站的相对位置计算观测方向，检查是否为零或非有限数值。
		const Vector3d relative_gcrs_norm = r_gcrs_norm - R_gcrs_norm.col(k);

		if (!relative_gcrs_norm.allFinite() || !(relative_gcrs_norm.norm() > 0.0))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__,
				"Cannot compute angular residual for observation " +
				std::to_string(k + 1) +
				": zero or non-finite observer-to-target position.",
				ERROR);
			return;
		}

		// 由目标与测站的相对位置反算本次观测的 RA 和 Dec。
		double equatorial[3]; // 输出为归一化rho、RA(rad)、Dec(rad)。
		RelativePosition_to_Equatorial(r_gcrs_norm(0), r_gcrs_norm(1), r_gcrs_norm(2),
			R_gcrs_norm(0, k), R_gcrs_norm(1, k), R_gcrs_norm(2, k),
			equatorial);

		// 赤经先作 O-C ，再乘cos(Dec_o)
		// 赤纬直接作 O-C；这两个量分别表示东西方向与南北方向的角偏差。
		const double res_ra = std::cos(job.obs.dec_rad(k)) *
							  std::remainder(job.obs.ra_rad(k) - equatorial[1], 2 * PI);
		const double res_dec = job.obs.dec_rad(k) - equatorial[2];

		if (!std::isfinite(res_ra) || !std::isfinite(res_dec))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite angular residual for observation " + std::to_string(k + 1), ERROR);
			return;
		}

		// 保存逐条残差，同时累计本轮 RMS 所需的平方和。
		job.solution.resRa(k) = res_ra;
		job.solution.resDec(k) = res_dec;

		sum_ra += res_ra * res_ra;
		sum_dec += res_dec * res_dec;
		sum_angle += res_ra * res_ra + res_dec * res_dec;
	}

	// 两个方向各有 N 个分量；合成 RMS 对全部 2N 个分量取均方根。
	rms_ra_cosdec_rad = std::sqrt(sum_ra / obs_count);
	rms_dec_rad = std::sqrt(sum_dec / obs_count);
	job.solution.RMS_angle = std::sqrt(sum_angle / (2.0 * obs_count));
	job.solution.RMS_sum = job.solution.RMS_angle;

	tracer.CheckOut(__func__);
}

/**
 * @brief 将一轮已保存的 IOD 数值记录格式化为报告文字。
 * @param[in] record 参考历元 GCRS 状态及改正量 [x,y,z,vx,vy,vz]，位置以 ae、速度以 ae/UNITTIME 归一化；
 *     读取轮次、最大归一化改正量和合成测角 RMS(rad)，起始记录的改正量为 NaN。
 * @param[in] state_tolerance 判断六个归一化改正分量最大绝对值是否收敛的正阈值。
 * @param[in] iteration_limit 最大迭代轮数，用于显示是否达到迭代上限。
 * @return 本轮改正量、Cartesian(km、km/s)、Keplerian(km、deg)、RMS(rad、arcsec) 和停止状态的文字。
 * @details 仅将已保存的状态转换为显示根数，不重新求解或计算残差；可在候选选优后按历史顺序调用。
 */
std::string iod_solver::Format_IterationRecord(const IodIterationRecord &record, double state_tolerance, int iteration_limit)
{
	// 数值历史包含生成报告所需的状态与改正量，按原显示顺序读取。
	const Eigen::Map<const Eigen::Matrix<double, 6, 1>> r0v0_norm(record.r0v0_norm.data());
	const Eigen::Map<const Eigen::Matrix<double, 6, 1>> dr0v0_norm(record.dr0v0_norm.data());
	const double max_dr0v0_norm = record.max_dr0v0_norm;

	// 每轮依次显示改正量、更新后的状态和拟合质量；显示单位与结果文件的单位选择分别处理。
	const int label_width = 40;
	const int value_width = 18;
	const std::string continuation(label_width + 2, ' ');
	std::ostringstream message;
	message << "\n------------------------------------- Iteration " << record.iteration << " -------------------------------------\n\n"
			<< std::left << std::setw(label_width) << "- Correction of Cartesian [norm]" << ": ";

	// 起始记录的六个改正分量均为 NaN，按 N/A 显示；每三个分量结束一行。
	for (int j = 0; j < 6; ++j)
	{
		if (j == 3) message << '\n' << continuation;

		if (j % 3 == 0) message << "[ ";

		if (std::isfinite(dr0v0_norm(j)))
		{
			message << std::right << std::scientific << std::setprecision(6) << std::setw(value_width) << dr0v0_norm(j);
		}
		else
		{
			message << std::right << std::setw(value_width) << "N/A";
		}

		message << (j % 3 == 2 ? " ]" : ", ");
	}

	message << "\n\n" << std::left << std::setw(label_width) << "- Cartesian [km, km/s]" << ": [ "
			<< std::right << std::fixed << std::setprecision(6)
			<< std::setw(value_width) << r0v0_norm(0) * UNIT_LENGTH_M / 1000.0 << ", "
			<< std::setw(value_width) << r0v0_norm(1) * UNIT_LENGTH_M / 1000.0 << ", "
			<< std::setw(value_width) << r0v0_norm(2) * UNIT_LENGTH_M / 1000.0 << " ]\n"
			<< continuation << "[ " << std::setw(value_width) << r0v0_norm(3) * UNIT_LENGTH_M / (1000.0 * UNIT_TIME_S) << ", "
			<< std::setw(value_width) << r0v0_norm(4) * UNIT_LENGTH_M / (1000.0 * UNIT_TIME_S) << ", "
			<< std::setw(value_width) << r0v0_norm(5) * UNIT_LENGTH_M / (1000.0 * UNIT_TIME_S) << " ]\n";

	// 当前根数转换适用于非退化椭圆轨道；其他中间状态用 N/A 表示，不改变求解器的状态。
	std::array<double, 6> display_kep;
	display_kep.fill(std::numeric_limits<double>::quiet_NaN());
	const Vector3d position = r0v0_norm.head<3>();
	const Vector3d velocity = r0v0_norm.tail<3>();
	const double inverse_a = 2.0 / position.norm() - velocity.squaredNorm();

	if (inverse_a > 0.0 && position.cross(velocity).squaredNorm() > 0.0)
	{
		CartesianToKeplerian(r0v0_norm(0), r0v0_norm(1), r0v0_norm(2), r0v0_norm(3), r0v0_norm(4), r0v0_norm(5), display_kep.data());
	}

	const bool has_keplerian = std::all_of(display_kep.begin(), display_kep.end(), [](double value) { return std::isfinite(value); });
	message << '\n' << std::left << std::setw(label_width) << "- Keplerian [km, deg]" << ": ";

	for (int j = 0; j < 6; ++j)
	{
		if (j == 3) message << '\n' << continuation;

		if (j % 3 == 0) message << "[ ";

		if (has_keplerian)
		{
			const double scale = j == 0 ? UNIT_LENGTH_M / 1000.0 : (j >= 2 ? constants::Conversion::RAD_TO_DEGREE : 1.0);
			message << std::right << std::setw(value_width) << display_kep[j] * scale;
		}
		else
		{
			message << std::right << std::setw(value_width) << "N/A";
		}

		message << (j % 3 == 2 ? " ]" : ", ");
	}

	message << "\n\n" << std::left << std::setw(label_width) << "- RMS [rad / arcsec]" << ": "
			<< std::scientific << std::setprecision(2) << record.rms_angle_rad << " / "
			<< std::fixed << std::setprecision(2) << record.rms_angle_rad / constants::Conversion::ARCSEC_TO_RAD
			<< '\n' << std::setw(label_width) << "- Max absolute correction [norm]" << ": ";

	if (std::isfinite(max_dr0v0_norm))
	{
		message << std::scientific << std::setprecision(6) << max_dr0v0_norm;
	}
	else
	{
		message << "N/A";
	}

	// 显示状态使用下方相同的停止条件，先判断收敛，再判断迭代上限。
	message << '\n' << std::setw(label_width) << "- Iteration Status" << ": ";

	if (std::isfinite(max_dr0v0_norm) && max_dr0v0_norm <= state_tolerance)
	{
		message << "CONVERGED";
	}
	else if (record.iteration >= iteration_limit)
	{
		message << (std::isfinite(max_dr0v0_norm) ? "NOT CONVERGED — ITERATION LIMIT" : "INITIALIZATION — ITERATION LIMIT");
	}
	else if (!std::isfinite(max_dr0v0_norm))
	{
		message << "INITIALIZATION";
	}
	else
	{
		message << "CONTINUE";
	}

	message << '\n';
	return message.str();
}

/**
 * @brief 使用 GL 方法迭代求解参考历元的轨道，并保存计算结果。
 * @param[in,out] tracer 记录计算事件和错误；调用方通过 Code() 检查是否发生错误。
 * @param[in,out] job 输入初值、动力学模型、约束类型、迭代设置和 N 条 RA/Dec 观测。
 *     观测时刻和参考历元采用 UTC MJD，测角采用 rad；每条观测对应测站位置矩阵的一列。
 *     job.obs.tau 应由调用方预先计算，保存 N 个相对参考历元的时间差，以 UNITTIME 归一化。
 *     测站位置和轨道状态均在 GCRS 中表示，位置以 ae 归一化，速度以 ae/UNITTIME 归一化。
 *     结果写入 job.solution，包括轨道、O-C 残差、RMS、求解状态、每轮迭代记录和诊断指标。
 * @param[in,out] iteration_output 可选的逐轮报告输出流；非空时写入与屏幕显示相同的迭代文字。
 *     传入 nullptr 时不额外写出报告；屏幕显示仍由 job.printIterInfo 控制。
 * @param[in] third_body_distances 可选的 N 条只读日月距离记录，按观测顺序对应；
 *     须由相同 dynamicsModel、epochMJDUTC 和 obs.tau 准备，传入 nullptr 时现场计算日月位置。
 *
 * @details
 * 先取得起始状态，再逐轮计算残差、保存迭代记录并检查停止条件。
 * saveIterationReport 为 true 时逐轮文字另存 solution.iteration_report；为 false 时保留数值历史，供选优后生成报告。
 * 需要继续迭代时，组装全部观测的方程，通过列主元 QR 求解并更新状态。
 * 简单迭代直接求下一轮状态；微分改正求状态改正量，并可加入选定的约束。
 * 残差和 RMS 的单位为 rad。solution.flag 为正数时表示收敛轮次，
 * 为 -1 时表示达到迭代上限，为 -2 时表示数值失败；具体错误由 tracer 记录。
 * job.obs.sigma_RaDec 为东西、南北独立高斯误差共同的标准差（rad）；J_likelihood 保存该模型的负对数似然。
 */
void iod_solver::InitialOrbitDetermination(Overwatch &tracer, IodJob &job, std::ostream* iteration_output,
	const std::vector<iod::IodThirdBodyDistances>* third_body_distances)
{
	// 旧残差检验停用后，每次任务都先标记为未计算，包括提前返回的任务。
	job.solution.resRa_Q = job.solution.resRa_p = -1.0;
	job.solution.resDec_Q = job.solution.resDec_p = -1.0;

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 复用数组必须与本次观测逐条对应；后续按观测下标直接读取。
	if (third_body_distances != nullptr && third_body_distances->size() != static_cast<std::size_t>(job.obs.number))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Third-body distance count does not match IOD observations.", ERROR);
		return;
	}

	// =================================== 01 初始化结果与检查输入 ===================================

	// 清空job上次计算的迭代记录，将结果先标记为数值失败，随后在完成求解时更新状态。
	job.solution.iterations.clear();
	job.solution.iteration_report.clear();
	job.solution.flag = -2;
	job.solution.RMS_angle = std::numeric_limits<double>::quiet_NaN();

	const int obs_count = job.obs.number;
	// 每列给出一条观测对应的测站位置，单位为 ae。
	const MatrixXd &R_gcrs_norm = job.observer.positions_gcrs_norm_emat;

	if (job.obs.type != OBS_TYPE_RADEC || obs_count < 3)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"GL requires at least three RA/Dec records; type=" +
			std::to_string(job.obs.type) + ", number=" + std::to_string(obs_count),
			ERROR);
		return;
	}

	// 观测时刻、归一化相对时间、测角和测站位置应逐条对应。
	if (job.obs.time_MJDUTC.size() != obs_count || job.obs.tau.size() != obs_count || job.obs.ra_rad.size() != obs_count ||
		job.obs.dec_rad.size() != obs_count || R_gcrs_norm.rows() != 3 || R_gcrs_norm.cols() != obs_count)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "GL needs N UTC times, N normalized time offsets, N RA/Dec pairs and a 3-by-N station-position matrix.", ERROR);
		return;
	}

	// 检查参考历元及各条观测的数值。
	if (!std::isfinite(job.epochMJDUTC))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Reference epoch must be finite UTC MJD.", ERROR);
		return;
	}

	for (int k = 0; k < obs_count; ++k)
	{
		if (!std::isfinite(job.obs.time_MJDUTC(k)) || !std::isfinite(job.obs.tau(k)) || !std::isfinite(job.obs.ra_rad(k)) ||
			!std::isfinite(job.obs.dec_rad(k)) || !R_gcrs_norm.col(k).allFinite())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__,
				"Non-finite UTC time, normalized time offset, RA/Dec or GCRS station position in observation " +
				std::to_string(k + 1),
				ERROR);
			return;
		}
	}

	// 迭代次数和容差取正值，动力学类型应属于已支持的模型。
	if (job.iterCountMax <= 0 || !std::isfinite(job.iterTOL) || job.iterTOL <= 0.0)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "GL iteration limit and state tolerance must be positive.", ERROR);
		return;
	}

	if (job.dynamicsModel < IOD_DYNAMICS_TWOBODY || job.dynamicsModel > IOD_DYNAMICS_J2_SOLAR_LUNAR)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown GL dynamics model: " + std::to_string(job.dynamicsModel), ERROR);
		return;
	}

	// 后面的似然评价包含标准差的倒数与对数，因此标准差应为有限正数。
	if (!std::isfinite(job.obs.sigma_RaDec) || job.obs.sigma_RaDec <= 0.0 ||
		!std::isfinite(2 * PI * job.obs.sigma_RaDec * job.obs.sigma_RaDec) ||
		job.obs.sigma_RaDec * job.obs.sigma_RaDec == 0.0)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Likelihood noise sigma must yield a finite positive variance in radian units.", ERROR);
		return;
	}

	// =================================== 02 构造各次观测的单位视线Ls ===================================

	// L 的每一列对应一次观测的 GCRS 单位视线，和 R_gcrs_norm 的列一一对应。
	MatrixXd L(3, obs_count);

	for (int k = 0; k < obs_count; ++k)
	{
		L.col(k) << std::cos(job.obs.dec_rad(k)) * std::cos(job.obs.ra_rad(k)),
			std::cos(job.obs.dec_rad(k)) * std::sin(job.obs.ra_rad(k)),
			std::sin(job.obs.dec_rad(k));
	}

	// =================================== 03 生成起始状态与准备迭代变量 ===================================

	// 取得参考历元的起始状态，作为后续轨道传播与方程线性化的展开点。
	VectorXd r0v0_norm;
	int iteration_count = 0;
	Solve_FirstIteration(tracer, job, L, r0v0_norm, iteration_count);
	if (tracer.Code() != ERR_NOERROR) return;

	// 起始记录尚无“本轮减上轮”的状态改正，使用 NaN 表示这一项未发生。
	VectorXd dr0v0_norm = VectorXd::Constant(6, std::numeric_limits<double>::quiet_NaN());
	double max_dr0v0_norm = std::numeric_limits<double>::quiet_NaN();
	MatrixXd matCoeff;
	VectorXd vecRHS;
	int flag = -1;

	// 每轮先评价当前轨道并检查停止条件，需要继续时再构造方程和更新状态。
	while (true)
	{
		// =================================== 04 计算当前轨道的残差并记录本轮结果 ===================================

		// r0v0_norm 保存当前参考历元的完整状态。
		// 先重新传播，取得同一条轨道的残差，再将状态和质量指标一起存入历史。
		double rms_ra_cosdec_rad, rms_dec_rad;
		Compute_Residuals(tracer, job, r0v0_norm, rms_ra_cosdec_rad, rms_dec_rad, third_body_distances);
		if (tracer.Code() != ERR_NOERROR) return;

		// 状态、改正量和残差指标均属于同一轮，作为一条完整记录保存。
		IodIterationRecord record;
		record.iteration = iteration_count;

		for (int j = 0; j < 6; ++j)
		{
			record.r0v0_norm[j] = r0v0_norm(j);
			record.dr0v0_norm[j] = dr0v0_norm(j);
		}

		record.max_dr0v0_norm = max_dr0v0_norm;
		record.rms_ra_cosdec_rad = rms_ra_cosdec_rad;
		record.rms_dec_rad = rms_dec_rad;
		record.rms_angle_rad = job.solution.RMS_angle;
		job.solution.iterations.push_back(record);

		// 需要屏幕、流或即时文字留存时才格式化；候选可在选优后从数值记录生成报告。
		if (job.saveIterationReport || job.printIterInfo || iteration_output != nullptr)
		{
			const std::string message = Format_IterationRecord(record, job.iterTOL, job.iterCountMax);

			if (job.printIterInfo) std::cerr << message;

			if (job.saveIterationReport) job.solution.iteration_report += message;

			if (iteration_output != nullptr) *iteration_output << message;
		}

		// =================================== 05 检查收敛条件与迭代上限 ===================================

		// 用六个归一化改正分量的最大绝对值判断状态是否收敛。
		// 起始记录的改正为 NaN，跳过收敛判断；达到次数上限时，结束迭代并保存当前解。
		if (std::isfinite(max_dr0v0_norm) && max_dr0v0_norm <= job.iterTOL)
		{
			flag = iteration_count;
			break;
		}

		if (iteration_count >= job.iterCountMax) break;

		// =================================== 06 构造方程并求解状态更新 ===================================

		// 如果未收敛且未达到迭代上限，则继续迭代。在内部会先reshape系数矩阵和RHS
		Build_LinearSystem(tracer, job, L, r0v0_norm, matCoeff, vecRHS, third_body_distances);
		if (tracer.Code() != ERR_NOERROR) return;

		// 使用列主元 QR 直接求解长矩阵的最小二乘问题。微分改正求 dr0v0；
		// 简单迭代固定本轮 f/g 求下一轮状态，再与旧状态作差以使用同一终止指标。
		if (job.iterMethod == IOD_ITERATION_NEWTON)
		{
			dr0v0_norm = matCoeff.colPivHouseholderQr().solve(vecRHS);
		}
		else
		{
			const VectorXd next_r0v0_norm = matCoeff.colPivHouseholderQr().solve(vecRHS);
			dr0v0_norm = next_r0v0_norm - r0v0_norm;
			r0v0_norm = next_r0v0_norm;
		}

		if (!dr0v0_norm.allFinite())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__,
				"QR returned a non-finite state correction in iteration " +
				std::to_string(iteration_count + 1),
				ERROR);
			return;
		}

		if (job.iterMethod == IOD_ITERATION_NEWTON) r0v0_norm += dr0v0_norm;

		max_dr0v0_norm = dr0v0_norm.cwiseAbs().maxCoeff();
		++iteration_count;
	}

	// =================================== 07 保存最终轨道与求解状态 ===================================

	// 将最终状态保存为参考历元的 Cartesian 状态与 Keplerian 根数，均使用归一化单位。
	job.solution.cart.Set(r0v0_norm(0), r0v0_norm(1), r0v0_norm(2), r0v0_norm(3), r0v0_norm(4), r0v0_norm(5), true);
	double keplerian[6];
	CartesianToKeplerian(r0v0_norm(0), r0v0_norm(1), r0v0_norm(2), r0v0_norm(3), r0v0_norm(4), r0v0_norm(5), keplerian);
	job.solution.kep.Set(keplerian[0], keplerian[1], keplerian[2], keplerian[3], keplerian[4], keplerian[5], true);
	job.solution.flag = flag;
	job.solution.RMS_length = 0.0; // 本入口只有测角资料，未计算测距残差。
	const std::string finish =
		flag > 0 ? "IOD iteration converged after " : "Iteration limit reached after ";
	const std::string finish_message = finish + std::to_string(iteration_count) + " iterations.\n";

	tracer.LogEvent(__FILE__, __func__, __LINE__, finish_message, MESSAGE, FORCE_ALERT_OFF);

	// =================================== 08 计算残差趋势与相关性 ===================================

	// 对两个残差序列分别进行一至三阶多项式拟合，评价残差随时间的趋势。
	// 以归一化相对时间 job.obs.tau 为自变量，COD 保存无量纲的决定系数。
	for (int degree = 1; degree <= 3; ++degree)
	{
		Polyfit_coeff(degree, job.obs.tau, job.solution.resRa, job.solution.COD[degree - 1]);
		Polyfit_coeff(degree, job.obs.tau, job.solution.resDec, job.solution.COD[degree + 2]);
	}
	job.solution.resRaCOD1d = job.solution.COD[0];
	job.solution.resRaCOD2d = job.solution.COD[1];
	job.solution.resDecCOD1d = job.solution.COD[3];
	job.solution.resDecCOD2d = job.solution.COD[4];

	// Ljung–Box 是采用最大似然选优前保留的残差检验，当前 GLML 不使用其 Q/p。
	// 调用代码留作追溯；四个旧字段在进入求解时已置为 -1。
	/*
	if (obs_count >= 10)
	{
		const int lag = std::min(obs_count / 5, 10);
		ljungBoxTest(job.solution.resRa, lag, job.solution.resRa_Q, job.solution.resRa_p);
		ljungBoxTest(job.solution.resDec, lag, job.solution.resDec_Q, job.solution.resDec_p);
	}
	*/

	// =================================== 09 计算似然评价与几何诊断 ===================================

	// 东西、南北方向的独立高斯误差共用标准差 sigma（rad）。
	// resRa 已乘观测赤纬的 cos，直接与 resDec 一起构成负对数似然，不再缩放标准差。
	const double sigma = job.obs.sigma_RaDec;
	job.solution.J_likelihood = obs_count * std::log(2 * PI * sigma * sigma) +
		(job.solution.resRa.square().sum() + job.solution.resDec.square().sum()) / (2 * sigma * sigma);

	// 用第一条观测所属测站在参考历元的位置，计算其与轨道平面的夹角（deg）。
	// 方程构造与残差始终使用逐条测站位置，angle_norm 仅用于结果诊断。
	Vector3d orbit_normal;
	Orbit_Normal(keplerian[2], keplerian[3], orbit_normal);
	job.solution.angle_norm =
		90.0 - Angle_Between_Vectors(job.observer.position_gcrs_epoch_norm, orbit_normal) * DR2D;

	tracer.CheckOut(__func__);
}
