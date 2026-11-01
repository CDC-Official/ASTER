// RKF78 的 NSG2、GCRS Cartesian 与 E-BCRS Cartesian 单步积分。
#include "integrator_rkf78.hpp"
#include "accel_complete.hpp"
#include <algorithm>
#include <limits>

// RKF78 是嵌入式单步法：在一个步长内依次计算 13 个阶段的状态导数，
// 再用不同权重组合出七阶、八阶近似，以两者之差估计本步的局部截断误差 TE。
// 自适应模式不满足容差就将步长减半重算；定步长模式仍计算 TE，但超限只警告。
// 两种模式使用相同积分公式；固定步长并不保证精度达标，仍需查看 TE 和步长收敛结果。
// 下方 NSG2／Cartesian 的流程相同，只是所积分的六个状态分量及其动力学导数不同。

void RKF78_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob,
	double t0_TT, double* X0, double h0_norm, double& h_return)
{
	double accepted_TE_norm;
	RKF78_NSG2(tracer, thisjob, t0_TT, X0, h0_norm, h_return, accepted_TE_norm);
}

void RKF78_NSG2(bastion::Overwatch& tracer,
	const OrbitPropagationJob& thisjob,
	double t0_TT,
	double* X0,
	double h0_norm, double& h_return, double& accepted_TE_norm)
{
	/*!
	 *
	* @param thisjob	当前预报任务的一系列设定，包含力模型、积分器模型等
	 * @param t0_TT		积分器当前这一单步的起始时刻MJDTT，千万注意不等于定义在thisjob里的开始全局积分的初始时刻！
	 * @param X0			积分器当前这一单步的起始6D状态量，在此基础上增加,为NSG2
	 * @param h0_norm	初始积分步长，归一化，其实可以从thisjob里获取，但是考虑到可能变步长还是单独传入一下
	 * @return			最终的这一步的步长，用于外部的时间累积
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);
	// 失败时不修改输入状态。自适应模式控制误差，定步长模式保留步长并监测误差。
	h_return = 0.0;
	accepted_TE_norm = std::numeric_limits<double>::quiet_NaN();
	if (tracer.Code() != ERR_NOERROR) return;
	if (!X0 || !std::isfinite(t0_TT) || !std::isfinite(h0_norm) || h0_norm == 0.0 ||
		!std::isfinite(thisjob.opn_setting.rkf.maxTE_norm) || thisjob.opn_setting.rkf.maxTE_norm <= 0.0 ||
		thisjob.opn_setting.rkf.iterCount_max < 1)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"RKF78 requires a finite epoch/state, a nonzero finite step, a positive finite tolerance and at least one trial.", bastion::ERROR);
		return;
	}
	for (int component = 0; component < 6; ++component)
		if (!std::isfinite(X0[component]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial state is not finite.", bastion::ERROR);
			return;
		}


	int count = 0;
	double h0 = h0_norm;

	int i, j, k;//for用
	int a;//调用各种函数用
	double initKep_array[6];
	thisjob.initKep_norm.ToArray(initKep_array);
	double temp_TE; //当前截断误差，当前时间
	double dt_mid[13] = {0};//13个中间点各自偏离t0时刻的时间长度（第一个为初始时刻为0）
	double f[13][6] = {0};//13个 dx/dt 的值
	double nsg2[13][6] = {0};//13个NSG2



	// 1. 在当前状态求第一个阶段导数。它不依赖试探步长，缩步重算时可以复用。
	nsg2[0][0] = X0[0];
	nsg2[0][1] = X0[1];
	nsg2[0][2] = X0[2];
	nsg2[0][3] = X0[3];
	nsg2[0][4] = X0[4];
	nsg2[0][5] = X0[5];
	Acceleration_NSG2(tracer, thisjob, t0_TT, 6, nsg2[0], f[0]);
	if (tracer.Code() != ERR_NOERROR) return;
	for (int component = 0; component < 6; ++component)
		if (!std::isfinite(f[0][component]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial derivative is not finite.", bastion::ERROR);
			return;
		}

	// 2. 依次构造其余阶段：用已知阶段导数预测中间状态，再在对应历元重新求导。
	// 这些是积分公式内部的试算节点，并不是要写入输出星历的观测／预报时刻。
	while (true)//用于变步长，第0步可以不写在循环内
	{
		count++;
		// 步长过小会使 MJD 无法前进；继续减半不能得到可靠的时变力模型。
		if (!std::isfinite(t0_TT + h0 * NORMALIZED_TIME_TO_DAY) ||
			t0_TT + h0 * NORMALIZED_TIME_TO_DAY == t0_TT)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 step cannot advance the epoch at double precision.", bastion::ERROR);
			return;
		}

		for (i = 1; i < 13; i++)
		{
			dt_mid[i] = A_rkf78[i] * h0;
			//变步长时，应当先把中间点的状态量清0！
			std::fill(std::begin(nsg2[i]), std::end(nsg2[i]), 0);

			//计算中间点的状态量 x = nsg2[i]
			for (j = 0; j < i; j++)//计算中间时刻的状态量，1.75式
				for (k = 0; k < 6; k++)//6个维度
					nsg2[i][k] += B_rkf78[i][j] * f[j][k];//先求sigma beta f

			for (k = 0; k < 6; k++)//6个维度
				nsg2[i][k] = nsg2[0][k] + h0 * nsg2[i][k];//1.75

			//计算中间点的f=dx/dt
			for (k = 0; k < 6; ++k)
				if (!std::isfinite(nsg2[i][k]))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 intermediate state is not finite.", bastion::ERROR);
					return;
				}
			Acceleration_NSG2(tracer, thisjob, t0_TT + dt_mid[i] * NORMALIZED_TIME_TO_DAY, 6, nsg2[i], f[i]);
			if (tracer.Code() != ERR_NOERROR) return;
			for (k = 0; k < 6; ++k)
				if (!std::isfinite(f[i][k]))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 intermediate derivative is not finite.", bastion::ERROR);
					return;
				}
		}

		// 3. 用嵌入公式之差估计本步误差；这里汇总六分量的绝对值，反向积分同样取 |h|。
		// TE 用于判断这一单步的数值积分质量，不代表从起点累积至此的轨道位置误差。
		temp_TE = 0;//当前截断误差 绝对值求和
		for (k = 0; k < 6; k++)
		{
			temp_TE += fabs(f[0][k] + f[10][k] - f[11][k] - f[12][k]);
		}
		temp_TE *= 41.0 / 840.0 * std::fabs(h0);
		if (!std::isfinite(temp_TE))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 truncation error is not finite.", bastion::ERROR);
			return;
		}

		// 定步长模式参考导师 RKF78 的处理：误差超限时警告，仍按给定步长推进。
		// 容差继续用于精度监测，不再参与是否缩步的判断。
		if (thisjob.opn_setting.rkf.fixed_step || temp_TE <= thisjob.opn_setting.rkf.maxTE_norm)
		{
			if (thisjob.opn_setting.rkf.fixed_step && temp_TE > thisjob.opn_setting.rkf.maxTE_norm)
			{
				std::ostringstream warning;
				warning << "RKF78 fixed-step TE " << std::scientific << std::setprecision(2) << temp_TE << " exceeds TOL "
					<< thisjob.opn_setting.rkf.maxTE_norm << "; fixed step retained: "
					<< std::fixed << std::setprecision(2) << h0 * NORMALIZED_TIME_TO_SEC << " s. Check step-size accuracy.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, warning.str(), bastion::WARNING);
			}

			if (thisjob.opn_setting.display_progress >= 2)
				fprintf(stderr, "RKF78 accepted | TE %.2e / TOL %.2e | step %+.2f s\n", temp_TE, thisjob.opn_setting.rkf.maxTE_norm, h0 * NORMALIZED_TIME_TO_SEC);

			break;
		}

		// 此轮被拒绝。达到试算上限就报错，不能用旧阶段导数和减半后的步长更新状态。
		if (count >= thisjob.opn_setting.rkf.iterCount_max)
		{
			std::ostringstream message;
			message << "RKF78 failed after " << count << " trial(s): TE = " << std::scientific << std::setprecision(2) << temp_TE
				<< ", TOL = " << thisjob.opn_setting.rkf.maxTE_norm << ", step [s] = " << std::fixed << std::setprecision(2) << h0 * NORMALIZED_TIME_TO_SEC
				<< ". Initial state preserved.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, message.str(), bastion::ERROR);
			return;
		}

		h0 /= 2.0;

		if (thisjob.opn_setting.display_progress >= 2)
			fprintf(stderr, "RKF78 rejected | TE %.2e / TOL %.2e | next step %+.2f s\n", temp_TE, thisjob.opn_setting.rkf.maxTE_norm, h0 * NORMALIZED_TIME_TO_SEC);
	}

	// 4. 步长接受后，按原有 C_rkf78 八阶权重组合阶段导数，得到这一步的末端状态。
	double sum_cf[6] = {0};
	for (i = 0; i < 13; i++)
		for (k = 0; k < 6; k++)
			sum_cf[k] += C_rkf78[i] * f[i][k];

	// 先完整检查候选终点，再一次性写回，避免失败时只更新了部分分量。
	double accepted_state[6];
	for (k = 0; k < 6; ++k)
	{
		accepted_state[k] = X0[k] + h0 * sum_cf[k];
		if (!std::isfinite(accepted_state[k]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 accepted state is not finite.", bastion::ERROR);
			return;
		}
	}
	std::copy(std::begin(accepted_state), std::end(accepted_state), X0);

	//	输出测试
	//	double temp_ele[6];
	//	CartesianToKeplerian(X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], temp_ele);
	//	cout << temp_ele[0] * ae / 1000 << '\t';
	//	cout << temp_ele[1] << '\t';
	//	cout << temp_ele[2] * DR2D << '\t';
	//	cout << temp_ele[3] * DR2D << '\t';
	//	cout << temp_ele[4] * DR2D << '\t';
	//	cout << temp_ele[5] * DR2D << '\t';
	h_return = h0;
	accepted_TE_norm = temp_TE;
	tracer.CheckOut(__func__);
}

void RKF78_Cart(bastion::Overwatch& tracer,
	const OrbitPropagationJob& thisjob,
	double t0_TT,
	double* X0,
	double h0_norm, double& h_return)
{
	// 保留原有六参数调用方式；不需要误差输出的调用方无需修改。
	double accepted_TE_norm;
	RKF78_Cart(tracer, thisjob, t0_TT, X0, h0_norm, h_return, accepted_TE_norm);
}

void RKF78_Cart(bastion::Overwatch& tracer,
	const OrbitPropagationJob& thisjob,
	double t0_TT,
	double* X0,
	double h0_norm, double& h_return, double& accepted_TE_norm)
{
	/*!
	 *
	* @param thisjob	当前预报任务的一系列设定，包含力模型、积分器模型等
	 * @param t0_TT		积分器当前这一单步的起始时刻MJDTT，千万注意不等于定义在thisjob里的开始全局积分的初始时刻！
	 * @param X0			积分器当前这一单步的起始6D状态量，在此基础上增加,为cart
	 * @param h0_norm	初始积分步长，归一化，其实可以从thisjob里获取，但是考虑到可能变步长还是单独传入一下
	 * @return			最终的这一步的步长，用于外部的时间累积
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);
	// 失败时不修改输入状态；定步长可接受误差超限的有限结果，自适应必须通过误差检验。
	h_return = 0.0;
	accepted_TE_norm = std::numeric_limits<double>::quiet_NaN();
	if (tracer.Code() != ERR_NOERROR) return;
	if (!X0 || !std::isfinite(t0_TT) || !std::isfinite(h0_norm) || h0_norm == 0.0 ||
		!std::isfinite(thisjob.opn_setting.rkf.maxTE_norm) || thisjob.opn_setting.rkf.maxTE_norm <= 0.0 ||
		thisjob.opn_setting.rkf.iterCount_max < 1)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"RKF78 requires a finite epoch/state, a nonzero finite step, a positive finite tolerance and at least one trial.", bastion::ERROR);
		return;
	}
	for (int component = 0; component < 6; ++component)
		if (!std::isfinite(X0[component]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial state is not finite.", bastion::ERROR);
			return;
		}


	int count = 0;
	double h0 = h0_norm;

	int i, j, k;//for用
	int a;//调用各种函数用
	double initKep_array[6];
	thisjob.initKep_norm.ToArray(initKep_array);
	double temp_TE; //当前截断误差，当前时间
	double dt_mid[13] = {0};//13个中间点各自偏离t0时刻的时间长度（第一个为初始时刻为0）
	double f[13][6] = {0};//13个 dx/dt 的值
	double cart[13][6] = {0};//13个cart

	//存在递推关系，必须要依次计算k[0]~k[12]
	//先计算t_mid[0]时刻的摄动加速度
	cart[0][0] = X0[0];
	cart[0][1] = X0[1];
	cart[0][2] = X0[2];
	cart[0][3] = X0[3];
	cart[0][4] = X0[4];
	cart[0][5] = X0[5];
	Acceleration_Cart(tracer, thisjob, t0_TT, 6, cart[0], f[0]);
	if (tracer.Code() != ERR_NOERROR) return;
	for (int component = 0; component < 6; ++component)
		if (!std::isfinite(f[0][component]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial derivative is not finite.", bastion::ERROR);
			return;
		}

	//递推计算f[1]~f[12]
	while (true)//用于变步长，第0步可以不写在循环内
	{
		count++;
		// 步长过小会使 MJD 无法前进；继续减半不能得到可靠的时变力模型。
		if (!std::isfinite(t0_TT + h0 * NORMALIZED_TIME_TO_DAY) ||
			t0_TT + h0 * NORMALIZED_TIME_TO_DAY == t0_TT)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 step cannot advance the epoch at double precision.", bastion::ERROR);
			return;
		}

		for (i = 1; i < 13; i++)
		{
			dt_mid[i] = A_rkf78[i] * h0;
			//变步长时，应当先把中间点的状态量清0！
			std::fill(std::begin(cart[i]), std::end(cart[i]), 0);

			//计算中间点的状态量 x = cart[i]
			for (j = 0; j < i; j++)//计算中间时刻的状态量，1.75式
				for (k = 0; k < 6; k++)//6个维度
					cart[i][k] += B_rkf78[i][j] * f[j][k];//先求sigma beta f

			for (k = 0; k < 6; k++)//6个维度
				cart[i][k] = cart[0][k] + h0 * cart[i][k];//1.75

			//计算中间点的f=dx/dt
			for (k = 0; k < 6; ++k)
				if (!std::isfinite(cart[i][k]))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 intermediate state is not finite.", bastion::ERROR);
					return;
				}
			Acceleration_Cart(tracer, thisjob, t0_TT + dt_mid[i] * NORMALIZED_TIME_TO_DAY, 6, cart[i], f[i]);
			if (tracer.Code() != ERR_NOERROR) return;
			for (k = 0; k < 6; ++k)
				if (!std::isfinite(f[i][k]))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 intermediate derivative is not finite.", bastion::ERROR);
					return;
				}
		}

		// 嵌入公式之差的六分量绝对值之和；反向积分也使用 |h|，误差不能为负。
		temp_TE = 0;//当前截断误差 绝对值求和
		for (k = 0; k < 6; k++)
		{
			temp_TE += fabs(f[0][k] + f[10][k] - f[11][k] - f[12][k]);
		}
		temp_TE *= 41.0 / 840.0 * std::fabs(h0);
		if (!std::isfinite(temp_TE))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 truncation error is not finite.", bastion::ERROR);
			return;
		}

		// 与 NSG2 共用 GCRS 积分设置；定步长时仍计算误差，超限仅警告。
		if (thisjob.opn_setting.rkf.fixed_step || temp_TE <= thisjob.opn_setting.rkf.maxTE_norm)
		{
			if (thisjob.opn_setting.rkf.fixed_step && temp_TE > thisjob.opn_setting.rkf.maxTE_norm)
			{
				std::ostringstream warning;
				warning << "RKF78 fixed-step TE " << std::scientific << std::setprecision(2) << temp_TE << " exceeds TOL "
					<< thisjob.opn_setting.rkf.maxTE_norm << "; fixed step retained: "
					<< std::fixed << std::setprecision(2) << h0 * NORMALIZED_TIME_TO_SEC << " s. Check step-size accuracy.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, warning.str(), bastion::WARNING);
			}

			if (thisjob.opn_setting.display_progress >= 2)
				fprintf(stderr, "RKF78 accepted | TE %.2e / TOL %.2e | step %+.2f s\n", temp_TE, thisjob.opn_setting.rkf.maxTE_norm, h0 * NORMALIZED_TIME_TO_SEC);

			break;
		}

		// 此轮被拒绝。达到试算上限就报错，不能用旧阶段导数和减半后的步长更新状态。
		if (count >= thisjob.opn_setting.rkf.iterCount_max)
		{
			std::ostringstream message;
			message << "RKF78 failed after " << count << " trial(s): TE = " << std::scientific << std::setprecision(2) << temp_TE
				<< ", TOL = " << thisjob.opn_setting.rkf.maxTE_norm << ", step [s] = " << std::fixed << std::setprecision(2) << h0 * NORMALIZED_TIME_TO_SEC
				<< ". Initial state preserved.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, message.str(), bastion::ERROR);
			return;
		}
		h0 /= 2.0;

		if (thisjob.opn_setting.display_progress >= 2)
			fprintf(stderr, "RKF78 rejected | TE %.2e / TOL %.2e | next step %+.2f s\n", temp_TE, thisjob.opn_setting.rkf.maxTE_norm, h0 * NORMALIZED_TIME_TO_SEC);
	}

	//完成。1.74第二个式子
	double sum_cf[6] = {0};
	for (i = 0; i < 13; i++)
		for (k = 0; k < 6; k++)
			sum_cf[k] += C_rkf78[i] * f[i][k];

	// 先完整检查候选终点，再一次性写回，避免失败时只更新了部分分量。
	double accepted_state[6];
	for (k = 0; k < 6; ++k)
	{
		accepted_state[k] = X0[k] + h0 * sum_cf[k];
		if (!std::isfinite(accepted_state[k]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 accepted state is not finite.", bastion::ERROR);
			return;
		}
	}
	std::copy(std::begin(accepted_state), std::end(accepted_state), X0);

	h_return = h0;
	accepted_TE_norm = temp_TE;
	tracer.CheckOut(__func__);
}

// =============================================================================
// 用于 EBCRS / TDB 下的单步积分
// =============================================================================

void RKF78_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& thisjob,
	double t0_TDB, double* X0, double h0_norm, double& h_return)
{
	double accepted_TE_norm;
	RKF78_Cart_EBCRS(tracer, thisjob, t0_TDB, X0, h0_norm, h_return, accepted_TE_norm);
}

void RKF78_Cart_EBCRS(bastion::Overwatch& tracer,
	const OrbitPropagationJob_EBCRS& thisjob,
	double t0_TDB,
	double* X0,
	double h0_norm, double& h_return, double& accepted_TE_norm)
{
	/*
	 * 函数用途：
	 *   在 EBCRS/TDB 框架下，对一个 6 维归一化 Cartesian 状态执行一次 RKF78 单步积分。
	 *   输入输出状态 X0 = [x, y, z, vx, vy, vz]，位置和速度均采用项目统一的地球归一化单位。
	 *   本函数只负责“当前这一单步”的推进；整段轨道预报中的循环、最后一步截断、星历记录等由外层函数负责。
	 *
	 * 时间含义：
	 *   t0_TDB 是当前 RKF78 单步的起始时刻，单位为 MJD(TDB)。
	 *   它不是整段轨道预报的全局初始时刻，外层传播函数每推进一步后都会更新下一步的 t0_TDB。
	 *   RKF78 的每个 stage 都有自己的中间时刻：
	 *       t_i_TDB = t0_TDB + dt_mid[i] * NORMALIZED_TIME_TO_DAY
	 *   其中 dt_mid[i] 是归一化时间长度，乘 NORMALIZED_TIME_TO_DAY 后转换为 MJD 的“天”。
	 *
	 * 状态递推：
	 *   cart[0] 保存当前单步起点状态，也就是输入的 X0。
	 *   cart[1] 到 cart[12] 是根据 RKF78 Butcher 表和已经算出的 f[j] 递推得到的中间状态。
	 *   每个 stage 调用 Acceleration_Cart_EBCRS 时，都会同时传入该 stage 的真实 MJD(TDB)
	 *   和该 stage 对应的中间 Cartesian 状态。因此加速度右函数看到的是“当前时刻 + 当前状态”，
	 *   不需要知道 RKF78 的 t0 或 dt 组织方式。
	 *
	 * 误差控制：
	 *   本函数沿用现有 RKF78 实现的局部截断误差估计：
	 *       TE = 41/840 * |h| * sum_k |f0[k] + f10[k] - f11[k] - f12[k]|
 *   自适应模式下，TE 超过 opn_setting.rkf.maxTE_norm 时减半步长并重算。
 *   定步长模式仍计算并显示 TE；超限时警告，但保持指定步长。
	 *   h_return 返回最终接受的归一化步长，外层传播函数应使用它累积时间。
	 *
	 * 参数说明：
	 *   tracer    运行日志和错误状态管理器。
	 *   thisjob   EBCRS 轨道预报任务设置，包括力模型、积分器参数、EOP 和卫星参数等。
	 *   t0_TDB    当前 RKF78 单步起点的 MJD(TDB)。
	 *   X0        输入/输出的 6 维归一化 EBCRS Cartesian 状态。
	 *   h0_norm   本单步的初始试探步长，归一化时间单位。
	 *   h_return  本单步最终接受的步长，归一化时间单位。
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);
	// 失败时不修改输入状态；只有通过误差检验的步长才能返回给外层累积时间。
	h_return = 0.0;
	accepted_TE_norm = std::numeric_limits<double>::quiet_NaN();
	if (tracer.Code() != ERR_NOERROR) return;
	if (!X0 || !std::isfinite(t0_TDB) || !std::isfinite(h0_norm) || h0_norm == 0.0 ||
		!std::isfinite(thisjob.opn_setting.rkf.maxTE_norm) || thisjob.opn_setting.rkf.maxTE_norm <= 0.0 ||
		thisjob.opn_setting.rkf.iterCount_max < 1)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"RKF78 requires a finite epoch/state, a nonzero finite step, a positive finite tolerance and at least one trial.", bastion::ERROR);
		return;
	}
	for (int component = 0; component < 6; ++component)
		if (!std::isfinite(X0[component]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial state is not finite.", bastion::ERROR);
			return;
		}


	int count = 0;
	double h0 = h0_norm;

	int i, j, k;//for用
	int a;//调用各种函数用

	double initKep_array[6];
	thisjob.initKep_norm.ToArray(initKep_array);

	double temp_TE; //当前截断误差，当前时间
	double dt_mid[13] = {0};//13个中间点各自偏离t0时刻的时间长度（第一个为初始时刻为0）
	double f[13][6] = {0};//13个 dx/dt 的值
	double cart[13][6] = {0};//13个cart



	//存在递推关系，必须要依次计算k[0]~k[12]
	//先计算t_mid[0]时刻的摄动加速度
	cart[0][0] = X0[0];
	cart[0][1] = X0[1];
	cart[0][2] = X0[2];
	cart[0][3] = X0[3];
	cart[0][4] = X0[4];
	cart[0][5] = X0[5];
	Acceleration_Cart_EBCRS(tracer, thisjob, t0_TDB, 6, cart[0], f[0]);
	if (tracer.Code() != ERR_NOERROR) return;
	for (int component = 0; component < 6; ++component)
		if (!std::isfinite(f[0][component]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial derivative is not finite.", bastion::ERROR);
			return;
		}

	//递推计算f[1]~f[12]
	while (true)//用于变步长，第0步可以不写在循环内
	{
		count++;
		// 步长过小会使 MJD 无法前进；继续减半不能得到可靠的时变力模型。
		if (!std::isfinite(t0_TDB + h0 * NORMALIZED_TIME_TO_DAY) ||
			t0_TDB + h0 * NORMALIZED_TIME_TO_DAY == t0_TDB)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 step cannot advance the epoch at double precision.", bastion::ERROR);
			return;
		}

		for (i = 1; i < 13; i++)
		{
			dt_mid[i] = A_rkf78[i] * h0;
			//变步长时，应当先把中间点的状态量清0！
			std::fill(std::begin(cart[i]), std::end(cart[i]), 0);

			//计算中间点的状态量 x = cart[i]
			for (j = 0; j < i; j++)//计算中间时刻的状态量，1.75式
				for (k = 0; k < 6; k++)//6个维度
					cart[i][k] += B_rkf78[i][j] * f[j][k];//先求sigma beta f

			for (k = 0; k < 6; k++)//6个维度
				cart[i][k] = cart[0][k] + h0 * cart[i][k];//1.75

			//计算中间点的f=dx/dt
			for (k = 0; k < 6; ++k)
				if (!std::isfinite(cart[i][k]))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 intermediate state is not finite.", bastion::ERROR);
					return;
				}
			Acceleration_Cart_EBCRS(tracer, thisjob, t0_TDB + dt_mid[i] * NORMALIZED_TIME_TO_DAY, 6, cart[i], f[i]);
			if (tracer.Code() != ERR_NOERROR) return;
			for (k = 0; k < 6; ++k)
				if (!std::isfinite(f[i][k]))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 intermediate derivative is not finite.", bastion::ERROR);
					return;
				}
		}

		// 嵌入公式之差的六分量绝对值之和；反向积分也使用 |h|，误差不能为负。
		temp_TE = 0;//当前截断误差 绝对值求和
		for (k = 0; k < 6; k++)
		{
			temp_TE += fabs(f[0][k] + f[10][k] - f[11][k] - f[12][k]);
		}
		temp_TE *= 41.0 / 840.0 * std::fabs(h0);
		if (!std::isfinite(temp_TE))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 truncation error is not finite.", bastion::ERROR);
			return;
		}

		// 固定步长只改变超限后的处理，阶段公式与截断误差估计均保持不变。
		if (temp_TE <= thisjob.opn_setting.rkf.maxTE_norm || thisjob.opn_setting.rkf.fixed_step)
		{
			if (temp_TE > thisjob.opn_setting.rkf.maxTE_norm)
			{
				std::ostringstream warning;
				warning << "RKF78 fixed-step TE " << std::scientific << std::setprecision(2) << temp_TE << " exceeds TOL "
					<< thisjob.opn_setting.rkf.maxTE_norm << "; fixed step retained: "
					<< std::fixed << std::setprecision(2) << h0 * NORMALIZED_TIME_TO_SEC << " s.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, warning.str(), bastion::WARNING);
			}

			if (thisjob.opn_setting.display_progress >= 2)
				fprintf(stderr, "RKF78 accepted | TE %.2e / TOL %.2e | step %+.2f s\n", temp_TE, thisjob.opn_setting.rkf.maxTE_norm, h0 * NORMALIZED_TIME_TO_SEC);

			break;
		}

		// 此轮被拒绝。达到试算上限就报错，不能用旧阶段导数和减半后的步长更新状态。
		if (count >= thisjob.opn_setting.rkf.iterCount_max)
		{
			std::ostringstream message;
			message << "RKF78 failed after " << count << " trial(s): TE = " << std::scientific << std::setprecision(2) << temp_TE
				<< ", TOL = " << thisjob.opn_setting.rkf.maxTE_norm << ", step [s] = " << std::fixed << std::setprecision(2) << h0 * NORMALIZED_TIME_TO_SEC
				<< ". Initial state preserved.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, message.str(), bastion::ERROR);
			return;
		}
		h0 /= 2.0;

		if (thisjob.opn_setting.display_progress >= 2)
			fprintf(stderr, "RKF78 rejected | TE %.2e / TOL %.2e | next step %+.2f s\n", temp_TE, thisjob.opn_setting.rkf.maxTE_norm, h0 * NORMALIZED_TIME_TO_SEC);
	}

	//完成。1.74第二个式子
	double sum_cf[6] = {0};
	for (i = 0; i < 13; i++)
		for (k = 0; k < 6; k++)
			sum_cf[k] += C_rkf78[i] * f[i][k];

	// 先完整检查候选终点，再一次性写回，避免失败时只更新了部分分量。
	double accepted_state[6];
	for (k = 0; k < 6; ++k)
	{
		accepted_state[k] = X0[k] + h0 * sum_cf[k];
		if (!std::isfinite(accepted_state[k]))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 accepted state is not finite.", bastion::ERROR);
			return;
		}
	}
	std::copy(std::begin(accepted_state), std::end(accepted_state), X0);

	h_return = h0;
	accepted_TE_norm = temp_TE;
	tracer.CheckOut(__func__);
}
