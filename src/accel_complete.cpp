//
// Created by 陈酊存 on 2024/10/16.
//

#include "orbit_geometry.hpp"
#include "accel_drag.hpp"
#include "accel_complete.hpp"
#include "core_propagation.hpp"//hpp中向前声明OrbitPropagationJob结构体，此处需要引用所在头文件补充其定义

namespace
{
	constexpr int PN_BODY_COUNT = 11;
	constexpr int PN_STATE_SIZE = 6;
	constexpr int PN_EARTH_INDEX = 2;

	double Norm3(const double x, const double y, const double z)
	{
		return std::sqrt(x * x + y * y + z * z);
	}

	bool HasInvalidMajorBodyDistance(const double major_body_bcrs_state_norm[PN_BODY_COUNT][PN_STATE_SIZE])
	{
		// 低层检查函数：只判断输入几何是否退化，不写日志。
		// 日志由上层 PostNewtonian_EIH_EBCRS 统一交给 Overwatch 处理，风格和项目里的其它模块更一致。
		for (int i = 0; i < PN_BODY_COUNT; ++i)
		{
			for (int j = i + 1; j < PN_BODY_COUNT; ++j)
			{
				const double dx = major_body_bcrs_state_norm[i][0] - major_body_bcrs_state_norm[j][0];
				const double dy = major_body_bcrs_state_norm[i][1] - major_body_bcrs_state_norm[j][1];
				const double dz = major_body_bcrs_state_norm[i][2] - major_body_bcrs_state_norm[j][2];
				if (Norm3(dx, dy, dz) <= 0.0) return true;
			}
		}
		return false;
	}

	void ComputeMajorBodyNewtonianAcceleration(const double major_body_bcrs_state_norm[PN_BODY_COUNT][PN_STATE_SIZE],
		double body_acc_norm[PN_BODY_COUNT][3],
		double body_potential_norm[PN_BODY_COUNT])
	{
		/*
		 * 函数用途：
		 *   为后牛顿 EIH 公式准备两个中间量：
		 *       body_acc_norm[i]       = 第 i 个大天体在其它大天体牛顿引力下的 BCRS 加速度；
		 *       body_potential_norm[i] = 第 i 个大天体处由其它大天体产生的 sum(mu_j / r_ij)。
		 *
		 * 为什么需要这个函数：
		 *   师兄 Fortran 的 FPNEWTON1/FPNEWTON2 里有 DDX 和 SUM2I 两组数组。
		 *   DDX 对应论文公式里天体 i 的 \ddot{r_i}；
		 *   SUM2I 对应论文公式里对 k != i 的 mu_k / Delta_ik 求和。
		 *   这两个量只依赖太阳系大天体彼此之间的 BCRS 状态，和卫星状态无关，
		 *   因此单独拆出来可以避免在 F_PN^S 和 F_PN^E 中重复写同一段双重循环。
		 *
		 * 单位说明：
		 *   输入位置、速度和输出加速度全部采用项目的地球归一化单位。
		 *   由于 mu_i 使用 GM_i / GM_Earth，并且距离用 ae 归一化，
		 *   -mu_j * (r_i-r_j)/|r_i-r_j|^3 得到的就是归一化加速度。
		 */

		for (int i = 0; i < PN_BODY_COUNT; ++i)
		{
			body_acc_norm[i][0] = 0.0;
			body_acc_norm[i][1] = 0.0;
			body_acc_norm[i][2] = 0.0;
			body_potential_norm[i] = 0.0;

			for (int j = 0; j < PN_BODY_COUNT; ++j)
			{
				if (j == i) continue;

				// r_i - r_j：从天体 j 指向天体 i 的向量。
				const double dx = major_body_bcrs_state_norm[i][0] - major_body_bcrs_state_norm[j][0];
				const double dy = major_body_bcrs_state_norm[i][1] - major_body_bcrs_state_norm[j][1];
				const double dz = major_body_bcrs_state_norm[i][2] - major_body_bcrs_state_norm[j][2];
				const double dr = Norm3(dx, dy, dz);
				const double dr3 = dr * dr * dr;

				// 经典牛顿加速度：a_i += -mu_j * (r_i - r_j) / |r_i - r_j|^3。
				body_acc_norm[i][0] -= MAJOR_BODY_EARTH_GM_RATIO[j] * dx / dr3;
				body_acc_norm[i][1] -= MAJOR_BODY_EARTH_GM_RATIO[j] * dy / dr3;
				body_acc_norm[i][2] -= MAJOR_BODY_EARTH_GM_RATIO[j] * dz / dr3;

				// EIH 公式中的势函数求和项 sum(mu_j / r_ij)。
				body_potential_norm[i] += MAJOR_BODY_EARTH_GM_RATIO[j] / dr;
			}
		}
	}

	void ComputePostNewtonianForMasslessBody(const double major_body_bcrs_state_norm[PN_BODY_COUNT][PN_STATE_SIZE],
		const double body_acc_norm[PN_BODY_COUNT][3],
		const double body_potential_norm[PN_BODY_COUNT],
		const double target_bcrs_state_norm[PN_STATE_SIZE],
		const double normalized_light_speed,
		double fpn_bcrs_norm[3])
	{
		/*
		 * 函数用途：
		 *   计算论文中“对于卫星”的 F_PN^S，也就是无质量航天器在 BCRS 中的 EIH 一阶后牛顿修正。
		 *   它是师兄 Fortran 子程序 FPNEWTON1 的 C++ 版本。
		 *
		 * 关键变量对应关系：
		 *   target_bcrs_state_norm  -> 论文中的 r 和 dot(r)，即卫星 BCRS 位置速度。
		 *   major_body_bcrs_state   -> 论文中的 r_i 和 dot(r_i)，即第 i 个大天体 BCRS 位置速度。
		 *   body_acc_norm[i]        -> 论文中的 ddot(r_i)，由其它大天体的牛顿引力给出。
		 *   normalized_light_speed  -> c / VUNIT，是归一化光速。
		 *
		 * 关于论文中的希腊字母 nu：
		 *   用户已经确认论文公式里有 OCR 错误，某些速度 v 被识别成希腊字母 nu。
		 *   本实现完全按速度向量处理，对应代码中的 target_bcrs_state_norm[3..5]。
		 */

		const double cl2 = normalized_light_speed * normalized_light_speed;
		double relative_state[PN_BODY_COUNT][PN_STATE_SIZE] = {{0.0}};
		double relative_distance[PN_BODY_COUNT] = {0.0};

		// 先准备卫星相对每个大天体的相对位置、相对速度和距离 Delta_i。
		for (int i = 0; i < PN_BODY_COUNT; ++i)
		{
			for (int k = 0; k < PN_STATE_SIZE; ++k)
			{
				relative_state[i][k] = target_bcrs_state_norm[k] - major_body_bcrs_state_norm[i][k];
			}
			relative_distance[i] = Norm3(relative_state[i][0], relative_state[i][1], relative_state[i][2]);
		}

		// SUM1 对应 sum_j(mu_j / |r - r_j|)，即卫星所在位置的主天体势函数求和。
		double spacecraft_potential_norm = 0.0;
		for (int j = 0; j < PN_BODY_COUNT; ++j)
		{
			spacecraft_potential_norm += MAJOR_BODY_EARTH_GM_RATIO[j] / relative_distance[j];
		}

		// v^2，注意这是卫星 BCRS 速度平方，不是希腊字母 nu。
		const double v2 =
			target_bcrs_state_norm[3] * target_bcrs_state_norm[3] +
			target_bcrs_state_norm[4] * target_bcrs_state_norm[4] +
			target_bcrs_state_norm[5] * target_bcrs_state_norm[5];

		fpn_bcrs_norm[0] = 0.0;
		fpn_bcrs_norm[1] = 0.0;
		fpn_bcrs_norm[2] = 0.0;

		for (int i = 0; i < PN_BODY_COUNT; ++i)
		{
			const double dr = relative_distance[i];
			const double dr3 = dr * dr * dr;

			// v_i^2：第 i 个大天体 BCRS 速度平方。
			const double vi2 =
				major_body_bcrs_state_norm[i][3] * major_body_bcrs_state_norm[i][3] +
				major_body_bcrs_state_norm[i][4] * major_body_bcrs_state_norm[i][4] +
				major_body_bcrs_state_norm[i][5] * major_body_bcrs_state_norm[i][5];

			// v dot v_i：卫星速度与第 i 个大天体速度的点积。
			const double vvi =
				target_bcrs_state_norm[3] * major_body_bcrs_state_norm[i][3] +
				target_bcrs_state_norm[4] * major_body_bcrs_state_norm[i][4] +
				target_bcrs_state_norm[5] * major_body_bcrs_state_norm[i][5];

			// [(r-r_i) dot v_i] / Delta_i。
			const double rvi =
				(relative_state[i][0] * major_body_bcrs_state_norm[i][3] +
				 relative_state[i][1] * major_body_bcrs_state_norm[i][4] +
				 relative_state[i][2] * major_body_bcrs_state_norm[i][5]) / dr;

			// (r-r_i) dot ddot(r_i)。
			const double rai =
				relative_state[i][0] * body_acc_norm[i][0] +
				relative_state[i][1] * body_acc_norm[i][1] +
				relative_state[i][2] * body_acc_norm[i][2];

			// 第一组项：乘在相对位置 (r-r_i) 前面的系数。
			double temp1 =
				-4.0 * spacecraft_potential_norm -
				body_potential_norm[i] +
				v2 +
				2.0 * vi2 -
				4.0 * vvi -
				1.5 * rvi * rvi -
				0.5 * rai;
			temp1 = -MAJOR_BODY_EARTH_GM_RATIO[i] * temp1 / dr3 / cl2;

			// 第二组项：乘在相对速度 (v-v_i) 前面的系数。
			const double sum3 =
				(4.0 * target_bcrs_state_norm[3] - 3.0 * major_body_bcrs_state_norm[i][3]) * relative_state[i][0] +
				(4.0 * target_bcrs_state_norm[4] - 3.0 * major_body_bcrs_state_norm[i][4]) * relative_state[i][1] +
				(4.0 * target_bcrs_state_norm[5] - 3.0 * major_body_bcrs_state_norm[i][5]) * relative_state[i][2];
			const double temp2 = MAJOR_BODY_EARTH_GM_RATIO[i] * sum3 / dr3 / cl2;

			// 第三组项：乘在大天体牛顿加速度 ddot(r_i) 前面的系数。
			const double temp3 = 3.5 * MAJOR_BODY_EARTH_GM_RATIO[i] / dr / cl2;

			fpn_bcrs_norm[0] += temp1 * relative_state[i][0] + temp2 * relative_state[i][3] + temp3 * body_acc_norm[i][0];
			fpn_bcrs_norm[1] += temp1 * relative_state[i][1] + temp2 * relative_state[i][4] + temp3 * body_acc_norm[i][1];
			fpn_bcrs_norm[2] += temp1 * relative_state[i][2] + temp2 * relative_state[i][5] + temp3 * body_acc_norm[i][2];
		}
	}

	void ComputePostNewtonianForMajorBody(const double major_body_bcrs_state_norm[PN_BODY_COUNT][PN_STATE_SIZE],
		const double body_acc_norm[PN_BODY_COUNT][3],
		const double body_potential_norm[PN_BODY_COUNT],
		const int target_body_index,
		const double normalized_light_speed,
		double fpn_bcrs_norm[3])
	{
		/*
		 * 函数用途：
		 *   计算论文中“对于地球”的 F_PN^E，更准确地说，是计算任意一个太阳系大天体在 BCRS 中的
		 *   EIH 一阶后牛顿修正。当前调用时 target_body_index 固定为地球。
		 *   它是师兄 Fortran 子程序 FPNEWTON2 的 C++ 版本。
		 *
		 * 为什么 EBCRS 要减去这个量：
		 *   EBCRS 是把 BCRS 原点平移到地心后的坐标系。
		 *   卫星在 EBCRS 中的相对加速度 = 卫星 BCRS 加速度 - 地球 BCRS 加速度。
		 *   所以后牛顿项也必须做同样的差分：
		 *       F_PN_EBCRS = F_PN^S - F_PN^E
		 */

		const double cl2 = normalized_light_speed * normalized_light_speed;
		double relative_state[PN_BODY_COUNT][PN_STATE_SIZE] = {{0.0}};
		double relative_distance[PN_BODY_COUNT] = {0.0};

		// 准备目标天体相对其它每个大天体的位置、速度和距离。
		for (int i = 0; i < PN_BODY_COUNT; ++i)
		{
			if (i == target_body_index) continue;

			for (int k = 0; k < PN_STATE_SIZE; ++k)
			{
				relative_state[i][k] =
					major_body_bcrs_state_norm[target_body_index][k] -
					major_body_bcrs_state_norm[i][k];
			}
			relative_distance[i] = Norm3(relative_state[i][0], relative_state[i][1], relative_state[i][2]);
		}

		// 目标天体 BCRS 速度平方。
		const double v2 =
			major_body_bcrs_state_norm[target_body_index][3] * major_body_bcrs_state_norm[target_body_index][3] +
			major_body_bcrs_state_norm[target_body_index][4] * major_body_bcrs_state_norm[target_body_index][4] +
			major_body_bcrs_state_norm[target_body_index][5] * major_body_bcrs_state_norm[target_body_index][5];

		// SUM1 对应目标天体处其它天体的 sum(mu_j / Delta_j)。
		double target_potential_norm = 0.0;
		for (int j = 0; j < PN_BODY_COUNT; ++j)
		{
			if (j == target_body_index) continue;
			target_potential_norm += MAJOR_BODY_EARTH_GM_RATIO[j] / relative_distance[j];
		}

		fpn_bcrs_norm[0] = 0.0;
		fpn_bcrs_norm[1] = 0.0;
		fpn_bcrs_norm[2] = 0.0;

		for (int i = 0; i < PN_BODY_COUNT; ++i)
		{
			if (i == target_body_index) continue;

			const double dr = relative_distance[i];
			const double dr3 = dr * dr * dr;

			const double vi2 =
				major_body_bcrs_state_norm[i][3] * major_body_bcrs_state_norm[i][3] +
				major_body_bcrs_state_norm[i][4] * major_body_bcrs_state_norm[i][4] +
				major_body_bcrs_state_norm[i][5] * major_body_bcrs_state_norm[i][5];

			const double vvi =
				major_body_bcrs_state_norm[target_body_index][3] * major_body_bcrs_state_norm[i][3] +
				major_body_bcrs_state_norm[target_body_index][4] * major_body_bcrs_state_norm[i][4] +
				major_body_bcrs_state_norm[target_body_index][5] * major_body_bcrs_state_norm[i][5];

			const double rvi =
				(relative_state[i][0] * major_body_bcrs_state_norm[i][3] +
				 relative_state[i][1] * major_body_bcrs_state_norm[i][4] +
				 relative_state[i][2] * major_body_bcrs_state_norm[i][5]) / dr;

			const double rai =
				relative_state[i][0] * body_acc_norm[i][0] +
				relative_state[i][1] * body_acc_norm[i][1] +
				relative_state[i][2] * body_acc_norm[i][2];

			double temp1 =
				-4.0 * target_potential_norm -
				body_potential_norm[i] +
				v2 +
				2.0 * vi2 -
				4.0 * vvi -
				1.5 * rvi * rvi -
				0.5 * rai;
			temp1 = -MAJOR_BODY_EARTH_GM_RATIO[i] * temp1 / dr3 / cl2;

			const double sum3 =
				(4.0 * major_body_bcrs_state_norm[target_body_index][3] - 3.0 * major_body_bcrs_state_norm[i][3]) * relative_state[i][0] +
				(4.0 * major_body_bcrs_state_norm[target_body_index][4] - 3.0 * major_body_bcrs_state_norm[i][4]) * relative_state[i][1] +
				(4.0 * major_body_bcrs_state_norm[target_body_index][5] - 3.0 * major_body_bcrs_state_norm[i][5]) * relative_state[i][2];
			const double temp2 = MAJOR_BODY_EARTH_GM_RATIO[i] * sum3 / dr3 / cl2;
			const double temp3 = 3.5 * MAJOR_BODY_EARTH_GM_RATIO[i] / dr / cl2;

			fpn_bcrs_norm[0] += temp1 * relative_state[i][0] + temp2 * relative_state[i][3] + temp3 * body_acc_norm[i][0];
			fpn_bcrs_norm[1] += temp1 * relative_state[i][1] + temp2 * relative_state[i][4] + temp3 * body_acc_norm[i][1];
			fpn_bcrs_norm[2] += temp1 * relative_state[i][2] + temp2 * relative_state[i][5] + temp3 * body_acc_norm[i][2];
		}
	}
}

void PostNewtonian_EIH_EBCRS(bastion::Overwatch& tracer,
	const double major_body_bcrs_state_norm[11][6],
	const double spacecraft_ebcrs_state_norm[6],
	double a_post_newtonian_ebcrs_norm[3])
{
	/*
	 * 函数用途：
	 *   计算 EBCRS/TDB 轨道预报中使用的一阶后牛顿 EIH 修正。
	 *   本函数对应师兄论文公式 2-21：
	 *       F_PN = F_PN^S - F_PN^E
	 *   其中 F_PN^S 是卫星在 BCRS 中的后牛顿修正，F_PN^E 是地球在 BCRS 中的后牛顿修正。
	 *
	 * 输入要求：
	 *   major_body_bcrs_state_norm 必须是 11 个大天体相对太阳系质心的 BCRS 状态。
	 *   spacecraft_ebcrs_state_norm 是卫星相对地心的 EBCRS 状态。
	 *   两者都必须已经使用项目统一归一化单位：
	 *       长度单位 = ae；
	 *       时间单位 = UNITTIME；
	 *       速度单位 = ae / UNITTIME。
	 *
	 * 计算步骤：
	 *   1. 先由 11 个大天体的 BCRS 状态计算每个天体的牛顿加速度和势函数求和；
	 *   2. 把卫星 EBCRS 状态加上地球 BCRS 状态，得到卫星 BCRS 状态；
	 *   3. 调用 FPNEWTON1 对应的 C++ 函数计算 F_PN^S；
	 *   4. 调用 FPNEWTON2 对应的 C++ 函数计算 F_PN^E；
	 *   5. 输出 F_PN^S - F_PN^E，也就是 EBCRS 中的相对论修正加速度。
	 *
	 * 注意：
	 *   本函数只实现师兄普通轨道预报分支实际加入的 EIH N 体项。
	 *   Fortran 中计算过的 Lense-Thirring 项和地球 J2 后牛顿项没有在 YHC_EBCRS_1 中加入，
	 *   因此这里也暂不加入，避免和参考程序不一致。
	 *
	 * Overwatch 追踪：
	 *   这是后牛顿模块对外暴露的主入口，不是最底层的数学循环，因此按照项目风格使用 tracer。
	 *   更低层的 ComputePostNewtonianForMasslessBody / ComputePostNewtonianForMajorBody 不使用 tracer，
	 *   避免在每个积分 stage 的内部双重循环里制造过多日志噪声。
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double body_acc_norm[PN_BODY_COUNT][3] = {{0.0}};
	double body_potential_norm[PN_BODY_COUNT] = {0.0};
	double spacecraft_bcrs_state_norm[PN_STATE_SIZE] = {0.0};
	double fpn_spacecraft_bcrs_norm[3] = {0.0};
	double fpn_earth_bcrs_norm[3] = {0.0};

	// 归一化光速 c_norm = c / (ae / UNITTIME)。
	const double normalized_light_speed = SPEED_OF_LIGHT_MPS / NORMALIZED_VELOCITY_TO_MPS;

	// 正常的 JPL 星历不可能让两个大天体位置完全重合；如果出现这种情况，通常意味着输入数组没有正确填充。
	if (HasInvalidMajorBodyDistance(major_body_bcrs_state_norm))
	{
		std::string msg = "Invalid major-body BCRS states for EBCRS post-Newtonian correction: at least two bodies have zero relative distance.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	ComputeMajorBodyNewtonianAcceleration(major_body_bcrs_state_norm, body_acc_norm, body_potential_norm);

	// EBCRS 只把 BCRS 原点纯几何平移到地心，不旋转轴向，保持同一 TDB 时刻。
	// 因此这里对位置和速度分别加回地球的 BCRS 分量，恢复卫星 BCRS 六维状态。
	for (int k = 0; k < PN_STATE_SIZE; ++k)
	{
		spacecraft_bcrs_state_norm[k] =
			spacecraft_ebcrs_state_norm[k] +
			major_body_bcrs_state_norm[PN_EARTH_INDEX][k];
	}

	ComputePostNewtonianForMasslessBody(major_body_bcrs_state_norm, body_acc_norm, body_potential_norm,
		spacecraft_bcrs_state_norm, normalized_light_speed, fpn_spacecraft_bcrs_norm);

	ComputePostNewtonianForMajorBody(major_body_bcrs_state_norm, body_acc_norm, body_potential_norm,
		PN_EARTH_INDEX, normalized_light_speed, fpn_earth_bcrs_norm);

	// EBCRS 中的相对加速度要减去地球原点本身的 BCRS 后牛顿加速度。
	a_post_newtonian_ebcrs_norm[0] = fpn_spacecraft_bcrs_norm[0] - fpn_earth_bcrs_norm[0];
	a_post_newtonian_ebcrs_norm[1] = fpn_spacecraft_bcrs_norm[1] - fpn_earth_bcrs_norm[1];
	a_post_newtonian_ebcrs_norm[2] = fpn_spacecraft_bcrs_norm[2] - fpn_earth_bcrs_norm[2];

	tracer.CheckOut(__func__);
}

void Acceleration_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob, const double MJD_TT, int n, const double* x, double* f)
{
	/*
	 * 函数用途：
	 *   计算 GCRS/TT 框架下 NonSingular2 状态量的右函数 f = dX/dt。
	 *   输入 x 是当前评估时刻的 6 维归一化 NonSingular2 状态，输出 f 是同一套变量的归一化导数。
	 *   本函数本身不负责推进时间，也不关心调用者使用的是 RKF78、固定步长积分器还是单点动力学检查。
	 *
	 * 时间含义：
	 *   MJD_TT 是当前右函数评估点的时刻，单位为 MJD(TT)。
	 *   如果本函数由 RKF78 调用，那么这个时刻应当已经由积分器根据 Butcher 表计算好：
	 *       MJD_TT = 当前单步起点 + 当前 stage 时间偏移
	 *   因此，这里不再传入 t0_TT 和 dt_norm，也不在函数内部拼接中间时刻。
	 *
	 * 状态量 x 的含义：
	 *   x 不是全局初始状态，也不一定是当前积分步起点状态。
	 *   对 RKF78 来说，x 是当前 stage 递推出来的中间状态。
	 *   本函数先把 NonSingular2 转换为 Cartesian，用 Cartesian 状态计算各种摄动加速度，
	 *   再把总摄动加速度投影到 RTN/STW 坐标系，并通过高斯行星方程转换为 NonSingular2 导数。
	 *
	 * 动力学项：
	 *   当前实现可按 thisjob.opn_model 的开关加入地球非球形引力、日月第三体、月球非球形引力、
	 *   光压和大气阻力。地球姿态、地固系转换和大气相关模型通过 MJD_TT -> UTC/UT1 的时间链获取 EOP。
	 *
	 * 参数说明：
	 *   tracer   运行日志和错误状态管理器。
	 *   thisjob  GCRS 轨道预报任务设置，包括力模型、EOP、卫星参数、大气数据等。
	 *   MJD_TT   当前评估点的 MJD(TT)。
	 *   n        状态维数，目前应为 6，暂未实际使用。
	 *   x        当前评估点的 6 维归一化 NonSingular2 状态。
	 *   f        输出的 6 维归一化 NonSingular2 状态导数。
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 先获取当前评估时刻的 EOP，用于计算各种时间系统和地固系转换。
	aEOP EOP_t;
	TimesMJD MJDs_t;
	double MJD_UTC = TT2UTC(MJD_TT);//当前时刻的UTC，用于获取EOP

	// 获取当前评估时刻的EOP，插值。
	thisjob.opn_EOP.Get_EOP(tracer, MJD_UTC, EOP_t, 1);
	if (tracer.Code() != ERR_NOERROR) return;

	int a; //用于部分函数的返回状态
	// 用当前评估时刻的UTC和插值得到的dut1生成其它时间尺度。
	MJDs_t = Utc2Else(EOP_t.dut1, MJD_UTC);

	//x[]定义/转换得到的三个等价的状态量
	NonSingular2 nsg2;
	Keplerian kep;
	Cartesian cart;
	double nsg2_array[6] = {0}, kep_array[6] = {0}, cart_array[6] = {0};
	nsg2.Set(tracer, x[0], x[1], x[2], x[3], x[4], x[5], true);
	if (tracer.Code() != ERR_NOERROR) return;

	nsg2.ToArray(nsg2_array);
	nsg2.ToKeplerian(tracer, kep_array);
	if (tracer.Code() != ERR_NOERROR) return;
	nsg2.ToCartesian(tracer, cart_array);
	if (tracer.Code() != ERR_NOERROR) return;
	double r3 =
		pow(cart_array[0] * cart_array[0] + cart_array[1] * cart_array[1] + cart_array[2] * cart_array[2], 1.5);//r^3
	double R[3] = {0};//卫星在地固系下的直角坐标
	//日月的位置速度矢量
	double rv_sun[6] = {0}, rv_moon[6] = {0};

	// 各种摄动加速度，都是xyz方向的
	double axyz_nsGravity_earth[3] = {0};//非球形引力
	double axyz_nsGravity_moon[3] = {0};//月球非球形引力
	double axyz_sunGravity[3] = {0};//太阳引力
	double axyz_moonGravity[3] = {0};//月球引力
	double axyz_srp[3] = {0};//光压加速度
	double rho_atm = 0, axyz_atmDrag[3] = {0};//大气密度和大气阻力加速度
	Vector3d axyz_sum(0.0, 0.0, 0.0);//总的xyz摄动加速度
	Vector3d artn_sum(0.0, 0.0, 0.0);;//总的rtn(STW)摄动加速度
	double artn_array[3] = {0};

	//处理太阳历表，避免重复生成
	if (!(thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_OFF && !thisjob.opn_model.use_drag &&
		  !thisjob.opn_model.use_SRP))
	{
		//光压、大气、太阳第三体只要有一个开启
		if (thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_JPL)//使用太阳数值历表，顺便算一下太阳第三体摄动
		{
			a = Get_JPLEPH(MJDs_t.TT, 11, 3, rv_sun);
		}
		else//使用分析历表，顺便算一下太阳第三体摄动，要用到UT1
		{
			Cart_MeanSun(MJDs_t.UT1, rv_sun);
		}
	}

//计算非球形引力
	if (thisjob.opn_model.use_earth_NSG)//开启非球形引力
	{
		Matrix3d gcrs_to_itrs_emat;    //	eigen matrix
		GCRS_to_ITRS_Matrix(MJDs_t.UTC, EOP_t, gcrs_to_itrs_emat);
		//	eigen vector
		Vector3d position_gcrs_norm_evec(cart_array[0], cart_array[1], cart_array[2]);
		Vector3d position_itrs_norm_evec = gcrs_to_itrs_emat * position_gcrs_norm_evec;//地固系
		R[0] = position_itrs_norm_evec[0];
		R[1] = position_itrs_norm_evec[1];
		R[2] = position_itrs_norm_evec[2];

		// 计算地固系下加速度
		double a_nsGravity_terrestrial[3] = {0};
		Earth_Nonspherical_Gravitational_Acceleration(a_nsGravity_terrestrial, R, 1, thisjob.opn_model.earth_zonal_order, thisjob.opn_model.earth_tesseral_degree);

		// 转换回惯性系
		Vector3d a_nsGravity_terrestrial_evec(a_nsGravity_terrestrial[0], a_nsGravity_terrestrial[1], a_nsGravity_terrestrial[2]);
		Vector3d a_nsGravity_inertial_evec = gcrs_to_itrs_emat.transpose() * a_nsGravity_terrestrial_evec;
		axyz_nsGravity_earth[0] = a_nsGravity_inertial_evec[0];
		axyz_nsGravity_earth[1] = a_nsGravity_inertial_evec[1];
		axyz_nsGravity_earth[2] = a_nsGravity_inertial_evec[2];
	}

	// 计算第三体摄动，太阳
	if (thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_JPL ||
		thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_ANALYTICAL)//1 调用数值历表 2分析
	{
		Thirdbody_Gravity_Acceleration(axyz_sunGravity, cart_array, rv_sun, SUN_EARTH_MASS_RATIO);
	}

	// 计算第三体摄动，月亮
	if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_JPL)//1 调用数值历表
	{
		a = Get_JPLEPH(MJDs_t.TT, 10, 3, rv_moon);
		Thirdbody_Gravity_Acceleration(axyz_moonGravity, cart_array, rv_moon, MOON_EARTH_MASS_RATIO);
	}
	if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_ANALYTICAL)//2 调用分析历表
	{
		Cart_MeanMoon(MJDs_t.UT1, rv_moon);
		Thirdbody_Gravity_Acceleration(axyz_moonGravity, cart_array, rv_moon, MOON_EARTH_MASS_RATIO);
	}

	// 260210 新增计算月球非球形引力，必须要在调用数值历表的前提下才可以，因为月球天平动的欧拉角需要从JPL数值历表中读.必须要在外部确保加载了
	if (thisjob.opn_model.use_lunar_NSG)
	{
		if (thisjob.opn_model.lunar_gravity_config != LUNAR_GRAVITY_JPL)
		{
			// 如果没有使用JPL数值历表，就无法获取月球天平动的欧拉角，无法计算月球非球形引力
			std::string msg = "To use lunar non-spherical gravity, you must use JPL numerical ephemeris for the Moon to obtain its Euler angles!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return;
		}
		else
		{
			if (!moon_gravity_field.is_ready)
			{
				// 如果月球非球形引力模型没有加载，就无法计算月球非球形引力
				std::string msg = "Lunar non-spherical gravity model not loaded!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return;
			}

			// 第一步，计算卫星在月心惯性系下的坐标，简单的平移，单位仍然是地球归一化
			Vector3d r_moon_inertial_earth_norm(
				cart_array[0] - rv_moon[0], cart_array[1] - rv_moon[1], cart_array[2] - rv_moon[2]);

			// 第二步，计算卫星在月心固连系下的坐标，需要用到月球天平动的欧拉角，单位仍然是地球归一化
			double Omega_p = 0, i_s = 0, Lambda = 0;//月球天平动的三个欧拉角，单位为弧度
			// 按理来说应该传入TDB 偷个懒只传TT 差别不大
			int ccc = Get_Libration(MJDs_t.TT, Omega_p, i_s, Lambda);
			Matrix3d moon_inertial_to_fixed_emat = Rz(Lambda) * Rx(i_s) * Rz(Omega_p);//月心惯性系到月心固连系的旋转矩阵
			Vector3d r_moon_fixed_earth_norm = moon_inertial_to_fixed_emat * r_moon_inertial_earth_norm;

			// 第三步，计算月球非球形引力加速度，得到月心固连系下的加速度，单位是月球归一化
			//卫星在月心固连系下的坐标，单位是月球归一化
			double r_moon_fixed_moon_norm[3] = {r_moon_fixed_earth_norm[0] * ae / moon_gravity_field.R_ref,
												r_moon_fixed_earth_norm[1] * ae / moon_gravity_field.R_ref,
												r_moon_fixed_earth_norm[2] * ae / moon_gravity_field.R_ref};

			double a_nsg_moon_fixed_moon_norm[3] = {0};
			Lunar_Nonspherical_Gravitational_Acceleration(a_nsg_moon_fixed_moon_norm, r_moon_fixed_moon_norm, thisjob.opn_model.lunar_zonal_order, thisjob.opn_model.lunar_tesseral_degree);

			// 第四步，把月固系下的月球归一化加速度转换回地心惯性系下的加速度，单位是地球归一化（加速度不涉及平移，旋转即可）
			Vector3d a_nsg_moon_fixed_moon_norm_evec(a_nsg_moon_fixed_moon_norm[0], a_nsg_moon_fixed_moon_norm[1], a_nsg_moon_fixed_moon_norm[2]);
			Vector3d a_nsg_moon_inertial_earth_norm_evec =
				moon_inertial_to_fixed_emat.transpose() * a_nsg_moon_fixed_moon_norm_evec *
				(moon_gravity_field.GM / GE) *
				pow(ae / moon_gravity_field.R_ref, 2);//最后乘上单位转换和归一化的反函数，得到在GRCS的地球归一化的加速度

			axyz_nsGravity_moon[0] = a_nsg_moon_inertial_earth_norm_evec[0];
			axyz_nsGravity_moon[1] = a_nsg_moon_inertial_earth_norm_evec[1];
			axyz_nsGravity_moon[2] = a_nsg_moon_inertial_earth_norm_evec[2];
		}
	}

	// 计算光压摄动
	if (thisjob.opn_model.use_SRP)// 开启光压
	{
		a = Solar_Radition_Acceleration(axyz_srp, cart_array, rv_sun, thisjob.opn_model.sat.AMR_SRP, thisjob.opn_model.sat.eta_SRP);
	}

	// 计算大气阻力摄动
	if (thisjob.opn_model.use_drag)
	{
		rho_atm = Get_AtmosphereDensity(cart_array, thisjob.refDensity, MJDs_t.UT1, thisjob.atm_list);
		AtmosphereDrag_Acceleration(axyz_atmDrag, cart_array, thisjob.opn_model.sat.Cd, thisjob.opn_model.sat.AMR_drag, rho_atm);
	}

	axyz_sum[0] =
		axyz_nsGravity_earth[0] + (axyz_sunGravity[0] + axyz_moonGravity[0] + axyz_nsGravity_moon[0]) +
		axyz_atmDrag[0] + axyz_srp[0];    //x\ddot
	axyz_sum[1] =
		axyz_nsGravity_earth[1] + (axyz_sunGravity[1] + axyz_moonGravity[1] + axyz_nsGravity_moon[1]) +
		axyz_atmDrag[1] + axyz_srp[1];    //y\ddot
	axyz_sum[2] =
		axyz_nsGravity_earth[2] + (axyz_sunGravity[2] + axyz_moonGravity[2] + axyz_nsGravity_moon[2]) +
		axyz_atmDrag[2] + axyz_srp[2];    //z\ddot

	Matrix3d inertial_to_rtn_emat;
	RTN_Matrix_from_Cart(cart_array[0], cart_array[1], cart_array[2], cart_array[3], cart_array[4], cart_array[5], inertial_to_rtn_emat);//使用cart计算RTN矩阵，应当与kep得到的结果一致
	artn_sum = inertial_to_rtn_emat * axyz_sum;

	//给出RTN/STW的摄动加速度
	artn_array[0] = artn_sum[0];// R/S
	artn_array[1] = artn_sum[1];// T
	artn_array[2] = artn_sum[2];// N/W

	//下一步，计算行星运动方程
	GaussPlanetaryEquation_RTNtoNSG2(tracer, nsg2, artn_array, f);
	if (tracer.Code() != ERR_NOERROR) return;

	tracer.CheckOut(__func__);
}

void Acceleration_Cart(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob, const double MJD_TT,  int n, const double* x, double* f)
{
	/*
	 * 函数用途：
	 *   计算 GCRS/TT 框架下 Cartesian 状态量的右函数 f = dX/dt。
	 *   输入 x 是当前评估时刻的 6 维归一化 Cartesian 状态：
	 *       X = [x, y, z, vx, vy, vz]
	 *   输出 f 的前三项是速度，后三项是归一化加速度：
	 *       f = [vx, vy, vz, ax, ay, az]
	 *
	 * 时间含义：
	 *   MJD_TT 是当前右函数评估点的时刻，单位为 MJD(TT)。
	 *   如果本函数由 RKF78 调用，那么积分器应当在调用前已经算好当前 stage 的 MJD(TT)。
	 *   因此，本函数只接收一个明确的当前时刻，不再接收 t0_TT 和 dt_norm。
	 *
	 * 状态量 x 的含义：
	 *   x 是当前评估点的状态，不一定是全局初值，也不一定是当前积分步起点状态。
	 *   对 RKF78 来说，x 是该 stage 根据前面斜率递推出来的中间状态。
	 *   这使得动力学模型只依赖“当前时刻 + 当前状态”，不绑定某一种积分器的内部组织方式。
	 *
	 * 动力学项：
	 *   本函数默认中心天体为地球，先计算归一化二体引力，再根据 thisjob.opn_model 的开关叠加
	 *   地球非球形引力、日月第三体、月球非球形引力、光压和大气阻力。
	 *   地球非球形引力和大气模型需要 EOP，因此会通过 MJD_TT -> UTC/UT1 的时间链获取地固系信息。
	 *
	 * 参数说明：
	 *   tracer   运行日志和错误状态管理器。
	 *   thisjob  GCRS 轨道预报任务设置，包括力模型、EOP、卫星参数、大气数据等。
	 *   MJD_TT   当前评估点的 MJD(TT)。
	 *   n        状态维数，目前应为 6，暂未实际使用。
	 *   x        当前评估点的 6 维归一化 Cartesian 状态。
	 *   f        输出的 6 维归一化 Cartesian 状态导数。
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	// 先获取当前评估时刻的EOP，用于计算各种时间系统和地固系转换。
	aEOP EOP_t;
	TimesMJD MJDs_t;

	double t_UTC = TT2UTC(MJD_TT);//当前时刻的UTC，用于获取EOP
	// 获取当前评估时刻的EOP，插值。
	thisjob.opn_EOP.Get_EOP(tracer, t_UTC, EOP_t, 1);
	if (tracer.Code() != ERR_NOERROR) return;

	int a; //用于部分函数的返回状态

	// 用当前评估时刻的UTC和插值得到的dut1生成其它时间尺度。
	MJDs_t = Utc2Else(EOP_t.dut1, t_UTC);

	//x[]定义/转换得到的三个等价的状态量
	NonSingular2 nsg2;
	Keplerian kep;
	Cartesian cart;
	double kep_array[6] = {0}, cart_array[6] = {0};
	cart.Set(x[0], x[1], x[2], x[3], x[4], x[5], true);
	cart.ToArray(cart_array);
	cart.ToKeplerian(tracer, kep_array);
	if (tracer.Code() != ERR_NOERROR) return;

	double r3 =
		pow(cart_array[0] * cart_array[0] + cart_array[1] * cart_array[1] + cart_array[2] * cart_array[2], 1.5);//r^3
	double R[3] = {0};//卫星在地固系下的直角坐标
	//日月的位置速度矢量
	double rv_sun[6] = {0}, rv_moon[6] = {0};

	// 各种摄动加速度，都是在惯性系且xyz方向的
	double axyz_tb[3] = {0};//质点二体引力
	double axyz_nsGravity_earth[3] = {0};//地球非球形引力
	double axyz_nsGravity_moon[3] = {0};//月球非球形引力
	double axyz_sunGravity[3] = {0};//太阳引力
	double axyz_moonGravity[3] = {0};//月球引力
	double axyz_srp[3] = {0};//光压加速度
	double rho_atm = 0, axyz_atmDrag[3] = {0};//大气密度和大气阻力加速度

	// 先计算二体加速度
	axyz_tb[0] = -cart_array[0] / r3;    //	-x/r^3
	axyz_tb[1] = -cart_array[1] / r3;    //	-y/r^3
	axyz_tb[2] = -cart_array[2] / r3;    //	-z/r^3

	// 处理太阳历表，避免重复生成
	if (!(thisjob.opn_model.solar_gravity_config == 0 && !thisjob.opn_model.use_drag && !thisjob.opn_model.use_SRP))
	{    // 光压、大气、太阳第三体只要有一个开启
		if (thisjob.opn_model.solar_gravity_config == 1)//使用太阳数值历表，顺便算一下太阳第三体摄动
		{
			// 按理来说应该传入TDB 偷个懒只传TT 差别不大
			a = Get_JPLEPH(MJDs_t.TT, 11, 3, rv_sun);
		}
		else//使用分析历表，顺便算一下太阳第三体摄动，要用到UT1
		{
			Cart_MeanSun(MJDs_t.UT1, rv_sun);
		}
	}

	//计算地球非球形引力
	if (thisjob.opn_model.use_earth_NSG)//开启非球形引力
	{
		Matrix3d gcrs_to_itrs_emat;//eigen matrix
		GCRS_to_ITRS_Matrix(MJDs_t.UTC, EOP_t, gcrs_to_itrs_emat);
		//eigen vector
		Vector3d position_gcrs_norm_evec(cart_array[0], cart_array[1], cart_array[2]);
		Vector3d position_itrs_norm_evec = gcrs_to_itrs_emat * position_gcrs_norm_evec;//地固系
		R[0] = position_itrs_norm_evec[0];
		R[1] = position_itrs_norm_evec[1];
		R[2] = position_itrs_norm_evec[2];

		// 计算地固系下加速度
		double a_nsGravity_terrestrial[3] = {0};
		Earth_Nonspherical_Gravitational_Acceleration(a_nsGravity_terrestrial, R, 1, thisjob.opn_model.earth_zonal_order, thisjob.opn_model.earth_tesseral_degree);

		// 转换回惯性系
		Vector3d a_nsGravity_terrestrial_evec(a_nsGravity_terrestrial[0], a_nsGravity_terrestrial[1], a_nsGravity_terrestrial[2]);
		Vector3d a_nsGravity_inertial_evec = gcrs_to_itrs_emat.transpose() * a_nsGravity_terrestrial_evec;
		axyz_nsGravity_earth[0] = a_nsGravity_inertial_evec[0];
		axyz_nsGravity_earth[1] = a_nsGravity_inertial_evec[1];
		axyz_nsGravity_earth[2] = a_nsGravity_inertial_evec[2];

	}

	// 计算第三体摄动，太阳
	if (thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_JPL ||
		thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_ANALYTICAL)//1 调用数值历表 2分析
	{
		Thirdbody_Gravity_Acceleration(axyz_sunGravity, cart_array, rv_sun, SUN_EARTH_MASS_RATIO);
	}

	// 计算第三体摄动，月亮
	if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_JPL)//1 调用数值历表
	{
		// 按理来说应该传入TDB 偷个懒只传TT 差别不大
		a = Get_JPLEPH(MJDs_t.TT, 10, 3, rv_moon);
		Thirdbody_Gravity_Acceleration(axyz_moonGravity, cart_array, rv_moon, MOON_EARTH_MASS_RATIO);
	}
	if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_ANALYTICAL)//2 调用分析历表
	{
		Cart_MeanMoon(MJDs_t.UT1, rv_moon);
		Thirdbody_Gravity_Acceleration(axyz_moonGravity, cart_array, rv_moon, MOON_EARTH_MASS_RATIO);
	}

	// 260210 新增计算月球非球形引力，必须要在调用数值历表的前提下才可以，因为月球天平动的欧拉角需要从JPL数值历表中读.必须要在外部确保加载了
	if (thisjob.opn_model.use_lunar_NSG)
	{
		if (thisjob.opn_model.lunar_gravity_config != LUNAR_GRAVITY_JPL)
		{
			// 如果没有使用JPL数值历表，就无法获取月球天平动的欧拉角，无法计算月球非球形引力
			std::string msg = "To use lunar non-spherical gravity, you must use JPL numerical ephemeris for the Moon to obtain its Euler angles!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return;
		}
		else
		{
			if (!moon_gravity_field.is_ready)
			{
				// 如果月球非球形引力模型没有加载，就无法计算月球非球形引力
				std::string msg = "Lunar non-spherical gravity model not loaded!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return;
			}

			// 第一步，计算卫星在月心惯性系下的坐标，简单的平移，单位仍然是地球归一化
			Vector3d r_moon_inertial_earth_norm(
				cart_array[0] - rv_moon[0], cart_array[1] - rv_moon[1], cart_array[2] - rv_moon[2]);

			// 第二步，计算卫星在月心固连系下的坐标，需要用到月球天平动的欧拉角，单位仍然是地球归一化
			double Omega_p = 0, i_s = 0, Lambda = 0;//月球天平动的三个欧拉角，单位为弧度
			// 按理来说应该传入TDB 偷个懒只传TT 差别不大
			int ccc = Get_Libration(MJDs_t.TT, Omega_p, i_s, Lambda);
			Matrix3d moon_inertial_to_fixed_emat = Rz(Lambda) * Rx(i_s) * Rz(Omega_p);//月心惯性系到月心固连系的旋转矩阵
			Vector3d r_moon_fixed_earth_norm = moon_inertial_to_fixed_emat * r_moon_inertial_earth_norm;

			// 第三步，计算月球非球形引力加速度，得到月心固连系下的加速度，单位是月球归一化
			//卫星在月心固连系下的坐标，单位是月球归一化
			double r_moon_fixed_moon_norm[3] = {r_moon_fixed_earth_norm[0] * ae / moon_gravity_field.R_ref,
												r_moon_fixed_earth_norm[1] * ae / moon_gravity_field.R_ref,
												r_moon_fixed_earth_norm[2] * ae / moon_gravity_field.R_ref};

			double a_nsg_moon_fixed_moon_norm[3] = {0};
			Lunar_Nonspherical_Gravitational_Acceleration(a_nsg_moon_fixed_moon_norm, r_moon_fixed_moon_norm, thisjob.opn_model.lunar_zonal_order, thisjob.opn_model.lunar_tesseral_degree);

			// 第四步，把月固系下的月球归一化加速度转换回地心惯性系下的加速度，单位是地球归一化（加速度不涉及平移，旋转即可）
			Vector3d a_nsg_moon_fixed_moon_norm_evec(a_nsg_moon_fixed_moon_norm[0], a_nsg_moon_fixed_moon_norm[1], a_nsg_moon_fixed_moon_norm[2]);
			Vector3d a_nsg_moon_inertial_earth_norm_evec =
				moon_inertial_to_fixed_emat.transpose() * a_nsg_moon_fixed_moon_norm_evec *
				(moon_gravity_field.GM / GE) *
				pow(ae / moon_gravity_field.R_ref, 2);//最后乘上单位转换和归一化的反函数，得到在GRCS的地球归一化的加速度

			axyz_nsGravity_moon[0] = a_nsg_moon_inertial_earth_norm_evec[0];
			axyz_nsGravity_moon[1] = a_nsg_moon_inertial_earth_norm_evec[1];
			axyz_nsGravity_moon[2] = a_nsg_moon_inertial_earth_norm_evec[2];
		}
	}


	// 计算光压摄动
	if (thisjob.opn_model.use_SRP)//1 开启光压
	{
		a = Solar_Radition_Acceleration(axyz_srp, cart_array, rv_sun, thisjob.opn_model.sat.AMR_SRP, thisjob.opn_model.sat.eta_SRP);
	}

	// 计算大气阻力摄动
	if (thisjob.opn_model.use_drag)
	{
		rho_atm = Get_AtmosphereDensity(cart_array, thisjob.refDensity, MJDs_t.UT1, thisjob.atm_list);
		AtmosphereDrag_Acceleration(axyz_atmDrag, cart_array, thisjob.opn_model.sat.Cd, thisjob.opn_model.sat.AMR_drag, rho_atm);
	}

	// FIXME:按理来说应该加入其他参数，但是先不管了，写定轨时候再加
	f[0] = cart_array[3];    //	x\dot=vx
	f[1] = cart_array[4];    //	y\dot=vy
	f[2] = cart_array[5];    //	z\dot=vz
	f[3] =
		(axyz_tb[0] + axyz_nsGravity_earth[0]) + (axyz_sunGravity[0] + axyz_moonGravity[0] + axyz_nsGravity_moon[0]) +
		axyz_atmDrag[0] + axyz_srp[0];    //x\ddot=ax
	f[4] =
		(axyz_tb[1] + axyz_nsGravity_earth[1]) + (axyz_sunGravity[1] + axyz_moonGravity[1] + axyz_nsGravity_moon[1]) +
		axyz_atmDrag[1] + axyz_srp[1];    //y\ddot=ay
	f[5] =
		(axyz_tb[2] + axyz_nsGravity_earth[2]) + (axyz_sunGravity[2] + axyz_moonGravity[2] + axyz_nsGravity_moon[2]) +
		axyz_atmDrag[2] + axyz_srp[2];    //z\ddot=az

	tracer.CheckOut(__func__);
}

void Acceleration_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& thisjob, double MJD_TDB, int n, const double* x, double* f)
{
	/*
	 * 函数用途：
	 *   计算 EBCRS/TDB 框架下 Cartesian 状态量的右函数 f = dX/dt。
	 *   这里的状态量 X = [x, y, z, vx, vy, vz]，位置和速度均为项目统一的地球归一化单位。
	 *   输出 f 的前三项是速度，后三项是归一化加速度，供 RKF78 积分器直接使用。
	 *
	 * 时间含义：
	 *   MJD_TDB 是当前右函数评估点的时刻，单位为 MJD(TDB)。
	 *   本函数不关心这个时刻是 RKF78 的单步起点、某个中间 stage，还是其它积分器给出的采样点。
	 *   积分器负责在调用本函数前计算好当前时刻；动力学右函数只使用“当前时刻 + 当前状态”。
	 *
	 * 状态量 x 的含义：
	 *   x 不是全局初值，也不是当前单步起点状态。
	 *   在 RKF78 的 13 个 stage 中，积分器会根据前面已经算出的斜率递推得到每个中间时刻的状态，
	 *   然后把该中间状态作为 x 传入本函数。因此，x[0..2] 就是当前 stage 对应时刻的 EBCRS 位置。
	 *   用 x[0..2] 参与 TDB -> TT 的位置项转换，是当前积分框架下应当采用的做法。
	 *
	 * 坐标与外部数据：
	 *   本函数假定 x 位于 EBCRS，也就是 BCRS 平移到地心后的坐标系。
	 *   日、月、行星等天体位置速度需要由 JPL DE 星历提供，因此调用本函数前必须已经 Read_JPLEPH。
	 *   地球非球形引力、大气阻力等需要地固系和 EOP 的模块，会通过 TDB -> TT -> UTC 的时间链获取 EOP。
	 *
	 * 参数说明：
	 *   thisjob   EBCRS 轨道预报任务设置，包括力模型、EOP、卫星参数和输出设置。
	 *   MJD_TDB   当前右函数评估点的 MJD(TDB)。
	 *   n         状态维数，目前应为 6，暂未实际使用。
	 *   x         当前评估点的 6 维 EBCRS Cartesian 归一化状态。
	 *   f         输出的 6 维归一化状态导数。
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	//用于部分函数的返回状态
	int a;

	// 先获取t时刻的EOP，用于计算各种时间系统
	aEOP EOP_t;
	TimesMJD MJDs_t;

	// x[0..2] 是当前评估点的 EBCRS 位置；转为米后用于 TDB -> TT 的位置项。
	double t_TT = TDB2TT(MJD_TDB, x[0] * ae, x[1] * ae, x[2] * ae);    // 当前评估点的TT
	double t_UTC = TT2UTC(t_TT);    //当前评估时刻的UTC，用于获取EOP

	// 获取t时刻的EOP，插值
	thisjob.opn_EOP.Get_EOP(tracer, t_UTC, EOP_t, 1);
	if (tracer.Code() != ERR_NOERROR) return;

	// 以 t_UTC 和插值得到的 dut1 生成当前评估点的其它时间尺度；TDB 保留由调用者传入的 MJD_TDB。
	MJDs_t = Utc2Else(EOP_t.dut1, t_UTC);
	MJDs_t.TDB = MJD_TDB;

	//x[]转为类
	Cartesian cart;
	double cart_array[6] = {0};
	cart.Set(x[0], x[1], x[2], x[3], x[4], x[5], true);
	cart.ToArray(cart_array);
	if (tracer.Code() != ERR_NOERROR) return;

	// 定义计算摄动需要的变量
	double r3 = pow(cart_array[0] * cart_array[0] + cart_array[1] * cart_array[1] + cart_array[2] * cart_array[2], 1.5);//r^3
	double R[3] = {0};//卫星在地固系下的直角坐标
	//日月的位置速度矢量
	double rv_sun[6] = {0}, rv_moon[6] = {0};

	// 各种摄动加速度，都是在惯性系且xyz方向的
	double axyz_tb[3] = {0};//质点二体引力
	double axyz_nsGravity_earth[3] = {0};//地球非球形引力
	double axyz_nsGravity_moon[3] = {0};//月球非球形引力
	double axyz_sunGravity[3] = {0};//太阳引力
	double axyz_moonGravity[3] = {0};//月球引力
	double axyz_planetGravity[3] = {0};//水金火木土天海冥的行星质点引力
	double axyz_postNewtonian[3] = {0};//后牛顿相对论修正加速度
	double axyz_srp[3] = {0};//光压加速度
	double rho_atm = 0, axyz_atmDrag[3] = {0};//大气密度和大气阻力加速度
	double major_body_bcrs_state_norm[11][6] = {{0.0}};//11个太阳系大天体相对太阳系质心的BCRS状态，用于后牛顿项

	// 先计算地球二体加速度
	axyz_tb[0] = -cart_array[0] / r3;    //	-x/r^3
	axyz_tb[1] = -cart_array[1] / r3;    //	-y/r^3
	axyz_tb[2] = -cart_array[2] / r3;    //	-z/r^3


	// 地球非球形引力采用如下近似：将当前 E-BCRS 位置视为 GCRS 位置，再旋转到 ITRS 计算摄动加速度。
	// 此处省略 E-BCRS 与 GCRS 之间的相对论坐标差异；其影响仅进入非球形摄动项，是坐标小差异对小摄动的修正。
	// 计算所得加速度用同一旋转矩阵的转置转回积分轴向，输入、传播和输出状态始终按 E-BCRS 解释。
	if (thisjob.opn_model.use_earth_NSG)
	{
		// 当前 TDB 和 E-BCRS 位置给出地球姿态计算所需的 UTC，旋转矩阵本身仍采用 GCRS 到 ITRS 的定义。
		Matrix3d gcrs_to_itrs_emat;
		GCRS_to_ITRS_Matrix(MJDs_t.UTC, EOP_t, gcrs_to_itrs_emat);

		// 将当前 E-BCRS 位置近似旋转到地固系，位置分量保持 ae 归一化。
		Vector3d position_ebcrs_norm_evec(cart_array[0], cart_array[1], cart_array[2]);
		Vector3d position_itrs_approx_norm_evec = gcrs_to_itrs_emat * position_ebcrs_norm_evec;
		R[0] = position_itrs_approx_norm_evec[0];
		R[1] = position_itrs_approx_norm_evec[1];
		R[2] = position_itrs_approx_norm_evec[2];

		// 计算地固系下加速度
		double a_nsGravity_terrestrial[3] = {0};
		Earth_Nonspherical_Gravitational_Acceleration(a_nsGravity_terrestrial, R, 1, thisjob.opn_model.earth_zonal_order, thisjob.opn_model.earth_tesseral_degree);

		// 转换回惯性系
		Vector3d a_nsGravity_terrestrial_evec(a_nsGravity_terrestrial[0], a_nsGravity_terrestrial[1], a_nsGravity_terrestrial[2]);
		Vector3d a_nsGravity_inertial_evec = gcrs_to_itrs_emat.transpose() * a_nsGravity_terrestrial_evec;
		axyz_nsGravity_earth[0] = a_nsGravity_inertial_evec[0];
		axyz_nsGravity_earth[1] = a_nsGravity_inertial_evec[1];
		axyz_nsGravity_earth[2] = a_nsGravity_inertial_evec[2];
	}


	//先获取日月相对地球的JPL历表，第三体摄动和光压接口仍沿用现有“目标相对地球”的输入格式。
	a = Get_JPLEPH(MJDs_t.TDB, 11, 3, rv_sun);
	a = Get_JPLEPH(MJDs_t.TDB, 10, 3, rv_moon);

	// 计算第三体摄动，太阳
	if (thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_JPL)	//1 调用数值历表
	{
		Thirdbody_Gravity_Acceleration(axyz_sunGravity, cart_array, rv_sun, SUN_EARTH_MASS_RATIO);
	}
	else if (thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_ANALYTICAL)	//2 调用分析历表，此处不支持，直接报错
	{
		// 不支持分析历表，报错
		std::string msg = "Analytical ephemeris for the Sun is not supported in EBCRS!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else if (thisjob.opn_model.solar_gravity_config == SOLAR_GRAVITY_OFF)
	{
		// 用户明确关闭太阳质点摄动。程序允许这样做，但地月空间高精度预报通常不建议关闭。
		std::string msg = "Solar point-mass gravity is OFF in EBCRS propagation.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		std::string msg = "Invalid solar_gravity_config in EBCRS propagation. Use 0 (off) or 1 (JPL).";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	// 计算第三体摄动，月亮，理论上地月空间是必加的
	if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_JPL)//1 调用数值历表
	{
		Thirdbody_Gravity_Acceleration(axyz_moonGravity, cart_array, rv_moon, MOON_EARTH_MASS_RATIO);
	}
	else if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_ANALYTICAL)	//2 调用分析历表，此处不支持，直接报错
	{
		// 不支持分析历表，报错
		std::string msg = "Analytical ephemeris for the Moon is not supported in EBCRS!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}
	else if (thisjob.opn_model.lunar_gravity_config == LUNAR_GRAVITY_OFF)
	{
		// 用户明确关闭月球质点摄动。程序允许这样做，但地月空间高精度预报通常不建议关闭。
		std::string msg = "Lunar point-mass gravity is OFF in EBCRS propagation.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
	}
	else
	{
		std::string msg = "Invalid lunar_gravity_config in EBCRS propagation. Use 0 (off) or 1 (JPL).";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	// 计算用户指定的行星质点摄动。
	// 注意：这个开关只控制经典牛顿第三体的行星项；后牛顿 EIH 修正只要开启，就会按公式自动使用全部11个大天体。
	// planetary_body_id 使用 JPL 天体编号，主程序已经排除了地球、月球、太阳和非法编号。
	if (thisjob.opn_model.use_planetary_gravity)
	{
		for (int body_id : thisjob.opn_model.planetary_body_id)
		{
			const int body_index = body_id - 1;
			double planet_ebcrs_state_norm[6] = {0.0};
			double axyz_planet_single[3] = {0.0};

			// 经典第三体摄动需要扰动天体相对地球的状态。
			// Get_JPLEPH 的 n_center = 3 直接返回行星相对地球的 EBCRS 状态，不需要先取 BCRS 再手动作差。
			a = Get_JPLEPH(MJDs_t.TDB, body_id, 3, planet_ebcrs_state_norm);

			Thirdbody_Gravity_Acceleration(axyz_planet_single, cart_array, planet_ebcrs_state_norm, MAJOR_BODY_EARTH_GM_RATIO[body_index]);
			axyz_planetGravity[0] += axyz_planet_single[0];
			axyz_planetGravity[1] += axyz_planet_single[1];
			axyz_planetGravity[2] += axyz_planet_single[2];
		}
	}


	// 计算一阶后牛顿 EIH N体修正。
	// 这里复刻师兄 Fortran 中 YHC_EBCRS_1 实际加入总加速度的 FPN1-FPN2：
	//   FPN1 = 卫星在BCRS中的EIH后牛顿修正；
	//   FPN2 = 地球在BCRS中的EIH后牛顿修正；
	//   EBCRS相对加速度 = FPN1 - FPN2。
	if (thisjob.opn_model.use_post_newtonian)
	{
		// 如果需要后牛顿 EIH 项，就一次性读取 11 个大天体相对太阳系质心的 BCRS 状态。
		// 经典行星第三体摄动不需要这一步，直接用 Get_JPLEPH(body, 3) 获取行星相对地球的 EBCRS 状态即可。
		// 这个数组的顺序与 JPL PLEPH 编号一致：1水星、2金星、3地球、4火星、5木星、6土星、
		// 7天王星、8海王星、9冥王星、10月球、11太阳。数组下标比JPL编号小1。
		for (int body_id = 1; body_id <= 11; ++body_id)
		{
			a = Get_JPLEPH(MJDs_t.TDB, body_id, 12, major_body_bcrs_state_norm[body_id - 1]);
		}

		PostNewtonian_EIH_EBCRS(tracer, major_body_bcrs_state_norm, cart_array, axyz_postNewtonian);
		if (tracer.Code() != ERR_NOERROR) return;
	}

	// 260210 新增计算月球非球形引力，必须要在外部确保加载了JPL
	if (thisjob.opn_model.use_lunar_NSG)
	{
		if (thisjob.opn_model.lunar_gravity_config != LUNAR_GRAVITY_JPL)
		{
			// 如果没有使用JPL数值历表，就无法获取月球天平动的欧拉角，无法计算月球非球形引力
			std::string msg = "To use lunar non-spherical gravity, you must use JPL numerical ephemeris for the Moon to obtain its Euler angles!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return;
		}
		else
		{
			if (!moon_gravity_field.is_ready)
			{
				// 如果月球非球形引力模型没有加载，就无法计算月球非球形引力
				std::string msg = "Lunar non-spherical gravity model not loaded!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				return;
			}

			// 第一步，计算卫星在月心惯性系下的坐标，简单的平移，单位仍然是地球归一化
			Vector3d r_moon_inertial_earth_norm(
				cart_array[0] - rv_moon[0], cart_array[1] - rv_moon[1], cart_array[2] - rv_moon[2]);

			// 第二步，计算卫星在月心固连系下的坐标，需要用到月球天平动的欧拉角，单位仍然是地球归一化
			double Omega_p = 0, i_s = 0, Lambda = 0;//月球天平动的三个欧拉角，单位为弧度
			a = Get_Libration(MJDs_t.TDB, Omega_p, i_s, Lambda); // 严格传入TDB
			Matrix3d moon_inertial_to_fixed_emat = Rz(Lambda) * Rx(i_s) * Rz(Omega_p);//月心惯性系到月心固连系的旋转矩阵
			Vector3d r_moon_fixed_earth_norm = moon_inertial_to_fixed_emat * r_moon_inertial_earth_norm;

			// 第三步，计算月球非球形引力加速度，得到月心固连系下的加速度，单位是月球归一化
			//卫星在月心固连系下的坐标，单位是月球归一化
			double r_moon_fixed_moon_norm[3] = {r_moon_fixed_earth_norm[0] * ae / moon_gravity_field.R_ref,
												r_moon_fixed_earth_norm[1] * ae / moon_gravity_field.R_ref,
												r_moon_fixed_earth_norm[2] * ae / moon_gravity_field.R_ref};

			double a_nsg_moon_fixed_moon_norm[3] = {0};
			Lunar_Nonspherical_Gravitational_Acceleration(a_nsg_moon_fixed_moon_norm, r_moon_fixed_moon_norm, thisjob.opn_model.lunar_zonal_order, thisjob.opn_model.lunar_tesseral_degree);

			// 第四步，扣除同一月球非球形场对地心的加速度，得到 E-BCRS 相对加速度。
			// 地心相对月心的位置为 -rv_moon；旋转与归一化均和卫星处直接项一致。
			Vector3d r_earth_moon_fixed_earth_norm = -moon_inertial_to_fixed_emat *
				Vector3d(rv_moon[0], rv_moon[1], rv_moon[2]);
			double r_earth_moon_fixed_moon_norm[3] = {
				r_earth_moon_fixed_earth_norm[0] * ae / moon_gravity_field.R_ref,
				r_earth_moon_fixed_earth_norm[1] * ae / moon_gravity_field.R_ref,
				r_earth_moon_fixed_earth_norm[2] * ae / moon_gravity_field.R_ref};
			double a_indirect_moon_fixed_moon_norm[3] = {0};
			Lunar_Nonspherical_Gravitational_Acceleration(a_indirect_moon_fixed_moon_norm,
				r_earth_moon_fixed_moon_norm, thisjob.opn_model.lunar_zonal_order, thisjob.opn_model.lunar_tesseral_degree);

			// 第五步，将直接项减间接项的结果旋回惯性系，单位转为 GE/ae^2。
			Vector3d a_nsg_moon_fixed_moon_norm_evec(
				a_nsg_moon_fixed_moon_norm[0] - a_indirect_moon_fixed_moon_norm[0],
				a_nsg_moon_fixed_moon_norm[1] - a_indirect_moon_fixed_moon_norm[1],
				a_nsg_moon_fixed_moon_norm[2] - a_indirect_moon_fixed_moon_norm[2]);
			Vector3d a_nsg_moon_inertial_earth_norm_evec =
				moon_inertial_to_fixed_emat.transpose() * a_nsg_moon_fixed_moon_norm_evec *
				(moon_gravity_field.GM / GE) *
				pow(ae / moon_gravity_field.R_ref, 2);

			axyz_nsGravity_moon[0] = a_nsg_moon_inertial_earth_norm_evec[0];
			axyz_nsGravity_moon[1] = a_nsg_moon_inertial_earth_norm_evec[1];
			axyz_nsGravity_moon[2] = a_nsg_moon_inertial_earth_norm_evec[2];
		}
	}


	// 计算光压摄动
	if (thisjob.opn_model.use_SRP)//1 开启光压
	{
		a = Solar_Radition_Acceleration(axyz_srp, cart_array, rv_sun, thisjob.opn_model.sat.AMR_SRP, thisjob.opn_model.sat.eta_SRP);
	}

	// 计算大气阻力摄动
	if (thisjob.opn_model.use_drag)
	{
		rho_atm = Get_AtmosphereDensity(cart_array, thisjob.refDensity, MJDs_t.UT1, thisjob.atm_list);
		AtmosphereDrag_Acceleration(axyz_atmDrag, cart_array, thisjob.opn_model.sat.Cd, thisjob.opn_model.sat.AMR_drag, rho_atm);
	}

	// FIXME:按理来说应该加入其他参数，但是先不管了，写定轨时候再加
	f[0] = cart_array[3];    //	x\dot=vx
	f[1] = cart_array[4];    //	y\dot=vy
	f[2] = cart_array[5];    //	z\dot=vz
	f[3] =
		(axyz_tb[0] + axyz_nsGravity_earth[0]) +
		(axyz_sunGravity[0] + axyz_moonGravity[0] + axyz_planetGravity[0] + axyz_nsGravity_moon[0]) +
		axyz_atmDrag[0] + axyz_srp[0] + axyz_postNewtonian[0];    //x\ddot=ax
	f[4] =
		(axyz_tb[1] + axyz_nsGravity_earth[1]) +
		(axyz_sunGravity[1] + axyz_moonGravity[1] + axyz_planetGravity[1] + axyz_nsGravity_moon[1]) +
		axyz_atmDrag[1] + axyz_srp[1] + axyz_postNewtonian[1];    //y\ddot=ay
	f[5] =
		(axyz_tb[2] + axyz_nsGravity_earth[2]) +
		(axyz_sunGravity[2] + axyz_moonGravity[2] + axyz_planetGravity[2] + axyz_nsGravity_moon[2]) +
		axyz_atmDrag[2] + axyz_srp[2] + axyz_postNewtonian[2];    //z\ddot=az

	tracer.CheckOut(__func__);
}
