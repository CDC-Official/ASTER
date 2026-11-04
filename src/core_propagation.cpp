/**
 * @file core_propagation.cpp
 * @brief 调度 RKF78、Adams 整段预报，保存轨道节点并提供 Cartesian 轨迹插值。
 * @details 地球轨道在 GCRS/TT 中积分，节点历元输出 UTC MJD；地月空间接口采用 EBCRS/TDB。
 * 完整节点可供观测插值及不同轨道表示的输出共用，插值函数不执行轨道预报。
 */

#include "core_propagation.hpp"
#include <limits>
#include <cstdio>
#include <algorithm>
#include <array>
#include <cmath>

using namespace constants::Conversion;

// =============================================================================
// 用于 GCRS / TT 下的结构和函数
// =============================================================================

void OrbitPropagation_RKF78_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob,
	vector<double>& record_MJD_UTC, vector<NonSingular2>& record_nsg2)
{
/*
 * 在外部先把初始时刻算好，initial_cart必须归一化
 * 修改在TT中积分，但是仍然传入UTC
 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	//参数初始化
	double opbt0_MJD_UTC = thisjob.opn_setting.t0_MJD_UTC;//必须是UTC op begin time t0
	if (!std::isfinite(opbt0_MJD_UTC))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial UTC epoch is not finite.", bastion::ERROR);
		return;
	}
	double opbt0_MJD_TT = UTC2TT(opbt0_MJD_UTC);
	double h0_norm = thisjob.opn_setting.h0_norm;
	double duration_norm = thisjob.opn_setting.duration_norm;
	double TOL_norm = thisjob.opn_setting.rkf.maxTE_norm;
	if (tracer.Code() != ERR_NOERROR) return;
	if (!std::isfinite(opbt0_MJD_TT) || !std::isfinite(duration_norm) ||
		!std::isfinite(h0_norm) || (duration_norm != 0.0 && h0_norm == 0.0) ||
		!std::isfinite(TOL_norm) || TOL_norm <= 0.0 || thisjob.opn_setting.rkf.iterCount_max < 1 ||
		!std::isfinite(opbt0_MJD_TT + duration_norm * NORMALIZED_TIME_TO_DAY))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid RKF78 propagation epoch, duration, step or error-control settings.", bastion::ERROR);
		return;
	}
	// 总时长的正负决定方向；步长只提供大小，避免反向预报时向相反方向循环。
	h0_norm = std::copysign(std::fabs(h0_norm), duration_norm);

	double X0[6];//初始的6D状态量，归一化的NSG2

	//记录初始值
	Keplerian initKep = thisjob.initKep_norm;
	Cartesian initCart = thisjob.initCart_norm;
	double initial_components[6];
	initCart.ToArray(initial_components);
	for (double component : initial_components)
		if (!std::isfinite(component))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial Cartesian state is not finite.", bastion::ERROR);
			return;
		}
	NonSingular2 initNsg2;

	initCart.ToNonSingular2(tracer, initNsg2);
	if (tracer.Code() != ERR_NOERROR) return;

	initCart.ToKeplerian(tracer, initKep);
	if (tracer.Code() != ERR_NOERROR) return;

	initNsg2.ToArray(X0);

	record_MJD_UTC.push_back(opbt0_MJD_UTC);
	record_nsg2.push_back(initNsg2);

	// 记录积分的中间值
	NonSingular2 temp_nsg2;

	// 在 TT 中积分，输出记录历元再转回 UTC。
	double tempdt_norm = 0;    // 从初始时刻开始，已经积分过的归一化时间长度，不是一个时刻，是一段时间！
	double tempt_MJD_TT = opbt0_MJD_TT;
	double tempt_MJD_UTC = opbt0_MJD_UTC;
	double temph_norm = h0_norm;    // 当前的步长，归一化单位

	double next_progress = 10.0;

	// 沿用 1e-5 s（10 微秒）的终点判据；最后一步按剩余时长截短，正反向相同。
	while (fabs((duration_norm - tempdt_norm) * NORMALIZED_TIME_TO_SEC) > 1E-5)
	{
		if (fabs((duration_norm - tempdt_norm) * NORMALIZED_TIME_TO_SEC) < std::fabs(h0_norm) * NORMALIZED_TIME_TO_SEC)
			h0_norm = duration_norm - tempdt_norm;

		double accepted_TE_norm; // 当前接受步的局部截断误差，不是累计轨道误差。
		RKF78_NSG2(tracer, thisjob, tempt_MJD_TT, X0, h0_norm, temph_norm, accepted_TE_norm);
		if (tracer.Code() != ERR_NOERROR) return;

		// 单步成功后才累积时间；防止浮点分辨率不足导致外层循环停滞。
		if (tempdt_norm + temph_norm == tempdt_norm)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 accepted step cannot advance elapsed time.", bastion::ERROR);
			return;
		}
		tempdt_norm += temph_norm;
		// 240925 找到每一步的MJD舍入误差始终较大问题了，是因为原本的temp_t_MJD_UTC不断累加，大数+小数：
		// tempt_MJD_UTC += temp_h * NORMALIZED_TIME_TO_DAY; 原本的错的
		tempt_MJD_TT = opbt0_MJD_TT + tempdt_norm * NORMALIZED_TIME_TO_DAY;
		tempt_MJD_UTC = TT2UTC(tempt_MJD_TT);

		temp_nsg2.Set(tracer, X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], true);
		if (tracer.Code() != ERR_NOERROR) return;

		record_MJD_UTC.push_back(tempt_MJD_UTC);
		record_nsg2.push_back(temp_nsg2);

		const double progress = tempdt_norm / duration_norm * 100;
		if (thisjob.opn_setting.display_progress >= 2 ||
			(thisjob.opn_setting.display_progress == 1 && progress + 1e-9 >= next_progress))
		{
			fprintf(stderr, "Progress: %6.2f%% | elapsed %+.2f s | step %+.2f s | TE %.2e / TOL %.2e | TT %.10f\n",
				progress, tempdt_norm * NORMALIZED_TIME_TO_SEC, temph_norm * NORMALIZED_TIME_TO_SEC,
				accepted_TE_norm, TOL_norm, tempt_MJD_TT);
			while (next_progress <= 100 && progress + 1e-9 >= next_progress) next_progress += 10;
		}

	}
	tracer.CheckOut(__func__);

}

void OrbitPropagation_RKF78_Cart(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob,
	vector<double>& record_MJD_UTC, vector<Cartesian>& record_cart)
{
/*
 * 在外部先把初始时刻算好，initial_cart必须归一化
 * 修改在TT中积分，但是仍然传入UTC
 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	//参数初始化
	double opbt0_MJD_UTC = thisjob.opn_setting.t0_MJD_UTC;//必须是UTC op begin time t0
	if (!std::isfinite(opbt0_MJD_UTC))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial UTC epoch is not finite.", bastion::ERROR);
		return;
	}
	double opbt0_MJD_TT = UTC2TT(opbt0_MJD_UTC);
	double h0_norm = thisjob.opn_setting.h0_norm;
	double duration_norm = thisjob.opn_setting.duration_norm;
	double TOL_norm = thisjob.opn_setting.rkf.maxTE_norm;
	if (tracer.Code() != ERR_NOERROR) return;
	if (!std::isfinite(opbt0_MJD_TT) || !std::isfinite(duration_norm) ||
		!std::isfinite(h0_norm) || (duration_norm != 0.0 && h0_norm == 0.0) ||
		!std::isfinite(TOL_norm) || TOL_norm <= 0.0 || thisjob.opn_setting.rkf.iterCount_max < 1 ||
		!std::isfinite(opbt0_MJD_TT + duration_norm * NORMALIZED_TIME_TO_DAY))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid RKF78 propagation epoch, duration, step or error-control settings.", bastion::ERROR);
		return;
	}
	// 总时长的正负决定方向；步长只提供大小，避免反向预报时向相反方向循环。
	h0_norm = std::copysign(std::fabs(h0_norm), duration_norm);

	double X0[6];//初始的6D状态量，归一化的cart

	//记录初始值
	Keplerian initKep = thisjob.initKep_norm;
	Cartesian initCart = thisjob.initCart_norm;
	double initial_components[6];
	initCart.ToArray(initial_components);
	for (double component : initial_components)
		if (!std::isfinite(component))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial Cartesian state is not finite.", bastion::ERROR);
			return;
		}
	initCart.ToArray(X0);

	record_MJD_UTC.push_back(opbt0_MJD_UTC);
	record_cart.push_back(initCart);

	// 记录积分的中间值
	Cartesian temp_cart;

	// 在 TT 中积分，输出记录历元再转回 UTC。
	double tempdt_norm = 0;    // 从初始时刻开始，已经积分过的归一化时间长度，不是一个时刻，是一段时间！
	double tempt_MJD_TT = opbt0_MJD_TT;
	double tempt_MJD_UTC = opbt0_MJD_UTC;
	int next_progress = 10; // 简洁模式最多输出约十条进度，不改变记录频率。
	double temph_norm = h0_norm;    // 当前的步长，归一化单位

	// 沿用 1e-5 s（10 微秒）的终点判据；最后一步按剩余时长截短，正反向相同。
	while (fabs((duration_norm - tempdt_norm) * NORMALIZED_TIME_TO_SEC) > 1E-5)
	{
		if (fabs((duration_norm - tempdt_norm) * NORMALIZED_TIME_TO_SEC) < std::fabs(h0_norm) * NORMALIZED_TIME_TO_SEC)
			h0_norm = duration_norm - tempdt_norm;

		double accepted_TE_norm; // 当前接受步的归一化局部截断误差估计，用于过程监测。
		RKF78_Cart(tracer, thisjob, tempt_MJD_TT, X0, h0_norm, temph_norm, accepted_TE_norm);
		if (tracer.Code() != ERR_NOERROR) return;

		// 单步成功后才累积时间；防止浮点分辨率不足导致外层循环停滞。
		if (tempdt_norm + temph_norm == tempdt_norm)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 accepted step cannot advance elapsed time.", bastion::ERROR);
			return;
		}
		tempdt_norm += temph_norm;
		// 240925 找到每一步的MJD舍入误差始终较大问题了，是因为原本的temp_t_MJD_UTC不断累加，大数+小数：
		// tempt_MJD_UTC += temp_h * NORMALIZED_TIME_TO_DAY; 原本的错的
		tempt_MJD_TT = opbt0_MJD_TT + tempdt_norm * NORMALIZED_TIME_TO_DAY;
		tempt_MJD_UTC = TT2UTC(tempt_MJD_TT);

		temp_cart.Set(X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], true);

		record_MJD_UTC.push_back(tempt_MJD_UTC);
		record_cart.push_back(temp_cart);
		const double progress = tempdt_norm / duration_norm * 100;
		if (thisjob.opn_setting.display_progress >= 2 ||
			(thisjob.opn_setting.display_progress == 1 && progress + 1e-9 >= next_progress))
		{
			fprintf(stderr, "Progress: %6.2f%% | elapsed %+.2f s | step %+.2f s | TE %.2e / TOL %.2e | TT %.10f\n",
				progress, tempdt_norm * NORMALIZED_TIME_TO_SEC, temph_norm * NORMALIZED_TIME_TO_SEC,
				accepted_TE_norm, TOL_norm, tempt_MJD_TT);
			while (next_progress <= 100 && progress + 1e-9 >= next_progress) next_progress += 10;
		}

	}

	tracer.CheckOut(__func__);
}

void OrbitPropagationKep_2Body(Overwatch& tracer, const double MJD_t0, const double MJD_t1, const Keplerian& kep0, Keplerian& kep1)
{
	// 使用纯二体模型，把MJD_t0时刻的归一化根数kep0，预报至MJD_t1时刻的归一化根数kep1
	// 仅考虑M就行
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double a0, e0, i0, OMEGA0, omega0, M0, M1;

	// 检查输入
	if (!kep0.IsNormalized())
	{
		std::string msg = "Normalize this Keplerian before 2body propagation!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
		return;
	}

	// 获取初始值
	a0 = kep0.GetValue(tracer, VAR_SEMIMAJORAXIS);
	e0 = kep0.GetValue(tracer, VAR_ECCENTRICITY);
	i0 = kep0.GetValue(tracer, VAR_INCLINATION);
	OMEGA0 = kep0.GetValue(tracer, VAR_RAAN);
	omega0 = kep0.GetValue(tracer, VAR_ARGPERIGEE);
	M0 = kep0.GetValue(tracer, VAR_MEANANOMALY);
	if (tracer.Code() != ERR_NOERROR) return;

	// 计算M1
	double dt_norm = (MJD_t1 - MJD_t0) * DAY_TO_NORMALIZED_TIME;
	double n0 = pow(a0, -1.5); // 归一化均角速度
	M1 = fmod(M0 + n0 * dt_norm, 2 * constants::PI);

	// 设置kep1
	kep1.Set(a0, e0, i0, OMEGA0, omega0, M1, true);

	tracer.CheckOut(__func__);
}

void OrbitPropagation_RKF78_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& thisjob,
	vector<double>& record_MJD_TDB, vector<Cartesian>& record_cart)
{
	/*
	 * 函数用途：
	 *   在 EBCRS/TDB 框架下，用 RKF78 变步长积分器完成一整段 Cartesian 轨道预报。
	 *   本函数是 EBCRS 主程序调用的高层预报入口；它负责循环推进、记录星历，
	 *   但不直接计算任何动力学项。真正的动力学右函数在 Acceleration_Cart_EBCRS 中。
	 *
	 * 输入任务 thisjob 的关键约定：
	 *   thisjob.opn_setting.t0_MJD_TDB 必须是积分起点的 MJD(TDB)。
	 *   thisjob.opn_setting.duration_norm 是总预报时长，单位是项目归一化时间。
	 *   thisjob.opn_setting.h0_norm 是初始试探步长，单位也是项目归一化时间。
	 *   thisjob.initCart_norm 必须已经是归一化的 EBCRS Cartesian 初值。
	 *
	 * 时间推进方式：
	 *   tempdt_norm 保存“已经从初始历元往前积分了多久”，它是时间长度，不是绝对时刻。
	 *   每接受一个 RKF78 单步，就用
	 *       tempt_MJD_TDB = opbt0_MJD_TDB + tempdt_norm * NORMALIZED_TIME_TO_DAY
	 *   重新从初始历元计算当前 MJD(TDB)。这样可以减少反复用大 MJD 加小步长造成的舍入误差。
	 *
	 * 与后牛顿补丁的关系：
	 *   本函数不需要知道后牛顿公式细节。只要 thisjob.opn_model.use_post_newtonian 为 true，
	 *   每个 RKF78 stage 在调用 Acceleration_Cart_EBCRS 时都会自动叠加 EIH 后牛顿修正。
	 *
	 * 输出说明：
	 *   record_MJD_TDB 会保存每个记录点的 MJD(TDB)，第一个点就是初始时刻。
	 *   record_cart 会保存对应的 EBCRS Cartesian 状态，单位保持归一化，主程序输出文件前再恢复单位。
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	//参数初始化
	double opbt0_MJD_TDB = thisjob.opn_setting.t0_MJD_TDB;//必须是TDB op begin time t0
	double h0_norm = thisjob.opn_setting.h0_norm;
	double duration_norm = thisjob.opn_setting.duration_norm;
	double TOL_norm = thisjob.opn_setting.rkf.maxTE_norm;
	if (tracer.Code() != ERR_NOERROR) return;
	if (!std::isfinite(opbt0_MJD_TDB) || !std::isfinite(duration_norm) ||
		!std::isfinite(h0_norm) || (duration_norm != 0.0 && h0_norm == 0.0) ||
		!std::isfinite(TOL_norm) || TOL_norm <= 0.0 || thisjob.opn_setting.rkf.iterCount_max < 1 ||
		!std::isfinite(opbt0_MJD_TDB + duration_norm * NORMALIZED_TIME_TO_DAY))
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid RKF78 propagation epoch, duration, step or error-control settings.", bastion::ERROR);
		return;
	}
	// 总时长的正负决定方向；步长只提供大小，避免反向预报时向相反方向循环。
	h0_norm = std::copysign(std::fabs(h0_norm), duration_norm);

	double X0[6];//初始的6D状态量，归一化的cart

	//记录初始值
	Cartesian initCart = thisjob.initCart_norm;
	double initial_components[6];
	initCart.ToArray(initial_components);
	for (double component : initial_components)
		if (!std::isfinite(component))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 initial Cartesian state is not finite.", bastion::ERROR);
			return;
		}
	initCart.ToArray(X0);

	record_MJD_TDB.push_back(opbt0_MJD_TDB);
	record_cart.push_back(initCart);

	// 记录积分的中间值
	Cartesian temp_cart;

	// 初始化用于积分的temp_t，分别是归一化单位以及MJD_TDB
	double tempdt_norm = 0;    // 从初始时刻开始，已经积分过的归一化时间长度，不是一个时刻，是一段时间！
	double tempt_MJD_TDB = opbt0_MJD_TDB;

	double temph_norm = h0_norm;    // 当前的步长，归一化单位

	double next_progress = 10.0;

	// 沿用 1e-5 s（10 微秒）的终点判据；最后一步按剩余时长截短，正反向相同。
	while (fabs((duration_norm - tempdt_norm) * NORMALIZED_TIME_TO_SEC) > 1E-5)
	{
		if (fabs((duration_norm - tempdt_norm) * NORMALIZED_TIME_TO_SEC) < std::fabs(h0_norm) * NORMALIZED_TIME_TO_SEC)
			h0_norm = duration_norm - tempdt_norm;

		double accepted_TE_norm; // 当前接受步的局部截断误差，不是累计轨道误差。
		RKF78_Cart_EBCRS(tracer, thisjob, tempt_MJD_TDB, X0, h0_norm, temph_norm, accepted_TE_norm);
		if (tracer.Code() != ERR_NOERROR) return;

		// 单步成功后才累积时间；防止浮点分辨率不足导致外层循环停滞。
		if (tempdt_norm + temph_norm == tempdt_norm)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RKF78 accepted step cannot advance elapsed time.", bastion::ERROR);
			return;
		}
		tempdt_norm += temph_norm;
		// 240925 找到每一步的MJD舍入误差始终较大问题了，是因为原本的temp_t_MJD_UTC不断累加，大数+小数：
		// tempt_MJD_UTC += temp_h * NORMALIZED_TIME_TO_DAY; 原本的错的
		tempt_MJD_TDB = opbt0_MJD_TDB + tempdt_norm * NORMALIZED_TIME_TO_DAY;

		temp_cart.Set(X0[0], X0[1], X0[2], X0[3], X0[4], X0[5], true);

		record_MJD_TDB.push_back(tempt_MJD_TDB);
		record_cart.push_back(temp_cart);
		const double progress = tempdt_norm / duration_norm * 100;
		if (thisjob.opn_setting.display_progress >= 2 ||
			(thisjob.opn_setting.display_progress == 1 && progress + 1e-9 >= next_progress))
		{
			fprintf(stderr, "Progress: %6.2f%% | elapsed %+.2f s | step %+.2f s | TE %.2e / TOL %.2e | TDB %.10f\n",
				progress, tempdt_norm * NORMALIZED_TIME_TO_SEC, temph_norm * NORMALIZED_TIME_TO_SEC,
				accepted_TE_norm, TOL_norm, tempt_MJD_TDB);
			while (next_progress <= 100 && progress + 1e-9 >= next_progress) next_progress += 10;
		}

	}

	tracer.CheckOut(__func__);
}


namespace {
/**
 * @brief 共用 GCRS 和 E-BCRS Cartesian 的固定网格 Adams 整段调度。
 * @tparam Job OrbitPropagationJob 或 OrbitPropagationJob_EBCRS。
 * @param[in,out] tracer 记录输入、积分和时间无法前进等错误。
 * @param[in] job 提供归一化 Cartesian 初态、动力学、传播时长、RKF 和 Adams 设置。
 * @param[in] epoch 动力学固定起点，GCRS 为 TT MJD，E-BCRS 为 TDB MJD。
 * @param[in] initial_record_epoch 输出初态历元，GCRS 为 UTC MJD，E-BCRS 为 TDB MJD。
 * @param[in] output_utc true 时将 TT 积分节点转为 UTC 输出；false 时保留 TDB。
 * @param[in,out] epochs 追加初态和接受节点的 UTC 或 TDB MJD，与 records 同下标配对。
 * @param[in,out] records 追加 [x,y,z,vx,vy,vz] 状态，位置尺度 ae，速度尺度 ae/UNITTIME。
 * @param[in] stepper 对应坐标系的 Adams 单步函数，维护同一份历史并返回接受步长和起步误差。
 * @param[in] rkf 对应坐标系的 RKF78 单步函数，处理不足固定步长的尾段。
 * @details duration_norm 决定方向，h0_norm 仅取大小，时间尺度 UNITTIME=sqrt(ae^3/GE)。
 * 剩余时长不超过 1e-5 s 时结束；失败时保留已追加节点，输入任务不变。
 */
template<class Job>
void Propagate_AdamsCartesian(bastion::Overwatch& tracer, const Job& job,
    double epoch, double initial_record_epoch, bool output_utc,
    vector<double>& epochs, vector<Cartesian>& records,
    void (*stepper)(bastion::Overwatch&, const Job&, double, double, double*, double,
                    AdamsHistory&, double&, double&),
    void (*rkf)(bastion::Overwatch&, const Job&, double, double*, double, double&, double&))
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    // 整段入口由 duration 决定方向，只取 h0_norm 的绝对值作为步长大小；
    // 因而负总时长配负 step，用户不必另外把 h0_norm 改为负数。
    // auto 在此只随 Job 选取现有的地球/地月设置类型。
    const auto& setting = job.opn_setting;
    const double duration = setting.duration_norm;
    const double step = std::copysign(std::fabs(setting.h0_norm), duration);

    // 初态六分量保存为局部数组，原 job 中的轨道不被积分过程改写。
    std::array<double, 6> state;
    job.initCart_norm.ToArray(state.data());

    // 六个 Cartesian 分量均须有限，避免坏初态进入动力学或写入输出。
    bool state_is_finite = true;

    for (double value : state)
    {
        if (!std::isfinite(value))
        {
            state_is_finite = false;
            break;
        }
    }

    if (!std::isfinite(epoch) || !std::isfinite(initial_record_epoch) || !std::isfinite(duration) ||
        !std::isfinite(step) || (duration != 0. && step == 0.) ||
        !std::isfinite(epoch + duration*NORMALIZED_TIME_TO_DAY) ||
        !std::isfinite(setting.rkf.maxTE_norm) || setting.rkf.maxTE_norm <= 0. ||
        setting.rkf.iterCount_max < 1 || setting.adams.steps < 2 || setting.adams.steps > 9 ||
        !state_is_finite)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid Cartesian Adams propagation settings or initial state.", bastion::ERROR);
        return;
    }

    // 向输出追加初态，随后只追加成功接受的终点；零时长只留下这条初态。
    // 中途失败直接返回，保留此前追加的有效记录，调用者仍须检查 tracer。
    epochs.push_back(initial_record_epoch);
    records.push_back(job.initCart_norm);

    // 同一份历史贯穿整段预报，不能在每次循环中重建。
    AdamsHistory history;
    double elapsed = 0.;

    // 累计时间从零开始，进度提示的下一阈值从 10% 开始。
    double next_progress = 10.;

    // 剩余时长按秒取绝对值判断，正反向使用同一终止条件；允许终点差不超过 1e-5 s。
    while (std::fabs((duration-elapsed)*NORMALIZED_TIME_TO_SEC) > 1e-5)
    {
        const double remaining = duration-elapsed;
        const bool tail = std::fabs(remaining) < std::fabs(step);

        // filled 未达到节点数时仍在起步，这个标记用于进度文字。
        const bool startup = history.filled < setting.adams.steps;

        // accepted 为本次实际前进量，error 仅由 RKF 起步或尾步提供。
        double accepted = 0.;
        double error = std::numeric_limits<double>::quiet_NaN();

        // 多步历史要求等间距；不足一个固定步长的尾段交给 RKF，不能换 h 后复用旧历史。
        // remaining 带传播方向；若自适应 RKF 只接受其中一部分，下次循环继续处理剩余尾段。
        if (tail)
        {
            rkf(tracer, job, epoch + elapsed*NORMALIZED_TIME_TO_DAY,
                state.data(), remaining, accepted, error);
        }
        else
        {
            // 整步继续固定网格；stepper 决定当前处于 RKF 起步还是 Adams 多步阶段。
            stepper(tracer, job, epoch, elapsed, state.data(), step, history, accepted, error);
        }
        if (tracer.Code() != ERR_NOERROR)
        {
            return;
        }

        // 接受步长仍须能够推进时间，避免极小步导致无限循环。
        if (elapsed + accepted == elapsed)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cartesian Adams propagation cannot advance elapsed time.", bastion::ERROR);
            return;
        }

        // 仅把成功接受的步长累加到时间，原 job 初态始终不变。
        elapsed += accepted;

        // 从整段起点加累计时长构造历元，避免反复累加大数 MJD。
        const double current_mjd = epoch + elapsed*NORMALIZED_TIME_TO_DAY;
        Cartesian orbit;
        orbit.Set(state[0], state[1], state[2], state[3], state[4], state[5], true);

        // GCRS 记录历元用 UTC，E-BCRS 则直接保留动力学使用的 TDB。
        double record_epoch = current_mjd;

        if (output_utc)
        {
            record_epoch = TT2UTC(current_mjd);
        }

        epochs.push_back(record_epoch);
        records.push_back(orbit);

        // 输出接受步的进度；Adams 正常步不伪造 RKF 的局部截断误差估计。
        const double progress = elapsed/duration*100.;

        if (setting.display_progress >= 2 || (setting.display_progress == 1 && progress+1e-9 >= next_progress))
        {
            // 依次判断尾步、起步和 PECE，其余完整步为 PE。
            const char* phase = "Adams PE";

            if (tail)
            {
                phase = "RKF78 tail";
            }
            else if (startup)
            {
                // 历史尚未积满时显示 RKF78 startup。
                phase = "RKF78 startup";
            }
            else if (setting.adams.use_correction)
            {
                // 历史充足且开启改正时显示 Adams PECE。
                phase = "Adams PECE";
            }

            // 两种 Cartesian 传播采用相同的精度和层次，动力学历元分别标明 TT 或 TDB。
            fprintf(stderr, "Progress: %6.2f%% | elapsed %+.2f s | step %+.2f s | %s | ",
                    progress, elapsed*NORMALIZED_TIME_TO_SEC, accepted*NORMALIZED_TIME_TO_SEC, phase);

            // RKF 显示局部误差；Adams 正常步明确没有相应的误差估计。
            if (std::isfinite(error))
            {
                fprintf(stderr, "RKF TE %.2e / TOL %.2e", error, setting.rkf.maxTE_norm);
            }
            else
            {
                fprintf(stderr, "TE unavailable (fixed-step Adams)");
            }

            fprintf(stderr, " | %s %.10f\n", output_utc ? "TT" : "TDB", current_mjd);

            // 跨过多个百分比节点时逐个推进，下一轮按 10% 间隔提示。
            while (next_progress <= 100 && progress+1e-9 >= next_progress)
            {
                next_progress += 10.;
            }
        }
    }

    tracer.CheckOut(__func__);
}
}

/**
 * @brief 在 GCRS/TT 下用 Adams 传播整段 Cartesian 轨道，保存初态和接受节点。
 * @param[in,out] tracer 记录设置、状态转换和积分错误；失败时保留错误调用链。
 * @param[in] job 提供归一化 initCart_norm、UTC 起点、传播时长、步长、动力学和 RKF/Adams 设置，调用中保持不变。
 * @param[in,out] epochs 追加各节点的 UTC MJD，与 records 同下标配对；负时长按历元递减顺序追加。
 * @param[in,out] records 追加归一化 GCRS Cartesian 状态 [x,y,z,vx,vy,vz]，位置尺度 ae，速度尺度 ae/UNITTIME。
 * @details duration_norm 的符号决定方向，h0_norm 只取大小，两者均以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * 同一份历史贯穿等间距起步及预测改正；不足整步的尾段用 RKF78 处理，剩余时长不超过 1e-5 s 时结束。
 * 零时长仅追加初态。失败时保留已追加的有效节点，不追加失败状态；调用者应检查 tracer。
 */
void OrbitPropagation_Adams_Cart(bastion::Overwatch& tracer, const OrbitPropagationJob& job,
    vector<double>& epochs, vector<Cartesian>& records)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    // 整段公共入口读取 UTC 起点，调用共同调度时才转为 TT。
    const double epoch_utc = job.opn_setting.t0_MJD_UTC;

    if (!std::isfinite(epoch_utc))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial Adams UTC epoch is not finite.", bastion::ERROR);
        return;
    }

    // 地球轨道在 TT 下积分；输出节点记录为 UTC 历元。
    Propagate_AdamsCartesian(tracer, job, UTC2TT(epoch_utc), epoch_utc, true,
                       epochs, records, Adams_Cart, RKF78_Cart);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 在 E-BCRS/TDB 下用 Adams 传播整段 Cartesian 轨道，保存初态和接受节点。
 * @param[in,out] tracer 记录设置、状态转换和积分错误；失败时保留错误调用链。
 * @param[in] job 提供归一化 initCart_norm、TDB 起点、传播时长、步长、动力学和 RKF/Adams 设置，调用中保持不变。
 * @param[in,out] epochs 追加各节点的 TDB MJD，与 records 同下标配对；负时长按历元递减顺序追加。
 * @param[in,out] records 追加归一化 E-BCRS Cartesian 状态 [x,y,z,vx,vy,vz]，位置尺度 ae，速度尺度 ae/UNITTIME。
 * @details duration_norm 的符号决定方向，h0_norm 只取大小，两者均以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * 同一份历史贯穿等间距起步及预测改正；不足整步的尾段用 RKF78 处理，剩余时长不超过 1e-5 s 时结束。
 * 零时长仅追加初态。失败时保留已追加的有效节点，不追加失败状态；调用者应检查 tracer。
 */
void OrbitPropagation_Adams_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& job,
    vector<double>& epochs, vector<Cartesian>& records)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    // 地月空间从输入、动力学到输出都使用 TDB，不作 UTC／TT 转换。
    const double epoch_tdb = job.opn_setting.t0_MJD_TDB;
    Propagate_AdamsCartesian(tracer, job, epoch_tdb, epoch_tdb, false,
                       epochs, records, Adams_Cart_EBCRS, RKF78_Cart_EBCRS);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 在 GCRS/TT 下用 Adams 传播整段 NSG2 轨道，保存初态和接受节点。
 * @param[in,out] tracer 记录设置、状态转换和积分错误；失败时保留错误调用链。
 * @param[in] job 提供归一化 initCart_norm、UTC 起点、传播时长、步长、动力学和 RKF/Adams 设置，调用中保持不变。
 * @param[in,out] epochs 追加各节点的 UTC MJD，与 records 同下标配对；负时长按历元递减顺序追加。
 * @param[in,out] records 追加归一化 GCRS NSG2 状态 [a,xi,eta,h,k,lambda]，半长轴尺度 ae，角度 rad。
 * @details duration_norm 的符号决定方向，h0_norm 只取大小，两者均以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * 同一份历史贯穿等间距起步及预测改正；不足整步的尾段用 RKF78 处理，剩余时长不超过 1e-5 s 时结束。
 * 零时长仅追加初态。失败时保留已追加的有效节点，不追加失败状态；调用者应检查 tracer。
 */
void OrbitPropagation_Adams_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& job,
    vector<double>& epochs, vector<NonSingular2>& records)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    // 这里是地球 NSG2 入口，设置类型明确为 OrbitPropagationSetting。
    const OrbitPropagationSetting& setting = job.opn_setting;
    const double epoch_utc = setting.t0_MJD_UTC;

    // SOFA 时间转换前先拒绝非有限输入，避免无效日期进入底层日历运算。
    if (!std::isfinite(epoch_utc))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial Adams UTC epoch is not finite.", bastion::ERROR);
        return;
    }

    // 绝对积分起点转成 TT，经过时间和步长仍为归一化量。
    const double epoch_tt = UTC2TT(epoch_utc);
    const double duration = setting.duration_norm;

    // 总时长给方向，h0_norm 只提供步长大小。
    const double step = std::copysign(std::fabs(setting.h0_norm), duration);
    double initial_cartesian[6];
    job.initCart_norm.ToArray(initial_cartesian);

    // 六个 Cartesian 分量均须有限，避免坏初态进入动力学或写入输出。
    bool state_is_finite = true;

    for (double value : initial_cartesian)
    {
        if (!std::isfinite(value))
        {
            state_is_finite = false;
            break;
        }
    }

    if (!std::isfinite(epoch_utc) || !std::isfinite(epoch_tt) || !std::isfinite(duration) ||
        !std::isfinite(step) || (duration != 0. && step == 0.) ||
        !std::isfinite(epoch_tt + duration*NORMALIZED_TIME_TO_DAY) ||
        !std::isfinite(setting.rkf.maxTE_norm) || setting.rkf.maxTE_norm <= 0. || setting.rkf.iterCount_max < 1 ||
        job.opn_setting.adams.steps < 2 || job.opn_setting.adams.steps > 9 ||
        !state_is_finite)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid Adams propagation settings or initial state.", bastion::ERROR);
        return;
    }

    // 将初始 Cartesian 转为 NSG2；后续积分直接更新六根数，记录也保存为 NSG2。
    NonSingular2 orbit;
    job.initCart_norm.ToNonSingular2(tracer, orbit);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    // NSG2 数组按 a、xi、eta、h、k、lambda 保存六个归一化定义量。
    std::array<double, 6> state;
    orbit.ToArray(state.data());

    for (double value : state)
    {
        if (!std::isfinite(value))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial NSG2 state is not finite.", bastion::ERROR);
            return;
        }
    }

    // 输出按调用顺序追加初态及已接受节点；负时长下历元递减，状态仍与同下标历元配对。
    epochs.push_back(epoch_utc);
    records.push_back(orbit);

    // 整段预报共用一份历史以连续积累固定网格导数，逐步重建会使积分始终停在 RKF 起步。
    AdamsHistory history;
    double elapsed = 0.;

    // 累计时间从零开始，进度提示的下一阈值从 10% 开始。
    double next_progress = 10.;

    // 零时长跳过循环；正反向都在剩余时长不超过 1e-5 s 时结束。
    while (std::fabs((duration-elapsed)*NORMALIZED_TIME_TO_SEC) > 1e-5)
    {
        const double remaining = duration-elapsed;
        const bool tail = std::fabs(remaining) < std::fabs(step);

        // accepted 为本次实际前进量，error 仅由 RKF 起步或尾步提供。
        double accepted = 0.;
        double error = std::numeric_limits<double>::quiet_NaN();

        // 本步开始前历史是否足够，只影响稍后显示的起步/Adams 阶段名称。
        const bool startup = history.filled < job.opn_setting.adams.steps;

        if (tail)
        {
            // 尾段不足整步时，直接改变 Adams 的 h 会破坏历史节点间距，改用 RKF78 收尾。
            // remaining 保留方向；若 RKF 缩步，下一轮继续处理尾段，不更新或复用 Adams 历史。
            // 正常主程序下尾段允许误差超限后缩步；终点仍按现有时间容差判定。
            RKF78_NSG2(tracer, job, epoch_tt + elapsed*NORMALIZED_TIME_TO_DAY,
                        state.data(), remaining, accepted, error);
        }
        else
        {
            // 完整固定步长继续交给 Adams，内部根据历史数量选择起步或多步公式。
            Adams_NSG2(tracer, job, epoch_tt, elapsed, state.data(), step, history, accepted, error);
        }

        // 本次失败时不追加记录；已经成功写入 vectors 的初态和前序节点继续保留。
        if (tracer.Code() != ERR_NOERROR)
        {
            return;
        }

        // 接受步长仍须能够推进时间，避免极小步导致无限循环。
        if (elapsed + accepted == elapsed)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Adams propagation cannot advance elapsed time.", bastion::ERROR);
            return;
        }

        // 仅累计成功接受的步长，再用临时状态保存此刻的轨道与 UTC 历元；原 job 初态不变。
        elapsed += accepted;
        const double current_tt = epoch_tt + elapsed*NORMALIZED_TIME_TO_DAY;
        orbit.Set(tracer, state[0], state[1], state[2], state[3], state[4], state[5], true);
        if (tracer.Code() != ERR_NOERROR)
        {
            return;
        }

        // 状态设置成功后才追加同下标配对的 UTC 历元与 NSG2 记录。
        epochs.push_back(TT2UTC(current_tt));
        records.push_back(orbit);

        // 输出接受步的进度；Adams 正常步不伪造 RKF 的局部截断误差估计。
        const double progress = elapsed/duration*100.;

        if (setting.display_progress >= 2 || (setting.display_progress == 1 && progress+1e-9 >= next_progress))
        {
            // 依次判断尾步、起步和 PECE，其余完整步为 PE。
            const char* phase = "Adams PE";

            if (tail)
            {
                phase = "RKF78 tail";
            }
            else if (startup)
            {
                // 历史尚未积满时显示 RKF78 startup。
                phase = "RKF78 startup";
            }
            else if (setting.adams.use_correction)
            {
                // 历史充足且开启改正时显示 Adams PECE。
                phase = "Adams PECE";
            }

            // 将每个接受步的已过时间和步长换成秒。
            fprintf(stderr, "Progress: %6.2f%% | elapsed %+.2f s | step %+.2f s | %s | ",
                    progress, elapsed*NORMALIZED_TIME_TO_SEC, accepted*NORMALIZED_TIME_TO_SEC, phase);

            // RKF 显示局部误差，普通 Adams 步说明未提供误差估计。
            if (std::isfinite(error))
            {
                fprintf(stderr, "RKF TE %.2e / TOL %.2e", error, setting.rkf.maxTE_norm);
            }
            else
            {
                fprintf(stderr, "TE unavailable (fixed-step Adams)");
            }

            // 进度显示动力学所用 TT，记录向量中保存的历元仍为 UTC。
            fprintf(stderr, " | TT %.10f\n", current_tt);

            // 跨过多个百分比节点时逐个推进，下一轮按 10% 间隔提示。
            while (next_progress <= 100 && progress+1e-9 >= next_progress)
            {
                next_progress += 10.;
            }
        }
    }

    tracer.CheckOut(__func__);
}


/**
 * @brief 按指定积分器传播 NSG2 轨道，保留完整节点及其 Cartesian 表示。
 * @param[in,out] tracer 记录积分器编号、传播、节点长度及轨道转换错误。
 * @param[in] input 提供归一化初态、UTC/TT 起点、传播时长、动力学及积分设置的任务；使用其副本。
 * @param[in] integrator 101/102 为自适应/固定 RKF78，201/202 为 Adams PE/PECE，据此设置积分模式。
 * @param[out] epochs_utc 覆盖初态及全部接受节点的 UTC MJD；覆盖传入数组。
 * @param[out] nsg2_norm 与 epochs_utc 对应的归一化 NSG2，顺序为 [a,xi,eta,h,k,lambda]，角度为 rad。
 * @param[out] cart_norm 与 epochs_utc 对应的归一化 GCRS Cartesian 状态，位置尺度 ae，速度尺度 ae/UNITTIME。
 * @details 只执行一次整弧传播，然后补齐 NSG2 近点角并转换为 Cartesian，供观测插值和不同星历输出共同使用。
 * 积分方向由传播时长决定；零时长只保留初态。失败时设置 tracer，输出可能包含已生成的部分节点。
 */
void OrbitPropagation_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& input, int integrator,
    std::vector<double>& epochs_utc, std::vector<orbit::NonSingular2>& nsg2_norm, std::vector<orbit::Cartesian>& cart_norm)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 积分器编号统一决定定步长和预测校正模式，输入任务保持不变。
    OrbitPropagationJob job = input;
    job.opn_setting.display_progress = 0;
    job.opn_setting.rkf.fixed_step = (integrator == 102);
    job.opn_setting.adams.use_correction = (integrator == 202);

    if (!(integrator == 101 || integrator == 102 || integrator == 201 || integrator == 202))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Unsupported integrator ID.", bastion::ERROR);
        return;
    }

    epochs_utc.clear();
    nsg2_norm.clear();
    cart_norm.clear();

    if (integrator == 101 || integrator == 102)
    {
        OrbitPropagation_RKF78_NSG2(tracer, job, epochs_utc, nsg2_norm);
    }
    else
    {
        OrbitPropagation_Adams_NSG2(tracer, job, epochs_utc, nsg2_norm);
    }

    if (tracer.Code() != bastion::ERR_NOERROR) return;

    if (epochs_utc.empty() || epochs_utc.size() != nsg2_norm.size())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Propagation returned empty or inconsistent orbit nodes.", bastion::ERROR);
        return;
    }

    // 同一节点保留两种表示；后续插值和文件写出均使用这些状态。
    cart_norm.resize(nsg2_norm.size());

    for (size_t i = 0; i < nsg2_norm.size(); ++i)
    {
        nsg2_norm[i].SyncAnomaly(tracer);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        nsg2_norm[i].ToCartesian(tracer, cart_norm[i]);
        if (tracer.Code() != bastion::ERR_NOERROR) return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 从已有积分节点插值得到指定历元的位置和速度。
 * @param[in,out] tracer 记录节点长度、时间顺序、非有限状态及请求历元越界错误。
 * @param[in] node_epochs_utc 积分节点的 UTC MJD，非空且严格递增。
 * @param[in] node_cart_norm 与节点历元一一对应的归一化 GCRS Cartesian 状态；位置尺度 ae，速度尺度 ae/UNITTIME。
 * @param[in] epochs 请求的 UTC MJD，非递减且位于节点覆盖范围内，允许重复历元。
 * @return 与 epochs 一一对应的归一化 GCRS Cartesian 状态；错误时返回空数组并设置 tracer。
 * @details 将 UTC 转为相对首节点的 TT 归一化时间，以附近最多四节点的位置和速度构造最高七次 Hermite 多项式。
 * 速度取同一位置多项式的导数；命中节点时直接取已有状态。本函数只插值，不调用轨道预报。
 */
std::vector<orbit::Cartesian> Interpolate_CartesianOrbit(bastion::Overwatch& tracer,
    const std::vector<double>& node_epochs_utc, const std::vector<orbit::Cartesian>& node_cart_norm,
    const std::vector<double>& epochs)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return {};

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (node_epochs_utc.empty() || node_epochs_utc.size() != node_cart_norm.size())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Orbit interpolation requires matching, non-empty nodes.", bastion::ERROR);
        return {};
    }

    // 插值要求正向且严格递增的节点；每个节点提供有限的位置和速度。
    for (size_t i = 0; i < node_epochs_utc.size(); ++i)
    {
        double values[6];
        node_cart_norm[i].ToArray(values);
        bool finite_state = true;

        for (double value : values) finite_state = finite_state && std::isfinite(value);

        if (!std::isfinite(node_epochs_utc[i]) || (i > 0 && node_epochs_utc[i] <= node_epochs_utc[i - 1]) || !finite_state)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Orbit interpolation nodes must be finite and strictly increasing.", bastion::ERROR);
            return {};
        }
    }

    const double epoch0_tt = UTC2TT(node_epochs_utc.front());

    // 以 TT 归一化时间做插值，避免在大数 MJD 上构建多项式。
    std::vector<double> t(node_epochs_utc.size());

    for (size_t i = 0; i < t.size(); ++i) t[i] = (UTC2TT(node_epochs_utc[i]) - epoch0_tt) * constants::SECONDS_PER_DAY / constants::Scale::UNIT_TIME_S;

    std::vector<orbit::Cartesian> output;

    // right 随非递减请求历元前移，不必为每条记录从头搜索积分网格。
    size_t right = 0;
    double previous = -std::numeric_limits<double>::infinity();

    for (double epoch : epochs)
    {
        if (!(std::isfinite(epoch) && epoch >= previous))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Requested epochs must be finite and ordered.", bastion::EventType::ERROR);
            return {};
        }

        // 先确认资料顺序，再转换为相对初始 TT 的归一化时间，供 Hermite 插值使用。
        previous = epoch;
        double x = (UTC2TT(epoch) - epoch0_tt) * constants::SECONDS_PER_DAY / constants::Scale::UNIT_TIME_S;

        if (!(x >= t.front() - 2e-8 && x <= t.back() + 2e-8))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Requested epoch outside propagated interval.", bastion::EventType::ERROR);
            return {};
        }

        // 找到请求历元附近的积分节点；恰好落在节点上时直接返回，免去插值误差。
        while (right + 1 < t.size() && t[right] < x) ++right;

        // 两种节点命中情况分别检查右节点与左邻节点，重复历元直接复用同一状态。
        if (std::abs(x - t[right]) < 1e-12)
        {
            output.push_back(node_cart_norm[right]);
            continue;
        }

        if (right > 0 && std::abs(x - t[right - 1]) < 1e-12)
        {
            output.push_back(node_cart_norm[right - 1]);
            continue;
        }

        if (t.size() < 2)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot interpolate a single-node trajectory.", bastion::EventType::ERROR);
            return {};
        }

        // 取附近最多四个节点，每个节点提供位置和速度两个约束，构造最高七次的 Hermite 多项式。
        // 插值位置和速度分别取该多项式的值与导数，使二者对应同一条局部轨迹。
        const size_t count = std::min<size_t>(4, t.size());
        const size_t begin = std::min(right > 1 ? right - 2 : 0, t.size() - count);
        std::array<double, 8> knots{};
        std::array<std::array<Eigen::Vector3d, 8>, 8> q;

        // 重复节点的位置值相同，一阶差商填入速度，以此建立 Hermite 差商表。
        for (size_t j = 0; j < count; ++j)
        {
            std::array<double, 6> y;
            node_cart_norm[begin + j].ToArray(y.data());

            // 将节点时间移到以本次请求历元为零点的局部时间轴，并填入重复的位置节点。
            knots[2 * j] = knots[2 * j + 1] = t[begin + j] - x;
            q[2 * j][0] = q[2 * j + 1][0] = Eigen::Vector3d(y[0], y[1], y[2]);

            // 重复节点的一阶差商直接用速度；不同节点之间仍按位置差商计算。
            q[2 * j + 1][1] = Eigen::Vector3d(y[3], y[4], y[5]);

            if (j) q[2 * j][1] = (q[2 * j][0] - q[2 * j - 1][0]) / (knots[2 * j] - knots[2 * j - 1]);
        }

        // 递推高阶差商；用请求时刻作为局部时间零点，随后同时求多项式的值与导数。
        for (size_t j = 2; j < 2 * count; ++j)
            for (size_t i = j; i < 2 * count; ++i) q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / (knots[i] - knots[i - j]);

        // 从最高阶差商反向嵌套求值，同时递推位置多项式及其对归一化时间的导数。
        Eigen::Vector3d position = q[2 * count - 1][2 * count - 1];
        Eigen::Vector3d velocity = Eigen::Vector3d::Zero();

        for (int i = static_cast<int>(2 * count) - 2; i >= 0; --i)
        {
            velocity = position - knots[i] * velocity;
            position = q[i][i] - knots[i] * position;
        }

        // Cartesian 构造器的 true 表示六个分量已经归一化，速度与位置使用同一多项式。
        output.emplace_back(position[0], position[1], position[2], velocity[0], velocity[1], velocity[2], true);
    }

    tracer.CheckOut(__func__);
    return output;
}
