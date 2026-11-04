/**
 * @file integrator_adams.cpp
 * @brief GCRS/TT 的 NSG2、Cartesian 及 E-BCRS/TDB Cartesian 定步长 Adams 单步积分。
 * @details 依次检查历史、用 RKF78 积累等间距节点、预测、可选改正，再保存接受状态的末端导数。
 * 长度尺度为 ae，时间尺度为 UNITTIME=sqrt(ae^3/GE)，角度为 rad；绝对历元为 TT 或 TDB MJD。
 * 每条轨道独立保存 AdamsHistory，单步失败时保留输入状态和历史。
 */
#include "integrator_adams.hpp"
#include "integrator_rkf78.hpp"
#include "accel_complete.hpp"
#include <array>
#include <cmath>
#include <limits>

namespace {
/**
 * @brief 将后向差分积分系数展开为 Adams 预测和改正的节点权重。
 * @param[in,out] history 读取 2..9 个节点数 steps，向已清零的 predictor、corrector 累加权重；其余成员不变。
 * @details 由 nabla^m f_n = sum_l (-1)^l C(m,l) f_(n-l) 展开 Bashforth/Moulton 权重。
 * predictor[l] 配 f_(n-l)，corrector[0] 配新节点导数，corrector[l] 配 f_(n+1-l)。
 * 对每个 l 按 m 递增求和，最后施加交替符号；权重无量纲。
 */
void Build_AdamsCoefficients(AdamsHistory& history)
{
    const double gex[9] = {1., .5, 5./12, .375, 251./720, 95./288,
                           19087./60480, 5257./17280, 1070017./3628800};

    // gim 为同一差分展开的隐式改正系数。
    const double gim[9] = {1., -.5, -1./12, -1./24, -19./720, -3./160,
                           -863./60480, -275./24192, -33953./3628800};

    // m 阶后向差分为 ∇^m f_n = Σ_l (-1)^l C(m,l) f_(n-l)。
    // 因而固定 l 后累加所有 m>=l 的贡献，得到直接乘第 l 个导数的权重。
    // 调用前由 Reset 将两组数组清零；同一节点数的后续步不再重复累加。
    for (int l = 0; l < history.steps; ++l)
    {
        for (int m = l; m < history.steps; ++m)
        {
            // 小整数二项式系数，先在整数域得到精确值，再转成 double。
            int binomial = 1;

            for (int j = 1; j <= l; ++j)
            {
                binomial = binomial * (m-j+1) / j;
            }

            // 按 m 从小到大累计预测与改正权重。
            history.predictor[l] += binomial * gex[m];
            history.corrector[l] += binomial * gim[m];
        }

        // 后向差分对第 l 个节点带交替符号，偶数为正、奇数为负。
        double sign = 1.;

        if (l % 2 != 0)
        {
            sign = -1.;
        }

        // 符号统一乘回本节点的两组权重。
        history.predictor[l] *= sign;
        history.corrector[l] *= sign;
    }
}

/**
 * @brief 共用 NSG2、GCRS Cartesian 和 E-BCRS Cartesian 的 Adams 单步算法。
 * @tparam Job OrbitPropagationJob 或 OrbitPropagationJob_EBCRS，决定任务及动力学函数类型。
 * @param[in,out] tracer 记录输入、起步和动力学错误。
 * @param[in] job 提供动力学、外部资料、adams 节点数及改正开关和 rkf 起步设置。
 * @param[in] epoch_MJD 固定积分起点，GCRS 使用 TT MJD，E-BCRS 使用 TDB MJD。
 * @param[in] elapsed_norm 相对起点的经过时间，以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * @param[in,out] X0 连续六个归一化分量；NSG2 为 [a,xi,eta,h,k,lambda]，Cartesian 为 [x,y,z,vx,vy,vz]。
 * @param[in] h_norm 带方向的固定步长，采用归一化时间。
 * @param[in,out] history 当前轨道的等间距导数、网格及系数，成功时与 X0 一起提交。
 * @param[out] h_return 成功返回 h_norm，失败返回零。
 * @param[out] startup_TE_norm RKF 起步的归一化六分量局部误差绝对值之和；普通 Adams 步或失败时为 NaN。
 * @param[in] right 对应状态的动力学函数，接收绝对 MJD、维数 6 和状态，写出对归一化时间的导数。
 * @param[in] rkf 对应状态的 RKF78 单步函数，用于按同一固定步长积累历史。
 * @details 长度尺度 ae、速度尺度 ae/UNITTIME、角度 rad；right 与 rkf 必须和输入状态及时间系统匹配。
 * 积累 k 个接受终点后执行预测和可选的一次改正，最后在接受状态重新求导并更新历史。
 * 所有试算在状态和历史副本上完成，失败时 X0、history 不变。
 */
template<class Job>
void Adams_Step(bastion::Overwatch& tracer, const Job& job,
    double epoch_MJD, double elapsed_norm, double* X0, double h_norm,
    AdamsHistory& history, double& h_return, double& startup_TE_norm,
    void (*right)(bastion::Overwatch&, const Job&, double, int, const double*, double*),
    void (*rkf)(bastion::Overwatch&, const Job&, double, double*, double, double&, double&))
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 先设为失败返回值，状态和历史一起提交后才返回实际步长。
    h_return = 0.;
    startup_TE_norm = std::numeric_limits<double>::quiet_NaN();

    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    // 检查节点数、状态与时间是否有效；一步必须能在 double 精度下推进经过时间和 MJD。
    // 单步调用直接采用 h_norm 的符号；epoch_MJD 固定不变，elapsed_norm 允许为负。
    const int k = job.opn_setting.adams.steps;
    const double current_mjd = epoch_MJD + elapsed_norm * constants::Conversion::NORMALIZED_TIME_TO_DAY;

    // 时间仍由固定起点加累计经过时间构造，避免每步直接累加大数 MJD。
    const double next_elapsed = elapsed_norm + h_norm;
    const double next_mjd = epoch_MJD + next_elapsed * constants::Conversion::NORMALIZED_TIME_TO_DAY;

    // 空指针不读取；逐分量检查初态，随后与时间及设置一并验收。
    bool state_is_finite = X0 != nullptr;

    if (X0)
    {
        for (int i = 0; i < 6; ++i)
        {
            if (!std::isfinite(X0[i]))
            {
                state_is_finite = false;
                break;
            }
        }
    }

    if (k < 2 || k > 9 || !X0 || !std::isfinite(epoch_MJD) || !std::isfinite(elapsed_norm) ||
        !std::isfinite(h_norm) || h_norm == 0. || !std::isfinite(next_elapsed) ||
        !std::isfinite(current_mjd) || !std::isfinite(next_mjd) || next_elapsed == elapsed_norm ||
        next_mjd == current_mjd || !state_is_finite ||
        !std::isfinite(job.opn_setting.rkf.maxTE_norm) || job.opn_setting.rkf.maxTE_norm <= 0. ||
        job.opn_setting.rkf.iterCount_max < 1)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__,
            "Invalid Adams state, epoch, step, history size (2..9) or RKF startup settings.", bastion::ERROR);
        return;
    }

    // ==================== 01 检查并准备等间距历史 ====================
    // 多步公式要求历史等间距；换步长（包括方向）、节点数、模式或不连续的时间／状态后须重新起步。
    // 这里按精确相等判断能否接续；力模型和外部资料不在比较项内，改变它们须由调用者显式 Reset。
    // 临时更新先落在副本上，失败时保留原状态及历史，避免后续接续一条已损坏的轨道。
    AdamsHistory next = history;
    bool same_state = true;

    // 六分量逐项精确比较，判断当前状态是否接续上次接受的节点。
    for (int i = 0; i < 6; ++i)
    {
        if (next.state_norm[i] != X0[i])
        {
            same_state = false;
            break;
        }
    }

    // 相同状态也必须属于同一起点、同一步长和同一模式才能接续历史。
    if (next.steps != k || next.use_correction != job.opn_setting.adams.use_correction ||
        next.step_norm != h_norm || next.epoch_MJD != epoch_MJD ||
        next.elapsed_norm != elapsed_norm || !same_state)
    {
        // 旧网格不能接续时清空节点和系数，重新积累等间距导数。
        next.Reset();
        next.steps = k;
        next.use_correction = job.opn_setting.adams.use_correction;

        // 保存新网格的步长和绝对起点，再为当前节点数重建系数。
        next.step_norm = h_norm;
        next.epoch_MJD = epoch_MJD;
        Build_AdamsCoefficients(next);
    }

    // 历史槽位数只允许 0..k，损坏计数不能作为数组下标。
    if (next.filled < 0 || next.filled > k)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid Adams history count.", bastion::ERROR);
        return;
    }

    // result 是待接受的临时状态；end_derivative 是该状态的末端导数，均为六分量数组。
    std::array<double, 6> result, end_derivative{};

    for (int i = 0; i < 6; ++i)
    {
        result[i] = X0[i];
    }

    // 起步误差只在 RKF 成功接受时得到数值，正常 Adams 步继续保持 NaN。
    double startup_error = std::numeric_limits<double>::quiet_NaN();

    if (next.filled < k)
    {
        // ==================== 02 用 RKF78 积累历史节点 ====================
        // 每次调用走一个完整固定步长，并存入接受终点的导数。
        // 把本次重新起步的历元记为 t_start；不缓存它本身的导数，而是依次缓存
        // t_start+h,...,t_start+k*h 的导数。连续 k 次成功后，下一次调用才进入 Adams。
        // 若 RKF 需要缩步，便不能维持这套等间距历史，要求减小 Adams 固定步长后重跑。
        double accepted = 0.;
        rkf(tracer, job, current_mjd, result.data(), h_norm, accepted, startup_error);
        if (tracer.Code() != ERR_NOERROR)
        {
            return;
        }

        // 即使库调用方启用了定步长 RKF，Adams 起步仍须满足原有的误差要求。
        if (accepted != h_norm || startup_error > job.opn_setting.rkf.maxTE_norm)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "RKF startup required a shorter step or exceeded tolerance. Reduce the fixed Adams step and restart; state/history preserved.", bastion::ERROR);
            return;
        }

        // RKF 最后一个阶段的中间状态不等于其最终接受状态，不能直接把阶段导数存进历史。
        right(tracer, job, next_mjd, 6, result.data(), end_derivative.data());
        if (tracer.Code() != ERR_NOERROR) return;

        // 检查六个导数分量，保证后续改正或历史更新使用有效数值。
        for (double value : end_derivative)
        {
            if (!std::isfinite(value))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Adams derivative is not finite.", bastion::ERROR);
                return;
            }
        }

        // 先写入下一个空历史槽，再增加已填节点数，两步分开便于核对下标。
        next.derivatives[next.filled] = end_derivative;
        ++next.filled;
    }
    else
    {
        // ==================== 03 Bashforth 预测 ====================
        // x_pred = x_n + h * Σ b_l f_(n-l)。
        // 复用已有导数预测下一状态，无须再计算 RKF78 的整组阶段导数。
        // derivatives[0..k-1] 按接受顺序存放，分别为 f_(n-k+1),...,f_n；
        // predictor[l] 配 f_(n-l)，所以从数组末端按 k-l-1 取值。
        // 反向传播也保持这个存储顺序：“最新”表示刚接受的节点，其绝对历元反而更小。
        for (int l = 0; l < k; ++l)
        {
            // 外层仍按 l 递增累加，内层依次更新六个状态分量。
            for (int i = 0; i < 6; ++i)
            {
                result[i] = result[i] + h_norm * next.predictor[l] * next.derivatives[k-l-1][i];
            }
        }

        // 预测状态有限才可送回动力学，失败不会提交临时状态与历史。
        for (double value : result)
        {
            if (!std::isfinite(value))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Adams predictor is not finite.", bastion::ERROR);
                return;
            }
        }

        // 仅 PECE 进入一次改正，PE 保留预测状态并计算末端导数。
        if (next.use_correction)
        {
            // ==================== 04 可选的一次 Moulton 改正 ====================
            // PECE 的第一次 E 在预测终点求导，为改正公式提供新节点信息。
            right(tracer, job, next_mjd, 6, result.data(), end_derivative.data());
            if (tracer.Code() != ERR_NOERROR) return;

            // 检查六个导数分量，保证后续改正或历史更新使用有效数值。
            for (double value : end_derivative)
            {
                if (!std::isfinite(value))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Adams derivative is not finite.", bastion::ERROR);
                    return;
                }
            }

            // 再做一次 C：由原始 x_n 出发，将历史导数与预测终点导数按 Moulton 权重求和。
            // corrector[0] 配预测的 f_(n+1)，corrector[l] 配 f_(n+1-l)。
            // 因此这里只用 k-1 个旧导数（数组下标 1..k-1），舍去最旧的 derivatives[0]。
            // 这是一次预测—改正；这里不反复迭代到收敛，也不由预测与改正之差控制步长。
            for (int i = 0; i < 6; ++i)
            {
                result[i] = X0[i];
            }

            // 先按 l 递增累加旧节点的贡献，再加入预测终点的导数贡献。
            for (int l = 1; l < k; ++l)
            {
                // 当前旧节点的同一权重，依次作用在六个状态分量上。
                for (int i = 0; i < 6; ++i)
                {
                    result[i] = result[i] + h_norm * next.corrector[l] * next.derivatives[k-l][i];
                }
            }

            // 最后这组权重对应新节点，不能提前放进前面的求和。
            for (int i = 0; i < 6; ++i)
            {
                result[i] = result[i] + h_norm * next.corrector[0] * end_derivative[i];
            }

            // 改正后的六个分量也须有限，失败时不提交本次状态和历史。
            for (double value : result)
            {
                if (!std::isfinite(value))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Adams corrector is not finite.", bastion::ERROR);
                    return;
                }
            }
        }

        // ==================== 05 接受状态求导并更新历史 ====================
        // PE 在预测状态、PECE 在改正状态计算导数，供下一步使用。
        // 因而正常 PE 每步求一次导数，PECE 求两次；随后左移并淘汰最旧节点，
        // 再把接受终点的导数放到末端，使下一步仍可按同一组下标取 k 个等间距导数。
        right(tracer, job, next_mjd, 6, result.data(), end_derivative.data());
        if (tracer.Code() != ERR_NOERROR) return;

        // 检查六个导数分量，保证后续改正或历史更新使用有效数值。
        for (double value : end_derivative)
        {
            if (!std::isfinite(value))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Adams derivative is not finite.", bastion::ERROR);
                return;
            }
        }

        // 按从旧到新的顺序左移，数组覆盖不会破坏尚未读取的下一槽位。
        for (int l = 0; l < k-1; ++l)
        {
            next.derivatives[l] = next.derivatives[l+1];
        }

        // 左移后最后一个槽位存新接受终点的导数。
        next.derivatives[k-1] = end_derivative;
    }

    // ==================== 06 提交终点及历史 ====================
    // 所有计算成功后才一起提交，供下一步接续同一条固定网格。
    // 上面任一检查、RKF 或动力学调用失败都会提前返回；对外仍为 h_return=0、误差 NaN，
    // X0 和传入 history 保持原值，错误原因由 tracer 传回调用层。
    next.state_norm = result;
    next.elapsed_norm = next_elapsed;

    // 此前均在副本中工作，到这里才逐项写回调用方的六维状态。
    for (int i = 0; i < 6; ++i)
    {
        X0[i] = result[i];
    }

    // 状态与历史一起接受，防止下一步拿到与状态不匹配的导数记录。
    history = next;
    h_return = h_norm;
    startup_TE_norm = startup_error; // 正常 Adams 步为 NaN。

    // 对外状态、历史和返回量均已更新完毕，正常退出当前调用层。
    tracer.CheckOut(__func__);
}

}

/**
 * @brief 在 GCRS/TT 下执行一个固定步长的 NSG2 Adams 积分，历史不足时使用 RKF78 起步。
 * @param[in,out] tracer 记录输入、起步及动力学错误，失败时保留错误调用链。
 * @param[in] job 提供 GCRS 动力学、外部资料以及 opn_setting.adams 和 rkf 起步设置。
 * @param[in] epoch_MJD_TT 固定积分起点，TT 的 MJD；同一条历史的起点保持不变。
 * @param[in] elapsed_norm 当前节点相对固定起点的经过时间，以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * @param[in,out] X0 连续六个归一化状态分量 [a,xi,eta,h,k,lambda]，半长轴尺度 ae，角度 rad；成功时改写为下一节点，失败时保持原值。
 * @param[in] h_norm 非零固定步长，以 UNITTIME 归一化；正负决定单步方向。
 * @param[in,out] history 当前轨道的等间距导数和系数；成功时与状态一起更新，失败时保持原值。
 * @param[out] h_return 成功时为 h_norm，失败时为零。
 * @param[out] startup_TE_norm RKF 起步成功时为六分量归一化局部误差绝对值之和；普通 Adams 步或失败时为 NaN。
 * @details 前 k 次成功调用用 RKF78 积累 k 个接受终点的导数；随后用 Bashforth 预测，按设置执行一次 Moulton 改正。
 * PE 在预测状态、PECE 在改正状态计算末端导数，供下一步使用；起步必须接受完整步长且满足 RKF 误差限。
 * 步长、节点数、模式、起点、经过时间或状态不连续时重建历史；改变力模型或外部资料后由调用者显式 Reset。
 */
void Adams_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& job,
    double epoch_MJD_TT, double elapsed_norm, double* X0, double h_norm,
    AdamsHistory& history, double& h_return, double& startup_TE_norm)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // GCRS 的 NSG2 动力学与 RKF 起步均接收 TT 历元，状态为归一化六根数。
    Adams_Step(tracer, job, epoch_MJD_TT, elapsed_norm, X0, h_norm, history,
              h_return, startup_TE_norm, Acceleration_NSG2, RKF78_NSG2);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 在 GCRS/TT 下执行一个固定步长的 Cartesian Adams 积分，历史不足时使用 RKF78 起步。
 * @param[in,out] tracer 记录输入、起步及动力学错误，失败时保留错误调用链。
 * @param[in] job 提供 GCRS 动力学、外部资料以及 opn_setting.adams 和 rkf 起步设置。
 * @param[in] epoch_MJD_TT 固定积分起点，TT 的 MJD；同一条历史的起点保持不变。
 * @param[in] elapsed_norm 当前节点相对固定起点的经过时间，以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * @param[in,out] X0 连续六个归一化状态分量 [x,y,z,vx,vy,vz]，位置尺度 ae，速度尺度 ae/UNITTIME；成功时改写为下一节点，失败时保持原值。
 * @param[in] h_norm 非零固定步长，以 UNITTIME 归一化；正负决定单步方向。
 * @param[in,out] history 当前轨道的等间距导数和系数；成功时与状态一起更新，失败时保持原值。
 * @param[out] h_return 成功时为 h_norm，失败时为零。
 * @param[out] startup_TE_norm RKF 起步成功时为六分量归一化局部误差绝对值之和；普通 Adams 步或失败时为 NaN。
 * @details 前 k 次成功调用用 RKF78 积累 k 个接受终点的导数；随后用 Bashforth 预测，按设置执行一次 Moulton 改正。
 * PE 在预测状态、PECE 在改正状态计算末端导数，供下一步使用；起步必须接受完整步长且满足 RKF 误差限。
 * 步长、节点数、模式、起点、经过时间或状态不连续时重建历史；改变力模型或外部资料后由调用者显式 Reset。
 */
void Adams_Cart(bastion::Overwatch& tracer, const OrbitPropagationJob& job,
    double epoch_MJD_TT, double elapsed_norm, double* X0, double h_norm,
    AdamsHistory& history, double& h_return, double& startup_TE_norm)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // GCRS 的 Cartesian 动力学与 RKF 起步均接收 TT 历元，状态为归一化位置和速度。
    Adams_Step(tracer, job, epoch_MJD_TT, elapsed_norm, X0, h_norm, history,
              h_return, startup_TE_norm, Acceleration_Cart, RKF78_Cart);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 在 E-BCRS/TDB 下执行一个固定步长的 Cartesian Adams 积分，历史不足时使用 RKF78 起步。
 * @param[in,out] tracer 记录输入、起步及动力学错误，失败时保留错误调用链。
 * @param[in] job 提供 E-BCRS 动力学、外部资料以及 opn_setting.adams 和 rkf 起步设置。
 * @param[in] epoch_MJD_TDB 固定积分起点，TDB 的 MJD；同一条历史的起点保持不变。
 * @param[in] elapsed_norm 当前节点相对固定起点的经过时间，以 UNITTIME=sqrt(ae^3/GE) 归一化。
 * @param[in,out] X0 连续六个归一化状态分量 [x,y,z,vx,vy,vz]，位置尺度 ae，速度尺度 ae/UNITTIME；成功时改写为下一节点，失败时保持原值。
 * @param[in] h_norm 非零固定步长，以 UNITTIME 归一化；正负决定单步方向。
 * @param[in,out] history 当前轨道的等间距导数和系数；成功时与状态一起更新，失败时保持原值。
 * @param[out] h_return 成功时为 h_norm，失败时为零。
 * @param[out] startup_TE_norm RKF 起步成功时为六分量归一化局部误差绝对值之和；普通 Adams 步或失败时为 NaN。
 * @details 前 k 次成功调用用 RKF78 积累 k 个接受终点的导数；随后用 Bashforth 预测，按设置执行一次 Moulton 改正。
 * PE 在预测状态、PECE 在改正状态计算末端导数，供下一步使用；起步必须接受完整步长且满足 RKF 误差限。
 * 步长、节点数、模式、起点、经过时间或状态不连续时重建历史；改变力模型或外部资料后由调用者显式 Reset。
 */
void Adams_Cart_EBCRS(bastion::Overwatch& tracer, const OrbitPropagationJob_EBCRS& job,
    double epoch_MJD_TDB, double elapsed_norm, double* X0, double h_norm,
    AdamsHistory& history, double& h_return, double& startup_TE_norm)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // E-BCRS 的动力学与 RKF 起步均接收 TDB 历元，状态为归一化位置和速度。
    Adams_Step(tracer, job, epoch_MJD_TDB, elapsed_norm, X0, h_norm, history,
              h_return, startup_TE_norm, Acceleration_Cart_EBCRS, RKF78_Cart_EBCRS);
    if (tracer.Code() != ERR_NOERROR)
    {
        return;
    }

    tracer.CheckOut(__func__);
}
