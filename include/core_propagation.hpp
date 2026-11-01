/**
 * @file core_propagation.hpp
 * @brief 整段轨道预报、完整积分节点保存及 Cartesian 轨迹插值接口。
 * @details 地球轨道采用 GCRS，输入输出历元为 UTC MJD，内部积分采用 TT；地月空间接口采用 EBCRS/TDB。
 * 轨道状态使用各传播任务规定的归一化尺度，文件写出由应用主程序完成。
 */
#ifndef CDC_CORE_PROPAGATION_HPP
#define CDC_CORE_PROPAGATION_HPP

#include "model_propagation.hpp"
#include "integrator_rkf78.hpp"
#include "integrator_adams.hpp"

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
	vector<double>& epochs, vector<NonSingular2>& records);

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
	vector<double>& epochs, vector<Cartesian>& records);


// 整段预报向输出 vector 追加初态及已接受状态；发生错误时保留先前的有效记录。
// duration_norm 决定方向，步长设置只取绝对值；零时长仅追加初态。
// GCRS 接口输入/输出历元为 UTC，内部在 TT 中积分；终点剩余时长允许 <= 1e-5 s。
void OrbitPropagation_RKF78_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob,
	vector<double>& record_MJD_UTC, vector<NonSingular2>& record_nsg2);

void OrbitPropagation_RKF78_Cart(bastion::Overwatch& tracer, const OrbitPropagationJob& thisjob,
	vector<double>& record_MJD_UTC, vector<Cartesian>& record_cart);

void OrbitPropagationKep_2Body(bastion::Overwatch& tracer, const double MJD_t0, const double MJD_t1, const Keplerian& kep0, Keplerian& kep1);


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
    std::vector<double>& epochs_utc, std::vector<orbit::NonSingular2>& nsg2_norm, std::vector<orbit::Cartesian>& cart_norm);

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
    const std::vector<double>& epochs);

// =============================================================================
// 用于 EBCRS / TDB 下的整段预报
// =============================================================================

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
	vector<double>& epochs, vector<Cartesian>& records);

/**
 * @brief Propagate an EBCRS/TDB Cartesian orbit with the RKF78 integrator.
 *
 * @param tracer Runtime event logger and error-state carrier.
 * @param thisjob EBCRS orbit propagation job. It must contain a normalized
 *                initial Cartesian state, a TDB start epoch, RKF78 settings,
 *                EOP data, and all requested force-model switches.
 * @param record_MJD_TDB Output vector of recorded epochs, expressed as MJD(TDB).
 *                       The first element is the initial epoch.
 * @param record_cart Output vector of recorded normalized EBCRS Cartesian states.
 *                    The first element is the initial state.
 *
 * @note 输出为追加语义；失败时保留初态和已经接受的记录，不追加被拒绝的状态。
 *       duration_norm 决定方向，步长只取大小；零时长仅记录初态，终点剩余时间 <= 1e-5 s。
 * @note This is the high-level propagation loop for the EBCRS executable. It
 *       repeatedly calls RKF78_Cart_EBCRS, clips the final step to the requested
 *       duration, and stores the accepted states.
 * @note The force model, including optional planetary gravity and the optional
 *       EIH post-Newtonian correction, is evaluated inside Acceleration_Cart_EBCRS.
 */
void OrbitPropagation_RKF78_Cart_EBCRS(bastion::Overwatch& tracer,
	const OrbitPropagationJob_EBCRS& thisjob,
	vector<double>& record_MJD_TDB,
	vector<Cartesian>& record_cart);
#endif // CDC_CORE_PROPAGATION_HPP
