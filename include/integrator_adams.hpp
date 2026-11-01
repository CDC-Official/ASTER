/**
 * @file integrator_adams.hpp
 * @brief Adams 定步长单步积分接口及单条轨道历史。
 * @details GCRS 状态使用 TT，E-BCRS 状态使用 TDB；长度尺度 ae，时间尺度 UNITTIME=sqrt(ae^3/GE)，角度 rad。
 * 单步负责起步、预测改正及历史更新，整段传播由 core_propagation 调度。
 */
#ifndef CDC_LIB_INCLUDE_INTEGRATOR_ADAMS_HPP_
#define CDC_LIB_INCLUDE_INTEGRATOR_ADAMS_HPP_

#include <array>
#include "model_propagation.hpp"

// 单条轨道的积分历史。一个目标使用一个实例；新任务应新建实例或调用 Reset。
// 不把这组数值积分状态放进物理力模型，也不在多个线程间共享。
struct AdamsHistory
{
	int steps = 0, filled = 0; // steps 为所需节点数，filled 为已保存的导数节点数。
	bool use_correction = false;
	double step_norm = 0, epoch_MJD = 0, elapsed_norm = 0; // GCRS 用 TT，E-BCRS 用 TDB；同一历史不可混用参考系或状态类型。
	std::array<double, 6> state_norm{};
	std::array<std::array<double, 6>, 9> derivatives{}; // 从最旧到最新；六个状态分量对归一化时间的导数。
	std::array<double, 9> predictor{}, corrector{}; // 历史导数的预测／改正权重，由节点数决定。

	/**
	 * @brief 清空当前轨道的网格、状态、导数历史及预测改正系数。
	 * @details 下一次单步调用将按任务设置重新建立系数，并用 RKF78 积累等间距导数节点。
	 */
	void Reset() { *this = AdamsHistory{}; }
};

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
	AdamsHistory& history, double& h_return, double& startup_TE_norm);

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
	AdamsHistory& history, double& h_return, double& startup_TE_norm);

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
	AdamsHistory& history, double& h_return, double& startup_TE_norm);

#endif // CDC_LIB_INCLUDE_INTEGRATOR_ADAMS_HPP_
