/**
 * @file core_iod_solver.hpp
 * @brief 将 GL 初轨确定拆为首轮求解、单次观测方程、约束和残差计算。
 *
 * 观测和测站的文件读取、单位转换及适用组合检查位于尝试性主程序中。
 * 核心接收同一目标的 RA/Dec，以及逐条对应的 GCRS 测站位置，支持多站。
 * job.obs.tau 由调用方预先计算，保存各次观测相对参考历元的时间差，以 UNITTIME 归一化。
 * 状态参考历元为 job.epochMJDUTC（UTC MJD）；长度单位为 ae，时间单位为
 * UNITTIME，速度单位为 ae/UNITTIME；测角和测角残差均为 rad。
 */
#ifndef CDC_EXPERIMENTAL_IOD_HPP
#define CDC_EXPERIMENTAL_IOD_HPP

#include "core_iod.hpp"

namespace iod_solver
{
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
 *
 * @details
 * 先取得起始状态，再逐轮计算残差、保存迭代记录并检查停止条件。
 * 逐轮文字另存 solution.iteration_report，供任务报告和 GLML 所选来源显示使用。
 * 需要继续迭代时，组装全部观测的方程，通过列主元 QR 求解并更新状态。
 * 简单迭代直接求下一轮状态；微分改正求状态改正量，并可加入选定的约束。
 * 残差和 RMS 的单位为 rad。solution.flag 为正数时表示收敛轮次，
 * 为 -1 时表示达到迭代上限，为 -2 时表示数值失败；具体错误由 tracer 记录。
 * job.obs.sigma_RaDec 为东西、南北独立高斯误差共同的标准差（rad）；J_likelihood 保存该模型的负对数似然。
 */
void InitialOrbitDetermination(bastion::Overwatch &tracer, iod::IodJob &job, std::ostream* iteration_output = nullptr);

/**
 * @brief 根据初值类型，取得参考历元的起始位置和速度。
 * @param[in,out] tracer 记录初值读取或线性求解中的错误。
 * @param[in] job 提供初值类型、参考历元、归一化相对时间和逐条测站位置。
 *     历元采用 UTC MJD，测站位置为 GCRS 坐标，已除以 ae。
 *     job.obs.tau 为 N 维数组，保存各次观测相对参考历元的时间差，以 UNITTIME 归一化。
 *     使用指定初值时，job.initValue.initCart 应保存参考历元的归一化状态。
 * @param[in] L 3×N 的 GCRS 单位视线矩阵，无量纲；第 k 列对应第 k 条观测，
 *     三个分量为 [cos(Dec)*cos(RA), cos(Dec)*sin(RA), sin(Dec)]。
 * @param[out] r0v0_norm 六维 GCRS 起始状态，顺序为 [x,y,z,vx,vy,vz]。
 *     前三项为位置除以 ae，后三项为速度除以 ae/UNITTIME。
 * @param[out] iteration 起始状态对应的轮次：默认策略求解一次线性方程，取 1；直接读取指定初值时取 0。
 *
 * @details
 * 默认策略令 f=fz=1、g=gz=tau，以 r(t)≈r0+tau*v0 近似目标运动，
 * 拼接各次观测的三个基本方程，通过最小二乘求出起始状态。
 * 初值策略由 job.initValue.type 选择，可在此处增加其他初值方法。
 */
void Solve_FirstIteration(bastion::Overwatch &tracer, const iod::IodJob &job,
                          const Eigen::MatrixXd &L, Eigen::VectorXd &r0v0_norm, int &iteration);

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
 *
 * @details
 * 简单迭代每条观测使用三个基本行。微分改正按约束组合使用 3、4、5 或 6 行。
 * 各行由对应的归一化方程构造，轨道约束随每条观测重复加入。
 * 输出是供 QR 最小二乘求解的长矩阵，不是 6×6 的法矩阵。
 */
void Build_LinearSystem(bastion::Overwatch &tracer, const iod::IodJob &job,
                        const Eigen::MatrixXd &L, const Eigen::VectorXd &r0v0_norm,
                        Eigen::MatrixXd &matCoeff, Eigen::VectorXd &vecRHS);

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
 *
 * @details
 * 先根据当前参考状态计算该观测时刻的传播系数，再构造三个基本方程。
 * 微分改正按约束类型追加相应行；简单迭代只构造三个基本行。
 * 每次使用第 k 条观测对应的测站位置，因此同一目标的不同测站观测可以共同组装求解。
 */
void Build_ObservationEquations(bastion::Overwatch &tracer, const iod::IodJob &job, int k,
                                const Eigen::Vector3d &los, const Eigen::VectorXd &r0v0_norm,
                                Eigen::MatrixXd &matCoeff, Eigen::VectorXd &vecRHS);

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
void Build_BasicEquations(bastion::Overwatch &tracer, int iteration_method, double tau,
                          const Eigen::Vector3d &los, const Eigen::Vector3d &Rk_gcrs_norm,
                          const Eigen::VectorXd &r0v0_norm, const double fg[4],
                          Eigen::MatrixXd &matCoeff, Eigen::VectorXd &vecRHS);

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
void Append_ConstraintEquations(bastion::Overwatch &tracer, int constraint_type,
                                double constraint_ecc, const Eigen::Vector3d &los,
                                const Eigen::Vector3d &Rk_gcrs_norm,
                                const Eigen::VectorXd &r0v0_norm, const double fg[4],
                                Eigen::MatrixXd &matCoeff, Eigen::VectorXd &vecRHS);

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
 *
 * @details
 * 每次调用都根据输入状态重新计算传播系数与观测量，使状态、残差和 RMS 对应同一条轨道。
 * 两个方向分别对 N 个残差分量计算 RMS，合成 RMS 对全部 2N 个分量计算：
 * RMS_angle=sqrt((sum(resRa²)+sum(resDec²))/(2N))，RMS_sum 保存相同数值。
 */
void Compute_Residuals(bastion::Overwatch &tracer, iod::IodJob &job,
                       const Eigen::VectorXd &r0v0_norm, double &rms_ra_cosdec_rad,
                       double &rms_dec_rad);
} // namespace iod_solver
#endif // CDC_EXPERIMENTAL_IOD_HPP
