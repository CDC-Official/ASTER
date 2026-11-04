/**
 * @file core_bodn.hpp
 * @brief 地球轨道六状态及可选 AMR 批处理定轨。GCRS，输入/输出历元 UTC，内部积分 TT。
 * @details 观测记录、理论观测值及几何偏导由 model_observation 提供；传播参数和动力学使用 OrbitPropagationJob。
 */
#ifndef CDC_CORE_BODN_HPP
#define CDC_CORE_BODN_HPP
#include "core_propagation.hpp"
#include "model_observation.hpp"
#include "orbit_geometry.hpp"
#include "orbit_stm.hpp"
#include <limits>
#include <string>
#include <iostream>

namespace determination
{
// 定轨策略与传播设置分开管理；传播步长、误差、力模型仍用 OrbitPropagationJob。
struct BatchSetting
{
    int integrator_id = 101; // 101 自适应 RKF78；102 固定 RKF78；201 Adams PE；202 PECE。
    int stm_method = 2; // 1 NSG2 解析近似；2 完整传播中心差分，包含已开启的摄动力。
    AnalyticStmSetting analytic_stm;
    int max_iterations = 12, min_iterations = 2;

    // 控制指标为所有采用分量的标准化残差 RMS（无量纲）；RMS 或 DRMS 任一达标即可。
    double rms_tolerance = 1e-6, drms_tolerance = 1e-8;
    bool drms_relative = false; // false: |RMS_k-RMS_(k-1)|；true: 再除以本轮 RMS。
    double fd_step = 1e-6;              // 中心差分：da=fd_step*max(1,|a|)，其余分量=fd_step。

    // 按观测类型分组统计剔除前 RMS，同类型的各测站共用一组统计量。
    double rejection_scale = 3; // 同组全部原始残差 preRMS 的倍数；0 关闭剔除。
    int rejection_start = 1; // 剔除开始轮次，从 1 计数。
    double rejection_angle_window_rad = 1; // 测角分量的硬窗口，单位 rad。
    double rejection_position_window_m = constants::Scale::UNIT_LENGTH_M; // 位置分量的硬窗口，单位 m。
    bool compute_covariance = true; // 返回已知测量 sigma 下的形式协方差，不自动按 RMS 放大。
    bool display_iterations = true;
    int estimate_parameter = 0; // 0 六状态；1 阻力 AMR；2 光压 AMR；3 自动选择（至少三圈，以 a*ae-ae=2000 km 为分界）。
    double parameter_fd_step = 1e-3; // AMR/先验AMR 比值的中心差分步长。

    /**
     * @brief 显示批处理定轨的估计、权重和迭代控制设置。
     * @param[in,out] tracer 读取已有错误状态，记录显示过程及输出错误。
     * @param[in] sigma_angle_rad 东西、南北方向共同使用的观测误差标准差，单位 rad，显示时转换为 arcsec。
     * @param[in] sigma_position_m GCRS XYZ 各分量使用的观测误差标准差，单位 m。
     * @param[in,out] output 接收 Estimation settings 文本的输出流；函数保持该流原有格式。
     * @param[in] precision 常规浮点数的有效数字位数，默认 14。
     * @details 读取待估参数模式、STM 方法及其差分或解析选项、剔除规则、迭代上下限、
     * RMS/DRMS 门限与协方差开关。相对 DRMS 模式只在本分区注明一次。
     * 本函数不修改观测权重、传播模型或迭代条件。
     */
    void Display(bastion::Overwatch& tracer, double sigma_angle_rad, double sigma_position_m, std::ostream& output = std::cerr, int precision = 14) const;
};

struct IterationRecord
{
    // iteration 从 1 开始；used_observations 是本轮接受的记录数，不是测角／XYZ 分量数。
    // rank 记录白化并作列尺度变换后设计矩阵的 QR 秩；仅满列秩时进入记录，六状态为 6，同估 AMR 为 7。
    int iteration = 0, used_observations = 0, rank = 0;

    // weighted_rms 是更新后残差除以各自 sigma 的 RMS，分母为接受的分量数，无量纲。
    // scaled_correction 为恢复列尺度后、尚未乘线搜索比例的六个归一化 NSG2 改正的最大绝对值。
    // 该统计量只包含轨道根数改正，不含 AMR。
    // step_fraction 为实际采用的线搜索比例；取 0 表示保留本轮更新前的状态。
    double weighted_rms = 0, scaled_correction = 0, step_fraction = 0;

    // parameter_amr_si 为本轮更新后的 AMR，单位 m²/kg；未估计附加参数时为 0。
    // parameter_correction 为 |delta(AMR/先验AMR)|，无量纲、无符号且尚未乘线搜索比例；未估计时为 0。
    double parameter_amr_si = 0, parameter_correction = 0;

    // 可比较时记录相邻轮更新后 WRMS 的绝对变化；相对模式再除以本轮 WRMS，均无量纲。
    double drms = std::numeric_limits<double>::quiet_NaN(); // 首轮或采用记录改变时不可比较。

    // drms_comparable 说明相邻轮采用集合相同且已有上一轮有限 WRMS；否则 drms 保持 NaN。
    // 后两项仅记录各自门限是否通过；实际停止还须满足最少轮次及粗差剔除启动条件。
    bool drms_comparable = false, rms_converged = false, drms_converged = false;
};

struct BatchResult
{
    bool converged = false;
    std::string convergence_reason; // RMS / DRMS_ABSOLUTE / DRMS_RELATIVE / ITERATION_LIMIT。
    std::string message;
    double fit_epoch_utc = 0; // 返回状态和协方差的参考历元；可在求解后传回先验历元。
    orbit::NonSingular2 fitted_nsg2_norm;
    orbit::Cartesian fitted_cart_norm;
    // 最终接受轨道的完整积分节点；三数组按节点一一对应，供插值和各种星历输出共用。
    std::vector<double> node_epochs_utc; // UTC MJD；更改输出参考历元时保持观测弧段节点不变。
    std::vector<orbit::NonSingular2> node_nsg2_norm; // 长度尺度 ae，角度 rad。
    std::vector<orbit::Cartesian> node_cart_norm; // GCRS，位置尺度 ae，速度尺度 ae/UNITTIME。
    std::vector<double> epochs_utc;
    std::vector<orbit::Cartesian> ephemeris_norm; // 与每条观测对应，包含同时刻的重复行。
    std::vector<Eigen::Vector3d> residuals; // O-C；角观测为 wrap(dRA)*cos(观测Dec)、dDec（rad）；XYZ 为 m。
    std::vector<bool> accepted; // 末轮更新前的筛选结果；最终输出时沿用，每条记录共用一个标记。
    std::vector<IterationRecord> iterations;

    // 逐轮显示及原始改正的报告文字，保留计算时的顺序，供主程序写入任务报告。
    // 数值迭代记录另存 iterations；这些文字不登记为 tracer 事件。
    std::string iteration_report;

    // 拟合历元的形式协方差：NSG2 六分量采用归一化单位；Cartesian 顺序为 r、v，单位为 m、m/s 的乘积。
    // 同估 AMR 时轨道块取联合协方差的左上角；未请求计算时保留 NaN。
    StateMatrix covariance_nsg2 = StateMatrix::Constant(std::numeric_limits<double>::quiet_NaN());
    StateMatrix covariance_cart_si = StateMatrix::Constant(std::numeric_limits<double>::quiet_NaN());

    int estimated_parameter = 0; // 最终实际选择：0 无，1 阻力，2 光压。
    double parameter_apriori_amr_si = 0, parameter_amr_si = 0;
    double parameter_sigma_amr_si = std::numeric_limits<double>::quiet_NaN();

    // 顺序为归一化 NSG2 六分量，启用附加参数时末项为 AMR(m²/kg)。含轨道与参数交叉协方差。
    Eigen::MatrixXd covariance_full;
    double weighted_rms = 0; // sqrt(sum((O-C)/sigma)^2 / 已接受分量数)，无量纲。

    /**
     * @brief 输出 BODN 完成状态、参考历元后验轨道及形式标准差。
     * @param[in,out] tracer 读取已有错误状态，记录显示过程及输出错误。
     * @param[in] setting 读取 compute_covariance，区分未计算协方差与无效方差。
     * @param[in] kep 主程序已准备的同历元后验 Keplerian，半长轴 km、角度 deg，末项为平近点角 M。
     * @param[in,out] output 接收 BODN summary 开头与后验结果；调用方在同一分区接续残差统计及资源信息。
     * @param[in] precision 状态值的有效数字位数，默认 14；标准差采用六位小数的科学记数法。
     * @details 读取 converged、convergence_reason、iterations、fit_epoch_utc、后验状态、协方差及可选 AMR。
     * Cartesian 从 ae、UNITTIME 尺度恢复为 km、km/s，标准差由 SI 协方差对角线开方后除以 1000。
     * NSG2 状态及其标准差保持归一化尺度；Keplerian 不显示尚未计算的形式标准差。
     * 形式标准差的标签分别注明 [km, km/s]、[norm] 或 [m²/kg]。
     * 本函数不传播轨道、不重新筛选观测或计算残差。要求至少有一轮迭代记录。
     */
    void PrintSummary(bastion::Overwatch& tracer, const BatchSetting& setting, const orbit::Keplerian& kep, std::ostream& output = std::cerr, int precision = 14) const;

    /**
     * @brief 写入科学报告中的结果概览及三种后验状态表。
     * @param[in,out] tracer 读取已有错误状态，记录显示过程及输出错误。
     * @param[in] setting 读取 RMS/DRMS 门限、相对 DRMS 开关及 compute_covariance。
     * @param[in] object_id 本次定轨目标编号，由主程序从选用观测中取得。
     * @param[in] initial_source 先验轨道来源说明，为 inline setting 或实际轨道文件路径。
     * @param[in] kep 已在主程序准备的同历元后验 Keplerian，半长轴为 km，角度为 deg，末项为平近点角 M。
     * @param[in,out] output 接收报告正文的输出流；函数保持该流的格式状态。
     * @param[in] precision 状态表使用的科学计数法小数位数，默认 16；UTC MJD 保留 12 位小数。
     * @details 写出第 01、02 节，读取状态、迭代和采用标记、WRMS 以及 AMR 参数类型。
     * Cartesian 表采用 m、m/s，Keplerian/NSG2 表采用 km、deg；报告的其余章节由主程序接续。
     * 本函数仅读取现有结果及已准备的 Keplerian，不执行传播或重复根数转换。要求至少有一轮迭代记录。
     */
    void PrintReport(bastion::Overwatch& tracer, const BatchSetting& setting, const std::string& object_id,
        const std::string& initial_source, const orbit::Keplerian& kep, std::ostream& output, int precision = 16) const;
};

/**
 * @brief 按观测类型的全体残差 RMS 和硬窗口，重新确定本轮采用记录。
 * @param[in,out] tracer 记录输入长度、观测类型和非有限残差错误。
 * @param[in] observations 原始顺序的计算观测，只读取 type；支持类型 1、6、7。
 * @param[in] residuals 与观测一一对应的 O-C；测角为东西、南北方向残差(rad)，惯性位置为 m。
 * @param[in] setting 读取 rejection_scale、rejection_start 及角度/位置硬窗口。
 * @param[in] iteration 当前迭代轮次，从 1 开始。
 * @return 每条记录的采用标记；任一有效分量超限则整条剔除。输入出错时返回空数组并设置 tracer。
 * @details 关闭剔除或未到启动轮次时全部采用。否则同类型跨站统计所有记录的标量 RMS，
 * 不减均值、不按 sigma 白化；分量须同时满足 scale*RMS 和硬窗口。每轮重新判断，允许重新纳入。
 */
std::vector<bool> RejectOutliersByRms(bastion::Overwatch& tracer, const std::vector<Observation>& observations,
    const std::vector<Eigen::Vector3d>& residuals, const BatchSetting& setting, int iteration);

/**
 * @brief 计算相邻两轮标准化残差 RMS 的绝对或相对变化。
 * @param[in] previous 上一轮 RMS，有限、非负且无量纲。
 * @param[in] current 本轮 RMS，有限、非负且无量纲。
 * @param[in] relative false 返回绝对变化；true 再除以本轮 RMS。
 * @return 无量纲 DRMS；相对模式双零返回 0，仅本轮零返回正无穷，非法输入返回 NaN。
 * @details 本函数只计算标量变化；两轮是否采用同一观测集合由调用方判断。
 */
double ResidualRmsChange(double previous, double current, bool relative);

/**
 * @brief 通过加权批处理最小二乘，估计首观测历元的六个 NSG2 根数及可选单个 AMR。
 * @param[in,out] tracer 记录迭代、传播及计算错误，库函数不退出进程。
 * @param[in] input 先验传播任务；initKep_norm/initCart_norm 已位于首观测 UTC 历元，
 * t0_MJD_UTC 与 t0_MJD_TT 对应同一时刻，动力学及积分参数已经准备好。
 * @param[in] observations 单目标、UTC 非递减的计算观测；类型 1/6 为 RA/Dec(rad)，类型 7 为 GCRS XYZ(m)。
 * sigma 为各分量独立的标准差，测角采用东西、南北方向的 rad 标准差；允许不同观测者具有相同历元。
 * @param[in] settings 控制积分器、STM/差分、RMS/DRMS、剔除、协方差及可选 AMR 估计。
 * @param[out] result 每次调用先重置，成功计算后保存拟合历元、归一化 NSG2/Cartesian、逐条 O-C、
 * 完整积分节点、观测历元状态、采用标记和迭代历史；可选联合协方差前六维为归一化 NSG2，末维为 AMR(m²/kg)，Cartesian 协方差为 SI。
 * @return 满足停止条件返回 true；迭代上限返回 false 并保留最近结果；计算失败返回 false 并设置 tracer。
 * @details 首次传播初始轨道，再通过线搜索传播候选轨道；接受候选时同时保留完整节点、观测历元状态和残差。
 * 下一轮复用已接受轨道，只为新的试轨和数值偏导执行必要传播；观测历元状态采用 Hermite 插值。
 * 解析偏导采用 H*G*Phi，数值模式对同一传播/插值链求中心差分；AMR 同估只支持数值模式。
 * 采用方程经列缩放 QR 求解，再用线搜索更新；不显式构造法方程。
 * 改正后的标准化 RMS 或可比较的 DRMS 任一达标，并满足最少/剔除启动轮次时停止。
 * 首轮或采用集合改变时不判断 DRMS。最终复用接受解的轨迹与残差，计算偏导及已知 sigma 下的形式协方差，
 * 不用 WRMS² 放大协方差；六状态或七参数问题均要求采用分量数多于参数数且满列秩。
 */
bool BatchOrbitDetermination_NSG2(bastion::Overwatch& tracer, const OrbitPropagationJob& input,
    const std::vector<Observation>& observations, const BatchSetting& settings, BatchResult& result);

/**
 * @brief 将后验轨道及协方差传到指定参考历元。
 * @param[in,out] tracer 记录历元检查、传播及协方差转换错误。
 * @param[in] input 提供动力学、积分设置及归一化先验 AMR 的传播任务。
 * @param[in] setting 读取 integrator_id、stm_method、analytic_stm、差分步长和 compute_covariance。
 * @param[in] reference_epoch_utc 目标 UTC MJD，允许位于观测弧段之前或之后。
 * @param[in,out] result 输入拟合历元、归一化状态及后验 AMR，更新 fit_epoch_utc 和两种轨道状态；
 * 启用协方差时同步更新联合、NSG2 和 SI Cartesian 协方差，完整积分节点、观测星历、残差及采用集合保持不变。
 * @details 轨道直接积分到目标历元；协方差采用 P_ref=J*P_fit*Jᵀ。
 * J 按当前设置取解析近似或中心差分，AMR 列恢复为对物理 m²/kg 的偏导。
 */
void MovePosteriorToEpoch(bastion::Overwatch& tracer, const OrbitPropagationJob& input,
    const BatchSetting& setting, double reference_epoch_utc, BatchResult& result);

/**
 * @brief 传播轨道并通过局部 Hermite 插值取得指定观测历元的状态。
 * @param[in,out] tracer 记录积分器编号、传播、转换、历元顺序及区间错误。
 * @param[in] input 提供初始归一化轨道、TT 起点、覆盖观测的正向传播时长、动力学及积分设置。
 * @param[in] integrator 101/102 为自适应/固定 RKF78，201/202 为 Adams PE/PECE。
 * @param[in] epochs UTC MJD 数组，非递减且位于传播区间内，允许重复历元。
 * @return 与 epochs 一一对应的归一化 GCRS Cartesian 状态；错误时返回空数组并保留 tracer 错误。
 * @details 先获得完整自由传播节点，再调用公用 Cartesian 插值函数。
 * 使用附近最多四节点的位置及速度约束最高七次多项式，速度由同一位置多项式求导；命中节点直接取值。
 */
std::vector<orbit::Cartesian> PropagateObservationEpochs(bastion::Overwatch& tracer,
    const OrbitPropagationJob& input, int integrator, const std::vector<double>& epochs);
} // namespace determination
#endif
