/**
 * @file core_bodn.cpp
 * @brief 六状态及可选阻力／光压面质比的加权批处理最小二乘。
 * @details 从先验任务准备候选轨道，在观测历元计算 O-C 和设计矩阵，按类型剔除，
 * 通过列缩放 QR 和线搜索更新初始 NSG2 及可选 AMR，最后保存同源轨迹、残差并计算形式协方差。
 * 观测历元和返回历元为 UTC MJD，内部传播使用 TT；状态使用 ae、UNITTIME 的归一化尺度。
 * 观测模型由 model_observation 提供；本文件负责估计流程、传播取样与后验历元传递。
 */
#include "orbit_geometry.hpp"
#include "core_bodn.hpp"
#include <algorithm>
#include <array>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <exception>
#include <utility>

namespace determination
{
namespace
{
/**
 * @brief 将归一化 NSG2 初态同步为传播任务中的两种轨道表示。
 * @param[in,out] tracer 记录根数定义域和状态转换错误。
 * @param[in,out] job 修改 initKep_norm、initCart_norm，其他传播设置和动力学参数不变。
 * @param[in] z 六维 [a,xi,eta,h,k,lambda]，长度尺度为 ae，角度为 rad。
 * @details 检查椭圆及 i<pi 的 NSG2 定义域，补齐近点角，再转换为 Keplerian 和 Cartesian。
 */
void Set_InitialOrbit(bastion::Overwatch& tracer, OrbitPropagationJob& job, const StateVector& z)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (!(z.allFinite() && z[0] > 0 && z[1] * z[1] + z[2] * z[2] < 1 && z[3] * z[3] + z[4] * z[4] < 1))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Correction left the supported elliptic NSG2 domain.", bastion::EventType::ERROR);
        return;
    }

    orbit::NonSingular2 n(tracer, z[0], z[1], z[2], z[3], z[4], z[5], true);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 同步偏近点角后生成两种状态表示；后续传播入口可读取一致的根数或 Cartesian。
    n.SyncAnomaly(tracer);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    n.ToKeplerian(tracer, job.initKep_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    n.ToCartesian(tracer, job.initCart_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckOut(__func__);
}

/**
 * @brief 用相对先验的比例更新所选力模型的面质比及面积。
 * @param[in,out] tracer 记录非正或非有限 AMR 比例错误。
 * @param[in,out] job 修改 opn_model.sat 中选定的 AMR 和 area，读取 mass，其他卫星参数不变。
 * @param[in] parameter 0 不修改，1 阻力，2 光压；调用方已完成自动类型选择。
 * @param[in] prior_amr_norm 选定 AMR 的归一化先验值，物理 AMR 乘 MASS_OF_EARTH/ae²。
 * @param[in] scale 当前物理 AMR 与先验 AMR 之比，无量纲且严格为正。
 * @details 同时更新 area=AMR*mass，使采用面积或面质比的力模型读取一致的参数。
 */
void Set_TrialAmr(bastion::Overwatch& tracer, OrbitPropagationJob& job, int parameter, double prior_amr_norm, double scale)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    if (parameter == 0) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (!(std::isfinite(scale) && scale > 0))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Trial AMR must be finite and positive.", bastion::EventType::ERROR);
        return;
    }

    // 按待估参数类型，同时更新对应力模型的面质比与面积。
    if (parameter == 1)
    {
        job.opn_model.sat.AMR_drag = prior_amr_norm * scale;
        job.opn_model.sat.area_drag = job.opn_model.sat.AMR_drag * job.opn_model.sat.mass;
    }
    else
    {
        job.opn_model.sat.AMR_SRP = prior_amr_norm * scale;
        job.opn_model.sat.area_SRP = job.opn_model.sat.AMR_SRP * job.opn_model.sat.mass;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 传播一条候选轨道，计算逐条 O-C 和白化残差。
 * @param[in,out] tracer 记录传播、观测计算及数值错误。
 * @param[in] input 先验传播任务，提供动力学、初始 TT 历元和归一化先验 AMR；本函数使用其副本。
 * @param[in] observations 按 UTC 非递减排序的计算观测；测角为 rad，位置为 m，sigma 同单位。
 * @param[in] epochs 与 observations 一一对应的 UTC MJD 数组。
 * @param[in] settings 读取 integrator_id 及已确定为 0/1/2 的 estimate_parameter。
 * @param[in] initial 初始归一化 NSG2 六向量 [a,xi,eta,h,k,lambda]，长度尺度 ae，时间尺度 UNITTIME。
 * @param[in] amr_scale 待估 AMR 与先验 AMR 之比；未估计参数时不使用。
 * @param[out] node_epochs_utc 本次传播完整节点的 UTC MJD。
 * @param[out] node_nsg2_norm 与节点历元对应的归一化 NSG2，长度尺度 ae，角度 rad。
 * @param[out] node_cart_norm 与节点历元对应的归一化 GCRS Cartesian，位置尺度 ae，速度尺度 ae/UNITTIME。
 * @param[out] states 逐条观测历元的归一化 GCRS Cartesian 状态，允许重复历元。
 * @param[out] residuals 已按观测条数分配的原始 O-C 数组；测角前两项为东西、南北方向残差(rad)，位置三项为 m。
 * @param[out] weighted_residual 按记录及分量展开的 (O-C)/sigma，无量纲；调用方已分配总行数。
 * @details 从 input 的副本设置候选初始根数及 AMR，传播到各观测时刻，再按 sigma 白化逐条 O-C。
 */
void Compute_OrbitResiduals(bastion::Overwatch& tracer, const OrbitPropagationJob& input,
    const std::vector<Observation>& observations, const std::vector<double>& epochs,
    const BatchSetting& settings, const StateVector& initial, double amr_scale,
    std::vector<double>& node_epochs_utc, std::vector<orbit::NonSingular2>& node_nsg2_norm,
    std::vector<orbit::Cartesian>& node_cart_norm, std::vector<orbit::Cartesian>& states,
    std::vector<Eigen::Vector3d>& residuals, Eigen::VectorXd& weighted_residual)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 从同一先验任务准备当前状态和 AMR，使每次试轨使用明确的初值与力参数。
    const int parameter = settings.estimate_parameter;
    const double prior_amr_norm = parameter == 1 ? input.opn_model.sat.AMR_drag : input.opn_model.sat.AMR_SRP;
    OrbitPropagationJob trial = input;
    Set_InitialOrbit(tracer, trial, initial);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    Set_TrialAmr(tracer, trial, parameter, prior_amr_norm, amr_scale);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    OrbitPropagation_NSG2(tracer, trial, settings.integrator_id, node_epochs_utc, node_nsg2_norm, node_cart_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    states = Interpolate_CartesianOrbit(tracer, node_epochs_utc, node_cart_norm, epochs);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 每条测角占两行，惯性位置占三行；所有记录均计算残差，采用集合由迭代流程决定。
    size_t row = 0;

    for (size_t i = 0; i < observations.size(); ++i)
    {
        // 在同一观测历元计算理论量 calc_obs，再由实测值构造 O-C。
        Eigen::Vector3d calc_obs;
        Eigen::Matrix<double, 3, 6> observation_partials;
        Compute_ObsAndPartials(tracer, observations[i], states[i], calc_obs, observation_partials);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        residuals[i].setZero();

        if (observations[i].type == ObservationType::InertialPosition)
        {
            residuals[i] = observations[i].value - calc_obs;
        }
        else
        {
            // 赤经差归入 [-pi,pi] 后乘 cos(观测Dec)，得到东西方向残差；赤纬差为南北方向残差。
            residuals[i][0] = std::remainder(observations[i].value[0] - calc_obs[0], 2 * pi) * std::cos(observations[i].value[1]);
            residuals[i][1] = observations[i].value[1] - calc_obs[1];
        }

        const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

        for (int component = 0; component < dimension; ++component)
            weighted_residual[row + component] = residuals[i][component] / observations[i].sigma[component];

        row += dimension;
    }

    if (!weighted_residual.allFinite())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite batch residuals.", bastion::EventType::ERROR);
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 构造全部观测对初始 NSG2 及可选 AMR 比例的白化设计矩阵。
 * @param[in,out] tracer 记录轨道转换、扰动传播和偏导错误。
 * @param[in] input 提供拟合初始 TT 历元、动力学及归一化先验 AMR 的传播任务。
 * @param[in] observations 按记录排列的观测及 sigma，测角采用东西、南北方向的 rad 标准差，位置采用 m。
 * @param[in] epochs 与 observations 一一对应的 UTC MJD 数组。
 * @param[in] settings 读取积分器、STM 模式、差分步长及已确定为 0/1/2 的待估参数类型。
 * @param[in] initial 当前线性化点的归一化 NSG2 六向量 [a,xi,eta,h,k,lambda]，a 以 ae 归一化，lambda 为 rad。
 * @param[in] amr_scale 当前 AMR/先验 AMR 比例，未同估时不使用。
 * @param[in] states 当前线性化点的观测历元归一化 GCRS Cartesian 状态。
 * @param[out] design_matrix 已分配的 N×6 或 N×7 矩阵；N 为全部标量观测分量数，行按记录及有效分量顺序排列。
 * 前六列对初始归一化 NSG2 求导，末列对 AMR 比例求导；每行均除以对应 sigma。
 * @details 解析模式计算 H*G*Phi；差分模式传播正负扰动轨道，再以名义轨道处的 H 映射到观测空间。
 * 设计矩阵覆盖全部记录，采用集合由批处理迭代流程另行确定。
 */
void Build_DesignMatrix(bastion::Overwatch& tracer, const OrbitPropagationJob& input,
    const std::vector<Observation>& observations, const std::vector<double>& epochs,
    const BatchSetting& settings, const StateVector& initial, double amr_scale,
    const std::vector<orbit::Cartesian>& states, Eigen::MatrixXd& design_matrix)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    const int parameter = settings.estimate_parameter;
    const int ncol = parameter ? 7 : 6;
    const double prior_amr_norm = parameter == 1 ? input.opn_model.sat.AMR_drag : input.opn_model.sat.AMR_SRP;

    if (settings.stm_method == 1)
    {
        size_t row = 0;

        for (size_t i = 0; i < observations.size(); ++i)
        {
            // H 为瞬时观测偏导（东西方向含 cos(观测Dec)），G=dx(t)/dz(t)，Phi=dz(t)/dz(t0)。
            Eigen::Vector3d calc_obs;
            Eigen::Matrix<double, 3, 6> observation_partials;
            Compute_ObsAndPartials(tracer, observations[i], states[i], calc_obs, observation_partials);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            orbit::NonSingular2 current;
            states[i].ToNonSingular2(tracer, current);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            const double elapsed_norm = (UTC2TT(epochs[i]) - input.opn_setting.t0_MJD_TT) * 86400 / UNITTIME;
            const StateMatrix geometry_partials = GeometryPartialsNSG2(tracer, current);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            const StateMatrix state_transition = AnalyticStmNSG2(initial, elapsed_norm, settings.analytic_stm);
            const Eigen::Matrix<double, 3, 6> current_element_partials = observation_partials * geometry_partials;
            const Eigen::Matrix<double, 3, 6> initial_element_partials = current_element_partials * state_transition;

            // 每行白化一次，矩阵行顺序与白化残差一致。
            const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

            for (int component = 0; component < dimension; ++component)
                design_matrix.row(row + component) = initial_element_partials.row(component) / observations[i].sigma[component];

            row += dimension;
        }
    }
    else
    {
        // 六状态需要十二条扰动轨道，同估 AMR 再加两条；对传播和插值的完整计算链求中心差分。
        for (int column = 0; column < ncol; ++column)
        {
            // AMR 的中心差分不允许跨过零；半长轴步长按 max(1,|a|) 缩放，其余根数用固定步长。
            double h;

            if (column == 6)
            {
                h = std::min(settings.parameter_fd_step * std::max(1.0, std::abs(amr_scale)), 0.5 * amr_scale);
            }
            else if (column == 0)
            {
                // 半长轴按自身尺度设置差分步长，其余五根数使用固定步长。
                h = settings.fd_step * std::max(1.0, std::abs(initial[0]));
            }
            else
            {
                h = settings.fd_step;
            }

            StateVector plus = initial;
            StateVector minus = initial;

            // 前六列扰动初始根数；末列保持初态不变，只在下面扰动 AMR 比例。
            if (column < 6)
            {
                plus[column] += h;
                minus[column] -= h;
            }

            // 从同一先验复制两份任务，分别施加本列参数的正、负扰动。
            OrbitPropagationJob plus_job = input;
            OrbitPropagationJob minus_job = input;
            Set_InitialOrbit(tracer, plus_job, plus);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            // 同步负扰动初态，再分别设置正、负轨道的 AMR 比例。
            Set_InitialOrbit(tracer, minus_job, minus);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            Set_TrialAmr(tracer, plus_job, parameter, prior_amr_norm, amr_scale + (column == 6 ? h : 0));
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            Set_TrialAmr(tracer, minus_job, parameter, prior_amr_norm, amr_scale - (column == 6 ? h : 0));
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            // 正负轨道沿相同积分器、观测历元和插值链传播，差值只反映本列待估分量。
            const std::vector<orbit::Cartesian> plus_states = PropagateObservationEpochs(tracer, plus_job, settings.integrator_id, epochs);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            const std::vector<orbit::Cartesian> minus_states = PropagateObservationEpochs(tracer, minus_job, settings.integrator_id, epochs);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            size_t difference_row = 0;

            // 先差分归一化 Cartesian，再由名义轨道的 H 映射到观测空间；不直接差分有角度绕回的 RA。
            for (size_t i = 0; i < observations.size(); ++i)
            {
                Eigen::Matrix<double, 3, 6> observation_partials;
                Eigen::Vector3d calc_obs;

                // 在名义轨道处计算几何偏导 H，用于将状态差分映射到观测空间。
                Compute_ObsAndPartials(tracer, observations[i], states[i], calc_obs, observation_partials);
                if (tracer.Code() != bastion::ERR_NOERROR) return;

                std::array<double, 6> plus_cartesian;
                plus_states[i].ToArray(plus_cartesian.data());
                std::array<double, 6> minus_cartesian;
                minus_states[i].ToArray(minus_cartesian.data());

                // state_derivative 的六行依次为位置、速度对当前待估分量的中心差分偏导。
                StateVector state_derivative;

                for (int component = 0; component < 6; ++component) state_derivative[component] = (plus_cartesian[component] - minus_cartesian[component]) / (2 * h);

                // H*dx/dp 得到 dC/dp，每行再除 sigma，填入当前待估参数对应的一列。
                const Eigen::Vector3d observation_derivative = observation_partials * state_derivative;
                const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

                for (int component = 0; component < dimension; ++component)
                    design_matrix(difference_row + component, column) = observation_derivative[component] / observations[i].sigma[component];

                difference_row += dimension;
            }
        }
    }

    if (!design_matrix.allFinite())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite batch design matrix.", bastion::EventType::ERROR);
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 将初始 NSG2 轨道向前或向后传播到指定历元。
 * @param[in,out] tracer 记录初值同步或传播错误。
 * @param[in] job 传播任务副本，提供起始 TT 历元、动力学和积分设置；调用方任务不变。
 * @param[in] initial 起始历元的归一化 NSG2 六向量 [a,xi,eta,h,k,lambda]，a 以 ae 归一化，lambda 为 rad。
 * @param[in] integrator 101/102 为自适应/固定 RKF78，201/202 为 Adams PE/PECE。
 * @param[in] target_epoch_utc 终点 UTC MJD，允许早于起点。
 * @param[out] endpoint 终点归一化 NSG2 六向量；出错时不写入。
 * @details 使用积分终点状态，不对观测弧段外的历元进行插值或外插。
 */
void Propagate_Endpoint(bastion::Overwatch& tracer, OrbitPropagationJob job, const StateVector& initial,
    int integrator, double target_epoch_utc, StateVector& endpoint)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    Set_InitialOrbit(tracer, job, initial);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    job.opn_setting.duration_norm = (UTC2TT(target_epoch_utc) - job.opn_setting.t0_MJD_TT) * 86400 / UNITTIME;
    job.opn_setting.display_progress = 0;
    job.opn_setting.rkf.fixed_step = integrator == 102;
    job.opn_setting.adams.use_correction = integrator == 202;
    std::vector<double> epochs;
    std::vector<orbit::NonSingular2> elements;

    if (integrator == 101 || integrator == 102)
    {
        OrbitPropagation_RKF78_NSG2(tracer, job, epochs, elements);
    }
    else
    {
        OrbitPropagation_Adams_NSG2(tracer, job, epochs, elements);
    }

    if (tracer.Code() != bastion::ERR_NOERROR) return;

    if (elements.empty())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Posterior epoch transfer failed.", bastion::EventType::ERROR);
        return;
    }

    elements.back().ToArray(endpoint.data());

    tracer.CheckOut(__func__);
}
} // namespace

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
    const std::vector<Eigen::Vector3d>& residuals, const BatchSetting& setting, int iteration)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return {};

    if (observations.size() != residuals.size())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Inconsistent outlier input sizes.", bastion::EventType::ERROR);
        return {};
    }

    std::vector<bool> accepted(observations.size(), true);

    if (setting.rejection_scale == 0 || iteration < setting.rejection_start) return accepted;

    std::array<double, 8> square_sum = {};
    std::array<std::size_t, 8> component_count = {};

    // 第一遍包含此前剔除的记录。RMS 不减均值，分母是同类型所有标量分量数。
    for (std::size_t i = 0; i < observations.size(); ++i)
    {
        const int type = static_cast<int>(observations[i].type);

        if (!(type == 1 || type == 6 || type == 7))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Unsupported rejection observation type.", bastion::EventType::ERROR);
            return {};
        }

        const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

        for (int j = 0; j < dimension; ++j)
        {
            if (!std::isfinite(residuals[i][j]))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite rejection residual.", bastion::EventType::ERROR);
                return {};
            }

            square_sum[type] += residuals[i][j] * residuals[i][j];
            ++component_count[type];
        }
    }

    // 第二遍按同组 RMS 和类型硬窗口判断；任一分量不合格时整条不用。
    // 位置窗口与残差均采用 m，测角窗口与东西、南北方向残差均采用 rad。
    for (std::size_t i = 0; i < observations.size(); ++i)
    {
        const int type = static_cast<int>(observations[i].type);
        const double pre_rms = std::sqrt(square_sum[type] / component_count[type]);
        const double window = type == 7 ? setting.rejection_position_window_m : setting.rejection_angle_window_rad;
        const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

        for (int j = 0; j < dimension; ++j)
            accepted[i] = accepted[i] && std::abs(residuals[i][j]) <= setting.rejection_scale * pre_rms &&
                          std::abs(residuals[i][j]) <= window;
    }

    return accepted;
}

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
    const BatchSetting& setting, double reference_epoch_utc, BatchResult& result)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    try
    {
        if (!std::isfinite(reference_epoch_utc))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite posterior reference epoch.", bastion::EventType::ERROR);
            return;
        }

        if (reference_epoch_utc == result.fit_epoch_utc)
        {
            tracer.CheckOut(__func__);
            return;
        }

        OrbitPropagationJob job = input;
        job.opn_setting.t0_MJD_UTC = result.fit_epoch_utc;
        job.opn_setting.t0_MJD_TT = UTC2TT(result.fit_epoch_utc);
        const int parameter = result.estimated_parameter;

        if (parameter)
        {
            const double prior_amr_norm = parameter == 1 ? input.opn_model.sat.AMR_drag : input.opn_model.sat.AMR_SRP;
            Set_TrialAmr(tracer, job, parameter, prior_amr_norm, result.parameter_amr_si / result.parameter_apriori_amr_si);
            if (tracer.Code() != bastion::ERR_NOERROR) return;
        }

        // 提取拟合历元的六个根数，再直接传播到指定参考历元。
        StateVector initial;
        result.fitted_nsg2_norm.ToArray(initial.data());

        StateVector endpoint;
        Propagate_Endpoint(tracer, job, initial, setting.integrator_id, reference_epoch_utc, endpoint);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        if (setting.compute_covariance)
        {
            const int dimension = parameter ? 7 : 6;
            Eigen::MatrixXd transition = Eigen::MatrixXd::Identity(dimension, dimension);

            if (setting.stm_method == 1)
            {
                const double elapsed = (UTC2TT(reference_epoch_utc) - job.opn_setting.t0_MJD_TT) * 86400 / UNITTIME;
                transition.topLeftCorner<6, 6>() = AnalyticStmNSG2(initial, elapsed, setting.analytic_stm);
            }
            else
            {
                // 对目标历元的 NSG2 求中心差分，构造协方差历元传递所需的 J。
                // 联合协方差的第七维为物理 AMR(m²/kg)，对应偏导也须恢复到这一参数尺度。
                for (int column = 0; column < dimension; ++column)
                {
                    StateVector plus = initial, minus = initial;
                    OrbitPropagationJob plus_job = job, minus_job = job;
                    double h = column == 0 ? setting.fd_step * std::max(1.0, std::abs(initial[0])) : setting.fd_step;

                    if (column == 6)
                    {
                        const double ratio = result.parameter_amr_si / result.parameter_apriori_amr_si;
                        h = std::min(setting.parameter_fd_step * std::max(1.0, std::abs(ratio)), 0.5 * ratio);
                        const double amr_norm = parameter == 1 ? input.opn_model.sat.AMR_drag : input.opn_model.sat.AMR_SRP;
                        Set_TrialAmr(tracer, plus_job, parameter, amr_norm, ratio + h);
                        if (tracer.Code() != bastion::ERR_NOERROR) return;

                        Set_TrialAmr(tracer, minus_job, parameter, amr_norm, ratio - h);
                        if (tracer.Code() != bastion::ERR_NOERROR) return;
                    }
                    else
                    {
                        plus[column] += h;
                        minus[column] -= h;
                    }

                    StateVector plus_endpoint, minus_endpoint;
                    Propagate_Endpoint(tracer, plus_job, plus, setting.integrator_id, reference_epoch_utc, plus_endpoint);
                    if (tracer.Code() != bastion::ERR_NOERROR) return;

                    Propagate_Endpoint(tracer, minus_job, minus, setting.integrator_id, reference_epoch_utc, minus_endpoint);
                    if (tracer.Code() != bastion::ERR_NOERROR) return;

                    StateVector difference = plus_endpoint - minus_endpoint;

                    // λ 是周期角；这里求局部导数，不把跨 2π 的绕回当成巨大变化。
                    difference[5] = std::remainder(difference[5], 2 * pi);
                    transition.block(0, column, 6, 1) = difference / (2 * h * (column == 6 ? result.parameter_apriori_amr_si : 1));
                }
            }

            result.covariance_full = (transition * result.covariance_full * transition.transpose()).eval();
            result.covariance_nsg2 = result.covariance_full.topLeftCorner<6, 6>();
        }

        Set_InitialOrbit(tracer, job, endpoint);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        job.initCart_norm.ToNonSingular2(tracer, result.fitted_nsg2_norm);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        result.fitted_cart_norm = job.initCart_norm;
        result.fit_epoch_utc = reference_epoch_utc;

        if (setting.compute_covariance)
        {
            StateMatrix geometry = GeometryPartialsNSG2(tracer, result.fitted_nsg2_norm);
            if (tracer.Code() != bastion::ERR_NOERROR) return;

            geometry.topRows<3>() *= ae;
            geometry.bottomRows<3>() *= ae / UNITTIME;
            result.covariance_cart_si = geometry * result.covariance_nsg2 * geometry.transpose();
        }

        tracer.LogEvent(__FILE__, __func__, __LINE__, "Posterior state and covariance moved to the apriori reference epoch.", bastion::EventType::MESSAGE, FORCE_ALERT_OFF);
        tracer.CheckOut(__func__);
    }
    catch (const std::exception& error)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), bastion::EventType::ERROR);
    }
}

/**
 * @brief 计算相邻两轮标准化残差 RMS 的绝对或相对变化。
 * @param[in] previous 上一轮 RMS，有限、非负且无量纲。
 * @param[in] current 本轮 RMS，有限、非负且无量纲。
 * @param[in] relative false 返回绝对变化；true 再除以本轮 RMS。
 * @return 无量纲 DRMS；相对模式双零返回 0，仅本轮零返回正无穷，非法输入返回 NaN。
 * @details 本函数只计算标量变化；两轮是否采用同一观测集合由调用方判断。
 */
double ResidualRmsChange(double previous, double current, bool relative)
{
    if (!std::isfinite(previous) || !std::isfinite(current) || previous < 0 || current < 0) return std::numeric_limits<double>::quiet_NaN();

    // 绝对 DRMS 比较两轮 RMS 差值；相对模式以本轮 RMS 为分母。
    const double change = std::abs(current - previous);

    if (!relative) return change;

    // 双零视为无变化，仅本轮为零且上一轮非零时返回无穷大。
    if (current > 0) return change / current;

    return previous == 0 ? 0 : std::numeric_limits<double>::infinity();
}

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
    const OrbitPropagationJob& input, int integrator, const std::vector<double>& epochs)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return {};

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    std::vector<double> node_epochs_utc;
    std::vector<orbit::NonSingular2> node_nsg2_norm;
    std::vector<orbit::Cartesian> node_cart_norm;
    OrbitPropagation_NSG2(tracer, input, integrator, node_epochs_utc, node_nsg2_norm, node_cart_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return {};

    std::vector<orbit::Cartesian> states = Interpolate_CartesianOrbit(tracer, node_epochs_utc, node_cart_norm, epochs);
    if (tracer.Code() != bastion::ERR_NOERROR) return {};

    tracer.CheckOut(__func__);
    return states;
}

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
    const std::vector<Observation>& observations, const BatchSetting& settings, BatchResult& result)
{
    // 每次调用从空结果开始；错误通过 tracer 和 result.message 返回，不退出进程。
    result = BatchResult{};
    if (tracer.Code() != bastion::ERR_NOERROR)
    {
        result.message = "Batch OD called with an existing error.";
        return false;
    }

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    try
    {
        // ==================== 01 检查迭代设置、历元和观测资料 ====================

        if (observations.empty())
        {
            result.message = "No observations.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        if (!(settings.max_iterations >= settings.min_iterations && settings.min_iterations >= 1 &&
                std::isfinite(settings.fd_step) && settings.fd_step > 0 && settings.fd_step < 1e-2 &&
                std::isfinite(settings.rms_tolerance) && settings.rms_tolerance >= 0 &&
                std::isfinite(settings.drms_tolerance) && settings.drms_tolerance >= 0 &&
                std::isfinite(settings.rejection_scale) && settings.rejection_scale >= 0 &&
                settings.rejection_start >= 1 &&
                std::isfinite(settings.rejection_angle_window_rad) && settings.rejection_angle_window_rad > 0 &&
                std::isfinite(settings.rejection_position_window_m) && settings.rejection_position_window_m > 0 && (settings.stm_method == 1 || settings.stm_method == 2)))
        {
            result.message = "Invalid batch iteration settings.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        // 初始轨道的 UTC 历元应与首条观测一致，批处理在这个历元估计初始根数。
        if (!(std::abs(observations.front().mjd_utc - input.opn_setting.t0_MJD_UTC) * 86400 < 1e-5))
        {
            result.message = "Initial orbit must be at the first observation epoch.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        // 传播使用 TT，输入任务中的 UTC 和 TT 起始历元应表示同一时刻。
        if (!(std::isfinite(input.opn_setting.t0_MJD_TT) &&
                std::abs(UTC2TT(input.opn_setting.t0_MJD_UTC) - input.opn_setting.t0_MJD_TT) * 86400 < 1e-6))
        {
            result.message = "Initial UTC/TT epochs are inconsistent.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        // 一条测角贡献两行方程，位置资料贡献三行；检查各分量的数值和标准差。
        std::vector<double> epochs;
        size_t nrow = 0;
        double previous = -1e100;

        for (const Observation& obs : observations)
        {
            if (!(obs.mjd_utc >= previous && std::isfinite(obs.mjd_utc)))
            {
                result.message = "Unordered/non-finite observation epoch.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            // 同时刻的多站观测允许重复历元，但整弧必须属于同一空间目标。
            previous = obs.mjd_utc;

            if (!(obs.object_id == observations.front().object_id))
            {
                result.message = "Batch OD supports one target per job.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            if (!(obs.value.allFinite() && obs.observer_position_gcrs_m.allFinite()))
            {
                result.message = "Non-finite observation value.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            // 限定为本观测模型支持的两类测角和惯性位置，随后才按类型检查分量。
            if (!(obs.type == ObservationType::GroundRaDec || obs.type == ObservationType::SpaceRaDec ||
                    obs.type == ObservationType::InertialPosition))
            {
                result.message = "Unsupported observation type.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            if (obs.type != ObservationType::InertialPosition)
            {
                if (!(std::abs(obs.value[1]) < pi / 2))
                {
                    result.message = "Declination must be between -90 and 90 deg.";
                    tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                    return false;
                }
            }

            // sigma 与该观测残差采用同一单位；严格为正才能对各行做白化。
            const int dimension = obs.type == ObservationType::InertialPosition ? 3 : 2;

            for (int j = 0; j < dimension; ++j)
            {
                if (!(std::isfinite(obs.sigma[j]) && obs.sigma[j] > 0))
                {
                    result.message = "Measurement sigma must be finite and positive.";
                    tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                    return false;
                }
            }

            // 保留逐条 UTC 历元，并累计展开后的方程行数；重复历元不能合并掉观测。
            epochs.push_back(obs.mjd_utc);
            nrow += dimension;
        }

        // AMR 差分步长是无量纲比例的步长；先检查设置，再决定是否实际增加第七列。
        if (!(settings.estimate_parameter >= 0 && settings.estimate_parameter <= 3 &&
                std::isfinite(settings.parameter_fd_step) && settings.parameter_fd_step > 0 &&
                settings.parameter_fd_step < 0.1))
        {
            result.message = "Invalid additional parameter settings.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        // ==================== 02 准备初始根数、AMR 与计算数组 ====================

        // 以首观测历元的 NSG2 为待估参数 z，传播区间覆盖全部资料。
        OrbitPropagationJob job = input;
        job.opn_setting.duration_norm = (UTC2TT(epochs.back()) - job.opn_setting.t0_MJD_TT) * 86400 / UNITTIME;

        // 将首观测历元的初轨转为归一化 NSG2，作为本次批处理的六个待估轨道参数 z。
        orbit::NonSingular2 initial;
        job.initKep_norm.ToNonSingular2(tracer, initial);
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            result.message = "Failed to convert the initial orbit to NSG2.";
            return false;
        }

        StateVector z;
        initial.ToArray(z.data());

        // 自动模式在弧长达到三圈后增加一个力参数，并以 a*ae-ae 表示的高度作选择。
        // 高度达到 2000 km 时选光压，较低时选已开启的阻力；显式选择则使用指定的参数类型。
        int parameter = settings.estimate_parameter;

        if (parameter == 3)
        {
            const double period = 2 * pi * std::pow(z[0], 1.5) * UNITTIME;
            parameter = 0;

            // 先要求至少覆盖三圈；满足后以平均高度选光压或已开启的阻力，短弧保留六状态。
            if ((epochs.back() - epochs.front()) * 86400 >= 3 * period) parameter = (z[0] - 1) * ae >= 2000000 ? 2 : (job.opn_model.use_drag ? 1 : 0);
        }

        const int ncol = parameter ? 7 : 6;

        if (nrow <= static_cast<size_t>(ncol))
        {
            result.message = "Need more observation components than estimated parameters.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        // 先验 AMR 保留归一化值供传播，同时恢复 SI 值供最终物理参数及协方差输出。
        const model::common::SphereSatellite& sat = job.opn_model.sat;
        const double prior_amr_norm = parameter == 1 ? sat.AMR_drag : sat.AMR_SRP;
        const double prior_amr_si = parameter ? prior_amr_norm * ae * ae / MASS_OF_EARTH : 0;

        // 同估力参数只支持完整传播差分，且必须启用对应力模型；解析近似没有该参数列。
        if (parameter)
        {
            if (settings.stm_method != 2)
            {
                result.message = "AMR estimation requires stm_method = 2 (full trajectory differences).";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            if (!(parameter == 1 ? job.opn_model.use_drag : job.opn_model.use_SRP))
            {
                result.message = "Estimated force parameter requires its force model to be enabled.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            // 比例参数要求正先验 AMR 与质量；固定的阻力/反射系数也必须处于当前模型定义域。
            if (!(std::isfinite(prior_amr_si) && prior_amr_si > 0 &&
                    std::isfinite(sat.mass) && sat.mass > 0 &&
                    (parameter == 1 ? std::isfinite(sat.Cd) && sat.Cd > 0 :
                    std::isfinite(sat.eta_SRP) && sat.eta_SRP >= 0)))
            {
                result.message = "AMR estimation requires positive prior AMR/mass and a valid fixed force coefficient.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }
        }

        // 保存本次实际估计的参数类型，以先验 AMR 为基准将初始比例设为 1。
        result.estimated_parameter = parameter;
        result.parameter_apriori_amr_si = prior_amr_si;
        double parameter_scale = 1;

        // 先验值保持在 job 中；试轨设置仅把自动选择结果固定为本弧实际使用的 0/1/2。
        // 后续残差和偏导计算均使用这一已确定的估计参数类型。
        BatchSetting trial_settings = settings;
        trial_settings.estimate_parameter = parameter;
        std::vector<bool> accepted(observations.size(), true);

        // 名义轨道与候选轨道使用不同的残差数组；设计矩阵只属于本轮线性化点。
        std::vector<double> node_epochs_utc, candidate_node_epochs_utc;
        std::vector<orbit::NonSingular2> node_nsg2_norm, candidate_node_nsg2_norm;
        std::vector<orbit::Cartesian> node_cart_norm, candidate_node_cart_norm;
        std::vector<orbit::Cartesian> states, candidate_states;
        std::vector<Eigen::Vector3d> residuals(observations.size()), candidate_residuals(observations.size());
        Eigen::VectorXd weighted_residual(nrow), candidate_weighted_residual(nrow);
        Eigen::MatrixXd design_matrix(nrow, ncol);

        // DRMS 必须比较相邻两轮完成改正后的同一组资料；首轮没有可比较的上一轮。
        double previous_rms = std::numeric_limits<double>::quiet_NaN();
        std::vector<bool> previous_accepted;

        // ==================== 03 迭代构造方程、求改正并判断收敛 ====================

        // 先传播初始轨道；每轮接受候选解后，轨迹和残差一并保留给下一轮。
        Compute_OrbitResiduals(tracer, job, observations, epochs, trial_settings, z, parameter_scale,
                               node_epochs_utc, node_nsg2_norm, node_cart_norm, states, residuals, weighted_residual);
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            result.message = "Failed to compute batch orbit residuals.";
            return false;
        }

        for (int iteration = 1; iteration <= settings.max_iterations; ++iteration)
        {
            // 当前轨迹与 z、AMR 同步；仅为构造设计矩阵计算必要的扰动轨道。
            Build_DesignMatrix(tracer, job, observations, epochs, trial_settings, z, parameter_scale, states, design_matrix);
            if (tracer.Code() != bastion::ERR_NOERROR)
            {
                result.message = "Failed to build the batch design matrix.";
                return false;
            }

            // 按观测类型，用当轮全部原始残差的 preRMS 和硬窗口重新筛选。
            // 筛选用原始 rad/m 残差，不使用除以输入 sigma 后的白化残差。
            accepted = RejectOutliersByRms(tracer, observations, residuals, settings, iteration);
            if (tracer.Code() != bastion::ERR_NOERROR)
            {
                result.message = "Failed to screen batch residuals.";
                return false;
            }

            Eigen::MatrixXd A;
            Eigen::VectorXd r;

            // 先数出接受记录贡献的方程行数，再按原顺序复制；整条剔除不会留下半条测角。
            size_t selected_rows = 0;

            for (size_t i = 0; i < accepted.size(); ++i)
                if (accepted[i]) selected_rows += (observations[i].type == ObservationType::InertialPosition ? 3 : 2);

            // 剩余分量必须多于未知数；否则无法继续当前估计问题。
            if (selected_rows <= static_cast<size_t>(ncol))
            {
                result.message = "Too few components remain after rejection.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            A.resize(selected_rows, ncol);
            r.resize(selected_rows);

            // from_row 指向全部观测方程中的源行，to_row 指向采用方程中的目标行。
            size_t from_row = 0, to_row = 0;

            for (size_t i = 0; i < accepted.size(); ++i)
            {
                const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

                // 采用一条记录时，设计矩阵与残差同步复制，保持各行的统计含义一致。
                if (accepted[i])
                {
                    A.middleRows(to_row, dimension) = design_matrix.middleRows(from_row, dimension);
                    r.segment(to_row, dimension) = weighted_residual.segment(from_row, dimension);

                    // 完成一条采用记录的复制后，将目标行移到下一条采用记录的位置。
                    to_row += dimension;
                }

                from_row += dimension;
            }

            // A 的每行对应一个已采用的标量观测，列对应初始根数及可选 AMR 比例；r 为同序的白化 O-C。
            // 每行已除以对应 sigma，本轮求解加权最小二乘问题 min ||A*dz-r||²。
            // 实际令 S=diag(scales)，用 QR 解 (A*S)*u=r，再以 dz=S*u 恢复改正，避免显式构造 A^T A。
            // S 将各列缩放到单位二范数，求出的 u 再通过 dz=S*u 恢复到原参数尺度。
            Eigen::VectorXd scales(ncol);

            for (int j = 0; j < ncol; ++j)
            {
                double norm = A.col(j).norm();

                if (!(norm > 0))
                {
                    result.message = "Zero sensitivity column.";
                    tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                    return false;
                }

                scales[j] = 1 / norm;
            }

            // 列尺度仅改善求解条件，不改变物理问题；QR 同时给出可观测秩。
            Eigen::MatrixXd scaled = A * scales.asDiagonal();
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(scaled);
            qr.setThreshold(1e-11);

            // 缺秩说明现有资料不能区分全部待估分量，不能继续给出唯一改正。
            if (qr.rank() != ncol)
            {
                result.message = "Rank deficient observation geometry; estimated parameters are not observable.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            // 先解列尺度归一后的改正，再恢复到初始根数和 AMR 比例的原坐标。
            const Eigen::VectorXd scaled_correction = qr.solve(r);
            const Eigen::VectorXd correction = scales.asDiagonal() * scaled_correction;

            // raw 只记录六根数改正的最大绝对分量；before 是本轮接受资料的白化残差平方和。
            const double raw = correction.head<6>().cwiseAbs().maxCoeff();
            const double before = r.squaredNorm();
            const double raw_parameter = parameter ? std::abs(correction[6]) : 0;

            // 只有相邻两轮采用集合相同且上一轮 RMS 有效时，才计算可比较的 DRMS。
            const bool comparable = std::isfinite(previous_rms) && previous_accepted == accepted;
            double fraction = 1, after = before;
            bool improved = false, evaluated_candidate = false;

            // 先试完整改正，若越出轨道定义域或残差增大则减半；最多试十二次。
            for (int trial = 0; trial < 12; ++trial)
            {
                StateVector candidate = z + fraction * correction.head<6>();
                const double candidate_scale = parameter_scale + (parameter ? fraction * correction[6] : 0);

                // 候选根数须保持椭圆 NSG2 定义域，候选 AMR 比例也必须有限且严格为正。
                if (candidate[0] <= 0 || candidate.segment<2>(1).squaredNorm() >= 1 ||
                    candidate.segment<2>(3).squaredNorm() >= 1 || !candidate.allFinite() ||
                    !std::isfinite(candidate_scale) || candidate_scale <= 0)
                {
                    // 越出定义域时不传播该候选轨道，直接减半改正后继续尝试。
                    fraction *= 0.5;
                    continue;
                }

                // 候选轨道单独保存，当前线性化点的残差、A 和设计矩阵在搜索期间保持不变。
                Compute_OrbitResiduals(tracer, job, observations, epochs, trial_settings, candidate, candidate_scale,
                                       candidate_node_epochs_utc, candidate_node_nsg2_norm, candidate_node_cart_norm,
                                       candidate_states, candidate_residuals, candidate_weighted_residual);
                if (tracer.Code() != bastion::ERR_NOERROR)
                {
                    result.message = "Failed to evaluate a line-search candidate.";
                    return false;
                }

                evaluated_candidate = true;
                after = 0;

                // 目标值只累计本轮已接受的记录；线搜索期间不重新剔除，保证每个候选解比较同一问题。
                size_t objective_row = 0;

                for (size_t i = 0; i < accepted.size(); ++i)
                {
                    const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

                    if (accepted[i]) after += candidate_weighted_residual.segment(objective_row, dimension).squaredNorm();

                    // 剔除记录仍占完整残差向量中的行，故源行必须继续前移。
                    objective_row += dimension;
                }

                // 允许 1e-12 相对舍入余量；一旦目标值不增，状态与 AMR 比例同时接受。
                if (after <= before * (1 + 1e-12))
                {
                    z = candidate;
                    parameter_scale = candidate_scale;

                    // 初态、力参数、完整轨迹和残差同时接受，下一轮及最终输出直接复用。
                    node_epochs_utc.swap(candidate_node_epochs_utc);
                    node_nsg2_norm.swap(candidate_node_nsg2_norm);
                    node_cart_norm.swap(candidate_node_cart_norm);
                    states.swap(candidate_states);
                    residuals.swap(candidate_residuals);
                    weighted_residual.swap(candidate_weighted_residual);
                    improved = true;
                    break;
                }

                // 该有效候选仍增大残差，缩短步长后重新传播。
                fraction *= 0.5;
            }

            // 搜索没有接受步时，只能在残差判据已达标的情况下保留当前解。
            // 全部试探均越出定义域不能被当作 DRMS 停滞；至少要有一次有效残差计算。
            if (!improved)
            {
                const double current_rms = std::sqrt(before / A.rows());
                const double change = comparable ? ResidualRmsChange(previous_rms, current_rms, settings.drms_relative)
                                                 : std::numeric_limits<double>::quiet_NaN();

                // 只有确实完成有效试轨，且现有解已满足 RMS/DRMS，才允许以零步长保留。
                if (!(evaluated_candidate && (current_rms <= settings.rms_tolerance ||
                        (comparable && change <= settings.drms_tolerance))))
                {
                    result.message = "Line search stalled before RMS/DRMS convergence.";
                    tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                    return false;
                }

                // 保留更新前轨道与目标值，后面的迭代记录据此报告实际采用的零改正。
                fraction = 0;
                after = before;
            }

            // 以本轮改正后的标准化 RMS 判断收敛；改正量仅记录，不参与停止条件。
            IterationRecord record;
            record.iteration = iteration;
            record.rank = ncol;

            // 记录资料数量和完成更新后的 WRMS，分母是接受的分量数而非记录数。
            record.used_observations = std::count(accepted.begin(), accepted.end(), true);
            record.weighted_rms = std::sqrt(after / A.rows());

            // 原始改正与线搜索比例分别保存，便于区分求解器建议步长和实际采用步长。
            record.scaled_correction = raw;
            record.step_fraction = fraction;

            // AMR 值恢复到 m²/kg；parameter_correction 仍是未缩短的无量纲比例改正。
            record.parameter_amr_si = parameter ? prior_amr_si * parameter_scale : 0;
            record.parameter_correction = raw_parameter;

            // 首轮或资料集合改变时不计算 DRMS；其缺省 NaN 明确表示不可比较。
            record.drms_comparable = comparable;

            if (comparable) record.drms = ResidualRmsChange(previous_rms, record.weighted_rms, settings.drms_relative);

            // RMS 或可比较的 DRMS 任一通过即可，但还必须满足最少轮次与剔除启动轮次。
            record.rms_converged = record.weighted_rms <= settings.rms_tolerance;
            record.drms_converged = comparable && record.drms <= settings.drms_tolerance;
            const bool criteria_met = record.rms_converged || record.drms_converged;

            // 到达剔除启动轮次前继续迭代，确保已启用的粗差策略至少真正执行一次。
            const bool stop_allowed = criteria_met && iteration >= settings.min_iterations &&
                (settings.rejection_scale == 0 || iteration >= settings.rejection_start);
            result.iterations.push_back(record);

            {
                // 本轮采用标记同时用于总计和分类计数，百分比各以本组完整记录数为分母。
                std::array<std::size_t, 8> used_by_type{}, rejected_by_type{};

                for (std::size_t i = 0; i < observations.size(); ++i)
                {
                    const int type = static_cast<int>(observations[i].type);

                    if (accepted[i])
                    {
                        ++used_by_type[type];
                    }
                    else
                    {
                        ++rejected_by_type[type];
                    }
                }

                const int label_width = 40, state_width = 20;
                const std::size_t used = record.used_observations, rejected = observations.size() - used;
                std::ostringstream progress;
                progress << "\n---------------------------- BODN iteration #" << iteration
                         << " ----------------------------\n\n"
                         << std::left << std::setw(label_width) << "- Used / rejected records" << ": "
                         << used << " (" << std::fixed << std::setprecision(2) << 100.0 * used / observations.size()
                         << "%) / " << rejected << " (" << 100.0 * rejected / observations.size() << "%)\n";

                for (int type : {1, 6, 7})
                {
                    const std::size_t total = used_by_type[type] + rejected_by_type[type];

                    if (total == 0) continue;

                    const std::string label = "  > " + std::string(observation::ObservationTypeName(type)) + " (" + std::to_string(type) + ")";
                    progress << std::setw(label_width) << label << ": " << used_by_type[type] << " ("
                             << 100.0 * used_by_type[type] / total << "%) / " << rejected_by_type[type] << " ("
                             << 100.0 * rejected_by_type[type] / total << "%)\n";
                }

                progress << '\n' << std::setw(label_width) << "- LS rank / parameters" << ": " << record.rank << " / " << ncol
                         << '\n' << std::setw(label_width) << "- Correction scale" << ": "
                         << std::defaultfloat << std::setprecision(14) << record.step_fraction << "\n\n";

                // 先列实际采用的 NSG2 改正，再列更新后的状态；两者均采用内部归一化尺度。
                for (int field = 0; field < 2; ++field)
                {
                    for (int row = 0; row < 2; ++row)
                    {
                        if (row == 0)
                        {
                            progress << std::left << std::setw(label_width)
                                     << (field == 0 ? "- Correction of NSG2 [norm]" : "- Value of updated NSG2 [norm]") << ": ";
                        }
                        else
                        {
                            progress << std::string(label_width + 2, ' ');
                        }

                        progress << "[ ";

                        for (int j = 0; j < 3; ++j)
                        {
                            if (j > 0) progress << ", ";

                            progress << std::right << std::setw(state_width);

                            if (field == 0)
                            {
                                progress << std::showpos << std::scientific << std::setprecision(6)
                                         << fraction * correction[row * 3 + j] << std::noshowpos;
                            }
                            else
                            {
                                progress << std::defaultfloat << std::setprecision(14) << z[row * 3 + j];
                            }
                        }

                        progress << " ]\n";
                    }
                }

                // 附加 AMR 使用物理单位；只在确实参与估计时列出实际改正和更新值。
                if (parameter)
                {
                    const std::string name = parameter == 1 ? "drag" : "SRP";
                    progress << '\n' << std::left << std::setw(label_width + 1) << "- Correction of " + name + " AMR [m²/kg]"
                             << ": " << std::defaultfloat << std::setprecision(14) << fraction * prior_amr_si * correction[6]
                             << '\n' << std::setw(label_width + 1) << "- Value of updated " + name + " AMR [m²/kg]"
                             << ": " << record.parameter_amr_si << '\n';
                }

                // 指标只列数值；继续或停止及其原因集中在本轮最后说明。
                progress << '\n' << std::left << std::setw(label_width) << "- WRMS pre / post" << ": "
                         << std::defaultfloat << std::setprecision(14) << std::sqrt(before / A.rows()) << " / " << record.weighted_rms
                         << '\n' << std::setw(label_width) << "- DRMS" << ": ";

                if (comparable)
                {
                    progress << record.drms;
                }
                else
                {
                    progress << "N/A";
                }

                progress << '\n' << std::setw(label_width) << "- Iteration result" << ": ";

                if (stop_allowed)
                {
                    progress << "CONVERGED — " << (record.rms_converged ? "WRMS" : "DRMS");
                }
                else if (iteration == settings.max_iterations)
                {
                    progress << "NOT CONVERGED — ITERATION LIMIT";
                }
                else if (criteria_met && iteration < settings.min_iterations)
                {
                    progress << "CONTINUE — MINIMUM ITERATIONS";
                }
                else if (criteria_met)
                {
                    progress << "CONTINUE — REJECTION START ROUND";
                }
                else
                {
                    progress << "CONTINUE";
                }

                progress << '\n';
                result.iteration_report += progress.str();

                if (settings.display_iterations) std::cerr << progress.str();

                // 原始改正追加到任务报告，控制台只显示实际采用的改正。
                std::ostringstream raw_correction;

                for (int row = 0; row < 2; ++row)
                {
                    if (row == 0)
                    {
                        raw_correction << std::left << std::setw(label_width) << "- Raw correction of NSG2 [norm]" << ": ";
                    }
                    else
                    {
                        raw_correction << std::string(label_width + 2, ' ');
                    }

                    raw_correction << "[ ";

                    for (int j = 0; j < 3; ++j)
                    {
                        if (j > 0) raw_correction << ", ";

                        raw_correction << std::right << std::setw(state_width) << std::defaultfloat << std::setprecision(14) << correction[row * 3 + j];
                    }

                    raw_correction << " ]\n";
                }

                if (parameter)
                {
                    raw_correction << std::left << std::setw(label_width + 1) << "- Raw correction of AMR [m²/kg]" << ": "
                                   << std::defaultfloat << std::setprecision(14) << prior_amr_si * correction[6] << '\n';
                }

                result.iteration_report += raw_correction.str();
            }

            // 保存已完成更新的 RMS 与资料集合，下一轮只能与这组同源结果作 DRMS 比较。
            previous_rms = record.weighted_rms;
            previous_accepted = accepted;

            // 通过所有停止条件后记录具体判据并退出；否则在已接受轨迹处继续线性化。
            if (stop_allowed)
            {
                result.converged = true;
                result.convergence_reason = record.rms_converged ? "RMS" :
                    settings.drms_relative ? "DRMS_RELATIVE" : "DRMS_ABSOLUTE";
                break;
            }
        }

        // ==================== 04 保存最终轨迹并计算最终偏导 ====================

        // 接受候选解时已经保存其完整轨迹与残差；零改正时保留当前轨迹。
        // 最终偏导仍在接受的 z 和 AMR 处重新计算，采用集合保持末轮的筛选结果。
        Build_DesignMatrix(tracer, job, observations, epochs, trial_settings, z, parameter_scale, states, design_matrix);
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            result.message = "Failed to build the batch design matrix.";
            return false;
        }

        // 将接受的 z 和 AMR 同步到最终任务；结果初态与保留的完整轨迹来自同一估计。
        result.fit_epoch_utc = epochs.front();
        Set_InitialOrbit(tracer, job, z);
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            result.message = "Failed to synchronize the final orbit or AMR.";
            return false;
        }

        Set_TrialAmr(tracer, job, parameter, prior_amr_norm, parameter_scale);
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            result.message = "Failed to synchronize the final orbit or AMR.";
            return false;
        }

        // 物理参数输出使用 SI；再从最终 Cartesian 得到最终 NSG2，保持两种表示同步。
        result.parameter_amr_si = parameter ? prior_amr_si * parameter_scale : 0;
        job.initCart_norm.ToNonSingular2(tracer, result.fitted_nsg2_norm);
        if (tracer.Code() != bastion::ERR_NOERROR)
        {
            result.message = "Failed to convert the final orbit to NSG2.";
            return false;
        }

        result.fitted_cart_norm = job.initCart_norm;

        // 星历按逐条观测历元输出，原始 O-C 保留角度/米制单位，不使用白化残差替代。
        result.node_epochs_utc = std::move(node_epochs_utc);
        result.node_nsg2_norm = std::move(node_nsg2_norm);
        result.node_cart_norm = std::move(node_cart_norm);
        result.epochs_utc = epochs;
        result.ephemeris_norm = states;
        result.residuals = residuals;

        // 接受标记继续使用末轮线性化时的选择；下面的协方差也使用完全相同的资料。
        result.accepted = accepted;
        Eigen::MatrixXd A;
        Eigen::VectorXd r;

        // 用末轮保留的采用集合提取最终方程，不在最终输出阶段改变已求解的资料集合。
        size_t selected_rows = 0;

        for (size_t i = 0; i < accepted.size(); ++i)
            if (accepted[i]) selected_rows += (observations[i].type == ObservationType::InertialPosition ? 3 : 2);

        // 最终采用分量数应多于待估参数数，再分配最终 QR 所需的矩阵和向量。
        if (selected_rows <= static_cast<size_t>(ncol))
        {
            result.message = "Too few components remain after rejection.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        A.resize(selected_rows, ncol);
        r.resize(selected_rows);

        // 同步提取残差和偏导；原始行包括剔除记录，目标行仅包括接受记录。
        size_t from_row = 0, to_row = 0;

        for (size_t i = 0; i < accepted.size(); ++i)
        {
            const int dimension = observations[i].type == ObservationType::InertialPosition ? 3 : 2;

            if (accepted[i])
            {
                // 整条记录的所有分量共用同一采用标记，不能只复制其中残差较小的分量。
                A.middleRows(to_row, dimension) = design_matrix.middleRows(from_row, dimension);
                r.segment(to_row, dimension) = weighted_residual.segment(from_row, dimension);
                to_row += dimension;
            }

            from_row += dimension;
        }

        // 最终 WRMS 只来自被采用的白化分量；这个统计量不用于额外放大形式协方差。
        result.weighted_rms = sqrt(r.squaredNorm() / r.size());

        // ==================== 05 检查最终秩并计算形式协方差 ====================

        // 最终轨道重新检查可观测秩，并准备与最终设计矩阵一致的协方差。
        Eigen::VectorXd scales(ncol);

        for (int j = 0; j < ncol; ++j)
        {
            if (!(A.col(j).norm() > 0))
            {
                result.message = "Final zero sensitivity column.";
                tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
                return false;
            }

            scales[j] = 1 / A.col(j).norm();
        }

        // 与逐轮求解一样先按列缩放，再由同一最终 A 检查秩并计算形式协方差。
        const Eigen::MatrixXd scaled = A * scales.asDiagonal();
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(scaled);
        qr.setThreshold(1e-11);

        // 检查最终设计矩阵的列秩，并用 NaN 初始化尚未计算的联合协方差。
        if (qr.rank() != ncol)
        {
            result.message = "Final solution has rank deficient geometry.";
            tracer.LogEvent(__FILE__, __func__, __LINE__, result.message, bastion::EventType::ERROR);
            return false;
        }

        result.covariance_full = Eigen::MatrixXd::Constant(ncol, ncol, std::numeric_limits<double>::quiet_NaN());

        if (settings.compute_covariance)
        {
            // A 已白化，已知观测方差下的形式协方差为 (A^T A)^(-1)，不按残差 WRMS² 放大。
            // 列主元 QR 满足 (A*S)*P=Q*R，先构造 P*R^(-1)*R^(-T)*P^T，再恢复 S 的尺度。
            const Eigen::MatrixXd R = qr.matrixR().topLeftCorner(ncol, ncol).triangularView<Eigen::Upper>();
            const Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(ncol, ncol);

            // 通过上三角回代计算 R 的逆，再按 QR 的列主元置换还原参数顺序。
            const Eigen::MatrixXd inverseR = R.triangularView<Eigen::Upper>().solve(identity);
            const Eigen::MatrixXd permuted_inverseR = qr.colsPermutation() * inverseR;
            const Eigen::MatrixXd permuted_covariance = permuted_inverseR * inverseR.transpose();

            // 最后右乘 P^T，得到列缩放坐标中的协方差；还未恢复根数与 AMR 的实际尺度。
            const Eigen::MatrixXd scaled_covariance = permuted_covariance * qr.colsPermutation().transpose();

            // 末维从无量纲 AMR 比值恢复到 m²/kg；轨道块取联合协方差的左上角，已包含 AMR 同估的影响。
            if (parameter) scales[6] *= prior_amr_si;

            const Eigen::MatrixXd left_scaled_covariance = scales.asDiagonal() * scaled_covariance;
            result.covariance_full = left_scaled_covariance * scales.asDiagonal();

            // 联合协方差左上六阶块包含与 AMR 同估的相关性，不能重新做仅六状态的独立逆。
            result.covariance_nsg2 = result.covariance_full.topLeftCorner<6, 6>();

            if (parameter) result.parameter_sigma_amr_si = std::sqrt(result.covariance_full(6, 6));

            // P_x = G P_z G^T；G 的位置行乘 ae、速度行乘 ae/UNITTIME，恢复到 SI。
            StateMatrix G = GeometryPartialsNSG2(tracer, result.fitted_nsg2_norm);
            if (tracer.Code() != bastion::ERR_NOERROR)
            {
                result.message = "Failed to transform the final covariance.";
                return false;
            }

            G.topRows<3>() *= ae;

            // 后三行对应速度，按 ae/UNITTIME 恢复到 m/s 的尺度。
            G.bottomRows<3>() *= ae / UNITTIME;

            // 先计算 G*P_z，再右乘 G^T，得到位置和速度的 SI 协方差。
            const StateMatrix state_element_covariance = G * result.covariance_nsg2;
            result.covariance_cart_si = state_element_covariance * G.transpose();
        }

        // 迭代上限仍保留最近解、星历与协方差，但返回 false 并明确标为未收敛。
        if (!result.converged) result.convergence_reason = "ITERATION_LIMIT";

        result.message = result.converged ? "Converged (" + result.convergence_reason + ")"
                                         : "Maximum iterations reached; not converged (ITERATION_LIMIT)";

        // 正常路径结束本次 tracer 调用层；返回值只表达是否真正满足停止条件。
        tracer.CheckOut(__func__);
        return result.converged;
    }
    catch (const std::exception& error)
    {
        // 异常文本写入结果并进入 tracer；调用方可继续处理失败，不会在库内退出整个程序。
        result.message = error.what();
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), bastion::EventType::ERROR);
        return false;
    }
}

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
void BatchSetting::Display(bastion::Overwatch& tracer, double sigma_angle_rad, double sigma_position_m, std::ostream& output, int precision) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 设置按待估量、观测精度、偏导、剔除及停止条件排列；模式名称对应当前配置。
    const int label_width = 40;
    std::ostringstream text;
    text << std::left << std::defaultfloat << std::setprecision(precision)
         << "\n================================= Estimation settings ==================================\n\n"
         << std::setw(label_width) << "- Estimated parameters" << ": Six NSG2 elements at first observation";

    if (estimate_parameter == 1) text << " + drag AMR";

    if (estimate_parameter == 2) text << " + SRP AMR";

    if (estimate_parameter == 3) text << "; automatic AMR selection";

    text << '\n' << std::setw(label_width) << "- Observation error standard deviation" << ": "
         << sigma_angle_rad * 180 * 3600 / pi << " arcsec (East / North); " << sigma_position_m << " m (XYZ)"
         << '\n' << std::setw(label_width) << "- STM computation method" << ": ";

    if (stm_method == 2)
    {
        text << "Full-trajectory central differences; step " << std::scientific << std::setprecision(2) << fd_step;
    }
    else
    {
        text << "Analytical approximation; J2 " << (analytic_stm.use_J2 ? "ON" : "OFF")
             << "; lunisolar " << (analytic_stm.use_lunisolar ? "ON" : "OFF");
    }

    // 剔除关闭时不展开倍数和启动轮次；显式收紧的硬窗口单独列出。
    text << "\n\n" << std::setw(label_width) << "- Outlier rejection rule" << ": " << std::defaultfloat << std::setprecision(precision);

    if (rejection_scale > 0)
    {
        text << rejection_scale << " × RMS; from round " << rejection_start;

        if (rejection_angle_window_rad != 1.0)
        {
            text << '\n' << std::setw(label_width) << "  > Angular hard window [rad]" << ": " << rejection_angle_window_rad;
        }

        if (rejection_position_window_m != ae)
        {
            text << '\n' << std::setw(label_width) << "  > Position hard window [m]" << ": " << rejection_position_window_m;
        }
    }
    else
    {
        text << "OFF";
    }

    text << '\n' << std::setw(label_width) << "- Min / Max of Iterations" << ": " << min_iterations << " / " << max_iterations
         << '\n' << std::setw(label_width) << "- WRMS / DRMS limits" << ": " << std::scientific << std::setprecision(2)
         << rms_tolerance << " / " << drms_tolerance;

    if (drms_relative) text << '\n' << std::setw(label_width) << "- DRMS mode" << ": Relative";

    text << '\n' << std::setw(label_width) << "- Formal covariance" << ": " << (compute_covariance ? "ON" : "OFF") << '\n';
    output << text.str();

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing estimation settings.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    tracer.CheckOut(__func__);
}

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
void BatchResult::PrintSummary(bastion::Overwatch& tracer, const BatchSetting& setting, const orbit::Keplerian& kep, std::ostream& output, int precision) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (iterations.empty())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Batch output requires at least one completed iteration.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    // 后验状态在同一参考历元准备；只对显示数组恢复单位，不修改结果中的归一化状态。
    double cart_values[6], kep_values[6], nsg2_values[6];
    double cart_sigma[6], nsg2_sigma[6];
    fitted_cart_norm.ToArray(cart_values);
    fitted_nsg2_norm.ToArray(nsg2_values);
    kep.ToArray(kep_values);

    for (int j = 0; j < 6; ++j)
    {
        cart_values[j] *= ae / 1000.0 / (j < 3 ? 1.0 : UNITTIME);
        const double cart_variance = covariance_cart_si(j, j);
        const double nsg2_variance = covariance_nsg2(j, j);
        cart_sigma[j] = std::isfinite(cart_variance) && cart_variance >= 0 ? std::sqrt(cart_variance) / 1000.0
                                                                        : std::numeric_limits<double>::quiet_NaN();
        nsg2_sigma[j] = std::isfinite(nsg2_variance) && nsg2_variance >= 0 ? std::sqrt(nsg2_variance)
                                                                        : std::numeric_limits<double>::quiet_NaN();
    }

    // 收敛原因使用控制指标的统一名称；内部结果标识仍保存绝对或相对 DRMS 的区别。
    const int label_width = 40, state_width = 20;
    std::ostringstream summary;
    summary << "\n===================================== BODN summary =====================================\n\n"
            << std::left << std::setw(label_width) << "- Status" << ": ";

    if (converged)
    {
        summary << "CONVERGED — " << (convergence_reason == "RMS" ? "WRMS" : "DRMS");
    }
    else
    {
        summary << "NOT CONVERGED — ITERATION LIMIT";
    }

    summary << '\n' << std::setw(label_width) << "- Iterations" << ": " << iterations.size()
            << '\n' << std::setw(label_width) << "- DRMS" << ": " << std::defaultfloat << std::setprecision(precision);

    if (iterations.back().drms_comparable)
    {
        summary << iterations.back().drms;
    }
    else
    {
        summary << "N/A";
    }

    summary << "\n\n" << std::setw(label_width) << "- Reference epoch [UTC]" << ": "
            << std::fixed << std::setprecision(12) << fit_epoch_utc << '\n';

    // 三种表示按相同的两行数组排列；Cartesian、NSG2 的形式标准差随相应状态展开。
    const std::array<const char*, 3> labels = {"- Posterior Cartesian [km, km/s]", "- Posterior Keplerian [km, deg]", "- Posterior NSG2 [norm]"};
    const std::array<const double*, 3> values = {cart_values, kep_values, nsg2_values};

    for (int representation = 0; representation < 3; ++representation)
    {
        if (representation > 0) summary << '\n';

        for (int row = 0; row < 2; ++row)
        {
            if (row == 0)
            {
                summary << std::left << std::setw(label_width) << labels[representation] << ": ";
            }
            else
            {
                summary << std::string(label_width + 2, ' ');
            }

            summary << "[ ";

            for (int j = 0; j < 3; ++j)
            {
                if (j > 0) summary << ", ";

                summary << std::right << std::setw(state_width) << std::defaultfloat << std::setprecision(precision)
                        << values[representation][row * 3 + j];
            }

            summary << " ]\n";
        }

        // Keplerian 仅显示已有根数，不借用其他表示的协方差。
        if (representation == 1) continue;

        summary << std::left << std::setw(label_width)
                << (representation == 0 ? "  > Formal standard deviation [km, km/s]" : "  > Formal standard deviation [norm]") << ": ";

        if (!setting.compute_covariance)
        {
            summary << "Not computed\n";
            continue;
        }

        const double* sigma = representation == 0 ? cart_sigma : nsg2_sigma;

        for (int row = 0; row < 2; ++row)
        {
            if (row > 0) summary << std::string(label_width + 2, ' ');

            summary << "[ ";

            for (int j = 0; j < 3; ++j)
            {
                if (j > 0) summary << ", ";

                const double value = sigma[row * 3 + j];
                summary << std::right << std::setw(state_width);

                if (std::isfinite(value))
                {
                    summary << std::scientific << std::setprecision(6) << value;
                }
                else
                {
                    summary << "N/A";
                }
            }

            summary << " ]\n";
        }
    }

    // 同估 AMR 时，附加后验、先验与改正、形式标准差，单位均为 m²/kg。
    if (estimated_parameter)
    {
        summary << '\n' << std::left << std::setw(label_width + 1)
                << (estimated_parameter == 1 ? "- Posterior drag AMR [m²/kg]" : "- Posterior SRP AMR [m²/kg]")
                << ": " << std::defaultfloat << std::setprecision(precision) << parameter_amr_si
                << '\n' << std::setw(label_width) << "  > Prior / correction" << ": " << parameter_apriori_amr_si
                << " / " << parameter_amr_si - parameter_apriori_amr_si
                << '\n' << std::setw(label_width + 1) << "  > Formal standard deviation [m²/kg]" << ": ";

        if (!setting.compute_covariance)
        {
            summary << "Not computed";
        }
        else if (std::isfinite(parameter_sigma_amr_si))
        {
            summary << std::scientific << std::setprecision(6) << parameter_sigma_amr_si;
        }
        else
        {
            summary << "N/A";
        }

        summary << '\n';
    }

    output << summary.str();

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing BODN summary.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    tracer.CheckOut(__func__);
}

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
void BatchResult::PrintReport(bastion::Overwatch& tracer, const BatchSetting& setting, const std::string& object_id,
    const std::string& initial_source, const orbit::Keplerian& kep, std::ostream& output, int precision) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (iterations.empty())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Batch output requires at least one completed iteration.", bastion::ERROR);
        return;
    }

    // 报告独立准备文字格式；状态数组直接从已有对象导出。
    std::ostringstream report;
    const std::string report_rule(100, '=');
    const std::string table_rule(100, '-');
    const auto accepted_count = std::count(accepted.begin(), accepted.end(), true);
    double final[6];
    fitted_cart_norm.ToArray(final);

    // 报告先写任务结果和停止条件，再说明目标、参考历元及残差统计约定。
    report << report_rule << "\n"
           << "                         BODN ORBIT DETERMINATION REPORT\n"
           << report_rule << "\n\n01  RESULT SUMMARY\n" << table_rule << "\n"
           << std::setprecision(17) << "Status: " << message
           << "\nStopping rule: standardized RMS OR DRMS; minimum rounds and rejection-start gate apply."
           << "\nRMS limit [1]: " << setting.rms_tolerance
           << "\nDRMS limit: " << setting.drms_tolerance << (setting.drms_relative ? " (relative)" : " (absolute)")

           // DRMS 只有在相邻两轮采用集合相同时才可比较，报告同时保留这一标记。
           << "\nFinal DRMS: " << iterations.back().drms
           << "; comparable: " << iterations.back().drms_comparable
           << "\nObject: " << object_id
           << "\nReference epoch MJD UTC: " << fit_epoch_utc

           // 记录初轨来源和状态约定，并列出迭代次数、采用记录数及最终 WRMS。
           << "\nInitial orbit source: " << initial_source
           << "\nFrame: GCRS; propagation time TT\nEstimated: six NSG2 elements"
           << (estimated_parameter == 1 ? " + drag AMR" : estimated_parameter == 2 ? " + SRP AMR" : "")
           << "\nIterations: " << iterations.size()
           << "\nUsed observations: " << accepted_count << " / " << accepted.size()
           << "\nWeighted RMS (dimensionless): " << std::scientific << std::setprecision(12) << weighted_rms

           // 解释 WRMS、东西方向角残差及形式协方差的定义。
           << "\nWRMS definition: sqrt(sum(((O-C)/sigma)^2) / accepted component count)."
           << "\nRA residual: wrap(O-C)*cos(observed declination); geometric measurements only."
           << "\nCovariance: "
           << (setting.compute_covariance ? "formal, prescribed independent sigmas; no RMS rescaling"
                                          : "not computed")

           // 以下状态表使用先验参考历元的后验轨道，Cartesian 以 m、m/s 显示。
           << "\n\n" << report_rule << "\n02  POSTERIOR STATE\n" << table_rule
           << "\n[Posterior state at reference epoch]\n"
           << "Epoch (MJD UTC): " << std::fixed << std::setprecision(12) << fit_epoch_utc
           << "\nFrame: GCRS; osculating orbit at the apriori reference epoch.\n\n"
           << "Cartesian state\n"
           << "  Component   Unit                       Value\n";

    // 六个 Cartesian 分量逐行列出，状态表使用 SI 单位，与形式协方差保持一致。
    const std::array<const char*, 6> cart_names{"x", "y", "z", "vx", "vy", "vz"};

    // 按六分量顺序写入后验 Cartesian，位置恢复为 m，速度恢复为 m/s。
    for (int j = 0; j < 6; ++j)
        report << "  " << std::left << std::setw(12) << cart_names[j]
               << std::setw(10) << (j < 3 ? "m" : "m/s") << std::right
               << std::scientific << std::setprecision(precision) << std::setw(26)
               << final[j] * ae / (j < 3 ? 1 : UNITTIME) << '\n';

    // 密切根数采用 km 与 deg；M 是平近点角，不能与真近点角混用。
    report << "\nKeplerian osculating elements (M = mean anomaly)\n"
           << "  Component   Unit                       Value\n";
    const std::array<const char*, 6> kep_names{"a", "e", "i", "OMEGA", "omega", "M"};

    // 逐项写入后验 Kepler 密切根数，最后一项为平近点角 M。
    double final_kep_values[6];
    kep.ToArray(final_kep_values);

    for (int j = 0; j < 6; ++j)
        report << "  " << std::left << std::setw(12) << kep_names[j]
               << std::setw(10) << (j == 0 ? "km" : j == 1 ? "1" : "deg") << std::right
               << std::scientific << std::setprecision(precision) << std::setw(26)
               << final_kep_values[j] << '\n';

    // NSG2 直接展示最终估计量，半长轴和平均经度恢复为 km、deg；其余定义量无量纲。
    report << "\nNSG2 osculating elements\n"
           << "  Component   Unit                       Value\n";
    const std::array<const char*, 6> nsg_names{"a", "xi", "eta", "h", "k", "lambda"};
    std::array<double, 6> final_nsg2;
    fitted_nsg2_norm.ToArray(final_nsg2.data());

    // 逐项恢复 NSG2 的半长轴与经度单位，其余四项保持无量纲。
    for (int j = 0; j < 6; ++j)
        report << "  " << std::left << std::setw(12) << nsg_names[j]
               << std::setw(10) << (j == 0 ? "km" : j == 5 ? "deg" : "1") << std::right
               << std::scientific << std::setprecision(precision) << std::setw(26)
               << final_nsg2[j] * (j == 0 ? ae / 1000 : j == 5 ? 180 / pi : 1) << '\n';

    output << report.str();

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing display output.", bastion::ERROR);
        return;
    }

    tracer.CheckOut(__func__);
}
} // namespace determination
