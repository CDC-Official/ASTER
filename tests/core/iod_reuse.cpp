/**
 * @file iod_reuse.cpp
 * @brief 核对任务级日月距离复用、完整数值轨迹和按需报告。
 *
 * 在五种 IOD 动力学中比较现场计算与只读距离数组；使用相同几何资料与初值，
 * 同时覆盖简单迭代、微分改正、负时间及重复历元。失败路径检查现有 tracer。
 */
#include "core_iod_solver.hpp"
#include <cstring>
#include <limits>
#include <sstream>

using namespace Eigen;
using namespace iod;
using namespace bastion;

int main()
{
    Overwatch tracer(__FILE__, __func__, __LINE__);
    std::ostringstream program_log;
    int comparisons = 0;

    // 中间参考历元使相对时间同时含负值、零和正值；重复历元仍逐条对应观测。
    const double epoch = 58406.887;
    ArrayXd tau(8);
    tau << -0.18, -0.12, -0.06, 0.0, 0.0, 0.06, 0.12, 0.18;
    const VectorXd truth = (VectorXd(6) << 1.2, 0.1, 0.4, -0.1, 0.7, 0.25).finished();

    for (int model = IOD_DYNAMICS_TWOBODY; model <= IOD_DYNAMICS_J2_SOLAR_LUNAR; ++model)
    {
        // 每种模型独立准备第三体量，随后与不传复用数组的传播路径逐位比较。
        std::vector<IodThirdBodyDistances> distances;
        Prepare_ThirdBodyDistances(tracer, model, epoch, tau, distances);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        for (int k = 0; k < tau.size(); ++k)
        {
            for (double displacement : {0.0, 0.15})
            {
                double state[6], direct[4], reused[4];
                Eigen::Map<VectorXd>(state, 6) = truth;
                state[0] += displacement;
                Lagrange_Coefficients(model, epoch, tau(k), state, direct);
                Lagrange_Coefficients(model, epoch, tau(k), state, reused, &distances[k]);

                if (std::memcmp(direct, reused, sizeof(direct)) != 0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Reused Lagrange coefficients changed, model=" + std::to_string(model), EventType::ERROR);
                    tracer.EmergencyExit(program_log);
                }

                ++comparisons;
            }
        }

        // 几何资料只用于实现等价检查；两条求解路径共用同一观测和同一扰动初值。
        IodJob base{};
        base.dynamicsModel = model;
        base.epochMJDUTC = epoch;
        base.iterCountMax = 8;
        base.iterTOL = 1e-10;
        base.constraintType = IOD_CONSTRAINT_ORIENTATION;
        base.printIterInfo = false;
        base.obs.type = OBS_TYPE_RADEC;
        base.obs.number = tau.size();
        base.obs.tau = tau;
        base.obs.time_MJDUTC = epoch + tau * constants::Conversion::NORMALIZED_TIME_TO_DAY;
        base.obs.ra_rad.resize(tau.size());
        base.obs.dec_rad.resize(tau.size());
        base.observer.positions_gcrs_norm_emat.resize(3, tau.size());
        base.observer.position_gcrs_epoch_norm = Vector3d(0.6, -0.5, 0.6);
        base.initValue.type = IOD_INITIAL_VALUE_SPECIFIED;
        base.initValue.initCart.Set(truth(0) + 1e-4, truth(1), truth(2), truth(3), truth(4), truth(5), true);

        for (int k = 0; k < tau.size(); ++k)
        {
            // 固定站址与真实状态的传播仅生成共同输入，不作为物理精度参考。
            base.observer.positions_gcrs_norm_emat.col(k) = base.observer.position_gcrs_epoch_norm;
            double state[6], fg[4];
            Eigen::Map<VectorXd>(state, 6) = truth;
            Lagrange_Coefficients(model, epoch, tau(k), state, fg);
            const Vector3d relative(fg[0] * truth(0) + fg[1] * truth(3) - 0.6,
                                    fg[0] * truth(1) + fg[1] * truth(4) + 0.5,
                                    fg[2] * truth(2) + fg[3] * truth(5) - 0.6);
            base.obs.ra_rad(k) = std::atan2(relative(1), relative(0));
            base.obs.dec_rad(k) = std::atan2(relative(2), std::hypot(relative(0), relative(1)));
        }

        for (int method : {IOD_ITERATION_SIMPLE, IOD_ITERATION_NEWTON})
        {
            IodJob direct = base, reused = base;
            direct.iterMethod = reused.iterMethod = method;
            reused.saveIterationReport = false;
            std::ostringstream streamed_report;
            iod_solver::InitialOrbitDetermination(tracer, direct);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            iod_solver::InitialOrbitDetermination(tracer, reused, &streamed_report, &distances);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 比较每一条历史的全部数值位，包括起始 NaN；延后格式化应重建同样的文字。
            std::string reconstructed;
            bool equal = direct.solution.flag == reused.solution.flag && reused.solution.iteration_report.empty() &&
                         direct.solution.iterations.size() == reused.solution.iterations.size();

            for (std::size_t k = 0; equal && k < direct.solution.iterations.size(); ++k)
            {
                const IodIterationRecord& a = direct.solution.iterations[k];
                const IodIterationRecord& b = reused.solution.iterations[k];
                equal = a.iteration == b.iteration &&
                        std::memcmp(a.r0v0_norm.data(), b.r0v0_norm.data(), 6 * sizeof(double)) == 0 &&
                        std::memcmp(a.dr0v0_norm.data(), b.dr0v0_norm.data(), 6 * sizeof(double)) == 0;
                const double values_a[] = {a.rms_ra_cosdec_rad, a.rms_dec_rad, a.rms_angle_rad, a.max_dr0v0_norm};
                const double values_b[] = {b.rms_ra_cosdec_rad, b.rms_dec_rad, b.rms_angle_rad, b.max_dr0v0_norm};
                equal = equal && std::memcmp(values_a, values_b, sizeof(values_a)) == 0;
                reconstructed += iod_solver::Format_IterationRecord(b, reused.iterTOL, reused.iterCountMax);
            }

            if (!equal || reconstructed != direct.solution.iteration_report || streamed_report.str() != reconstructed)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "IOD history or report changed after reuse/deferred formatting.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 无即时输出流时也应完整保存历史，并跳过文字生成。
            IodJob deferred = base;
            deferred.iterMethod = method;
            deferred.saveIterationReport = false;
            iod_solver::InitialOrbitDetermination(tracer, deferred, nullptr, &distances);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            if (!deferred.solution.iteration_report.empty() || deferred.solution.iterations.size() != direct.solution.iterations.size())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Deferred report lost numeric history or generated unwanted text.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }
        }
    }

    // 错误输入必须进入 tracer，且不访问长度错误的复用数组。
    Overwatch invalid_model, invalid_times, invalid_count;
    std::vector<IodThirdBodyDistances> distances;
    Prepare_ThirdBodyDistances(invalid_model, 6, epoch, tau, distances);
    ArrayXd bad_tau = tau;
    bad_tau(0) = std::numeric_limits<double>::quiet_NaN();
    Prepare_ThirdBodyDistances(invalid_times, IOD_DYNAMICS_J2_SOLAR_LUNAR, epoch, bad_tau, distances);
    IodJob invalid_job{};
    invalid_job.obs.number = 8;
    iod_solver::InitialOrbitDetermination(invalid_count, invalid_job, nullptr, &distances);

    if (invalid_model.Code() == ERR_NOERROR || invalid_times.Code() == ERR_NOERROR || invalid_count.Code() == ERR_NOERROR)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid reuse input was accepted.", EventType::ERROR);
        tracer.EmergencyExit(program_log);
    }

    std::cout << "Exact coefficient comparisons: " << comparisons << "; 5 models x 2 iteration methods, history/report and invalid-input checks passed.\n";
    return 0;
}
