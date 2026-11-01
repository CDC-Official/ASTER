#include <array>
// 附加力参数验证：Cartesian 独立步进生成资料，NSG2 批处理同时恢复轨道和面质比。
#include "core_bodn.hpp"
#include <Eigen/Cholesky>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace determination;
namespace {
void require(bool pass, const std::string& message) {
    if (!pass) throw std::runtime_error(message);
}
OrbitPropagationJob make_job(bastion::Overwatch& t, bool drag) {
    OrbitPropagationJob j{};
    std::string eop = "data/eop/EOP-All_250516.txt";
    j.opn_EOP.Read_EOP(t, eop);
    j.initKep_norm.Set(drag ? 6800 : 42164, drag ? 0.005 : 0.03, drag ? 55 : 15, 20, 40, 10, false);
    j.initKep_norm.Normalize(t, ae / 1000);
    j.initKep_norm.SyncAnomaly(t);
    j.initKep_norm.ToCartesian(t, j.initCart_norm);
    j.opn_setting.t0_MJD_UTC = 59215;
    j.opn_setting.t0_MJD_TT = UTC2TT(59215);
    j.opn_setting.h0_norm = (drag ? 15 : 120) / UNITTIME;
    j.opn_setting.rkf.maxTE_norm = 1e-13;
    j.opn_setting.display_progress = 0;
    j.opn_model.use_earth_NSG = true;
    j.opn_model.earth_zonal_order = 2;
    j.opn_model.earth_tesseral_degree = 2;
    j.opn_model.use_drag = drag;
    j.opn_model.use_SRP = !drag;
    auto& s = j.opn_model.sat;
    s.mass = 100 / MASS_OF_EARTH;
    s.Cd = 2.2;
    s.eta_SRP = 0.3;
    s.AMR_drag = 0.02 * MASS_OF_EARTH / (ae * ae);
    s.AMR_SRP = 0.05 * MASS_OF_EARTH / (ae * ae);
    s.area_drag = s.AMR_drag * s.mass;
    s.area_SRP = s.AMR_SRP * s.mass;
    if (drag) {
        double spacing = 0;
        require(Read_AtmosphereDensity("data/atmosphere/atmosphere.txt", j.atm_list, spacing) > 0,
                "Missing density fixture");
        Set_ReferenceDensity(j.atm_list, j.initKep_norm, j.refDensity);
    }
    return j;
}
std::vector<Observation> truth(bastion::Overwatch& t, OrbitPropagationJob j, bool drag) {
    // 每一步直接到达资料时刻，未调用 NSG2 传播与 Hermite 插值；同一力模型下检查接口闭合。
    j.opn_setting.rkf.fixed_step = true;
    std::array<double, 6> exported_values_1; j.initKep_norm.ToArray(exported_values_1.data());
    const double period = 2 * pi * std::pow(exported_values_1[0], 1.5) * UNITTIME;
    std::array<double, 6> exported_values_2; j.initCart_norm.ToArray(exported_values_2.data());
    auto x = exported_values_2;
    double elapsed = 0;
    std::vector<Observation> obs;
    for (int i = 0; i <= 96; ++i) {
        Observation o;
        o.mjd_utc = 59215 + (3.2 * period * i / 96 + (i % 3) * 0.31) / 86400;
        o.type = ObservationType::InertialPosition;
        o.object_id = drag ? "DRAG" : "SRP";
        o.observer_id = "XYZ";
        o.sigma.setConstant(1);
        double target = (UTC2TT(o.mjd_utc) - UTC2TT(59215)) * 86400 / UNITTIME;
        while (target - elapsed > 1e-13) {
            double accepted = 0, step = std::min((drag ? 2 : 20) / UNITTIME, target - elapsed);
            RKF78_Cart(t, j, UTC2TT(59215) + elapsed * UNITTIME / 86400, x.data(), step, accepted);
            require(accepted > 0 && t.Code() == bastion::ERR_NOERROR, "Cartesian truth propagation failed");
            elapsed += accepted;
        }
        for (int k = 0; k < 3; ++k) o.value[k] = x[k] * ae;
        obs.push_back(o);
    }
    return obs;
}
OrbitPropagationJob perturbed(bastion::Overwatch& t, OrbitPropagationJob j, bool drag) {
    orbit::NonSingular2 n;
    j.initKep_norm.ToNonSingular2(t, n);
    std::array<double, 6> exported_values_3; n.ToArray(exported_values_3.data());
    auto z = exported_values_3;
    z[0] += 1e-5; z[1] += 2e-6; z[2] -= 3e-6;
    z[3] += 4e-6; z[4] -= 2e-6; z[5] += 1e-5;
    orbit::NonSingular2 p(t, z[0], z[1], z[2], z[3], z[4], z[5], true);
    p.SyncAnomaly(t); p.ToKeplerian(t, j.initKep_norm); p.ToCartesian(t, j.initCart_norm);
    auto& s = j.opn_model.sat;
    if (drag) { s.AMR_drag *= 0.7; s.area_drag *= 0.7; }
    else { s.AMR_SRP *= 0.7; s.area_SRP *= 0.7; }
    return j;
}
}
int main(int argc, char** argv) {
    try {
        bastion::Overwatch t;
        require(Read_JPLEPH("data/ephemeris/JPLEPH_C.430") > 0, "Missing JPL fixture");
        for (bool drag : {true, false}) {
            auto job = make_job(t, drag);
            auto obs = truth(t, job, drag);
            if (argc >= 2) {
                std::filesystem::create_directories(argv[1]);
                auto prefix = std::filesystem::path(argv[1]) / (drag ? "drag" : "srp");
                std::ofstream data(prefix.string() + ".obs");
                data << std::setprecision(17);
                for (const auto& o : obs)
                    data << o.mjd_utc << " 7 " << o.object_id << " XYZ " << o.value.transpose() << '\n';
                // 将仿真初轨直接写成 BODN 使用的 SimpleInput 文件，单位为 km、km/s、deg。
                orbit::Cartesian output_cart = job.initCart_norm;
                orbit::Keplerian output_kep;
                output_cart.ToKeplerian(t, output_kep);
                require(t.Code() == bastion::ERR_NOERROR, "Cannot convert the orbit fixture.");

                output_cart.Restore(t, ae / 1000, UNITTIME);
                output_kep.Restore(t, ae / 1000);
                require(t.Code() == bastion::ERR_NOERROR, "Cannot restore orbit fixture units.");

                double cart_values[6], kep_values[6];
                output_cart.ToArray(cart_values);
                output_kep.ToArray(kep_values);

                std::ofstream orbit_file(prefix.string() + ".si");
                orbit_file << std::setprecision(17) << "object_id = " << obs[0].object_id
                           << "\nstate_units = 2\nreference_epoch = 59215\nreference_Keplerian =";

                for (double value : kep_values) orbit_file << ' ' << value;

                orbit_file << "\nreference_Cartesian =";

                for (double value : cart_values) orbit_file << ' ' << value;

                orbit_file << '\n';
                orbit_file.close();
                require(bool(orbit_file), "Cannot write the orbit fixture.");
            }
            if (argc == 3 && std::string(argv[2]) == "--fixtures-only") continue;
            BatchResult reference;
            for (int variant = 0; variant < 3; ++variant) {
                const double fd = variant == 1 ? 3e-4 : 1e-3;
                BatchSetting s;
                s.estimate_parameter = 3;
                s.parameter_fd_step = fd;
                s.rms_tolerance = 1e-6;
                s.display_iterations = false;
                auto fitjob = perturbed(t, job, drag);
                if (variant == 2) {
                    // 固定半步积分另作一遍，检查估计值不依赖自适应节点与初始步长。
                    s.integrator_id = 102;
                    fitjob.opn_setting.h0_norm *= 0.5;
                }
                BatchResult r;
                bool ok = BatchOrbitDetermination_NSG2(t, fitjob, obs, s, r);
                std::cout << std::setprecision(12) << (drag ? "drag" : "srp") << " integrator " << s.integrator_id << " fd " << fd
                          << " converged " << ok << " AMR " << r.parameter_amr_si << " sigma " << r.parameter_sigma_amr_si
                          << " WRMS " << r.weighted_rms << " message " << r.message << '\n';
                require(ok, "Additional parameter fit failed: " + r.message);
                require(r.estimated_parameter == (drag ? 1 : 2), "Teacher auto chose the wrong parameter");
                const double expected = drag ? 0.02 : 0.05;
                require(std::abs(r.parameter_amr_si / expected - 1) < 2e-4, "AMR recovery inaccurate");
                require(r.weighted_rms < 0.01, "Parameter fit residual too large");
                std::array<double, 6> exported_values_4; job.initCart_norm.ToArray(exported_values_4.data());
                std::array<double, 6> exported_values_5; r.fitted_cart_norm.ToArray(exported_values_5.data());
                auto x = exported_values_5, target = exported_values_4;
                double distance = 0;
                for (int k = 0; k < 3; ++k) distance += std::pow((x[k] - target[k]) * ae, 2);
                std::cout << "position_error_m " << std::sqrt(distance) << " rounds " << r.iterations.size() << '\n';
                require(std::sqrt(distance) < 0.02, "Joint orbit recovery inaccurate");
                require(r.covariance_full.rows() == 7 && r.covariance_full.allFinite(), "Missing joint covariance");
                Eigen::VectorXd d = r.covariance_full.diagonal().cwiseSqrt().cwiseInverse();
                Eigen::MatrixXd correlation = d.asDiagonal() * r.covariance_full * d.asDiagonal();
                require(Eigen::LLT<Eigen::MatrixXd>(correlation).info() == Eigen::Success, "Joint covariance is not positive definite");
                require((r.covariance_full.topLeftCorner<6,6>() - r.covariance_nsg2).norm() == 0, "Wrong marginal covariance");
                require(r.parameter_sigma_amr_si > 0 &&
                        std::abs(r.parameter_sigma_amr_si * r.parameter_sigma_amr_si / r.covariance_full(6,6) - 1) < 1e-14,
                        "Wrong AMR sigma units");
                if (variant == 0)
                {
                    // 第七项为物理 AMR；历元转换必须携带完整交叉协方差，不能只传六阶块。
                    BatchResult moved = r;
                    const double epoch = 59215 + 30.0/86400;
                    MovePosteriorToEpoch(t, fitjob, s, epoch, moved);
                    require(t.Code() == bastion::ERR_NOERROR && moved.fit_epoch_utc == epoch,
                            "Joint posterior epoch transfer failed");
                    require(moved.covariance_full.rows() == 7 && moved.covariance_full.allFinite() &&
                            moved.covariance_full(6,6) == r.covariance_full(6,6), "Joint covariance lost AMR variance");
                    MovePosteriorToEpoch(t, fitjob, s, 59215, moved);
                    require((moved.covariance_full-r.covariance_full).norm()/r.covariance_full.norm()<5e-4,
                            "Joint covariance epoch round-trip failed");
                }
                if (reference.converged)
                    require(std::abs(reference.parameter_amr_si - r.parameter_amr_si) < 1e-6, "FD/integration sensitivity too large");
                reference = r;
            }
            // 非法选择在传播前拒绝；短弧自动模式保留六状态。
            auto invalid = [&](OrbitPropagationJob j, BatchSetting s, std::vector<Observation> o, const std::string& why) {
                bastion::Overwatch local;
                BatchResult r;
                require(!BatchOrbitDetermination_NSG2(local, j, o, s, r) && r.message.find(why) != std::string::npos,
                        "Invalid parameter setup accepted: " + why + "; " + r.message);
            };
            BatchSetting s;
            s.display_iterations = false;
            s.estimate_parameter = drag ? 1 : 2;
            s.stm_method = 1;
            invalid(job, s, obs, "stm_method");
            s.stm_method = 2;
            auto disabled = job;
            disabled.opn_model.use_drag = disabled.opn_model.use_SRP = false;
            invalid(disabled, s, obs, "force model");
            auto zero = job;
            if (drag) zero.opn_model.sat.AMR_drag = 0; else zero.opn_model.sat.AMR_SRP = 0;
            invalid(zero, s, obs, "positive prior");
            s.estimate_parameter = 4;
            invalid(job, s, obs, "additional parameter settings");
            // 短于三圈时不增加参数；六状态协方差仍应正常计算。
            s.estimate_parameter = 3;
            auto short_obs = std::vector<Observation>(obs.begin(), obs.begin() + 15);
            BatchResult short_result;
            require(BatchOrbitDetermination_NSG2(t, job, short_obs, s, short_result), "Short auto fit failed");
            require(short_result.estimated_parameter == 0 && short_result.covariance_full.rows() == 6,
                    "Auto estimated a seventh parameter on a short arc");
        }
        require(t.Code() == bastion::ERR_NOERROR, "Unexpected model error");
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
