#include <array>
// 2026-10-21：批处理定轨科学验证；真值使用独立二体解析推进，不由求解器反推自身。
#include "orbit_geometry.hpp"
#include "core_bodn.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace determination;
namespace
{
void require(bool pass, const std::string &message)
{
    if (!pass)
        throw std::runtime_error(message);
}
orbit::NonSingular2 nsg(bastion::Overwatch &t, const StateVector &z)
{
    return orbit::NonSingular2(t, z[0], z[1], z[2], z[3], z[4], z[5], true);
}
StateVector statevector(const orbit::NonSingular2 &z)
{
    std::array<double, 6> exported_values_1; z.ToArray(exported_values_1.data());
    auto a = exported_values_1;
    return Eigen::Map<const StateVector>(a.data());
}
OrbitPropagationJob makejob(bastion::Overwatch &t, double a = 7200, double e = 0.05, double inc = 45)
{
    OrbitPropagationJob j{};
    std::string path = "data/eop/EOP-All_250516.txt";
    j.opn_EOP.Read_EOP(t, path);
    j.initKep_norm.Set(a, e, inc, 20, 40, 10, false);
    j.initKep_norm.Normalize(t, ae / 1000);
    j.initKep_norm.SyncAnomaly(t);
    j.initKep_norm.ToCartesian(t, j.initCart_norm);
    j.opn_setting.t0_MJD_UTC = 59215;
    j.opn_setting.t0_MJD_TT = UTC2TT(59215);
    j.opn_setting.h0_norm = 10 / UNITTIME;
    j.opn_setting.rkf.maxTE_norm = 1e-13;
    j.opn_setting.display_progress = 0;
    return j;
}
std::vector<Observation> truth(bastion::Overwatch &t, const OrbitPropagationJob &j, int type,
                               bool noise = false)
{
    orbit::NonSingular2 initial;
    j.initKep_norm.ToNonSingular2(t, initial);
    StateVector z = statevector(initial);
    std::vector<Observation> obs;
    std::mt19937 random(20261021);
    std::normal_distribution<double> normal;
    for (int i = 0; i < 91; ++i)
    {
        double seconds = i * 40.0 + (i % 3) * 0.37;
        double epoch = 59215 + seconds / 86400;
        StateVector zi = z;
        zi[5] += pow(z[0], -1.5) * (UTC2TT(epoch) - UTC2TT(59215)) * 86400 / UNITTIME;
        auto state = nsg(t, zi);
        state.SyncAnomaly(t);
        orbit::Cartesian cart;
        state.ToCartesian(t, cart);
        std::array<double, 6> exported_values_2; cart.ToArray(exported_values_2.data());
        auto a = exported_values_2;
        for (int station = 0; station < (type == 7 ? 1 : 2); ++station)
        {
            Observation o;
            o.mjd_utc = epoch;
            o.object_id = "SYNTHETIC";
            o.observer_id = station ? "B" : "A";
            o.type = static_cast<ObservationType>(type);
            Eigen::Vector3d position(a[0] * ae, a[1] * ae, a[2] * ae);
            if (type == 7)
            {
                o.value = position;
                o.sigma.setConstant(5);
                if (noise)
                    for (int k = 0; k < 3; ++k)
                        o.value[k] += 5 * normal(random);
            }
            else
            {
                o.observer_position_gcrs_m =
                    station ? Eigen::Vector3d(-3e6, 5e6, 2e6) : Eigen::Vector3d(6e6, 1e6, -1e6);
                Eigen::Vector3d dr_gcrs_m = position - o.observer_position_gcrs_m;
                o.value << atan2(dr_gcrs_m[1], dr_gcrs_m[0]), atan2(dr_gcrs_m[2], hypot(dr_gcrs_m[0], dr_gcrs_m[1])), 0;
                o.sigma.setConstant(pi / (180 * 3600));
                if (noise)
                {
                    o.value[0] += normal(random) * o.sigma[0] / cos(o.value[1]);
                    o.value[1] += normal(random) * o.sigma[1];
                }
            }
            obs.push_back(o);
        }
    }
    return obs;
}
OrbitPropagationJob perturbed(bastion::Overwatch &t, OrbitPropagationJob j)
{
    orbit::NonSingular2 initial;
    j.initKep_norm.ToNonSingular2(t, initial);
    StateVector z = statevector(initial);
    z[0] += 1e-4;
    z[1] += 2e-5;
    z[2] -= 3e-5;
    z[3] += 4e-5;
    z[4] -= 2e-5;
    z[5] += 1e-4;
    auto state = nsg(t, z);
    state.SyncAnomaly(t);
    state.ToKeplerian(t, j.initKep_norm);
    state.ToCartesian(t, j.initCart_norm);
    return j;
}
} // namespace
int main(int argc, char **argv)
{
    try
    {
        bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
        // 用独立三次多项式检验公用插值的位置及速度，不调用动力学或积分器。
        std::vector<double> polynomial_epochs;
        std::vector<orbit::Cartesian> polynomial_nodes;
        const double polynomial_epoch0_tt = UTC2TT(59215);

        for (double seconds : {0.0, 20.0, 57.0, 100.0})
        {
            const double epoch = 59215 + seconds / 86400;
            const double t = (UTC2TT(epoch) - polynomial_epoch0_tt) * 86400 / UNITTIME;
            polynomial_epochs.push_back(epoch);
            polynomial_nodes.emplace_back(1 + t + 0.5*t*t + 0.03*t*t*t, 2 - 0.2*t + 0.4*t*t, -1 + 0.4*t - 0.07*t*t*t,
                                          1 + t + 0.09*t*t, -0.2 + 0.8*t, 0.4 - 0.21*t*t, true);
        }

        const std::vector<double> sample_epochs = {59215, 59215 + 7.0/86400, 59215 + 7.0/86400,
                                                   59215 + 49.0/86400, polynomial_epochs.back()};
        const auto polynomial_samples = Interpolate_CartesianOrbit(tracer, polynomial_epochs, polynomial_nodes, sample_epochs);
        require(tracer.Code() == bastion::ERR_NOERROR && polynomial_samples.size() == sample_epochs.size(),
                "Polynomial orbit interpolation failed.");

        for (std::size_t i = 0; i < sample_epochs.size(); ++i)
        {
            const double t = (UTC2TT(sample_epochs[i]) - polynomial_epoch0_tt) * 86400 / UNITTIME;
            const double expected[6] = {1 + t + 0.5*t*t + 0.03*t*t*t, 2 - 0.2*t + 0.4*t*t, -1 + 0.4*t - 0.07*t*t*t,
                                        1 + t + 0.09*t*t, -0.2 + 0.8*t, 0.4 - 0.21*t*t};
            double actual[6];
            polynomial_samples[i].ToArray(actual);

            for (int j = 0; j < 6; ++j)
                require(std::abs(actual[j] - expected[j]) < 1e-12, "Hermite position/velocity polynomial mismatch.");
        }

        // UTC→TT→UTC 应保留 MJD 的亚微秒分辨率，不能先折叠为完整 JD。
        for (double utc : {55256.123456789, 59215.00000428241, 61334.999876543})
            require(std::abs(TT2UTC(UTC2TT(utc)) - utc) * 86400 < 1e-6, "UTC/TT round trip lost precision.");
        // 独立已知几何：相对视线 (-3, 3, sqrt(6)) 的距离 sqrt(24)，RA=135°、Dec=30°。
        const Eigen::Vector3d line_of_sight(-3, 3, std::sqrt(6.0));
        const auto equatorial = LOS_to_Equatorial(line_of_sight);
        require(std::abs(equatorial.rho - std::sqrt(24.0)) < 1e-14 &&
                    std::abs(equatorial.right_ascension - 3 * pi / 4) < 1e-14 &&
                    std::abs(equatorial.declination - pi / 6) < 1e-14,
                "Relative line-of-sight geometry mismatch.");
        double equator[3];
        RelativePosition_to_Equatorial(7, 23, 30 + std::sqrt(6.0), 10, 20, 30, equator);
        require(std::abs(equator[0] - equatorial.rho) < 1e-14 &&
                    std::abs(equator[1] - equatorial.right_ascension) < 1e-14 &&
                    std::abs(equator[2] - equatorial.declination) < 1e-14,
                "RelativePosition_to_Equatorial observer subtraction mismatch.");
        const auto zero_line = LOS_to_Equatorial(Eigen::Vector3d::Zero());
        const auto pole_line = LOS_to_Equatorial(Eigen::Vector3d(0, 0, -2));
        require(zero_line.rho == 0 && std::isnan(zero_line.right_ascension) &&
                    std::isnan(zero_line.declination), "Zero line of sight must have undefined angles.");
        require(pole_line.rho == 2 && pole_line.right_ascension == 0 &&
                    pole_line.declination == -pi / 2, "Polar line-of-sight convention mismatch.");
        const double high_ra = 359.999 * pi / 180, low_ra = 0.001 * pi / 180;
        require(std::abs(std::remainder(high_ra - low_ra, 2 * pi) + 0.002 * pi / 180) < 1e-14 &&
                    std::abs(std::remainder(low_ra - high_ra, 2 * pi) - 0.002 * pi / 180) < 1e-14,
                "Signed RA difference across zero changed sign or magnitude.");
        auto job = makejob(tracer);
        orbit::NonSingular2 state;
        job.initKep_norm.ToNonSingular2(tracer, state);
        const StateVector z = statevector(state);
        const auto G = GeometryPartialsNSG2(tracer, state);
        double geometry_error = 0;
        for (double a : {7200.0, 42164.0})
            for (double e : {0.0, 0.6})
                for (double inc : {0.0, 45.0, 120.0})
                {
                    auto j = makejob(tracer, a, e, inc);
                    orbit::NonSingular2 s;
                    j.initKep_norm.ToNonSingular2(tracer, s);
                    StateVector x = statevector(s);
                    const StateMatrix analytic = GeometryPartialsNSG2(tracer, s);
                    for (int k = 0; k < 6; ++k)
                    {
                        StateVector p = x, m = x;
                        p[k] += 1e-6;
                        m[k] -= 1e-6;
                        auto sp = nsg(tracer, p), sm = nsg(tracer, m);
                        sp.SyncAnomaly(tracer);
                        sm.SyncAnomaly(tracer);
                        orbit::Cartesian cp, cm;
                        sp.ToCartesian(tracer, cp);
                        sm.ToCartesian(tracer, cm);
                        std::array<double, 6> exported_values_3; cm.ToArray(exported_values_3.data());
                        std::array<double, 6> exported_values_4; cp.ToArray(exported_values_4.data());
                        auto ap = exported_values_4, am = exported_values_3;
                        StateVector fd;
                        for (int q = 0; q < 6; ++q)
                            fd[q] = (ap[q] - am[q]) / 2e-6;
                        double err = (fd - analytic.col(k)).norm() / std::max(1.0, analytic.col(k).norm());
                        if (err > 2e-7)
                            std::cerr << "geometry a/e/i/column/error " << a << " " << e << " " << inc << " "
                                      << k << " " << err << "\nFD " << fd.transpose() << "\nG "
                                      << analytic.col(k).transpose() << "\n";
                        geometry_error = std::max(geometry_error, err);
                    }
                }
        require(geometry_error < 2e-7, "Geometry Jacobian finite-difference mismatch.");
        std::cout << std::setprecision(12) << "geometry_max_relative_error " << geometry_error << '\n';
        auto observations = truth(tracer, job, 6);
        // 对理论观测量 calc_obs 做中心差分，东西方向使用观测 Dec 的余弦因子，再与 H 比较。
        Eigen::Vector3d calc_obs;
        Eigen::Matrix<double, 3, 6> H;
        Compute_ObsAndPartials(tracer, observations[0], job.initCart_norm, calc_obs, H);
        double observation_error = 0;
        std::array<double, 6> cart;
        job.initCart_norm.ToArray(cart.data());

        for (int k = 0; k < 3; ++k)
        {
            auto p = cart, m = cart;
            p[k] += 1e-7;
            m[k] -= 1e-7;
            Eigen::Vector3d calc_obs_plus, calc_obs_minus;
            Eigen::Matrix<double, 3, 6> unused;
            Compute_ObsAndPartials(tracer, observations[0], orbit::Cartesian(p[0], p[1], p[2], p[3], p[4], p[5], true), calc_obs_plus, unused);
            Compute_ObsAndPartials(tracer, observations[0], orbit::Cartesian(m[0], m[1], m[2], m[3], m[4], m[5], true), calc_obs_minus, unused);

            Eigen::Vector3d difference = (calc_obs_plus - calc_obs_minus) / 2e-7;
            difference[0] = std::remainder(calc_obs_plus[0] - calc_obs_minus[0], 2 * pi) * std::cos(observations[0].value[1]) / 2e-7;
            observation_error = std::max(observation_error, (difference - H.col(k)).norm());
        }

        require(observation_error < 1e-7, "Observation Jacobian mismatch.");

        // 理论 RA 应按 [0,2pi) 返回，不随输入实测 RA 改变；零点两侧均用已知几何检查。
        Observation wrap;
        wrap.type = ObservationType::SpaceRaDec;
        wrap.value << high_ra, 0, 0;
        Compute_ObsAndPartials(tracer, wrap, orbit::Cartesian(cos(low_ra), sin(low_ra), 0, 0, 0, 0, true), calc_obs, H);
        require(std::abs(calc_obs[0] - low_ra) < 1e-12 && std::abs(calc_obs[1]) < 1e-12, "Computed RA/Dec near zero mismatch.");
        wrap.value[0] = low_ra;
        Compute_ObsAndPartials(tracer, wrap, orbit::Cartesian(cos(high_ra), sin(high_ra), 0, 0, 0, 0, true), calc_obs, H);
        require(std::abs(calc_obs[0] - high_ra) < 1e-12 && std::abs(calc_obs[1]) < 1e-12, "Computed RA near 2pi mismatch.");

        // 观测模型的定义域错误通过 tracer 返回，不抛出异常或进入后续求解。
        for (int geometry : {0, 1, 2})
        {
            bastion::Overwatch invalid_geometry;
            Observation invalid = wrap;
            orbit::Cartesian test_state(1, 0, 0, 0, 0, 0, true);
            if (geometry == 0)
            {
                invalid.type = static_cast<ObservationType>(99);
            }
            else if (geometry == 1)
            {
                invalid.observer_position_gcrs_m = Eigen::Vector3d(ae, 0, 0);
            }
            else
            {
                test_state.Set(0, 0, 1, 0, 0, 0, true);
            }

            Compute_ObsAndPartials(invalid_geometry, invalid, test_state, calc_obs, H);
            require(invalid_geometry.Code() != bastion::ERR_NOERROR, "Invalid observation geometry was accepted.");
        }

        std::cout << "observation_jacobian_error " << observation_error << '\n';
        // 非网格历元用独立二体真值验证 Hermite 插值；四种积分器均参与。
        auto fitjob = perturbed(tracer, job);
        for (int id : {101, 102, 201, 202})
        {
            for (int method : {1, 2})
            {
                BatchSetting settings;
                settings.integrator_id = id;
                settings.stm_method = method;
                settings.display_iterations = false;
                settings.rms_tolerance = 1e-6;
                BatchResult result;
                bool ok = BatchOrbitDetermination_NSG2(tracer, fitjob, observations, settings, result);
                require(ok, "Noise-free optical fit failed: " + result.message);
                require(!result.node_epochs_utc.empty() && result.node_epochs_utc.size() == result.node_nsg2_norm.size() &&
                            result.node_epochs_utc.size() == result.node_cart_norm.size(), "Final integration nodes were not retained.");
                const auto sampled = Interpolate_CartesianOrbit(tracer, result.node_epochs_utc, result.node_cart_norm, result.epochs_utc);
                require(tracer.Code() == bastion::ERR_NOERROR && sampled.size() == result.ephemeris_norm.size(),
                        "Final nodes do not cover observation epochs.");

                // 保留节点必须对应最终接受解，且插值后与求残差时的状态逐分量一致。
                for (std::size_t i = 0; i < sampled.size(); ++i)
                {
                    double a[6], b[6];
                    sampled[i].ToArray(a);
                    result.ephemeris_norm[i].ToArray(b);

                    for (int j = 0; j < 6; ++j) require(a[j] == b[j], "Final trajectory and residual states are inconsistent.");
                }

                double first_node[6], final_initial[6];
                result.node_cart_norm.front().ToArray(first_node);
                result.fitted_cart_norm.ToArray(final_initial);

                for (int j = 0; j < 6; ++j)
                    require(std::abs(first_node[j] - final_initial[j]) < 1e-13, "Retained trajectory belongs to a different initial orbit.");

                std::array<double, 6> exported_values_6; job.initCart_norm.ToArray(exported_values_6.data());
                std::array<double, 6> exported_values_7; result.fitted_cart_norm.ToArray(exported_values_7.data());
                auto x = exported_values_7, y = exported_values_6;
                double pos = 0, vel = 0;
                for (int k = 0; k < 3; ++k)
                {
                    pos += pow((x[k] - y[k]) * ae, 2);
                    vel += pow((x[k + 3] - y[k + 3]) * ae / UNITTIME, 2);
                }
                require(sqrt(pos) < 0.002 && sqrt(vel) < 2e-6, "Noise-free recovered orbit inaccurate.");
                require(result.covariance_cart_si.allFinite() &&
                            result.covariance_cart_si.diagonal().minCoeff() > 0,
                        "Invalid formal covariance.");
                Eigen::Vector3d calc_obs, checkres = Eigen::Vector3d::Zero();
                Eigen::Matrix<double, 3, 6> unused;
                Compute_ObsAndPartials(tracer, observations[0], result.fitted_cart_norm, calc_obs, unused);
                checkres[0] = std::remainder(observations[0].value[0] - calc_obs[0], 2 * pi) * std::cos(observations[0].value[1]);
                checkres[1] = observations[0].value[1] - calc_obs[1];
                require((checkres - result.residuals[0]).norm() < 1e-12,
                        "Final residual/state not from same orbit.");
                std::cout << "optical id " << id << " stm " << method << " position_m " << sqrt(pos)
                          << " velocity_mps " << sqrt(vel) << " wrms " << result.weighted_rms << " rounds "
                          << result.iterations.size() << '\n';
            }
        }
        // 更高半长轴和偏心率：独立二体真值检查 GEO 与 HEO 的收敛。
        for (bool high_e : {false, true})
        {
            auto j = makejob(tracer, high_e ? 26000 : 42164, high_e ? 0.6 : 0.001, high_e ? 63.4 : 0.1);
            auto obs = truth(tracer, j, 6);
            for (int id : {101, 102, 201, 202})
            {
                BatchSetting s;
                s.integrator_id = id;
                s.display_iterations = false;
                s.rms_tolerance = 1e-6;
                BatchResult result;
                require(BatchOrbitDetermination_NSG2(tracer, perturbed(tracer, j), obs, s, result),
                        "GEO/HEO fit failed: " + result.message);
                std::array<double, 6> exported_values_8; j.initCart_norm.ToArray(exported_values_8.data());
                std::array<double, 6> exported_values_9; result.fitted_cart_norm.ToArray(exported_values_9.data());
                auto x = exported_values_9, y = exported_values_8;
                double pos = 0;
                for (int k = 0; k < 3; ++k)
                    pos += pow((x[k] - y[k]) * ae, 2);
                require(sqrt(pos) < 0.02, "GEO/HEO truth recovery inaccurate.");
                std::cout << (high_e ? "HEO" : "GEO") << " id " << id << " position_m " << sqrt(pos)
                          << " WRMS " << result.weighted_rms << '\n';
            }
        }
        // 固定种子噪声和少量粗差：检验统计量及延迟剔除，不要求噪声轨道等于真值。
        auto noisy = truth(tracer, job, 6, true);
        for (int i : {10, 70, 130})
            noisy[i].value[0] += 200 * noisy[i].sigma[0] / cos(noisy[i].value[1]);
        BatchSetting settings;
        settings.stm_method = 2;
        settings.rejection_scale = 5;
        settings.display_iterations = false;
        BatchResult noisy_result;
        require(BatchOrbitDetermination_NSG2(tracer, fitjob, noisy, settings, noisy_result),
                "Noisy fit failed.");
        require(noisy_result.weighted_rms > 0.7 && noisy_result.weighted_rms < 1.3,
                "Noise WRMS implausible.");
        require(std::count(noisy_result.accepted.begin(), noisy_result.accepted.end(), false) == 3,
                "Outlier rejection mismatch.");
        std::cout << "noisy_WRMS " << noisy_result.weighted_rms << " rejected 3\n";
        auto xyz = truth(tracer, job, 7);
        BatchResult xyz_result;
        settings.rejection_scale = 0;
        require(BatchOrbitDetermination_NSG2(tracer, fitjob, xyz, settings, xyz_result), "XYZ fit failed.");
        std::cout << "XYZ_WRMS " << xyz_result.weighted_rms << '\n';
        // 输出历元与先验一致：检查前后 30s 的独立二体真值，并检查协方差往返传递。
        for (int method : {1, 2})
            for (double seconds : {-30.0, 30.0})
            {
                BatchResult moved = xyz_result;
                BatchSetting transfer = settings;
                transfer.stm_method = method;
                const double target_epoch = 59215 + seconds / 86400;
                MovePosteriorToEpoch(tracer, fitjob, transfer, target_epoch, moved);
                require(tracer.Code() == bastion::ERR_NOERROR && moved.fit_epoch_utc == target_epoch, "Posterior epoch label/transfer failed.");
                StateVector expected_z = statevector(state);
                expected_z[5] += std::pow(expected_z[0], -1.5) * (UTC2TT(target_epoch)-UTC2TT(59215)) * 86400 / UNITTIME;
                auto expected = nsg(tracer, expected_z);
                expected.SyncAnomaly(tracer);
                orbit::Cartesian expected_cart;
                expected.ToCartesian(tracer, expected_cart);
                std::array<double, 6> a, b;
                expected_cart.ToArray(a.data()); moved.fitted_cart_norm.ToArray(b.data());
                double error = 0;
                for (int j=0; j<3; ++j) error += std::pow((a[j]-b[j])*ae,2);
                require(std::sqrt(error)<0.002, "Reference-epoch posterior disagrees with independent truth.");
                require(moved.covariance_full.allFinite() && (moved.covariance_full.diagonal().array()>0).all(), "Transferred covariance invalid.");
                require(moved.accepted == xyz_result.accepted && moved.residuals == xyz_result.residuals, "Epoch transfer changed observations.");
                require(moved.node_epochs_utc == xyz_result.node_epochs_utc, "Reference epoch transfer changed the stored arc.");

                for (std::size_t i = 0; i < moved.node_epochs_utc.size(); ++i)
                {
                    double a[6], b[6];
                    moved.node_cart_norm[i].ToArray(a);
                    xyz_result.node_cart_norm[i].ToArray(b);

                    for (int j = 0; j < 6; ++j) require(a[j] == b[j], "Reference epoch transfer changed integration nodes.");
                }

                MovePosteriorToEpoch(tracer, fitjob, transfer, 59215, moved);
                require((moved.covariance_full-xyz_result.covariance_full).norm()/xyz_result.covariance_full.norm()<5e-4,
                        "Covariance epoch round-trip failed.");
            }
        // 摄动真值走 Cartesian 单步积分，逐个到达观测历元，不使用被测 NSG2 的网格/插值。
        // 保持相同力模型，专门检查状态表示、传播及定轨链连接是否一致。
        auto forcejob = makejob(tracer, 7000, 0.01, 55);
        forcejob.opn_model.use_earth_NSG = true;
        forcejob.opn_model.earth_zonal_order = 4;
        forcejob.opn_model.earth_tesseral_degree = 4;
        forcejob.opn_model.solar_gravity_config = 2;
        forcejob.opn_model.lunar_gravity_config = 2;
        forcejob.opn_model.use_drag = true;
        auto &sat = forcejob.opn_model.sat;
        sat.mass = 100 / MASS_OF_EARTH;
        sat.area_drag = 1 / (ae * ae);
        sat.AMR_drag = sat.area_drag / sat.mass;
        sat.Cd = 2.2;
        double density_step = 0;
        require(Read_AtmosphereDensity("data/atmosphere/atmosphere.txt",
                                       forcejob.atm_list, density_step) > 0,
                "Atmosphere fixture missing.");
        Set_ReferenceDensity(forcejob.atm_list, forcejob.initKep_norm, forcejob.refDensity);
        auto truthjob = forcejob;
        truthjob.opn_setting.rkf.fixed_step = true;
        std::array<double, 6> exported_values_10; truthjob.initCart_norm.ToArray(exported_values_10.data());
        auto x = exported_values_10;
        double elapsed = 0;
        std::vector<Observation> forceobs;
        for (int i = 0; i < 61; ++i)
        {
            Observation obs;
            obs.type = ObservationType::InertialPosition;
            obs.object_id = "PERTURBED";
            obs.observer_id = "XYZ";
            obs.mjd_utc = 59215 + (i * 60 + (i % 3) * 0.31) / 86400.;
            obs.sigma.setConstant(1);
            const double target = (UTC2TT(obs.mjd_utc) - UTC2TT(59215)) * 86400 / UNITTIME;
            while (target - elapsed > 1e-13)
            {
                double accepted = 0, step = std::min(2 / UNITTIME, target - elapsed);
                RKF78_Cart(tracer, truthjob, UTC2TT(59215) + elapsed * UNITTIME / 86400, x.data(), step,
                           accepted);
                require(accepted > 0 && tracer.Code() == bastion::ERR_NOERROR,
                        "Independent Cartesian truth failed.");
                elapsed += accepted;
            }
            for (int k = 0; k < 3; ++k)
                obs.value[k] = x[k] * ae;
            forceobs.push_back(obs);
        }
        for (int id : {101, 102, 201, 202})
            for (double fd : {1e-6, 3e-7})
            {
                BatchSetting s;
                s.integrator_id = id;
                s.fd_step = fd;
                s.display_iterations = false;
                s.rms_tolerance = 1e-6;
                BatchResult result;
                require(
                    BatchOrbitDetermination_NSG2(tracer, perturbed(tracer, forcejob), forceobs, s, result),
                    "Perturbed fit failed: " + result.message);
                std::array<double, 6> exported_values_11; forcejob.initCart_norm.ToArray(exported_values_11.data());
                std::array<double, 6> exported_values_12; result.fitted_cart_norm.ToArray(exported_values_12.data());
                auto fitted = exported_values_12, reference = exported_values_11;
                double distance = 0;
                for (int k = 0; k < 3; ++k)
                    distance += pow((fitted[k] - reference[k]) * ae, 2);
                require(sqrt(distance) < 0.003 && result.weighted_rms < 0.005,
                        "Perturbed orbit recovery inaccurate.");
                std::cout << "perturbed id " << id << " fd " << fd << " position_m " << sqrt(distance)
                          << " WRMS " << result.weighted_rms << '\n';
            }
        // 高偏心率的广义开普勒方程必须真正收敛，不能只达到迭代次数便返回。
        for (double eccentricity : {0., 0.6, 0.9, 0.99})
            for (double longitude : {0.01, 1., 4., 6.2})
            {
                double xi = eccentricity * cos(0.4), eta = eccentricity * sin(0.4);
                double u = SolveGeneralizedKepler(xi, eta, longitude, 1e-18, 30);
                require(std::abs(std::remainder(u - xi * sin(u) + eta * cos(u) - longitude, 2 * pi)) < 3e-14,
                        "Generalized Kepler did not converge.");
            }
        // 直接测试残差统计的数学边界，不依赖轨道转换能否生成机器精确的零残差。
        require(ResidualRmsChange(2, 1, false) == 1, "Absolute DRMS definition failed");
        require(ResidualRmsChange(1, 2, true) == 0.5, "Relative DRMS must divide by current RMS");
        require(ResidualRmsChange(0, 0, true) == 0, "Relative zero/zero must be zero");
        require(std::isinf(ResidualRmsChange(1, 0, true)), "Relative nonzero/zero must be infinite");
        require(std::isnan(ResidualRmsChange(-1, 0, false)), "Negative RMS was accepted");
        require(std::isnan(ResidualRmsChange(std::numeric_limits<double>::quiet_NaN(), 0, true)), "NaN RMS was accepted");
        std::cout << "Absolute/relative DRMS, zero/zero, nonzero/zero and invalid-domain checks PASS\n";
        // 数值秩不足、非法配置和乱序资料必须拒绝；使用独立 tracer 不污染后续成功案例。
        auto invalid = [&](std::vector<Observation> obs, BatchSetting s)
        {
            bastion::Overwatch t;
            BatchResult result;
            return !BatchOrbitDetermination_NSG2(t, fitjob, obs, s, result);
        };
        auto same = std::vector<Observation>(10, observations[0]);
        require(invalid(same, settings), "Rank-deficient case accepted.");
        auto unordered = observations;
        std::swap(unordered[0], unordered.back());
        require(invalid(unordered, settings), "Unordered data accepted.");
        auto bad = settings;
        bad.stm_method = 0;
        require(invalid(observations, bad), "Invalid STM accepted.");
        auto poles = observations;
        poles[0].value[1] = pi / 2;
        require(invalid(poles, settings), "Polar RA singularity accepted.");
        // 可供真实 CLI 重跑的输入，保留东西、南北方向的 sigma 约定。
        if (argc > 1)
        {
            std::filesystem::create_directories(argv[1]);
            std::ofstream f(std::filesystem::path(argv[1]) / "observations.txt");
            f << std::setprecision(17);
            for (const auto &o : observations)
                f << o.mjd_utc << " 6 " << o.object_id << ' ' << o.observer_id << ' ' << o.value[0] * 180 / pi
                  << ' ' << o.value[1] * 180 / pi << ' ' << o.observer_position_gcrs_m.transpose() << '\n';
            std::ofstream tf(std::filesystem::path(argv[1]) / "truth.txt");
            tf << std::setprecision(17);
            std::array<double, 6> exported_values_13; job.initCart_norm.ToArray(exported_values_13.data());
            auto a = exported_values_13;
            for (int j = 0; j < 6; ++j)
                tf << a[j] * ae / (j < 3 ? 1 : UNITTIME) << ' ';
            tf << '\n';
            std::ofstream kf(std::filesystem::path(argv[1]) / "initial_keplerian.txt");
            auto kep = fitjob.initKep_norm;
            kep.Restore(tracer, ae / 1000);
            kf << std::setprecision(17);
            std::array<double, 6> exported_values_14; kep.ToArray(exported_values_14.data());
            for (double x : exported_values_14)
                kf << x << ' ';
            kf << '\n';
            // 导出含摄动的独立 Cartesian 真值，用于验证 CLI 的力模型配置确实生效。
            std::ofstream force_file(std::filesystem::path(argv[1]) / "perturbed_observations.txt"),
                force_initial(std::filesystem::path(argv[1]) / "perturbed_initial.txt"),
                force_truth(std::filesystem::path(argv[1]) / "perturbed_truth.txt");
            force_file << std::setprecision(17);
            for (const auto &o : forceobs)
                force_file << o.mjd_utc << " 7 " << o.object_id << ' ' << o.observer_id << ' '
                           << o.value.transpose() << '\n';
            auto perturbed_kep = perturbed(tracer, forcejob).initKep_norm;
            perturbed_kep.Restore(tracer, ae / 1000);
            force_initial << std::setprecision(17);
            std::array<double, 6> exported_values_15; perturbed_kep.ToArray(exported_values_15.data());
            for (double value : exported_values_15)
                force_initial << value << ' ';
            force_truth << std::setprecision(17);
            std::array<double, 6> exported_values_16; forcejob.initCart_norm.ToArray(exported_values_16.data());
            auto reference = exported_values_16;
            for (int k = 0; k < 6; ++k)
                force_truth << reference[k] * ae / (k < 3 ? 1 : UNITTIME) << ' ';
            std::ofstream ground(std::filesystem::path(argv[1]) / "ground.txt"),
                space(std::filesystem::path(argv[1]) / "ground_as_space.txt"),
                stations(std::filesystem::path(argv[1]) / "stations.txt");
            ground << std::setprecision(17);
            space << std::setprecision(17);
            stations << std::setprecision(17) << "A 1 30 120 100\nB 2 0 " << ae << " 0\n";
            for (auto o : observations)
            {
                aEOP eop;
                job.opn_EOP.Get_EOP(tracer, o.mjd_utc, eop);
                Eigen::Matrix3d rotation;
                ITRS_to_GCRS_Matrix(o.mjd_utc, eop, rotation);
                // 独立椭球公式给出 A 站位置，用非零纬经高核对 CLI 的单位和三列顺序。
                double latitude = 30 * pi / 180, longitude = 120 * pi / 180,
                       N = ae / sqrt(1 - E_SQUARED * pow(sin(latitude), 2));
                Eigen::Vector3d fixed = o.observer_id == "A"
                                            ? Eigen::Vector3d((N + 100) * cos(latitude) * cos(longitude),
                                                              (N + 100) * cos(latitude) * sin(longitude),
                                                              (N * (1 - E_SQUARED) + 100) * sin(latitude))
                                            : Eigen::Vector3d(0, ae, 0);
                StateVector zi = z;
                zi[5] += pow(z[0], -1.5) * (UTC2TT(o.mjd_utc) - UTC2TT(59215)) * 86400 / UNITTIME;
                auto ns = nsg(tracer, zi);
                ns.SyncAnomaly(tracer);
                orbit::Cartesian c;
                ns.ToCartesian(tracer, c);
                std::array<double, 6> exported_values_17; c.ToArray(exported_values_17.data());
                auto a = exported_values_17;
                Eigen::Vector3d position = Eigen::Vector3d(a[0], a[1], a[2]) * ae;
                o.observer_position_gcrs_m = rotation * fixed;
                Eigen::Vector3d dr_gcrs_m = position - o.observer_position_gcrs_m;
                double ra = atan2(dr_gcrs_m[1], dr_gcrs_m[0]) * 180 / pi,
                       dec = atan2(dr_gcrs_m[2], hypot(dr_gcrs_m[0], dr_gcrs_m[1])) * 180 / pi;
                ground << o.mjd_utc << " 1 " << o.object_id << ' ' << o.observer_id << ' ' << ra << ' ' << dec
                       << '\n';
                space << o.mjd_utc << " 6 " << o.object_id << ' ' << o.observer_id << ' ' << ra << ' ' << dec
                      << ' ' << o.observer_position_gcrs_m.transpose() << '\n';
            }
        }
        require(tracer.Code() == bastion::ERR_NOERROR, "Unexpected tracer error.");
        std::cout << "PASS batch determination\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL " << e.what() << '\n';
        return 1;
    }
}
