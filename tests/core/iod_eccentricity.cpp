// 偏心率偏导修正的专项测试：对照 Mathematica 参考值，并保存同一批定轨算例的结果。
#include "core_iod_solver.hpp"
#include "frame_definition.hpp"
#include <array>
#include <filesystem>
#include <random>
#include <sstream>
#include <stdexcept>

using namespace Eigen;
using namespace iod;

namespace
{
// 通过独立的偏心率矢量计算偏心率和约束差额；状态使用 mu=1 的归一化单位。
void Compute_Constraint(const VectorXd &state, double target, double &eccentricity, double &value)
{
    const Vector3d r = state.head<3>();
    const Vector3d v = state.tail<3>();
    const Vector3d evec = v.cross(r.cross(v)) - r / r.norm();
    eccentricity = evec.norm();
    value = evec.squaredNorm() - target * target;
}

void Require(bool condition, const std::string &message)
{
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        Require(argc == 3 || (argc == 4 && std::string(argv[3]) == "--baseline"),
                "Usage: aster_test_iod_eccentricity reference.tsv output_directory [--baseline]");
        const bool baseline = argc == 4;
        const std::filesystem::path folder = argv[2];
        std::filesystem::create_directories(folder);
        std::ifstream reference(argv[1]);
        Require(reference.good(), "Cannot read Mathematica reference table");
        std::ofstream derivative_file(folder / "derivatives.csv");
        derivative_file << "case,component,mathematica,cpp,finite_difference\n"
                        << std::setprecision(17);
        int reference_count = 0, derivative_failures = 0;
        double maximum_reference_error = 0.0, maximum_fd_error = 0.0;
        std::string line;
        while (std::getline(reference, line))
        {
            if (line.empty())
                continue;
            std::istringstream input(line);
            VectorXd state(6), analytic(6);
            for (int j = 0; j < 6; ++j)
                input >> state(j);
            double target, constraint;
            input >> target >> constraint;
            for (int j = 0; j < 6; ++j)
                input >> analytic(j);
            Require(!input.fail(), "Malformed Mathematica reference row");
            ++reference_count;
            MatrixXd matrix = MatrixXd::Zero(3, 6);
            VectorXd rhs = VectorXd::Zero(3);
            const Vector3d los(0.4, 0.5, std::sqrt(0.59));
            const Vector3d station(0.8, 0.1, 0.5);
            const double fg[4] = {1.0, 0.0, 1.0, 0.0};
            bastion::Overwatch tracer;
            iod_solver::Append_ConstraintEquations(
                tracer, IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY, target, los, station, state, fg, matrix, rhs);
            if (tracer.Code() != bastion::ERR_NOERROR)
            {
                ++derivative_failures;
                Require(baseline, "Corrected gradient rejected a valid reference state");
                continue;
            }
            Require(std::abs(rhs(3) + constraint) < 1e-10 * (1 + std::abs(constraint)),
                    "Constraint RHS changed");
            for (int j = 0; j < 6; ++j)
            {
                // 五点差分直接使用偏心率矢量，不复用待检验的解析偏导。
                const double h = 1e-5 * std::max(1.0, std::abs(state(j)));
                double values[4];
                const int offsets[4] = {-2, -1, 1, 2};
                for (int p = 0; p < 4; ++p)
                {
                    VectorXd shifted = state;
                    shifted(j) += offsets[p] * h;
                    double ecc;
                    Compute_Constraint(shifted, target, ecc, values[p]);
                }
                const double difference =
                    (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * h);
                const double reference_error = std::abs(matrix(3, j) - analytic(j));
                const double fd_error = std::abs(matrix(3, j) - difference);
                maximum_reference_error = std::max(maximum_reference_error, reference_error);
                maximum_fd_error = std::max(maximum_fd_error, fd_error);
                derivative_file << reference_count << ',' << j << ',' << analytic(j) << ','
                                << matrix(3, j) << ',' << difference << '\n';
                if (!baseline)
                {
                    Require(reference_error < 1e-11 * (1 + std::abs(analytic(j))),
                            "C++ gradient differs from Mathematica");
                    Require(fd_error < 1e-7 * (1 + std::abs(analytic(j))),
                            "C++ gradient differs from eccentricity-vector difference");
                }
            }
        }
        Require(reference_count == 41, "Expected 41 Mathematica reference states");
        derivative_file.close();

        std::ofstream results(folder / "orbits.csv");
        std::ofstream history(folder / "iterations.csv");
        results << "dataset,truth_e,span_seconds,stations,noise_arcsec,constraint,target_e,flag,"
                   "iterations,actual_e,ecc_error,rms_rad,position_error_m,velocity_error_mps\n"
                << std::setprecision(17);
        history << "dataset,target_e,iteration,actual_e,ecc_error,rms_rad,max_correction\n"
                << std::setprecision(17);
        int dataset = 0, completed = 0, failed = 0, limited = 0, runs = 0;
        for (double truth_e : {0.01, 0.1, 0.4})
            for (double duration : {180.0, 540.0})
                for (int stations : {1, 2})
                    for (double noise : {0.0, 1.0})
                    {
                        ++dataset;
                        const int n = 20;
                        bastion::Overwatch geometry;
                        // 各轨道近地点半径均为 7000 km，e 改变时同步调整半长轴。
                        const double a = 7000.0 / (1.0 - truth_e) / (constants::Scale::UNIT_LENGTH_M / 1000.0);
                        const double inclination = 50.0 * constants::Conversion::DEGREE_TO_RAD,
                                     node = 40.0 * constants::Conversion::DEGREE_TO_RAD;
                        const double argument = 20.0 * constants::Conversion::DEGREE_TO_RAD, mean0 = 30.0 * constants::Conversion::DEGREE_TO_RAD;
                        orbit::Keplerian truth_kep(a, truth_e, inclination, node, argument, mean0,
                                                   true);
                        double truth_array[6];
                        truth_kep.ToCartesian(geometry, truth_array);
                        const VectorXd truth = Map<VectorXd>(truth_array, 6);
                        IodJob base{};
                        base.obs.type = OBS_TYPE_RADEC;
                        base.obs.number = n;
                        base.obs.time_MJDUTC.resize(n);
                        base.obs.tau.resize(n);
                        base.obs.ra_rad.resize(n);
                        base.obs.dec_rad.resize(n);
                        base.observer.positions_gcrs_norm_emat.resize(3, n);
                        base.dynamicsModel = IOD_DYNAMICS_TWOBODY;
                        base.iterMethod = IOD_ITERATION_NEWTON;
                        base.iterTOL = 1e-10;
                        base.iterCountMax = 40;
                        base.epochMJDUTC = 60000.25;
                        base.printIterInfo = false;
                        base.initValue.type = IOD_INITIAL_VALUE_SPECIFIED;
                        VectorXd initial = truth;
                        initial(0) += 1000.0 / constants::Scale::UNIT_LENGTH_M;
                        initial(1) -= 700.0 / constants::Scale::UNIT_LENGTH_M;
                        initial(2) += 400.0 / constants::Scale::UNIT_LENGTH_M;
                        initial(3) += 1.0 / (constants::Scale::UNIT_LENGTH_M / constants::Scale::UNIT_TIME_S);
                        initial(4) -= 0.7 / (constants::Scale::UNIT_LENGTH_M / constants::Scale::UNIT_TIME_S);
                        initial(5) += 0.4 / (constants::Scale::UNIT_LENGTH_M / constants::Scale::UNIT_TIME_S);
                        base.initValue.initCart.Set(initial(0), initial(1), initial(2), initial(3),
                                                    initial(4), initial(5), true);
                        std::mt19937 generator(20261028 + dataset);
                        std::normal_distribution<double> standard_normal(0.0, 1.0);
                        aEOP eop;
                        eop.dat = 37;
                        for (int k = 0; k < n; ++k)
                        {
                            base.obs.time_MJDUTC(k) =
                                base.epochMJDUTC + duration * k / (n - 1) / DAYSEC;
                            base.obs.tau(k) = (base.obs.time_MJDUTC(k) - base.epochMJDUTC) *
                                              constants::Conversion::DAY_TO_NORMALIZED_TIME;
                            const bool second = stations == 2 && k >= n / 2;
                            frame::Geodetic site(second ? 500.0 : 943.0, second ? 100.0 : 86.174633,
                                                 second ? 35.0 : 41.725892, false);
                            site.Normalize(geometry, constants::Scale::UNIT_LENGTH_M);
                            frame::GCRS gcrs;
                            site.ToGCRS(geometry, base.obs.time_MJDUTC(k), eop, gcrs);
                            double station[3];
                            gcrs.ToArray(station);
                            base.observer.positions_gcrs_norm_emat.col(k) =
                                Vector3d(station[0], station[1], station[2]);
                            // 独立二体真值以 M=M0+n*tau 传播，不调用被检验求解器的拉格朗日系数。
                            orbit::Keplerian propagated(
                                a, truth_e, inclination, node, argument,
                                mean0 + std::sqrt(1.0 / (a * a * a)) * base.obs.tau(k), true);
                            double target[6];
                            propagated.ToCartesian(geometry, target);
                            const Vector3d relative = Vector3d(target[0], target[1], target[2]) -
                                                      base.observer.positions_gcrs_norm_emat.col(k);
                            base.obs.ra_rad(k) = std::atan2(relative(1), relative(0)) +
                                                 noise * constants::Conversion::ARCSEC_TO_RAD * standard_normal(generator);
                            base.obs.dec_rad(k) =
                                std::atan2(relative(2), std::hypot(relative(0), relative(1))) +
                                noise * constants::Conversion::ARCSEC_TO_RAD * standard_normal(generator);
                        }
                        base.observer.position_gcrs_epoch_norm =
                            base.observer.positions_gcrs_norm_emat.col(0);
                        Require(geometry.Code() == bastion::ERR_NOERROR,
                                "Geometry preparation failed");
                        // 三个固定偏心率及不使用该约束的 GL-O 对照，起点、观测和阈值保持一致。
                        const std::array<double, 4> targets = {std::max(0.0, truth_e - 0.03),
                                                               truth_e, truth_e + 0.03, -1.0};
                        for (double target : targets)
                        {
                            IodJob job = base;
                            job.constraintType = target < 0 ? IOD_CONSTRAINT_ORIENTATION
                                                            : IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY;
                            job.constraintEcc = target;
                            bastion::Overwatch tracer;
                            iod_solver::InitialOrbitDetermination(tracer, job);
                            ++runs;
                            const double nan = std::numeric_limits<double>::quiet_NaN();
                            double actual_e = nan, value = nan, position_error = nan,
                                   velocity_error = nan;
                            if (tracer.Code() == bastion::ERR_NOERROR && job.solution.flag != -2)
                            {
                                double result[6];
                                job.solution.cart.ToArray(result);
                                const VectorXd state = Map<VectorXd>(result, 6);
                                Compute_Constraint(state, target, actual_e, value);
                                position_error = (state.head<3>() - truth.head<3>()).norm() * constants::Scale::UNIT_LENGTH_M;
                                velocity_error =
                                    (state.tail<3>() - truth.tail<3>()).norm() * constants::Scale::UNIT_LENGTH_M / constants::Scale::UNIT_TIME_S;
                                if (job.solution.flag > 0)
                                    ++completed;
                                else
                                    ++limited;
                            }
                            else
                                ++failed;
                            results << dataset << ',' << truth_e << ',' << duration << ','
                                    << stations << ',' << noise << ',' << job.constraintType << ','
                                    << target << ',' << job.solution.flag << ','
                                    << (job.solution.iterations.empty()
                                            ? 0
                                            : job.solution.iterations.back().iteration)
                                    << ',' << actual_e << ','
                                    << (target >= 0 ? std::abs(actual_e - target) : nan) << ','
                                    << job.solution.RMS_angle << ',' << position_error << ','
                                    << velocity_error << '\n';
                            for (const IodIterationRecord &record : job.solution.iterations)
                            {
                                const VectorXd state =
                                    Map<const VectorXd>(record.r0v0_norm.data(), 6);
                                double ecc, g;
                                Compute_Constraint(state, target, ecc, g);
                                history << dataset << ',' << target << ',' << record.iteration
                                        << ',' << ecc << ','
                                        << (target >= 0 ? std::abs(ecc - target) : nan) << ','
                                        << record.rms_angle_rad << ',' << record.max_dr0v0_norm
                                        << '\n';
                            }
                        }
                    }
        results.close();
        history.close();
        Require(results.good() && history.good(), "Cannot save test results");
        std::cout << std::scientific << std::setprecision(12) << (baseline ? "BASELINE" : "PASS")
                  << " reference_states=" << reference_count
                  << " rejected_reference_states=" << derivative_failures
                  << " max_mathematica_error=" << maximum_reference_error
                  << " max_fd_error=" << maximum_fd_error << " runs=" << runs
                  << " converged=" << completed << " limited=" << limited << " failed=" << failed << '\n';
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "FAIL " << error.what() << '\n';
        return 1;
    }
}
