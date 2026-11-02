/**
 * @file main_generate_observations.cpp
 * @brief 根据地心惯性系星历和地面测站位置，生成用于初轨确定的光学观测弧段。
 * @author 陈酊存
 *
 * 输入每行：ID、MJD(UTC)、六维 Cartesian 或 Keplerian 状态；两种文件只选一个。
 * 未归一化输入使用 km、km/s、deg；观测输出为 MJD(UTC)、类型、目标号、测站号、RA、Dec。
 * 测站从 Observers 读入原始站表；Observations 写出标准观测，角度保持 deg。
 * main 按“读配置 → 读星历与转换单位 → 计算测站 → 筛选与测角 → 分弧输出”展开。
 * 运行：aster_generate_observations -I <配置文件>。
 */
#include "orbit_geometry.hpp"
#include "io_settings.hpp"
#include "frame_definition.hpp"
#include "model_observation.hpp"
#include "math_backend.hpp"
#include "io_observation_arc.hpp"
#include "math_basic.hpp"
#include "model_ephemeris_analytic.hpp"
#include "model_ephemeris_jpl.hpp"
#include "model_propagation.hpp"
#include <sofa.h>
#include <cmath>
#include <array>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace io;
using namespace bastion;
using orbit::Cartesian;
using orbit::Keplerian;
using observation::SimulatedObservationRecord;
using observation::WriteObservationArc;

int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_generate_observations " << ASTER_PROGRAM_VERSION << '\n';
        std::cout << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_generate_observations -I <setting.txt>\n";
        return 2;
    }

    Overwatch tracer(__FILE__, __func__, __LINE__);
    std::ofstream program_log;
    try
    {
        // =================================== 01 选择配置文件 ===================================
        namespace fs = std::filesystem;
        string setting_file;
        fs::path output_dir;

        for (int i = 1; i < argc; ++i)
        {
            const string arg = argv[i];

            // 帮助说明命令行接口；正常运行使用显式给定的配置。
            if (arg == "-h" || arg == "--help")
            {
                cout << "Usage: aster_generate_observations [-I|-input|--input <setting.txt>]" << endl;
                tracer.CheckOut(__func__);
                return 0;
            }

            // 三种写法均只改变 setting 路径，输出位置由其中的设置决定。
            if (arg != "-I" && arg != "-input" && arg != "--input")
                throw runtime_error("Unknown argument: " + arg);

            if (i + 1 >= argc || argv[i + 1][0] == '-' || string(argv[i + 1]).empty())
                throw runtime_error("Missing value for " + arg);

            setting_file = argv[++i];
        }

        fs::path config_path = setting_file;
        const auto start = chrono::steady_clock::now();
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Observation simulation started. Config: " + config_path.string(), EventType::MESSAGE, FORCE_ALERT_ON);

        // =================================== 02 读取星历、测站与输出设置 ===================================
        SimpleInput settings;
        settings.ReadSimpleInput(tracer, setting_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 可选运行日志：省略或留空时不建文件，日志路径相对于 setting 目录。
        std::string fn_programlog;

        if (settings.HasKeyword("file_programlog"))
            settings.GetString(tracer, "file_programlog", fn_programlog);

        if (!fn_programlog.empty())
        {
            const std::filesystem::path log_config = std::filesystem::weakly_canonical(std::filesystem::absolute(setting_file));
            std::filesystem::path log_path = fn_programlog;

            if (log_path.is_relative()) log_path = log_config.parent_path() / log_path;

            log_path = std::filesystem::weakly_canonical(log_path);

            // 日志不能覆盖已有资料，也不能与本次计划写出的其他文件重名。
            if (std::filesystem::exists(log_path))
                throw std::runtime_error("Program log already exists: " + log_path.string());

            for (std::size_t record_index = 1; record_index <= settings.Count(); ++record_index)
            {
                io::SimpleRecord path_record;
                settings.GetRecord(tracer, record_index, path_record);

                if (path_record.kw.compare(0, 5, "file_") != 0 ||
                    path_record.kw == "file_programlog" || path_record.va.empty()) continue;
                std::filesystem::path other_path = path_record.va;

                if (other_path.is_relative()) other_path = log_config.parent_path() / other_path;

                if (std::filesystem::weakly_canonical(other_path) == log_path)
                    throw std::runtime_error("Program log conflicts with " + path_record.kw);
            }

            program_log.open(log_path);

            if (!program_log) throw std::runtime_error("Cannot open program log: " + log_path.string());

            const std::time_t log_started = std::time(nullptr);
            program_log << "Program log\nStarted: " << std::ctime(&log_started)
                        << "Executable: " << argv[0] << "\nConfig: " << log_config.string()
                        << "\nLog: " << log_path.string() << "\n\n[Input settings]\n";
            std::ifstream log_input(log_config);

            if (!log_input) throw std::runtime_error("Cannot read settings for program log.");

            program_log << log_input.rdbuf() << "\n[Resolved file paths]\n";

            for (std::size_t record_index = 1; record_index <= settings.Count(); ++record_index)
            {
                io::SimpleRecord path_record;
                settings.GetRecord(tracer, record_index, path_record);

                if (path_record.kw.compare(0, 5, "file_") != 0 || path_record.va.empty()) continue;

                std::filesystem::path recorded_path = path_record.va;

                if (recorded_path.is_relative()) recorded_path = log_config.parent_path() / recorded_path;

                program_log << path_record.kw << " = "
                            << std::filesystem::weakly_canonical(recorded_path).string() << '\n';
            }

            program_log.flush();

            if (!program_log) throw std::runtime_error("Failed while writing program log.");
        }

        config_path = fs::weakly_canonical(fs::absolute(config_path));
        string fn_cart, fn_kep, fn_eop, fn_observer, configured_output, object_id, observer_id;
        // 【必选】object_id：目标号，用于 OBS 正文及输出文件名；不得留空或含路径分隔符。
        settings.GetString(tracer, "object_id", object_id);

        // 【二选一】两种星历路径：Cartesian 为 km、km/s，Keplerian 为 km、deg（归一化开关为 0 时）。
        settings.GetString(tracer, "file_eph_Cartesian", fn_cart);
        settings.GetString(tracer, "file_eph_Keplerian", fn_kep);

        // 【必选】EOP、测站表、站号和现成输出目录；相对路径均从本配置所在目录解释。
        settings.GetString(tracer, "file_EOP", fn_eop);
        settings.GetString(tracer, "file_observer", fn_observer);
        settings.GetString(tracer, "observer_id", observer_id);
        settings.GetString(tracer, "dir_output", configured_output);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        const bool has_cart = !fn_cart.empty();

        if (has_cart == !fn_kep.empty())
            throw runtime_error("Provide exactly one of file_eph_Cartesian and file_eph_Keplerian.");

        if (object_id.empty() || object_id.find_first_of("/\\ \t\r\n") != string::npos || object_id == "." || object_id == "..")
            throw runtime_error("object_id must be a non-empty identifier without path separators or whitespace.");

        int normalized = 0, output_kep = 1, output_cart = 0, output_all = 0, observation_type = 1;
        double min_elevation_deg = 5.0;
        // 【取值】is_input_normalized：0 为 km、km/s、deg；1 为 ae、UNITTIME、rad，历元始终为 UTC MJD。
        // 普通数值读取器在缺项或无法读取时警告并给零；下方数值关键词遵守同一规则。
        settings.GetInteger(tracer, "is_input_normalized", normalized);

        // 随机误差由配套 Python 脚本生成；旧加噪设置不能留在理论观测配置中。
        if (settings.HasKeyword("noise_sigma_as") || settings.HasKeyword("random_seed"))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                "Noise settings moved to tools/add_angular_noise.py; remove noise_sigma_as and random_seed from this setting.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 【可选】三个输出开关：0 关闭，1 分别写出 Keplerian 真值、Cartesian 真值及全部历元；缺省为 0。
        settings.GetInteger(tracer, "need_output_Keplerian", output_kep);
        settings.GetInteger(tracer, "need_output_Cartesian", output_cart);
        settings.GetInteger(tracer, "need_output_all", output_all);

        // 【必选】observation_type：当前仅 1（地基 GCRS RA/Dec）。
        // 【取值】min_elevation_deg：最低高度角，单位 deg，范围 [0,90)；缺项警告后取 0。
        settings.GetInteger(tracer, "observation_type", observation_type);
        settings.GetDouble(tracer, "min_elevation_deg", min_elevation_deg);

        int use_station_night = 0, use_earth_shadow = 0, solar_ephemeris = SOLAR_GRAVITY_ANALYTICAL;
        double max_sun_elevation_deg = 0.0;
        string fn_jpl;

        // 【可选】夜间、地影开关独立，0 关闭、1 开启；缺省为 0。
        settings.GetInteger(tracer, "use_station_night", use_station_night);
        settings.GetInteger(tracer, "use_earth_shadow", use_earth_shadow);

        // 【条件必选】启用任一筛选时 solar_ephemeris 为 1（JPL/TDB）或 2（解析/UT1）。
        // 夜间筛选另需 max_sun_elevation_deg（太阳高度角上限，deg，[-90,0)）；
        // 选择 JPL 时另需 file_JPLEPH，格式为本库可读的数值历表二进制。
        if (use_station_night || use_earth_shadow) settings.GetInteger(tracer, "solar_ephemeris", solar_ephemeris);

        if (use_station_night) settings.GetDouble(tracer, "max_sun_elevation_deg", max_sun_elevation_deg);

        if ((use_station_night || use_earth_shadow) && solar_ephemeris == SOLAR_GRAVITY_JPL) settings.GetString(tracer, "file_JPLEPH", fn_jpl);

        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        if ((normalized != 0 && normalized != 1) || (output_kep != 0 && output_kep != 1) ||
            (output_cart != 0 && output_cart != 1) || (output_all != 0 && output_all != 1))
            throw runtime_error("is_input_normalized and output switches must be 0 or 1.");

        if (observer_id.empty() || observer_id.find_first_of("/\\ \t\r\n#") != string::npos || observer_id == "." || observer_id == "..")
            throw runtime_error("observer_id must be a non-empty identifier without path separators, whitespace or #.");

        if (min_elevation_deg < 0.0 || min_elevation_deg >= 90.0) throw runtime_error("min_elevation_deg must be in [0, 90).");

        if ((use_station_night != 0 && use_station_night != 1) || (use_earth_shadow != 0 && use_earth_shadow != 1))
            throw runtime_error("use_station_night and use_earth_shadow must be 0 or 1.");

        if (solar_ephemeris != SOLAR_GRAVITY_JPL && solar_ephemeris != SOLAR_GRAVITY_ANALYTICAL)
            throw runtime_error("solar_ephemeris must be 1 (JPL) or 2 (analytical).");

        if (use_station_night && (max_sun_elevation_deg < -90.0 || max_sun_elevation_deg >= 0.0))
            throw runtime_error("Station-night filtering requires max_sun_elevation_deg in [-90, 0); choose it for the observation conditions.");

        const bool need_sun = use_station_night || use_earth_shadow;

        if (observation_type != 1)
            throw runtime_error("Only observation_type=1 (ground-based GCRS RA/Dec) is implemented; types 2..8 are reserved.");

        // 相对输入路径以配置文件所在目录为基准，输出路径同样以配置目录为基准。
        for (auto* filename : {&fn_cart, &fn_kep, &fn_eop, &fn_observer})
        {
            if (filename->empty()) continue;

            const fs::path path(*filename);
            *filename = fs::weakly_canonical(path.is_absolute() ? path : config_path.parent_path() / path).string();

            if (!fs::is_regular_file(*filename)) throw runtime_error("Input file not found: " + *filename);
        }

        if (fn_eop.empty()) throw runtime_error("file_EOP is required.");

        if (fn_observer.empty()) throw runtime_error("file_observer is required.");

        if (need_sun && solar_ephemeris == SOLAR_GRAVITY_JPL)
        {
            if (fn_jpl.empty()) throw runtime_error("file_JPLEPH is required for solar_ephemeris=1 when an optical filter is enabled.");

            const fs::path path(fn_jpl);
            fn_jpl = fs::weakly_canonical(path.is_absolute() ? path : config_path.parent_path() / path).string();

            if (!fs::is_regular_file(fn_jpl)) throw runtime_error("JPL ephemeris not found: " + fn_jpl);

            if (Read_JPLEPH(fn_jpl) < 0) throw runtime_error("Cannot load the configured JPL ephemeris.");

            if (fs::file_size(fn_jpl) < static_cast<uintmax_t>(3 * RECSIZE))
                throw runtime_error("JPL file contains no complete ephemeris record after its headers.");
        }

        if (configured_output.empty()) throw runtime_error("dir_output is required.");

        const fs::path output_path(configured_output);
        output_dir = fs::weakly_canonical(output_path.is_absolute() ? output_path : config_path.parent_path() / output_path);

        // 多弧段沿用目录内的科学命名；目录由外部准备，不自动创建或改名。
        if (!fs::is_directory(output_dir)) throw runtime_error("Prepare output directory before running: " + output_dir.string());

        // =================================== 03 读星历并转换单位 ===================================

        // 内部 Cartesian 使用归一化单位；Keplerian 真值保留 km、deg。
        vector<int> ID;
        vector<double> MJD;
        vector<array<double, 6>> cartesian, keplerian;
        const string& fn_eph = has_cart ? fn_cart : fn_kep;
        ifstream input(fn_eph);

        if (!input) throw runtime_error("Cannot open ephemeris: " + fn_eph);

        string line;
        size_t line_number = 0;

        while (getline(input, line))
        {
            ++line_number;

            if (line.find_first_not_of(" \t\r") == string::npos) continue;

            istringstream row(line);
            int id;
            double mjd;
            array<double, 6> values{}, cart{}, kep{};

            if (!(row >> id >> mjd >> values[0] >> values[1] >> values[2] >> values[3] >> values[4] >> values[5]))
                throw runtime_error("Expected ID, MJD UTC and six state components at ephemeris line " + to_string(line_number));

            string extra;

            if (row >> extra) throw runtime_error("Unexpected extra column at ephemeris line " + to_string(line_number));

            if (!isfinite(mjd)) throw runtime_error("Non-finite ephemeris epoch at line " + to_string(line_number));

            for (double value : values)

                if (!isfinite(value)) throw runtime_error("Non-finite ephemeris state at line " + to_string(line_number));

            if (!MJD.empty() && mjd <= MJD.back()) throw runtime_error("Ephemeris epochs must be strictly increasing.");

            if (has_cart)
            {
                cart = values;

                if (!normalized)
                {
                    for (int j = 0; j < 3; ++j) cart[j] *= KM_TO_NORMALIZED_LENGTH;

                    for (int j = 3; j < 6; ++j) cart[j] *= KMPS_TO_NORMALIZED_VELOCITY;
                }

                if (output_kep)
                {
                    Cartesian state(cart[0], cart[1], cart[2], cart[3], cart[4], cart[5], true);
                    state.ToKeplerian(tracer, kep.data());
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
                    kep[0] *= NORMALIZED_LENGTH_TO_KM;

                    for (int j = 2; j < 6; ++j) kep[j] *= DR2D;
                }
            }
            else
            {
                if (values[0] <= 0.0 || values[1] < 0.0 || values[1] >= 1.0)
                    throw runtime_error("Keplerian input requires a > 0 and 0 <= e < 1 at line " + to_string(line_number));

                Keplerian state(values[0], values[1], values[2], values[3], values[4], values[5], normalized != 0);

                if (!normalized) state.Normalize(tracer, ae / 1000.0);

                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

                state.ToCartesian(tracer, cart.data());
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
                kep = values;

                if (normalized)
                {
                    kep[0] *= NORMALIZED_LENGTH_TO_KM;

                    for (int j = 2; j < 6; ++j) kep[j] *= DR2D;
                }
            }

            for (double value : cart)

                if (!isfinite(value)) throw runtime_error("State conversion produced a non-finite value at line " + to_string(line_number));

            if (Module_3D(cart[0], cart[1], cart[2]) <= 0.0) throw runtime_error("Target geocentric distance must be positive.");

            ID.push_back(id); MJD.push_back(mjd);
            cartesian.push_back(cart); keplerian.push_back(kep);
        }

        if (input.bad()) throw runtime_error("Failed while reading ephemeris: " + fn_eph);

        if (MJD.empty()) throw runtime_error("Input ephemeris is empty.");

        const size_t num = MJD.size();

        // =================================== 04 计算同历元测站与太阳位置 ===================================

        // 完整读取原始站表，再按 observer_id 选择一个站。
        observation::Observers simulation_observers;
        simulation_observers.ReadFile(tracer, fn_observer);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        observation::ObserverRecord ground_observer;
        bool observer_found = false;

        for (const observation::ObserverRecord& candidate : simulation_observers.GetRecords())
        {
            if (candidate.id == observer_id)
            {
                ground_observer = candidate;
                observer_found = true;
                break;
            }
        }

        if (!observer_found) throw runtime_error("Unknown ground observer: " + observer_id);

        // 读取类不做转换；在使用阶段把纬经高或 ITRS XYZ 统一变成归一化 ITRS。
        frame::ITRS site_itrs;

        if (ground_observer.coordinate_type == 1)
        {
            if (ground_observer.geodetic_phi_deg < -90.0 || ground_observer.geodetic_phi_deg > 90.0)
                throw runtime_error("Ground observer latitude must be between -90 and 90 deg.");

            frame::Geodetic site_geodetic(ground_observer.geodetic_height_m, ground_observer.geodetic_lambda_deg,
                                 ground_observer.geodetic_phi_deg, false);
            site_geodetic.Normalize(tracer, ae); // 高度 m→ae，经纬度 deg→rad。
            site_geodetic.ToITRS(tracer, site_itrs);
        }
        else
        {
            site_itrs.Set(ground_observer.position_itrs_m[0], ground_observer.position_itrs_m[1],
                          ground_observer.position_itrs_m[2], false);
            site_itrs.Normalize(tracer, ae); // 地固 XYZ 从 m 变成 ae。
        }

        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        EOP eop;
        eop.Read_EOP(tracer, fn_eop);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        vector<Vector3d> site_positions(num);
        vector<Vector3d> sun_directions;

        if (need_sun) sun_directions.resize(num);

        for (size_t k = 0; k < num; ++k)
        {
            // 每个观测历元从内存表插值 EOP，再计算同历元的测站位置。
            aEOP eop_at_time;
            eop.Get_EOP(tracer, MJD[k], eop_at_time, 1);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            frame::GCRS site_gcrs;
            site_itrs.ToGCRS(tracer, MJD[k], eop_at_time, site_gcrs);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            site_positions[k][0] = site_gcrs.GetValue(tracer, VAR_CARTX);
            site_positions[k][1] = site_gcrs.GetValue(tracer, VAR_CARTY);
            site_positions[k][2] = site_gcrs.GetValue(tracer, VAR_CARTZ);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            if (need_sun)
            {
                // 太阳均取同历元的地心几何方向；沿用 JPL/解析星历轴向与 GCRS 近似一致的处理。
                const double ut1 = MJD[k] + eop_at_time.dut1 / 86400.0;
                double sun[6] = {};

                if (solar_ephemeris == SOLAR_GRAVITY_ANALYTICAL)
                {
                    Cart_MeanSun(ut1, sun);  // 沿用本库分析历表的 UT1 接口，长度归一化。
                }
                else
                {
                    // 当前 EOP 的 TAI-UTC 给出 TT；SOFA 地心 TDB-TT 项用于 JPL 查询。
                    double tt1, tt2, tdb1, tdb2;
                    iauTaitt(START_OF_MJD, MJD[k] + eop_at_time.dat / 86400.0, &tt1, &tt2);
                    const double dtr = iauDtdb(tt1, tt2, ut1 - floor(ut1), 0.0, 0.0, 0.0);
                    iauTttdb(tt1, tt2, dtr, &tdb1, &tdb2);
                    const double tdb = (tdb1 - START_OF_MJD) + tdb2;
                    const double jd = tdb + START_OF_MJD;

                    if (jd < JPLEPH_R1.r1.ss[0] || jd > JPLEPH_R1.r1.ss[1])
                        throw runtime_error("Observation epoch is outside JPL coverage at input row " + to_string(k + 1));

                    if (Get_JPLEPH(tdb, 11, 3, sun) < 0) throw runtime_error("Failed to obtain the geocentric Sun position from JPL.");
                }

                const Vector3d direction(sun[0], sun[1], sun[2]);
                const double distance = direction.norm();

                if (!isfinite(distance) || distance <= 0.0) throw runtime_error("Invalid geocentric Sun position at input row " + to_string(k + 1));

                sun_directions[k] = direction / distance;
            }
        }

        const double station_radius = site_positions.front().norm();

        if (!isfinite(station_radius) || station_radius <= 0.0) throw runtime_error("Station geocentric distance must be positive.");

        // 数字站号仅在文件名中补足四位；正文始终保留站表中的原始编号。
        ostringstream station_label;

        if (observer_id.find_first_not_of("0123456789") == string::npos)
        {
            station_label << setw(4) << setfill('0') << observer_id;
        }
        else
        {
            station_label << observer_id;
        }

        // =================================== 05 准备弧段与完整观测输出 ===================================
        fs::path all_filename;
        ofstream all_observations;

        if (output_all)
        {
            ostringstream name;
            name << "all_" << object_id << '_' << station_label.str() << '_'
                 << fixed << setprecision(10) << MJD.front() << ".txt";
            all_filename = output_dir / name.str();
            all_observations.open(all_filename);

            if (!all_observations) throw runtime_error("Cannot open complete observation file: " + all_filename.string());
        }

        ostringstream setup;
        setup << "\n================ Optical observation simulation ================\n"
              << "Input frame / time : GCRS (J2000 approximation) / MJD UTC\n"
              << "Input samples      : " << num << '\n'
              << "Epoch range [UTC]  : " << fixed << setprecision(10) << MJD.front() << " .. " << MJD.back() << '\n'
              << "Observation type   : 1 (ground-based GCRS RA/Dec)\n"
              << "Observer file / ID : " << fn_observer << " / " << observer_id << '\n'
              << "Minimum elevation  : " << min_elevation_deg << " deg\n"
              << "EOP                : interpolated at every observation epoch\n"
              << "Night / shadow     : " << use_station_night << " / " << use_earth_shadow << '\n'
              << "Sun ephemeris      : " << (need_sun ? (solar_ephemeris == 1 ? "1 = JPL / TDB" : "2 = analytical / UT1") : "unused; filters off") << '\n'
              << "Observation values : noise-free geometric angles\n"
              << "Truth output       : Kep=" << output_kep << ", Cart=" << output_cart << '\n'
              << "All-epoch output   : " << (output_all ? all_filename.string() : "OFF") << '\n'
              << "Output directory   : " << output_dir.string();

        if (use_station_night) setup << "\nSun elevation limit: " << max_sun_elevation_deg << " deg";

        std::cerr << setup.str() << '\n';

        vector<SimulatedObservationRecord> current_arc;
        vector<string> generated_obs_files;
        size_t arc_count = 0, observation_count = 0;
        size_t elevation_pass = 0, night_rejected = 0, shadow_rejected = 0;

        // =================================== 06 筛选、计算理论测角并分弧 ===================================

        // 按教材第 1.4 节筛选；可见点进入弧段，不可见点仅在完整输出时保留。
        int last_visible_id = -1;

        for (size_t k = 0; k < num; ++k)
        {
            const auto& cart = cartesian[k];
            const auto& station = site_positions[k];
            const double target_radius = Module_3D(cart[0], cart[1], cart[2]);

            // theta0 是测站与目标地心向量夹角的可见上限；输入长度为归一化单位，角度为 rad。
            const double max_central_angle = Visibility_Para_Theta0(
                station_radius, target_radius, min_elevation_deg * DD2R);
            bool visible = Visibility(station[0], station[1], station[2], cart[0], cart[1], cart[2], max_central_angle);

            if (visible) ++elevation_pass;

            if (visible && use_station_night)
            {
                // 式（1.74）：忽略太阳视差，用地心太阳方向与测站地心径向判断黑夜。
                if (!(sun_directions[k].dot(station) < station_radius * sin(max_sun_elevation_deg * DD2R)))
                {
                    ++night_rejected;
                    visible = false;
                }
            }

            if (visible && use_earth_shadow)
            {
                // 式（1.75）：球形地球、平行光柱形影；地球半径 ae 在此单位下为 1。
                // 影外要求 cos(psi) > -sqrt(1 - 1/r^2)，边界按教材严格不等式处理。
                const Vector3d target_position(cart[0], cart[1], cart[2]);

                if (target_radius <= 1.0 || !(sun_directions[k].dot(target_position) / target_radius > -sqrt(max(0.0, 1.0 - 1.0 / (target_radius * target_radius)))))
                {
                    ++shadow_rejected;
                    visible = false;
                }
            }

            if (!visible && !output_all) continue;

            // 遇到下一个可见点时，若其 ID 未接续上一可见点，则先结束上一弧段。
            if (visible && last_visible_id != -1 && static_cast<long long>(ID[k]) != static_cast<long long>(last_visible_id) + 1)
            {
                // 空弧不生成文件；原实现对单点、双点弧也写出，不增加最少记录数门限。
                if (!current_arc.empty())
                {
                    const string arc_file = WriteObservationArc(tracer, current_arc, output_dir, object_id, station_label.str(),
                                                               arc_count + 1, output_kep, output_cart);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

                    // 仅在全部文件成功关闭后推进弧号、登记路径，再清空缓存接收新弧。
                    ++arc_count;
                    generated_obs_files.push_back(arc_file);
                    current_arc.clear();
                }
            }

            // 完整输出还计算不可见点的几何方向；每条有效星历记录恰好输出一行。
            // 两位置同处地心惯性轴，位置已归一化；先作目标减测站的站心视线，再得距离及赤经赤纬(rad)。
            const Vector3d line_of_sight = Vector3d(cart[0], cart[1], cart[2]) - station;
            const EquatorialCoordinates equatorial = LOS_to_Equatorial(line_of_sight);
            SimulatedObservationRecord sample;
            sample.input_id = ID[k];
            observation::ObservationRecord& record = sample.observation;
            record.mjd_utc = MJD[k];
            record.type = observation_type;
            record.object_id = object_id;
            record.observer_id = observer_id;
            record.ra_deg = equatorial.right_ascension * DR2D;
            record.dec_deg = equatorial.declination * DR2D;

            // 理论角度直接写入 OBS；随机误差由 Python 在分弧结果上加入。
            if (!isfinite(record.ra_deg) || !isfinite(record.dec_deg))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Undefined theoretical RA/Dec at input row " + to_string(k + 1), EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            if (output_all)
            {
                // 前六列与正式 OBS 相同；之后追加可见标记和输入 ID。
                // 显式分隔编号，避免目标号较长时与相邻数值粘连。
                all_observations << setprecision(17) << record.mjd_utc << ' ' << record.type
                    << ' ' << record.object_id << ' ' << record.observer_id
                    << ' ' << record.ra_deg << ' ' << record.dec_deg
                    << ' ' << (visible ? 1 : 0) << ' ' << sample.input_id << '\n';
            }

            if (!visible) continue;

            sample.kep = keplerian[k];
            sample.cart = cart;

            for (int j = 0; j < 3; ++j) sample.cart[j] *= NORMALIZED_LENGTH_TO_KM;

            for (int j = 3; j < 6; ++j) sample.cart[j] *= NORMALIZED_VELOCITY_TO_KMPS;

            current_arc.push_back(sample);
            last_visible_id = ID[k];
            ++observation_count;
        }

        // 输入结束时写出尚未被后续不连续可见点触发的尾弧；与弧中断处使用同一文件格式。
        if (!current_arc.empty())
        {
            const string arc_file = WriteObservationArc(tracer, current_arc, output_dir, object_id, station_label.str(),
                                                       arc_count + 1, output_kep, output_cart);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 尾弧成功后同样先登记完成，再清空缓存；无可见点时不生成空文件。
            ++arc_count;
            generated_obs_files.push_back(arc_file);
            current_arc.clear();
        }

        if (output_all)
        {
            all_observations.close();

            if (!all_observations) throw runtime_error("Failed while writing complete observations: " + all_filename.string());
        }

        // =================================== 07 汇总与关闭日志 ===================================

        ostringstream completed;
        completed << "Observation simulation finished successfully. Arcs: " << arc_count << "; observations: " << observation_count
                  << "; elevation pass: " << elevation_pass << "; night rejected: " << night_rejected << "; shadow rejected: " << shadow_rejected
                  << "; all-epoch rows: " << (output_all ? num : 0)
                  << "; wall time: " << fixed << setprecision(3) << chrono::duration<double>(chrono::steady_clock::now() - start).count() << " s.";
        std::cerr << completed.str() << '\n';
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Optical observation simulation completed successfully.", EventType::MESSAGE, FORCE_ALERT_ON);
        cout << generated_obs_files.size() << '\n';

        for (const auto& filename : generated_obs_files) cout << filename << '\n';

        // 正常结束时集中保存事件；失败出口会写出相同事件和当时的调用链。
        if (program_log.is_open())
        {
            tracer.ReportEvents(program_log);
            program_log << "\nExit status: 0\n";
            program_log.close();

            if (!program_log) throw std::runtime_error("Failed while writing program log.");
        }

        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const exception& error)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR, FORCE_ALERT_ON);
        tracer.EmergencyExit(program_log);
    }

    return 1;
}
