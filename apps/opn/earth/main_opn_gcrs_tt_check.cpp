/**
 * @file main_opn_gcrs_tt_check.cpp
 * @brief 为程序间核对提供 GCRS/TT Cartesian 轨道预报、有效 GM 表和运行摘要。
 * @details 初值采用 UTC MJD、GCRS 位置 km 和速度 km/s；积分自变量为 TT，RKF78 定步长。
 * 地球质点和 GGM03C 球谐项共用 GM_TT=398600.4415 km³/s²、参考半径 6378.1363 km。
 * 日月及额外行星采用 constants::DE440::GM_TDB × (1-L_G)/(1-L_B)，单位换算后参与 TT 方程。
 * 月球球谐文件提供参考半径、归一化 C/S；本入口读取文件 GM 作记录，再将有效 GM 覆盖为 DE440 月球 TT 值。
 * 此处保留原 C/S，因此改变球谐加速度幅度，属于本次双方核对采用的模型约定。
 * 归一化以地球 GM、半径定义单位；第三体 GM/GM_E 在右函数中使用，输出时恢复 km、km/s。
 * 地球 GM 的取值来自 GGM03C/IERS，日月来自 DE440；实际地月 GM 比由有效值相除，另列 JPL EMRAT。
 *
 * 历表使用 DE440；JPL 查询暂以 TT 近似 TDB，直接采用地心差矢量；积分包含牛顿引力与所选球谐项。
 * 本入口的相对论处理限于 GM 尺度约定；完整的 1PN 动力学和 BCRS→GCRS 状态变换见文献 [1]。
 * [1] Turyshev, Bar-Sever & Bertiger (2026), Astronomical Journal, 172(1), 58.
 *     DOI: 10.3847/1538-3881/ae7436；II.4 节讨论尺度，III.2 节讨论地心动力学。
 * [2] Klioner (2008), Astronomy & Astrophysics, 478(3), 951–958.
 *     Relativistic scaling of astronomical quantities and the system of astronomical units.
 *     DOI: 10.1051/0004-6361:20077786；时间、坐标和 GM 的尺度变换须配套。
 * [3] Park et al. (2021), Astronomical Journal, 161, 105. DOI: 10.3847/1538-3881/abd414.
 *     DE440/DE441 的 TDB 相容 GM 见表 2；源码集中定义在 utils_constants.hpp。
 *
 * 使用方法、Windows SYSTEM 后端构建命令及本算例的具体约定见 doc/18_gcrs-tt-orbit-check.md。
 */
#include "ASTER.h"
#include "io_settings.hpp"
#include "math_backend.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

using namespace bastion;
using namespace constants::Scale;

/**
 * @brief 读取对比设置，以固定 TT 步长传播六维 Cartesian 状态并保存 CSV 和摘要。
 * @param[in] argc 命令行参数个数。
 * @param[in] argv 参数数组；支持 -I setting.txt、--help、--version。
 * @return 0 表示完成；配置、模型、积分或输出错误通过 Overwatch 记录并以非零状态退出。
 */
int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_opn_gcrs_tt_check " << ASTER_PROGRAM_VERSION << '\n'
                  << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--help"))
    {
        std::cout << "Usage: aster_opn_gcrs_tt_check -I setting.txt\n"
                  << "GCRS/TT Cartesian; DE440; RKF78 fixed step; SYSTEM math backend.\n";
        return argc == 1 ? 2 : 0;
    }

    Overwatch tracer(__FILE__, __func__, __LINE__);
    std::ofstream program_log;

    try
    {
        // ==================== 1. 读取任务配置与输出位置 ====================
        if (argc != 3 || std::string(argv[1]) != "-I")
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected -I setting.txt.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 本交付入口要求工具链原生数学库，Windows 与本机测试使用相同选项。
#if ASTER_MATH_OPENLIBM
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Configure this comparison with -DASTER_MATH_BACKEND=SYSTEM.", EventType::ERROR);
        tracer.EmergencyExit(program_log);
#endif
        namespace fs = std::filesystem;
        const fs::path setting_path = fs::absolute(argv[2]);
        const fs::path base = setting_path.parent_path();
        std::string setting_file = setting_path.string();
        io::SimpleInput input;
        input.ReadSimpleInput(tracer, setting_file, true);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 沿用原 Cartesian 主程序的关键词；各字段的单位、开关和路径含义保持原定义。
        const std::set<std::string> keywords = {
            "reference_epoch", "reference_Keplerian", "reference_Cartesian", "propagation_length",
            "file_EOP", "file_JPLEPH", "file_atmosphere_density", "file_lunar_gravity_field",
            "earth_zonal_order", "earth_tesseral_degree", "solar_gravity", "lunar_gravity",
            "use_lunar_NSG", "lunar_zonal_order", "lunar_tesseral_degree", "use_drag_perturb",
            "use_SRP_perturb", "use_solidtide_perturb", "satellite_id", "satellite_mass",
            "satellite_drag_cross_area", "satellite_srp_cross_area", "satellite_drag_Cd", "satellite_srp_eta",
            "integrator_id", "integrator_adams_steps", "integrator_stepsize", "integrator_maxTE",
            "file_eph_Cartesian", "file_eph_Keplerian", "file_eph_NonSingular2", "file_report",
            "file_programlog", "output_state_units", "display_model", "display_setting", "display_progress"
        };

        for (size_t i = 0; i < input.Count(); ++i)
        {
            io::SimpleRecord record;
            input.GetRecord(tracer, i + 1, record);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            if (keywords.count(record.kw) == 0)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown keyword: " + record.kw, EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 输出文件分别采用原有 file_* 关键词，相对路径以 setting 目录为基准。
        const std::string cartesian_name = input.GetString(tracer, "file_eph_Cartesian");
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        if (cartesian_name.empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "file_eph_Cartesian must be nonempty for the Cartesian comparison.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        const fs::path cartesian_path = base / cartesian_name;
        fs::create_directories(cartesian_path.parent_path());
        std::string report_name, log_name;

        if (input.HasKeyword("file_report")) input.GetString(tracer, "file_report", report_name);

        if (input.HasKeyword("file_programlog")) input.GetString(tracer, "file_programlog", log_name);

        if (!log_name.empty())
        {
            fs::create_directories((base / log_name).parent_path());
            program_log.open(base / log_name);

            if (!program_log)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open file_programlog: " + log_name, EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 原始配置随日志保存，核对实际输入及模型开关。
            std::ifstream setting_stream(setting_path);
            program_log << "Input settings:\n" << setting_stream.rdbuf() << "\nEvents:\n";
        }

        OrbitPropagationJob job;
        double epoch_utc, duration_s, step_s, max_te_setting;
        double state_km[6];
        const std::string epoch_text = input.GetString(tracer, "reference_epoch");
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // reference_epoch 保留 UTC MJD、JD 和六项年月日时分秒三种原有形式。
        std::vector<std::string> epoch_parts;
        split(epoch_text, epoch_parts);

        if (epoch_parts.size() == 1)
        {
            input.GetDoubleStrict(tracer, "reference_epoch", epoch_utc);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            if (epoch_utc > constants::MJD_OFFSET) epoch_utc -= constants::MJD_OFFSET;
        }
        else if (epoch_parts.size() == 6)
        {
            int year, month, day, hour, minute;
            double second;
            std::istringstream epoch_stream(epoch_text);
            epoch_stream >> year >> month >> day >> hour >> minute >> second;

            if (!epoch_stream || !std::isfinite(second) || hour < 0 || hour > 23 || minute < 0 || minute > 59 || second < 0 || second >= 60)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid reference_epoch calendar fields.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            const int status = Calendar2MJD(year, month, day, hour, minute, second, epoch_utc);
            if (status != 0)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid reference_epoch calendar date.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }
        }
        else
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "reference_epoch requires UTC MJD, JD, or six calendar fields.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 原 reference_Cartesian 一行接收 x y z vx vy vz，单位 km、km/s。
        const std::string state_text = input.GetString(tracer, "reference_Cartesian");
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        std::istringstream state_stream(state_text);

        for (double& component : state_km)
        {
            state_stream >> component;

            if (!state_stream || !std::isfinite(component))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "reference_Cartesian requires six finite values in km and km/s.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        std::string extra_state;
        state_stream >> extra_state;

        std::string keplerian_text;

        if (input.HasKeyword("reference_Keplerian"))
        {
            input.GetString(tracer, "reference_Keplerian", keplerian_text);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        if (!extra_state.empty() || !keplerian_text.empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Use six reference_Cartesian values and leave reference_Keplerian empty for this comparison.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 原积分设置仍以秒给出；102 为 RKF78 定步长，maxTE 为归一化局部误差告警阈值。
        input.GetDoubleStrict(tracer, "propagation_length", duration_s);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetDoubleStrict(tracer, "integrator_stepsize", step_s);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetDoubleStrict(tracer, "integrator_maxTE", max_te_setting);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 原有地月带谐与田谐分别设置；use_lunar_NSG 独立控制月球非球形项。
        int earth_order, earth_degree, moon_order, moon_degree;
        int solar_gravity, lunar_gravity, use_lunar_nsg, integrator_id, output_units;
        int use_drag, use_srp, display_model, display_setting, display_progress;
        input.GetIntegerStrict(tracer, "earth_zonal_order", earth_order);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "earth_tesseral_degree", earth_degree);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "lunar_zonal_order", moon_order);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "lunar_tesseral_degree", moon_degree);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "solar_gravity", solar_gravity);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "lunar_gravity", lunar_gravity);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "use_lunar_NSG", use_lunar_nsg);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "integrator_id", integrator_id);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "output_state_units", output_units);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "use_drag_perturb", use_drag);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "use_SRP_perturb", use_srp);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "display_model", display_model);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "display_setting", display_setting);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        input.GetIntegerStrict(tracer, "display_progress", display_progress);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        int use_solid_tide = 0;

        if (input.HasKeyword("use_solidtide_perturb"))
        {
            input.GetIntegerStrict(tracer, "use_solidtide_perturb", use_solid_tide);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 核对入口运行本次约定的定步长、日月 DE440 引力模型，检查误设的其他模型开关。
        if (duration_s < 0 || step_s <= 0 || duration_s / step_s > 1e8 || max_te_setting <= 0
            || earth_order < 0 || earth_order > 80 || earth_degree < 0 || earth_degree > 80
            || moon_order < 0 || moon_order > 80 || moon_degree < 0 || moon_degree > 80
            || integrator_id != 102 || output_units != 2 || use_drag != 0 || use_srp != 0 || use_solid_tide != 0
            || solar_gravity < 0 || solar_gravity > 1 || lunar_gravity < 0 || lunar_gravity > 1
            || use_lunar_nsg < 0 || use_lunar_nsg > 1 || (use_lunar_nsg && lunar_gravity != 1)
            || display_model < 0 || display_model > 1 || display_setting < 0 || display_setting > 1
            || display_progress < 0 || display_progress > 2
            || std::hypot(state_km[0], state_km[1], state_km[2]) == 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Check propagation_length>=0, integrator_stepsize>0, maxTE>0, gravity degrees 0..80, integrator_id=102, output_state_units=2, solar/lunar_gravity=0 or 1, and disabled drag/SRP/tides.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // ==================== 2. 历表、EOP 与力模型 ====================
        const double epoch_tt = UTC2TT(epoch_utc);
        const double end_tt = epoch_tt + duration_s / constants::SECONDS_PER_DAY;
        const std::string jpl_name = input.GetString(tracer, "file_JPLEPH");
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        const int jpl_status = Read_JPLEPH((base / jpl_name).string());
        if (jpl_name.empty() || jpl_status <= 0 || JPLEPH_R1.r1.numde != 440
            || epoch_tt < JPLEPH_R1.r1.ss[0] - constants::MJD_OFFSET
            || end_tt > JPLEPH_R1.r1.ss[1] - constants::MJD_OFFSET)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Require native binary DE440 covering the complete TT interval.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 有资料时按日插值；超出覆盖范围时定向量归零，DAT 由 SOFA 单独取得。
        job.opn_EOP.SetZeroFallback(true);
        const std::string eop_name = input.GetString(tracer, "file_EOP");
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        if (!eop_name.empty())
        {
            std::string eop_path = (base / eop_name).string();
            job.opn_EOP.Read_EOP(tracer, eop_path);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        job.opn_model.use_earth_NSG = earth_order > 0 || earth_degree > 0;
        job.opn_model.earth_zonal_order = earth_order;
        job.opn_model.earth_tesseral_degree = earth_degree;
        job.opn_model.solar_gravity_config = solar_gravity;
        job.opn_model.lunar_gravity_config = lunar_gravity;
        job.opn_model.use_lunar_NSG = use_lunar_nsg != 0;
        job.opn_model.lunar_zonal_order = moon_order;
        job.opn_model.lunar_tesseral_degree = moon_degree;

        // 月球场的参考半径与 C/S 来自文件；本次核对将 GM 显式覆盖为与月球质点相同的 TT 值。
        double lunar_file_gm = 0;
        const double moon_gm_tt = constants::DE440::GM_TDB[9] * constants::TDB_TO_TT_SCALE;
        const std::string moon_name = input.GetString(tracer, "file_lunar_gravity_field");
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        if (use_lunar_nsg)
        {
            const bool loaded = moon_gravity_field.Load((base / moon_name).string());
            if (!loaded || std::max(moon_order, moon_degree) > moon_gravity_field.max_degree)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot load lunar field at requested degree.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            lunar_file_gm = moon_gravity_field.GM;
            moon_gravity_field.GM = moon_gm_tt;
        }

        // ==================== 3. 打印实际 GM、初态与积分约定 ====================
        // 地球 GM 同时定义归一化时间，恢复有量纲加速度后仍是该物理 GM。
        // 全部第三体统一乘同一个 TDB→TT 因子；球谐 C/S 保持无量纲原值。
        const std::array<const char*, 11> body_names = {
            "Mercury", "Venus", "Earth", "Mars_system", "Jupiter_system", "Saturn_system",
            "Uranus_system", "Neptune_system", "Pluto_system", "Moon", "Sun"
        };
        std::ostringstream summary;
        summary << std::setprecision(17)
                << "aster_opn_gcrs_tt_check " << ASTER_PROGRAM_VERSION << '\n'
                << aster::math_backend::BuildDescription() << '\n'
                << "Frame=GCRS; independent_time=TT; ephemeris=DE440\n"
                << "JPL_argument=TT_approx_TDB; vectors=geocentric_ephemeris_difference; 1PN=off\n"
                << "Integrator=RKF78_fixed; step_s=" << step_s << "; duration_s=" << duration_s << '\n'
                << "epoch_UTC_MJD=" << epoch_utc << "; epoch_TT_MJD=" << epoch_tt << '\n'
                << "TT_minus_UTC_s=" << (epoch_tt - epoch_utc) * constants::SECONDS_PER_DAY << '\n'
                << "EOP=interpolate_if_covered_else_zero_orientation; DAT=SOFA_leap_table\n"
                << "Earth_order_degree=" << earth_order << "," << earth_degree << "; Moon_NSG=" << use_lunar_nsg << "; Moon_order_degree=" << moon_order << "," << moon_degree << '\n'
                << "Sun=" << solar_gravity << "; Moon=" << lunar_gravity << "; other_planets=0" << '\n'
                << "SRP=off; drag=off; tides=off\n"
                << "Earth_GM_TT_km3_s2=" << UNIT_GM_M3_S2 / 1e9 << '\n'
                << "Earth_radius_m=" << constants::GGM03C::REFERENCE_RADIUS << '\n'
                << "TDB_to_TT_GM_factor=" << constants::TDB_TO_TT_SCALE << '\n'
                << "Earth_over_Moon_effective=" << UNIT_GM_M3_S2 / moon_gm_tt << '\n'
                << "Moon_over_Earth_effective=" << moon_gm_tt / UNIT_GM_M3_S2 << '\n'
                << "DE440_header_EMRAT=" << JPLEPH_R1.r1.emrat << '\n'
                << "GM table: body,DE440_TDB_km3_s2,DE440_scaled_TT_km3_s2\n";

        for (int body = 0; body < 11; ++body)
        {
            summary << body_names[body] << ',' << constants::DE440::GM_TDB[body] / 1e9 << ','
                    << constants::DE440::GM_TDB[body] * constants::TDB_TO_TT_SCALE / 1e9 << '\n';
        }

        summary << "Earth force uses Earth_GM_TT above; DE440 Earth row is for comparison.\n";

        if (use_lunar_nsg)
        {
            summary << "Moon_field_file_GM_km3_s2=" << lunar_file_gm / 1e9 << '\n'
                    << "Moon_field_effective_GM_TT_km3_s2=" << moon_gravity_field.GM / 1e9 << '\n'
                    << "Moon_field_radius_m=" << moon_gravity_field.R_ref << '\n'
                    << "Moon_field_amplitude_factor=" << moon_gravity_field.GM / lunar_file_gm << '\n';
        }

        summary << "Initial_cartesian_km_km_s=";

        for (const double component : state_km) summary << ' ' << component;

        summary << '\n';

        if (display_model || display_setting) std::cout << summary.str() << std::flush;

        std::ofstream report;

        if (!report_name.empty())
        {
            fs::create_directories((base / report_name).parent_path());
            report.open(base / report_name);
        }

        std::ofstream ephemeris(cartesian_path);

        if ((!report_name.empty() && !report) || !ephemeris)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot create file_report or file_eph_Cartesian.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        if (report.is_open()) report << summary.str();

        ephemeris << std::setprecision(17) << "elapsed_TT_s,MJD_UTC,MJD_TT,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s\n";

        // ==================== 4. 定步长积分与 Cartesian 输出 ====================
        job.opn_setting.t0_MJD_UTC = epoch_utc;
        job.opn_setting.t0_MJD_TT = epoch_tt;
        job.opn_setting.rkf.fixed_step = true;
        job.opn_setting.rkf.maxTE_norm = max_te_setting;
        job.opn_setting.display_progress = display_progress == 2 ? 2 : 0;
        double state_norm[6];

        for (int component = 0; component < 6; ++component)
        {
            const double scale = component < 3 ? 1000.0 / UNIT_LENGTH_M : 1000.0 * UNIT_TIME_S / UNIT_LENGTH_M;
            state_norm[component] = state_km[component] * scale;
        }

        double elapsed_s = 0, max_te = 0;
        size_t rows = 0;
        int next_progress = 10;

        while (true)
        {
            // 先写初态或本次接受步的状态。历元从起点和累计秒数计算，避免 MJD 反复相加。
            const double current_tt = epoch_tt + elapsed_s / constants::SECONDS_PER_DAY;
            ephemeris << elapsed_s << ',' << TT2UTC(current_tt) << ',' << current_tt;

            for (int component = 0; component < 6; ++component)
            {
                const double scale = component < 3 ? UNIT_LENGTH_M / 1000.0 : UNIT_LENGTH_M / UNIT_TIME_S / 1000.0;
                ephemeris << ',' << state_norm[component] * scale;
            }

            ephemeris << '\n';
            ++rows;

            if (elapsed_s >= duration_s) break;

            // RKF78 的 13 级右函数始终使用 TT；定步长模式仅记录嵌入误差估计。
            const double next_elapsed = std::min(duration_s, static_cast<double>(rows) * step_s);
            const double h_norm = (next_elapsed - elapsed_s) / UNIT_TIME_S;
            double accepted_step, accepted_te;
            RKF78_Cart(tracer, job, current_tt, state_norm, h_norm, accepted_step, accepted_te);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            if (accepted_step != h_norm || !std::isfinite(accepted_te))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Fixed-step RKF78 returned an invalid accepted step or error estimate.", EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            elapsed_s = next_elapsed;
            max_te = std::max(max_te, accepted_te);

            if (display_progress > 0 && 100.0 * elapsed_s / duration_s >= next_progress)
            {
                std::cout << "Propagation " << static_cast<int>(100.0 * elapsed_s / duration_s) << "%; elapsed_TT_s=" << elapsed_s << std::endl;
                next_progress += 10;
            }
        }

        // ==================== 5. 保存末状态、统计量与运行事件 ====================
        std::ostringstream results;
        results << std::setprecision(17) << "rows=" << rows << '\n'
                << "final_TT_MJD=" << end_tt << "; final_UTC_MJD=" << TT2UTC(end_tt) << '\n'
                << "max_local_TE_normalized=" << max_te << '\n'
                << "zero_EOP_queries=" << job.opn_EOP.ZeroFallbackCount() << '\n'
                << "Final_cartesian_km_km_s=";

        for (int component = 0; component < 6; ++component)
        {
            const double scale = component < 3 ? UNIT_LENGTH_M / 1000.0 : UNIT_LENGTH_M / UNIT_TIME_S / 1000.0;
            results << ' ' << state_norm[component] * scale;
        }

        results << '\n';

        if (report.is_open())
        {
            report << results.str();
            report.close();
        }

        ephemeris.close();

        if (!ephemeris || (!report_name.empty() && !report))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing propagation output.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        tracer.CheckOut(__func__);

        if (program_log.is_open())
        {
            tracer.ReportEvents(program_log);
            program_log.close();
        }

        std::cout << results.str() << "Output: " << cartesian_path.lexically_normal().string() << '\n';
    }
    catch (const std::exception& error)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
        tracer.EmergencyExit(program_log);
    }

    return 0;
}
