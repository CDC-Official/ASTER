/**
 * @file main_opn_cislunar.cpp
 * @brief 在 E-BCRS/TDB 中传播 Cartesian，输出同坐标系的轨道星历、任务摘要和运行事件。
 * @details 输入状态直接解释为 E-BCRS，积分和输出沿用该坐标系，不增加 E-BCRS 与 GCRS 的坐标转换。
 * 主程序依次读取设置、准备 EOP 与力模型、显示任务参数、调用积分并写出结果。
 * 输入与输出历元均为 TDB MJD；状态计算采用 ae、UNITTIME 归一化，控制台显示 km、km/s。
 * 控制台初末状态只显示 Cartesian，两行分别为位置和速度；星历单位独立由配置控制。
 * 相对路径以 setting 所在目录为基准；tracer 记录文件读写路径、警告和错误。
 * file_report 保存任务设置与结果，file_programlog 保存原始设置、事件及失败调用链，均可单独关闭。
 * 用法：aster_opn_cislunar -I <配置>。地月引力、行星引力和后牛顿修正按现有配置计算。
 */
#include "math_backend.hpp"
#include "io_settings.hpp"
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ctime>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

#include "ASTER.h"

using namespace std;
using namespace Eigen;
using namespace orbit;
using namespace bastion;
using orbit::Keplerian;
using orbit::Cartesian;
using namespace io;

int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_opn_cislunar " << ASTER_PROGRAM_VERSION << '\n';
        std::cout << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_opn_cislunar -I <setting.txt>\n";
        return 2;
    }

	Overwatch tracer(__FILE__, __func__, __LINE__);
	std::ofstream program_log;

	try
	{
		// 从进入任务开始计实际耗时；传播调用本身另行计时。
		const auto wall_start = std::chrono::steady_clock::now();

		// 说明任务信息
		std::time_t timeNow = std::time(nullptr);    // 获取当前时间
		std::tm* tmNow = std::localtime(&timeNow);
		char timeStr[20];    // 创建用于存储时间格式化结果的字符串
		char jobId[15];    // 任务编号
		std::strftime(jobId, sizeof(jobId), "%Y%m%d%H%M%S", tmNow);// 生成 jobId 格式：YYYYMMDDHHMMSS
		std::strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", tmNow);    // 生成日志消息格式：YYYY-MM-DD HH:MM:SS
		std::stringstream jobMsg;    // 构建日志消息
		jobMsg << "E-BCRS orbit propagation job started at " << timeStr << " (jobId = " << jobId
			   << ")";
		tracer.LogEvent(__FILE__, __func__, __LINE__, jobMsg.str(), bastion::EventType::MESSAGE, FORCE_ALERT_ON);

		string msg;

		namespace fs = std::filesystem;

		// 程序名取实际执行文件名，版本由 CMake 统一提供。
		const string program_name = fs::path(argv[0]).filename().string();
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Program: " + program_name + " | Version: v" + ASTER_PROGRAM_VERSION, EventType::MESSAGE, FORCE_ALERT_ON);

		// 默认读取随入口提供的设置，命令行可指定其他配置文件。
		std::string fn_EBCRS_input = "";

		// 按顺序处理参数；缺值和未知选项均报错，避免误用默认设置。
		for (int i = 1; i < argc; ++i)
		{
			const std::string arg = argv[i];
			if (arg == "-h" || arg == "--help")
			{
				std::cout << "Usage: aster_opn_cislunar [-I|-input|--input <setting>]\n"
					<< "Relative paths use the setting directory. Existing outputs are never overwritten.\n";
				tracer.CheckOut(__func__);
				return 0;
			}
			// 三种输入选项都直接改写同一个文件名；后出现的设置优先。
			if (arg != "-I" && arg != "-input" && arg != "--input")
				throw std::runtime_error("Unknown argument: " + arg);
			// 输入选项后必须紧跟路径，不能把下一个选项误当成文件名。
			if (i + 1 >= argc || std::string(argv[i + 1]).empty() || argv[i + 1][0] == '-')
				throw std::runtime_error("Missing value for " + arg);
			fn_EBCRS_input = argv[++i];
		}

		/* MARK: --------------------------------- read input ------------------------------------ */

		msg = "Reading and analyzing input from: " + fs::absolute(fn_EBCRS_input).lexically_normal().string();
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		SimpleInput si_EBCRSInput;
		OrbitPropagationNModel_EBCRS EBCRS_model;
		OrbitPropagationSetting_EBCRS EBCRS_setting;
		OrbitPropagationJob_EBCRS EBCRS_job;

		si_EBCRSInput.ReadSimpleInput(tracer, fn_EBCRS_input);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_programlog：保存原始配置、运行事件与失败调用链；省略或留空时关闭。
		// 相对路径以 setting 所在目录为基准，与 file_report 的任务报告分别输出。
		std::string fn_programlog;

		if (si_EBCRSInput.HasKeyword("file_programlog"))
		{
			si_EBCRSInput.GetString(tracer, "file_programlog", fn_programlog);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		if (!fn_programlog.empty())
		{
		    const std::filesystem::path log_config = std::filesystem::weakly_canonical(std::filesystem::absolute(fn_EBCRS_input));
		    std::filesystem::path log_path = fn_programlog;
		    if (log_path.is_relative()) log_path = log_config.parent_path() / log_path;
		    log_path = std::filesystem::weakly_canonical(log_path);
		    fn_programlog = log_path.string();
		    // 日志不能覆盖已有资料，也不能与本次计划写出的其他文件重名。
		    if (std::filesystem::exists(log_path))
		        throw std::runtime_error("Program log already exists: " + log_path.string());
		    for (std::size_t record_index = 1; record_index <= si_EBCRSInput.Count(); ++record_index)
		    {
		        io::SimpleRecord path_record;
		        si_EBCRSInput.GetRecord(tracer, record_index, path_record);
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
		    for (std::size_t record_index = 1; record_index <= si_EBCRSInput.Count(); ++record_index)
		    {
		        io::SimpleRecord path_record;
		        si_EBCRSInput.GetRecord(tracer, record_index, path_record);
		        if (path_record.kw.compare(0, 5, "file_") != 0 || path_record.va.empty()) continue;
		        std::filesystem::path recorded_path = path_record.va;
		        if (recorded_path.is_relative()) recorded_path = log_config.parent_path() / recorded_path;
		        program_log << path_record.kw << " = "
		                    << std::filesystem::weakly_canonical(recorded_path).string() << '\n';
		    }
		    program_log.flush();
		    if (!program_log) throw std::runtime_error("Failed while writing program log.");
		}

		// 设置读取成功后确定资料基准目录，符号链接沿用实际设置文件的位置。
		const fs::path config_path = fs::weakly_canonical(fs::absolute(fn_EBCRS_input));

		// MARK: >>>>>>>>>>>>>>>>>>>>>>>> compulsory keywords
		double t0_mjd_tdb = 0.0, propagation_length_sec = 0.0, integrator_step_sec = 0.0, integrator_maxTE = 0.0;
		string str_epoch, str_cart, str_kep, str_planetary_body;
		string fn_EOP, fn_JPL, fn_ephCart, fn_report, fn_atm_density, fn_lunar_gravity_field;
		string sat_id;
		double sat_mass = 0.0, sat_area_drag = 0.0, sat_area_SRP = 0.0, sat_eta_SRP = 0.0, sat_Cd_drag = 0.0;
		int earth_zonal_order = 0, earth_tesseral_degree = 0, solar_gravity_config = 0, lunar_gravity_config = 0;
		int use_drag = 0, use_SRP = 0, use_post_newtonian = 0, use_planetary_gravity = 0;
		int use_lunar_NSG = 0, lunar_zonal_order = 0, lunar_tesseral_degree = 0;
		int eph_normalize = STATE_UNIT_KM_SEC_DEGREE, display_summary = 0;
		int display_progress = 1; // 0 静默；1 约每 10%；2 每步及试算误差。
		int integrator_id = 101, adams_steps = 9;
		int status = 0;
		bool isKepReady = false, isCartReady = false;
		vector<int> planetary_body_id;
		vector<AtmosphereData> atm_list;
		Keplerian initial_kep;
		Cartesian initial_cart;

		// initial epoch, interpreted as TDB.
		si_EBCRSInput.GetString(tracer, "reference_epoch", str_epoch);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (str_epoch.empty())
		{
			msg = "Parameter reference_epoch is missing!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		vector<string> vec_epoch;
		split(str_epoch, vec_epoch);
		switch (vec_epoch.size())
		{
			case 1:
				// 输入为 JD 或 MJD，统一转成 MJD(TDB) 存储。
				t0_mjd_tdb = stod(vec_epoch[0]) > 2400000.5 ? stod(vec_epoch[0]) - 2400000.5 : stod(vec_epoch[0]);
				break;
			case 6:
				// 输入为年月日时分秒，解释为 TDB 日历时间。
				Calendar2MJD(stoi(vec_epoch[0]), stoi(vec_epoch[1]), stoi(vec_epoch[2]), stoi(vec_epoch[3]), stoi(vec_epoch[4]), stod(vec_epoch[5]), t0_mjd_tdb);
				break;
			default:
				msg = "Unknown format for parameter reference_epoch! Use JD, MJD, or yy mm dd hh mm ss.s in TDB.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
				break;
		}

		// try initial Keplerian
		si_EBCRSInput.GetString(tracer, "reference_Keplerian", str_kep);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!str_kep.empty())
		{
			vector<string> vec_kep;
			split(str_kep, vec_kep);
			if (vec_kep.size() != 6)
			{
				msg = "Invalid format of parameter reference_Keplerian! Use a e i Omega omega M in km and deg.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 输入为 EBCRS 下、以地球为中心定义的 Keplerian osculating elements：半长轴 km，角度 deg。
			initial_kep.Set(stod(vec_kep[0]), stod(vec_kep[1]), stod(vec_kep[2]), stod(vec_kep[3]), stod(vec_kep[4]), stod(vec_kep[5]), false);
			initial_kep.Normalize(tracer, ae / 1000.0);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 后续地月空间积分入口仍以 Cartesian 为主；Keplerian 只作为输入便利格式。
			initial_kep.ToCartesian(tracer, initial_cart);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 根数转换函数以归一化 GM=1 计算速度；E-BCRS 初速度按 DE440 地球 GM 缩放。
			double state_norm[6];
			initial_cart.ToArray(state_norm);
			const double velocity_scale = sqrt(MAJOR_BODY_GM_NORMALIZED[2]);

			for (int i = 3; i < 6; ++i) state_norm[i] *= velocity_scale;

			initial_cart.Set(state_norm[0], state_norm[1], state_norm[2], state_norm[3], state_norm[4], state_norm[5], true);
			isKepReady = true;
		}

		// try initial Cartesian
		si_EBCRSInput.GetString(tracer, "reference_Cartesian", str_cart);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!str_cart.empty())
		{
			vector<string> vec_cart;
			split(str_cart, vec_cart);
			if (vec_cart.size() != 6)
			{
				msg = "Invalid format of parameter reference_Cartesian! Use x y z vx vy vz in km and km/s.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 输入为 EBCRS Cartesian 状态：位置 km，速度 km/s。
			initial_cart.Set(stod(vec_cart[0]), stod(vec_cart[1]), stod(vec_cart[2]), stod(vec_cart[3]), stod(vec_cart[4]), stod(vec_cart[5]), false);
			initial_cart.Normalize(tracer, ae / 1000.0, UNITTIME);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
			isCartReady = true;
		}

		// 指定初始状态量：Cartesian 最常用，但 Keplerian 也支持；两者必须且只能存在一个。
		if (isCartReady && isKepReady)
		{
			msg = "Both Cartesian and Keplerian are fed as initial states. Select only one type!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		else if (!isCartReady && !isKepReady)
		{
			msg = "Neither Cartesian nor Keplerian is fed as initial states!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		if (isCartReady)
		{
			// 用 DE440 地球 GM 求参考根数，先将副本的速度换成转换函数要求的 GM=1 尺度。
			double state_norm[6];
			initial_cart.ToArray(state_norm);
			const double velocity_scale = sqrt(MAJOR_BODY_GM_NORMALIZED[2]);

			for (int i = 3; i < 6; ++i) state_norm[i] /= velocity_scale;

			Cartesian state_for_elements;
			state_for_elements.Set(state_norm[0], state_norm[1], state_norm[2], state_norm[3], state_norm[4], state_norm[5], true);
			state_for_elements.ToKeplerian(tracer, initial_kep);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// propagation length
		si_EBCRSInput.GetDouble(tracer, "propagation_length", propagation_length_sec);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		// 与地球轨道入口一致：正时长前向，负时长反向，零时长只输出初态。
		if (!std::isfinite(propagation_length_sec))
		{
			msg = "Invalid propagation_length. It must be finite seconds.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// EOP file
		si_EBCRSInput.GetString(tracer, "file_EOP", fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// JPLEPH file
		si_EBCRSInput.GetString(tracer, "file_JPLEPH", fn_JPL);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (fn_JPL.empty())
		{
			msg = "Parameter file_JPLEPH is missing. EBCRS propagation needs a JPL ephemeris.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// Atmosphere density file. It is only required when atmospheric drag is enabled.
		si_EBCRSInput.GetString(tracer, "file_atmosphere_density", fn_atm_density);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 与地球轨道入口使用相同编号；自适应 RKF78 统一使用 101。
		si_EBCRSInput.GetInteger(tracer, "integrator_id", integrator_id);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (integrator_id != 101 && integrator_id != 102 && integrator_id != 201 && integrator_id != 202)
			throw std::runtime_error("Supported E-BCRS integrator_id values: 101=adaptive RKF78, 102=fixed-step RKF78, 201=Adams PE, 202=Adams PECE.");
		if (integrator_id == 201 || integrator_id == 202)
		{
			si_EBCRSInput.GetInteger(tracer, "integrator_adams_steps", adams_steps);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
			if (adams_steps < 2 || adams_steps > 9)
				throw std::runtime_error("integrator_adams_steps must be an integer in 2..9.");
		}

		// integrator stepsize
		si_EBCRSInput.GetDouble(tracer, "integrator_stepsize", integrator_step_sec);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!std::isfinite(integrator_step_sec) || integrator_step_sec <= 0.0)
		{
			msg = "Invalid integrator_stepsize. It must be positive seconds.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// integrator maxTE
		si_EBCRSInput.GetDouble(tracer, "integrator_maxTE", integrator_maxTE);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!std::isfinite(integrator_maxTE) || integrator_maxTE <= 0.0)
		{
			msg = "Invalid integrator_maxTE. It must be positive.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// Earth gravity order and degree
		si_EBCRSInput.GetInteger(tracer, "earth_zonal_order", earth_zonal_order);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetInteger(tracer, "earth_tesseral_degree", earth_tesseral_degree);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// Solar and lunar point-mass gravity. EBCRS propagation should use JPL.
		si_EBCRSInput.GetInteger(tracer, "solar_gravity", solar_gravity_config);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetInteger(tracer, "lunar_gravity", lunar_gravity_config);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (solar_gravity_config != SOLAR_GRAVITY_OFF && solar_gravity_config != SOLAR_GRAVITY_JPL)
		{
			msg = "Invalid solar_gravity for EBCRS propagation. Use 0 (off) or 1 (JPL); analytical ephemeris is not supported.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		if (lunar_gravity_config != LUNAR_GRAVITY_OFF && lunar_gravity_config != LUNAR_GRAVITY_JPL)
		{
			msg = "Invalid lunar_gravity for EBCRS propagation. Use 0 (off) or 1 (JPL); analytical ephemeris is not supported.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// Planetary gravity: planetary_gravity is the global switch, planetary_body lists selected JPL body IDs.
		si_EBCRSInput.GetInteger(tracer, "planetary_gravity", use_planetary_gravity);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetString(tracer, "planetary_body", str_planetary_body);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		vector<string> vec_planetary_body;
		split(str_planetary_body, vec_planetary_body);
		planetary_body_id.clear();
		for (const string& body_token : vec_planetary_body)
		{
			int body_id = 0;
			try
			{
				size_t parsed_pos = 0;
				body_id = stoi(body_token, &parsed_pos);
				if (parsed_pos != body_token.size())
				{
					msg = "Invalid planetary_body entry: " + body_token + ". Use pure integer JPL body IDs separated by spaces.";
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}
			}
			catch (...)
			{
				msg = "Invalid planetary_body entry: " + body_token + ". Use JPL integer body IDs separated by spaces.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			if (body_id == 3)
			{
				msg = "Invalid planetary_body entry: BODY_EARTH(3). Earth is the EBCRS origin and cannot be used as a third-body perturbation.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			if (body_id == 10 || body_id == 11)
			{
				msg = "Invalid planetary_body entry: BODY_MOON(10) and BODY_SUN(11) are controlled by lunar_gravity and solar_gravity.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			if (body_id < 1 || body_id > 9)
			{
				msg = "Invalid planetary_body entry: " + body_token + ". Allowed IDs are 1, 2, 4, 5, 6, 7, 8, and 9.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			if (find(planetary_body_id.begin(), planetary_body_id.end(), body_id) != planetary_body_id.end())
			{
				msg = "Duplicate planetary_body entry: " + body_token + ". Each perturbing planet should appear only once.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			planetary_body_id.push_back(body_id);
		}
		if (use_planetary_gravity != 0 && planetary_body_id.empty())
		{
			msg = "planetary_gravity is enabled, but planetary_body is empty. Specify at least one allowed JPL body ID.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// Lunar non-spherical gravity
		si_EBCRSInput.GetInteger(tracer, "use_lunar_NSG", use_lunar_NSG);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetInteger(tracer, "lunar_zonal_order", lunar_zonal_order);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetInteger(tracer, "lunar_tesseral_degree", lunar_tesseral_degree);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetString(tracer, "file_lunar_gravity_field", fn_lunar_gravity_field);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// Drag, SRP, and post-Newtonian switches
		si_EBCRSInput.GetInteger(tracer, "use_drag_perturb", use_drag);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetInteger(tracer, "use_SRP_perturb", use_SRP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetInteger(tracer, "use_post_newtonian", use_post_newtonian);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// Cartesian eph file
		si_EBCRSInput.GetString(tracer, "file_eph_Cartesian", fn_ephCart);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_report：保存任务设置、动力学模型和传播结果；省略或留空时关闭。
		if (si_EBCRSInput.HasKeyword("file_report"))
		{
			si_EBCRSInput.GetString(tracer, "file_report", fn_report);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// output state normalization
		si_EBCRSInput.GetInteger(tracer, "output_state_units", eph_normalize);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 是否把本次 EBCRS 预报任务的完整配置总览输出到控制台。报告文件始终会写入完整总览。
		si_EBCRSInput.GetInteger(tracer, "display_summary", display_summary);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (display_summary != 0 && display_summary != 1)
		{
			msg = "Invalid display_summary. Use 0 (hide model and integrator details) or 1 (show details).";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// satellite id and physical parameters
		si_EBCRSInput.GetString(tracer, "satellite_id", sat_id);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetDouble(tracer, "satellite_mass", sat_mass);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetDouble(tracer, "satellite_drag_cross_area", sat_area_drag);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetDouble(tracer, "satellite_drag_Cd", sat_Cd_drag);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetDouble(tracer, "satellite_srp_cross_area", sat_area_SRP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_EBCRSInput.GetDouble(tracer, "satellite_srp_eta", sat_eta_SRP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		si_EBCRSInput.GetInteger(tracer, "display_progress", display_progress);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (display_progress < 0 || display_progress > 2)
			throw std::runtime_error("display_progress must be 0 (quiet), 1 (progress) or 2 (all steps).");
		// 核对设置与文件。可选输出留空表示关闭；停用的摄动力不要求相关资料存在。
		if (!std::isfinite(t0_mjd_tdb)) throw std::runtime_error("reference_epoch must be finite.");
		if (!std::isfinite(propagation_length_sec) || !std::isfinite(integrator_step_sec) ||
			(propagation_length_sec != 0 && integrator_step_sec == 0) || integrator_maxTE <= 0)
			throw std::runtime_error("Invalid propagation_length, integrator_stepsize or integrator_maxTE.");
		if (eph_normalize < 1 || eph_normalize > 3) throw std::runtime_error("output_state_units must be 1, 2 or 3.");
		if (solar_gravity_config < 0 || solar_gravity_config > 2 || lunar_gravity_config < 0 || lunar_gravity_config > 2)
			throw std::runtime_error("solar_gravity and lunar_gravity must be 0, 1 or 2.");
		for (int flag : {use_lunar_NSG, use_drag, use_SRP, use_post_newtonian, use_planetary_gravity, display_summary})
			if (flag != 0 && flag != 1) throw std::runtime_error("Force and display switches must be 0 or 1.");
		if (earth_zonal_order < 0 || earth_zonal_order > 80 || earth_tesseral_degree < 0 || earth_tesseral_degree > 80)
			throw std::runtime_error("Earth gravity settings must be within the built-in GGM03C range 0..80.");
		if ((use_drag || use_SRP) && (sat_mass <= 0 || sat_area_drag < 0 || sat_area_SRP < 0))
			throw std::runtime_error("Drag/SRP requires positive satellite_mass and nonnegative areas.");
		// 六个分量依次对应位置和速度，仍采用当前归一化单位检查有限性。
		double initial_components[6];
		initial_cart.ToArray(initial_components);
		for (int i = 0; i < 6; i++)
			if (!std::isfinite(initial_components[i])) throw std::runtime_error("Initial Cartesian state is not finite.");
		// 把非空资料路径相对于设置目录展开；停用模型的资料不强制加载。
		// 这里只整理路径，实际读取检查仍由下面各资料读取器完成。
		for (auto* filename : {&fn_EOP, &fn_JPL, &fn_atm_density, &fn_lunar_gravity_field})
		{
			if (filename->empty()) continue;
			const fs::path path(*filename);
			*filename = (path.is_absolute() ? path : config_path.parent_path() / path).lexically_normal().string();
		}

		/* MARK: ------------------------- set EBCRS model and setting ---------------------------- */

		// EOP 供地球非球形引力和大气阻力内部的地固系计算使用。
		EOP thisEOP;
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading EOP from: " + fn_EOP, EventType::MESSAGE, FORCE_ALERT_ON);
		thisEOP.Read_EOP(tracer, fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// JPL 星历提供当前 TDB 时刻的天体状态，供引力及后牛顿修正使用。
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading JPL ephemeris from: " + fn_JPL, EventType::MESSAGE, FORCE_ALERT_ON);
		status = Read_JPLEPH(fn_JPL);
		if (status <= 0)
		{
			msg = "Failed to read JPL ephemeris: " + fn_JPL;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		else if (status == 2)
		{
			msg = "JPL ephemeris has already been loaded; reuse the existing global JPL handle.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE);
		}

		// 设置地球非球形引力
		if (earth_zonal_order < 0 || earth_tesseral_degree < 0)
		{
			msg = "Invalid earth_zonal_order or earth_tesseral_degree!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		else if (earth_zonal_order == 0 && earth_tesseral_degree == 0)
		{
			EBCRS_model.use_earth_NSG = false;
			EBCRS_model.earth_zonal_order = 0;
			EBCRS_model.earth_tesseral_degree = 0;
		}
		else
		{
			EBCRS_model.use_earth_NSG = true;
			EBCRS_model.earth_zonal_order = earth_zonal_order;
			EBCRS_model.earth_tesseral_degree = earth_tesseral_degree;
		}

		// 设置第三体和行星摄动
		EBCRS_model.solar_gravity_config = solar_gravity_config;
		EBCRS_model.lunar_gravity_config = lunar_gravity_config;
		EBCRS_model.use_planetary_gravity = use_planetary_gravity != 0;
		if (EBCRS_model.use_planetary_gravity)
		{
			EBCRS_model.planetary_body_id = planetary_body_id;
		}

		// 设置月球非球形引力
		if (use_lunar_NSG)
		{
			if (lunar_zonal_order < 0 || lunar_tesseral_degree < 0)
			{
				msg = "Invalid lunar earth_zonal_order or earth_tesseral_degree!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			else if (lunar_zonal_order == 0 && lunar_tesseral_degree == 0)
			{
				msg = "No lunar NSG perturbation!";
				EBCRS_model.use_lunar_NSG = false;
				EBCRS_model.lunar_zonal_order = 0;
				EBCRS_model.lunar_tesseral_degree = 0;
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE);
			}
			else
			{
				if (fn_lunar_gravity_field.empty())
				{
					msg = "file_lunar_gravity_field is missing while use_lunar_NSG = 1.";
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}

				msg = "Reading lunar gravity field from: " + fn_lunar_gravity_field;
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
				if (!moon_gravity_field.Load(fn_lunar_gravity_field))
				{
					msg = "Failed to load lunar gravity field file for lunar NSG perturbation!";
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}

				if (lunar_zonal_order > moon_gravity_field.max_degree || lunar_tesseral_degree > moon_gravity_field.max_degree)
				{
					msg = "Lunar earth_zonal_order or earth_tesseral_degree exceeds the maximum order and degree of the loaded lunar gravity field!";
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}

				EBCRS_model.use_lunar_NSG = true;
				EBCRS_model.lunar_zonal_order = lunar_zonal_order;
				EBCRS_model.lunar_tesseral_degree = lunar_tesseral_degree;
			}
		}

		// 设置大气阻力。地月空间通常不建议打开，但接口沿用现有指数大气模型，近地 EBCRS 测试可用。
		EBCRS_model.use_drag = use_drag != 0;
		if (EBCRS_model.use_drag)
		{
			if (fn_atm_density.empty())
			{
				msg = "file_atmosphere_density is missing while use_drag_perturb = 1.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			double atm_step = 0.0;
			msg = "Reading atmospheric density from: " + fn_atm_density;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
			// 密度数组供动力学使用，表头条件供动力学 Display 显示。
			status = Read_AtmosphereDensity(fn_atm_density, atm_list, atm_step, &EBCRS_model.atmosphere_Ap, &EBCRS_model.atmosphere_F107);
			if (status <= 0 || atm_list.empty() || !std::isfinite(atm_step) || atm_step <= 0)
				throw std::runtime_error("Failed to load a usable atmosphere density table: " + fn_atm_density);
			if (status != 1)
			{
				msg = "Failed to read atmosphere density file for EBCRS drag perturbation!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 设置光压
		EBCRS_model.use_SRP = use_SRP != 0;

		// 设置后牛顿开关
		EBCRS_model.use_post_newtonian = use_post_newtonian != 0;
		if (EBCRS_model.use_post_newtonian)
		{
			msg = "Post-Newtonian EIH correction is enabled. The EBCRS acceleration will use all 11 JPL major bodies for F_PN = F_PN^S - F_PN^E.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
		}

		// 读入卫星模型
		EBCRS_model.sat.Id = sat_id;
		if (sat_mass > 0.0)
		{
			EBCRS_model.sat.mass = sat_mass / MASS_OF_EARTH;
		}
		if (sat_area_drag > 0.0)
		{
			EBCRS_model.sat.area_drag = sat_area_drag / ae / ae;
		}
		if (EBCRS_model.sat.mass > 0.0 && EBCRS_model.sat.area_drag > 0.0)
		{
			EBCRS_model.sat.AMR_drag = EBCRS_model.sat.area_drag / EBCRS_model.sat.mass;
		}
		if (sat_area_SRP > 0.0)
		{
			EBCRS_model.sat.area_SRP = sat_area_SRP / ae / ae;
		}
		if (EBCRS_model.sat.mass > 0.0 && EBCRS_model.sat.area_SRP > 0.0)
		{
			EBCRS_model.sat.AMR_SRP = EBCRS_model.sat.area_SRP / EBCRS_model.sat.mass;
		}
		EBCRS_model.sat.eta_SRP = sat_eta_SRP;
		EBCRS_model.sat.Cd = sat_Cd_drag;

		// 共同预报参数与积分器专用参数分开；时间仍为 TDB，状态仍为归一化 Cartesian。
		EBCRS_setting.t0_MJD_TDB = t0_mjd_tdb;
		EBCRS_setting.display_progress = display_progress;
		EBCRS_setting.duration_norm = propagation_length_sec * SEC_TO_NORMALIZED_TIME;
		EBCRS_setting.h0_norm = integrator_step_sec * SEC_TO_NORMALIZED_TIME;
		EBCRS_setting.rkf.maxTE_norm = integrator_maxTE;
		EBCRS_setting.rkf.fixed_step = integrator_id == 102;
		EBCRS_setting.adams.steps = adams_steps;
		EBCRS_setting.adams.use_correction = integrator_id == 202;

		// 最后配置当前任务
		EBCRS_job.opn_model = EBCRS_model;
		EBCRS_job.opn_setting = EBCRS_setting;
		EBCRS_job.initCart_norm = initial_cart;
		EBCRS_job.initKep_norm = initial_kep;
		EBCRS_job.opn_EOP = thisEOP;
		EBCRS_job.atm_list = atm_list;
		if (EBCRS_model.use_drag)
		{
			Set_ReferenceDensity(EBCRS_job.atm_list, EBCRS_job.initKep_norm, EBCRS_job.refDensity);
		}
		// 四、展开各项输出路径；空文件名关闭对应输出，实际写入时直接打开文件。
		if (fn_ephCart.empty())
			throw std::runtime_error("Enable at least one ephemeris output.");
		// 相对文件名从设置文件所在目录展开，各星历与报告分别保留其指定位置。
		if (!fn_ephCart.empty() && fs::path(fn_ephCart).is_relative())
			fn_ephCart = (config_path.parent_path() / fn_ephCart).lexically_normal().string();
		if (!fn_report.empty() && fs::path(fn_report).is_relative())
			fn_report = (config_path.parent_path() / fn_report).lexically_normal().string();


		// =========================== 05 显示任务、动力学与积分设置 ===========================

		tracer.LogEvent(__FILE__, __func__, __LINE__, "Input checks and model preparation completed.", EventType::MESSAGE, FORCE_ALERT_ON);

		// 显示只恢复当前 E-BCRS 状态副本的单位，原归一化初态留给积分使用。
		Cartesian initial_cart_km = initial_cart;
		initial_cart_km.Restore(tracer, ae / 1000.0, UNITTIME);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		double initial_cart_values[6];
		initial_cart_km.ToArray(initial_cart_values);

		// TDB 历元保留数值日，时长保留正负号；地月距离使用足够的列宽对齐位置和速度。
		const int label_width = 31, state_width = 16;
		std::ostringstream setup;
		setup << "\n=================================== General settings ===================================\n\n"
			  << "- Propagation Frame / Time     : E-BCRS / TDB"
			  << "\n- Input / Output Time          : TDB"
			  << "\n- Target ID                    : " << sat_id
			  << "\n- Mass [kg]                    : " << std::defaultfloat << std::setprecision(14) << sat_mass
			  << "\n\n- Start epoch [TDB]            : " << std::fixed << std::setprecision(12) << t0_mjd_tdb
			  << "\n- End epoch [TDB]              : " << t0_mjd_tdb + propagation_length_sec / 86400.0
			  << "\n- Duration [s]                 : " << std::showpos << std::setprecision(2) << propagation_length_sec << std::noshowpos
			  << "\n\n";

		// 两行依次为 [x,y,z]、[vx,vy,vz]，共用列宽以对齐括号和逗号。
		for (int row = 0; row < 2; ++row)
		{
			if (row == 0)
			{
				setup << std::left << std::setw(label_width) << "- Initial Cartesian [km, km/s]" << ": ";
			}
			else
			{
				setup << std::string(label_width + 2, ' ');
			}

			setup << "[ ";

			for (int j = 0; j < 3; ++j)
			{
				if (j > 0) setup << ", ";

				setup << std::fixed << std::setprecision(6) << std::right << std::setw(state_width) << initial_cart_values[row * 3 + j];
			}

			setup << " ]\n";
		}

		// 控制台详情由现有 display_summary 开关控制；报告始终保存完整模型和积分设置。
		const string general_settings = setup.str();

		if (display_summary)
		{
			EBCRS_model.Display(tracer, setup);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			EBCRS_setting.Display(tracer, integrator_id, setup);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		std::cerr << setup.str() << '\n';

		// =========================== 06 执行传播并保存最终状态副本 ===========================

		vector<double> record_MJD_TDB;
		vector<Cartesian> record_cart;

		const string integrator_name = integrator_id == 101 ? "adaptive RKF78" : integrator_id == 102 ? "fixed-step RKF78"
			: integrator_id == 201 ? "Adams PE" : "Adams PECE";
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Starting E-BCRS propagation with " + integrator_name + ".", EventType::MESSAGE, FORCE_ALERT_ON);

		if (display_progress > 0) std::cerr << "\n================================= Propagation progress =================================\n\n";

		// 只对实际传播调用计时；输入、积分和保存节点均使用 E-BCRS/TDB。
		const auto propagation_start = std::chrono::steady_clock::now();

		if (integrator_id == 201 || integrator_id == 202)
		{
			OrbitPropagation_Adams_Cart_EBCRS(tracer, EBCRS_job, record_MJD_TDB, record_cart);
		}
		else
		{
			OrbitPropagation_RKF78_Cart_EBCRS(tracer, EBCRS_job, record_MJD_TDB, record_cart);
		}

		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		const double propagation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - propagation_start).count();

		if (display_progress > 0 && record_cart.size() > 1) std::cerr << '\n';

		tracer.LogEvent(__FILE__, __func__, __LINE__, "Requested final epoch reached.", EventType::MESSAGE, FORCE_ALERT_ON);

		// 摘要复用已保存的末状态，只在副本上恢复 km、km/s，不再传播或变换坐标系。
		Cartesian final_cart_km = record_cart.back();
		final_cart_km.Restore(tracer, ae / 1000.0, UNITTIME);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		double final_cart_values[6];
		final_cart_km.ToArray(final_cart_values);

		// =========================== 07 按配置恢复单位并写出星历 ===========================

		// 写出事件在对应位置登记；路径消息集中显示于完成消息之前。
		std::ostringstream file_messages;

		if (!fn_ephCart.empty())
		{
			ofstream ofs_ephCart(fn_ephCart.c_str());
			if (!ofs_ephCart)
			{
				msg = "Cannot open file_eph_Cartesian: " + fn_ephCart + " to write!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			else
			{
				msg = "Writing Cartesian ephemeris to: " + fn_ephCart;
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
				file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';

				vector<Cartesian> output_cart;
				for (int k = 0; k < record_cart.size(); k++)
				{
					Cartesian temp_cart = record_cart[k];

					switch (eph_normalize)
					{
						case STATE_UNIT_NORMALIZE:
							break;
						case STATE_UNIT_KM_SEC_DEGREE:
							temp_cart.Restore(tracer, ae / 1000.0, UNITTIME);
							if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
							break;
						case STATE_UNIT_M_SEC_DEGREE:
							temp_cart.Restore(tracer, ae, UNITTIME);
							if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
							break;
						default:
							msg = "Unclear ephemeris unit! Use normalized unit in default!";
							tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
							break;
					}

					output_cart.push_back(temp_cart);
				}

				double cartArr[6];
				for (int k = 0; k < output_cart.size(); k++)
				{
					output_cart[k].ToArray(cartArr);
					writeDouble(ofs_ephCart, record_MJD_TDB[k]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[0]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[1]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[2]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[3]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[4]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[5]); ofs_ephCart << ' ';
					ofs_ephCart << endl;
				}
				ofs_ephCart.close();
				if (!ofs_ephCart) throw std::runtime_error("Failed while writing EBCRS ephemeris.");
			}
		}
		else
		{
			msg = "No EBCRS Cartesian ephemeris file is specified.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
		}

		// =========================== 08 完成摘要、报告与程序日志 ===========================

		// 记录整个进程截至此处的峰值常驻内存；macOS 返回 byte，Linux 返回 KiB，统一换算为 MiB。
		double peak_memory_mib = std::numeric_limits<double>::quiet_NaN();

#if defined(__APPLE__) || defined(__linux__)
		struct rusage usage{};

		if (getrusage(RUSAGE_SELF, &usage) == 0)
		{
#if defined(__APPLE__)
			peak_memory_mib = static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
			peak_memory_mib = static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
		}
		else
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot obtain peak resident memory usage.", EventType::WARNING);
		}
#endif

		// 摘要固定显示物理状态，独立于星历文件的表示和单位；时长由实际保存的首末历元计算。
		std::ostringstream summary;
		summary << "\n================================= Propagation summary ==================================\n\n"
				<< "- Status                       : COMPLETED"
				<< "\n- Final epoch [TDB]            : " << std::fixed << std::setprecision(12) << record_MJD_TDB.back()
				<< "\n- Propagated duration [s]      : " << std::showpos << std::setprecision(2)
				<< (record_MJD_TDB.back() - EBCRS_setting.t0_MJD_TDB) * 86400.0 << std::noshowpos
				<< "\n- Ephemeris records            : " << record_cart.size() << "\n\n";

		for (int row = 0; row < 2; ++row)
		{
			if (row == 0)
			{
				summary << std::left << std::setw(label_width) << "- Final Cartesian [km, km/s]" << ": ";
			}
			else
			{
				summary << std::string(label_width + 2, ' ');
			}

			summary << "[ ";

			for (int j = 0; j < 3; ++j)
			{
				if (j > 0) summary << ", ";

				summary << std::fixed << std::setprecision(6) << std::right << std::setw(state_width) << final_cart_values[row * 3 + j];
			}

			summary << " ]\n";
		}

		summary << "\n- Propagation time [s]         : " << std::fixed << std::setprecision(2) << propagation_seconds
				<< "\n- Total elapsed time [s]       : " << std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count()
				<< "\n- Peak memory [MiB]            : ";

		if (std::isfinite(peak_memory_mib))
		{
			summary << std::fixed << std::setprecision(1) << peak_memory_mib;
		}
		else
		{
			summary << "Not available";
		}

		summary << '\n';

		// 任务报告直接写设置和结果；运行事件与源码位置集中保存在独立的程序日志中。
		if (!fn_report.empty())
		{
			msg = "Writing propagation report to: " + fn_report;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
			std::ofstream report(fn_report);

			if (!report)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open propagation report: " + fn_report, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			report << "E-BCRS orbit propagation report\n"
				   << "Program       : " << program_name << " | Version: v" << ASTER_PROGRAM_VERSION
				   << "\nJob ID        : " << jobId << "\nStarted at    : " << timeStr
				   << "\nInput setting : " << config_path.string() << '\n'
				   << general_settings;
			EBCRS_model.Display(tracer, report);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			EBCRS_setting.Display(tracer, integrator_id, report);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			report << summary.str();
			report.close();

			if (!report)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing report: " + fn_report, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 程序日志的路径也作为独立事件记录，并与其他写出路径一起在结束消息前显示。
		if (program_log.is_open())
		{
			msg = "Writing program log to: " + fn_programlog;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
		}

		// 成功事件先记入日志，只有日志也保存成功后才在控制台显示摘要和结束消息。
		const string completed = "E-BCRS propagation job completed successfully.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, completed, EventType::MESSAGE, FORCE_ALERT_OFF);

		if (program_log.is_open())
		{
			tracer.ReportEvents(program_log);
			program_log << "\nExit status: 0\n";
			program_log.close();

			if (!program_log)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing program log.", EventType::ERROR);
				tracer.EmergencyExit();
			}
		}

		std::cerr << summary.str() << '\n' << file_messages.str();
		std::cerr << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << completed << '\n' << std::setfill(' ');

		tracer.CheckOut(__func__);
		return 0;
	}
	catch (const std::exception& error)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
		tracer.EmergencyExit(program_log);
		return 1;
	}
}
