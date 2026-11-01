/**
 * @file main_opn_nonsingular2.cpp
 * @brief 在 GCRS/TT 中传播 NSG2，输出 UTC 历元的轨道星历、任务摘要和运行事件。
 * @details 主程序依次读取设置、准备 EOP 与力模型、显示任务参数、调用积分、转换并写出星历。
 * 输入初态为 GCRS，历元为 UTC MJD，单位 km、km/s、deg；计算状态使用本库归一化尺度。
 * 控制台初末轨道固定显示 Keplerian 和 Cartesian，两行分别保存前三项与后三项。
 * 星历单位由配置控制；相对路径以 setting 所在目录为基准，文件读取与写出路径由 tracer 记录。
 * file_report 保存任务设置与结果，file_programlog 保存原始配置、运行事件及失败调用链，均可单独关闭。
 * 用法：aster_opn_nonsingular2 -I <配置>。动力学与积分设置由已有模型和设置类的 Display 显示。
 */
#include "io_settings.hpp"
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
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

// 输出星历的单位
#define STATE_UNIT_NORMALIZE 1
#define STATE_UNIT_KM_SEC_DEGREE 2
#define STATE_UNIT_M_SEC_DEGREE 3

int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_opn_nonsingular2 " << ASTER_PROGRAM_VERSION << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_opn_nonsingular2 -I <setting.txt>\n";
        return 2;
    }

	// 登记 数值法轨道预报主程序开始
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
		jobMsg << "NSG2 orbit propagation job started at " << timeStr << " (jobId = " << jobId
			   << ")";
		tracer.LogEvent(__FILE__, __func__, __LINE__, jobMsg.str(), bastion::EventType::MESSAGE, FORCE_ALERT_ON);

		// 一些程序内通用定义
		string msg;    // 用于记录事件

		namespace fs = std::filesystem;

		// 程序名取实际执行文件名，版本字符串由 CMake 统一提供。
		const string program_name = fs::path(argv[0]).filename().string();
		tracer.LogEvent(__FILE__, __func__, __LINE__,
			"Program: " + program_name + " | Version: v" + ASTER_PROGRAM_VERSION, EventType::MESSAGE, FORCE_ALERT_ON);

		// 默认读取随入口提供的设置，命令行可指定其他配置文件。
		std::string fn_opn_input = "";

		// 按顺序处理参数；缺值和未知选项均报错，避免误用默认设置。
		for (int i = 1; i < argc; ++i)
		{
			const std::string arg = argv[i];
			if (arg == "-h" || arg == "--help")
			{
				std::cout << "Usage: aster_opn_nonsingular2 [-I|-input|--input <setting>]\n"
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
			fn_opn_input = argv[++i];
		}

		/* MARK: --------------------------------- read input ------------------------------------ */

		msg = "Reading and analyzing input from: " + fs::absolute(fn_opn_input).lexically_normal().string();
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		// SimpleInput 读入 key = value，再按下面的科学计算步骤逐项取值。
		SimpleInput si_opnInput;
		OrbitPropagationNModel opn_model;    //	动力学模型
		OrbitPropagationSetting opn_setting; // 共同预报设置及各积分器参数。
		OrbitPropagationJob opn_job;

		// 二、读取初态、摄动力及输出设置。初态转换后统一使用归一化单位。
		si_opnInput.ReadSimpleInput(tracer, fn_opn_input);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_programlog：保存原始配置、运行事件与失败调用链；省略或留空时关闭。
		// 相对路径以 setting 所在目录为基准，与 file_report 的任务报告分别输出。
		std::string fn_programlog;

		if (si_opnInput.HasKeyword("file_programlog"))
		{
			si_opnInput.GetString(tracer, "file_programlog", fn_programlog);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		if (!fn_programlog.empty())
		{
		    const std::filesystem::path log_config = std::filesystem::weakly_canonical(std::filesystem::absolute(fn_opn_input));
		    std::filesystem::path log_path = fn_programlog;
		    if (log_path.is_relative()) log_path = log_config.parent_path() / log_path;
		    log_path = std::filesystem::weakly_canonical(log_path);
		    fn_programlog = log_path.string();
		    // 日志不能覆盖已有资料，也不能与本次计划写出的其他文件重名。
		    if (std::filesystem::exists(log_path))
		        throw std::runtime_error("Program log already exists: " + log_path.string());
		    for (std::size_t record_index = 1; record_index <= si_opnInput.Count(); ++record_index)
		    {
		        io::SimpleRecord path_record;
		        si_opnInput.GetRecord(tracer, record_index, path_record);
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
		    for (std::size_t record_index = 1; record_index <= si_opnInput.Count(); ++record_index)
		    {
		        io::SimpleRecord path_record;
		        si_opnInput.GetRecord(tracer, record_index, path_record);
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
		const fs::path config_path = fs::weakly_canonical(fs::absolute(fn_opn_input));

		// 初态只选一种表示；积分器内部统一采用归一化 NSG2。
		bool isKepReady = false, isCartReady = false;
		double t0_mjd_utc = 0, opn_length_sec = 0, integrator_step_sec = 0, integrator_maxTE = 0;
		string str_epoch, str_kep, str_cart;
		string fn_EOP, fn_JPL, fn_ephCart, fn_ephKep, fn_ephNSG2, fn_report, fn_atm_density;
		string sat_id;
		double sat_mass, sat_area_drag, sat_area_SRP, sat_eta_SRP, sat_Cd_drag;
		int adams_steps = 9;
		int integrator_id, earth_zonal_order, earth_tesseral_degree, solar_gravity_config, lunar_gravity_config, use_drag = 0, use_SRP = 0;
		// 地球与月球非球形引力分别配置，月球项还需要外部重力场。
		int use_lunar_NSG = 0, lunar_zonal_order = 0, lunar_tesseral_degree = 0;
		string fn_lunar_gravity_field;

		int eph_normalize, display_model, display_setting;
		int display_progress = 1; // 0 静默；1 每 10% 输出进度；2 输出每步及试算误差。
		int status;
		vector<string> vec_epoch, vec_kep, vec_cart;
		Keplerian initial_kep;
		Cartesian initial_cart;

		// initial epoch
		si_opnInput.GetString(tracer, "reference_epoch", str_epoch);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!str_epoch.empty())    // 在ReadSimpleInput就已经strip过，因此只考虑是不是空的""就行
		{
			split(str_epoch, vec_epoch);
			switch (vec_epoch.size())
			{
				case 1:
					// 原本的str_epoch仅为一个数字，spilt完被原封不动加进vec_epoch，作为第一个也是唯一一个元素
					// 对应的日期格式为JD或MJD，统一转换存储为MJD
					t0_mjd_utc = stod(vec_epoch[0]) > 2400000.5 ? stod(vec_epoch[0]) - 2400000.5 : stod(vec_epoch[0]);
					break;
				case 6:
					// 原本的str_epoch仅为年月日时分秒 2025 05 25 12 40 0.322
					Calendar2MJD(stoi(vec_epoch[0]), stoi(vec_epoch[1]), stoi(vec_epoch[2]), stoi(vec_epoch[3]), stoi(vec_epoch[4]), stod(vec_epoch[5]), t0_mjd_utc);
					break;
				default:
					msg = "Unknown format for parameter reference_epoch!";
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
					tracer.EmergencyExit(program_log);
					break;
			}
		}
		else
		{
			msg = "Parameter reference_epoch is missing!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// try initial Keplerian
		si_opnInput.GetString(tracer, "reference_Keplerian", str_kep);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!str_kep.empty())
		{
			split(str_kep, vec_kep);
			if (vec_kep.size() == 6)
			{
				// 输入只接受 km deg！
				initial_kep.Set(stod(vec_kep[0]), stod(vec_kep[1]), stod(vec_kep[2]), stod(vec_kep[3]), stod(vec_kep[4]), stod(vec_kep[5]), false);
				initial_kep.Normalize(tracer, ae / 1000.0);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
				// 转换一下
				initial_kep.ToCartesian(tracer, initial_cart);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

				isKepReady = true;
			}
			else
			{
				msg = "Invalid format of parameter reference_Keplerian!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// try initial Cartesian
		si_opnInput.GetString(tracer, "reference_Cartesian", str_cart);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (!str_cart.empty())
		{
			split(str_cart, vec_cart);
			if (vec_cart.size() == 6)
			{
				// 输入只接受 km km/s！
				initial_cart.Set(stod(vec_cart[0]), stod(vec_cart[1]), stod(vec_cart[2]), stod(vec_cart[3]), stod(vec_cart[4]), stod(vec_cart[5]), false);
				initial_cart.Normalize(tracer, ae / 1000.0, UNITTIME);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
				// 转换一下
				initial_cart.ToKeplerian(tracer, initial_kep);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

				isCartReady = true;
			}
			else
			{
				msg = "Invalid format of parameter reference_Cartesian!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 初态必须且只能指定一种；两种表示同时给出或都留空均报错。
		if (isCartReady && isKepReady)
		{
			msg = "Both Cartesian and Keplerian are fed as initial states. Select only one type!";
			// 同时指定了直角坐标和Kepler根数，不对！...........
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		else if (!isCartReady && !isKepReady)
		{
			// 直角坐标和Kepler根数都没有，也不对！...........
			msg = "Neither Cartesian or Keplerian is fed as initial states!";
			// 同时指定了直角坐标和Kepler根数，不对！...........
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// propagation length
		si_opnInput.GetDouble(tracer, "propagation_length", opn_length_sec);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// EOP file
		si_opnInput.GetString(tracer, "file_EOP", fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// JPLEPH file
		si_opnInput.GetString(tracer, "file_JPLEPH", fn_JPL);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// atmosphere density file
		si_opnInput.GetString(tracer, "file_atmosphere_density", fn_atm_density);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 101 为自适应 RKF78，102 为定步长 RKF78；201／202 为 Adams PE／PECE。
		si_opnInput.GetInteger(tracer, "integrator_id", integrator_id);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (integrator_id != 101 && integrator_id != 102 && integrator_id != 201 && integrator_id != 202)
			throw std::runtime_error("Supported NSG2 integrator_id values: 101=adaptive RKF78, 102=fixed-step RKF78, 201=Adams PE, 202=Adams PECE.");
		if (integrator_id == 201 || integrator_id == 202)
		{
			si_opnInput.GetInteger(tracer, "integrator_adams_steps", adams_steps);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
			if (adams_steps < 2 || adams_steps > 9)
				throw std::runtime_error("integrator_adams_steps must be an integer in 2..9.");
		}

		// integrator stepsize
		si_opnInput.GetDouble(tracer, "integrator_stepsize", integrator_step_sec);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// integrator maxTE
		si_opnInput.GetDouble(tracer, "integrator_maxTE", integrator_maxTE);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		//Earth gravity order and degree
		si_opnInput.GetInteger(tracer, "earth_zonal_order", earth_zonal_order);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_opnInput.GetInteger(tracer, "earth_tesseral_degree", earth_tesseral_degree);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// solar gravity perturbation (central body)
		si_opnInput.GetInteger(tracer, "solar_gravity", solar_gravity_config);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// lunar gravity perturbation (central body)
		si_opnInput.GetInteger(tracer, "lunar_gravity", lunar_gravity_config);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 260211 lunar NSG perturbation
		si_opnInput.GetInteger(tracer, "use_lunar_NSG", use_lunar_NSG);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_opnInput.GetInteger(tracer, "lunar_zonal_order", lunar_zonal_order);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		si_opnInput.GetInteger(tracer, "lunar_tesseral_degree", lunar_tesseral_degree);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		// lunar gravity field file for lunar NSG
		si_opnInput.GetString(tracer, "file_lunar_gravity_field", fn_lunar_gravity_field);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// drag perturbation
		si_opnInput.GetInteger(tracer, "use_drag_perturb", use_drag);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// SRP perturbation
		si_opnInput.GetInteger(tracer, "use_SRP_perturb", use_SRP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// Cartesian eph file
		si_opnInput.GetString(tracer, "file_eph_Cartesian", fn_ephCart);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// Keplerian eph file
		si_opnInput.GetString(tracer, "file_eph_Keplerian", fn_ephKep);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// NonSingular2 eph file
		si_opnInput.GetString(tracer, "file_eph_NonSingular2", fn_ephNSG2);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_report：保存任务设置、动力学模型、积分参数和最终结果；省略或留空时关闭。
		if (si_opnInput.HasKeyword("file_report"))
		{
			si_opnInput.GetString(tracer, "file_report", fn_report);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// output state normalization 读取输出的单位要求 跟导师的有区别
		si_opnInput.GetInteger(tracer, "output_state_units", eph_normalize);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// display model
		si_opnInput.GetInteger(tracer, "display_model", display_model);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// display setting
		si_opnInput.GetInteger(tracer, "display_setting", display_setting);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 进度显示：0 静默，1 简洁进度，2 全部接受步与截断误差。
		si_opnInput.GetInteger(tracer, "display_progress", display_progress);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		if (display_progress < 0 || display_progress > 2)
			throw std::runtime_error("display_progress must be 0 (quiet), 1 (progress) or 2 (all steps).");

		// satellite id
		si_opnInput.GetString(tracer, "satellite_id", sat_id);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// satellite mass kg
		si_opnInput.GetDouble(tracer, "satellite_mass", sat_mass);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// satellite cross-section area of drag m^2
		si_opnInput.GetDouble(tracer, "satellite_drag_cross_area", sat_area_drag);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// satellite Cd coefficient of drag
		si_opnInput.GetDouble(tracer, "satellite_drag_Cd", sat_Cd_drag);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// satellite cross-section area of SRP m^2
		si_opnInput.GetDouble(tracer, "satellite_srp_cross_area", sat_area_SRP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// satellite eta of SRP
		si_opnInput.GetDouble(tracer, "satellite_srp_eta", sat_eta_SRP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 核对设置与文件。可选输出留空表示关闭；停用的摄动力不要求相关资料存在。
		if (!std::isfinite(t0_mjd_utc)) throw std::runtime_error("reference_epoch must be finite.");
		if (!std::isfinite(opn_length_sec) || !std::isfinite(integrator_step_sec) ||
			(opn_length_sec != 0 && integrator_step_sec == 0) || integrator_maxTE <= 0)
			throw std::runtime_error("Invalid propagation_length, integrator_stepsize or integrator_maxTE.");
		if (eph_normalize < 1 || eph_normalize > 3) throw std::runtime_error("output_state_units must be 1, 2 or 3.");
		if (solar_gravity_config < 0 || solar_gravity_config > 2 || lunar_gravity_config < 0 || lunar_gravity_config > 2)
			throw std::runtime_error("solar_gravity and lunar_gravity must be 0, 1 or 2.");
		for (int flag : {use_lunar_NSG, use_drag, use_SRP, display_model, display_setting})
			if (flag != 0 && flag != 1) throw std::runtime_error("Force and model/setting display switches must be 0 or 1.");
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

		// 三、加载实际启用的摄动力资料，配置目标物理参数与积分设置。
		/* MARK: ------------------------- set OPN model and setting ---------------------------- */

		// 先加载 EOP，再设置各项摄动力。

		// 加载EOP数据，Read_EOP函数内部保证同一个实例不会重复读取
		EOP thisEOP;
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading EOP from: " + fn_EOP, EventType::MESSAGE, FORCE_ALERT_ON);
		thisEOP.Read_EOP(tracer, fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 分别设置地球带谐项和田谐项的截断；2 / 0 只开启 J2 摄动。
		if (earth_zonal_order < 0 || earth_tesseral_degree < 0)
		{
			// 非法阶次
			msg = "Invalid earth_zonal_order or earth_tesseral_degree!";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}
		else if (earth_zonal_order == 0 && earth_tesseral_degree == 0)
		{
			// 都为0 质点引力
			msg = "No NSG perturbation!";
			opn_model.use_earth_NSG = false;
			opn_model.earth_zonal_order = 0;
			opn_model.earth_tesseral_degree = 0;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE);
		}
		else
		{
			// 开启非球形引力
			opn_model.use_earth_NSG = true;
			opn_model.earth_zonal_order = earth_zonal_order;
			opn_model.earth_tesseral_degree = earth_tesseral_degree;
		}

		// 设置第三体摄动
		opn_model.solar_gravity_config = solar_gravity_config;
		opn_model.lunar_gravity_config = lunar_gravity_config;
		if (solar_gravity_config == SOLAR_GRAVITY_JPL || lunar_gravity_config == LUNAR_GRAVITY_JPL)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading JPL ephemeris from: " + fn_JPL, EventType::MESSAGE, FORCE_ALERT_ON);
			status = Read_JPLEPH(fn_JPL);    // 有需要再读JPL
			if (status <= 0)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed to read JPL ephemeris: " + fn_JPL, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 设置月球非球形引力：月固坐标转换要求使用 JPL 月球星历。
		if (use_lunar_NSG)
		{
			// 使用月球非球形需要lunar_gravity_config == LUNAR_GRAVITY_JPL
			if (lunar_gravity_config != LUNAR_GRAVITY_JPL)
			{
				msg = "Lunar gravity model must be JPL when using lunar NSG perturbation!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			else
			{
				// 先排除有效性问题，再开启月球非球形引力
				if (lunar_zonal_order < 0 || lunar_tesseral_degree < 0)
				{
					// 非法阶次
					msg = "Invalid lunar earth_zonal_order or earth_tesseral_degree!";
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}
				else if (lunar_zonal_order == 0 && lunar_tesseral_degree == 0)
				{
					// 都为0 质点引力
					msg = "No lunar NSG perturbation!";
					opn_model.use_lunar_NSG = false;
					opn_model.lunar_zonal_order = 0;
					opn_model.lunar_tesseral_degree = 0;
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE);
				}
				else
				{
					// 启用月球非球形引力，并加载归一化谐系数。
					opn_model.use_lunar_NSG = true;
					// 由于月球非球形引力的计算需要用到月球的重力场数据，因此在此加载月球重力场文件
					msg = "Reading lunar gravity field from: " + fn_lunar_gravity_field;
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
					if (!moon_gravity_field.Load(fn_lunar_gravity_field))
					{
						msg = "Failed to load lunar gravity field file for lunar NSG perturbation!";
						tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
						tracer.EmergencyExit(program_log);
					}
					// 检查月球非球形引力的阶次是否超过了加载的月球重力场的最大阶次
					if (lunar_zonal_order > moon_gravity_field.max_degree || lunar_tesseral_degree > moon_gravity_field.max_degree)
					{
						msg = "Lunar earth_zonal_order or earth_tesseral_degree exceeds the maximum order and degree of the loaded lunar gravity field!";
						tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
						tracer.EmergencyExit(program_log);
					}
					//没有问题了，设置月球非球形引力的阶次和度数
					msg = "Lunar gravity field loaded successfully. Max degree: " + to_string(moon_gravity_field.max_degree);
					tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

					opn_model.lunar_zonal_order = lunar_zonal_order;
					opn_model.lunar_tesseral_degree = lunar_tesseral_degree;

				}
			}

		}

		// 设置大气阻力
		opn_model.use_drag = use_drag;
		vector<AtmosphereData> atm_list;
		double atm_step = 0;
		if (opn_model.use_drag)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading atmospheric density from: " + fn_atm_density, EventType::MESSAGE, FORCE_ALERT_ON);
			status = Read_AtmosphereDensity(fn_atm_density, atm_list, atm_step, &opn_model.atmosphere_Ap, &opn_model.atmosphere_F107);    // 用了再读
			if (status <= 0 || atm_list.empty() || !std::isfinite(atm_step) || atm_step <= 0)
				throw std::runtime_error("Failed to load a usable atmosphere density table: " + fn_atm_density);
		}

		// 设置光压
		opn_model.use_SRP = use_SRP;

		// 配置目标质量、面积与系数；质量除以地球质量，面积除以 ae²。

		// 读入卫星编号
		opn_model.sat.Id = sat_id;

		// 读入卫星质量，同时归一化...........
		if (sat_mass > 0.0)
		{
			opn_model.sat.mass = sat_mass / MASS_OF_EARTH;
		}
		// 读入迎风向面积，同时归一化...........
		if (sat_area_drag > 0.0)
		{
			opn_model.sat.area_drag = sat_area_drag / ae / ae;
		}
		// 质量和迎风向面积都就位时，计算归一化大气阻力面质比...........
		if (opn_model.sat.mass > 0.0 && opn_model.sat.area_drag > 0.0)
		{
			opn_model.sat.AMR_drag = opn_model.sat.area_drag / opn_model.sat.mass;
		}
		// 读入光压面积，同时归一化...........
		if (sat_area_SRP > 0.0)
		{
			opn_model.sat.area_SRP = sat_area_SRP / ae / ae;
		}
		// 质量和光压面积都就位时，计算归一化光压面质比...........
		if (opn_model.sat.mass > 0.0 && opn_model.sat.area_SRP > 0.0)
		{
			opn_model.sat.AMR_SRP = opn_model.sat.area_SRP / opn_model.sat.mass;
		}
		// 设置光压eta
		opn_model.sat.eta_SRP = sat_eta_SRP;
		// 设置Cd
		opn_model.sat.Cd = sat_Cd_drag;

		// 设置共同的积分起点、归一化时长与步长，再分别填写 RKF 和 Adams 参数。

		vector<double> record_MJD_UTC;
		vector<NonSingular2> record_nsg2;

		opn_setting.t0_MJD_UTC = t0_mjd_utc;
		opn_setting.t0_MJD_TT = UTC2TT(t0_mjd_utc);
		opn_setting.rkf.maxTE_norm = integrator_maxTE;
		opn_setting.duration_norm = opn_length_sec * SEC_TO_NORMALIZED_TIME;
		opn_setting.h0_norm = integrator_step_sec * SEC_TO_NORMALIZED_TIME;
		opn_setting.display_progress = display_progress;
		opn_setting.rkf.fixed_step = integrator_id == 102;
		opn_setting.adams.steps = adams_steps;
		opn_setting.adams.use_correction = integrator_id == 202;

		// 汇总为预报任务，保留 main 中的模型和数据准备顺序。

		// 最后配置当前任务
		opn_job.opn_model = opn_model;
		opn_job.opn_setting = opn_setting;
		opn_job.initCart_norm = initial_cart;
		opn_job.initKep_norm = initial_kep;
		opn_job.opn_EOP = thisEOP;
		opn_job.atm_list = atm_list;
		if (opn_model.use_drag)
		{
			Set_ReferenceDensity(opn_job.atm_list, opn_job.initKep_norm, opn_job.refDensity);
		}

		// 四、展开各项输出路径；空文件名关闭对应输出，实际写入时直接打开文件。
		if (fn_ephCart.empty() && fn_ephKep.empty() && fn_ephNSG2.empty())
			throw std::runtime_error("Enable at least one ephemeris output.");
		// 相对文件名从设置文件所在目录展开，各星历与报告分别保留其指定位置。
		if (!fn_ephCart.empty() && fs::path(fn_ephCart).is_relative())
			fn_ephCart = (config_path.parent_path() / fn_ephCart).lexically_normal().string();
		if (!fn_ephKep.empty() && fs::path(fn_ephKep).is_relative())
			fn_ephKep = (config_path.parent_path() / fn_ephKep).lexically_normal().string();
		if (!fn_ephNSG2.empty() && fs::path(fn_ephNSG2).is_relative())
			fn_ephNSG2 = (config_path.parent_path() / fn_ephNSG2).lexically_normal().string();
		if (!fn_report.empty() && fs::path(fn_report).is_relative())
			fn_report = (config_path.parent_path() / fn_report).lexically_normal().string();

		// =========================== 05 显示任务、动力学与积分设置 ===========================

		tracer.LogEvent(__FILE__, __func__, __LINE__, "Input checks and model preparation completed.", EventType::MESSAGE, FORCE_ALERT_ON);

		// 初态显示统一为 km、km/s、deg；只恢复副本的单位，积分仍使用原归一化状态。
		Cartesian initial_cart_km = initial_cart;
		Keplerian initial_kep_km = initial_kep;
		initial_cart_km.Restore(tracer, ae / 1000.0, UNITTIME);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		initial_kep_km.Restore(tracer, ae / 1000.0);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		double initial_cart_values[6], initial_kep_values[6];
		initial_cart_km.ToArray(initial_cart_values);
		initial_kep_km.ToArray(initial_kep_values);

		// 起终点以 UTC 的 MJD 数值显示，时长保留符号；每种状态分两行、每行三个分量。
		const int label_width = 31;
		const int state_width = 12;
		std::ostringstream setup;
		setup << "\n=================================== General settings ===================================\n\n"
			  << "- Propagation Frame / Time     : GCRS / TT"
			  << "\n- Input / Output Time          : UTC"
			  << "\n- Target ID                    : " << sat_id
			  << "\n- Mass [kg]                    : " << std::defaultfloat << std::setprecision(14) << sat_mass
			  << "\n\n- Start epoch [UTC]            : " << std::fixed << std::setprecision(12) << t0_mjd_utc
			  << "\n- End epoch [UTC]              : " << TT2UTC(opn_setting.t0_MJD_TT + opn_length_sec / 86400.0)
			  << "\n- Duration [s]                 : " << std::showpos << std::setprecision(2) << opn_length_sec << std::noshowpos
			  << "\n\n";

		// Keplerian 两行为 [a,e,i]、[Ω,ω,M]；Cartesian 两行为位置、速度。
		// 两行共用标签宽度和数值列宽，使左右方括号、逗号及各分量严格对齐。
		for (int row = 0; row < 2; ++row)
		{
			if (row == 0)
			{
				setup << std::left << std::setw(label_width) << "- Initial Keplerian [km, deg]" << ": ";
			}
			else
			{
				setup << std::string(label_width + 2, ' ');
			}

			setup << "[ ";

			for (int j = 0; j < 3; ++j)
			{
				if (j > 0) setup << ", ";

				setup << std::fixed << std::setprecision(6) << std::right << std::setw(state_width) << initial_kep_values[row * 3 + j];
			}

			setup << " ]\n";
		}

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

		// 报告保留完整设置；控制台开关只决定模型和积分分区是否在屏幕上显示。
		const string general_settings = setup.str();

		if (display_model)
		{
			opn_model.Display(tracer, setup);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		if (display_setting)
		{
			opn_setting.Display(tracer, integrator_id, setup);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// 设置文字直接显示到控制台，任务报告使用相同的内容。
		std::cerr << setup.str() << '\n';

		// =========================== 06 执行传播并保存最终状态副本 ===========================

		const string integrator_name = integrator_id == 101 ? "adaptive RKF78"                               : integrator_id == 102 ? "fixed-step RKF78"
			: integrator_id == 201 ? "Adams PE"                               : "Adams PECE";
		tracer.LogEvent(__FILE__, __func__, __LINE__, "Starting NSG2 propagation with " + integrator_name + ".", EventType::MESSAGE, FORCE_ALERT_ON);

		if (display_progress > 0) std::cerr << "\n================================= Propagation progress =================================\n\n";

		// 传播函数负责进度和积分警告；计时区间只覆盖本次实际传播调用。
		const auto propagation_start = std::chrono::steady_clock::now();

		if (integrator_id == 201 || integrator_id == 202)
		{
			OrbitPropagation_Adams_NSG2(tracer, opn_job, record_MJD_UTC, record_nsg2);
		}
		else
		{
			OrbitPropagation_RKF78_NSG2(tracer, opn_job, record_MJD_UTC, record_nsg2);
		}

		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		const double propagation_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - propagation_start).count();

		// 进度逐行连续显示；最后一条进度之后留一行，再显示传播完成事件。
		if (display_progress > 0 && record_nsg2.size() > 1) std::cerr << '\n';

		tracer.LogEvent(__FILE__, __func__, __LINE__, "Requested final epoch reached.", EventType::MESSAGE, FORCE_ALERT_ON);

		// 科学文件输出随后可能原地恢复 record_nsg2 的单位，因此先从归一化末状态取得摘要副本。
		const NonSingular2 final_nsg2_norm = record_nsg2.back();
		Cartesian final_cart_km;
		Keplerian final_kep_km;
		final_nsg2_norm.ToCartesian(tracer, final_cart_km);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		final_nsg2_norm.ToKeplerian(tracer, final_kep_km);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		final_cart_km.Restore(tracer, ae / 1000.0, UNITTIME);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		final_kep_km.Restore(tracer, ae / 1000.0);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		double final_cart_values[6], final_kep_values[6];
		final_cart_km.ToArray(final_cart_values);
		final_kep_km.ToArray(final_kep_values);

		// =========================== 07 按配置转换并写出科学星历 ===========================
		// 写文件事件在发生时记录；路径消息待全部输出成功后集中显示。
		std::ostringstream file_messages;

		vector<Cartesian> record_cart;
		vector<Keplerian> record_kep;
		// 每列显式留空格；固定最小列宽不足时也不会把相邻数字连在一起。
		// 如需要cart
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


				// 转换状态量以及单位
				for (int k = 0; k < record_nsg2.size(); k++)
				{
					Cartesian temp_cart;
					record_nsg2[k].ToCartesian(tracer, temp_cart);
					if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
					// 转换单位
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
					record_cart.push_back(temp_cart);
				}

				// 最终输出
				double cartArr[6];
				for (int k = 0; k < record_cart.size(); k++)
				{
					record_cart[k].ToArray(cartArr);
					writeDouble(ofs_ephCart, record_MJD_UTC[k]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[0]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[1]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[2]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[3]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[4]); ofs_ephCart << ' ';
					writeDouble(ofs_ephCart, cartArr[5]); ofs_ephCart << ' ';
					ofs_ephCart << endl;
				}
				ofs_ephCart.close();
				if (!ofs_ephCart) throw std::runtime_error("Failed while writing an ephemeris output.");
			}


		}

		// 如需要kep
		if (!fn_ephKep.empty())
		{
			ofstream ofs_ephKep(fn_ephKep.c_str());
			if (!ofs_ephKep)
			{
				msg = "Cannot open file_eph_Keplerian: " + fn_ephKep + " to write!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			else
			{
				msg = "Writing Keplerian ephemeris to: " + fn_ephKep;
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
				file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';

				// 转换状态量以及单位
				for (int k = 0; k < record_nsg2.size(); k++)
				{
					Keplerian temp_kep;
					record_nsg2[k].ToKeplerian(tracer, temp_kep);
					if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
					// 转换单位
					switch (eph_normalize)
					{
						case STATE_UNIT_NORMALIZE:
							break;
						case STATE_UNIT_KM_SEC_DEGREE:
							temp_kep.Restore(tracer, ae / 1000.0);
							if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
							break;
						case STATE_UNIT_M_SEC_DEGREE:
							temp_kep.Restore(tracer, ae);
							if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
							break;
						default:
							msg = "Unclear ephemeris unit! Use normalized unit in default!";
							tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
							break;

					}
					record_kep.push_back(temp_kep);
				}

				// 最终输出
				double KepArr[6];
				for (int k = 0; k < record_kep.size(); k++)
				{
					record_kep[k].ToArray(KepArr);
					writeDouble(ofs_ephKep, record_MJD_UTC[k]); ofs_ephKep << ' ';
					writeDouble(ofs_ephKep, KepArr[0]); ofs_ephKep << ' ';
					writeDouble(ofs_ephKep, KepArr[1]); ofs_ephKep << ' ';
					writeDouble(ofs_ephKep, KepArr[2]); ofs_ephKep << ' ';
					writeDouble(ofs_ephKep, KepArr[3]); ofs_ephKep << ' ';
					writeDouble(ofs_ephKep, KepArr[4]); ofs_ephKep << ' ';
					writeDouble(ofs_ephKep, KepArr[5]); ofs_ephKep << ' ';
					ofs_ephKep << endl;
				}
				ofs_ephKep.close();
				if (!ofs_ephKep) throw std::runtime_error("Failed while writing an ephemeris output.");
			}
		}

		// 如需要nsg2
		if (!fn_ephNSG2.empty())
		{
			ofstream ofs_ephNsg2(fn_ephNSG2.c_str());
			if (!ofs_ephNsg2)
			{
				msg = "Cannot open file_eph_NonSingular2: " + fn_ephNSG2 + " to write!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
			else
			{
				msg = "Writing NonSingular2 ephemeris to: " + fn_ephNSG2;
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
				file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';

				// 转换状态量以及单位
				for (int k = 0; k < record_nsg2.size(); k++)
				{
					// 转换单位
					switch (eph_normalize)
					{
						case STATE_UNIT_NORMALIZE:
							break;
						case STATE_UNIT_KM_SEC_DEGREE:
							record_nsg2[k].Restore(tracer, ae / 1000.0);
							if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
							break;
						case STATE_UNIT_M_SEC_DEGREE:
							record_nsg2[k].Restore(tracer, ae);
							if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
							break;
						default:
							msg = "Unclear ephemeris unit! Use normalized unit in default!";
							tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::WARNING);
							break;

					}
				}

				// 最终输出
				double nsg2Arr[6];
				for (int k = 0; k < record_nsg2.size(); k++)
				{
					record_nsg2[k].ToArray(nsg2Arr);
					writeDouble(ofs_ephNsg2, record_MJD_UTC[k]); ofs_ephNsg2 << ' ';
					writeDouble(ofs_ephNsg2, nsg2Arr[0]); ofs_ephNsg2 << ' ';
					writeDouble(ofs_ephNsg2, nsg2Arr[1]); ofs_ephNsg2 << ' ';
					writeDouble(ofs_ephNsg2, nsg2Arr[2]); ofs_ephNsg2 << ' ';
					writeDouble(ofs_ephNsg2, nsg2Arr[3]); ofs_ephNsg2 << ' ';
					writeDouble(ofs_ephNsg2, nsg2Arr[4]); ofs_ephNsg2 << ' ';
					writeDouble(ofs_ephNsg2, nsg2Arr[5]); ofs_ephNsg2 << ' ';
					ofs_ephNsg2 << endl;
				}
				ofs_ephNsg2.close();
				if (!ofs_ephNsg2) throw std::runtime_error("Failed while writing an ephemeris output.");
			}
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
				<< "\n- Final epoch [UTC]            : " << std::fixed << std::setprecision(12) << record_MJD_UTC.back()
				<< "\n- Propagated duration [s]      : " << std::showpos << std::setprecision(2)
				<< (UTC2TT(record_MJD_UTC.back()) - opn_setting.t0_MJD_TT) * 86400.0 << std::noshowpos
				<< "\n- Ephemeris records            : " << record_nsg2.size() << "\n\n";

		for (int row = 0; row < 2; ++row)
		{
			if (row == 0)
			{
				summary << std::left << std::setw(label_width) << "- Final Keplerian [km, deg]" << ": ";
			}
			else
			{
				summary << std::string(label_width + 2, ' ');
			}

			summary << "[ ";

			for (int j = 0; j < 3; ++j)
			{
				if (j > 0) summary << ", ";

				summary << std::fixed << std::setprecision(6) << std::right << std::setw(state_width) << final_kep_values[row * 3 + j];
			}

			summary << " ]\n";
		}

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

			report << "NSG2 orbit propagation report\n"
				   << "Program       : " << program_name << " | Version: v" << ASTER_PROGRAM_VERSION
				   << "\nJob ID        : " << jobId << "\nStarted at    : " << timeStr
				   << "\nInput setting : " << config_path.string() << '\n'
				   << general_settings;
			opn_model.Display(tracer, report);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			opn_setting.Display(tracer, integrator_id, report);
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
		const string completed = "NSG2 propagation job completed successfully.";
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
