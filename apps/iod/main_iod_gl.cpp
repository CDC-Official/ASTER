/**
 * @file main_iod_gl.cpp
 * @brief 根据同一目标的单站或多站 RA/Dec 观测，使用 GL 方法确定参考历元的轨道。
 *
 * @details
 * 程序先读取配置、观测和测站资料，检查本次选用的观测，再计算各观测时刻的测站位置。
 * 随后调用 GL 定轨核心，求解参考历元的位置和速度，并输出轨道解及选定的结果文件。
 * 简单迭代求解无约束的基本方程；微分改正支持方向、圆轨道及方向与偏心率等约束组合。
 *
 * @par 输入资料
 * 配置文件采用 SimpleInput 格式，可通过 -I、-input 或 --input 指定；省略配置参数时显示用法。
 * 观测文件每行依次为 MJD_UTC、观测类型、目标号、站号、RA(deg)、Dec(deg)，本程序使用类型 1。
 * 测站文件每行以站号和坐标类型开头：类型 1 后接大地纬度 phi(deg)、经度 lambda(deg)、高度(m)；
 * 类型 2 后接 ITRS 直角坐标 X(m)、Y(m)、Z(m)。观测通过站号与测站资料关联。
 * EOP 文件提供测站坐标转换所需的地球定向参数。配置中的相对路径以配置文件所在目录为基准。
 * 可选 file_initial_orbit 使用 GL 轨道解格式，提供同一目标在参考历元的 GCRS 初值；
 * 省略时采用默认线性初值，提供时直接从指定状态开始迭代。
 *
 * @par 时间、坐标与单位
 * 观测时刻和参考历元均采用 UTC MJD；定轨使用 GCRS 下的地心位置与速度。
 * 角度在计算前转换为 rad，位置、速度和相对时间在核心计算中使用项目定义的归一化单位。
 * 每条观测的测站位置按该条观测时刻计算，EOP 统一取参考历元的值。
 *
 * @par 输出结果
 * 轨道解文件采用 SimpleInput 格式，保存参考历元、轨道状态、测角残差 RMS 和求解状态标识。
 * 轨道状态的输出单位由 output_state_units 指定；残差文件、运行报告和运行日志按配置选择输出。
 * 残差文件每行依次为 MJD_UTC、观测类型、目标号、站号、东西方向残差和南北方向残差；
 * 两个残差均取 O-C，单位为 rad，文件不含表头。
 * 任务报告保存设置、逐轮迭代和结果；程序日志保存配置原文、文件路径、事件及错误调用链。
 */
#include "math_backend.hpp"
#include "core_iod_solver.hpp"

#include "io_settings.hpp"
#include "model_observation.hpp"
#include "frame_definition.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <ctime>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <cstdio>

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

#include "core_iod.hpp"
#include "time_transform.hpp"

namespace fs = std::filesystem;

using namespace std;
using namespace Eigen;
using namespace iod;
using namespace bastion;
using namespace orbit;
using namespace io;

// 观测记录、观测集合和测站集合用于保存输入资料。
using observation::ObservationRecord;
using observation::Observations;
using observation::Observers;

/**
 * @brief 读取单目标测角资料，执行 GL 初轨确定并输出结果和任务报告。
 * @param[in] argc 命令行参数个数，包含程序路径。
 * @param[in] argv 命令行参数；-I、-input 或 --input 后接 SimpleInput 配置文件路径。
 * @return 正常完成返回 0；文件、输入或数值错误由 tracer 报告并以失败状态退出。
 * @details General settings 和 Observation summary 在本函数组织，IOD 设置和最终轨道
 * 分别调用现有任务、结果对象的显示成员函数。科学结果文件的字段和精度直接在本函数中设置。
 */
int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_iod_gl " << ASTER_PROGRAM_VERSION << '\n';
        std::cout << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_iod_gl -I <setting.txt>\n";
        return 2;
    }

	// =================================== 01 初始化任务与事件追踪 ===================================

	Overwatch tracer(__FILE__, __func__, __LINE__);
	std::ofstream program_log;

	try
	{
		const auto wall_start = std::chrono::steady_clock::now();
		const std::string program_name = fs::path(argv[0]).filename().string();

		// 用当前本地时间生成任务编号和启动时间，便于在输出信息中识别本次运行。
		std::time_t timeNow = std::time(nullptr); // 获取当前时间
		std::tm* tmNow = std::localtime(&timeNow);

		char timeStr[20]; // 创建用于存储时间格式化结果的字符串
		char jobId[15];

		std::strftime(jobId, sizeof(jobId), "%Y%m%d%H%M%S", tmNow); // 生成 jobId 格式：YYYYMMDDHHMMSS
		std::strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", tmNow);     // 生成日志消息格式：YYYY-MM-DD HH:MM:SS

		std::stringstream jobMsg; // 构建日志消息
		jobMsg << "GL initial orbit determination job started at " << timeStr
			   << " (jobId = " << jobId << ")";
		tracer.LogEvent(__FILE__, __func__, __LINE__, jobMsg.str(), bastion::EventType::MESSAGE, FORCE_ALERT_ON);

		const std::string version_message = "Program: " + program_name + " | Version: v" + ASTER_PROGRAM_VERSION;
		tracer.LogEvent(__FILE__, __func__, __LINE__, version_message, EventType::MESSAGE, FORCE_ALERT_ON);

		// =================================== 02 解析命令行参数 ===================================

		// 设置默认配置文件；使用 -I、-input 或 --input 可以指定其他配置文件。
		string fn_iod_input = "";

		// 逐个读取命令行选项，处理帮助请求或取得配置文件路径。
		for (int i = 1; i < argc; ++i)
		{
			string argument = argv[i];

			// 收到帮助选项时，显示命令行用法并结束程序。
			if (argument == "-h" || argument == "--help")
			{
				cout << "Usage: aster_iod_gl [-I|-input|--input <setting>]\n"
					 << "Output paths come from setting; prepare their parent directories before "
						"running.\n";
				tracer.CheckOut(__func__);
				return 0;
			}

			// 其余选项用于指定配置文件，检查选项名称是否有效。
			if (argument != "-I" && argument != "-input" && argument != "--input")
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown argument: " + argument, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// -I、-input 和 --input 后面需要填写配置文件路径。
			if (i + 1 >= argc || argv[i + 1][0] == '\0' || argv[i + 1][0] == '-')
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Missing value for " + argument, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			fn_iod_input = argv[++i];
		}

		// =================================== 03 读取配置文件并准备运行日志 ===================================

		string msg;              // 用于记录事件
		IodJob thisJob{};        // 保存本次定轨的设置、观测资料和计算结果。
		SimpleInput si_iodInput; // 读取输入文件

		fn_iod_input = fs::absolute(fn_iod_input).lexically_normal().string();
		msg = "Reading and analyzing input from: " + fn_iod_input;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		// 先读入配置记录；后续按用途逐项读取文件路径和求解参数。
		si_iodInput.ReadSimpleInput(tracer, fn_iod_input);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_programlog：运行日志路径，保存配置、事件及出错时的调用链。
		// 填写尚不存在的文件名后启用；省略或留空时关闭文件日志。
		std::string fn_programlog;

		if (si_iodInput.HasKeyword("file_programlog"))
		{
			si_iodInput.GetString(tracer, "file_programlog", fn_programlog);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// 指定了日志路径时，创建日志并记录配置内容及各文件的完整路径。
		if (!fn_programlog.empty())
		{
			const std::filesystem::path log_config = std::filesystem::absolute(fn_iod_input);
			std::filesystem::path log_path = fn_programlog;

			if (log_path.is_relative()) log_path = log_config.parent_path() / log_path;

			log_path = std::filesystem::weakly_canonical(log_path);
			fn_programlog = log_path.string();

			// 检查日志文件是否已经存在；已有文件需要另取文件名。
			if (std::filesystem::exists(log_path))
			{
				msg = "Program log already exists: " + log_path.string();
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR, FORCE_ALERT_ON);
				tracer.EmergencyExit(program_log);
			}

			// 将日志路径与配置中的其他文件路径逐一比较，检查是否重名。
			for (std::size_t record_index = 1; record_index <= si_iodInput.Count(); ++record_index)
			{
				io::SimpleRecord path_record;
				si_iodInput.GetRecord(tracer, record_index, path_record);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

				if (path_record.kw.compare(0, 5, "file_") != 0 || path_record.kw == "file_programlog" || path_record.va.empty()) continue;

				std::filesystem::path other_path = path_record.va;

				if (other_path.is_relative()) other_path = log_config.parent_path() / other_path;

				if (std::filesystem::weakly_canonical(other_path) == log_path)
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "Program log conflicts with " + path_record.kw, EventType::ERROR, FORCE_ALERT_ON);
					tracer.EmergencyExit(program_log);
				}
			}

			// 创建日志文件，记录启动时间、程序路径和配置文件路径。
			program_log.open(log_path);
			if (!program_log)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open program log: " + log_path.string(), EventType::ERROR, FORCE_ALERT_ON);
				tracer.EmergencyExit(program_log);
			}

			const std::time_t log_started = std::time(nullptr);
			program_log << "Program log\nStarted: " << std::ctime(&log_started)
						<< "Executable: " << argv[0] << "\nConfig: " << log_config.string()
						<< "\nLog: " << log_path.string() << "\n\n[Input settings]\n";

			// 将本次使用的配置原文保存到日志，供日后核对算例设置。
			std::ifstream log_input(log_config);

			if (!log_input)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot read settings for program log.", EventType::ERROR, FORCE_ALERT_ON);
				tracer.EmergencyExit(program_log);
			}

			program_log << log_input.rdbuf() << "\n[Resolved file paths]\n";

			// 把配置中的文件路径展开为绝对路径后写入日志。
			for (std::size_t record_index = 1; record_index <= si_iodInput.Count(); ++record_index)
			{
				io::SimpleRecord path_record;
				si_iodInput.GetRecord(tracer, record_index, path_record);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

				if (path_record.kw.compare(0, 5, "file_") != 0 || path_record.va.empty()) continue;

				std::filesystem::path recorded_path = path_record.va;

				if (recorded_path.is_relative()) recorded_path = log_config.parent_path() / recorded_path;

				program_log << path_record.kw << " = "
							<< std::filesystem::weakly_canonical(recorded_path).string() << '\n';
			}

			// 将日志缓冲区写入文件，并检查写入是否成功。
			program_log.flush();
			if (!program_log)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing program log.", EventType::ERROR, FORCE_ALERT_ON);
				tracer.EmergencyExit(program_log);
			}
		}

		// =================================== 04 读取求解与输出设置 ===================================

		// 保存配置文件所在目录，用于解析输入和输出文件的相对路径。
		const fs::path configdir = fs::absolute(fn_iod_input).parent_path();

		string epochSetting; // 参考历元设定：First、Mid、Last 或明确的 UTC MJD。

		double noiseLevelAs; // 东西、南北两个方向共同的误差标准差，arcsec；用于似然评价。

		// 【必选】file_observation：地面测站 RA/Dec 观测文件路径，不能留空。
		si_iodInput.GetString(tracer, "file_observation", thisJob.fn_obs);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】file_EOP：地球定向参数文件路径，用于将测站位置从 ITRS 转换到 GCRS。
		si_iodInput.GetString(tracer, "file_EOP", thisJob.fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】file_orbit_solution_GL：参考历元轨道解的输出路径，不能留空。
		si_iodInput.GetString(tracer, "file_orbit_solution_GL", thisJob.fn_solution);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_report：任务报告路径，保存设置、观测概览、逐轮迭代和最终结果；省略或留空时关闭输出。
		si_iodInput.GetString(tracer, "file_report", thisJob.fn_report);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_residual：逐条观测残差的输出路径；省略或留空时关闭输出。
		// 残差为东西和南北方向的 O-C，单位 rad。
		si_iodInput.GetString(tracer, "file_residual", thisJob.fn_residual);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_initial_orbit：参考历元的 GCRS 初值文件，采用 GL 轨道解的 SimpleInput 格式。
		// 省略或留空时使用默认线性初值；填写时使用文件中指定的 Cartesian 或 Keplerian 状态。
		std::string fn_initial_orbit;

		if (si_iodInput.HasKeyword("file_initial_orbit"))
		{
			si_iodInput.GetString(tracer, "file_initial_orbit", fn_initial_orbit);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// 【必选】iteration_method_GL：1 为简单迭代，2 为微分改正。
		// 简单迭代用于无约束的基本方程；加入约束时选 2。
		si_iodInput.GetInteger(tracer, "iteration_method_GL", thisJob.iterMethod);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】constraint_type_GL：0 为无约束，1 为方向，2 为圆轨道，3 为方向＋圆轨道，4 为方向＋偏心率。
		// 当前读取器在缺项时返回 0，程序按无约束处理。
		si_iodInput.GetInteger(tracer, "constraint_type_GL", thisJob.constraintType);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 检查迭代方法的有效性
		if (thisJob.iterMethod != IOD_ITERATION_SIMPLE && thisJob.iterMethod != IOD_ITERATION_NEWTON)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "iteration_method_GL must be 1 or 2.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 检查约束类型的有效性
		if (thisJob.constraintType < IOD_CONSTRAINT_NONE || thisJob.constraintType > IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown constraint_type_GL.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 简单迭代求解 3N 个基本方程；加入约束时需要选择微分改正。
		if ((thisJob.iterMethod == IOD_ITERATION_SIMPLE) && (thisJob.constraintType != IOD_CONSTRAINT_NONE))
		{
			tracer.LogEvent( __FILE__, __func__, __LINE__, "Simple iteration requires constraint_type_GL=0; use method 2 for constraints.",
				EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 【条件可选】constraint_eccentricity：constraint_type_GL 为 4 时读取的指定偏心率，无量纲，范围为 [0,1)。
		// 当前读取器在缺项时返回 0，对应圆轨道的偏心率；其他约束类型下不读取此项。
		if (thisJob.constraintType == IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY)
		{
			si_iodInput.GetDouble(tracer, "constraint_eccentricity", thisJob.constraintEcc);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			if (!std::isfinite(thisJob.constraintEcc) || thisJob.constraintEcc < 0.0 || thisJob.constraintEcc >= 1.0)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "constraint_eccentricity must be in [0,1).", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 【必选】iteration_maxround：最大迭代次数，取正整数。
		si_iodInput.GetInteger(tracer, "iteration_maxround", thisJob.iterCountMax);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】iteration_state_tolerance：六个归一化状态改正分量的最大绝对值阈值，取有限正数。
		si_iodInput.GetDouble(tracer, "iteration_state_tolerance", thisJob.iterTOL);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】display_iterations：0 为关闭屏幕迭代信息，1 为显示；缺项时读取器返回 0。
		// 每轮计算记录始终保存在任务结果中。
		int print_iterInfo;
		si_iodInput.GetInteger(tracer, "display_iterations", print_iterInfo);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		if (print_iterInfo == IOD_PRINT_ITERATION_INFO_OFF)
		{
			thisJob.printIterInfo = false;
		}
		else if (print_iterInfo == IOD_PRINT_ITERATION_INFO_ON)
		{
			thisJob.printIterInfo = true;
		}
		else
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "display_iterations must be 0 or 1.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 【必选】reference_epoch：待求轨道所属的参考历元，可填 First、Mid、Last 或明确的 UTC MJD。
		// First、Mid、Last 分别取选用观测的首条、中间一条和末条记录的时刻。
		si_iodInput.GetString(tracer, "reference_epoch", epochSetting);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】dynamics_model_GL：1 为二体，2 为 J2，3 为 J2＋太阳，4 为 J2＋月球，5 为 J2＋日月。
		si_iodInput.GetInteger(tracer, "dynamics_model_GL", thisJob.dynamicsModel);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		if (thisJob.dynamicsModel < IOD_DYNAMICS_TWOBODY || thisJob.dynamicsModel > IOD_DYNAMICS_J2_SOLAR_LUNAR)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "dynamics_model_GL must be between 1 and 5.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 【必选】output_state_units：1 为归一化单位/rad，2 为 km、km/s、deg，3 为 m、m/s、deg。
		si_iodInput.GetInteger(tracer, "output_state_units", thisJob.solutionUnit);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】observations_to_use：All 表示使用全部观测，整数 N 表示使用文件中的前 N 条；N 至少为 3。
		std::string observationCount;
		si_iodInput.GetString(tracer, "observations_to_use", observationCount);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】noise_sigma_arcsecond：东西、南北独立高斯误差共同的标准差，单位 arcsec，取有限正数。
		// 读入后换算为 rad，保存到本次任务的观测设置中。
		si_iodInput.GetDouble(tracer, "noise_sigma_arcsecond", noiseLevelAs);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		if (!isfinite(noiseLevelAs) || noiseLevelAs <= 0.0)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "noise_sigma_arcsecond must be positive and finite.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		thisJob.obs.sigma_RaDec = noiseLevelAs * ARCSEC_TO_RAD; // 转化为弧度单位

		// 最大迭代次数和状态改正量的收敛阈值均应为正数。
		if (thisJob.iterCountMax <= 0 || !isfinite(thisJob.iterTOL) || thisJob.iterTOL <= 0.0)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "iteration_maxround and iteration_state_tolerance must be positive.",
				EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 输出单位可选归一化单位、km/deg 或 m/deg。
		if (thisJob.solutionUnit < IOD_SOLUTION_NORMALIZE || thisJob.solutionUnit > IOD_SOLUTION_MDEG)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "output_state_units must be 1, 2 or 3.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 检查轨道解的输出路径和参考历元是否已填写。
		if (thisJob.fn_solution.empty())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "file_orbit_solution_GL is required.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		if (epochSetting.empty())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "reference_epoch is required.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// =================================== 05 解析文件路径与观测选用范围 ===================================

		// 【必选】file_observer：地面测站表路径，每条选用观测的站号都应能在表中找到。
		// 各文件的相对路径均以配置文件所在目录为基准。
		string stationFile;
		si_iodInput.GetString(tracer, "file_observer", stationFile);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 将各文件路径相对于配置目录展开，并消去路径中的 . 和 ..。
		if (!thisJob.fn_obs.empty()) thisJob.fn_obs = (configdir / thisJob.fn_obs).lexically_normal().string();

		if (!thisJob.fn_EOP.empty()) thisJob.fn_EOP = (configdir / thisJob.fn_EOP).lexically_normal().string();

		if (!stationFile.empty()) stationFile = (configdir / stationFile).lexically_normal().string();

		if (!thisJob.fn_solution.empty()) thisJob.fn_solution = (configdir / thisJob.fn_solution).lexically_normal().string();

		if (!thisJob.fn_report.empty()) thisJob.fn_report = (configdir / thisJob.fn_report).lexically_normal().string();

		if (!thisJob.fn_residual.empty()) thisJob.fn_residual = (configdir / thisJob.fn_residual).lexically_normal().string();

		if (!fn_initial_orbit.empty()) fn_initial_orbit = (configdir / fn_initial_orbit).lexically_normal().string();

		// All 使用全部记录；整数 N 使用文件顺序中的前 N 条，且 N 至少为三。
		size_t readCount = 0;

		if (observationCount != "All")
		{
			if (observationCount.empty() || observationCount.find_first_not_of("0123456789") != string::npos)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "observations_to_use must be All or an integer >= 3.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 将条数转换为整数；转换失败表示数值超出了 size_t 的取值范围。
			istringstream countInput(observationCount);
			if (!(countInput >> readCount))
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "observations_to_use must be All or an integer >= 3.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// GL 初轨确定至少需要三条观测。
			if (readCount < 3)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "observations_to_use must be All or an integer >= 3.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// =================================== 06 读取原始观测并关联测站 ===================================

		// 完整读入观测和测站资料，并按站号为每条观测关联相应的测站。
		// 读入的资料保留文件中的单位，供后续检查和坐标转换使用。

		// 在读取事件中给出完整路径，随后检查读取器记录的错误。
		msg = "Reading observations from: " + thisJob.fn_obs;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
		Observations iod_observations;
		iod_observations.ReadFile(tracer, thisJob.fn_obs);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 读取地面测站表，原始站址在后续使用阶段转换为 GCRS。
		msg = "Reading observers from: " + stationFile;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
		Observers iod_stations;
		iod_stations.ReadFile(tracer, stationFile);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 按观测记录中的站号查找测站，将对应的站址资料关联到每条观测。
		iod_observations.AttachObservers(tracer, iod_stations);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		const std::vector<ObservationRecord>& records = iod_observations.GetRecords();

		// =================================== 07 检查本次选用的观测记录 ===================================

		// readCount 为零时使用全部记录，否则使用前 readCount 条记录。
		// 本次使用的条数应至少为三，且应在文件实际提供的条数以内。
		const std::size_t used_count = readCount == 0 ? records.size() : readCount;

		if (used_count < 3 || used_count > records.size())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Too few observations for observations_to_use.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 由于后续算法使用 int 作为观测点数，若超过 int 最大值则直接报错退出。
		if (used_count > static_cast<std::size_t>(numeric_limits<int>::max()))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Too many observation records.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 以第一条观测的目标号确定本次定轨的目标。
		thisJob.objectId = records[0].object_id;
		std::vector<std::string> station_ids; // 按首次出现顺序保存本次使用的站号。

		// 逐条检查选用观测的目标号、类型、角度和时间顺序，并统计参与观测的测站。
		for (std::size_t i = 0; i < used_count; ++i)
		{
			const ObservationRecord& record = records[i];
			const std::string location = "Observation record " + std::to_string(i + 1) + ": ";

			// 目前仅支持地面测站的 RA/Dec 观测，类型 1
			if (record.type != 1)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "GL supports ground RA/Dec (type 1) only.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 各测站提供的观测应属于同一个目标。
			if (record.object_id != thisJob.objectId)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "GL requires one target; ground station IDs may differ.",
					EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 将首次出现的站号加入列表，用于汇总本次参与观测的测站。
			if (std::find(station_ids.begin(), station_ids.end(), record.observer_id) == station_ids.end())
				station_ids.push_back(record.observer_id);

			// 赤纬的有效范围为 −90° 到 90°。
			if (std::abs(record.dec_deg) > 90.0)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "Declination outside [-90,90] deg.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 观测应按时间先后排列；同一时刻可以有多条观测。
			if (i > 0 && record.mjd_utc < records[i - 1].mjd_utc)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "Observation epochs must be nondecreasing UTC MJD.",
					EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 每条观测都需要对应的站址，才能计算观测者在 GCRS 中的位置。
			if (!record.is_observer_attached)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "Unknown ground station: " + record.observer_id, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// =================================== 08 准备观测数组并转换测角单位 ===================================

		// 选用记录检查通过后，为 N 个观测历元、角度、测站位置和残差分配空间。
		// 后续将选用记录中的时间、角度和测站位置填入这些数组。
		const int num = static_cast<int>(used_count);
		thisJob.obs.type = OBS_TYPE_RADEC;
		thisJob.obs.number = num;
		thisJob.obs.time_MJDUTC.resize(num);
		thisJob.obs.tau.resize(num);
		thisJob.obs.ra_rad.resize(num);
		thisJob.obs.dec_rad.resize(num);

		thisJob.observer.positions_gcrs_norm_emat.resize(3, num);

		thisJob.solution.resRa.resize(num);
		thisJob.solution.resDec.resize(num);

		// 取前 num 条记录，将 UTC MJD 与弧度制测角填入算法数组。
		for (int k = 0; k < num; ++k)
		{
			const ObservationRecord& record = records[k];
			thisJob.obs.time_MJDUTC(k) = record.mjd_utc;

			// 将赤经、赤纬从 deg 换算为 rad；赤经为负值时加上 2π。
			double ra_rad = record.ra_deg * DEGREE_TO_RAD;

			if (ra_rad < 0.0) ra_rad += 2 * pi;

			thisJob.obs.ra_rad(k) = ra_rad;
			thisJob.obs.dec_rad(k) = record.dec_deg * DEGREE_TO_RAD;
		}

		// =================================== 09 确定参考历元并读取 EOP ===================================

		// 先确定待估状态所属的参考历元，再据此计算各次观测的相对时间。
		// Mid 取中间一条观测的时刻，不是首尾时刻的平均值。
		if (epochSetting == "Mid")
		{
			thisJob.epochMJDUTC = thisJob.obs.time_MJDUTC(num / 2);
		}
		else if (epochSetting == "First")
		{
			thisJob.epochMJDUTC = thisJob.obs.time_MJDUTC(0);
		}
		else if (epochSetting == "Last")
		{
			thisJob.epochMJDUTC = thisJob.obs.time_MJDUTC(num - 1);
		}
		else
		{
			si_iodInput.GetDouble(tracer, "reference_epoch", thisJob.epochMJDUTC);
		}
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 参考历元必须是有限数值，才能计算相对时间和测站位置。
		if (!isfinite(thisJob.epochMJDUTC))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "reference_epoch must be finite.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 观测时刻减去定轨历元，由日换算为归一化时间。
		thisJob.obs.tau = (thisJob.obs.time_MJDUTC - thisJob.epochMJDUTC) * DAY_TO_NORMALIZED_TIME;

		// 读取地球定向参数（EOP），并取得参考历元的参数用于测站坐标转换。
		msg = "Reading EOP from: " + thisJob.fn_EOP;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
		EOP thisEOP;
		thisEOP.Read_EOP(tracer, thisJob.fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		aEOP EOP_at_epoch;
		thisEOP.Get_EOP(tracer, thisJob.epochMJDUTC, EOP_at_epoch, 1);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// =================================== 10 逐条计算观测时刻的测站位置 ===================================

		// 根据每条观测关联的测站，先计算归一化 ITRS 站址，再转换到该次观测时刻的 GCRS。
		// 各次转换使用各自的观测时刻，EOP 统一取参考历元的值。
		for (int k = 0; k < num; ++k)
		{
			const observation::ObserverRecord& ground_observer = records[k].ground_observer;
			frame::ITRS site_itrs;

			if (ground_observer.coordinate_type == 1)
			{
				// 大地站址为高度(m)、lambda(deg)、phi(deg)；先转换为归一化 ITRS XYZ。
				frame::Geodetic site_geodetic(ground_observer.geodetic_height_m, ground_observer.geodetic_lambda_deg,
					ground_observer.geodetic_phi_deg, false);
				site_geodetic.Normalize(tracer, ae);
				site_geodetic.ToITRS(tracer, site_itrs);
			}
			else
			{
				// 类型 2 直接提供 ITRS XYZ(m)，以 ae 归一化后用于地固到天球的变换。
				site_itrs.Set(ground_observer.position_itrs_m[0], ground_observer.position_itrs_m[1],
					ground_observer.position_itrs_m[2], false);
				site_itrs.Normalize(tracer, ae);
			}

			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 在本次观测时刻将站址从 ITRS 转换到 GCRS，位置仍以地球赤道半径 ae 归一化。
			frame::GCRS site_gcrs;
			site_itrs.ToGCRS(tracer, thisJob.obs.time_MJDUTC(k), EOP_at_epoch, site_gcrs);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 第 k 列保存第 k 条观测对应的测站位置，供该条观测的方程使用。
			double position[3];
			site_gcrs.ToArray(position);
			thisJob.observer.positions_gcrs_norm_emat.col(k) = Vector3d(position[0], position[1], position[2]);

			// angle_norm 的计算需要第一条观测所属测站在参考历元的位置，因此另算并保存该位置。
			if (k == 0)
			{
				site_itrs.ToGCRS(tracer, thisJob.epochMJDUTC, EOP_at_epoch, site_gcrs);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

				site_gcrs.ToArray(position);
				thisJob.observer.position_gcrs_epoch_norm = Vector3d(position[0], position[1], position[2]);
			}
		}

		// =================================== 11 设置默认初值或读取指定初轨 ===================================

		// 未提供初轨文件时，由核心的首轮线性求解生成起始状态。
		thisJob.initValue.type = IOD_INITIAL_VALUE_DEFAULT;

		// 指定初轨采用与 GL 轨道解相同的字段，状态必须属于本次目标和参考历元。
		if (!fn_initial_orbit.empty())
		{
			msg = "Reading initial orbit from: " + fn_initial_orbit;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
			SimpleInput initial_orbit;
			initial_orbit.ReadSimpleInput(tracer, fn_initial_orbit);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 【必选】object_id、state_units、reference_epoch：目标号、状态单位和 UTC 参考历元。
			std::string initial_object_id;
			initial_orbit.GetString(tracer, "object_id", initial_object_id);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			int initial_units;
			initial_orbit.GetIntegerStrict(tracer, "state_units", initial_units);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			double initial_epoch;
			initial_orbit.GetDoubleStrict(tracer, "reference_epoch", initial_epoch);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 输入初轨直接用于本次参考历元，允许 1e-10 日的文件小数舍入差。
			if (initial_object_id != thisJob.objectId || std::abs(initial_epoch - thisJob.epochMJDUTC) > 1e-10)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial orbit must match the selected target and reference epoch: " + fn_initial_orbit, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			if (initial_units < IOD_SOLUTION_NORMALIZE || initial_units > IOD_SOLUTION_MDEG)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial orbit state_units must be 1, 2 or 3: " + fn_initial_orbit, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 【必选其一】reference_Cartesian / reference_Keplerian：六个分量，同时提供时采用 Cartesian。
			// 单位码 1 为归一化/rad，2 为 km、km/s、deg，3 为 m、m/s、deg。
			const bool initial_is_cartesian = initial_orbit.HasKeyword("reference_Cartesian");
			const std::string state_keyword = initial_is_cartesian ? "reference_Cartesian" : "reference_Keplerian";

			if (!initial_orbit.HasKeyword(state_keyword))
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial orbit requires reference_Cartesian or reference_Keplerian: " + fn_initial_orbit, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			std::string state_text;
			initial_orbit.GetString(tracer, state_keyword, state_text);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			// 按文件顺序解析恰好六个有限分量，再交给已有状态类进行单位转换。
			std::array<double, 6> initial_values{};
			std::istringstream state_input(state_text);

			for (double& value : initial_values)
			{
				if (!(state_input >> value) || !std::isfinite(value))
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial orbit requires six finite state values: " + fn_initial_orbit, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}
			}

			std::string trailing;

			if (state_input >> trailing)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial orbit requires exactly six state values: " + fn_initial_orbit, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			const bool is_normalized = initial_units == IOD_SOLUTION_NORMALIZE;
			const double length_unit = initial_units == IOD_SOLUTION_KMDEG ? ae / 1000.0 : ae;

			if (initial_is_cartesian)
			{
				thisJob.initValue.initCart.Set(initial_values[0], initial_values[1], initial_values[2],
					initial_values[3], initial_values[4], initial_values[5], is_normalized);

				if (!is_normalized)
				{
					thisJob.initValue.initCart.Normalize(tracer, length_unit, UNITTIME);
					if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
				}
			}
			else
			{
				// 现有 Keplerian 转换采用椭圆轨道，先检查半长轴和偏心率的适用范围。
				if (!(initial_values[0] > 0.0) || initial_values[1] < 0.0 || initial_values[1] >= 1.0)
				{
					tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial Keplerian orbit requires a > 0 and 0 <= e < 1: " + fn_initial_orbit, EventType::ERROR);
					tracer.EmergencyExit(program_log);
				}

				thisJob.initValue.initKep.Set(initial_values[0], initial_values[1], initial_values[2],
					initial_values[3], initial_values[4], initial_values[5], is_normalized);

				if (!is_normalized)
				{
					thisJob.initValue.initKep.Normalize(tracer, length_unit);
					if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
				}

				thisJob.initValue.initKep.ToCartesian(tracer, thisJob.initValue.initCart);
				if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
			}

			thisJob.initValue.type = IOD_INITIAL_VALUE_SPECIFIED;
		}

		// =================================== 12 显示任务设置并执行 GL 定轨 ===================================

		msg = "Input checks and model preparation completed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		// 通用信息由主程序组织；历元保持输入的 UTC 日数表示，状态所属坐标系为 GCRS。
		const int label_width = 40;
		std::ostringstream general_settings;
		general_settings << "\n=================================== General settings ===================================\n\n"
						 << std::left << std::setw(label_width) << "- Reference frame" << ": GCRS\n"
						 << std::setw(label_width) << "- Input / Output Time" << ": UTC\n"
						 << std::setw(label_width) << "- Target ID" << ": " << thisJob.objectId << '\n'
						 << std::setw(label_width) << "- Reference epoch [UTC]" << ": "
						 << std::fixed << std::setprecision(12) << thisJob.epochMJDUTC << '\n';
		std::cerr << general_settings.str();

		// 统计实际参与本次定轨的记录，文件中未选用的部分不计入历元范围和类型分项。
		std::array<std::size_t, 9> type_counts{};
		std::size_t unique_epochs = 0;

		for (std::size_t k = 0; k < used_count; ++k)
		{
			++type_counts[records[k].type];

			if (k == 0 || records[k].mjd_utc != records[k - 1].mjd_utc) ++unique_epochs;
		}

		std::ostringstream observation_summary;
		observation_summary << "\n================================= Observation summary ==================================\n\n"
							<< std::left << std::setw(label_width) << "- First / Last observation [UTC]" << ": "
							<< std::fixed << std::setprecision(12) << records[0].mjd_utc << " / " << records[used_count - 1].mjd_utc
							<< '\n' << std::setw(label_width) << "- Arc span [s]" << ": " << std::setprecision(2)
							<< (records[used_count - 1].mjd_utc - records[0].mjd_utc) * 86400.0
							<< '\n' << std::setw(label_width) << "- Unique epochs" << ": " << unique_epochs
							<< "\n\n" << std::setw(label_width) << "- Observers" << ": " << station_ids.size() << " (";

		for (std::size_t k = 0; k < station_ids.size(); ++k)
		{
			if (k > 0) observation_summary << ", ";

			observation_summary << station_ids[k];
		}

		observation_summary << ")\n\n";

		// 仅在选用了部分资料时附带读入条数，使选用数量与原始文件数量清楚区分。
		if (used_count != records.size())
		{
			observation_summary << std::setw(label_width) << "- Loaded observation records" << ": " << records.size() << '\n';
		}

		observation_summary << std::setw(label_width) << "- Observation records" << ": " << used_count << '\n';

		for (int type = 1; type <= 8; ++type)
		{
			if (type_counts[type] == 0) continue;

			const std::string type_label = "  > " + std::string(observation::ObservationTypeName(type)) + " (" + std::to_string(type) + ")";
			observation_summary << std::setw(label_width) << type_label << ": " << type_counts[type] << '\n';
		}

		observation_summary << std::setw(label_width) << "- Scalar components" << ": " << 2 * used_count << '\n';
		std::cerr << observation_summary.str();

		// 求解设置由任务对象显示，同一份文字用于屏幕和任务报告。
		std::ostringstream iod_settings;
		thisJob.Display(iod_settings);
		std::cerr << iod_settings.str();

		// 执行初轨确定，迭代求解参考历元的位置和速度。
		std::ostringstream iteration_report;
		std::ostream* iteration_output = thisJob.fn_report.empty() ? nullptr : &iteration_report;
		iod_solver::InitialOrbitDetermination(tracer, thisJob, iteration_output);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// =================================== 13 检查定轨结果是否可以写出 ===================================

		// 先检查求解状态与 RMS，再检查轨道分量和逐条残差是否为有限数。
		// 数值失败时通过 tracer 报告状态并结束本次任务。
		if (thisJob.solution.flag == -2 || !std::isfinite(thisJob.solution.RMS_angle))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "IOD status: NUMERICAL FAILURE; no solution written.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 检查直角坐标状态和开普勒根数的六个分量是否均为有限数值。
		double checkedCart[6], checkedKep[6];
		thisJob.solution.cart.ToArray(checkedCart);
		thisJob.solution.kep.ToArray(checkedKep);

		for (int i = 0; i < 6; ++i)
		{
			if (!isfinite(checkedCart[i]) || !isfinite(checkedKep[i]))
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite IOD orbit.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 检查每条观测的两个残差分量是否均为有限数值。
		if (!thisJob.solution.resRa.isFinite().all() || !thisJob.solution.resDec.isFinite().all())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite IOD residuals.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// =================================== 14 写出轨道解与观测残差 ===================================

		// 写出事件先保存到 tracer，控制台统一在摘要之后显示这些路径。
		std::ostringstream file_messages;

		// 按指定单位写出参考历元的轨道，同时保存测角残差 RMS 和求解状态标识。
		// 从最终 Cartesian 状态准备两种输出表示，保留求解器中的归一化结果。
		orbit::Cartesian output_cart = thisJob.solution.cart;
		orbit::Keplerian output_kep;
		output_cart.ToKeplerian(tracer, output_kep);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 单位码 1 保留归一化量和 rad；2 恢复为 km、km/s、deg；3 恢复为 m、m/s、deg。
		if (thisJob.solutionUnit != IOD_SOLUTION_NORMALIZE)
		{
			const double length_scale = thisJob.solutionUnit == IOD_SOLUTION_KMDEG ? ae / 1000 : ae;
			output_cart.Restore(tracer, length_scale, UNITTIME);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

			output_kep.Restore(tracer, length_scale);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		double output_cart_values[6], output_kep_values[6];
		output_cart.ToArray(output_cart_values);
		output_kep.ToArray(output_kep_values);

		msg = "Writing orbit solution to: " + thisJob.fn_solution;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
		file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
		std::ofstream ofs_solution(thisJob.fn_solution);
		if (!ofs_solution)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open orbit solution file: " + thisJob.fn_solution, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// SimpleInput 轨道解：历元为 UTC MJD，状态单位由 state_units 指明，RMS 为 rad。
		// defaultfloat 下 setprecision(17) 表示 17 位有效数字；字段和精度直接在这里调整。
		ofs_solution << std::defaultfloat << std::setprecision(17)
					 << "object_id = " << thisJob.objectId
					 << "\nstate_units = " << thisJob.solutionUnit
					 << "\nreference_epoch = " << thisJob.epochMJDUTC
					 << "\nreference_Keplerian =";

		for (double value : output_kep_values) ofs_solution << ' ' << value;

		ofs_solution << "\nreference_Cartesian =";

		for (double value : output_cart_values) ofs_solution << ' ' << value;

		ofs_solution << "\nrms_angle_rad = " << thisJob.solution.RMS_angle;

		if (thisJob.solution.flag != 0) ofs_solution << "\nsolution_flag = " << thisJob.solution.flag;

		ofs_solution << '\n';
		ofs_solution.close();
		if (!ofs_solution)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing orbit solution file: " + thisJob.fn_solution, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 配置了残差文件路径时，按观测记录的顺序写出残差。
		if (thisJob.fn_residual.empty())
		{
			// 不输出残差
			tracer.LogEvent(__FILE__, __func__, __LINE__, "No residual output.", EventType::MESSAGE, FORCE_ALERT_OFF);
		}
		else
		{
			ofstream ofs_residual(thisJob.fn_residual.c_str());
			if (!ofs_residual)
			{
				msg = "Cannot open single IOD residual file: " + thisJob.fn_residual;
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			msg = "Writing residuals to: " + thisJob.fn_residual;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';

			// 每行写入历元、观测类型、目标号、站号及东西和南北方向的 O-C 残差（rad）。
			for (int k = 0; k < num; k++)
			{
				const ObservationRecord& record = records[k];
				ofs_residual << std::fixed << std::setprecision(16) << record.mjd_utc << ' '
							 << record.type << ' ' << record.object_id << ' ' << record.observer_id
							 << ' ' << std::scientific << std::setprecision(16)
							 << thisJob.solution.resRa(k) << ' ' << thisJob.solution.resDec(k)
							 << '\n';
			}
			ofs_residual.close();
			if (!ofs_residual)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing residuals.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// =================================== 15 组织最终摘要并保存报告与日志 ===================================

		// 峰值常驻内存在 macOS 中以 byte 返回，在 Linux 中以 KiB 返回，显示时统一为 MiB。
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

		// 结果对象显示科学结果，主程序在同一摘要末尾补充任务耗时和内存。
		std::ostringstream summary;
		summary << "\n====================================== GL summary ======================================\n\n";
		thisJob.solution.PrintSummary(summary);
		summary << '\n' << std::left << std::setw(label_width) << "- Total elapsed time [s]" << ": "
				<< std::fixed << std::setprecision(2)
				<< std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count()
				<< '\n' << std::setw(label_width) << "- Peak memory [MiB]" << ": ";

		if (std::isfinite(peak_memory_mib))
		{
			summary << std::setprecision(1) << peak_memory_mib;
		}
		else
		{
			summary << "Not available";
		}

		summary << '\n';

		// 任务报告保存可阅读的设置、逐轮记录和结果；事件序号与源码位置留在独立程序日志中。
		if (!thisJob.fn_report.empty())
		{
			msg = "Writing GL report to: " + thisJob.fn_report;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
			std::ofstream report(thisJob.fn_report);
			if (!report)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open GL report: " + thisJob.fn_report, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			report << "GL initial orbit determination report\n"
				   << "Program       : " << program_name << " | Version: v" << ASTER_PROGRAM_VERSION
				   << "\nJob ID        : " << jobId << "\nStarted at    : " << timeStr
				   << "\nInput setting : " << fn_iod_input << '\n'
				   << general_settings.str() << observation_summary.str() << iod_settings.str()
				   << iteration_report.str() << summary.str();
			report.close();
			if (!report)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing report: " + thisJob.fn_report, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 将日志路径纳入同组写出消息，成功事件也保存到本次程序日志。
		if (program_log.is_open())
		{
			msg = "Writing program log to: " + fn_programlog;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
		}

		const std::string completed = "GL job completed successfully.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, completed, EventType::MESSAGE, FORCE_ALERT_OFF);

		// 所有文件保存成功后再显示最终摘要和完成消息，实际写入错误仍由 tracer 立即报告。
		if (program_log.is_open())
		{
			tracer.ReportEvents(program_log);
			program_log << "\nExit status: 0\n";
			program_log.close();
			if (!program_log)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing program log.", EventType::ERROR, FORCE_ALERT_ON);
				tracer.EmergencyExit();
			}
		}

		std::cerr << summary.str() << '\n' << file_messages.str();
		std::cerr << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << completed << '\n' << std::setfill(' ');

		tracer.CheckOut(__func__);
		return 0;
	}

	// =================================== 16 处理标准库异常 ===================================
	catch (const std::exception& error)
	{
		// 将文件系统操作、内存分配等标准库异常交给 tracer，显示错误原因并结束任务。
		tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
		tracer.EmergencyExit(program_log);
		return 1;
	}
}
