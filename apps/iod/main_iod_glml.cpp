/**
 * @file main_iod_glml.cpp
 * @brief 根据同一目标的单站或多站 RA/Dec 观测，使用 GLML 方法选择参考历元的轨道。
 *
 * @details
 * 程序先读取配置、观测和测站资料，检查本次选用的观测，再计算各观测时刻的测站位置。
 * 随后求基准 GL/GL-O 解，再逐个指定偏心率，调用GL 核心生成候选轨道。
 * 候选经过轨道边界处理后重新计算残差，以东西、南北等标准差高斯模型的负对数似然 J 选择轨道。
 * 采样、边界处理、评价及输出在 main 中顺序展开；候选轨道不再做额外的无约束迭代。
 *
 * @par 输入资料
 * 配置文件采用 SimpleInput 格式，可通过 -I、-input 或 --input 指定；省略配置参数时显示用法。
 * 观测文件每行依次为 MJD_UTC、观测类型、目标号、站号、RA(deg)、Dec(deg)，本程序使用类型 1。
 * 测站文件每行以站号和坐标类型开头：类型 1 后接大地纬度 phi(deg)、经度 lambda(deg)、高度(m)；
 * 类型 2 后接 ITRS 直角坐标 X(m)、Y(m)、Z(m)。观测通过站号与测站资料关联。
 * EOP 文件提供测站坐标转换所需的地球定向参数。配置中的相对路径以配置文件所在目录为基准。
 *
 * @par 时间、坐标与单位
 * 观测时刻和参考历元均采用 UTC MJD；定轨使用 GCRS 下的地心位置与速度。
 * 角度在计算前转换为 rad，位置、速度和相对时间在核心计算中使用项目定义的归一化单位。
 * 每条观测的测站位置按该条观测时刻计算，EOP 统一取参考历元的值。
 *
 * @par 输出结果
 * 轨道解文件采用 SimpleInput 格式，保存参考历元、轨道状态、测角残差 RMS、候选编号及 J。
 * source_solution_flag 保存边界处理前的来源 IOD 求解状态；它不代表处理后的轨道重新收敛。
 * 轨道状态的输出单位由 output_state_units 指定；残差文件、运行报告和运行日志按配置选择输出。
 * 残差文件每行依次为 MJD_UTC、观测类型、目标号、站号、东西方向残差和南北方向残差；
 * 两个残差均取 O-C，单位为 rad，文件不含表头。
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
#include <omp.h>
#include <random>

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

#include "core_iod.hpp"
#include "time_transform.hpp"

using namespace constants::Scale;
using namespace constants::Conversion;
using constants::PI;

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
 * @brief 读取测角资料，生成 GLML 候选，并按负对数似然选择和输出轨道。
 * @param[in] argc 命令行参数个数，包含可执行程序名。
 * @param[in] argv 命令行参数；-I、-input 或 --input 后指定 SimpleInput 配置文件路径。
 * @return 成功完成候选选择和文件输出时返回 0；错误由 tracer 记录并以非零状态结束。
 * @details 观测与通用设置、候选生成和评价统计在主程序组织；共同的 IOD 设置、
 * 逐轮记录及轨道摘要使用现有 GL 接口。科学显示直接写入控制台与报告，事件另存日志。
 */
int main(int argc, char* argv[])
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_iod_glml " << ASTER_PROGRAM_VERSION << '\n';
        std::cout << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_iod_glml -I <setting.txt>\n";
        return 2;
    }

	// =================================== 01 初始化任务与事件追踪 ===================================

	Overwatch tracer(__FILE__, __func__, __LINE__);
	std::ofstream program_log;

	try
	{
		const double start = omp_get_wtime();
		const std::string program_name = fs::path(argv[0]).filename().string();

		// 用当前本地时间生成任务编号和启动时间，便于在输出信息中识别本次运行。
		std::time_t timeNow = std::time(nullptr);
		std::tm* tmNow = std::localtime(&timeNow);
		char timeStr[20];
		char jobId[15];

		std::strftime(jobId, sizeof(jobId), "%Y%m%d%H%M%S", tmNow);
		std::strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", tmNow);

		std::stringstream jobMsg;
		jobMsg << "GLML initial orbit determination job started at " << timeStr
			   << " (jobId = " << jobId << ")";
		tracer.LogEvent(__FILE__, __func__, __LINE__, jobMsg.str(), bastion::EventType::MESSAGE, FORCE_ALERT_ON);

		const std::string version_message = "Program: " + program_name + " | Version: v" + ASTER_PROGRAM_VERSION
			+ " | " + aster::math_backend::BuildDescription();
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
				cout << "Usage: aster_iod_glml [-I|-input|--input <setting>]\n"
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

		string msg;             // 当前阶段的事件说明。
		IodJob thisJob{};        // 本次定轨的设置、观测资料和计算结果。
		SimpleInput si_iodInput; // 配置记录及按关键词读取的接口。

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

		string epochSetting;       // First、Mid、Last 或明确的 UTC MJD。
		string fn_candidates_GLML; // 可选的候选表路径。
		double noiseLevelAs;       // 东西、南北共同的高斯误差标准差，arcsec。

		// 【必选】file_observation：地面测站 RA/Dec 观测文件路径，不能留空。
		si_iodInput.GetString(tracer, "file_observation", thisJob.fn_obs);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】file_EOP：地球定向参数文件路径，用于将测站位置从 ITRS 转换到 GCRS。
		si_iodInput.GetString(tracer, "file_EOP", thisJob.fn_EOP);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】file_orbit_solution_GLML：参考历元轨道解的输出路径，不能留空。
		si_iodInput.GetString(tracer, "file_orbit_solution_GLML", thisJob.fn_solution);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_report：任务报告路径，保存基准轨道、最终选解及逐轮记录；省略或留空时关闭输出。
		si_iodInput.GetString(tracer, "file_report", thisJob.fn_report);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_residual：逐条观测残差的输出路径；省略或留空时关闭输出。
		// 残差为东西和南北方向的 O-C，单位 rad。
		si_iodInput.GetString(tracer, "file_residual", thisJob.fn_residual);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】file_candidates_GLML：全部候选的 CSV 路径；省略或留空时关闭。
		// 根数固定使用 km、deg，Cartesian 使用 km、km/s；RMS 为 rad，J 为负对数似然。
		si_iodInput.GetString(tracer, "file_candidates_GLML", fn_candidates_GLML);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】iteration_method_GL：GLML 的固定偏心率候选需要微分改正，取 2。
		si_iodInput.GetInteger(tracer, "iteration_method_GL", thisJob.iterMethod);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】constraint_type_GL：基准解的约束；0 为无约束，1 为方向，2 为圆轨道，3 为方向＋圆轨道，4 为方向＋偏心率。
		// 当前读取器在缺项时返回 0，程序按无约束处理。
		si_iodInput.GetInteger(tracer, "constraint_type_GL", thisJob.constraintType);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 固定偏心率约束需要微分改正，检查本次配置是否选择了方法 2。
		if (thisJob.iterMethod != IOD_ITERATION_NEWTON)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "GLML eccentricity sampling requires iteration_method_GL=2.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 基准解可使用已实现的五种约束组合。
		if (thisJob.constraintType < IOD_CONSTRAINT_NONE || thisJob.constraintType > IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown constraint_type_GL.", EventType::ERROR);
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

		// 【必选】iteration_maxround：基准解的最大迭代次数，取正整数；固定偏心率候选各自使用 40 轮上限。
		si_iodInput.GetInteger(tracer, "iteration_maxround", thisJob.iterCountMax);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【必选】iteration_state_tolerance：六个归一化状态改正分量的最大绝对值阈值，取有限正数。
		si_iodInput.GetDouble(tracer, "iteration_state_tolerance", thisJob.iterTOL);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 【可选】display_iterations：0 为关闭屏幕迭代信息，1 为显示；缺项时读取器返回 0。
		// 此项控制基准解及所选来源记录的屏幕显示；启用 file_report 后写出相应迭代记录。
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
		// 读入后换算为 rad，用于评价文件中观测的似然。
		si_iodInput.GetDouble(tracer, "noise_sigma_arcsecond", noiseLevelAs);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		if (!isfinite(noiseLevelAs) || noiseLevelAs <= 0.0)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "noise_sigma_arcsecond must be positive and finite.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 似然公式需要计算方差的倒数和对数，换算后的 rad 方差应能表示为有限正数。
		thisJob.obs.sigma_RaDec = noiseLevelAs * ARCSEC_TO_RAD;

		if (!std::isfinite(2 * PI * thisJob.obs.sigma_RaDec * thisJob.obs.sigma_RaDec) ||
			thisJob.obs.sigma_RaDec * thisJob.obs.sigma_RaDec == 0.0)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "noise_sigma_arcsecond must yield a finite positive radian variance.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 【可选】random_seed：随机回退的基础种子，整数范围 [0,2147483647]，省略时取 0。
		// 候选按编号建立独立序列，使抽样不受 OpenMP 线程数和调度顺序影响。
		int random_seed = 0;

		if (si_iodInput.HasKeyword("random_seed"))
		{
			si_iodInput.GetIntegerStrict(tracer, "random_seed", random_seed);
			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		if (random_seed < 0)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "random_seed must be an integer in [0,2147483647].", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

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
			tracer.LogEvent(__FILE__, __func__, __LINE__, "file_orbit_solution_GLML is required.", EventType::ERROR);
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

		if (!fn_candidates_GLML.empty()) fn_candidates_GLML = (configdir / fn_candidates_GLML).lexically_normal().string();

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

			// 基准轨道和固定偏心率候选均至少需要三条观测。
			if (readCount < 3)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "observations_to_use must be All or an integer >= 3.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// =================================== 06 读取原始观测并关联测站 ===================================

		// 完整读入观测和测站资料，并按站号为每条观测关联相应的测站。
		// 读入的资料保留文件中的单位，供后续检查和坐标转换使用。

		// 读入观测记录，保留 UTC 历元、原始身份列及 deg 制测角。
		msg = "Reading observations from: " + thisJob.fn_obs;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
		Observations iod_observations;
		iod_observations.ReadFile(tracer, thisJob.fn_obs);
		if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

		// 读入各测站的大地坐标或 ITRS 直角坐标，供逐条观测关联。
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

			// 当前入口使用地面测站的 RA/Dec 观测，对应类型 1。
			if (record.type != 1)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "GLML supports ground RA/Dec (type 1) only.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 各测站提供的观测应属于同一个目标。
			if (record.object_id != thisJob.objectId)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, location + "GLML requires one target; ground station IDs may differ.",
					EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 将首次出现的站号加入列表，用于汇总本次参与观测的测站。
			if (std::find(station_ids.begin(), station_ids.end(), record.observer_id) == station_ids.end()) station_ids.push_back(record.observer_id);

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

		// =================================== 08 准备算法使用的时间与角度数组 ===================================

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

			// 将赤经、赤纬从 deg 换算为 rad；赤经归入 [0,2π)。
			double ra_rad = std::fmod(record.ra_deg * DEGREE_TO_RAD, 2 * PI);

			if (ra_rad < 0.0) ra_rad += 2 * PI;

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
				site_geodetic.Normalize(tracer, UNIT_LENGTH_M);
				site_geodetic.ToITRS(tracer, site_itrs);
			}
			else
			{
				// 类型 2 直接提供 ITRS XYZ(m)，以 ae 归一化后用于地固到天球的变换。
				site_itrs.Set(ground_observer.position_itrs_m[0], ground_observer.position_itrs_m[1],
					ground_observer.position_itrs_m[2], false);
				site_itrs.Normalize(tracer, UNIT_LENGTH_M);
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

		// =================================== 11 求基准轨道并准备候选初值 ===================================

		// 在 [0,0.8) 均匀取 10000 个偏心率；候选 0 对应基准解，1..10000 对应指定偏心率解。
		const int eccentricity_samples = 10000;
		const int candidate_count = eccentricity_samples + 1;
		const double ecc_min = 0.0, ecc_max = 0.8;

		// 边界处理允许的偏心率为 [0,0.9]，近地点高度不低于 170 km，远地点高度不高于 52884.148 km。
		// 高度先换算为以 ae 归一化的长度，随后用于确定半长轴的可行区间。
		const double e_min = 0.0, e_max = 0.9;
		const double hp_min = 170.0 * KM_TO_NORMALIZED_LENGTH;
		const double ha_max = 52884.148 * KM_TO_NORMALIZED_LENGTH;

		// 显示准备好的任务资料；同一份文字随后写入任务报告。
		thisJob.initValue.type = IOD_INITIAL_VALUE_DEFAULT;
		const int label_width = 40;
		msg = "Input checks and model preparation completed.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		std::ostringstream general_settings;
		general_settings << "\n=================================== General settings ===================================\n\n"
						 << std::left << std::setw(label_width) << "- Reference frame" << ": GCRS\n"
						 << std::setw(label_width) << "- Input / Output Time" << ": UTC\n"
						 << std::setw(label_width) << "- Target ID" << ": " << thisJob.objectId << '\n'
						 << std::setw(label_width) << "- Reference epoch [UTC]" << ": "
						 << std::fixed << std::setprecision(12) << thisJob.epochMJDUTC << '\n';
		std::cerr << general_settings.str();

		// 观测概览只统计参与求解的记录，按时间范围、观测者、记录数的顺序显示。
		std::size_t unique_epochs = 0;

		for (std::size_t k = 0; k < used_count; ++k)
		{
			if (k == 0 || records[k].mjd_utc != records[k - 1].mjd_utc) ++unique_epochs;
		}

		std::ostringstream observation_summary;
		observation_summary << "\n================================= Observation summary ==================================\n\n"
							<< std::left << std::setw(label_width) << "- First / Last observation [UTC]" << ": "
							<< std::fixed << std::setprecision(12) << records[0].mjd_utc << " / " << records[used_count - 1].mjd_utc
							<< '\n' << std::setw(label_width) << "- Arc span [s]" << ": " << std::setprecision(2)
							<< (records[used_count - 1].mjd_utc - records[0].mjd_utc) * constants::SECONDS_PER_DAY
							<< '\n' << std::setw(label_width) << "- Unique epochs" << ": " << unique_epochs
							<< "\n\n" << std::setw(label_width) << "- Observers" << ": " << station_ids.size() << " (";

		for (std::size_t k = 0; k < station_ids.size(); ++k)
		{
			if (k > 0) observation_summary << ", ";

			observation_summary << station_ids[k];
		}

		observation_summary << ")\n\n";

		if (used_count != records.size())
		{
			observation_summary << std::setw(label_width) << "- Loaded observation records" << ": " << records.size() << '\n';
		}

		observation_summary << std::setw(label_width) << "- Observation records" << ": " << used_count << '\n'
							<< std::setw(label_width) << "  > Ground-based RA/Dec (1)" << ": " << used_count << '\n'
							<< std::setw(label_width) << "- Scalar components" << ": " << 2 * used_count << '\n';
		std::cerr << observation_summary.str();

		// 基准 GL 的设置、迭代和摘要连续显示，完成后再进入 GLML 候选阶段。
		std::cerr << '\n';
		msg = "Starting baseline GL orbit determination.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
		std::ostringstream iod_settings;
		thisJob.Display(iod_settings);
		std::cerr << iod_settings.str();

		// 基准 GL/GL-O 解使用配置给出的约束、迭代上限和容差。
		// 求解事件单独保存，数值失败时仍可按下面的随机距离策略准备候选初值。
		std::cerr << '\n';
		Overwatch baseline_tracer(__FILE__, __func__, __LINE__);
		iod_solver::InitialOrbitDetermination(baseline_tracer, thisJob);
		if (baseline_tracer.Code() != ERR_NOERROR) thisJob.solution.flag = -2;

		// 立即汇总基准结果；数值失败时不把尚未完成的状态当作轨道解显示。
		std::ostringstream baseline_summary;
		baseline_summary << "\n================================== Baseline GL summary ==================================\n\n";

		if (thisJob.solution.flag == -2)
		{
			baseline_summary << std::left << std::setw(label_width) << "- Iteration Status" << ": NUMERICAL FAILURE\n";
		}
		else
		{
			thisJob.solution.PrintSummary(baseline_summary);
		}

		std::cerr << baseline_summary.str() << '\n';
		msg = "Starting GLML candidate generation and selection.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		std::ostringstream sampling;
		sampling << "\n==================================== GLML settings =====================================\n\n"
				 << std::left << std::setw(label_width) << "- Random seed" << ": " << random_seed
				 << '\n' << std::setw(label_width) << "- Eccentricity samples" << ": " << eccentricity_samples
				 << '\n' << std::setw(label_width) << "- Sampling eccentricity" << ": [" << ecc_min << ", " << ecc_max << ")"
				 << '\n' << std::setw(label_width) << "- Total candidates" << ": " << candidate_count << " (including baseline)"
				 << '\n' << std::setw(label_width) << "- Candidate constraints" << ": Orientation + Eccentricity"
				 << '\n' << std::setw(label_width) << "- Candidate initial state" << ": Specified"
				 << '\n' << std::setw(label_width) << "- Candidate max iterations" << ": 40"
				 << "\n\n" << std::setw(label_width) << "- Eccentricity bounds" << ": [" << e_min << ", " << e_max << "]"
				 << '\n' << std::setw(label_width) << "- Perigee / Apogee height [km]" << ": >= "
				 << std::fixed << std::setprecision(3) << hp_min * NORMALIZED_LENGTH_TO_KM
				 << " / <= " << ha_max * NORMALIZED_LENGTH_TO_KM
				 << "\n\n" << std::setw(label_width) << "- Selection criterion" << ": Minimum J (negative log likelihood)\n";
		std::cerr << sampling.str() << '\n';

		// 候选共用参考历元的六维 GCRS 初值，位置以 ae、速度以 ae/UNITTIME 归一化。
		Cartesian candidate_initial_cart;

		if (thisJob.solution.flag == -2)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__,
				"Baseline IOD failed numerically; initialize candidates with a random distance along the reference observation.",
				EventType::WARNING, FORCE_ALERT_ON);

			// 随机距离初始化需要参考历元处的观测方向，First/Mid/Last 均直接取相应记录。
			int epoch_index = 0;

			if (epochSetting == "Mid")
			{
				epoch_index = num / 2;
			}
			else if (epochSetting == "Last")
			{
				epoch_index = num - 1;
			}
			else if (epochSetting != "First")
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__,
					"Failed-baseline random-distance initialization requires First, Mid or Last.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 用参考观测的赤经、赤纬构造 GCRS 单位视线，并取得同历元的测站位置。
			const double alpha = thisJob.obs.ra_rad(epoch_index), delta = thisJob.obs.dec_rad(epoch_index);
			const Vector3d los(std::cos(delta) * std::cos(alpha), std::cos(delta) * std::sin(alpha), std::sin(delta));
			const Vector3d station = thisJob.observer.positions_gcrs_norm_emat.col(epoch_index);

			// 在归一化距离区间 [hp_min,ha_max) 抽样，以 r=R+rho*los 得到位置，初速度取零。
			std::seed_seq seed{static_cast<unsigned int>(random_seed), 0u};
			std::mt19937 generator(seed);
			std::uniform_real_distribution<double> range_distribution(hp_min, ha_max);
			const Vector3d r_guess = station + range_distribution(generator) * los;
			candidate_initial_cart.Set(r_guess(0), r_guess(1), r_guess(2), 0.0, 0.0, 0.0, true);
		}
		else
		{
			candidate_initial_cart = thisJob.solution.cart;
		}

		// 一条候选把边界处理后的轨道、残差和评价值放在一起，供选择及输出共用。
		// 来源 IOD 的状态和迭代事件另外保存，避免与边界处理后的结果混淆。
		struct OrbitCandidate
		{
			Cartesian cart;        // 参考历元的 GCRS 六维状态，归一化。
			Keplerian kep;         // 对应的 Keplerian 根数，长度归一化，角度为 rad。
			ArrayXd resRa, resDec;  // 每条观测的东西、南北 O-C 残差，rad。
			double rms = std::numeric_limits<double>::quiet_NaN(); // 全部 2N 个分量的合成 RMS，rad。
			double j = std::numeric_limits<double>::quiet_NaN();   // 等标准差高斯模型的负对数似然。
		};

		// 处理后的候选、来源求解任务与来源事件按相同编号对应。
		vector<OrbitCandidate> candidates(candidate_count);
		vector<IodJob> candidate_jobs(candidate_count);
		vector<Overwatch> candidate_tracers(candidate_count);
		candidate_jobs[0] = thisJob;

		// =================================== 12 逐个指定偏心率并行生成候选轨道 ===================================

		std::cerr << '\n';
		msg = "Generating " + std::to_string(eccentricity_samples) + " eccentricity-constrained candidates.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);

		// 记录本次实际线程数；每个线程仅修改当前编号的候选任务及其事件记录。
		int actual_threads = 1;
		int completed_candidates = 0;
		const int progress_step = std::max(1, eccentricity_samples / 10);
		const auto generation_start = std::chrono::steady_clock::now();
		double last_progress_time = 0.0;

		// 进度只直接显示到控制台，耗时从本次候选生成开始计算。
		std::fprintf(stderr, "\nProgress: %6.2f%% | candidates %d / %d | elapsed %.2f s\n", 0.0, 0, eccentricity_samples, 0.0);
		std::fflush(stderr);

#pragma omp parallel
		{
#pragma omp master
			{
				actual_threads = omp_get_num_threads();
			}

#pragma omp for schedule(dynamic)
			for (int k = 1; k < candidate_count; ++k)
			{
				candidate_tracers[k].Set(__FILE__, __func__, __LINE__, k);
				IodJob& candidate_job = candidate_jobs[k];
				candidate_job = thisJob;
				candidate_job.printIterInfo = false;

				// 从共用初值出发，指定当前候选的偏心率；每条观测使用三个基本行、一个偏心率行和一个方向行。
				candidate_job.initValue.type = IOD_INITIAL_VALUE_SPECIFIED;
				candidate_job.initValue.initCart = candidate_initial_cart;

				candidate_job.constraintEcc = ecc_min + (ecc_max - ecc_min) * (k - 1) / eccentricity_samples;
				candidate_job.constraintType = IOD_CONSTRAINT_ORIENTATION_ECCENTRICITY;
				candidate_job.iterCountMax = 40;

				// 求解后保存来源状态：正数为收敛轮次，-1 为达到迭代上限，-2 为数值失败。
				// 本层保留求解事件和任务中留存的轨道，后续统一进行边界处理与残差评价。
				iod_solver::InitialOrbitDetermination(candidate_tracers[k], candidate_job);
				if (candidate_tracers[k].Code() != ERR_NOERROR) candidate_job.solution.flag = -2;

				// 按实际完成数更新进度；每完成 10% 或距上次显示达到一秒时输出。
				// 计数和整行输出共用临界区，使并行完成的候选按递增顺序显示。
#pragma omp critical(glml_candidate_progress)
				{
					++completed_candidates;
					const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - generation_start).count();

					if (completed_candidates % progress_step == 0 || elapsed - last_progress_time >= 1.0 ||
						completed_candidates == eccentricity_samples)
					{
						std::fprintf(stderr, "Progress: %6.2f%% | candidates %d / %d | elapsed %.2f s\n",
							100.0 * completed_candidates / eccentricity_samples, completed_candidates, eccentricity_samples, elapsed);
						std::fflush(stderr);
						last_progress_time = elapsed;
					}
				}
			}
		}

		std::cerr << '\n';

		// 按来源 IOD 的 flag 汇总指定偏心率候选，基准解不计入这三个计数。
		int converged_count = 0, limited_count = 0, failed_count = 0;

		for (int k = 1; k < candidate_count; ++k)
		{
			if (candidate_jobs[k].solution.flag > 0)
			{
				++converged_count;
			}
			else if (candidate_jobs[k].solution.flag == -1)
			{
				++limited_count;
			}
			else
			{
				++failed_count;
			}
		}

		// =================================== 13 处理轨道边界并计算东西、南北残差与 J ===================================

		// 用边界处理后的状态反算各次观测，残差及 RMS 由 Compute_Residuals 保存。
		// sigma 为两个方向共同的标准差，单位 rad；resRa 已乘观测赤纬的 cos，表示东西分量。
		msg = "Evaluating " + std::to_string(candidate_count) + " candidates after orbit bounds handling.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_ON);
		const double sigma = thisJob.obs.sigma_RaDec;
		vector<int> invalid_bounds(candidate_count, 0);
		vector<Overwatch> evaluation_tracers(candidate_count);

#pragma omp parallel for schedule(dynamic)
		for (int k = 0; k < candidate_count; ++k)
		{
			// 阶段编号 1 与基准失败时的距离抽样分开；候选编号确定本候选的序列。
			std::seed_seq seed{static_cast<unsigned int>(random_seed), 1u, static_cast<unsigned int>(k)};
			std::mt19937 generator(seed);

			double kep_norm[6];
			candidate_jobs[k].solution.kep.ToArray(kep_norm);

			// 偏心率先补全 NaN，再限制到 [0,0.9]；边界长度采用归一化单位。
			if (std::isnan(kep_norm[1]))
			{
				std::uniform_real_distribution<double> distribution(e_min, e_max);
				kep_norm[1] = distribution(generator);
			}

			kep_norm[1] = std::clamp(kep_norm[1], e_min, e_max);

			// 由 a(1-e)>=1+hp_min、a(1+e)<=1+ha_max 确定半长轴范围；1 为归一化地球半径。
			const double a_min = (hp_min + 1.0) / (1.0 - kep_norm[1]);
			const double a_max = (ha_max + 1.0) / (1.0 + kep_norm[1]);

			// 两个高度限制没有交集时，该候选无法满足边界；不对反序区间调用 clamp。
			if (a_min > a_max)
			{
				invalid_bounds[k] = 1;
				continue;
			}

			// 半长轴为 NaN 时先在可行区间内补全，再将其限制到同一可行区间。
			if (std::isnan(kep_norm[0]))
			{
				std::uniform_real_distribution<double> distribution(a_min, a_max);
				kep_norm[0] = distribution(generator);
			}

			kep_norm[0] = std::clamp(kep_norm[0], a_min, a_max);

			// 倾角为 NaN 时在 [0,pi) rad 内均匀抽样补全。
			if (std::isnan(kep_norm[2]))
			{
				std::uniform_real_distribution<double> distribution(0.0, PI);
				kep_norm[2] = distribution(generator);
			}

			// 升交点赤经为 NaN 时在 [0,2pi) rad 内均匀抽样补全。
			if (std::isnan(kep_norm[3]))
			{
				std::uniform_real_distribution<double> distribution(0.0, 2 * PI);
				kep_norm[3] = distribution(generator);
			}

			// 近地点幅角或平近点角存在 NaN 时，分别在 [0,2pi) rad 内均匀抽样补全两者。
			if (std::isnan(kep_norm[4]) || std::isnan(kep_norm[5]))
			{
				std::uniform_real_distribution<double> distribution(0.0, 2 * PI);
				kep_norm[4] = distribution(generator);
				kep_norm[5] = distribution(generator);
			}

			// 评价阶段单独记录事件，来源 IOD 的诊断保持可查。
			Overwatch& evaluation_tracer = evaluation_tracers[k];
			evaluation_tracer.Set(__FILE__, __func__, __LINE__, k);

			// 将处理后的根数转换为同一参考历元的 GCRS Cartesian 状态，位置和速度均归一化。
			Keplerian bounded_kep(kep_norm[0], kep_norm[1], kep_norm[2], kep_norm[3], kep_norm[4], kep_norm[5], true);
			bounded_kep.ToCartesian(evaluation_tracer, candidates[k].cart);
			if (evaluation_tracer.Code() != ERR_NOERROR) continue;

			// 将 [x,y,z,vx,vy,vz] 数组转为六维 Eigen 向量，复算全部观测的两个残差分量及 RMS。
			double cart_norm[6];
			candidates[k].cart.ToArray(cart_norm);
			const VectorXd state = Eigen::Map<VectorXd>(cart_norm, 6);
			double rms_east, rms_north;

			iod_solver::Compute_Residuals(evaluation_tracer, candidate_jobs[k], state, rms_east, rms_north);
			if (evaluation_tracer.Code() != ERR_NOERROR) continue;

			// 从处理后的 Cartesian 状态取得对应根数，保证两种输出表示属于同一条轨道。
			candidates[k].cart.ToKeplerian(evaluation_tracer, candidates[k].kep);
			if (evaluation_tracer.Code() != ERR_NOERROR) continue;

			candidates[k].resRa = candidate_jobs[k].solution.resRa;
			candidates[k].resDec = candidate_jobs[k].solution.resDec;
			candidates[k].rms = candidate_jobs[k].solution.RMS_angle;

			// N 条观测共含 2N 个角度分量，J=N*ln(2*pi*sigma^2)+sum(r_E^2+r_N^2)/(2*sigma^2)。
			// 两个分量直接共用 sigma，resRa 不再乘 cos，标准差也不另乘 cos。
			candidates[k].j = num * std::log(2 * PI * sigma * sigma) +
				(candidates[k].resRa.square().sum() + candidates[k].resDec.square().sum()) / (2 * sigma * sigma);
		}

		// =================================== 14 按最小 J 选择候选并检查最终结果 ===================================

		// 顺序扫描候选编号，相等时保留先出现的候选；非有限 J 不参与选择。
		int best_index = -1;
		double j_best = std::numeric_limits<double>::infinity();

		for (int k = 0; k < candidate_count; ++k)
		{
			if (std::isfinite(candidates[k].j) && candidates[k].j < j_best)
			{
				best_index = k;
				j_best = candidates[k].j;
			}
		}

		if (best_index < 0)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "No finite GLML candidate was selected.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 输出共用已经保存的最优轨道与残差；先检查其状态、根数及质量指标是否有效。
		const OrbitCandidate& selected = candidates[best_index];
		double selected_cart[6], selected_kep[6];
		selected.cart.ToArray(selected_cart);
		selected.kep.ToArray(selected_kep);

		for (int j = 0; j < 6; ++j)
		{
			if (!std::isfinite(selected_cart[j]) || !std::isfinite(selected_kep[j]))
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite selected GLML orbit.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		if (!std::isfinite(selected.rms) || !selected.resRa.isFinite().all() || !selected.resDec.isFinite().all())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite selected GLML residuals or RMS.", EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 生成阶段的三种状态只统计 10000 个指定偏心率来源；评价计数包含基准候选 0。
		const int invalid_bound_count = std::count(invalid_bounds.begin(), invalid_bounds.end(), 1);
		int finite_candidate_count = 0;

		for (const OrbitCandidate& candidate : candidates)
		{
			if (std::isfinite(candidate.j)) ++finite_candidate_count;
		}

		std::ostringstream candidate_status;
		candidate_status << "\n================================== Candidate summary ===================================\n\n"
						 << std::left << std::setw(label_width) << "- Constrained IOD candidates" << ": " << eccentricity_samples
						 << '\n' << std::setw(label_width) << "  > Converged" << ": " << converged_count
						 << '\n' << std::setw(label_width) << "  > Iteration limit" << ": " << limited_count
						 << '\n' << std::setw(label_width) << "  > Numerical failure" << ": " << failed_count
						 << "\n\n" << std::setw(label_width) << "- Evaluated candidates" << ": " << candidate_count
						 << '\n' << std::setw(label_width) << "  > Finite J" << ": " << finite_candidate_count
						 << '\n' << std::setw(label_width) << "  > Infeasible orbit bounds" << ": " << invalid_bound_count
						 << '\n' << std::setw(label_width) << "  > Other evaluation failure" << ": "
						 << candidate_count - finite_candidate_count - invalid_bound_count << '\n';

		std::cerr << candidate_status.str();

		// 此处重放的是所选候选留存的来源迭代，属于边界处理前的轨道。
		// 标出候选编号，使它与前面的基准 GL 迭代清楚区分。
		std::ostringstream selected_source_report;

		if (best_index > 0)
		{
			selected_source_report << "\n======================= Selected source IOD (before orbit bounds) =======================\n\n"
								   << std::left << std::setw(label_width) << "- Candidate ID" << ": " << best_index << '\n'
								   << candidate_jobs[best_index].solution.iteration_report;
		}

		if (thisJob.printIterInfo) std::cerr << selected_source_report.str();

		// =================================== 15 写出候选表、最终轨道和逐条残差 ===================================

		std::ostringstream file_messages;

		// 候选 CSV 写出编号、六个根数、RMS、J、六维 Cartesian 状态及来源 IOD 的 flag。
		if (!fn_candidates_GLML.empty())
		{
			msg = "Writing GLML candidates to: " + fn_candidates_GLML;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
			ofstream output(fn_candidates_GLML);

			if (!output)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open candidate file: " + fn_candidates_GLML, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			output << "id,a(km),e,i(deg),RAAN(deg),omega(deg),M(deg),RMS,J,x(km),y(km),z(km),vx(km/s),vy(km/s),vz(km/s),source_flag\n";
			output << std::defaultfloat << std::setprecision(17);

			for (int k = 0; k < candidate_count; ++k)
			{
				double kep[6], cart[6];
				candidates[k].kep.ToArray(kep);
				candidates[k].cart.ToArray(cart);

				// 未能完成边界处理或残差评价的候选用 NaN 标记，避免默认零状态被误读为轨道。
				if (!std::isfinite(candidates[k].j))
				{
					std::fill(kep, kep + 6, std::numeric_limits<double>::quiet_NaN());
					std::fill(cart, cart + 6, std::numeric_limits<double>::quiet_NaN());
				}

				// 根数按 a,e,i,OMEGA,omega,M 写出，长度恢复为 km，角度恢复为 deg。
				output << k << ',' << kep[0] * NORMALIZED_LENGTH_TO_KM << ',' << kep[1];

				for (int j = 2; j < 6; ++j) output << ',' << kep[j] * RAD_TO_DEGREE;

				output << ',' << candidates[k].rms << ',' << candidates[k].j;

				// Cartesian 按位置、速度分组恢复为 km、km/s；来源 flag 属于边界处理前的求解。
				for (int j = 0; j < 3; ++j) output << ',' << cart[j] * NORMALIZED_LENGTH_TO_KM;

				for (int j = 3; j < 6; ++j) output << ',' << cart[j] * NORMALIZED_VELOCITY_TO_KMPS;

				output << ',' << candidate_jobs[k].solution.flag << '\n';
			}

			output.close();

			if (!output)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing candidate file: " + fn_candidates_GLML, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// 仅对输出副本恢复单位；候选和来源轨道始终保持归一化，供报告和诊断使用。
		Cartesian output_cart = selected.cart;
		Keplerian output_kep = selected.kep;

		// 单位码 1 保留归一化状态；2、3 分别以 km、m 恢复长度，根数角度均恢复为 deg。
		if (thisJob.solutionUnit != IOD_SOLUTION_NORMALIZE)
		{
			const double length_scale = thisJob.solutionUnit == IOD_SOLUTION_KMDEG ? UNIT_LENGTH_M / 1000.0 : UNIT_LENGTH_M;
			output_cart.Restore(tracer, length_scale, UNIT_TIME_S);
			output_kep.Restore(tracer, length_scale);

			if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
		}

		// 最终选解写入轨道文件，原始基准解仅用于后面的对照报告。
		double output_cart_values[6], output_kep_values[6];
		output_cart.ToArray(output_cart_values);
		output_kep.ToArray(output_kep_values);

		msg = "Writing orbit solution to: " + thisJob.fn_solution;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
		file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
		ofstream solution(thisJob.fn_solution);

		if (!solution)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open orbit solution file: " + thisJob.fn_solution, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 轨道文件保存参考历元、两种状态、RMS、J、候选编号及来源状态；数值使用 17 位有效数字。
		solution << std::defaultfloat << std::setprecision(17)
				 << "object_id = " << thisJob.objectId
				 << "\nstate_units = " << thisJob.solutionUnit
				 << "\nreference_epoch = " << thisJob.epochMJDUTC
				 << "\nreference_Keplerian =";

		for (double value : output_kep_values) solution << ' ' << value;

		solution << "\nreference_Cartesian =";

		for (double value : output_cart_values) solution << ' ' << value;

		solution << "\nrms_angle_rad = " << selected.rms
				 << "\nnegative_log_likelihood = " << j_best
				 << "\nselected_candidate = " << best_index
				 << "\nrandom_seed = " << random_seed
				 << "\nsource_solution_flag = " << candidate_jobs[best_index].solution.flag << '\n';
		solution.close();

		if (!solution)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing orbit solution file: " + thisJob.fn_solution, EventType::ERROR);
			tracer.EmergencyExit(program_log);
		}

		// 按原观测的顺序输出已经算好的残差，文件中的身份列取自原始记录。
		if (!thisJob.fn_residual.empty())
		{
			msg = "Writing residuals to: " + thisJob.fn_residual;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
			ofstream residuals(thisJob.fn_residual);

			if (!residuals)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open residual file: " + thisJob.fn_residual, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			// 六列残差带上原始身份；两个分量分别是东西、南北方向的 O-C，单位 rad。
			for (int k = 0; k < num; ++k)
			{
				const ObservationRecord& record = records[k];
				residuals << std::fixed << std::setprecision(16) << record.mjd_utc << ' '
						  << record.type << ' ' << record.object_id << ' ' << record.observer_id
						  << ' ' << std::scientific << std::setprecision(16)
						  << selected.resRa(k) << ' ' << selected.resDec(k) << '\n';
			}

			residuals.close();

			if (!residuals)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing residual file: " + thisJob.fn_residual, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		// =================================== 16 汇总结果并保存报告与运行日志 ===================================

		// 显示来源状态时明确其属于边界处理前的 IOD，不把它作为最终轨道的收敛结论。
		const IodSolution& selected_source = candidate_jobs[best_index].solution;
		const char* source_status = "NUMERICAL FAILURE";

		if (selected_source.flag > 0)
		{
			source_status = "CONVERGED";
		}
		else if (selected_source.flag == -1)
		{
			source_status = "NOT CONVERGED — ITERATION LIMIT";
		}

		// 最终摘要只使用保存的选解和残差；共同的轨道与 RMS 排版交给 IodSolution。
		IodSolution selected_solution{};
		selected_solution.cart = selected.cart;
		selected_solution.kep = selected.kep;
		selected_solution.resRa = selected.resRa;
		selected_solution.resDec = selected.resDec;
		selected_solution.RMS_angle = selected.rms;

		// 峰值常驻内存在 macOS 中以 byte 返回，在 Linux 中以 KiB 返回，统一换算为 MiB。
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

		std::ostringstream summary;
		summary << "\n===================================== GLML summary =====================================\n\n"
				<< std::left << std::setw(label_width) << "- Selected candidate" << ": " << best_index
				<< (best_index == 0 ? " (baseline)" : " (eccentricity-constrained)")
				<< '\n' << std::setw(label_width) << "- Source IOD Iteration Status" << ": " << source_status << " (before orbit bounds)"
				<< '\n' << std::setw(label_width) << "- Source IOD iterations" << ": ";

		if (!selected_source.iterations.empty())
		{
			summary << selected_source.iterations.back().iteration;
		}
		else
		{
			summary << "N/A";
		}

		summary << "\n\n" << std::setw(label_width) << "- J baseline / selected" << ": "
				<< std::scientific << std::setprecision(6) << thisJob.solution.J_likelihood << " / " << j_best
				<< '\n' << std::setw(label_width) << "- RMS baseline / selected [arcsec]" << ": "
				<< std::fixed << std::setprecision(2) << thisJob.solution.RMS_angle / ARCSEC_TO_RAD
				<< " / " << selected.rms / ARCSEC_TO_RAD << "\n\n";
		selected_solution.PrintSummary(summary, false);
		summary << '\n' << std::left << std::setw(label_width) << "- Total elapsed time [s]" << ": "
				<< std::fixed << std::setprecision(2) << omp_get_wtime() - start
				<< '\n' << std::setw(label_width) << "- Peak memory [MiB]" << ": ";

		if (std::isfinite(peak_memory_mib))
		{
			summary << std::setprecision(1) << peak_memory_mib;
		}
		else
		{
			summary << "Not available";
		}

		summary << '\n' << std::setw(label_width) << "- OpenMP threads" << ": " << actual_threads << '\n';

		// 报告按设置、基准求解、候选统计、所选来源、最终轨道的顺序保存科学输出。
		if (!thisJob.fn_report.empty())
		{
			msg = "Writing GLML report to: " + thisJob.fn_report;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
			ofstream report(thisJob.fn_report);

			if (!report)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open report: " + thisJob.fn_report, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}

			report << "GLML initial orbit determination report\n"
				   << "Program       : " << program_name << " | Version: v" << ASTER_PROGRAM_VERSION
				   << "\nBuild         : " << aster::math_backend::BuildDescription()
				   << "\nJob ID        : " << jobId << "\nStarted at    : " << timeStr
				   << "\nInput setting : " << fn_iod_input << '\n'
				   << general_settings.str() << observation_summary.str() << iod_settings.str()
				   << "\n================================= Baseline GL iterations ================================\n"
				   << thisJob.solution.iteration_report << baseline_summary.str()
				   << sampling.str() << candidate_status.str() << selected_source_report.str() << summary.str();
			report.close();

			if (!report)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing report: " + thisJob.fn_report, EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		if (program_log.is_open())
		{
			msg = "Writing program log to: " + fn_programlog;
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::MESSAGE, FORCE_ALERT_OFF);
			file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << msg << '\n';
		}

		const std::string completed = "GLML job completed successfully.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, completed, EventType::MESSAGE, FORCE_ALERT_OFF);

		// 启用运行日志时，保存任务事件及需要保留的各候选诊断。
		if (program_log.is_open())
		{
			tracer.ReportEvents(program_log);
			baseline_tracer.ReportEvents(program_log);

			// 来源求解的错误与警告保留在各自的事件记录中，日志同时保存所选候选的完整过程。
			for (int k = 1; k < candidate_count; ++k)
			{
				if (candidate_tracers[k].Code() == ERR_NOERROR && candidate_tracers[k].WarningCount() == 0 &&
					evaluation_tracers[k].Code() == ERR_NOERROR && k != best_index)
					continue;

				program_log << "\n[Candidate " << k << "]\n";
				candidate_tracers[k].ReportEvents(program_log);
				evaluation_tracers[k].ReportEvents(program_log);
			}

			program_log << "\nExit status: 0\n";
			program_log.close();

			if (!program_log)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing program log.", EventType::ERROR);
				tracer.EmergencyExit(program_log);
			}
		}

		std::cerr << summary.str() << '\n' << file_messages.str();
		std::cerr << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << completed << '\n' << std::setfill(' ');

		tracer.CheckOut(__func__);
		return 0;
	}

	// 文件系统、分配和其他标准库异常统一交给 tracer，保留诊断并退出本次任务。
	catch (const std::exception& error)
	{
		tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
		tracer.EmergencyExit(program_log);
		return 1;
	}
}
