/**
 * @file main_odn_batch.cpp
 * @brief 读取轨道与观测资料，执行数值法批处理定轨，并写出后验结果。
 * @details 在主程序中读取设置、观测及测站资料，准备初轨和动力学模型，
 * 将初轨传到首观测历元后调用批处理核心，再将后验轨道及协方差传回先验参考历元。
 * 轨道文件的字段读取、单位转换、输出格式和精度均在相应使用位置直接写明。
 *
 * @par 输入资料
 * 通过 -I、-input 或 --input 指定 SimpleInput 设置文件；相对路径以该文件所在目录为基准。
 * 观测类型 1 为地基 RA/Dec，6 为自带平台 GCRS 位置的天基 RA/Dec，7 为目标 GCRS XYZ。
 * 测角输入为 deg，位置输入为 m；测站表保存原始大地坐标或 ITRS XYZ。
 * 先验轨道可在设置中按 km、km/s、deg 填写，或通过 file_apriori_orbit 指定轨道文件。
 * 轨道文件的 state_units=1/2/3 分别表示归一化/rad、km与deg、m与deg。
 *
 * @par 时间、坐标与单位
 * 观测、先验及后验参考历元采用 UTC MJD，地球轨道传播使用 GCRS/TT。
 * 核心计算采用归一化 NSG2，长度尺度为 ae，时间尺度为 UNITTIME，角度为 rad。
 * 测角残差采用东西方向 ΔRA*cos(观测Dec) 与南北方向 ΔDec，单位 rad；位置残差为 m。
 *
 * @par 输出结果
 * 报告、日志、逐条残差、三种后验星历及 SimpleInput 后验轨道解分别由设置指定路径。
 * 空路径关闭对应输出；输出父目录由使用者准备，同名已有文件会被拒绝。
 */
#include "math_backend.hpp"
#include "model_observation.hpp"
#include "io_settings.hpp"
#include "frame_definition.hpp"
#include "frame_earth.hpp"
#include "core_bodn.hpp"
#include "io_common.hpp"
#include <algorithm>
#include <cmath>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>
#include <map>
#include <set>
#include <sstream>

#if defined(__APPLE__) || defined(__linux__)
#include <sys/resource.h>
#endif

using namespace determination;
using namespace bastion;

int main(int argc, char **argv)
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_odn_batch " << ASTER_PROGRAM_VERSION << '\n';
        std::cout << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_odn_batch -I <setting.txt>\n";
        return 2;
    }


    // =================================== 01 初始化 ===================================

    // 先创建 tracer，记录 main 的调用过程，以及后续读取、计算中发生的事件和错误。
    Overwatch tracer(__FILE__, __func__, __LINE__);

    // 日志流在整个 main 中有效，任何阶段失败都先保存事件和调用链。
    std::ofstream program_log;

    // 预期错误由 tracer.Code() 就地控制；末尾 catch 只接收外部接口或标准库异常。
    try
    {
        namespace fs = std::filesystem;
        std::error_code filesystem_error;

        // 用 steady_clock 统计运行耗时，再用当前日期和时间生成日志时间及任务编号。
        const auto wall_start = std::chrono::steady_clock::now();
        const std::time_t now = std::time(nullptr);
        char started[20], job_id[15];

        // 日志中的时间写成 YYYY-MM-DD HH:MM:SS，便于直接查看任务何时开始。
        std::strftime(started, sizeof(started), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

        // 再将同一时刻写成 YYYYMMDDHHMMSS，作为本次运行的任务编号。
        std::strftime(job_id, sizeof(job_id), "%Y%m%d%H%M%S", std::localtime(&now));

        // 记录任务的开始时间和编号，并通过 FORCE_ALERT_ON 立即显示在控制台。
        tracer.LogEvent(__FILE__, __func__, __LINE__,
                        std::string("BODN job started at ") + started + " (jobId = " + job_id + ")",
                        EventType::MESSAGE, FORCE_ALERT_ON);

        // 程序名使用实际执行文件名，版本由 CMake 提供，并与预报程序保持一致。
        tracer.LogEvent(__FILE__, __func__, __LINE__,
                        "Program: " + fs::path(argv[0]).filename().string() + " | Version: v" + ASTER_PROGRAM_VERSION,
                        EventType::MESSAGE, FORCE_ALERT_ON);

        // 默认使用工程算例；普通命令行循环只替换 setting 文件名。
        std::string filename = "";
        for (int i = 1; i < argc; ++i)
        {
            std::string argument = argv[i];
            // 帮助不读配置；未知开关直接报告，避免误跑默认算例。
            if (argument == "-h" || argument == "--help")
            {
                std::cout << "aster_odn_batch [-I|-input|--input <setting>]\n"
                             "Each output file has its own path in setting. Prepare parent directories before running.\n";
                tracer.CheckOut(__func__);
                return 0;
            }
            if (argument != "-I" && argument != "-input" && argument != "--input")
                throw std::runtime_error("Unknown argument: " + argument + ". Use aster_odn_batch -I <setting>.");
            // 配置开关后必须跟随文件名，文件读取错误仍交给 SimpleInput 处理。
            if (i + 1 >= argc || argv[i + 1][0] == '\0' || argv[i + 1][0] == '-')
                throw std::runtime_error("Missing value for " + argument);
            filename = argv[++i];
        }

        // 显示本次实际读取的 setting 路径，运行时可以据此确认是否选对了算例。
        tracer.LogEvent(__FILE__, __func__, __LINE__,
                        "Reading and analyzing input from: " + fs::absolute(filename).lexically_normal().string(),
                        EventType::MESSAGE, FORCE_ALERT_ON);

        // 先载入 setting 文件，随后逐项读取配置。
        io::SimpleInput input;
        input.ReadSimpleInput(tracer, filename);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
        // 设置读取成功后确定路径基准，配置内相对文件名都从这里开始解析。
        const fs::path config = fs::weakly_canonical(fs::absolute(filename));
        // 关键词缺失时按原有开关警告并返回零或空字符串；字符串显式留空不警告。

        // =========================== 03 配置关键词集中读取 ===========================

        // 各输出文件独立指定路径；缺失或留空时关闭该项，相对路径从 setting 所在目录解析。
        // 输出目录由外部准备；程序只打开指定文件，不创建目录或自动改名。

        // ----------------------------- 3.1 选择输出文件 ------------------------------

        // 【可选】file_aposteriori_eph_Cartesian_BODN：后验 Cartesian 星历。
        // 缺失或留空时不输出；填写文件名后启用。
        std::string fn_eph_cartesian;
        input.GetString(tracer, "file_aposteriori_eph_Cartesian_BODN", fn_eph_cartesian);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选】file_aposteriori_eph_Keplerian_BODN：后验 Kepler 星历，缺失或留空时不输出。
        std::string fn_eph_keplerian;
        input.GetString(tracer, "file_aposteriori_eph_Keplerian_BODN", fn_eph_keplerian);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选】file_aposteriori_eph_NonSingular2_BODN：第二类无奇点根数星历，缺失或留空时不输出。
        // 与另外两种星历共用采样历元；a 为 km、lambda 为 deg，其余四项无量纲。
        std::string fn_eph_nsg2;
        input.GetString(tracer, "file_aposteriori_eph_NonSingular2_BODN", fn_eph_nsg2);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选】file_residual：逐条观测残差。
        // 填写文件名后启用，缺失或留空时关闭输出。
        std::string fn_residual;
        input.GetString(tracer, "file_residual", fn_residual);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选】file_report：科学结果报告，包含状态、迭代、统计和形式协方差。
        // 填写文件名后启用，缺失或留空时关闭输出。
        std::string fn_report;
        input.GetString(tracer, "file_report", fn_report);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选】file_programlog：运行日志，保存原始设置、事件和失败时的调用链。
        // 填写文件名后启用，缺失或留空时关闭输出。
        std::string fn_programlog;
        if (input.HasKeyword("file_programlog"))
            input.GetString(tracer, "file_programlog", fn_programlog);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选】file_orbit_solution_BODN：后验轨道解，含参考历元、GCRS 状态及已启用力项的 AMR。
        // 缺失或留空时不输出；可填写 bodn_orbit_solution.si，供下一次 BODN 读取。
        std::string fn_orbit_solution;
        input.GetString(tracer, "file_orbit_solution_BODN", fn_orbit_solution);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // 【可选覆盖】ephemeris_interval：后验星历采样间隔，单位 s。
        // 0 输出积分节点；正值按 UTC 秒均匀采样且须不小于 1e-5，填写时不能留空。
        double ephemeris_interval = 0;
        if (input.HasKeyword("ephemeris_interval"))
            input.GetDouble(tracer, "ephemeris_interval", ephemeris_interval);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // ----------------------------- 3.2 检查输出路径 ------------------------------

        // 两个数组按同一下标配对：前者是 setting 关键词，后者是对应输出路径。
        // 下标名称只用于后面打开各类文件；列举顺序沿用报告中的关键词排序。
        enum { CartesianEphemeris, KeplerianEphemeris, Nsg2Ephemeris, OrbitSolution,
               ProgramLog, Report, Residual, OutputFileCount };
        const char *output_keywords[OutputFileCount] = {
            "file_aposteriori_eph_Cartesian_BODN", "file_aposteriori_eph_Keplerian_BODN", "file_aposteriori_eph_NonSingular2_BODN", "file_orbit_solution_BODN",
            "file_programlog", "file_report", "file_residual"};

        // 路径数组与上面的关键词数组按下标对应，后面统一检查并打开这些输出。
        fs::path output_files[OutputFileCount] = {
            fn_eph_cartesian, fn_eph_keplerian, fn_eph_nsg2, fn_orbit_solution,
            fn_programlog, fn_report, fn_residual};
        std::array<std::string, OutputFileCount> write_messages;

        // 逐项规范化输出路径；空值表示关闭输出，不参与后面的冲突检查。
        for (int i = 0; i < OutputFileCount; ++i)
        {
            if (output_files[i].empty()) continue;
            output_files[i] = fs::weakly_canonical(output_files[i].is_absolute() ? output_files[i] :
                                                  config.parent_path() / output_files[i]);

            // 与此前已经解析的路径逐一比较，两个输出不能写到同一文件。
            for (int j = 0; j < i; ++j)
                if (output_files[i] == output_files[j])
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__,
                        "Output paths must be distinct: " + output_files[i].string(), ERROR);
                    tracer.EmergencyExit();
                }

            // 输出流按覆盖方式写入已有文件；目录不能作为文件目标。
            if (fs::is_directory(output_files[i]))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Output path is a directory: " + output_files[i].string(), ERROR);
                tracer.EmergencyExit();
            }
        }

        // 文件路径不能同时用作目录；逐层检查父目录，排除输出之间的嵌套冲突。
        for (int i = 0; i < OutputFileCount; ++i)
        {
            if (output_files[i].empty()) continue;
            for (fs::path parent = output_files[i].parent_path(); parent != parent.root_path(); parent = parent.parent_path())
                for (int j = 0; j < OutputFileCount; ++j)
                    // 任一父目录与其他输出文件同名时，文件和目录的用途发生冲突。
                    if (parent == output_files[j])
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__,
                            "Output file is also used as a directory: " + parent.string(), ERROR);
                        tracer.EmergencyExit();
                    }
        }

        // 星历间隔允许为 0（使用积分节点）；负值或不足 1e-5 s 的正间隔均不接受。
        if (ephemeris_interval < 0 ||
            (ephemeris_interval > 0 && ephemeris_interval < 1e-5))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "ephemeris_interval must be 0 or at least 1e-5 seconds.", ERROR);
            tracer.EmergencyExit();
        }
        // ----------------------------- 3.3 创建程序日志 ------------------------------

        // EmergencyExit 可直接接收日志流并负责写事件、写调用链和关闭文件。
        // 日志启用时，退出前另外在控制台显示调用链；没有日志时直接报告到 stderr。
        if (!output_files[ProgramLog].empty())
        {
            const fs::path &destination = output_files[ProgramLog];
            program_log.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!program_log)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[ProgramLog] = "Writing program log to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[ProgramLog], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        // -------------------------- 3.4 集中读取其余关键词 ---------------------------

        // 已按设置打开程序日志，以下读取失败也能保存完整事件和调用链。
        // 在本节读完所有配置，后续阶段使用这些值进行校验、单位转换和科学计算。
        OrbitPropagationJob job{};
        BatchSetting setting;
        auto &model = job.opn_model;

        // 各类资料路径先按字符串接收，读取完成后统一以 setting 目录为基准解析。
        std::string observation_file;
        std::string observer_file;
        std::string eop_file;
        std::string jpl_file;
        std::string density_file;
        std::string initial_orbit_file;

        // 观测选取范围与内联初轨保留为原始文本，后面分别检查格式和数值。
        std::string observation_count;
        std::string epoch_text, kep_text, cart_text;
        double apriori_epoch;
        int use_apriori_orbit_amr = 0;

        // 布尔开关先用整数接收；可选项从现有设置取得缺省值，再检查是否为 0 或 1。
        int stm_use_J2_value = setting.analytic_stm.use_J2 ? 1 : 0;
        int stm_use_lunisolar_value = setting.analytic_stm.use_lunisolar ? 1 : 0;
        int drms_relative;
        int compute_covariance_value = setting.compute_covariance ? 1 : 0;
        int iteration_display_value;

        // 卫星摄动力开关、质量、面积和观测标准差在读取后统一校验并转换单位。
        int use_drag_perturb_value;
        int use_SRP_perturb_value;
        double mass, drag_area, srp_area;
        double sigma_angle, sigma_xyz;

        // ------------------------ 资料路径与读取范围 ------------------------

        // 【必选】file_observation：标准 OBS 观测文件路径，不能留空。
        input.GetString(tracer, "file_observation", observation_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【条件必选】file_observer：含地基赤经赤纬观测时，必须填写测站表路径。
        input.GetString(tracer, "file_observer", observer_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】file_EOP：地球定向参数文件路径，不能留空，须覆盖初轨及观测历元。
        input.GetString(tracer, "file_EOP", eop_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【条件必选】file_JPLEPH：solar_gravity 或 lunar_gravity 为 1，或启用 SRP 时必须填写。
        input.GetString(tracer, "file_JPLEPH", jpl_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【条件必选】file_atmosphere_density：use_drag_perturb 为 1 时必须填写大气密度表路径。
        input.GetString(tracer, "file_atmosphere_density", density_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选】file_apriori_orbit：SI 格式初轨文件路径，缺失或留空时改用 setting 内的初轨。
        input.GetString(tracer, "file_apriori_orbit", initial_orbit_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】observations_to_use：读取全部记录，或只读取前 N 条观测。
        input.GetString(tracer, "observations_to_use", observation_count);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ------------------------ 初始轨道与物理参数导入 ------------------------

        // 先保留 apriori_epoch 的原始文本，用于检查文件初轨与内联初轨是否同时填写。
        input.GetString(tracer, "apriori_epoch", epoch_text);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【条件必选】使用内联初轨时，将 apriori_epoch 读为 MJD(UTC)；文件初轨的历元从 SI 读取。
        if (initial_orbit_file.empty())
        {
            input.GetDouble(tracer, "apriori_epoch", apriori_epoch);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 【条件必选】apriori_Keplerian：内联根数，km、deg，末项为平近点角 M。
        input.GetString(tracer, "apriori_Keplerian", kep_text);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【条件必选】apriori_Cartesian：内联 GCRS 状态，km、km/s，与根数二选一。
        input.GetString(tracer, "apriori_Cartesian", cart_text);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】use_apriori_orbit_amr：默认 0；设为 1 时从先验轨道文件导入可用 AMR。
        if (input.HasKeyword("use_apriori_orbit_amr"))
            input.GetInteger(tracer, "use_apriori_orbit_amr", use_apriori_orbit_amr);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ------------------------ 数值预报与偏导设置 ------------------------

        // 【必选】integrator_id：积分器编号，101/102 为自适应/定步长 RKF78，201/202 为 Adams PE/PECE。
        input.GetInteger(tracer, "integrator_id", setting.integrator_id);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】integrator_stepsize：积分步长，s；配置传播时再归一化。
        input.GetDouble(tracer, "integrator_stepsize", job.opn_setting.h0_norm);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】integrator_maxTE：RKF78 的归一化局部截断误差限。
        input.GetDouble(tracer, "integrator_maxTE", job.opn_setting.rkf.maxTE_norm);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】integrator_adams_steps：Adams 多步法使用的历史节点数。
        input.GetInteger(tracer, "integrator_adams_steps", job.opn_setting.adams.steps);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】stm_method：状态转移矩阵的求法，1 为解析近似，2 为完整传播轨迹中心差分。
        if (input.HasKeyword("stm_method"))
            input.GetInteger(tracer, "stm_method", setting.stm_method);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】stm_use_J2：解析 STM 的 J2 项开关，不改变实际传播的力模型。
        if (setting.stm_method == 1 && input.HasKeyword("stm_use_J2"))
            input.GetInteger(tracer, "stm_use_J2", stm_use_J2_value);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】stm_use_lunisolar：解析 STM 的日月平均项开关。
        if (setting.stm_method == 1 && input.HasKeyword("stm_use_lunisolar"))
            input.GetInteger(tracer, "stm_use_lunisolar", stm_use_lunisolar_value);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】stm_fd_step：中心差分的扰动尺度，采用内部归一化 NSG2 单位。
        if (setting.stm_method == 2 && input.HasKeyword("stm_fd_step"))
            input.GetDouble(tracer, "stm_fd_step", setting.fd_step);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ------------------------ 定轨与迭代设置 ------------------------

        // 【必选】iteration_maxround：允许完成的最大迭代轮数。
        input.GetInteger(tracer, "iteration_maxround", setting.max_iterations);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】iteration_minround：允许停止前至少需要完成的迭代轮数。
        input.GetInteger(tracer, "iteration_minround", setting.min_iterations);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】iteration_rms_tolerance：标准化残差 RMS 门限，无量纲。
        input.GetDouble(tracer, "iteration_rms_tolerance", setting.rms_tolerance);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】iteration_drms_tolerance：相邻轮 RMS 变化门限。
        input.GetDouble(tracer, "iteration_drms_tolerance", setting.drms_tolerance);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】iteration_drms_relative：0 为绝对差，1 为相对本轮 RMS 的变化。
        input.GetInteger(tracer, "iteration_drms_relative", drms_relative);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】rejection_scale：同类型剔除前 RMS 的倍数门限，0 表示不剔除。
        if (input.HasKeyword("rejection_scale"))
            input.GetDouble(tracer, "rejection_scale", setting.rejection_scale);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】rejection_start：从这一轮起允许剔除异常观测。
        if (setting.rejection_scale > 0 && input.HasKeyword("rejection_start"))
            input.GetInteger(tracer, "rejection_start", setting.rejection_start);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 硬窗口按类型使用内部默认；仅显式覆盖时读取，文件单位直接体现在名称中。
        if (setting.rejection_scale > 0 && input.HasKeyword("rejection_angle_window_rad"))
            input.GetDouble(tracer, "rejection_angle_window_rad", setting.rejection_angle_window_rad);
        if (setting.rejection_scale > 0 && input.HasKeyword("rejection_position_window_m"))
            input.GetDouble(tracer, "rejection_position_window_m", setting.rejection_position_window_m);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】compute_covariance：是否计算最终形式协方差，0 关闭，1 开启。
        if (input.HasKeyword("compute_covariance"))
            input.GetInteger(tracer, "compute_covariance", compute_covariance_value);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】display_iterations：是否显示逐轮过程，0 关闭，1 开启；任务报告仍保留逐轮记录。
        input.GetInteger(tracer, "display_iterations", iteration_display_value);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】estimate_parameter：附加参数模式，0 只估计轨道，1 阻力 AMR，2 光压 AMR，3 自动选择。
        if (input.HasKeyword("estimate_parameter"))
            input.GetInteger(tracer, "estimate_parameter", setting.estimate_parameter);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选覆盖】parameter_fd_step：计算 AMR 偏导时，对 AMR 与先验之比施加的差分步长。
        if (setting.stm_method == 2 && setting.estimate_parameter != 0 && input.HasKeyword("parameter_fd_step"))
            input.GetDouble(tracer, "parameter_fd_step", setting.parameter_fd_step);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ------------------------ 动力学模型与卫星参数 ------------------------

        // 【必选】earth_zonal_order：地球重力场带谐项的截断阶数。
        input.GetInteger(tracer, "earth_zonal_order", model.earth_zonal_order);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】earth_tesseral_degree：地球重力场田谐项的截断阶数。
        input.GetInteger(tracer, "earth_tesseral_degree", model.earth_tesseral_degree);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】solar_gravity：太阳引力模式，0 关闭，1 使用 JPL，2 使用解析星历。
        input.GetInteger(tracer, "solar_gravity", model.solar_gravity_config);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】lunar_gravity：月球引力模式，0 关闭，1 使用 JPL，2 使用解析星历。
        input.GetInteger(tracer, "lunar_gravity", model.lunar_gravity_config);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】use_drag_perturb：大气阻力开关，0 关闭，1 开启。
        input.GetInteger(tracer, "use_drag_perturb", use_drag_perturb_value);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】use_SRP_perturb：太阳光压开关，0 关闭，1 开启。
        input.GetInteger(tracer, "use_SRP_perturb", use_SRP_perturb_value);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】satellite_mass：卫星质量，单位 kg。
        input.GetDouble(tracer, "satellite_mass", mass);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】satellite_drag_cross_area：迎风面积，单位 m²。
        input.GetDouble(tracer, "satellite_drag_cross_area", drag_area);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】satellite_srp_cross_area：受光面积，单位 m²。
        input.GetDouble(tracer, "satellite_srp_cross_area", srp_area);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】satellite_drag_Cd：固定阻力系数 Cd，无量纲。
        input.GetDouble(tracer, "satellite_drag_Cd", model.sat.Cd);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】satellite_srp_eta：固定光压反射参数 eta，无量纲。
        input.GetDouble(tracer, "satellite_srp_eta", model.sat.eta_SRP);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ------------------------ 观测标准差 ------------------------

        // 【必选】noise_sigma_arcsecond：测角标准差，arcsec；后续转换为 rad。
        input.GetDouble(tracer, "noise_sigma_arcsecond", sigma_angle);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】noise_sigma_position_m：位置分量标准差，m。
        input.GetDouble(tracer, "noise_sigma_position_m", sigma_xyz);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ----------------------------- 3.5 解析资料路径 ------------------------------

        // 资料路径已集中读入；这里从 setting 目录解析相对路径。
        if (!observation_file.empty())
        {
            const fs::path input_path = observation_file;
            observation_file = fs::weakly_canonical(input_path.is_absolute() ? input_path :
                config.parent_path()/input_path, filesystem_error).string();
            // 资料路径解析失败时报告文件系统错误，不继续使用未确认的路径。
            if (filesystem_error)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot resolve file_observation: " + filesystem_error.message(), ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 仅含天基赤经赤纬或 GCRS 位置观测时，可留空。
        if (!observer_file.empty())
        {
            const fs::path input_path = observer_file;
            observer_file = fs::weakly_canonical(input_path.is_absolute() ? input_path :
                config.parent_path()/input_path, filesystem_error).string();
            // 资料路径解析失败时报告文件系统错误，不继续使用未确认的路径。
            if (filesystem_error)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot resolve file_observer: " + filesystem_error.message(), ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 填写了 EOP 路径时解析其相对位置；后续模型加载仍负责检查资料是否可用。
        if (!eop_file.empty())
        {
            const fs::path input_path = eop_file;
            eop_file = fs::weakly_canonical(input_path.is_absolute() ? input_path :
                config.parent_path()/input_path, filesystem_error).string();
            // 资料路径解析失败时报告文件系统错误，不继续使用未确认的路径。
            if (filesystem_error)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot resolve file_EOP: " + filesystem_error.message(), ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 其余情况下可留空。
        if (!jpl_file.empty())
        {
            const fs::path input_path = jpl_file;
            jpl_file = fs::weakly_canonical(input_path.is_absolute() ? input_path :
                config.parent_path()/input_path, filesystem_error).string();
            // 资料路径解析失败时报告文件系统错误，不继续使用未确认的路径。
            if (filesystem_error)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot resolve file_JPLEPH: " + filesystem_error.message(), ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 未启用大气阻力时，可留空。
        if (!density_file.empty())
        {
            const fs::path input_path = density_file;
            density_file = fs::weakly_canonical(input_path.is_absolute() ? input_path :
                config.parent_path()/input_path, filesystem_error).string();
            // 资料路径解析失败时报告文件系统错误，不继续使用未确认的路径。
            if (filesystem_error)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot resolve file_atmosphere_density: " + filesystem_error.message(), ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 填写路径时，apriori_epoch、apriori_Keplerian 和 apriori_Cartesian 均须留空。
        if (!initial_orbit_file.empty())
        {
            const fs::path input_path = initial_orbit_file;
            initial_orbit_file = fs::weakly_canonical(input_path.is_absolute() ? input_path :
                config.parent_path()/input_path, filesystem_error).string();
            // 资料路径解析失败时报告文件系统错误，不继续使用未确认的路径。
            if (filesystem_error)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot resolve file_apriori_orbit: " + filesystem_error.message(), ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // ------------------------ 3.6 保存配置原文与实际路径 -------------------------
        if (program_log.is_open())
        {
            program_log << "BODN program log\nStarted: " << started << "\nConfig: " << config.string()
                        << "\n\n[Original setting]\n";

            // 保存 setting 文件的原始文本，便于核对本次运行使用的完整输入。
            std::ifstream original_input(config);
            program_log << original_input.rdbuf() << "\n[Resolved paths]\n";

            // 随后列出解析后的实际资料路径，将相对文件名对应到本次读取的位置。
            program_log << "file_observation = " << observation_file << '\n';
            program_log << "file_observer = " << observer_file << '\n';
            program_log << "file_EOP = " << eop_file << '\n';
            program_log << "file_JPLEPH = " << jpl_file << '\n';
            program_log << "file_atmosphere_density = " << density_file << '\n';
            program_log << "file_apriori_orbit = " << initial_orbit_file << '\n';

            // 逐项记录所有输出路径；空路径表示对应文件没有启用。
            for (int i = 0; i < OutputFileCount; ++i)
                program_log << output_keywords[i] << " = " << output_files[i].string() << '\n';
        }
        // ============================ 04 初轨读取与归一化 ============================

        // 读取初轨及其历元，按输入单位准备 km、km/s、deg 状态，再转为程序内部的归一化量。
        // 初轨、积分器和力模型放在 job 中，定轨所需的求导方法与迭代设置放在 setting 中。
        const fs::path initial_file = initial_orbit_file;

        std::string apriori_object_id;
        double apriori_drag_amr_si = std::numeric_limits<double>::quiet_NaN();
        double apriori_srp_amr_si = std::numeric_limits<double>::quiet_NaN();

        // AMR 导入开关只允许 0 或 1，其他数值没有定义。
        if (use_apriori_orbit_amr != 0 && use_apriori_orbit_amr != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "use_apriori_orbit_amr must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        const bool import_amr = use_apriori_orbit_amr == 1;
        // 启用 AMR 导入时必须提供 SI 初轨文件，因为 AMR 将从该文件读取。
        if (import_amr && initial_file.empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "use_apriori_orbit_amr requires file_apriori_orbit.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        bool initial_is_keplerian;
        std::array<double, 6> v{};
        // 填写 SI 路径时从文件读取初轨，否则使用 setting 中的内联初轨。
        if (!initial_file.empty())
        {
            // 文件初轨与内联初轨只能选用一种；若同时填写，就提示冲突并停止。
            if (!epoch_text.empty() || !kep_text.empty() || !cart_text.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "With file_apriori_orbit, leave all inline apriori_* entries empty.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 轨道文件采用 SimpleInput，读取本次定轨需要的身份、历元、状态及可选 AMR。
            io::SimpleInput orbit_input;
            std::string orbit_filename = initial_file.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading initial orbit from: " + orbit_filename, EventType::MESSAGE, FORCE_ALERT_ON);
            orbit_input.ReadSimpleInput(tracer, orbit_filename);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 【必选】object_id、state_units、reference_epoch：目标编号、状态单位码和 UTC MJD。
            if (!orbit_input.HasKeyword("object_id") || !orbit_input.HasKeyword("state_units") ||
                !orbit_input.HasKeyword("reference_epoch"))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Orbit file requires object_id, state_units and reference_epoch.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            orbit_input.GetString(tracer, "object_id", apriori_object_id);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 单位码和历元应各含一个完整数值，额外字段或解析失败均作为输入错误。
            std::string unit_text, reference_epoch_text, trailing;
            orbit_input.GetString(tracer, "state_units", unit_text);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 参考历元单独读取为文本，下一步检查它是否只包含一个 UTC MJD 数值。
            orbit_input.GetString(tracer, "reference_epoch", reference_epoch_text);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 将单位码和历元分别解析为整数、浮点数，并拒绝尾随的额外字段。
            int orbit_state_units = 0;
            std::istringstream unit_input(unit_text), reference_epoch_input(reference_epoch_text);
            if (!(unit_input >> orbit_state_units) || (unit_input >> trailing) ||
                !(reference_epoch_input >> apriori_epoch) || (reference_epoch_input >> trailing))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid orbit units or reference_epoch.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 【必选其一】reference_Cartesian / reference_Keplerian：六个状态分量。
            // 同时提供两种表示时，使用 Cartesian。单位码 1 为归一化/rad，2 为 km、km/s、deg，3 为 m、m/s、deg。
            initial_is_keplerian = !orbit_input.HasKeyword("reference_Cartesian");
            const std::string state_keyword = initial_is_keplerian ? "reference_Keplerian" : "reference_Cartesian";

            // 检查目标编号非空、历元有限、单位码有效，并确认选中的状态字段存在。
            if (apriori_object_id.empty() || !std::isfinite(apriori_epoch) || orbit_state_units < 1 || orbit_state_units > 3 ||
                !orbit_input.HasKeyword(state_keyword))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid orbit identity, epoch, units or missing state.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 取出选定的六分量状态文本；此时仍保持轨道文件声明的单位。
            std::string state_text;
            orbit_input.GetString(tracer, state_keyword, state_text);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 用文本流逐个解析状态值，结果按文件顺序存入六元素数组 v。
            std::istringstream state_input(state_text);

            for (double& value : v)
            {
                // 任一分量缺失、无法解析或不是有限数值时，报告轨道状态输入错误。
                if (!(state_input >> value) || !std::isfinite(value))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Orbit state requires six finite values.", ERROR);
                    tracer.EmergencyExit(program_log);
                }
            }

            // 读满六个分量后再尝试读取；若仍有内容，说明状态列数多于六列。
            if (state_input >> trailing)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Orbit state requires exactly six values.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 将文件中的状态换为本段初轨转换使用的 km、km/s、deg。
            if (orbit_state_units == 1)
            {
                // 归一化 Keplerian 的半长轴恢复为 km，四个角度由 rad 转为 deg。
                if (initial_is_keplerian)
                {
                    v[0] *= ae / 1000;

                    for (int j = 2; j < 6; ++j) v[j] *= 180 / pi;
                }
                else
                {
                    // Cartesian 的位置乘 ae/1000，速度还需除以 UNITTIME，分别得到 km、km/s。
                    for (int j = 0; j < 6; ++j) v[j] *= ae / 1000 / (j < 3 ? 1 : UNITTIME);
                }
            }
            else if (orbit_state_units == 3)
            {
                // 米制 Keplerian 只需将半长轴除以 1000；其角度已经是 deg。
                if (initial_is_keplerian)
                {
                    v[0] /= 1000;
                }
                else
                {
                    // 米制 Cartesian 的位置和速度同时除以 1000，得到 km、km/s。
                    for (double& value : v) value /= 1000;
                }
            }

            // 【可选】drag_amr_m2_kg、srp_amr_m2_kg：物理面质比，单位 m²/kg；缺项时保持 NaN。
            // 后续仅在 use_apriori_orbit_amr=1 时采用文件中的有效值。
            if (orbit_input.HasKeyword("drag_amr_m2_kg"))
            {
                orbit_input.GetDouble(tracer, "drag_amr_m2_kg", apriori_drag_amr_si);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            }

            // 光压面质比与阻力面质比分别保存，文件缺少这一项时保留未提供状态。
            if (orbit_input.HasKeyword("srp_amr_m2_kg"))
            {
                orbit_input.GetDouble(tracer, "srp_amr_m2_kg", apriori_srp_amr_si);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            }

            // 日志保存本次实际采用的原始轨道文件；方法诊断不参与数值求解。
            if (program_log.is_open())
            {
                std::ifstream source(initial_file);
                program_log << "\n[Original apriori orbit]\n" << source.rdbuf() << '\n';
            }
        }
        // 未提供 SI 文件时，从 setting 内联状态解析初轨，仍采用原有 km、km/s、deg 单位。
        else
        {
            // 内联历元必须是有限数值，才能转换时标并作为传播起点。
            if (!std::isfinite(apriori_epoch))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial epoch must be finite.", ERROR);
                tracer.EmergencyExit(program_log);
            }
            // Kepler 根数与 Cartesian 状态必须二选一，不能同时填写，也不能同时留空。
            if (kep_text.empty() == cart_text.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Supply exactly one initial state representation.", ERROR);
                tracer.EmergencyExit(program_log);
            }
            initial_is_keplerian = !kep_text.empty();
            std::istringstream state_input(initial_is_keplerian ? kep_text : cart_text);
            // 依次解析六个状态分量，任何一项缺失或非法都不能作为初轨。
            for (double &x : v)
                // 若分量无法解析或不是有限数值，就停止，避免非法初态进入轨道转换。
                if (!(state_input >> x) || !std::isfinite(x))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial state needs six finite values.", ERROR);
                    tracer.EmergencyExit(program_log);
                }
            std::string extra;
            // 读完六个分量后若还有内容，就按列数错误处理，提示检查初轨格式。
            if (state_input >> extra)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Too many initial state columns.", ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // -------------------------- 4.1 初轨单位与状态转换 ---------------------------

        // 如果填写的是开普勒根数，就先读入根数，再由它得到后面需要的 Cartesian 状态。
        if (initial_is_keplerian)
        {
            // 这里使用的 NSG2 表示要求初轨满足 a>0、0≤e<1、0≤i<180°，先检查是否处于适用范围。
            if (v[0] <= 0 || v[1] < 0 || v[1] >= 1 || v[2] < 0 || v[2] >= 180)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Require elliptical orbit, 0<=i<180 deg.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 六个输入值依次为 a(km)、e、i、升交点赤经、近地点幅角和平近点角 M（角度均为 deg）。
            job.initKep_norm.Set(v[0], v[1], v[2], v[3], v[4], v[5], false);

            // 将 a 从 km 转为以地球半径为尺度的归一化长度，角度从 deg 转为 rad。
            job.initKep_norm.Normalize(tracer, ae / 1000);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 根据平近点角补齐其他近点角，然后转换为归一化的 Cartesian 位置和速度。
            job.initKep_norm.SyncAnomaly(tracer);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            job.initKep_norm.ToCartesian(tracer, job.initCart_norm);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }
        // 输入采用 Cartesian 时，先归一化位置速度，再转换为传播所需的 Kepler 根数。
        else
        {
            // 位置和速度分别按地球半径、时间尺度归一化。
            job.initCart_norm.Set(v[0], v[1], v[2], v[3], v[4], v[5], false);
            job.initCart_norm.Normalize(tracer, ae / 1000, UNITTIME);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            job.initCart_norm.ToKeplerian(tracer, job.initKep_norm);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ============================= 05 传播与定轨策略 =============================

        // 输入历元采用 UTC，而轨道积分使用 TT，因此同时保存两种时标的历元，并将步长从秒归一化。
        job.opn_setting.t0_MJD_UTC = apriori_epoch;
        job.opn_setting.t0_MJD_TT = UTC2TT(apriori_epoch);

        job.opn_setting.h0_norm /= UNITTIME;

        // 定轨过程中要多次传播轨道，这里先关闭传播进度条。每轮定轨的信息由后面的开关单独控制。
        job.opn_setting.display_progress = 0;

        // 检查积分器编号：101 为自适应 RKF78，102 为定步长 RKF78，201/202 分别为 Adams PE/PECE。
        if (setting.integrator_id != 101 && setting.integrator_id != 102 && setting.integrator_id != 201 &&
            setting.integrator_id != 202)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Unsupported integrator_id.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 步长和误差限都应大于零，Adams 节点数应在 2～9 之间；超出这些范围时就报告输入错误。
        if (job.opn_setting.h0_norm <= 0 || job.opn_setting.rkf.maxTE_norm <= 0 ||
            job.opn_setting.adams.steps < 2 || job.opn_setting.adams.steps > 9)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid integrator settings.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 根据编号设置积分方式：选择 102 时采用固定步长，选择 202 时采用校正后再评估的 PECE。
        job.opn_setting.rkf.fixed_step = setting.integrator_id == 102;
        job.opn_setting.adams.use_correction = setting.integrator_id == 202;

        // ------------------------ 5.1 状态转移矩阵与迭代控制 -------------------------

        // 下面检查求导和迭代设置。状态转移矩阵描述初态变化如何传递到其他历元，是定轨求导的一部分。
        // stm_method 为 1 时使用解析近似，为 2 时对完整传播轨迹做中心差分，得到初态扰动的影响。

        // 这两个开关只控制解析 STM 中的 J2、日月平均项，不改变实际轨道传播的摄动力。
        if (stm_use_J2_value != 0 && stm_use_J2_value != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "stm_use_J2 must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        setting.analytic_stm.use_J2 = stm_use_J2_value == 1;

        // 解析 STM 的日月项开关只能取 0 或 1，不接受其他整数。
        if (stm_use_lunisolar_value != 0 && stm_use_lunisolar_value != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "stm_use_lunisolar must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        setting.analytic_stm.use_lunisolar = stm_use_lunisolar_value == 1;

        // DRMS 只支持绝对差与相对差两种模式，对应开关 0 和 1。
        if (drms_relative != 0 && drms_relative != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "iteration_drms_relative must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        setting.drms_relative = drms_relative == 1;

        // 这两个开关分别决定是否计算最终形式协方差，以及是否显示每轮定轨的信息。
        if (compute_covariance_value != 0 && compute_covariance_value != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "compute_covariance must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        setting.compute_covariance = compute_covariance_value == 1;

        // 逐轮显示开关只接受 0 或 1，其他数值作为配置错误处理。
        if (iteration_display_value != 0 && iteration_display_value != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "display_iterations must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        setting.display_iterations = iteration_display_value == 1;

        // ============================= 06 动力学模型配置 =============================

        // 准备动力学模型。质量、Cd、eta 固定，可按配置增加一个 AMR 待估量。
        // 加载地球定向参数 EOP 表，用于地固系与惯性系转换及相关动力学计算。
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading EOP from: " + eop_file, EventType::MESSAGE, FORCE_ALERT_ON);
        job.opn_EOP.Read_EOP(tracer, eop_file);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // ---------------------------- 6.1 地球及日月引力 -----------------------------

        // 先检查已读入的地球带谐项和田谐项截断值，后面再设置日月引力。

        // 两个截断值都应在 0～80 之间，分别检查，才能独立选择带谐项和田谐项。
        if (model.earth_zonal_order < 0 || model.earth_zonal_order > 80 || model.earth_tesseral_degree < 0 ||
            model.earth_tesseral_degree > 80)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Earth gravity cutoffs must be 0..80.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 只要其中一项达到 2，就开启地球非球形引力；如果两项都小于 2，则只保留地球中心引力。
        model.use_earth_NSG = model.earth_zonal_order >= 2 || model.earth_tesseral_degree >= 2;

        // 日、月点质量引力各自选择：0 关闭，1 使用 JPL 数值历表，2 使用分析历表。

        // 日、月引力的模式都应在 0～2 之间。其他编号没有对应的位置计算方法，需要报错。
        if (model.solar_gravity_config < 0 || model.solar_gravity_config > 2 ||
            model.lunar_gravity_config < 0 || model.lunar_gravity_config > 2)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid Sun/Moon gravity mode.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 再检查大气阻力和太阳光压的开关：0 表示关闭，1 表示开启。
        if (use_drag_perturb_value != 0 && use_drag_perturb_value != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "use_drag_perturb must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        model.use_drag = use_drag_perturb_value == 1;

        // 太阳光压开关只能取 0 或 1，随后据此决定是否启用该摄动力。
        if (use_SRP_perturb_value != 0 && use_SRP_perturb_value != 1)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "use_SRP_perturb must be 0 or 1.", ERROR);
            tracer.EmergencyExit(program_log);
        }
        model.use_SRP = use_SRP_perturb_value == 1;

        // 如果日或月选择了 JPL 历表，或开启了太阳光压，就需要在这里加载 JPL 数据。
        if (model.solar_gravity_config == 1 || model.lunar_gravity_config == 1 || model.use_SRP)
        {
            // 使用已解析的 JPL 文件路径加载历表。
            const auto& fn = jpl_file;

            tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading JPL ephemeris from: " + fn, EventType::MESSAGE, FORCE_ALERT_ON);

            // JPL 星历加载失败时停止计算，保留用户选定的模型，不自动切换历表。
            if (Read_JPLEPH(fn) <= 0)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "JPL loading failed.", ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // ----------------------------- 6.2 卫星物理参数 ------------------------------

        // 接下来配置已读入的卫星物理参数。输入的质量用 kg，迎风面积和受光面积用 m²。
        model.sat.Id = "batch-target";

        // 阻力迎风面积，m²；允许为零。

        // 光压受光面积，m²；与阻力面积分别配置。

        // 质量应大于零，两种面积都不能为负。面积允许填 0，此时对应的阻力或光压作用为零。
        if (!std::isfinite(mass) || !std::isfinite(drag_area) || !std::isfinite(srp_area) ||
            mass <= 0 || drag_area < 0 || srp_area < 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid mass/area.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 质量以地球质量归一化，面积以地球半径平方归一化，供内部动力学使用。
        model.sat.mass = mass / MASS_OF_EARTH;
        model.sat.area_drag = drag_area / (ae * ae);
        model.sat.area_SRP = srp_area / (ae * ae);

        // 分别计算阻力和光压使用的归一化面积质量比。
        model.sat.AMR_drag = model.sat.area_drag / model.sat.mass;
        model.sat.AMR_SRP = model.sat.area_SRP / model.sat.mass;

        // 检查阻力系数 Cd 和光压系数 eta，它们都是无量纲参数，本次不随轨道一起估计。

        // 光压反射参数 eta，无量纲，与阻力 Cd 无关。

        // 这两个系数在当前模型中不能为负，若填入负值就报告输入错误。
        if (!std::isfinite(model.sat.Cd) || !std::isfinite(model.sat.eta_SRP) ||
            model.sat.Cd < 0 || model.sat.eta_SRP < 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Negative Cd/eta.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 只有显式启用导入时，才用 SI 中有效的 AMR 替换 setting 给出的面积质量比。
        // 质量和力系数仍取 setting；同时按导入的 AMR 更新面积，保持 AMR=A/m。
        if (import_amr)
        {
            // 启用 AMR 导入后，SI 中至少应有一个有效 AMR；否则提示缺少参数并停止。
            if (!std::isfinite(apriori_drag_amr_si) && !std::isfinite(apriori_srp_amr_si))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "SI contains no AMR to import.", ERROR);
                tracer.EmergencyExit(program_log);
            }
            // SI 提供有效阻力 AMR 时，只更新阻力面积质量比及相应面积。
            if (std::isfinite(apriori_drag_amr_si))
            {
                model.sat.AMR_drag = apriori_drag_amr_si * MASS_OF_EARTH / (ae * ae);
                model.sat.area_drag = model.sat.AMR_drag * model.sat.mass;
            }
            // SI 提供有效光压 AMR 时，只更新光压面积质量比及相应面积。
            if (std::isfinite(apriori_srp_amr_si))
            {
                model.sat.AMR_SRP = apriori_srp_amr_si * MASS_OF_EARTH / (ae * ae);
                model.sat.area_SRP = model.sat.AMR_SRP * model.sat.mass;
            }
            tracer.LogEvent(__FILE__, __func__, __LINE__,
                            "SI AMR imported (m^2/kg); mass, Cd and eta still supplied by setting.",
                            EventType::MESSAGE, FORCE_ALERT_ON);
        }

        // ----------------------------- 6.3 大气密度资料 ------------------------------

        // 如果开启了大气阻力，还需要读取密度表，并根据初轨选取参考密度。
        if (model.use_drag)
        {
            // 使用已解析的密度表路径，同时用 spacing 接收读取函数返回的表格间距。
            double spacing = 0;
            const auto& fn = density_file;

            tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading atmospheric density from: " + fn, EventType::MESSAGE, FORCE_ALERT_ON);

            // 密度表必须读取成功且至少含一条记录，否则无法按选定的模型计算阻力。
            if (Read_AtmosphereDensity(fn, job.atm_list, spacing, &model.atmosphere_Ap, &model.atmosphere_F107) <= 0 || job.atm_list.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Atmosphere loading failed.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 根据初始轨道和密度表选择参考密度，供当前阻力模型使用。
            Set_ReferenceDensity(job.atm_list, job.initKep_norm, job.refDensity);
        }
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // =========================== 07 观测资料与测站几何 ===========================

        // 读取观测资料。文件中的历元、角度和位置分别采用 MJD(UTC)、deg 和 m。
        // model_observation 中的原始资料接口解释每行的观测类型和各列内容。
        // 资料路径和读取范围已在输入区读入：observations_to_use 可以选择全部记录，也可以只取前 N 条。
        const fs::path observation_path = observation_file;
        tracer.LogEvent(__FILE__, __func__, __LINE__,
                        "Reading observation records from: " + observation_path.string(),
                        EventType::MESSAGE, FORCE_ALERT_ON);

        // 原始观测类负责读入各类型的数据列，记录仍保留文件中的角度和位置单位。
        observation::Observations bodn_observations;
        bodn_observations.ReadFile(tracer, observation_path.string());
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 取得原始记录，并先按全部资料初始化选取条数。
        const std::vector<observation::ObservationRecord> &records = bodn_observations.GetRecords();
        std::size_t selected_count = records.size();

        // 指定前 N 条时，检查文本只含十进制数字，且 N 不小于三。
        if (observation_count != "All")
        {
            if (observation_count.empty() || observation_count.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("observations_to_use must be All or an integer >= 3.");
            selected_count = std::stoull(observation_count);
            if (selected_count < 3) throw std::runtime_error("observations_to_use must be All or an integer >= 3.");
        }

        // 选取条数还须落在实际文件范围内，随后为计算观测预留空间。
        if (selected_count < 3 || selected_count > records.size())
            throw std::runtime_error("Too few observations for observations_to_use.");
        std::vector<Observation> observations;
        observations.reserve(selected_count);
        // 读取已结束；以下限制仅属于本次 BODN，而不限制原始观测类支持的文件类型。
        for (std::size_t i = 0; i < selected_count; ++i)
        {
            const observation::ObservationRecord &record = records[i];

            // 本次 BODN 仅接收地基测角、天基测角和 GCRS 位置观测，且全部记录属于同一目标。
            if (record.type != 1 && record.type != 6 && record.type != 7)
                throw std::runtime_error("BODN supports observation types 1, 6 and 7 only.");
            if (record.object_id != records[0].object_id)
                throw std::runtime_error("Batch OD supports one target per job.");

            // 观测历元须非递减；测角记录的赤纬须严格位于两极之间。
            if (i > 0 && record.mjd_utc < records[i - 1].mjd_utc)
                throw std::runtime_error("Observation epochs must be nondecreasing UTC MJD.");
            if (record.type != 7 && std::abs(record.dec_deg) >= 90)
                throw std::runtime_error("Declination must be between -90 and 90 deg.");

            // 把历元、观测类型、目标号和观测者号复制到本次定轨使用的计算记录。
            Observation obs;
            obs.mjd_utc = record.mjd_utc;
            obs.type = static_cast<ObservationType>(record.type);
            obs.object_id = record.object_id;
            obs.observer_id = record.observer_id;

            // 位置观测直接使用 GCRS XYZ(m)；测角观测在下面转成 rad。
            if (record.type == 7)
                for (int j = 0; j < 3; ++j) obs.value[j] = record.position_gcrs_m[j];
            else
            {
                obs.value[0] = record.ra_deg * pi / 180;
                obs.value[1] = record.dec_deg * pi / 180;

                // 天基测角同时复制该条记录提供的平台 GCRS 位置，单位为 m。
                if (record.type == 6)
                    for (int j = 0; j < 3; ++j) obs.observer_position_gcrs_m[j] = record.observer_position_gcrs_m[j];
            }

            // 保留选取记录的原始顺序，保证后续残差和采用标记能按下标对应。
            observations.push_back(obs);
        }

        // 准备加权计算所需的观测标准差：测角标准差从 arcsec 转为 rad，位置标准差保持为 m。
        sigma_angle = sigma_angle * pi / (180 * 3600);

        // 加权时需要用残差除以标准差，因此标准差必须大于零。
        if (sigma_angle <= 0 || sigma_xyz <= 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Observation sigmas must be positive.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // ------------------------- 7.1 观测权重与观测者位置 --------------------------

        // 为每条观测指定标准差，并准备观测方程需要的测站位置。
        // 地球姿态随时间变化，因此每条地基观测都使用自身历元的 EOP 计算测站惯性位置。
        // 先按站号保存 ITRS 下的站址（m）。同一测站的观测共用这个固定的地固系位置。
        std::map<std::string, Eigen::Vector3d> station_positions_itrs_m;

        // 先检查是否包含地基测角；只有这类观测需要读取地面站表。
        bool has_ground = false;
        for (const Observation &obs : observations)
            if (obs.type == ObservationType::GroundRaDec) has_ground = true;

        // 如果存在地基测角，就读取站表，后面通过每条观测的站号查找相应站址。
        if (has_ground)
        {
            // 配置只提供文件路径；相对路径仍以 setting 所在目录为基准。
            std::string stationFile;
            input.GetString(tracer, "file_observer", stationFile);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 出现地基观测时，站表路径必须非空。
            if (stationFile.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Ground observations require file_observer (station catalogue).", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 将站表相对路径转换为基于 setting 目录的路径。
            if (fs::path(stationFile).is_relative()) stationFile = (config.parent_path() / stationFile).string();

            // 读取站表，再按观测者编号将地基记录与站表关联。
            observation::Observers bodn_stations;
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Reading observers from: " + stationFile, EventType::MESSAGE, FORCE_ALERT_ON);
            bodn_stations.ReadFile(tracer, stationFile);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            bodn_observations.AttachObservers(tracer, bodn_stations);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 逐站准备 ITRS 位置，按站号保存，供不同观测历元重复使用。
            for (const observation::ObserverRecord &observer : bodn_stations.GetRecords())
            {
                Eigen::Vector3d position_itrs_m_evec;
                if (observer.coordinate_type == 1)
                {
                    // 纬度、经度(deg)、高度(m)转为归一化 ITRS，再恢复为 BODN 使用的 m。
                    frame::Geodetic site_geodetic(observer.geodetic_height_m, observer.geodetic_lambda_deg, observer.geodetic_phi_deg, false);
                    site_geodetic.Normalize(tracer, ae);
                    frame::ITRS site_itrs;
                    site_geodetic.ToITRS(tracer, site_itrs);
                    site_itrs.Restore(tracer, ae);

                    // 将恢复为 m 的 ITRS 坐标导出为三维向量。
                    double position_itrs_m_array[3];
                    site_itrs.ToArray(position_itrs_m_array);
                    position_itrs_m_evec = Eigen::Vector3d(position_itrs_m_array[0], position_itrs_m_array[1], position_itrs_m_array[2]);
                }
                else
                {
                    // 类型 2 已经提供 ITRS XYZ(m)，直接复制三个坐标分量。
                    position_itrs_m_evec = Eigen::Vector3d(observer.position_itrs_m[0], observer.position_itrs_m[1], observer.position_itrs_m[2]);
                }
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
                station_positions_itrs_m.emplace(observer.id, position_itrs_m_evec);
            }
        }

        // 接着逐条处理观测，根据资料类型填写标准差；地基测角还要计算相应历元的 GCRS 站址。
        for (Observation &obs : observations)
        {
            // 如果观测直接给出了惯性位置，就为 XYZ 三个分量分别赋予相同的标准差，单位为 m。
            if (obs.type == ObservationType::InertialPosition)
            {
                obs.sigma.setConstant(sigma_xyz); // XYZ 独立分量，m。
            }
            // 其余支持的观测属于测角，按角标准差赋权；地基观测还需计算惯性站址。
            else
            {
                // 东西、南北两个测角分量采用相同的角标准差，单位为 rad。
                obs.sigma.setConstant(sigma_angle); // 天球切平面两方向，rad。
                // 如果是地基测角，还要把测站转到惯性系；天基平台的惯性位置已经在观测中给出。
                if (obs.type == ObservationType::GroundRaDec)
                {
                    // 用当前观测的 observer_id 查找站址，这样同一份资料中也能包含多个测站。
                    const auto station = station_positions_itrs_m.find(obs.observer_id);
                    // 如果找不到这个站号，就无法确定观测者位置，需要报告站号并停止计算。
                    if (station == station_positions_itrs_m.end())
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown ground station: " + obs.observer_id, ERROR);
                        tracer.EmergencyExit(program_log);
                    }
                    // 按当前观测的 UTC 历元插取 EOP，得到这一时刻的地球定向参数。
                    aEOP eop;
                    job.opn_EOP.Get_EOP(tracer, obs.mjd_utc, eop);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
                    // 计算这一历元从地固系到惯性系的旋转矩阵，将 ITRS 站址转为 GCRS 位置，单位仍为 m。
                    Eigen::Matrix3d itrs_to_gcrs_emat;
                    ITRS_to_GCRS_Matrix(obs.mjd_utc, eop, itrs_to_gcrs_emat);
                    obs.observer_position_gcrs_m = itrs_to_gcrs_emat * station->second;
                }
            }
        }
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 使用 SI 初值时核对目标编号，避免把其他目标的轨道用于本组观测。
        if (!initial_file.empty() && apriori_object_id != observations.front().object_id)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Orbit object_id does not match the observation object.", ERROR);
            tracer.EmergencyExit(program_log);
        }

        // =========================== 08 初轨传递与估计历元 ===========================

        // 保留输入参考历元的先验，最终报告在同一时刻比较先验和后验。
        const orbit::Cartesian apriori_cart_norm = job.initCart_norm;
        const orbit::Keplerian apriori_kep_norm = job.initKep_norm;
        // 将初轨传递到首条观测历元，后面的轨道估计和协方差计算都以这个历元为准。
        // 输入初轨可能在首条观测之前或之后，需要用已选定的动力学模型向前或向后传播。
        // fit_epoch 仅是内部迭代历元；最终结果传回 apriori_epoch。
        const double fit_epoch = observations.front().mjd_utc;

        // 如果两个历元相差超过 10 微秒，就先传播初轨；差值更小时，直接使用给定初态。
        if (std::abs(fit_epoch - apriori_epoch) * 86400 > 1e-5)
        {
            // 先把两个历元转为 TT，再将时间差归一化。正值表示向后续时刻传播，负值表示向更早时刻传播。
            job.opn_setting.duration_norm = (UTC2TT(fit_epoch) - UTC2TT(apriori_epoch)) * 86400 / UNITTIME;

            // 用 times 和 states 接收传播得到的时间节点与 NSG2 状态，随后取最后一个状态。
            std::vector<double> times;
            std::vector<orbit::NonSingular2> states;

            // 按所选积分器调用 RKF78 或 Adams，两种情况都使用前面准备好的动力学模型。
            if (setting.integrator_id == 101 || setting.integrator_id == 102)
                OrbitPropagation_RKF78_NSG2(tracer, job, times, states);
            // 选择 Adams 时调用对应传播入口，与 RKF 分支采用相同的目标历元。
            else
                OrbitPropagation_Adams_NSG2(tracer, job, times, states);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 如果传播没有返回状态，就无法得到估计历元的初轨，需要在这里停止并报告错误。
            if (states.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Initial orbit epoch transfer failed.", ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 取最后一个状态，补齐近点角，再更新 job 中的开普勒根数和 Cartesian 初态。
            states.back().SyncAnomaly(tracer);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            states.back().ToKeplerian(tracer, job.initKep_norm);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            states.back().ToCartesian(tracer, job.initCart_norm);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 此时初态已经对应首条观测历元，接着将 job 中的 UTC、TT 起始历元也改到这一时刻。
        job.opn_setting.t0_MJD_UTC = fit_epoch;
        job.opn_setting.t0_MJD_TT = UTC2TT(fit_epoch);

        // ========================== 09 任务概览与批处理求解 ==========================

        tracer.LogEvent(__FILE__, __func__, __LINE__, "Input checks and model preparation completed.", EventType::MESSAGE, FORCE_ALERT_ON);

        // 选用资料的类型、观测者和独立历元用于本次任务概览。
        std::map<int, std::size_t> type_counts;
        std::set<std::string> observer_ids;
        std::set<double> observation_epochs;
        std::size_t scalar_components = 0;

        for (const Observation& obs : observations)
        {
            ++type_counts[static_cast<int>(obs.type)];
            observer_ids.insert(obs.observer_id);
            observation_epochs.insert(obs.mjd_utc);
            scalar_components += obs.type == ObservationType::InertialPosition ? 3 : 2;
        }

        // 保存归一化先验数组供同历元比较；显示副本恢复为 km、km/s、deg。
        std::array<double, 6> prior_cart, prior_kep, prior_cart_display, prior_kep_display;
        apriori_cart_norm.ToArray(prior_cart.data());
        apriori_kep_norm.ToArray(prior_kep.data());

        for (int j = 0; j < 6; ++j)
        {
            prior_cart_display[j] = prior_cart[j] * ae / 1000.0 / (j < 3 ? 1.0 : UNITTIME);
            prior_kep_display[j] = prior_kep[j] * (j == 0 ? ae / 1000.0 : j == 1 ? 1.0 : 180.0 / pi);
        }

        // General settings 保留本次输入的参考历元及初态；每种状态两行、每行三个分量。
        const int label_width = 40, state_width = 20;
        std::ostringstream setup;
        setup << "\n=================================== General settings ===================================\n\n"
              << std::left << std::setw(label_width) << "- Propagation Frame / Time" << ": GCRS / TT\n"
              << std::setw(label_width) << "- Input / Output Time" << ": UTC\n"
              << std::setw(label_width) << "- Target ID" << ": " << observations.front().object_id << '\n'
              << std::setw(label_width) << "- Mass [kg]" << ": " << std::defaultfloat << std::setprecision(14) << mass
              << "\n\n" << std::setw(label_width) << "- Reference epoch [UTC]" << ": "
              << std::fixed << std::setprecision(12) << apriori_epoch << '\n';

        const std::array<const char*, 2> initial_labels = {"- Initial Keplerian [km, deg]", "- Initial Cartesian [km, km/s]"};
        const std::array<const double*, 2> initial_values = {prior_kep_display.data(), prior_cart_display.data()};

        for (int representation = 0; representation < 2; ++representation)
        {
            for (int row = 0; row < 2; ++row)
            {
                if (row == 0)
                {
                    setup << std::left << std::setw(label_width) << initial_labels[representation] << ": ";
                }
                else
                {
                    setup << std::string(label_width + 2, ' ');
                }

                setup << "[ ";

                for (int j = 0; j < 3; ++j)
                {
                    if (j > 0) setup << ", ";

                    setup << std::right << std::setw(state_width) << std::defaultfloat << std::setprecision(14)
                          << initial_values[representation][row * 3 + j];
                }

                setup << " ]\n";
            }
        }

        // 观测概览按时间覆盖、观测者、记录总量及类型组成排列；计数均对应实际选用的资料。
        setup << "\n================================= Observation summary ==================================\n\n"
              << std::left << std::setw(label_width) << "- First / Last observation [UTC]" << ": "
              << std::fixed << std::setprecision(12) << fit_epoch << " / " << observations.back().mjd_utc
              << '\n' << std::setw(label_width) << "- Arc span [s]" << ": " << std::setprecision(2)
              << (observations.back().mjd_utc - fit_epoch) * 86400
              << '\n' << std::setw(label_width) << "- Unique epochs" << ": " << observation_epochs.size()
              << "\n\n" << std::setw(label_width) << "- Observers" << ": " << observer_ids.size() << " (";
        bool first_observer = true;

        for (const std::string& id : observer_ids)
        {
            if (!first_observer) setup << ' ';

            setup << id;
            first_observer = false;
        }

        setup << ")\n\n";

        // 选用部分资料时补充读取总数；全选时不重复列出同一个数量。
        if (records.size() != observations.size()) setup << std::setw(label_width) << "- Loaded records" << ": " << records.size() << '\n';

        setup << std::setw(label_width) << "- Observation records" << ": " << observations.size() << '\n';

        for (const auto& entry : type_counts)
        {
            const std::string label = "  > " + std::string(observation::ObservationTypeName(entry.first)) + " (" + std::to_string(entry.first) + ")";
            setup << std::setw(label_width) << label << ": " << entry.second << '\n';
        }

        setup << std::setw(label_width) << "- Scalar components" << ": " << scalar_components << '\n';

        // 动力学、积分器与估计设置各由已有成员函数显示，使用同一个字段列宽。
        model.Display(tracer, setup, 14, label_width);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        job.opn_setting.Display(tracer, setting.integrator_id, setup, 14, label_width);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        setting.Display(tracer, sigma_angle, sigma_xyz, setup);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        std::cerr << setup.str() << '\n';

        // ---------------------------- 9.1 执行批处理迭代 -----------------------------

        // 调用批处理接口，估计首观测历元的轨道与可选 AMR，接收同源残差、星历和联合协方差。
        BatchResult result;
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Starting batch orbit determination.", EventType::MESSAGE, FORCE_ALERT_ON);

        bool success = BatchOrbitDetermination_NSG2(tracer, job, observations, setting, result);

        // 求解错误立即经 tracer 退出；达到迭代上限不属于错误，仍保存已有结果。
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 内部迭代结束后，将后验状态及其协方差传回输入先验的参考历元。
        // 星历、观测残差和采用标记保留原始观测历元，不因输出参考历元改变而重新拟合。
        MovePosteriorToEpoch(tracer, job, setting, apriori_epoch, result);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // =========================== 10 后验结果与科学文件 ===========================

        // 先准备与先验同历元的后验状态，再分别写科学结果、运行报告与程序日志。

        // --------------------------- 10.1 创建科学结果文件 ---------------------------

        // 报告集中保存状态、迭代与协方差；星历和 SI 后验解仅在设置了文件名时写出。
        std::ofstream report_file;
        // 只有配置了科学报告路径时才创建该文件，留空表示不输出。
        if (!output_files[Report].empty())
        {
            const fs::path &destination = output_files[Report];
            report_file.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!report_file)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[Report] = "Writing BODN report to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[Report], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        std::ofstream residual_file;
        // 只有配置了逐条残差文件路径时才创建该文件，留空表示不输出。
        if (!output_files[Residual].empty())
        {
            const fs::path &destination = output_files[Residual];
            residual_file.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!residual_file)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[Residual] = "Writing residuals to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[Residual], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        std::ofstream eph;
        // 只有配置了Cartesian 星历路径时才创建该文件，留空表示不输出。
        if (!output_files[CartesianEphemeris].empty())
        {
            const fs::path &destination = output_files[CartesianEphemeris];
            eph.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!eph)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[CartesianEphemeris] = "Writing Cartesian ephemeris to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[CartesianEphemeris], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        std::ofstream eph_kep;
        // 只有配置了Kepler 星历路径时才创建该文件，留空表示不输出。
        if (!output_files[KeplerianEphemeris].empty())
        {
            const fs::path &destination = output_files[KeplerianEphemeris];
            eph_kep.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!eph_kep)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[KeplerianEphemeris] = "Writing Keplerian ephemeris to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[KeplerianEphemeris], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        std::ofstream eph_nsg2;
        // 只有配置了NSG2 星历路径时才创建该文件，留空表示不输出。
        if (!output_files[Nsg2Ephemeris].empty())
        {
            const fs::path &destination = output_files[Nsg2Ephemeris];
            eph_nsg2.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!eph_nsg2)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[Nsg2Ephemeris] = "Writing NSG2 ephemeris to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[Nsg2Ephemeris], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        std::ofstream orbit_solution;
        // 只有配置了SI 后验轨道解路径时才创建该文件，留空表示不输出。
        if (!output_files[OrbitSolution].empty())
        {
            const fs::path &destination = output_files[OrbitSolution];
            orbit_solution.open(destination);
            // 父目录不存在或文件无法打开时，报告实际目标路径并停止。
            if (!orbit_solution)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Cannot create output file: " + destination.string(), ERROR);
                tracer.EmergencyExit(program_log);
            }

            write_messages[OrbitSolution] = "Writing orbit solution to: " + destination.string();
            tracer.LogEvent(__FILE__, __func__, __LINE__, write_messages[OrbitSolution], EventType::MESSAGE, FORCE_ALERT_OFF);
        }

        // -------------------------- 10.2 后验状态与计算设置 --------------------------
        std::ostringstream report, covariance_text, history_text;
        // 报告使用等宽英文表格；分隔线只划分章节，表格内部依靠空行和数值列对齐。
        const std::string report_rule(100, '=');
        const std::string table_rule(100, '-');

        // 导出参考历元的归一化 Cartesian 后验状态，供先验与后验的差值统计使用。
        std::array<double, 6> final;
        result.fitted_cart_norm.ToArray(final.data());

        // 由同一后验状态计算 Keplerian，恢复为 km、deg，供报告和轨道解文件共同使用。
        orbit::Keplerian final_kep;
        result.fitted_cart_norm.ToKeplerian(tracer, final_kep);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        final_kep.Restore(tracer, ae / 1000);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 写入报告第 01、02 节；沿用已准备的后验 Keplerian，不为显示重复转换。
        result.PrintReport(tracer, setting, observations.front().object_id,
            initial_file.empty() ? "inline setting" : initial_file.string(), final_kep, report);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 记录本次实际采用的积分、求导和观测标准差设置。
        report << std::scientific << std::setprecision(8)
               << "\n" << report_rule << "\n03  MODELS AND SOLVER SETTINGS\n" << table_rule
               << "\n[Computation settings]\nSTM method: " << setting.stm_method
               << "\nIntegrator: " << setting.integrator_id
               << "\nStepsize (s): " << job.opn_setting.h0_norm * UNITTIME
               << "\nRKF maxTE (normalized): " << job.opn_setting.rkf.maxTE_norm
               << "\nAdams nodes: " << job.opn_setting.adams.steps

               // 测角标准差显示为 arcsec，位置标准差显示为 m；随后列出收敛与剔除门限。
               << "\nAngular sigma (arcsec): " << sigma_angle * 180 * 3600 / pi
               << "\nXYZ sigma (m): " << sigma_xyz
               << "\nConvergence: standardized RMS <= " << setting.rms_tolerance
               << " OR " << (setting.drms_relative ? "relative DRMS" : "absolute DRMS") << " <= " << setting.drms_tolerance
               << "; minimum rounds " << setting.min_iterations
               << "\nRejection RMS scale / first round: " << setting.rejection_scale << " / " << setting.rejection_start

               // 记录摄动力开关、重力场截断阶次，以及卫星的质量、面积和固定系数。
               << "\nEarth zonal/tesseral cutoffs: " << model.earth_zonal_order << ' ' << model.earth_tesseral_degree
               << "\nSun/Moon modes (0 off, 1 JPL, 2 analytical): " << model.solar_gravity_config << ' '
               << model.lunar_gravity_config << "\nDrag/SRP enabled: " << model.use_drag << ' ' << model.use_SRP
               << "\nMass (kg), prior drag/SRP area (m^2), Cd, eta: " << mass << ' '
               << model.sat.area_drag * ae * ae << ' ' << model.sat.area_SRP * ae * ae
               << ' ' << model.sat.Cd << ' ' << model.sat.eta_SRP

               // 解析 STM 的摄动开关仅在 stm_method=1 时参与偏导计算。
               << "\nAnalytic STM J2/lunisolar flags: " << setting.analytic_stm.use_J2 << ' '
               << setting.analytic_stm.use_lunisolar << " (only used for stm_method=1)\n";

        // 标明初始 AMR 的来源，以及附加参数模式最终解析出的待估参数。
        report << "Initial AMR source: " << (import_amr ? "SI AMR entries; remaining entries from setting" : "setting areas/mass")
               << "\nParameter selection (0 none, 1 drag, 2 SRP, 3 teacher auto): " << setting.estimate_parameter
               << "\nResolved additional parameter: " << result.estimated_parameter << '\n';
        // 启用光压时附上当前模型的遮影约定，便于解释实际使用的摄动力。
        if (model.use_SRP)
            report << "SRP model: eclipse disabled in existing force model; reference distance 23455 Earth radii.\n";

        // ------------------------------ 10.3 后验力参数 ------------------------------

        // 输出任务副本保存最终 AMR；轨道文件中的力参数应与后验轨道对应。
        report << "\n" << report_rule << "\n04  ADDITIONAL FORCE PARAMETER\n" << table_rule << '\n';
        auto posterior = job;
        // 联合估计了 AMR 时，将后验值写入输出任务副本，供轨道文件保存并报告参数结果。
        if (result.estimated_parameter)
        {
            const bool drag = result.estimated_parameter == 1;
            const double amr_norm = result.parameter_amr_si * MASS_OF_EARTH / (ae * ae);
            // 待估参数为阻力 AMR 时，同步输出任务副本的阻力面质比和面积。
            if (drag)
            {
                posterior.opn_model.sat.AMR_drag = amr_norm;
                posterior.opn_model.sat.area_drag = amr_norm * model.sat.mass;
            }
            // 待估参数为光压 AMR 时，同步输出任务副本的光压面质比和面积。
            else
            {
                posterior.opn_model.sat.AMR_SRP = amr_norm;
                posterior.opn_model.sat.area_SRP = amr_norm * model.sat.mass;
            }

            // 将面质比换算为报告中的阻力系数或光压加速度，并用相同因子换算形式标准差。
            const double factor = drag ? model.sat.Cd * job.refDensity.rho * MASS_OF_EARTH / (ae * ae * ae) :
                                         (1 + model.sat.eta_SRP) * 4.5605e-6;

            // 同时列出 AMR 的先验、后验和形式标准差，单位均为 m²/kg。
            report << "[Additional force parameter]\n"
                   << "Estimated AMR: " << (drag ? "drag" : "SRP") << " [m^2/kg]\n"
                   << "Prior / posterior / formal sigma: " << std::setprecision(17)
                   << result.parameter_apriori_amr_si << ' ' << result.parameter_amr_si << ' '
                   << result.parameter_sigma_amr_si

                   // 报告中注明固定力参数及差分步长，说明本次联合估计的参数范围。
                   << "\nMass, Cd, eta and density model fixed; coefficient and AMR not independently estimated.\n"
                   << "AMR correction is diagnostic; stopping uses the joint fit RMS/DRMS."
                   << " Parameter FD step (ratio) " << setting.parameter_fd_step
                   << (drag ? "\nEffective drag coefficient Cd*rho_reference*AMR [1/m] (not TLE BSTAR): " :
                              "\nSRP acceleration at model reference distance [m/s^2]: ")
                   << factor * result.parameter_amr_si << " ; formal sigma " << factor * result.parameter_sigma_amr_si << '\n';

        }

        // 只估计六维轨道时明确说明未估计附加力参数，避免把固定值误认作估计结果。
        if (!result.estimated_parameter)
            report << "No additional force parameter estimated; six-state orbit fit only.\n";

        // ----------------------------- 10.4 分站残差统计 -----------------------------

        // 按观测者编号和观测类型分组，组内只统计最终采用记录的分量 RMS。
        // 测角分别统计东西、南北方向，位置分别统计 GCRS XYZ；计数单位为观测记录。
        std::set<std::pair<std::string, int>> residual_groups;

        for (const Observation& obs : observations)
            residual_groups.emplace(obs.observer_id, static_cast<int>(obs.type));

        report << "\n" << report_rule << "\n05  RESIDUAL STATISTICS\n" << table_rule
               << "\n[Residual statistics: accepted records only]\n"
               << "Angular components: dRA*cos(observed Dec), dDec [arcsec]. XYZ: dX,dY,dZ [m].\n"

               // 正式报告每组保留三列 RMS；测角第三列和全组剔除时的数值稍后写为 NaN。
               << "  " << std::left << std::setw(16) << "Observer" << std::right
               << std::setw(5) << "Type" << std::setw(9) << "Used" << std::setw(9) << "Rejected"
               << ' ' << std::setw(22) << "RMS_1" << ' ' << std::setw(22) << "RMS_2"
               << ' ' << std::setw(22) << "RMS_3" << " Unit\n"
               << std::scientific << std::setprecision(12);

        // set 按编号、类型排序；每组按原始记录顺序累加平方和，不减均值、不按 sigma 白化。
        for (const auto& group : residual_groups)
        {
            // 取出当前组的观测者编号和类型，确定这一组是两个测角分量还是三个位置分量。
            const std::string& observer_id = group.first;
            const int type = group.second;
            const bool xyz = type == static_cast<int>(ObservationType::InertialPosition);
            const int components = xyz ? 3 : 2;

            // 每组独立清零采用数、剔除数和残差平方和。
            std::size_t used = 0, rejected = 0;
            Eigen::Vector3d square_sum = Eigen::Vector3d::Zero();

            // 遍历观测，只处理编号和类型都属于当前组的记录。
            for (std::size_t i = 0; i < observations.size(); ++i)
            {
                if (observations[i].observer_id != observer_id || static_cast<int>(observations[i].type) != type) continue;

                // 采用记录进入平方和与 RMS 分母；剔除记录只增加剔除计数。
                if (result.accepted[i])
                {
                    ++used;
                    square_sum += result.residuals[i].cwiseAbs2();
                }
                else
                {
                    ++rejected;
                }
            }

            // RMS 以采用记录数为分母；全组剔除和测角的第三分量保留 NaN。
            Eigen::Vector3d component_rms = Eigen::Vector3d::Constant(std::numeric_limits<double>::quiet_NaN());

            if (used > 0)
            {
                for (int component = 0; component < components; ++component)
                    component_rms[component] = std::sqrt(square_sum[component] / used);
            }

            // 计算结果保持 rad/m，显示测角 RMS 时换算为 arcsec。
            const double scale = xyz ? 1 : 180 * 3600 / pi;
            const char* unit = xyz ? "m" : "arcsec";

            // 先写当前组的编号、类型与记录计数，再填写各分量 RMS。
            report << "  " << std::left << std::setw(16) << observer_id << std::right
                   << std::setw(5) << type << std::setw(9) << used << std::setw(9) << rejected;

            // 有定义的分量恢复显示单位；未定义或没有采用记录的分量显示 NaN。
            for (int component = 0; component < 3; ++component)
            {
                if (component >= components || used == 0)
                {
                    report << ' ' << std::setw(22) << "NaN";
                }
                else
                {
                    report << ' ' << std::setw(22) << component_rms[component] * scale;
                }
            }

            report << ' ' << unit << '\n';

        }

        // 控制台按观测类型跨观测者汇总；直接使用最终残差、标准差和采用标记，不重新传播。
        std::array<std::size_t, 8> used_by_type{}, rejected_by_type{};
        std::array<Eigen::Vector3d, 8> residual_square_sum;
        std::array<double, 8> weighted_square_sum{};
        std::size_t used_records = 0, rejected_records = 0, used_components = 0, rejected_components = 0;

        for (Eigen::Vector3d& sum : residual_square_sum) sum.setZero();

        for (std::size_t i = 0; i < observations.size(); ++i)
        {
            const Observation& obs = observations[i];
            const int type = static_cast<int>(obs.type);
            const int components = obs.type == ObservationType::InertialPosition ? 3 : 2;

            if (result.accepted[i])
            {
                ++used_by_type[type];
                ++used_records;
                used_components += components;

                // 物理分量平方和供 RMS 使用；除以各自 sigma 的平方和供分类 WRMS 使用。
                for (int j = 0; j < components; ++j)
                {
                    const double residual = result.residuals[i][j];
                    residual_square_sum[type][j] += residual * residual;
                    const double weighted_residual = residual / obs.sigma[j];
                    weighted_square_sum[type] += weighted_residual * weighted_residual;
                }
            }
            else
            {
                ++rejected_by_type[type];
                ++rejected_records;
                rejected_components += components;
            }
        }

        // 总 WRMS 复用求解结果；记录与分量的百分比分别使用各自的完整数量。
        std::ostringstream residual_summary;
        residual_summary << '\n' << std::left << std::setw(label_width) << "- All observation types" << ": WRMS "
                         << std::defaultfloat << std::setprecision(14) << result.weighted_rms
                         << '\n' << std::setw(label_width) << "  > Used / rejected records" << ": "
                         << used_records << " (" << std::fixed << std::setprecision(2) << 100.0 * used_records / observations.size()
                         << "%) / " << rejected_records << " (" << 100.0 * rejected_records / observations.size() << "%)"
                         << '\n' << std::setw(label_width) << "  > Used / rejected components" << ": "
                         << used_components << " (" << 100.0 * used_components / scalar_components
                         << "%) / " << rejected_components << " (" << 100.0 * rejected_components / scalar_components << "%)\n";

        for (const auto& entry : type_counts)
        {
            const int type = entry.first;
            const std::size_t used = used_by_type[type], rejected = rejected_by_type[type];
            const int components = type == static_cast<int>(ObservationType::InertialPosition) ? 3 : 2;
            const std::string label = "- " + std::string(observation::ObservationTypeName(type)) + " (" + std::to_string(type) + ")";
            residual_summary << '\n' << std::setw(label_width) << label << ": WRMS " << std::defaultfloat << std::setprecision(14);

            if (used > 0)
            {
                residual_summary << std::sqrt(weighted_square_sum[type] / (used * components));
            }
            else
            {
                residual_summary << "N/A";
            }

            residual_summary << '\n' << std::setw(label_width) << "  > Used / rejected records" << ": "
                             << used << " (" << std::fixed << std::setprecision(2) << 100.0 * used / entry.second
                             << "%) / " << rejected << " (" << 100.0 * rejected / entry.second << "%)"
                             << '\n' << std::setw(label_width)
                             << (components == 2 ? "  > RMS East / North [arcsec]" : "  > RMS X / Y / Z [m]") << ": ";

            // 仅显示真实存在的分量，测角从 rad 换算为 arcsec；全组剔除时没有 RMS。
            for (int j = 0; j < components; ++j)
            {
                if (j > 0) residual_summary << " / ";

                if (used > 0)
                {
                    const double scale = components == 2 ? 180.0 * 3600.0 / pi : 1.0;
                    residual_summary << std::defaultfloat << std::setprecision(10) << std::sqrt(residual_square_sum[type][j] / used) * scale;
                }
                else
                {
                    residual_summary << "N/A";
                }
            }

            residual_summary << '\n';
        }

        // ----------------------------- 10.5 逐条观测残差 -----------------------------

        // 按原始记录顺序输出 O-C 和采用标记；测角为 rad，位置按 ae 归一化。
        if (residual_file.is_open())
        {
            for (std::size_t i = 0; i < observations.size(); ++i)
            {
                const Observation& obs = observations[i];
                const bool xyz = obs.type == ObservationType::InertialPosition;

                // 每行先写序号、UTC MJD、观测类型和采用标记。
                residual_file << std::setw(8) << i + 1 << ' ' << std::fixed << std::setprecision(12)
                              << std::setw(20) << obs.mjd_utc << ' ' << std::setw(4) << static_cast<int>(obs.type)
                              << ' ' << result.accepted[i] << std::scientific << std::setprecision(16);

                // 测角第三分量写 NaN；其余分量按测角 rad、位置除以 ae 的约定写出。
                for (int component = 0; component < 3; ++component)
                {
                    if (!xyz && component == 2)
                    {
                        residual_file << ' ' << std::setw(24) << "NaN";
                    }
                    else
                    {
                        residual_file << ' ' << std::setw(24) << result.residuals[i][component] * (xyz ? 1 / ae : 1);
                    }
                }

                // 最后附上目标号与观测者号，便于将残差定位到原始记录。
                residual_file << ' ' << obs.object_id << ' ' << obs.observer_id << '\n';
            }
        }

        // ----------------------------- 10.6 完整后验星历 -----------------------------

        // 最终接受轨道的完整积分节点已保存在 result 中，三种表示共用同一条轨迹。
        // 正间隔只在这些节点上插值；输出阶段不再执行整弧轨道预报。
        std::size_t ephemeris_rows = 0;

        // 仅在至少开启一种星历文件时准备输出历元和状态。
        if (eph.is_open() || eph_kep.is_open() || eph_nsg2.is_open())
        {
            std::vector<double> sample_epochs;
            std::vector<orbit::Cartesian> sample_cart;

            // 正间隔采用 UTC 等间隔采样；间隔为零时直接使用求解器保留的原生节点。
            if (ephemeris_interval > 0)
            {
                const double first_epoch = observations.front().mjd_utc;
                const double last_epoch = observations.back().mjd_utc;
                const double span = (last_epoch - first_epoch) * 86400;

                // 采样区间数超过一千万时要求增大间隔。
                if (span / ephemeris_interval > 10000000)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Ephemeris exceeds 10000000 rows; increase ephemeris_interval.", ERROR);
                    tracer.EmergencyExit(program_log);
                }

                // 生成 UTC 等间隔采样历元，并保留观测弧段的两个端点。
                sample_epochs.push_back(first_epoch);

                for (std::size_t i = 1; i * ephemeris_interval < span; ++i)
                {
                    const double epoch = first_epoch + i * ephemeris_interval / 86400;

                    // 到达终点时停止生成内部点，终点在循环后单独补入。
                    if (epoch >= last_epoch) break;

                    // 新历元必须严格大于上一历元；间隔低于 MJD 浮点分辨率时报告错误。
                    if (epoch <= sample_epochs.back())
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__, "ephemeris_interval is below epoch resolution.", ERROR);
                        tracer.EmergencyExit(program_log);
                    }

                    sample_epochs.push_back(epoch);
                }

                // 若等间隔网格尚未包含最后观测历元，就补入终点。
                if (last_epoch > sample_epochs.back()) sample_epochs.push_back(last_epoch);

                // 在最终积分轨迹上插值一次，结果同时供三种星历表示使用。
                sample_cart = Interpolate_CartesianOrbit(tracer, result.node_epochs_utc, result.node_cart_norm, sample_epochs);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
            }

            // 间隔为零时直接读取原生节点；两种模式都不复制整份节点数组。
            const std::vector<double>& output_epochs = ephemeris_interval == 0 ? result.node_epochs_utc : sample_epochs;
            const std::vector<orbit::Cartesian>& output_cart = ephemeris_interval == 0 ? result.node_cart_norm : sample_cart;
            ephemeris_rows = output_epochs.size();

            for (std::size_t i = 0; i < output_epochs.size(); ++i)
            {
                // Cartesian 输出使用 km、km/s；单位恢复在副本上完成。
                if (eph.is_open())
                {
                    orbit::Cartesian cart = output_cart[i];
                    cart.Restore(tracer, ae / 1000, UNITTIME);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

                    // Cartesian 行先写序号和 UTC MJD，再按 XYZ、VXYZ 顺序写六个物理量。
                    double values[6];
                    cart.ToArray(values);
                    eph << std::setw(8) << i + 1 << ' ' << std::fixed << std::setprecision(12)
                        << std::setw(20) << output_epochs[i] << std::scientific << std::setprecision(16);

                    for (double value : values) eph << ' ' << std::setw(24) << value;

                    eph << '\n';
                }

                // 原生节点已有 NSG2；均匀采样时从同一插值状态转换，a、lambda 恢复为 km、deg。
                if (eph_nsg2.is_open())
                {
                    orbit::NonSingular2 elements;

                    // 原生节点直接取已存的 NSG2；均匀采样点由对应 Cartesian 状态转换。
                    if (ephemeris_interval == 0)
                    {
                        elements = result.node_nsg2_norm[i];
                    }
                    else
                    {
                        output_cart[i].ToNonSingular2(tracer, elements);
                        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
                    }

                    // 恢复 NSG2 的半长轴和平均经度单位，保留共享的归一化节点不变。
                    elements.Restore(tracer, ae / 1000);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

                    // 将平均经度归入 [0,360) deg，统一各星历行中的周期角表示。
                    double values[6];
                    elements.ToArray(values);
                    values[5] = std::fmod(values[5], 360.0);

                    if (values[5] < 0) values[5] += 360.0;

                    // NSG2 星历按目标号、UTC MJD、a、xi、eta、h、k、lambda 的顺序写出。
                    eph_nsg2 << observations.front().object_id << ' ' << std::fixed << std::setprecision(12)
                             << std::setw(20) << output_epochs[i] << std::scientific << std::setprecision(16);

                    for (double value : values) eph_nsg2 << ' ' << std::setw(24) << value;

                    eph_nsg2 << '\n';
                }

                // Keplerian 从同一 Cartesian 状态转换，采用 km、deg 和平近点角。
                if (eph_kep.is_open())
                {
                    orbit::Keplerian elements;
                    output_cart[i].ToKeplerian(tracer, elements);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

                    // Keplerian 的半长轴恢复为 km，倾角、升交点、近地点角和平近点角恢复为 deg。
                    elements.Restore(tracer, ae / 1000);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

                    // 按目标号、UTC MJD 和六个 Kepler 根数写出，所有表示使用相同采样历元。
                    double values[6];
                    elements.ToArray(values);
                    eph_kep << observations.front().object_id << ' ' << std::fixed << std::setprecision(12)
                            << std::setw(20) << output_epochs[i] << std::scientific << std::setprecision(16);

                    for (double value : values) eph_kep << ' ' << std::setw(24) << value;

                    eph_kep << '\n';
                }
            }
        }

        // 在报告中登记实际星历行数和采样方式，未开启星历时行数为零。
        report << "\n" << report_rule << "\n06  POSTERIOR EPHEMERIS\n" << table_rule
               << "\n[Posterior ephemeris]\nRecords: " << ephemeris_rows
               << "\nSampling: " << (ephemeris_rows == 0 ? "not written" :
                                     ephemeris_interval == 0 ? "integration nodes" : "uniform UTC seconds, with both endpoints")
               << "\nRequested interval (s): " << ephemeris_interval << '\n';

        // ---------------------------- 10.7 后验轨道解文件 ----------------------------

        // 后验解使用与 IOD 相同的精简 SimpleInput 格式；仅在填写路径时输出，供后续按路径传入。
        if (orbit_solution.is_open())
        {
            // 后验状态与输入先验采用同一 UTC 历元，轨道文件使用 km、km/s、deg。
            orbit::Cartesian output_cart = result.fitted_cart_norm;
            output_cart.Restore(tracer, ae / 1000, UNITTIME);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 导出两种状态表示的六个分量，准备写入同一份 SimpleInput 轨道文件。
            double output_cart_values[6], output_kep_values[6];
            output_cart.ToArray(output_cart_values);
            final_kep.ToArray(output_kep_values);

            // 分别统计已采用测角和位置观测的分量 RMS，单位为 rad、m。
            double angle_square_sum = 0, position_square_sum = 0;
            std::size_t angle_count = 0, position_count = 0;

            // 只累计最终采用的记录；位置计三个标量分量，测角计两个标量分量。
            for (std::size_t i = 0; i < observations.size(); ++i)
            {
                if (!result.accepted[i]) continue;

                if (observations[i].type == ObservationType::InertialPosition)
                {
                    position_square_sum += result.residuals[i].squaredNorm();
                    position_count += 3;
                }
                else
                {
                    // 测角使用已算好的东西、南北方向残差，其平方和不再乘 cos(Dec)。
                    angle_square_sum += result.residuals[i].head<2>().squaredNorm();
                    angle_count += 2;
                }
            }

            // SimpleInput 字段、输出顺序和 17 位有效数字在本段直接指定。
            orbit_solution << std::defaultfloat << std::setprecision(17)
                           << "object_id = " << observations.front().object_id
                           << "\nstate_units = 2\nreference_epoch = " << result.fit_epoch_utc
                           << "\nreference_Keplerian =";

            // 先写入以平近点角 M 结尾的 Keplerian 六分量，再写入 Cartesian 六分量。
            for (double value : output_kep_values) orbit_solution << ' ' << value;

            orbit_solution << "\nreference_Cartesian =";

            for (double value : output_cart_values) orbit_solution << ' ' << value;

            // 观测中实际包含相应类型时，才写入其物理量单位下的合成 RMS。
            if (angle_count) orbit_solution << "\nrms_angle_rad = " << std::sqrt(angle_square_sum / angle_count);

            if (position_count) orbit_solution << "\nrms_position_m = " << std::sqrt(position_square_sum / position_count);

            // 保存实际启用力项的后验面质比，包括固定的力参数，单位 m²/kg。
            if (model.use_drag)
            {
                const double drag_amr_si = posterior.opn_model.sat.AMR_drag * ae * ae / MASS_OF_EARTH;

                if (std::isfinite(drag_amr_si)) orbit_solution << "\ndrag_amr_m2_kg = " << drag_amr_si;
            }

            // 光压面质比独立写出；启用但未估计时保存设置中的固定值。
            if (model.use_SRP)
            {
                const double srp_amr_si = posterior.opn_model.sat.AMR_SRP * ae * ae / MASS_OF_EARTH;

                if (std::isfinite(srp_amr_si)) orbit_solution << "\nsrp_amr_m2_kg = " << srp_amr_si;
            }

            orbit_solution << '\n';

            // 检查轨道文件写出状态，流错误由 tracer 记录并结束运行。
            if (!orbit_solution)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing orbit file.", ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // ------------------------- 10.8 迭代记录与形式协方差 -------------------------

        // 矩阵按空格分列、换行分行，保留 NSG2 归一化协方差和 Cartesian SI 协方差。
        const Eigen::IOFormat matrix_format(Eigen::StreamPrecision, 0, " ", "\n", "  ", "");
        covariance_text << "Reference epoch MJD UTC " << std::fixed << std::setprecision(12) << result.fit_epoch_utc
                        << "\n\nNSG2 normalized covariance\n" << std::scientific << std::setprecision(16)
                        << result.covariance_nsg2.format(matrix_format)
                        << "\n\nGCRS Cartesian SI covariance (m, m/s)\n"
                        << result.covariance_cart_si.format(matrix_format) << '\n';
        history_text << std::scientific << std::setprecision(16);
        // 按实际迭代顺序记录 WRMS、改正量和收敛诊断，便于复核停止过程。
        for (const auto &record : result.iterations)
        {
            history_text << std::setw(5) << record.iteration << std::setw(9) << record.used_observations
                         << std::setw(6) << record.rank << ' ' << std::setw(24) << record.weighted_rms
                         << ' ' << std::setw(24) << record.scaled_correction
                         << ' ' << std::setw(24) << record.step_fraction;
            // 联合估计 AMR 时增加参数列；只估计轨道时保留六状态记录格式。
            if (result.estimated_parameter)
                history_text << ' ' << std::setw(24) << record.parameter_amr_si
                             << ' ' << std::setw(24) << record.parameter_correction;

            // 每轮最后保存 DRMS 及两个收敛标记；这些标记仍需结合最少轮次和剔除启动条件判断。
            history_text << ' ' << std::setw(24) << record.drms << ' ' << record.rms_converged
                         << ' ' << record.drms_converged << '\n';
        }

        // 先写迭代记录的列定义和数值，再说明各诊断量的含义。
        report << "\n" << report_rule << "\n07  ITERATION HISTORY\n" << table_rule
               << "\n" << result.iteration_report
               << "\n[Iteration history]\n"
               << std::setw(5) << "round" << std::setw(9) << "used" << std::setw(6) << "rank"
               << ' ' << std::setw(24) << "WRMS" << ' ' << std::setw(24) << "max_dz_normalized"
               << ' ' << std::setw(24) << "step_fraction"
               << (result.estimated_parameter ? "                AMR_m2_kg     abs_dAMR_over_prior" : "")
               << "                     DRMS RMS_met DRMS_met" << '\n' << history_text.str()

               // WRMS 为改正后的无量纲统计量；原始改正量与实际采用的线搜索比例分别记录。
               << "\nWRMS: post-correction, dimensionless. max_dz_normalized: unshortened six-state correction.\n"
               << "step_fraction: applied line-search fraction; AMR correction: unshortened abs(dAMR)/prior_AMR.\n"
               << "DRMS: change between completed rounds; NaN if first round or accepted records changed.\n"
               << "RMS_met/DRMS_met: diagnostic flags before minimum-round/rejection-start gates.\n"

               // 协方差说明明确分量顺序和 SI 各子块的单位，避免把归一化方差当作物理方差。
               << "\n" << report_rule << "\n08  FORMAL COVARIANCE\n" << table_rule
               << "\n[Formal covariance at reference epoch]\n"
               << "NSG2 order: a,xi,eta,h,k,lambda (normalized). Cartesian order: x,y,z,vx,vy,vz.\n"
               << "SI covariance blocks: position m^2, cross m^2/s, velocity m^2/s^2.\n";
        // 启用形式协方差计算时写出已准备的 NSG2 与 Cartesian 协方差。
        if (setting.compute_covariance)
            report << covariance_text.str();
        // 未计算协方差时注明这一点，不输出默认矩阵。
        else
            report << "Not computed.\n";
        // 计算了轨道与 AMR 的联合协方差时，再输出完整矩阵及参数顺序。
        if (setting.compute_covariance && result.estimated_parameter)
            report << "\nJoint covariance order: normalized a,xi,eta,h,k,lambda, AMR[m^2/kg].\n"
                   << "Orbit blocks above are marginal blocks from this joint fit.\n"
                   << std::scientific << std::setprecision(16)
                   << result.covariance_full.format(matrix_format) << '\n';

        // 集中说明各数值文件的列顺序、时间系统和坐标系，具体文件按是否启用列出。
        report << "\n" << report_rule << "\n09  DATA FORMATS\n" << table_rule
               << "\n[Numerical file layouts]\n"
               << "Whitespace-separated data, without headers. Epochs: MJD UTC; frame: GCRS.\n";
        // 生成 Cartesian 星历时，在报告中列出对应的列序和单位。
        if (eph.is_open())
            report << "file_aposteriori_eph_Cartesian_BODN: index MJD_UTC x_km y_km z_km vx_km_s vy_km_s vz_km_s\n";
        // 开启 Kepler 星历时说明其列序、单位及平近点角约定。
        if (eph_kep.is_open())
            report << "file_aposteriori_eph_Keplerian_BODN: object MJD_UTC a_km e i_deg OMEGA_deg omega_deg M_deg\n";
        // 开启残差文件时说明分量定义、采用标记和观测类型对应关系。
        if (residual_file.is_open())
            report << "file_residual: index MJD_UTC type accepted r1 r2 r3 object observer\n"
                   << "  O-C; types 1/6: wrapped dRA*cos(observed Dec), dDec, NaN [rad]; type 7: dX/ae,dY/ae,dZ/ae [1].\n"
                   << "  accepted: 1 used, 0 rejected; original observation order retained.\n";
        // 开启 NSG2 星历时说明六分量及其几何定义，避免与其他根数混淆。
        if (eph_nsg2.is_open())
            report << "file_aposteriori_eph_NonSingular2_BODN: object MJD_UTC a_km xi eta h k lambda_deg\n"
                   << "  xi=e*cos(OMEGA+omega), eta=e*sin(OMEGA+omega).\n"
                   << "  h=sin(i/2)*cos(OMEGA), k=sin(i/2)*sin(OMEGA), lambda=OMEGA+omega+M in [0,360) deg.\n";
        // 开启 SI 后验轨道解时说明文件协议、坐标系、时标和参数单位。
        if (orbit_solution.is_open())
            report << "file_orbit_solution_BODN: SimpleInput; UTC MJD, GCRS, state_units=2 (km/deg); optional AMR in m^2/kg.\n";
        report << "\n" << report_rule << "\n10  OUTPUT FILES\n" << table_rule
               << "\n[Output files]\n";
        // 遍历所有可配置输出，报告中只列出实际启用的文件路径。
        for (int i = 0; i < OutputFileCount; ++i)
            if (!output_files[i].empty())
                report << output_keywords[i] << " = " << output_files[i].string() << '\n';
        report << '\n' << report_rule << "\nEND OF BODN REPORT\n" << report_rule << '\n';
        // 报告文件已打开时写入整理好的文本；未启用报告时跳过这一步。
        if (report_file.is_open())
            report_file << report.str();

        // ------------------------ 10.9 写出结果并检查文件状态 ------------------------

        // 所有启用的科学文件都完成刷新和关闭后，才允许显示正常结束。
        for (auto* stream : {&report_file, &residual_file, &eph, &eph_kep, &eph_nsg2, &orbit_solution})
        {
            if (!stream->is_open()) continue;

            stream->flush();
            stream->close();

            if (!*stream)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Result write failed.", ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // =========================== 11 完成摘要与日志结束 ===========================

        // 结果成员输出状态、后验轨道及形式标准差；主程序接续同一分区中的统计和资源信息。
        std::ostringstream summary;
        result.PrintSummary(tracer, setting, final_kep, summary);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        double position_change2 = 0, velocity_change2 = 0;

        for (int j = 0; j < 6; ++j)
        {
            const double delta = (final[j] - prior_cart[j]) * ae / (j < 3 ? 1 : UNITTIME);

            if (j < 3)
            {
                position_change2 += delta * delta;
            }
            else
            {
                velocity_change2 += delta * delta;
            }
        }

        // 同历元的轨道改变量用一行表示，随后追加已从保存残差取得的总体及分类统计。
        // UTF-8 的减号和两个 Δ 比 ASCII 多四个字节，字段列宽相应补足。
        summary << '\n' << std::left << std::setw(label_width + 4) << "- Posterior − prior |Δr| / |Δv|" << ": "
                << std::defaultfloat << std::setprecision(14) << std::sqrt(position_change2) << " m / "
                << std::sqrt(velocity_change2) << " m/s\n" << residual_summary.str();

        // 峰值常驻内存采用系统给出的进程统计：macOS 为 byte，Linux 为 KiB。
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

        summary << '\n' << std::setw(label_width) << "- Total elapsed time [s]" << ": " << std::fixed << std::setprecision(2)
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

        // 路径按星历、残差、轨道解、报告、日志排序；事件此前已在真实写入位置记录。
        std::ostringstream file_messages;

        for (int index : {CartesianEphemeris, KeplerianEphemeris, Nsg2Ephemeris, Residual, OrbitSolution, Report, ProgramLog})
        {
            if (write_messages[index].empty()) continue;

            file_messages << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << write_messages[index] << '\n';
        }

        const std::string completed = success ? "BODN job completed successfully." : "BODN job finished without convergence.";
        tracer.LogEvent(__FILE__, __func__, __LINE__, completed, EventType::MESSAGE, FORCE_ALERT_OFF);

        if (program_log.is_open())
        {
            tracer.ReportEvents(program_log);
            program_log << "\nExit status: " << (success ? 0 : 2) << '\n';
            program_log.close();

            if (!program_log)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Program log write failed: " + output_files[ProgramLog].string(), ERROR);
                tracer.EmergencyExit();
            }
        }

        std::cerr << summary.str() << '\n' << file_messages.str();
        std::cerr << std::setfill('0') << std::setw(5) << tracer.GetId() << "> " << completed << '\n' << std::setfill(' ');

        // 收敛时返回 0，迭代上限时返回 2；已有计算错误仍由原错误出口处理。
        tracer.CheckOut(__func__);
        return success ? 0 : 2;
    }

    // ============================= 12 异常边界与错误退出 =============================

    // 将标准库或外部接口异常交给 tracer；普通输入、模型和输出校验不走异常跳转。
    catch (const std::exception &error)
    {
        // 若异常尚未由 tracer 记录，就补充一条事件；底层已记录的错误不再重复报告。
        if (tracer.Code() == ERR_NOERROR)
            tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
        tracer.EmergencyExit(program_log);
    }
}
