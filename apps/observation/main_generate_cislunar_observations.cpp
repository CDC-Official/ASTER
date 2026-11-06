/**
 * @file main_generate_cislunar_observations.cpp
 * @brief 读取两颗星的状态序列，生成正反方向的单程测距观测。
 * @details 依次读取设置、状态及接收事件，准备 JPL，计算理论伪距，再加入系统差和白噪声。
 * 状态列为年、月、日、时、分、秒、x、y、z、vx、vy、vz、钟差，采用 TDB/E-BCRS 和 m、m/s、s。
 * 接收事件列为年、月、日、时、分、秒、发射星编号、接收星编号；历元为真实 TDB 接收时刻。
 * 输出 observations.txt 保存接收钟读数和伪距，components.csv 保存对应的物理分项。
 * 用法：aster_generate_cislunar_observations -I <配置>；相对路径以配置目录为基准。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include "io_settings.hpp"
#include "math_backend.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

using namespace std;
using namespace bastion;

/**
 * @brief 读取双星状态与接收事件，生成两种发收方向的单程测距及物理分项。
 * @param[in] argc 命令行参数个数，包含程序名。
 * @param[in] argv -I后接配置路径；--version显示完整版本和构建信息，--help显示用法。
 * @return 正常完成返回0，缺少参数返回2，输入与计算错误通过tracer记录并失败退出。
 * @details 状态采用TDB/E-BCRS及m、m/s、s；输出接收钟读数、伪距和各物理分项。
 */
int main(int argc, char* argv[])
{
    // 建立主程序的调用追踪，运行事件和错误写入 program_log。
    Overwatch tracer(__FILE__, __func__, __LINE__);
    ofstream program_log;

    // 请求版本号时，输出与其他入口一致的版本、数学后端及源码提交。
    if (argc == 2 && string(argv[1]) == "--version")
    {
        cout << "aster_generate_cislunar_observations " << ASTER_PROGRAM_VERSION << '\n';
        cout << aster::math_backend::BuildDescription() << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }

    // 未给参数或请求帮助时，显示配置文件的传入方式。
    if (argc == 1 || (argc == 2 && (string(argv[1]) == "--help" || string(argv[1]) == "-h")))
    {
        cout << "Usage: aster_generate_cislunar_observations -I <settings>\n";
        tracer.CheckOut(__func__);
        return argc == 1 ? 2 : 0;
    }

    try
    {
        // =================================== 01 读取设置 ===================================
        // 正常运行使用 -I 指定一份配置文件。
        if (argc != 3 || string(argv[1]) != "-I")
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected -I <settings>.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 以配置文件所在目录解释各项相对路径。
        namespace fs = filesystem;
        const fs::path config = fs::absolute(argv[2]);
        // 读取关键词设置，后续按输入资料、模型和输出分别取值。
        io::SimpleInput settings;
        string config_name = config.string();
        settings.ReadSimpleInput(tracer, config_name);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】两份状态文件、接收事件文件、JPL 二进制历表和输出目录。
        // 所有相对路径以配置目录为基准，输出目录内生成观测、分项及日志三个文件。
        const array<string, 5> keywords = {"file_navigation_state", "file_user_state", "file_receive_events", "file_JPLEPH", "directory_output"};
        array<fs::path, 5> paths;

        for (size_t i = 0; i < keywords.size(); ++i)
        {
            // 取出当前关键词对应的路径字符串。
            string value;
            settings.GetString(tracer, keywords[i], value);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 这五项路径均为必填项，缺项时报告关键词名称。
            if (value.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Missing path: " + keywords[i], EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 绝对路径直接使用，相对路径接到配置目录下。
            paths[i] = fs::path(value).is_absolute() ? fs::path(value) : config.parent_path() / value;
        }

        // 【必选】navigation_id、user_id：两颗星的正整数编号，必须互不相同。
        array<int, 2> ids{};
        settings.GetIntegerStrict(tracer, "navigation_id", ids[0]);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 第二个编号对应用户星，两个编号用于确定每条观测的发收方向。
        settings.GetIntegerStrict(tracer, "user_id", ids[1]);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【可选】gravity_delay：0 关闭、1 开启；缺省为 1，天体为除地球外的十个主要天体。
        // 【可选】systematic_bias_m、noise_sigma_m：系统差和白噪声标准差，单位 m，缺省均为 0。
        // 【可选】noise_seed：非负整数随机种子，缺省为 1。
        int gravity_delay = 1, noise_seed = 1;
        double systematic_bias_m = 0.0, noise_sigma_m = 0.0;

        // 配置给出引力时延开关时，更新其默认值。
        if (settings.HasKeyword("gravity_delay"))
        {
            settings.GetIntegerStrict(tracer, "gravity_delay", gravity_delay);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 读取统一加到各条伪距上的系统差。
        if (settings.HasKeyword("systematic_bias_m"))
        {
            settings.GetDoubleStrict(tracer, "systematic_bias_m", systematic_bias_m);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 读取白噪声标准差，零表示本次仿真不加白噪声。
        if (settings.HasKeyword("noise_sigma_m"))
        {
            settings.GetDoubleStrict(tracer, "noise_sigma_m", noise_sigma_m);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 固定随机种子，使同一平台上的重复运行得到相同噪声序列。
        if (settings.HasKeyword("noise_seed"))
        {
            settings.GetIntegerStrict(tracer, "noise_seed", noise_seed);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);
        }

        // 检查两颗星的编号、模型开关和噪声参数的取值。
        if (ids[0] <= 0 || ids[1] <= 0 || ids[0] == ids[1] || (gravity_delay != 0 && gravity_delay != 1)
            || !isfinite(systematic_bias_m) || !isfinite(noise_sigma_m) || noise_sigma_m < 0.0 || noise_seed < 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid satellite IDs, gravity switch or noise settings.", EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 创建输出目录，随后检查三个输出文件是否已经存在。
        fs::create_directories(paths[4]);

        for (const string name : {"observations.txt", "components.csv", "run.log"})
        {
            // 已有同名结果时报告路径，由配置指定新的输出目录。
            if (fs::exists(paths[4] / name))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Output already exists: " + (paths[4] / name).string(), EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }
        }

        // 打开运行日志，接下来的读取与计算事件统一写入该文件。
        program_log.open(paths[4] / "run.log");
        if (!program_log)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot create log: " + paths[4].string(), EventType::ERROR);
            tracer.EmergencyExit();
        }

        // =================================== 02 读取两颗星的状态 ===================================
        // states[0] 保存导航星，states[1] 保存用户星，两者采用相同的状态记录。
        array<vector<Cislunar::StateRecord>, 2> states;

        for (size_t satellite = 0; satellite < states.size(); ++satellite)
        {
            // 按配置中的顺序打开当前卫星的状态文件。
            ifstream input(paths[satellite]);
            if (!input)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot read state file: " + paths[satellite].string(), EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 保存当前文本行和行号，读取错误可定位到具体文件行。
            string line;
            size_t line_number = 0;

            while (getline(input, line))
            {
                // 删除 # 后的说明文字，再处理这一行的数据。
                ++line_number;
                const size_t comment = line.find('#');

                if (comment != string::npos) line.erase(comment);

                // 跳过行首空白，空行和纯注释行直接进入下一轮。
                istringstream row(line);
                row >> ws;

                if (row.eof()) continue;

                // 按固定列序读取绝对 TDB 历元、六维状态和总钟差。
                Cislunar::StateRecord state;
                string extra;
                row >> state.year >> state.month >> state.day >> state.hour >> state.minute >> state.second
                    >> state.position_ebcrs_m[0] >> state.position_ebcrs_m[1] >> state.position_ebcrs_m[2]
                    >> state.velocity_ebcrs_mps[0] >> state.velocity_ebcrs_mps[1] >> state.velocity_ebcrs_mps[2] >> state.clock_offset_s;

                // 一行应恰好包含十三列，缺列或多列时报告所在行。
                if (!row || (row >> extra))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected 13 state columns: " + paths[satellite].string()
                        + ":" + to_string(line_number), EventType::ERROR);
                    tracer.EmergencyExit(program_log);
                }

                // 将读出的状态按文件顺序加入该星的采样序列。
                states[satellite].push_back(state);
            }

            // 确认读取正常结束，并且已有至少两个插值节点。
            if (!input.eof() || states[satellite].size() < 2)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Incomplete state file: " + paths[satellite].string(), EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 在日志中记下当前状态资料的路径。
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Loaded state records: " + paths[satellite].string(), EventType::MESSAGE);
        }

        // =================================== 03 准备传播模型及输出 ===================================
        // 载入 JPL 历表，供地球位置平移和引力时延计算共用。
        const int jpl_status = Read_JPLEPH(paths[3].string());
        if (jpl_status < 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot load JPL ephemeris: " + paths[3].string(), EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 引力开关开启时，计入除地球外的十个天体；关闭时使用空列表。
        vector<int> gravity_bodies;

        if (gravity_delay == 1) gravity_bodies = {1, 2, 4, 5, 6, 7, 8, 9, 10, 11};

        // 先生成标准正态随机数，再在每条观测中乘以指定的标准差。
        mt19937 generator(static_cast<unsigned>(noise_seed));
        normal_distribution<double> standard_normal(0.0, 1.0);
        // 打开接收事件输入，以及观测量和计算分项两个输出文件。
        ifstream events(paths[2]);
        ofstream observations(paths[4] / "observations.txt"), components(paths[4] / "components.csv");

        if (!events || !observations || !components)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open receive events or output files: "
                + paths[2].string() + ", " + paths[4].string(), EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 写明观测列的含义，并在分项文件中记录实际使用的历表和引力天体。
        observations << "# receiver clock: year month day hour minute second; transmitter receiver pseudorange_m\n";
        components << "# jpl_de=" << JPLEPH_R1.r1.numde << "; GM=DE440; gravity_bodies=";

        for (int body : gravity_bodies) components << body << ' ';

        // 分项列依次排列事件历元、传播与钟差贡献、伪距和迭代结果。
        components << "\nrecord,tx,rx,rx_year,rx_month,rx_day,rx_hour,rx_minute,rx_second,"
            << "tx_year,tx_month,tx_day,tx_hour,tx_minute,tx_second,light_time_s,geometric_m,gravity_s,"
            << "tx_clock_s,rx_clock_s,clock_m,theoretical_pseudorange_m,systematic_bias_m,noise_m,pseudorange_m,iterations,update_s\n";
        // 分项使用十七位有效数字输出，文本行号和观测序号分别计数。
        components << setprecision(17);
        string line;
        size_t line_number = 0, record_number = 0;

        // =================================== 04 逐条计算观测 ===================================
        while (getline(events, line))
        {
            // 去掉事件行中的说明文字，空行和纯注释行均跳过。
            ++line_number;
            const size_t comment = line.find('#');

            if (comment != string::npos) line.erase(comment);

            istringstream row(line);
            row >> ws;

            if (row.eof()) continue;

            // 输入时刻是真实接收 TDB，接收端位置和钟差随后从状态表插值得到。
            Cislunar::StateRecord receive;
            int transmitter_id = 0, receiver_id = 0;
            string extra;
            row >> receive.year >> receive.month >> receive.day >> receive.hour >> receive.minute >> receive.second
                >> transmitter_id >> receiver_id;

            // 核对八列事件数据，并检查发收编号是否构成这两颗星之间的一个方向。
            if (!row || (row >> extra) || !((transmitter_id == ids[0] && receiver_id == ids[1])
                || (transmitter_id == ids[1] && receiver_id == ids[0])))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid receive event: " + paths[2].string()
                    + ":" + to_string(line_number), EventType::ERROR);
                tracer.EmergencyExit(program_log);
            }

            // 根据发射编号选择状态序列，另一颗星就是本次接收端。
            const size_t transmit_index = transmitter_id == ids[0] ? 0 : 1;
            const size_t receive_index = 1 - transmit_index;
            // 逐条记下观测序号和方向，计算出错时可回到这条输入记录。
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Computing event " + to_string(++record_number)
                + ": " + to_string(transmitter_id) + " -> " + to_string(receiver_id), EventType::MESSAGE);

            // 在真实接收时刻取得接收端状态；发射端状态由核心迭代按需查询。
            Cislunar::InterpolateState(tracer, states[receive_index], receive);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 将本次发射星的状态序列交给查询函数，供光行时迭代反复调用。
            const auto query_transmit = [&](Overwatch& query_tracer, Cislunar::StateRecord& state)
            {
                Cislunar::InterpolateState(query_tracer, states[transmit_index], state);
            };
            // 从已知接收事件出发，求解发射事件并计算理论伪距。
            Cislunar::OneWayRangeRecord record;
            Cislunar::ComputeOneWayRange(tracer, receive, query_transmit, gravity_bodies, record);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            // 设备系统差和白噪声以米为单位叠加到理论伪距上。
            record.transmitter_id = transmitter_id;
            record.receiver_id = receiver_id;
            record.systematic_bias_m = systematic_bias_m;

            // 启用白噪声时，为当前观测抽取一个随机距离误差。
            if (noise_sigma_m > 0.0) record.noise_m = noise_sigma_m * standard_normal(generator);

            // 把系统差和本次随机误差加到理论伪距上，得到最终仿真观测。
            record.pseudorange_m = record.theoretical_pseudorange_m + record.systematic_bias_m + record.noise_m;

            // 观测文件保存接收钟读数；分项文件分别保存两个事件的真实 TDB 历元。
            for (int field : record.receive_clock_calendar) observations << field << ' ';

            // 接收钟读数的小数秒与伪距均以十二位小数输出。
            observations << fixed << setprecision(12) << record.receive_clock_second << ' '
                << record.transmitter_id << ' ' << record.receiver_id << ' ' << record.pseudorange_m << '\n';
            // 分项行先写观测序号和发收编号，便于与观测文件逐行对应。
            components << record_number << ',' << record.transmitter_id << ',' << record.receiver_id << ',';

            // 依次写出真实接收历元和真实发射历元，均使用 TDB 年月日时分秒。
            for (const Cislunar::StateRecord* state : {&record.receive_state, &record.transmit_state})
            {
                components << state->year << ',' << state->month << ',' << state->day << ',' << state->hour << ','
                    << state->minute << ',' << state->second << ',';
            }

            // 写出传播、钟差和误差分项，最后记录迭代次数与末次光行时变化。
            components << record.light_time_s << ',' << record.geometric_m << ',' << record.gravity_delay_s << ','
                << record.transmit_state.clock_offset_s << ',' << record.receive_state.clock_offset_s << ','
                << record.clock_contribution_m << ',' << record.theoretical_pseudorange_m << ','
                << record.systematic_bias_m << ',' << record.noise_m << ',' << record.pseudorange_m << ','
                << record.iterations << ',' << record.light_time_update_s << '\n';
        }

        // =================================== 05 完成输出并保存日志 ===================================
        // 确认事件文件正常读完，且本次至少生成了一条观测。
        if (!events.eof() || record_number == 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Incomplete or empty receive events: " + paths[2].string(), EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 关闭两个结果文件，将输出缓冲写入磁盘。
        observations.close();
        components.close();

        // 检查最终写入状态，失败时报告输出目录。
        if (!observations || !components)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing results: " + paths[4].string(), EventType::ERROR);
            tracer.EmergencyExit(program_log);
        }

        // 报告生成条数和结果目录，再把完整事件记录保存到日志。
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Generated " + to_string(record_number) + " observations in "
            + paths[4].string(), EventType::MESSAGE, FORCE_ALERT_ON);
        tracer.ReportEvents(program_log);
        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const exception& error)
    {
        // 文件系统等标准库操作抛出异常时，将错误纳入同一条 Tracer 调用追踪。
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
        tracer.EmergencyExit(program_log);
        return 1;
    }
}
