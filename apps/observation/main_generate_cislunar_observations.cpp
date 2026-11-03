/**
 * @file main_generate_cislunar_observations.cpp
 * @brief 读取地月空间轨道、钟差和时隙表，分别生成两个方向的单程理论伪距。
 * @details 主流程依次读取配置与资料、按接收钟时准备插值、求解收发事件，并按原时隙顺序输出。
 * 轨道为 E-BCRS/TDB，位置 m、速度 m/s；钟差为星钟减 TDB，单位 s。
 * 输出 CSV 保留接收钟时、收发 TDB、传播距离、引力时延、两端钟差和伪距。
 * 用法：aster_generate_cislunar_observations -I <配置>；相对路径以配置文件或清单所在目录为基准。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include "io_settings.hpp"
#include "utils_constants.hpp"
#include <sofam.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

using namespace std;
using namespace observation;
using namespace bastion;

int main(int argc, char* argv[])
{
    if (argc == 2 && string(argv[1]) == "--version")
    {
        cout << "aster_generate_cislunar_observations " << ASTER_PROGRAM_VERSION << '\n';
        return 0;
    }

    if (argc == 1 || (argc == 2 && (string(argv[1]) == "--help" || string(argv[1]) == "-h")))
    {
        cout << "Usage: aster_generate_cislunar_observations -I <setting.txt>\n";
        return argc == 1 ? 2 : 0;
    }

    Overwatch tracer(__FILE__, __func__, __LINE__);
    ofstream program_log;
    size_t current_record = 0;

    try
    {
        // =================================== 01 读取配置与输入路径 ===================================
        namespace fs = filesystem;

        if (argc != 3 || string(argv[1]) != "-I") throw runtime_error("Expected -I <setting.txt>.");

        const fs::path config = fs::absolute(argv[2]);
        io::SimpleInput settings;
        string config_name = config.string();
        settings.ReadSimpleInput(tracer, config_name);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

        // 【必选】四项均为文件路径，相对路径从配置目录解释；file_output 为新建的 CSV 文件。
        // file_satellite_list 每行两个带双引号的路径，依次为轨道表和钟差表，从清单目录解释。
        // file_slots 为 LHH 时隙表；file_JPLEPH 为 ASTER 可直接读取的 JPL 二进制历表。
        const array<string, 4> keywords = {"file_satellite_list", "file_slots", "file_JPLEPH", "file_output"};
        array<fs::path, 4> paths;

        for (size_t i = 0; i < keywords.size(); ++i)
        {
            if (!settings.HasKeyword(keywords[i])) throw runtime_error("Missing keyword: " + keywords[i]);

            string value;
            settings.GetString(tracer, keywords[i], value);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit(program_log);

            if (value.empty()) throw runtime_error("Empty path: " + keywords[i]);

            paths[i] = value;

            if (paths[i].is_relative()) paths[i] = config.parent_path() / paths[i];
        }

        if (fs::exists(paths[3])) throw runtime_error("Output already exists: " + paths[3].string());

        array<double, 2> reference_mjd;
        const vector<CislunarSlot> slots = ReadCislunarSlots(paths[1].string(), reference_mjd);
        const int jpl_status = Read_JPLEPH(paths[2].string());
        if (jpl_status < 0) throw runtime_error("Cannot load JPL ephemeris: " + paths[2].string());

        // =================================== 02 读取各星的轨道和钟差 ===================================
        ifstream input_list(paths[0]);
        if (!input_list) throw runtime_error("Cannot read satellite list: " + paths[0].string());

        map<int, CislunarSeries> orbits, clocks;
        string line;

        while (getline(input_list, line))
        {
            istringstream row(line);
            string orbit_name, clock_name;

            if (!(row >> quoted(orbit_name)) || orbit_name[0] == '#') continue;

            if (!(row >> quoted(clock_name))) throw runtime_error("Missing clock path in satellite list: " + line);

            fs::path orbit_path = orbit_name, clock_path = clock_name;

            if (orbit_path.is_relative()) orbit_path = paths[0].parent_path() / orbit_path;

            if (clock_path.is_relative()) clock_path = paths[0].parent_path() / clock_path;

            CislunarSeries orbit = ReadCislunarSeries(orbit_path.string(), reference_mjd, false);
            CislunarSeries clock = ReadCislunarSeries(clock_path.string(), reference_mjd, true);
            const int id = orbit.satellite_id;

            if (id != clock.satellite_id || orbits.count(id))
                throw runtime_error("Mismatched or repeated satellite ID in: " + line);

            orbits.emplace(id, std::move(orbit));
            clocks.emplace(id, std::move(clock));
        }

        // =================================== 03 按接收钟时计算各条观测 ===================================
        vector<size_t> order(slots.size());
        iota(order.begin(), order.end(), 0);
        stable_sort(order.begin(), order.end(), [&slots](size_t a, size_t b)
        {
            return slots[a].clock_tag_day < slots[b].clock_tag_day;
        });
        vector<OneWayRange> observations(slots.size());

        for (size_t index : order)
        {
            current_record = index + 1;
            const CislunarSlot& slot = slots[index];

            if (!orbits.count(slot.transmitter) || !orbits.count(slot.receiver))
                throw runtime_error("Missing satellite table for link " + to_string(slot.transmitter) + " -> " + to_string(slot.receiver));

            // 本条观测的四个窗口按接收钟时固定，随后迭代收发事件。
            const CislunarInterpolant transmit_orbit = PrepareCislunarInterpolant(orbits.at(slot.transmitter), slot.clock_tag_day);
            const CislunarInterpolant receive_orbit = PrepareCislunarInterpolant(orbits.at(slot.receiver), slot.clock_tag_day);
            const CislunarInterpolant transmit_clock = PrepareCislunarInterpolant(clocks.at(slot.transmitter), slot.clock_tag_day);
            const CislunarInterpolant receive_clock = PrepareCislunarInterpolant(clocks.at(slot.receiver), slot.clock_tag_day);

            observations[index] = ComputeCislunarOneWayRange(reference_mjd, slot.clock_tag_day,
                transmit_orbit, receive_orbit, transmit_clock, receive_clock);
        }

        // =================================== 04 按时隙表原顺序写出分项 ===================================
        ofstream output(paths[3]);
        if (!output) throw runtime_error("Cannot create output: " + paths[3].string());

        output << setprecision(17)
               << "# reference_mjd_tdb_integer=" << reference_mjd[0] << '\n'
               << "# reference_mjd_tdb_fraction=" << reference_mjd[1] << '\n'
               << "# jpl_de=" << JPLEPH_R1.r1.numde << '\n'
               << "# gravity_bodies=1,2,4,5,6,7,8,9,10,11; device_delay=0; antenna_offset=0; noise=0\n"
               << "record,type,tx,tx_antenna,rx,rx_antenna,clock_tag_s,tx_tdb_s,rx_tdb_s,geometric_m,gravity_s,gravity_m,"
               << "propagation_m,tx_physical_s,tx_relativistic_s,rx_physical_s,rx_relativistic_s,clock_m,pseudorange_m,"
               << "rx_iterations,tx_iterations,rx_update_s,tx_update_s\n";

        for (size_t i = 0; i < slots.size(); ++i)
        {
            const CislunarSlot& slot = slots[i];
            const OneWayRange& result = observations[i];
            output << i + 1 << ',' << slot.type << ',' << slot.transmitter << ',' << slot.transmit_antenna << ','
                   << slot.receiver << ',' << slot.receive_antenna << ',' << slot.clock_tag_day * DAYSEC << ','
                   << result.transmit_day * DAYSEC << ',' << result.receive_day * DAYSEC << ','
                   << result.geometric_m << ',' << result.gravity_m / SPEED_OF_LIGHT_MPS << ',' << result.gravity_m << ','
                   << result.propagation_m << ',' << result.transmit_physical_s << ',' << result.transmit_relativistic_s << ','
                   << result.receive_physical_s << ',' << result.receive_relativistic_s << ',' << result.clock_m << ','
                   << result.pseudorange_m << ',' << result.receive_iterations << ',' << result.transmit_iterations << ','
                   << result.receive_update_s << ',' << result.transmit_update_s << '\n';
        }

        output.close();
        if (!output) throw runtime_error("Failed while writing: " + paths[3].string());

        cout << "Generated " << observations.size() << " one-way observations with DE" << JPLEPH_R1.r1.numde
             << ".\nOutput: " << paths[3] << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const exception& error)
    {
        cerr << "Cislunar observation error";

        if (current_record) cerr << " at slot " << current_record;

        cerr << ": " << error.what() << '\n';
        return 1;
    }
}
