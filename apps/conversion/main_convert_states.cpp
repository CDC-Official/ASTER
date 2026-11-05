/**
 * @file main_convert_states.cpp
 * @brief 批量转换 Cartesian、Keplerian 和 NSG2 三种地球轨道状态表示。
 * @details 每行输入 UTC MJD 和六分量状态，输出保留历元及参考系。
 * --units 仅支持 deg-m（默认）和 deg-km，输入输出共用一套单位。
 * Cartesian 的位置、速度为 m、m/s 或 km、km/s；Keplerian 的 a 为长度，e 无量纲，其余角度为 deg。
 * NSG2 顺序为 a、xi、eta、h、k、lambda，a 为长度，lambda 为 deg，其余无量纲。
 * 使用 constants::Scale 的地球 GM 和数值尺度，六个方向均调用现有轨道类。
 */
#include "orbit_statedef.hpp"
#include "math_backend.hpp"
#include "utils_constants.hpp"
#include <array>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace bastion;
using constants::Scale::UNIT_LENGTH_M;
using constants::Scale::UNIT_TIME_S;

int main(int argc, char* argv[])
{
    Overwatch tracer(__FILE__, __func__, __LINE__);

    // 版本与帮助用于脚本检查入口，无需读取输入资料。
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_convert_states " << ASTER_PROGRAM_VERSION << '\n'
                  << aster::math_backend::BuildDescription() << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }

    const bool help_requested = argc == 2 &&
        (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h");

    if (argc == 1 || help_requested)
    {
        std::cout << "Usage: aster_convert_states <from> <to> <input> <output> [--units deg-m|deg-km]\n"
                  << "States: cart, kep, nsg2 (six distinct directions).\n"
                  << "Each row: MJD_UTC followed by six state components.\n"
                  << "cart: x y z vx vy vz; kep: a e i OMEGA omega M.\n"
                  << "nsg2: a xi eta h k lambda; xi, eta, h, k are dimensionless.\n"
                  << "Units: deg-m (default) or deg-km; velocities are m/s or km/s.\n"
                  << "Earth elliptical orbits; GM uses constants::Scale::UNIT_GM_M3_S2.\n"
                  << "Blank lines and # comments are accepted. Existing output is overwritten.\n";
        tracer.CheckOut(__func__);
        return help_requested ? 0 : 2;
    }

    try
    {
        // =================================== 01 确定状态表示与单位 ===================================
        // 单位选项放在命令末尾，输入与输出采用同一套度、长度单位。
        int positional_argc = argc;
        std::string units = "deg-m";

        if (argc >= 3 && std::string(argv[argc - 2]) == "--units")
        {
            units = argv[argc - 1];
            positional_argc -= 2;
        }

        if (units != "deg-m" && units != "deg-km")
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Units must be deg-m or deg-km.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        // 类接口的归一化仍以同一物理长度为基准，在这里选择米制或千米制的数值尺度。
        const double length_unit = units == "deg-m" ? UNIT_LENGTH_M : UNIT_LENGTH_M / 1000.0;
        const std::string length_symbol = units == "deg-m" ? "m" : "km";

        if (positional_argc != 5)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected from, to, input, output and optional --units deg-m|deg-km.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        const std::string source = argv[1], target = argv[2];
        const bool valid_source = source == "cart" || source == "kep" || source == "nsg2";
        const bool valid_target = target == "cart" || target == "kep" || target == "nsg2";

        if (!valid_source || !valid_target || source == target)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected two distinct states: cart, kep or nsg2.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        // =================================== 02 逐行读取并转换 ===================================
        std::ifstream input_file(argv[3]);

        if (!input_file)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open input: " + std::string(argv[3]), EventType::ERROR);
            tracer.EmergencyExit();
        }

        // 七列依次保存历元和六维状态；转换完全部记录后再写输出文件。
        std::vector<std::array<double, 7>> output;
        std::string line;
        size_t line_number = 0;

        while (std::getline(input_file, line))
        {
            ++line_number;
            const size_t comment = line.find('#');

            if (comment != std::string::npos) line.resize(comment);

            std::istringstream row(line);
            row >> std::ws;

            if (row.eof()) continue;

            // 每条记录应恰有七个有限数值，时间和坐标的单位由入口约定。
            std::array<double, 7> record{};
            bool valid_record = true;

            for (double& value : record)
            {
                if (!(row >> value) || !std::isfinite(value)) valid_record = false;
            }

            row >> std::ws;

            if (!valid_record || !row.eof())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected seven finite numbers at input line " + std::to_string(line_number), EventType::ERROR);
                tracer.EmergencyExit();
            }

            // 三种表示共用现有轨道对象，转换前统一进入地球二体问题的归一化尺度。
            orbit::Cartesian cartesian;
            orbit::Keplerian keplerian;
            orbit::NonSingular2 nsg2;

            if (source == "cart")
            {
                cartesian.Set(record[1], record[2], record[3], record[4], record[5], record[6], false);
                cartesian.Normalize(tracer, length_unit, UNIT_TIME_S);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                // 当前转换采用椭圆轨道公式，归一化后二体引力参数为 1。
                double state[6];
                cartesian.ToArray(state);
                const double radius = std::hypot(state[0], state[1], state[2]);
                const double speed_squared = state[3] * state[3] + state[4] * state[4] + state[5] * state[5];

                if (!(radius > 0.0) || speed_squared >= 2.0 / radius)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected an Earth elliptical Cartesian state at input line " + std::to_string(line_number), EventType::ERROR);
                    tracer.EmergencyExit();
                }

                // 直接调用目标表示的现有转换函数。
                if (target == "kep")
                {
                    cartesian.ToKeplerian(tracer, keplerian);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
                }
                else
                {
                    cartesian.ToNonSingular2(tracer, nsg2);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
                }
            }
            else if (source == "kep")
            {
                // 六根数依次为 a、e、i、升交点赤经、近地点幅角、平近点角。
                if (!(record[1] > 0.0) || record[2] < 0.0 || record[2] >= 1.0 || record[3] < 0.0 || record[3] > 180.0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected a>0, 0<=e<1 and 0<=i<=180 deg at input line " + std::to_string(line_number), EventType::ERROR);
                    tracer.EmergencyExit();
                }

                keplerian.Set(record[1], record[2], record[3], record[4], record[5], record[6], false);
                keplerian.Normalize(tracer, length_unit);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                if (target == "cart")
                {
                    keplerian.ToCartesian(tracer, cartesian);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
                }
                else
                {
                    keplerian.ToNonSingular2(tracer, nsg2);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
                }
            }
            else
            {
                // xi²+eta²=e²，h²+k²=sin²(i/2)；只有 a 和 lambda 带长度、角度单位。
                if (!(record[1] > 0.0) || std::hypot(record[2], record[3]) >= 1.0 || std::hypot(record[4], record[5]) > 1.0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected a>0, xi^2+eta^2<1 and h^2+k^2<=1 at input line " + std::to_string(line_number), EventType::ERROR);
                    tracer.EmergencyExit();
                }

                nsg2.Set(tracer, record[1], record[2], record[3], record[4], record[5], record[6], false);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                // Normalize 同时同步 NSG2 的近点角等导出量。
                nsg2.Normalize(tracer, length_unit);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                if (target == "cart")
                {
                    nsg2.ToCartesian(tracer, cartesian);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
                }
                else
                {
                    nsg2.ToKeplerian(tracer, keplerian);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
                }
            }

            // 输出统一恢复为所选长度单位和度；Cartesian 速度同时恢复为每秒。
            double converted[6];

            if (target == "cart")
            {
                cartesian.Restore(tracer, length_unit, UNIT_TIME_S);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                cartesian.ToArray(converted);
            }
            else if (target == "kep")
            {
                keplerian.Restore(tracer, length_unit);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                keplerian.ToArray(converted);
            }
            else
            {
                nsg2.Restore(tracer, length_unit);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                nsg2.ToArray(converted);
            }

            // 记录保留输入历元，转换结果按目标表示的六分量顺序保存。
            for (int axis = 0; axis < 6; ++axis)
            {
                if (!std::isfinite(converted[axis]))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Nonfinite converted state at input line " + std::to_string(line_number), EventType::ERROR);
                    tracer.EmergencyExit();
                }

                record[axis + 1] = converted[axis];
            }

            output.push_back(record);
        }

        if (input_file.bad() || output.empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Input is unreadable or contains no state records: " + std::string(argv[3]), EventType::ERROR);
            tracer.EmergencyExit();
        }

        input_file.close();

        // =================================== 03 写出坐标 ===================================
        std::ofstream output_file(argv[4]);

        if (!output_file)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open output: " + std::string(argv[4]), EventType::ERROR);
            tracer.EmergencyExit();
        }

        // 表头说明各分量单位和地球 GM，数据用 17 位有效数字输出。
        output_file << std::setprecision(17);

        if (target == "cart")
        {
            output_file << "# MJD_UTC x_" << length_symbol << " y_" << length_symbol << " z_" << length_symbol
                        << " vx_" << length_symbol << "_s vy_" << length_symbol << "_s vz_" << length_symbol << "_s\n";
        }
        else if (target == "kep")
        {
            output_file << "# MJD_UTC a_" << length_symbol << " e i_deg OMEGA_deg omega_deg M_deg\n";
        }
        else
        {
            output_file << "# MJD_UTC a_" << length_symbol << " xi eta h k lambda_deg\n";
        }

        output_file << "# state=" << target << " units=" << units
                    << " GM_m3_s2=" << constants::Scale::UNIT_GM_M3_S2 << '\n';

        for (const std::array<double, 7>& record : output)
        {
            output_file << record[0];

            for (int axis = 1; axis < 7; ++axis) output_file << ' ' << record[axis];

            output_file << '\n';
        }

        output_file.close();

        if (!output_file)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing output: " + std::string(argv[4]), EventType::ERROR);
            tracer.EmergencyExit();
        }

        std::cerr << "Converted " << output.size() << " records to " << argv[4] << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const std::exception& error)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), EventType::ERROR);
        tracer.EmergencyExit();
    }

    return 1;
}
