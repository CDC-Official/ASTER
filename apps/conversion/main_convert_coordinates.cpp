/**
 * @file main_convert_coordinates.cpp
 * @brief 在大地经纬高、ITRS 和 GCRS 之间批量转换位置坐标。
 * @details 每行输入 UTC MJD 和三个位置分量，按原顺序写出相同时刻的目标坐标。
 * 大地坐标按经度、大地纬度、椭球高排列，直角坐标为 X、Y、Z；输入输出共用所选单位。
 * --units 仅支持 deg-m（默认）和 deg-km，经纬度为 deg，长度分别为 m 和 km。
 * 椭球采用 WGS84；涉及 GCRS 时读取 EOP，调用公共位置转换接口。
 * 用法：aster_convert_coordinates <源坐标> <目标坐标> <输入> <输出> [EOP文件] [--units deg-m|deg-km]。
 * 坐标名称为 geodetic、itrs、gcrs；文件读写及批量流程直接在主程序中组织。
 */
#include "frame_definition.hpp"
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

int main(int argc, char* argv[])
{
    Overwatch tracer(__FILE__, __func__, __LINE__);

    // 版本与帮助用于脚本检查入口，无需读取输入资料。
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_convert_coordinates " << ASTER_PROGRAM_VERSION << '\n'
                  << aster::math_backend::BuildDescription() << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }

    const bool help_requested = argc == 2 &&
        (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h");

    if (argc == 1 || help_requested)
    {
        std::cout << "Usage: aster_convert_coordinates <from> <to> <input> <output> [EOP] [--units deg-m|deg-km]\n"
                  << "Frames: geodetic, itrs, gcrs (six distinct directions).\n"
                  << "Geodetic row: MJD_UTC longitude_deg latitude_deg height (WGS84).\n"
                  << "ITRS/GCRS row: MJD_UTC X Y Z. Positions only.\n"
                  << "Units: deg-m (default) or deg-km; one choice applies to input and output.\n"
                  << "EOP is required for every direction involving GCRS.\n"
                  << "Blank lines and # comments are accepted. Existing output is overwritten.\n";
        tracer.CheckOut(__func__);
        return help_requested ? 0 : 2;
    }

    try
    {
        // =================================== 01 确定转换方向与单位 ===================================
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

        if (positional_argc != 5 && positional_argc != 6)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected from, to, input, output and optional EOP path.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        const std::string source = argv[1], target = argv[2];
        const bool valid_source = source == "geodetic" || source == "itrs" || source == "gcrs";
        const bool valid_target = target == "geodetic" || target == "itrs" || target == "gcrs";

        if (!valid_source || !valid_target || source == target)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected two distinct frames: geodetic, itrs or gcrs.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        // 纯椭球转换只用 WGS84；地固与天球之间的旋转还需要观测时刻及 EOP。
        const bool needs_eop = source == "gcrs" || target == "gcrs";

        if (positional_argc != (needs_eop ? 6 : 5))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Supply an EOP file exactly when converting to or from GCRS.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        // 复用现有 EOP 读取器，逐条转换时再按 UTC 插值。
        EOP eop_table;

        if (needs_eop)
        {
            std::string eop_path = argv[5];
            eop_table.Read_EOP(tracer, eop_path);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
        }

        // =================================== 02 逐行读取并转换 ===================================
        std::ifstream input_file(argv[3]);

        if (!input_file)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open input: " + std::string(argv[3]), EventType::ERROR);
            tracer.EmergencyExit();
        }

        // 四列依次保存历元和位置；转换完全部记录后再写输出文件。
        std::vector<std::array<double, 4>> output;
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

            // 每条记录应恰有四个有限数值，时间和坐标的单位由入口约定。
            std::array<double, 4> record{};
            bool valid_record = true;

            for (double& value : record)
            {
                if (!(row >> value) || !std::isfinite(value)) valid_record = false;
            }

            row >> std::ws;

            if (!valid_record || !row.eof())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected four finite numbers at input line " + std::to_string(line_number), EventType::ERROR);
                tracer.EmergencyExit();
            }

            // 按本条记录的 UTC 时刻提取 EOP，供地固与天球坐标之间的旋转使用。
            aEOP eop;

            if (needs_eop)
            {
                eop_table.Get_EOP(tracer, record[0], eop, 1);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
            }

            // 大地坐标对象按高度、经度、纬度构造，文件则按经度、纬度、高度排列。
            double converted[3];

            if (source == "geodetic")
            {
                frame::Geodetic geodetic(record[3], record[1], record[2], false);
                geodetic.Normalize(tracer, length_unit);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                if (target == "itrs")
                {
                    frame::ITRS itrs;
                    geodetic.ToITRS(tracer, itrs);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    // 科学接口输出归一化位置，文件恢复为所选的 m 或 km。
                    itrs.Restore(tracer, length_unit);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    itrs.ToArray(converted);
                }
                else
                {
                    frame::GCRS gcrs;
                    geodetic.ToGCRS(tracer, record[0], eop, gcrs);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    gcrs.Restore(tracer, length_unit);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    gcrs.ToArray(converted);
                }
            }
            else if (source == "itrs")
            {
                frame::ITRS itrs(record[1], record[2], record[3], false);

                if (target == "gcrs")
                {
                    // 位置旋转保持所选长度单位，直接得到地心天球位置。
                    frame::GCRS gcrs;
                    itrs.ToGCRS(tracer, record[0], eop, gcrs);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    gcrs.ToArray(converted);
                }
                else
                {
                    itrs.Normalize(tracer, length_unit);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    frame::Geodetic geodetic;
                    itrs.ToGeodetic(tracer, geodetic);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    // 先恢复所选长度单位和 deg，再按经度、纬度、高度顺序保存。
                    geodetic.Restore(tracer, length_unit);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    double values[3];
                    geodetic.ToArray(values);
                    converted[0] = values[1];
                    converted[1] = values[2];
                    converted[2] = values[0];
                }
            }
            else
            {
                frame::GCRS gcrs(record[1], record[2], record[3], false);

                if (target == "itrs")
                {
                    frame::ITRS itrs;
                    gcrs.ToITRS(tracer, record[0], eop, itrs);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    itrs.ToArray(converted);
                }
                else
                {
                    // 一步接口在内部依次进行天球到地固、地固到大地坐标的转换。
                    gcrs.Normalize(tracer, length_unit);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    frame::Geodetic geodetic;
                    gcrs.ToGeodetic(tracer, record[0], eop, geodetic);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    geodetic.Restore(tracer, length_unit);
                    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

                    double values[3];
                    geodetic.ToArray(values);
                    converted[0] = values[1];
                    converted[1] = values[2];
                    converted[2] = values[0];
                }
            }

            output.push_back({record[0], converted[0], converted[1], converted[2]});
        }

        if (input_file.bad() || output.empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Input is unreadable or contains no coordinate records: " + std::string(argv[3]), EventType::ERROR);
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

        // 表头明确坐标含义，17 位有效数字保留 double 的数值精度。
        if (target == "geodetic")
        {
            output_file << "# MJD_UTC longitude_deg latitude_deg height_" << length_symbol << '\n';
        }
        else
        {
            output_file << "# MJD_UTC X_" << length_symbol << " Y_" << length_symbol << " Z_" << length_symbol << '\n';
        }

        output_file << "# frame=" << target;

        if (source == "geodetic" || target == "geodetic") output_file << " ellipsoid=WGS84";

        if (needs_eop) output_file << " rotation=IAU2000B";

        output_file << '\n' << std::setprecision(17);

        for (const std::array<double, 4>& record : output)
        {
            output_file << record[0] << ' ' << record[1] << ' ' << record[2] << ' ' << record[3] << '\n';
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
