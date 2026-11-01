/**
 * @file main_convert_cartesian_to_keplerian.cpp
 * @brief 将同一历元的位置、速度转换为地心轨道六根数。
 *
 * 输入每行：MJD_UTC x y z vx vy vz，支持 m、m/s 或 km、km/s。
 * 输出每行：MJD_UTC a[km] e i[deg] OMEGA[deg] omega[deg] M[deg]，M 为平近点角。
 * 本入口使用现有椭圆轨道转换公式；各记录独立处理，不改变参考系和历元。
 *
 * 用法：aster_convert_cartesian_to_keplerian <input.txt> <output.txt> [m|km]
 * 省略单位参数时，保留旧程序的交互选择；CLion 示例显式使用 m，避免等待输入。
 * 输出若已存在会覆盖，建议写入 build/ 或单独的结果目录。
 *
 * @par 本次整理维护记录
 * 2026-10-20：恢复 Cartesian 转根数的独立入口、CMake 目标及 CLion 配置，补齐与反向转换配套的程序。
 * 显式指定输入、输出路径，保留 m／km 与交互式单位选择；增加逐行格式及有限数值检查。
 * 已对比历史算例并验证双向转换，结果在浮点精度范围内吻合；详见 doc/cart_to_kep_restoration_20261020.md。
 */
#include "utils_constants.hpp"
#include "io_common.hpp"
#include "orbit_statedef.hpp"
#include "utils_overwatch.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
// 历元与对应的状态一起保存；values 保持输入所选的 m 或 km 单位。
struct CartesianRecord {
    double mjd_utc;
    std::array<double, 6> values;
};

// 转换后的六根数，依次为 a、e、i、OMEGA、omega、M，采用 km、deg。
struct KeplerianRecord {
    double mjd_utc;
    std::array<double, 6> km_deg;
};
}

int main(int argc, const char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_convert_cartesian_to_keplerian " << ASTER_PROGRAM_VERSION << '\n';
        return 0;
    }

    // 一、确定输入、输出文件和状态量单位；不再使用写在源码中的私人数据路径。
    const bool help_requested = argc == 2 &&
        (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h");
    if (help_requested || (argc != 3 && argc != 4)) {
        std::cerr << "Usage: " << argv[0] << " <input.txt> <output.txt> [m|km]\n"
                  << "Input: MJD_UTC x y z vx vy vz (m, m/s or km, km/s)\n"
                  << "Output: MJD_UTC a[km] e i[deg] OMEGA[deg] omega[deg] M[deg]\n"
                  << "Omitting the unit prompts for 1 (m) or 2 (km). Existing output is overwritten.\n";
        return help_requested ? 0 : 1;
    }
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    try {
        const std::string input_path = argv[1];
        const std::string output_path = argv[2];
        std::string unit;
        if (argc == 4) {
            unit = argv[3];
        } else {
            std::cerr << "Please select the input unit:\n1. m, m/s\n2. km, km/s\n";
            int choice = 0;
            if (!(std::cin >> choice) || (choice != 1 && choice != 2))
                throw std::runtime_error("Input unit must be 1 (m) or 2 (km)");
            unit = choice == 1 ? "m" : "km";
        }
        if (unit != "m" && unit != "km")
            throw std::runtime_error("Input unit must be m or km");
        const double input_length = unit == "m" ? ae : ae / 1000.0;

        // 二、按行读取历元及六维状态。允许空行，每个非空行须为七个有限数值。
        std::ifstream input_file(input_path);
        if (!input_file) throw std::runtime_error("Cannot open input: " + input_path);
        std::vector<CartesianRecord> input;
        std::string line;
        size_t line_number = 0;
        while (std::getline(input_file, line)) {
            ++line_number;
            std::istringstream row(line);
            row >> std::ws;

            if (row.eof()) continue;

            CartesianRecord record{};
            // 每行读入历元和六个有限数值，行末只允许空白。
            bool valid_record = static_cast<bool>(row >> record.mjd_utc) && std::isfinite(record.mjd_utc);

            for (double& value : record.values)
            {
                if (!(row >> value) || !std::isfinite(value)) valid_record = false;
            }

            row >> std::ws;

            if (!valid_record || !row.eof())
                throw std::runtime_error("Expected seven finite numbers at line " + std::to_string(line_number));

            input.push_back(record);
        }
        if (input_file.bad()) throw std::runtime_error("Failed while reading: " + input_path);
        if (input.empty()) throw std::runtime_error("Input contains no records: " + input_path);

        // 三、逐历元归一化并转换；根数恢复为 km、deg，输入历元保持不变。
        std::vector<KeplerianRecord> output;
        output.reserve(input.size());
        for (const auto& record : input) {
            const auto& x = record.values;
            orbit::Cartesian state(x[0], x[1], x[2], x[3], x[4], x[5], false);
            orbit::Keplerian elements;
            // 位置除以所选长度尺度，速度除以该尺度/UNITTIME，转为归一化状态。
            state.Normalize(tracer, input_length, UNITTIME);
            if (tracer.Code() != ERR_NOERROR) throw std::runtime_error("Cartesian normalization failed");

            // 转换为 Keplerian 六根数，结果仍为归一化量和 rad。
            state.ToKeplerian(tracer, elements);
            if (tracer.Code() != ERR_NOERROR) throw std::runtime_error("Cartesian conversion failed");

            // 半长轴恢复为 km，轨道角度由 rad 转为 deg。
            elements.Restore(tracer, ae / 1000.0);
            if (tracer.Code() != ERR_NOERROR) throw std::runtime_error("Keplerian unit restoration failed");
            std::array<double, 6> values;
            elements.ToArray(values.data());
            for (double value : values)
                if (!std::isfinite(value))
                    throw std::runtime_error("Conversion produced non-finite elements; check for a degenerate or unsupported orbit");
            output.push_back({record.mjd_utc, values});
        }

        // 四、全部转换成功后再打开输出文件，避免坏输入覆盖已有结果。
        std::ofstream output_file(output_path);
        if (!output_file) throw std::runtime_error("Cannot open output: " + output_path);
        for (const auto& record : output) {
            writeDouble(output_file, record.mjd_utc, 22, 16);
            for (size_t i = 0; i < record.km_deg.size(); ++i)
                writeDouble(output_file, record.km_deg[i], 22, 16, i + 1 == record.km_deg.size());
        }
        output_file.close();
        if (!output_file) throw std::runtime_error("Failed while writing: " + output_path);
        std::cerr << "Converted " << output.size() << " records to " << output_path << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const std::exception& error)
    {
        if (tracer.Code() != ERR_NOERROR)
        {
            tracer.ReportFailure();
        }
        else
        {
            std::cerr << "Conversion failed: " << error.what() << '\n';
        }

        return 1;
    }
}
