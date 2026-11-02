/**
 * @file main_convert_keplerian_to_cartesian.cpp
 * @brief 已知某一时刻的轨道六根数，计算同一时刻的位置和速度。
 *
 * 输入每行：MJD_UTC a[km] e i[deg] OMEGA[deg] omega[deg] M[deg]。
 * M 为平近点角；本入口要求椭圆轨道（a > 0、0 <= e < 1）。
 * 输出每行：MJD_UTC x[m] y[m] z[m] vx[m/s] vy[m/s] vz[m/s]。
 * 各条记录独立转换，保留输入历元和参考系。
 *
 * 运行示例（在项目目录下）：
 *   build/managed-release/bin/aster_convert_keplerian_to_cartesian kep2cart/test_kep.txt build/cart_result.txt
 * CLion 中可选择共享配置“Keplerian to Cartesian - example”。更换计算资料时，
 * 修改运行参数中的输入、输出路径即可。输出文件若已存在会被覆盖，请指定合适路径。
 *
 * @par 本次整理维护记录
 * 2026-10-20：迁入 apps/conversion，统一 CMake 目标及 CLion 运行配置，移除私人路径。
 * 按读取、归一化、根数转 Cartesian、恢复单位和写出的顺序展开 main；增加逐行输入检查及中文功能注释。
 * 原算例与多行样例的输出已与整理前逐字节对比一致；转换公式、输出列和单位未改变。
 */
#include "math_backend.hpp"
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
// 一条输入记录：历元及对应的六根数。
struct KeplerianRecord {
    double mjd_utc; ///< 这组根数对应的 UTC 历元，单位为日（MJD）。
    std::array<double, 6> km_deg; ///< 依次为 a、e、i、OMEGA、omega、M，单位见文件头。
    size_t source_line; ///< 保留原始行号，以定位转换结果异常的记录。
};

// 一条输出记录：原历元及地心位置、速度，单位分别为 m、m/s。
struct CartesianRecord {
    double mjd_utc;
    std::array<double, 6> meters;
};

} // namespace

int main(int argc, const char* argv[]) {
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_convert_keplerian_to_cartesian " << ASTER_PROGRAM_VERSION << '\n';
        std::cout << aster::math_backend::BuildDescription() << '\n';
        return 0;
    }

    // 一、读取运行参数：输入根数文件、输出状态文件。
    const bool help_requested = argc == 2 &&
        (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h");
    if (help_requested || argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input.txt> <output.txt>\n"
                  << "Input: MJD_UTC a[km] e i[deg] OMEGA[deg] omega[deg] M[deg]\n"
                  << "Output: MJD_UTC x[m] y[m] z[m] vx[m/s] vy[m/s] vz[m/s]\n"
                  << "Existing output files are overwritten.\n";
        return help_requested ? 0 : 1;
    }
    // 记录库函数的错误信息。
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    try {
        const std::string input_path = argv[1];
        const std::string output_path = argv[2];
        // 二、按行读取各历元的根数，半长轴为 km，角度为 deg。
        std::ifstream input_file(input_path);
        if (!input_file) throw std::runtime_error("Cannot open input: " + input_path);
        std::vector<KeplerianRecord> input;
        std::string line;
        size_t line_number = 0;
        while (std::getline(input_file, line)) {
            ++line_number;
            std::istringstream row(line);
            row >> std::ws;

            if (row.eof()) continue;

            KeplerianRecord record{};
            // 每行读入历元和六个有限数值，行末只允许空白。
            bool valid_record = static_cast<bool>(row >> record.mjd_utc) && std::isfinite(record.mjd_utc);

            for (double& value : record.km_deg)
            {
                if (!(row >> value) || !std::isfinite(value)) valid_record = false;
            }

            row >> std::ws;

            if (!valid_record || !row.eof())
                throw std::runtime_error("Expected seven finite numbers at line " + std::to_string(line_number));

            // 现有转换公式仅用于椭圆轨道；在进入归一化和求解前检查 a、e。
            if (!(record.km_deg[0] > 0.0 && record.km_deg[1] >= 0.0 && record.km_deg[1] < 1.0))
                throw std::runtime_error("Expected an elliptical orbit (a > 0, 0 <= e < 1) at line " +
                                         std::to_string(line_number));
            record.source_line = line_number;
            input.push_back(record);
        }
        if (input_file.bad()) throw std::runtime_error("Failed while reading: " + input_path);
        if (input.empty()) throw std::runtime_error("Input contains no records: " + input_path);

        // 三、逐条将六根数转换为位置、速度。
        std::vector<CartesianRecord> output;
        output.reserve(input.size());
        for (const auto& record : input) {
            // 当前记录的六根数数组：[a, e, i, OMEGA, omega, M]。
            const auto& kepler = record.km_deg;
            // false 表示输入根数尚未归一化。
            orbit::Keplerian elements(kepler[0], kepler[1], kepler[2], kepler[3], kepler[4], kepler[5], false);
            orbit::Cartesian state;
            // 根数归一化：半长轴由 km 除以 ae/1000（km），角度由 deg 转为 rad。
            elements.Normalize(tracer, ae / 1000.0);
            if (tracer.Code() != ERR_NOERROR) throw std::runtime_error("Keplerian normalization failed");

            // 转换为 Cartesian 位置、速度，结果仍为归一化量。
            elements.ToCartesian(tracer, state);
            if (tracer.Code() != ERR_NOERROR) throw std::runtime_error("Keplerian conversion failed");

            // 恢复单位：位置乘 ae 得到 m，速度乘 ae/UNITTIME 得到 m/s。
            state.Restore(tracer, ae, UNITTIME);
            if (tracer.Code() != ERR_NOERROR) throw std::runtime_error("Cartesian unit restoration failed");
            std::array<double, 6> values;
            state.ToArray(values.data());
            for (double value : values)
                if (!std::isfinite(value))
                    throw std::runtime_error("Conversion produced a non-finite Cartesian state at line " +
                                             std::to_string(record.source_line));

            // 保存原历元及转换后的位置、速度。
            output.push_back({record.mjd_utc, values});
        }

        // 四、按输入顺序输出：MJD_UTC、位置（m）、速度（m/s）。
        // 全部转换完成后再打开输出文件；同名文件会被覆盖，写入失败可能留下部分结果。
        std::ofstream output_file(output_path);
        if (!output_file) throw std::runtime_error("Cannot open output: " + output_path);
        for (const auto& record : output) {
            writeDouble(output_file, record.mjd_utc, 22, 16);
            for (size_t i = 0; i < record.meters.size(); ++i)
                writeDouble(output_file, record.meters[i], 22, 16, i + 1 == record.meters.size());
        }
        output_file.close();
        if (!output_file) throw std::runtime_error("Failed while writing: " + output_path);
        std::cerr << "Converted " << output.size() << " records to " << output_path << '\n';
        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const std::exception& error)
    {
        // 库函数错误已由 tracer 提示，失败时补充调用位置；文件异常直接说明路径。
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
