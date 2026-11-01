/**
 * @file main_convert_osculating_to_mean.cpp
 * @brief 逐历元读取 GCRS 密切 Kepler 根数，转换为同历元均根数并写出七列星历。
 * @details 流程为读取配置、检查路径、读取全部轨道、准备模型、逐条瞬转平、写出星历和日志。
 * 输入输出列为 [MJD_UTC,a,e,i,Omega,omega,M]，a 为 km，四个角为 deg；库内用 ae、UNITTIME 和 rad 归一化。
 * 配置采用严格 SimpleInput，所有相对路径以配置目录为基准。结果独占创建，全部计算成功后才写出。
 */
#include "analyt_averaging.hpp"
#include "io_common.hpp"
#include "io_settings.hpp"
#include "utils_constants.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using bastion::Overwatch;
using orbit::Keplerian;
using orbit::NonSingular2;

/**
 * @brief 执行一次逐历元瞬转平任务。
 * @param[in] argc 命令行参数个数。
 * @param[in] argv 参数数组，支持 -I/-input/--input/--setting 指定配置，-h/--help 显示帮助。
 * @return 成功为 0，任务失败为 1，命令行格式错误为 2。
 * @details 科学文件在本函数中直接读取和写出，模型摘要由 AveragingOptions::Display 显示。
 * 主动错误写入 tracer，跳出当前任务后统一关闭句柄、撤回本轮不完整结果并保留失败日志。
 */
int main(int argc, char** argv)
{
    if (argc == 2 && std::string(argv[1]) == "--version")
    {
        std::cout << "aster_convert_osculating_to_mean " << ASTER_PROGRAM_VERSION << '\n';
        return 0;
    }

    if (argc == 1)
    {
        std::cerr << "Usage: aster_convert_osculating_to_mean -I <setting.txt>\n";
        return 2;
    }

    // =================================== 01 读取命令行 ===================================

    Overwatch tracer(__FILE__, __func__, __LINE__);
    std::string setting_file = "";

    for (int i = 1; i < argc; ++i)
    {
        const std::string arg = argv[i];

        if (arg == "-h" || arg == "--help")
        {
            std::cout << "Usage: aster_convert_osculating_to_mean [-I|-input|--input|--setting <setting.txt>]" << std::endl;
            return 0;
        }

        if (arg != "-I" && arg != "-input" && arg != "--input" && arg != "--setting")
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Unknown argument: " + arg, bastion::ERROR);
            return 2;
        }

        if (i + 1 >= argc || argv[i + 1][0] == '-' || std::string(argv[i + 1]).empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Missing setting path for " + arg, bastion::ERROR);
            return 2;
        }

        setting_file = argv[++i];
    }

    // 文件句柄及结果所有权在任务块外保存，失败时可在同一出口清理。
    FILE* log_file = nullptr;
    FILE* output_file = nullptr;
    fs::path created_output;
    bool output_created = false;

    try
    {
        // 单次任务块：错误从各阶段跳出，块后统一清理；正常结束直接返回成功。
        do
        {
            // =================================== 02 读取文件路径及模型设置 ===================================

            const fs::path setting_path = fs::absolute(setting_file).lexically_normal();
            const fs::path config_dir = setting_path.parent_path();
            std::string config_name = setting_path.string();
            io::SimpleInput input;

            // 严格模式拒绝缺少等号、空键和重复键，# 可用于行尾注释。
            input.ReadSimpleInput(tracer, config_name, true);
            if (tracer.Code() != bastion::ERR_NOERROR) break;

            // 只接受本程序定义的关键词，拼写错误不作为默认配置继续计算。
            const std::set<std::string> keywords{
                "input_file", "output_file", "log_file", "gravity_degree", "iteration_count",
                "hansen_expansion_order", "use_J2", "use_J22", "use_frame", "use_JL", "use_JLM",
                "use_lunar", "use_solar", "use_SRP", "display_records", "SRP_AMR_m2_kg",
                "SRP_AMR_normalized", "drag_AMR_normalized", "reflectivity", "lunar_short_period"
            };

            for (std::size_t index = 1; index <= input.Count(); ++index)
            {
                io::SimpleRecord record;
                input.GetRecord(tracer, index, record);
                if (tracer.Code() != bastion::ERR_NOERROR) break;

                if (keywords.count(record.kw) == 0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "未知 setting 项: " + record.kw + " = " + record.va, bastion::ERROR);
                    break;
                }
            }

            if (tracer.Code() != bastion::ERR_NOERROR) break;

            // 【必选】input_file：七列密切 Kepler 根数文件，不能省略或留空。
            std::string input_name;
            input.GetString(tracer, "input_file", input_name);
            if (tracer.Code() != bastion::ERR_NOERROR) break;

            if (input_name.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "setting 缺少必填项: input_file", bastion::ERROR);
                break;
            }

            fs::path input_path = input_name;

            if (input_path.is_relative()) input_path = (config_dir / input_path).lexically_normal();

            // 【必选】output_file：同历元七列均根数文件，不能留空；已有文件拒绝覆盖。
            std::string output_name;
            input.GetString(tracer, "output_file", output_name);
            if (tracer.Code() != bastion::ERR_NOERROR) break;

            if (output_name.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "setting 缺少必填项: output_file", bastion::ERROR);
                break;
            }

            fs::path output_path = output_name;

            if (output_path.is_relative()) output_path = (config_dir / output_path).lexically_normal();

            // 【可选】log_file：模型、警告和任务状态日志；省略或留空时不创建日志。
            fs::path log_path;

            if (input.HasKeyword("log_file"))
            {
                std::string log_name;
                input.GetString(tracer, "log_file", log_name);
                if (tracer.Code() != bastion::ERR_NOERROR) break;

                if (!log_name.empty())
                {
                    log_path = log_name;

                    if (log_path.is_relative()) log_path = (config_dir / log_path).lexically_normal();
                }
            }

            // 模型参数直接填入库内选项；缺项时保留构造默认值。
            orbit::AveragingOptions model;
            bool display_records = false;
            double srp_amr_norm = 1468568181.9309881;
            double drag_amr_norm = 1468568181.9309881;

            // 【可选】gravity_degree：JL/JLM 最高重力次数，整数 2..10，默认 10；J2/J22 各有独立开关。
            if (input.HasKeyword("gravity_degree"))
            {
                input.GetIntegerStrict(tracer, "gravity_degree", model.gravity_degree);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】iteration_count：固定点迭代次数，整数 1..100，默认 3，不按残差提前终止。
            if (input.HasKeyword("iteration_count"))
            {
                input.GetIntegerStrict(tracer, "iteration_count", model.iteration_count);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】hansen_expansion_order：Hansen 偏心率展开阶数及 q 截断，整数 0..4，默认 4。
            if (input.HasKeyword("hansen_expansion_order"))
            {
                input.GetIntegerStrict(tracer, "hansen_expansion_order", model.hansen_order);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_J2：是否扣除 J2 短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_J2"))
            {
                input.GetBoolean(tracer, "use_J2", model.use_j2);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_J22：是否扣除 J22 短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_J22"))
            {
                input.GetBoolean(tracer, "use_J22", model.use_j22);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_frame：是否扣除岁差章动与 J2 的耦合项，默认开启，需同时启用 J2。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_frame"))
            {
                input.GetBoolean(tracer, "use_frame", model.use_frame);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_JL：是否计算 l>=3、m=0 的带谐短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_JL"))
            {
                input.GetBoolean(tracer, "use_JL", model.use_jl);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_JLM：是否计算 l>=3、m>0 的田谐短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_JLM"))
            {
                input.GetBoolean(tracer, "use_JLM", model.use_jlm);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_lunar：是否扣除月球引力短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_lunar"))
            {
                input.GetBoolean(tracer, "use_lunar", model.use_lunar);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_solar：是否扣除太阳引力短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_solar"))
            {
                input.GetBoolean(tracer, "use_solar", model.use_solar);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】use_SRP：是否扣除无地影光压短周期项，默认开启。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("use_SRP"))
            {
                input.GetBoolean(tracer, "use_SRP", model.use_srp);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】display_records：是否在控制台逐条显示七列结果，默认关闭。布尔值支持 true/false、1/0、on/off。
            if (input.HasKeyword("display_records"))
            {
                input.GetBoolean(tracer, "display_records", display_records);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【二选一可选】SRP_AMR_m2_kg 与 SRP_AMR_normalized：非负光压面积质量比。
            // 两项都省略时使用归一化默认值 1468568181.9309881，约为 0.01 m²/kg。
            if (input.HasKeyword("SRP_AMR_m2_kg") && input.HasKeyword("SRP_AMR_normalized"))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "SRP_AMR_m2_kg 与 SRP_AMR_normalized 只能配置一个", bastion::ERROR);
                break;
            }

            if (input.HasKeyword("SRP_AMR_m2_kg"))
            {
                double srp_amr_m2_kg;
                input.GetDoubleStrict(tracer, "SRP_AMR_m2_kg", srp_amr_m2_kg);
                if (tracer.Code() != bastion::ERR_NOERROR) break;

                srp_amr_norm = srp_amr_m2_kg / NORMALIZED_AM_TO_IAM;
            }
            else if (input.HasKeyword("SRP_AMR_normalized"))
            {
                input.GetDoubleStrict(tracer, "SRP_AMR_normalized", srp_amr_norm);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】drag_AMR_normalized：非负归一化阻力面积质量比，默认 1468568181.9309881；只在日志中留存。
            if (input.HasKeyword("drag_AMR_normalized"))
            {
                input.GetDoubleStrict(tracer, "drag_AMR_normalized", drag_amr_norm);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】reflectivity：无量纲 Cr，范围 0..2，默认 1.44。
            if (input.HasKeyword("reflectivity"))
            {
                input.GetDoubleStrict(tracer, "reflectivity", model.radiation_coefficient);
                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            // 【可选】lunar_short_period：single_fastarg（默认）或 double_fastarg，只在启用月球项时参与计算。
            if (input.HasKeyword("lunar_short_period"))
            {
                std::string lunar_mode;
                input.GetString(tracer, "lunar_short_period", lunar_mode);
                if (tracer.Code() != bastion::ERR_NOERROR) break;

                if (lunar_mode == "single_fastarg")
                {
                    model.lunar_double_average = false;
                }
                else if (lunar_mode == "double_fastarg")
                {
                    model.lunar_double_average = true;
                }
                else
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "lunar_short_period 仅支持 single_fastarg 或 double_fastarg", bastion::ERROR);
                    break;
                }
            }

            // 检查各模型的支持范围，面积质量比在单位转换后再次检查有限性。
            if (model.gravity_degree < 2 || model.gravity_degree > 10)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "当前内置 GGM03C 数据支持 gravity_degree 2 到 10", bastion::ERROR);
                break;
            }

            if (model.iteration_count < 1 || model.iteration_count > 100)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "iteration_count 必须在 1 到 100 之间", bastion::ERROR);
                break;
            }

            if (model.hansen_order < 0 || model.hansen_order > 4)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "hansen_expansion_order 必须在 0 到 4 之间", bastion::ERROR);
                break;
            }

            if (!std::isfinite(srp_amr_norm) || srp_amr_norm < 0 || drag_amr_norm < 0 ||
                model.radiation_coefficient < 0 || model.radiation_coefficient > 2)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "AMR 必须非负，reflectivity 必须在 0 到 2 之间", bastion::ERROR);
                break;
            }

            // 库接口接收物理 m²/kg；原始归一化数值只留作日志显示。
            model.srp_area_mass_m2_kg = srp_amr_norm * NORMALIZED_AM_TO_IAM;

            // =================================== 03 检查路径并创建日志 ===================================

            // 前两项为输入，后两项为产品；检查产品与此前各项的相等、包含及硬链接关系。
            const std::array<fs::path, 4> paths{setting_path, input_path, output_path, log_path};

            for (std::size_t i = 2; i < paths.size(); ++i)
            {
                if (paths[i].empty()) continue;

                if (fs::exists(paths[i]) && !fs::is_regular_file(paths[i]))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "输出路径已存在且不是普通文件: " + paths[i].string(), bastion::ERROR);
                    break;
                }

                const fs::path product = fs::weakly_canonical(paths[i]);

                for (std::size_t j = 0; j < i; ++j)
                {
                    const fs::path other = fs::weakly_canonical(paths[j]);
                    const auto prefix = std::mismatch(product.begin(), product.end(), other.begin(), other.end());

                    // 逐段比较路径，任意一方先结束表示相等或父子关系；equivalent 另外识别硬链接。
                    if (prefix.first == product.end() || prefix.second == other.end() ||
                        (fs::exists(paths[i]) && fs::exists(paths[j]) && fs::equivalent(paths[i], paths[j])))
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__, "输出路径与输入资料或其他产品冲突: " + paths[i].string(), bastion::ERROR);
                        break;
                    }
                }

                if (tracer.Code() != bastion::ERR_NOERROR) break;
            }

            if (tracer.Code() != bastion::ERR_NOERROR) break;

            // 日志先独占创建，后续输入或计算错误可写入这份日志；空路径关闭日志。
            if (!log_path.empty())
            {
                log_file = std::fopen(log_path.string().c_str(), "wx");

                if (!log_file)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "无法创建日志文件: " + log_path.string(), bastion::ERROR);
                    break;
                }
            }

            // =================================== 04 读取并检查七列密切根数 ===================================

            std::ifstream ephemeris(input_path);

            if (!ephemeris)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "无法打开输入星历: " + input_path.string(), bastion::ERROR);
                break;
            }

            // 数组保存原始七列及原始行号，计算阶段据此报告具体记录的警告或错误。
            std::vector<std::array<double, 7>> records;
            std::vector<std::size_t> input_lines;
            std::string line;
            std::size_t line_number = 0;
            double last_epoch = -std::numeric_limits<double>::infinity();

            while (std::getline(ephemeris, line))
            {
                ++line_number;
                const std::size_t comment = line.find('#');

                if (comment != std::string::npos) line.erase(comment);

                if (strip(line).empty()) continue;

                // 一行只读一个 UTC 历元及六根数，缺列和多列都拒绝。
                std::istringstream row(line);
                std::array<double, 7> values{};

                if (!(row >> values[0]))
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "输入第 " + std::to_string(line_number) + " 行历元不是有效数值", bastion::ERROR);
                    break;
                }

                for (std::size_t i = 1; i < values.size(); ++i)
                {
                    if (!(row >> values[i]))
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__, "输入第 " + std::to_string(line_number) + " 行不足 7 列", bastion::ERROR);
                        break;
                    }
                }

                if (tracer.Code() != bastion::ERR_NOERROR) break;

                std::string trailing;

                if (row >> trailing)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "输入第 " + std::to_string(line_number) + " 行多于 7 列", bastion::ERROR);
                    break;
                }

                // 原始单位为 km、deg；确认有限性后按椭圆定义域和非递减历元验收。
                for (double value : values)
                {
                    if (!std::isfinite(value))
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__, "输入包含非有限数值，行 " + std::to_string(line_number), bastion::ERROR);
                        break;
                    }
                }

                if (tracer.Code() != bastion::ERR_NOERROR) break;

                if (values[1] <= 0 || values[2] < 0 || values[2] >= 1 || values[3] < 0 || values[3] >= 180)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "输入根数要求 a>0、0<=e<1、0<=i<180 deg，行 " + std::to_string(line_number), bastion::ERROR);
                    break;
                }

                if (values[0] < last_epoch)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "输入历元必须非递减，行 " + std::to_string(line_number), bastion::ERROR);
                    break;
                }

                last_epoch = values[0];
                records.push_back(values);
                input_lines.push_back(line_number);
            }

            if (tracer.Code() != bastion::ERR_NOERROR) break;

            if (ephemeris.bad() || records.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "读取输入星历失败或没有有效数据行: " + input_path.string(), bastion::ERROR);
                break;
            }

            // =================================== 05 显示任务资料及模型设置 ===================================

            if (log_file)
            {
                // 路径、列定义及记录数由主程序写出，模型部分由已有选项对象显示。
                std::ostringstream header;
                header << "ASTER osculating-to-mean conversion\n"
                    << "Implementation: ASTER analytical averaging\n"
                    << "Setting: " << setting_path << "\nInput: " << input_path << "\nOutput: " << output_path << "\nRecords: " << records.size()
                    << "\nLog: " << log_path << "\nEOP: none (DUT1=0, advisor default)"
                    << "\nEpoch: MJD UTC input; TT used internally by ASTER\n"
                    << "Elements: osculating/mean Keplerian; a in km, e dimensionless, angles in deg; last angle is mean anomaly M\n"
                    << "Frame: GCRS input convention; no input frame transformation\n";

                model.Display(tracer, srp_amr_norm, drag_amr_norm, header);
                if (tracer.Code() != bastion::ERR_NOERROR) break;

                if (std::fputs(header.str().c_str(), log_file) < 0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "写出模型日志失败: " + log_path.string(), bastion::ERROR);
                    break;
                }
            }

            std::cout << "O2M: " << records.size() << " records\nOutput: " << output_path << '\n';

            if (!log_path.empty()) std::cout << "Log: " << log_path << '\n';

            // =================================== 06 逐历元密切根数转均根数 ===================================

            std::vector<std::array<double, 7>> results;
            results.reserve(records.size());

            for (std::size_t index = 0; index < records.size(); ++index)
            {
                const std::array<double, 7>& values = records[index];
                const std::string location = "，输入行 " + std::to_string(input_lines[index]);
                const int previous_warnings = tracer.WarningCount();

                // 输入采用 km、deg；Normalize 将半长轴换为 a/ae、角度换为 rad。
                Keplerian osculating(values[1], values[2], values[3], values[4], values[5], values[6], false);
                osculating.Normalize(tracer, ae / 1000.);
                if (tracer.Code() != bastion::ERR_NOERROR)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "根数归一化失败" + location, bastion::ERROR);
                    break;
                }

                // 先转换为 NSG2，再在同一 UTC 标签下扣除所选模型的解析短周期项。
                NonSingular2 osculating_nsg2, mean_nsg2;
                osculating.ToNonSingular2(tracer, osculating_nsg2);
                if (tracer.Code() != bastion::ERR_NOERROR)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Keplerian 到 NSG2 转换失败" + location, bastion::ERROR);
                    break;
                }

                orbit::OsculatingToMean(tracer, values[0], osculating_nsg2, model, mean_nsg2);
                if (tracer.Code() != bastion::ERR_NOERROR)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "O2M 失败" + location, bastion::ERROR);
                    break;
                }

                // 将均 NSG2 转回 Kepler，随后恢复与输入相同的 km、deg。
                Keplerian mean;
                mean_nsg2.ToKeplerian(tracer, mean);
                if (tracer.Code() != bastion::ERR_NOERROR)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "NSG2 到 Keplerian 转换失败" + location, bastion::ERROR);
                    break;
                }

                mean.Restore(tracer, ae / 1000.);
                if (tracer.Code() != bastion::ERR_NOERROR)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "均根数单位恢复失败" + location, bastion::ERROR);
                    break;
                }

                std::array<double, 6> elements;
                mean.ToArray(elements.data());
                std::array<double, 7> row{};
                row[0] = values[0];

                // 保留原 UTC 历元，检查并填入 [a,e,i,Omega,omega,M] 六列。
                for (std::size_t i = 0; i < elements.size(); ++i)
                {
                    if (!std::isfinite(elements[i]))
                    {
                        tracer.LogEvent(__FILE__, __func__, __LINE__, "O2M 产生非有限结果" + location, bastion::ERROR);
                        break;
                    }

                    row[i + 1] = elements[i];
                }

                if (tracer.Code() != bastion::ERR_NOERROR) break;

                if (row[1] <= 0 || row[2] < 0 || row[2] >= 1)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "O2M 结果超出椭圆轨道范围" + location, bastion::ERROR);
                    break;
                }

                results.push_back(row);

                // 【display_records】仅控制终端显示；七列顺序与文件一致，采用 17 位有效数字。
                if (display_records)
                {
                    std::cout << std::setprecision(17) << row[0];

                    for (std::size_t i = 1; i < row.size(); ++i) std::cout << ' ' << row[i];

                    std::cout << '\n';
                }

                // 将本条记录新增的模型警告关联到原始输入行，事件正文在日志末尾集中写出。
                if (log_file && tracer.WarningCount() != previous_warnings)
                {
                    std::fprintf(log_file, "Input line %zu: %d new warning(s)\n", input_lines[index], tracer.WarningCount() - previous_warnings);
                }
            }

            if (tracer.Code() != bastion::ERR_NOERROR) break;

            // =================================== 07 写出七列均根数 ===================================

            created_output = output_path;
            output_file = std::fopen(output_path.string().c_str(), "wx");

            if (!output_file)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "无法创建输出文件: " + output_path.string(), bastion::ERROR);
                break;
            }

            // 只有本轮独占创建成功的结果，才允许在后续失败时撤回。
            output_created = true;

            for (const std::array<double, 7>& row : results)
            {
                // 无表头的七列纯数值，MJD UTC、km、deg 均以 17 位有效数字写出。
                if (std::fprintf(output_file, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    row[0], row[1], row[2], row[3], row[4], row[5], row[6]) < 0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "写出结果失败: " + output_path.string(), bastion::ERROR);
                    break;
                }
            }

            if (tracer.Code() != bastion::ERR_NOERROR) break;

            // fclose 已处理句柄，先置空再检查，避免失败清理时重复关闭。
            const int output_close_status = std::fclose(output_file);
            output_file = nullptr;

            if (output_close_status != 0)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "关闭结果文件失败: " + output_path.string(), bastion::ERROR);
                break;
            }

            // =================================== 08 保存事件和完成状态 ===================================

            if (log_file)
            {
                std::ostringstream events;
                tracer.ReportEvents(events);

                if (std::fputs("\nPer-orbit model warnings and events:\n", log_file) < 0 ||
                    std::fputs(events.str().c_str(), log_file) < 0 ||
                    std::fprintf(log_file, "\nStatus: success\nWarnings: %d\n", tracer.WarningCount()) < 0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "写出事件日志失败: " + log_path.string(), bastion::ERROR);
                    break;
                }

                const int log_close_status = std::fclose(log_file);
                log_file = nullptr;

                if (log_close_status != 0)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "写出日志失败: " + log_path.string(), bastion::ERROR);
                    break;
                }
            }

            std::cout << "O2M completed: " << results.size() << " records, " << tracer.WarningCount() << " warning(s)\n";

            if (!log_path.empty() && tracer.WarningCount() > 0) std::cout << "See log for model warnings.\n";

            // 模型警告已由 tracer 实时提示，完成时保留数量，不重复重放事件。

            return 0;
        }
        while (false);
    }
    catch (const std::exception& error)
    {
        // 文件系统、容器等标准库异常统一转为任务错误，再进入同一清理出口。
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), bastion::ERROR);
    }
    catch (...)
    {
        // 未知异常交回运行环境，先关闭已打开的文件句柄。
        if (output_file) std::fclose(output_file);

        if (log_file) std::fclose(log_file);

        throw;
    }

    // =================================== 09 失败清理与错误留存 ===================================

    if (output_file) std::fclose(output_file);

    if (output_created)
    {
        std::error_code ignored;
        fs::remove(created_output, ignored);
    }

    if (log_file)
    {
        // 失败日志保留事件及调用链；撤回的只是本轮不完整结果。
        std::fputs("\nStatus: failed\n", log_file);
        std::ostringstream events;
        tracer.ReportEvents(events);
        tracer.ReportTrace(events);
        std::fputs(events.str().c_str(), log_file);
        std::fclose(log_file);
    }

    tracer.ReportFailure();
    return 1;
}
