/**
 * @file io_cislunar_observation.cpp
 * @brief 读取 LHH 星间时隙、E-BCRS 轨道和 TDB 钟差文本资料。
 * @details 日期拆分为整数 MJD 和日内小数；采样时刻保存为相对参考历元的天数。
 */
#include "model_cislunar_observation.hpp"
#include <sofa.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace observation
{

/**
 * @brief 将日历日期转换成拆分 MJD。
 * @param[in] text YYYY-MM-DDThh:mm:ss 格式的 TDB 参考日期或星钟日期。
 * @return 整数 MJD 和日内小数，两者分别保存以保留亚秒分辨率。
 */
static std::array<double, 2> ParseCislunarDate(const std::string& text)
{
    int year, month, day, hour, minute;
    double second, zero, mjd;
    const int fields = std::sscanf(text.c_str(), "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &minute, &second);

    if (fields != 6) throw std::runtime_error("Invalid TDB calendar: " + text);

    const int status = iauCal2jd(year, month, day, &zero, &mjd);
    if (status || hour < 0 || hour > 23 || minute < 0 || minute > 59 || !std::isfinite(second) || second < 0 || second >= 60)
        throw std::runtime_error("Invalid TDB calendar: " + text);

    return {mjd, (hour * 3600.0 + minute * 60.0 + second) / 86400.0};
}

/**
 * @brief 读取收发时隙表，并保存接收端钟时标签。
 * @param[in] path 时隙表路径；T 行给出参考日期，ISLL 行依次为安排发射秒、接收秒、收发星与天线编号。
 * @param[out] reference_mjd 首个 T 行的 TDB 参考 MJD，分别保存整数日和小数日。
 * @return 文件顺序的单程记录；clock_tag_day 为相对于 reference_mjd 的天数。
 * @details 接收秒是接收星钟的时间标签；发射事件由传播方程求解。
 */
std::vector<CislunarSlot> ReadCislunarSlots(const std::string& path, std::array<double, 2>& reference_mjd)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read slot table: " + path);

    std::vector<CislunarSlot> slots;
    std::array<double, 2> current = {};
    bool have_reference = false;
    std::string line, type;
    int line_number = 0;

    while (std::getline(input, line))
    {
        ++line_number;
        std::istringstream row(line);

        if (!(row >> type) || type[0] == '#') continue;

        if (type == "T")
        {
            int year, month, day, hour, minute;
            double second;

            if (!(row >> year >> month >> day >> hour >> minute >> second))
                throw std::runtime_error(path + ": invalid T row " + std::to_string(line_number));

            char date[80];
            std::snprintf(date, sizeof(date), "%04d-%02d-%02dT%02d:%02d:%.9f", year, month, day, hour, minute, second);
            current = ParseCislunarDate(date);

            if (!have_reference) reference_mjd = current;

            have_reference = true;
        }
        else
        {
            CislunarSlot slot;
            double scheduled_transmit_s, receive_s;
            slot.type = type;

            if (!have_reference || type != "ISLL" || !(row >> scheduled_transmit_s >> receive_s >> slot.transmitter
                >> slot.transmit_antenna >> slot.receiver >> slot.receive_antenna) || !std::isfinite(receive_s))
                throw std::runtime_error(path + ": invalid ISLL row " + std::to_string(line_number));

            slot.clock_tag_day = (current[0] - reference_mjd[0]) + (current[1] + receive_s / 86400.0) - reference_mjd[1];
            slots.push_back(slot);
        }
    }

    if (slots.empty()) throw std::runtime_error("Empty slot table: " + path);

    return slots;
}

/**
 * @brief 读取一颗卫星的轨道表或钟差表。
 * @param[in] path 含 OBJECT_ID、DATA_START 和 DATA_STOP 的 LHH 文本文件路径。
 * @param[in] reference_mjd 观测参考 TDB MJD，整数日与小数日分别保存。
 * @param[in] is_clock true 读取相对论钟差和物理钟偏两列，单位 s；false 读取 E-BCRS 位置和速度六列，单位 m、m/s。
 * @return 卫星编号、相对参考历元的天数和等间隔采样值，滑动窗口从第一段开始。
 */
CislunarSeries ReadCislunarSeries(const std::string& path, const std::array<double, 2>& reference_mjd, bool is_clock)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot read satellite table: " + path);

    CislunarSeries series;
    series.is_clock = is_clock;
    const int components = is_clock ? 2 : 6;
    std::vector<Eigen::VectorXd> samples;
    std::string line, word;
    bool in_data = false, data_stop = false;
    int line_number = 0;

    while (std::getline(input, line))
    {
        ++line_number;
        std::istringstream row(line);

        if (!(row >> word)) continue;

        if (word == "OBJECT_ID") row >> word >> series.satellite_id;

        if (word == "DATA_START")
        {
            in_data = true;
            continue;
        }

        if (word == "DATA_STOP")
        {
            data_stop = true;
            break;
        }

        if (!in_data) continue;

        const std::array<double, 2> epoch = ParseCislunarDate(word);
        const double day = (epoch[0] - reference_mjd[0]) + epoch[1] - reference_mjd[1];
        Eigen::VectorXd value(components);

        for (int j = 0; j < components; ++j)
        {
            if (!(row >> value(j)) || !std::isfinite(value(j)))
                throw std::runtime_error(path + ": invalid sample at line " + std::to_string(line_number));
        }

        if (!series.days.empty() && day <= series.days.back())
            throw std::runtime_error(path + ": non-increasing sample at line " + std::to_string(line_number));

        series.days.push_back(day);
        samples.push_back(value);
    }

    if (!data_stop || series.satellite_id <= 0 || samples.size() < 10)
        throw std::runtime_error("Incomplete satellite table: " + path);

    // Neville 节点按等间隔排列，读入时核对采样间隔。
    const double step = series.days[1] - series.days[0];
    series.values.resize(samples.size(), components);

    for (std::size_t i = 0; i < samples.size(); ++i)
    {
        if (i > 0 && std::abs((series.days[i] - series.days[i - 1]) - step) * 86400.0 > 1e-6)
            throw std::runtime_error("Nonuniform sampling: " + path);

        series.values.row(i) = samples[i].transpose();
    }

    return series;
}

}
