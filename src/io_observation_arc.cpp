/**
 * @file io_observation_arc.cpp
 * @brief 按相同历元顺序写出六列 OBS 和可选 Keplerian、Cartesian 真值。
 */
#include "io_observation_arc.hpp"
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace bastion;

namespace observation
{
/**
 * @brief 写出一条理论测角弧段及其同历元轨道真值。
 * @param[in,out] tracer 记录文件打开、写入和关闭时的错误。
 * @param[in] arc 非空弧段；观测时刻为 UTC MJD，GCRS 赤经、赤纬为 deg。
 *     kep 顺序为 [a,e,i,OMEGA,omega,M]，单位 km、deg；cart 为 [x,y,z,vx,vy,vz]，单位 km、km/s。
 * @param[in] output_dir 已存在的结果目录。
 * @param[in] object_id 目标标识，用于构造文件名。
 * @param[in] station_label 用于文件名的测站标识；观测正文保留原始站号。
 * @param[in] arc_number 本次输出的弧段编号，从 1 开始。
 * @param[in] output_kep 1 表示写出同历元 Keplerian 真值，0 表示关闭。
 * @param[in] output_cart 1 表示写出同历元 Cartesian 真值，0 表示关闭。
 * @return 成功时返回 OBS 文件路径，发生错误时返回空字符串并设置 tracer。
 * @details 观测使用标准六列格式；三类文件共用后缀，各行历元和顺序一一对应。
 */
string WriteObservationArc(Overwatch& tracer, const vector<SimulatedObservationRecord>& arc,
                           const std::filesystem::path& output_dir,
                           const string& object_id, const string& station_label,
                           size_t arc_number, int output_kep, int output_cart)
{
    // 文件名取首条观测的 UTC 历元，弧号为五位，MJD 保留十位小数。
    ostringstream suffix;
    suffix << object_id << '_' << station_label << '_' << setw(5) << setfill('0') << arc_number
           << '_' << fixed << setprecision(10) << arc.front().observation.mjd_utc << ".txt";

    // 观测文件总是写出，另外两类真值按开关打开；三者共用同一后缀保证可对应。
    const std::filesystem::path obs_file = output_dir / ("obs_" + suffix.str());
    observation::Observations arc_observations;
    ofstream kep_file, cart_file;

    // 按开关打开轨道真值文件，失败时由 tracer 报告具体输出路径。
    if (output_kep) kep_file.open(output_dir / ("kep_" + suffix.str()));

    if (output_cart) cart_file.open(output_dir / ("cart_" + suffix.str()));

    if ((output_kep && !kep_file) || (output_cart && !cart_file))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open arc truth output: " + obs_file.string(), EventType::ERROR);
        return "";
    }

    // 每条记录写一行观测，并按相同顺序写可选真值，不重新计算或修改角度和轨道数据。
    for (const SimulatedObservationRecord& sample : arc)
    {
        arc_observations.AddRecord(sample.observation);

        // Keplerian 真值保持原六列顺序 a,e,i,OMEGA,omega,M，历元与观测行一致。
        if (output_kep)
        {
            kep_file << setiosflags(ios::left | ios::showpoint) << setw(22) << setprecision(15) << sample.observation.mjd_utc;

            for (double value : sample.kep) kep_file << setw(22) << value;

            kep_file << '\n';
        }

        // Cartesian 真值按位置、速度顺序写出，采用十五位有效数字。
        if (output_cart)
        {
            cart_file << setiosflags(ios::left | ios::showpoint) << setw(22) << setprecision(15) << sample.observation.mjd_utc;

            for (double value : sample.cart) cart_file << setw(22) << value;

            cart_file << '\n';
        }
    }

    // 正式 OBS 的公共前四列及各类型数据列统一由资料类写出，真值仍由仿真程序管理。
    arc_observations.WriteFile(tracer, obs_file.string());
    if (tracer.Code() != ERR_NOERROR) return "";

    // 显式关闭后再检查流状态，确保写入或关闭错误不会被当作成功登记到输出文件清单。
    if (output_kep) kep_file.close();

    if (output_cart) cart_file.close();

    // 这里不递增弧号、不清空输入；main 只在本函数成功返回后推进弧段生命周期。
    if ((output_kep && !kep_file) || (output_cart && !cart_file))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing arc: " + obs_file.string(), EventType::ERROR);
        return "";
    }

    return obs_file.string();
}

} // namespace observation
