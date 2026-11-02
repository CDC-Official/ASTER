/**
 * @file io_observation_arc.hpp
 * @brief 写出一条地基理论测角弧段及其同历元轨道真值。
 * @details 观测和轨道真值采用同一后缀与记录顺序，供后续定轨和误差评价读取。
 */
#pragma once
#include "model_observation.hpp"
#include <array>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace observation
{
// 仿真时把原始观测与同历元真值放在一起，避免分弧后丢失对应关系。
// 当前生成地基 GCRS 赤经、赤纬的理论值，观测历元采用 UTC MJD。
struct SimulatedObservationRecord
{
    int input_id = 0;
    observation::ObservationRecord observation;
    std::array<double, 6> kep{};       // a,e,i,OMEGA,omega,M；单位 km、deg。
    std::array<double, 6> cart{};      // x,y,z,vx,vy,vz；单位 km、km/s。
};

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
std::string WriteObservationArc(bastion::Overwatch& tracer, const std::vector<SimulatedObservationRecord>& arc,
                                const std::filesystem::path& output_dir,
                                const std::string& object_id, const std::string& station_label,
                                std::size_t arc_number, int output_kep, int output_cart);
} // namespace observation
