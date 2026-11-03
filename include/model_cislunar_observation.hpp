/**
 * @file model_cislunar_observation.hpp
 * @brief 地月空间状态记录、单程测距记录及科学计算接口。
 * @details 状态历元为 TDB，位置与速度为 E-BCRS，采用 m、m/s、s。
 * 钟差定义为卫星钟读数减 TDB；接收钟读数与真实接收历元分别保存。
 */
#ifndef ASTER_MODEL_CISLUNAR_OBSERVATION_HPP
#define ASTER_MODEL_CISLUNAR_OBSERVATION_HPP

#include "utils_overwatch.hpp"
#include <array>
#include <functional>
#include <vector>

namespace Cislunar
{

// 一个 TDB 历元的状态；秒字段保留小数，位置、速度和钟差对应同一事件。
struct StateRecord
{
    // 绝对 TDB 历元：整数保存年月日时分，double 保存带小数的秒。
    int year = 0, month = 0, day = 0;
    int hour = 0, minute = 0;
    double second = 0.0;

    // 同一历元的 E-BCRS 位置和速度，分量顺序均为 x、y、z。
    std::array<double, 3> position_ebcrs_m{};
    std::array<double, 3> velocity_ebcrs_mps{};
    double clock_offset_s = 0.0; // 卫星钟读数减 TDB，单位 s。
};

// 一条观测及其计算分项。发收状态的历元为 TDB，接收钟读数包含接收钟差。
struct OneWayRangeRecord
{
    // 观测文件使用的发收编号、接收钟读数和最终伪距。
    int transmitter_id = 0, receiver_id = 0;
    std::array<int, 5> receive_clock_calendar{}; // 年、月、日、时、分。
    double receive_clock_second = 0.0;
    double pseudorange_m = 0.0; // 最终观测伪距，包含调用方加入的误差。

    // 光行时求解得到的发射状态，以及给定的真实接收状态。
    StateRecord transmit_state, receive_state;
    // 理论伪距的三个贡献：几何距离、引力时延对应的距离、两端钟差贡献。
    double geometric_m = 0.0; // 两个发收事件的 BCRS 位置之间的直线距离。
    double gravity_delay_s = 0.0; // 引力传播时延，乘光速后得到距离贡献。
    double clock_contribution_m = 0.0; // c ×（接收钟差 − 发射钟差）。
    double theoretical_pseudorange_m = 0.0; // 加入系统差和白噪声之前的理论伪距。
    // 主程序加入的系统差和本条观测抽取的白噪声，单位 m。
    double systematic_bias_m = 0.0, noise_m = 0.0;

    // 光行时取自当前发收事件；update 保存末次更新量的绝对值，单位 s。
    double light_time_s = 0.0, light_time_update_s = 0.0;
    int iterations = 0;
};

/**
 * @brief 从离散状态中插值指定 TDB 时刻的位置、速度和钟差。
 * @param[in,out] tracer 记录节点不足、时间乱序、越界及非有限状态错误。
 * @param[in] records 按 TDB 递增的状态序列，至少两点；位置 m、速度 m/s、钟差 s。
 * @param[in,out] state 输入年月日时分秒，输出该历元的 E-BCRS 位置、速度和总钟差。
 * @details 选取查询时刻附近最多十个节点，逐分量使用 Neville 多项式插值，最高九次。
 * 节点间隔可不等；位置、速度和钟差分别插值。时间差由整日差、整分钟差和秒差组成。
 */
void InterpolateState(bastion::Overwatch& tracer, const std::vector<StateRecord>& records, StateRecord& state);

/**
 * @brief 计算多个天体单极场产生的单程引力传播时延。
 * @param[in,out] tracer 记录天体编号、重复天体及无效传播几何错误。
 * @param[in] transmit_m 发射事件的 BCRS 位置，单位 m。
 * @param[in] receive_m 接收事件的 BCRS 位置，单位 m。
 * @param[in] body_positions_m 十一个天体在接收时刻的 BCRS 位置，按 JPL 编号减一索引，单位 m。
 * @param[in] gravity_bodies 参与求和的 JPL 天体编号，取 1—11；空列表对应零时延。
 * @param[out] delay_s 引力传播时延，单位 s。
 * @details 各项为 2GM/c³ ln[(r发+r收+距离+修正)/(r发+r收−距离+修正)]。
 * GM 使用公共 DE440 参数；太阳项的修正为 2GM/c²，其余天体为零。
 */
void ComputeGravityDelay(bastion::Overwatch& tracer, const std::array<double, 3>& transmit_m,
    const std::array<double, 3>& receive_m, const std::array<std::array<double, 3>, 11>& body_positions_m,
    const std::vector<int>& gravity_bodies, double& delay_s);

/**
 * @brief 从真实接收事件计算单程理论伪距及发射事件。
 * @param[in,out] tracer 记录状态查询、JPL 覆盖范围及光行时收敛错误。
 * @param[in] receive 接收事件的 TDB 历元、E-BCRS 位置和钟差，单位 m、s。
 * @param[in] query_transmit 根据传入状态的 TDB 历元填写发射端位置、速度和总钟差，并通过 tracer 报错。
 * @param[in] gravity_bodies 参与引力时延的 JPL 天体编号 1—11；空列表关闭引力时延。
 * @param[out] record 发收状态、接收钟读数、几何距离、引力时延、钟差贡献和理论伪距；编号由调用方填写。
 * @details 两端 E-BCRS 位置分别加各自事件的地球 BCRS 位置。引力天体位置取接收 TDB。
 * 以秒为单位迭代光行时，末次变化不超过 1e-12 s 时收敛，最多二十轮。
 * 伪距为几何距离加 c 倍引力时延，再加 c 倍接收减发射钟差；接收钟读数为接收 TDB 加接收钟差。
 */
void ComputeOneWayRange(bastion::Overwatch& tracer, const StateRecord& receive,
    const std::function<void(bastion::Overwatch&, StateRecord&)>& query_transmit,
    const std::vector<int>& gravity_bodies, OneWayRangeRecord& record);

}
#endif
