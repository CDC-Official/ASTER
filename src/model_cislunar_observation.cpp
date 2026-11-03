/**
 * @file model_cislunar_observation.cpp
 * @brief 插值地月空间状态，计算单程光行时、引力时延和伪距。
 * @details 状态历元采用 TDB 年月日时分秒，位置为 E-BCRS，计算量使用 SI 单位。
 * 传播端点分别转换为各自事件的 BCRS 位置，再求解发射事件。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include "utils_constants.hpp"
#include <sofa.h>
#include <sofam.h>
#include <Eigen/Core>
#include <algorithm>
#include <cmath>

namespace Cislunar
{

/**
 * @brief 将状态历元分为午夜儒略日和日内小数日。
 * @param[in,out] tracer 记录无效日期或时分秒。
 * @param[in] state 读取 TDB 年月日时分秒；秒数应在 [0,60) 内。
 * @param[out] jd 两部分 TDB 儒略日，供 JPL 查询及日差计算使用。
 */
static void EpochParts(bastion::Overwatch& tracer, const StateRecord& state, double jd[2])
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);
    // SOFA 返回儒略日基准和该日期的 MJD，二者相加得到当天午夜的儒略日。
    double origin, day;
    const int status = iauCal2jd(state.year, state.month, state.day, &origin, &day);

    // 检查日期和日内时间；TDB 的秒字段取 0 到 60 之间的数值。
    if (status != 0 || state.hour < 0 || state.hour > 23 || state.minute < 0 || state.minute > 59
        || !std::isfinite(state.second) || state.second < 0.0 || state.second >= 60.0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid TDB calendar epoch.", bastion::EventType::ERROR);
        return;
    }

    // 午夜日期与日内时间分别保存，后续查询将两部分一起传给 JPL。
    jd[0] = origin + day;
    jd[1] = (state.hour * 3600.0 + state.minute * 60.0 + state.second) / DAYSEC;
    tracer.CheckOut(__func__);
}

/**
 * @brief 将历元移动指定秒数，保留秒的小数部分。
 * @param[in,out] tracer 记录日期转换错误。
 * @param[in] epoch 起始 TDB 历元，只读取年月日时分秒。
 * @param[in] offset_s 时间增量，单位 s，可为负数。
 * @param[in,out] shifted 只更新年月日时分秒，位置、速度和钟差保持原值。
 * @details 先在秒字段上加增量，再将整分钟进位到日期，跨日时使用 SOFA 日历转换。
 */
static void ShiftEpoch(bastion::Overwatch& tracer, const StateRecord& epoch, double offset_s, StateRecord& shifted)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);
    double jd[2];
    EpochParts(tracer, epoch, jd);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    if (!std::isfinite(offset_s))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite epoch offset.", bastion::EventType::ERROR);
        return;
    }

    // 先把时间增量加到秒字段，再用整分钟数表示进位或借位。
    const double seconds = epoch.second + offset_s;
    const double minute_shift = std::floor(seconds / 60.0);
    // 将进位计入时分，分出整日偏移和目标日期内的分钟数。
    const double minutes = epoch.hour * 60.0 + epoch.minute + minute_shift;
    const double day_shift = std::floor(minutes / 1440.0);
    const int minute_of_day = static_cast<int>(minutes - day_shift * 1440.0);
    // 用整日偏移更新年月日；秒的小数部分留在 seconds 中。
    double fraction;
    const int status = iauJd2cal(jd[0], day_shift, &shifted.year, &shifted.month, &shifted.day, &fraction);

    // 确认目标日期和累加后的秒数有效。
    if (status != 0 || !std::isfinite(seconds))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot shift TDB epoch.", bastion::EventType::ERROR);
        return;
    }

    // 把日内分钟拆回时分，余下的秒连同小数写入目标历元。
    shifted.hour = minute_of_day / 60;
    shifted.minute = minute_of_day % 60;
    shifted.second = seconds - minute_shift * 60.0;
    tracer.CheckOut(__func__);
}

/**
 * @brief 按 TDB 历元查询天体的 BCRS 位置。
 * @param[in,out] tracer 记录历表覆盖范围错误。
 * @param[in] epoch 读取 TDB 年月日时分秒。
 * @param[in] body JPL 天体编号，1—11，其中 3 为地球、10 为月球、11 为太阳。
 * @param[out] position_m 天体的 BCRS 三维位置，单位 m。
 * @details 使用主程序已经载入的 JPL 历表，以两部分儒略日查询。
 */
static void BodyPosition(bastion::Overwatch& tracer, const StateRecord& epoch, int body, Eigen::Vector3d& position_m)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);
    double jd[2];
    EpochParts(tracer, epoch, jd);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 将历表起止日期减去查询日期的午夜值，再与日内小数日比较。
    if (jd[1] < JPLEPH_R1.r1.ss[0] - jd[0] || jd[1] > JPLEPH_R1.r1.ss[1] - jd[0])
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "TDB epoch outside loaded JPL ephemeris.", bastion::EventType::ERROR);
        return;
    }

    // JPL 的中心编号 12 表示太阳系质心，返回状态的前三项是位置。
    double state_km[6];
    pleph(jd[0], jd[1], body, 12, state_km);
    // 将历表返回的位置从 km 换成 m。
    position_m = Eigen::Vector3d(state_km[0], state_km[1], state_km[2]) * 1000.0;
    tracer.CheckOut(__func__);
}

/**
 * @brief 从离散状态中插值指定 TDB 时刻的位置、速度和钟差。
 * @param[in,out] tracer 记录节点不足、时间乱序、越界及非有限状态错误。
 * @param[in] records 按 TDB 递增的状态序列，至少两点；位置 m、速度 m/s、钟差 s。
 * @param[in,out] state 输入年月日时分秒，输出该历元的 E-BCRS 位置、速度和总钟差。
 * @details 选取查询时刻附近最多十个节点，逐分量使用 Neville 多项式插值，最高九次。
 * 节点间隔可不等；位置、速度和钟差分别插值。时间差由整日差、整分钟差和秒差组成。
 */
void InterpolateState(bastion::Overwatch& tracer, const std::vector<StateRecord>& records, StateRecord& state)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 本函数至少使用两个节点进行状态插值。
    if (records.size() < 2)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "State interpolation needs at least two records.", bastion::EventType::ERROR);
        return;
    }

    // 以查询时刻为零点，分别计算整日、整分钟及秒的差值。
    double query_jd[2];
    EpochParts(tracer, state, query_jd);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // times 保存各采样点相对查询时刻的秒数，负值表示采样点在查询时刻之前。
    std::vector<double> times(records.size());

    for (std::size_t i = 0; i < records.size(); ++i)
    {
        // 逐点检查日历历元，并取得该采样日期的午夜儒略日。
        const StateRecord& record = records[i];
        double jd[2];
        EpochParts(tracer, record, jd);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        // 整日差与时分秒差分别计算，最后合成以秒为单位的时间间隔。
        times[i] = (jd[0] - query_jd[0]) * DAYSEC + (record.hour - state.hour) * 3600.0
            + (record.minute - state.minute) * 60.0 + (record.second - state.second);

        // 采样时刻必须严格递增，重复或倒序时记录具体节点。
        if (i > 0 && times[i] <= times[i - 1])
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "State epochs must increase: record " + std::to_string(i + 1), bastion::EventType::ERROR);
            return;
        }
    }

    // 查询时刻应落在首末采样点之间，包含两个端点。
    if (times.front() > 0.0 || times.back() < 0.0)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Requested TDB epoch outside state records.", bastion::EventType::ERROR);
        return;
    }

    // 窗口围绕查询时刻选取；靠近两端时使用端部节点。
    const int count = static_cast<int>(std::min<std::size_t>(10, records.size()));
    // 找到不早于查询时刻的首个节点，用它确定局部窗口的位置。
    const int right = static_cast<int>(std::lower_bound(times.begin(), times.end(), 0.0) - times.begin());
    const int start = std::clamp(right - count / 2, 0, static_cast<int>(records.size()) - count);
    // 七个插值分量依次为三维位置、三维速度和总钟差。
    std::array<double, 7> value{};

    for (int component = 0; component < 7; ++component)
    {
        // 每次取一个分量的节点值，work 同时保存 Neville 递推的中间结果。
        double work[10];

        for (int i = 0; i < count; ++i)
        {
            // 按分量编号从当前节点取出位置、速度或钟差。
            const StateRecord& record = records[start + i];
            work[i] = component < 3 ? record.position_ebcrs_m[component]
                : (component < 6 ? record.velocity_ebcrs_mps[component - 3] : record.clock_offset_s);

            // 记录当前窗口内非有限数值所在的原始采样行。
            if (!std::isfinite(work[i]))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite state at record " + std::to_string(start + i + 1), bastion::EventType::ERROR);
                return;
            }
        }

        // Neville 递推在零点求值，节点时刻保留实际采样间隔。
        for (int order = 1; order < count; ++order)
        {
            for (int i = 0; i < count - order; ++i)
            {
                // 用左右两个低一阶插值值，构造跨越这两个端点的高一阶插值值。
                const double left = times[start + i];
                const double right_time = times[start + i + order];
                work[i] = (right_time * work[i] - left * work[i + 1]) / (right_time - left);
            }
        }

        // 完成全部阶次后，首项就是该分量在查询时刻的插值结果。
        value[component] = work[0];
    }

    // 将七个插值结果写回状态记录，历元保持为传入的查询时刻。
    for (int i = 0; i < 3; ++i)
    {
        state.position_ebcrs_m[i] = value[i];
        state.velocity_ebcrs_mps[i] = value[i + 3];
    }

    state.clock_offset_s = value[6];
    tracer.CheckOut(__func__);
}

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
    const std::vector<int>& gravity_bodies, double& delay_s)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);
    // 把发收位置作为三维向量使用，求两个事件之间的直线距离。
    const Eigen::Map<const Eigen::Vector3d> transmit(transmit_m.data()), receive(receive_m.data());
    const double distance = (receive - transmit).norm();
    const double c = SPEED_OF_LIGHT_MPS;
    // 记录已经计入的天体，引力时延从零开始逐项累加。
    std::array<bool, 11> included{};
    delay_s = 0.0;

    for (int body : gravity_bodies)
    {
        // 参与求和的天体编号取 1—11，同一天体只计入一次。
        if (body < 1 || body > 11 || included[body - 1])
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid or repeated gravity body: " + std::to_string(body), bastion::EventType::ERROR);
            return;
        }

        // 求发射端、接收端到当前天体的距离之和。
        included[body - 1] = true;
        const Eigen::Map<const Eigen::Vector3d> position(body_positions_m[body - 1].data());
        const double radius_sum = (transmit - position).norm() + (receive - position).norm();
        // 取当前天体的 GM；太阳项在对数分子、分母中各加 2GM/c²。
        const double gm = MAJOR_BODY_GM[body - 1];
        const double correction = body == 11 ? 2.0 * gm / (c * c) : 0.0;
        const double denominator = radius_sum - distance + correction;

        // 对数的分母应为有限正数，几何异常时报告对应天体。
        if (!std::isfinite(denominator) || denominator <= 0.0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid gravitational delay geometry for body " + std::to_string(body), bastion::EventType::ERROR);
            return;
        }

        // log1p 在距离远小于端点至天体距离时仍能保留有效数字。
        delay_s += 2.0 * gm / (c * c * c) * std::log1p(2.0 * distance / denominator);
    }

    tracer.CheckOut(__func__);
}

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
    const std::vector<int>& gravity_bodies, OneWayRangeRecord& record)
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);
    // 初始化本次计算记录，保存给定的真实接收事件。
    record = OneWayRangeRecord{};
    record.receive_state = receive;
    const double c = SPEED_OF_LIGHT_MPS;

    // =================================== 01 准备接收事件及引力天体 ===================================
    // 查询接收时刻的地球位置，用于把接收端坐标平移到太阳系质心。
    Eigen::Vector3d earth_receive;
    BodyPosition(tracer, receive, 3, earth_receive);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 接收事件在整轮光行时迭代中固定，位置只需转换一次。
    const Eigen::Vector3d receive_bcrs = Eigen::Map<const Eigen::Vector3d>(receive.position_ebcrs_m.data()) + earth_receive;
    std::array<double, 3> receive_m{receive_bcrs[0], receive_bcrs[1], receive_bcrs[2]};
    // 引力天体的位置按 JPL 编号减一存放，供每轮时延计算使用。
    std::array<std::array<double, 3>, 11> body_positions{};

    for (int body : gravity_bodies)
    {
        // 查询前检查天体编号，数组下标与 JPL 编号一一对应。
        if (body < 1 || body > 11)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid gravity body: " + std::to_string(body), bastion::EventType::ERROR);
            return;
        }

        // 统一在接收时刻取得引力天体位置，并在本次迭代中固定使用。
        Eigen::Vector3d position;
        BodyPosition(tracer, receive, body, position);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        body_positions[body - 1] = {position[0], position[1], position[2]};
    }

    // =================================== 02 迭代发射事件与光行时 ===================================
    // 从零光行时开始，首轮在接收时刻查询发射端，随后逐步回推发射事件。
    double light_time_s = 0.0;
    bool converged = false;

    for (int iteration = 1; iteration <= 20; ++iteration)
    {
        // 本轮发射时刻由真实接收时刻减去当前光行时得到。
        ShiftEpoch(tracer, receive, -light_time_s, record.transmit_state);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        // 查询函数按这个发射时刻填写发射端的位置、速度和钟差。
        query_transmit(tracer, record.transmit_state);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        // 发射端的坐标平移使用发射时刻的地球位置。
        Eigen::Vector3d earth_transmit;
        BodyPosition(tracer, record.transmit_state, 3, earth_transmit);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        // 两端都表达为 BCRS 位置后，计算不同事件之间的几何距离。
        const Eigen::Vector3d transmit_bcrs = Eigen::Map<const Eigen::Vector3d>(record.transmit_state.position_ebcrs_m.data()) + earth_transmit;
        const std::array<double, 3> transmit_m{transmit_bcrs[0], transmit_bcrs[1], transmit_bcrs[2]};
        record.geometric_m = (receive_bcrs - transmit_bcrs).norm();

        // 用当前发收几何计算引力传播时延，并计入下一轮光行时。
        ComputeGravityDelay(tracer, transmit_m, receive_m, body_positions, gravity_bodies, record.gravity_delay_s);
        if (tracer.Code() != bastion::ERR_NOERROR) return;

        // 新光行时由几何距离除以光速，再加引力时延得到。
        const double next_light_time_s = record.geometric_m / c + record.gravity_delay_s;
        if (!std::isfinite(next_light_time_s))
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite one-way light time.", bastion::EventType::ERROR);
            return;
        }

        // 记录新旧光行时的差值，以及产生当前发射状态所用的光行时。
        record.light_time_update_s = std::abs(next_light_time_s - light_time_s);
        record.light_time_s = light_time_s;
        record.iterations = iteration;

        // 光行时变化达到收敛阈值时，保存当前事件和分项并结束迭代。
        if (record.light_time_update_s <= 1e-12)
        {
            converged = true;
            break;
        }

        // 尚未收敛时，用新光行时重新求下一轮的发射事件。
        light_time_s = next_light_time_s;
    }

    // 二十轮结束仍未达到阈值时，通过 Tracer 报告求解失败。
    if (!converged)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "One-way light time did not converge in 20 iterations.", bastion::EventType::ERROR);
        return;
    }

    // =================================== 03 合成伪距与接收钟读数 ===================================
    // 钟差贡献取接收钟差减发射钟差，乘光速后换成距离。
    record.clock_contribution_m = c * (receive.clock_offset_s - record.transmit_state.clock_offset_s);
    // 将几何距离、引力时延对应的距离和钟差贡献相加，得到理论伪距。
    record.theoretical_pseudorange_m = record.geometric_m + c * record.gravity_delay_s + record.clock_contribution_m;
    record.pseudorange_m = record.theoretical_pseudorange_m;

    // 检查合成后的伪距是否为有限数值。
    if (!std::isfinite(record.pseudorange_m))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Non-finite one-way pseudorange.", bastion::EventType::ERROR);
        return;
    }

    // 在真实接收 TDB 上加接收钟差，得到接收端实际记录的钟读数。
    StateRecord clock_epoch;
    ShiftEpoch(tracer, receive, receive.clock_offset_s, clock_epoch);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 将接收钟读数的年月日时分与小数秒分别写入观测记录。
    record.receive_clock_calendar = {clock_epoch.year, clock_epoch.month, clock_epoch.day, clock_epoch.hour, clock_epoch.minute};
    record.receive_clock_second = clock_epoch.second;
    tracer.CheckOut(__func__);
}

}
