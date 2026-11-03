/**
 * @file model_cislunar_observation.hpp
 * @brief 地月星间单程测距的轨道、钟差插值和收发事件计算。
 * @details 输入轨道为 E-BCRS/TDB，位置 m、速度 m/s；钟差为星钟减 TDB，单位 s。
 * 时刻以参考 MJD 的整数日、小数日和相对天数保存，传播几何使用 BCRS。
 * 模型采用 10 点 Neville 插值生成 16 阶切比雪夫系数，包含除地球外的十个天体引力时延。
 */
#ifndef ASTER_MODEL_CISLUNAR_OBSERVATION_HPP
#define ASTER_MODEL_CISLUNAR_OBSERVATION_HPP

#include <Eigen/Core>
#include <array>
#include <string>
#include <vector>

namespace observation
{

// 一份等间隔轨道或钟差表；values 每行对应一个历元，轨道六列，钟差两列。
struct CislunarSeries
{
    int satellite_id = 0;
    std::vector<double> days;
    Eigen::MatrixXd values;
    bool is_clock = false;
    int window_number = 1;
    int window_start = 0;
};

// 一段切比雪夫插值，coefficients 的行对应 0—16 阶，列对应输入分量。
struct CislunarInterpolant
{
    double start_day = 0.0;
    double span_days = 0.0;
    Eigen::MatrixXd coefficients;
};

struct CislunarSlot
{
    std::string type;
    int transmitter = 0, receiver = 0;
    int transmit_antenna = 0, receive_antenna = 0;
    double clock_tag_day = 0.0;
};

struct OneWayRange
{
    double transmit_day = 0.0, receive_day = 0.0;
    double geometric_m = 0.0, gravity_m = 0.0, propagation_m = 0.0;
    double transmit_physical_s = 0.0, transmit_relativistic_s = 0.0;
    double receive_physical_s = 0.0, receive_relativistic_s = 0.0;
    double clock_m = 0.0, pseudorange_m = 0.0;
    double receive_update_s = 0.0, transmit_update_s = 0.0;
    int receive_iterations = 0, transmit_iterations = 0;
};

/**
 * @brief 读取收发时隙表，并保存接收端钟时标签。
 * @param[in] path 时隙表路径；T 行给出参考日期，ISLL 行依次为安排发射秒、接收秒、收发星与天线编号。
 * @param[out] reference_mjd 首个 T 行的 TDB 参考 MJD，分别保存整数日和小数日。
 * @return 文件顺序的单程记录；clock_tag_day 为相对于 reference_mjd 的天数。
 * @details 接收秒是接收星钟的时间标签；发射事件由传播方程求解。
 */
std::vector<CislunarSlot> ReadCislunarSlots(const std::string& path, std::array<double, 2>& reference_mjd);

/**
 * @brief 读取一颗卫星的轨道表或钟差表。
 * @param[in] path 含 OBJECT_ID、DATA_START 和 DATA_STOP 的 LHH 文本文件路径。
 * @param[in] reference_mjd 观测参考 TDB MJD，整数日与小数日分别保存。
 * @param[in] is_clock true 读取相对论钟差和物理钟偏两列，单位 s；false 读取 E-BCRS 位置和速度六列，单位 m、m/s。
 * @return 卫星编号、相对参考历元的天数和等间隔采样值，滑动窗口从第一段开始。
 */
CislunarSeries ReadCislunarSeries(const std::string& path, const std::array<double, 2>& reference_mjd, bool is_clock);

/**
 * @brief 为当前接收钟时选取采样窗口，并计算切比雪夫系数。
 * @param[in,out] series 读取 days、values、is_clock，向后更新 window_number 和 window_start；观测须按钟时递增调用。
 * @param[in] clock_tag_day 接收端钟时距参考历元的天数，用于固定本条观测的插值窗口。
 * @return 当前窗口的起始时刻、跨度及 17 行系数；分量和单位与 series 相同。
 * @details 轨道窗口每次前进 9 点；钟差从第二段起按 7 点步长、提前一行取 10 点。
 * 先以 Neville 方法求 17 个 Chebyshev–Gauss–Lobatto 节点值，再计算 16 阶展开。
 */
CislunarInterpolant PrepareCislunarInterpolant(CislunarSeries& series, double clock_tag_day);

/**
 * @brief 在固定的切比雪夫窗口内计算轨道或钟差。
 * @param[in] interpolant 起始时刻、跨度和 0—16 阶系数，单位与原始采样值相同。
 * @param[in] day 待求时刻距参考 TDB 历元的天数。
 * @return 按原始列顺序排列的六维轨道或二维钟差；单位为 m、m/s 或 s。
 * @details 使用 T0=1、T1=x 的递推，零阶系数取一半；收发事件共用接收钟时选定的窗口。
 */
Eigen::VectorXd EvaluateCislunarInterpolant(const CislunarInterpolant& interpolant, double day);

/**
 * @brief 求解一条星间单程观测的收发 TDB 时刻和理论伪距。
 * @param[in] reference_mjd 参考 TDB MJD 的整数日和小数日，用于查询已载入的 JPL 二进制历表。
 * @param[in] clock_tag_day 接收端钟时距参考历元的天数。
 * @param[in] transmit_orbit 发射星 E-BCRS 轨道插值，位置 m、速度 m/s。
 * @param[in] receive_orbit 接收星 E-BCRS 轨道插值，位置 m、速度 m/s。
 * @param[in] transmit_clock 发射星钟差插值，两列为相对论钟差和物理钟偏，单位 s。
 * @param[in] receive_clock 接收星钟差插值，定义及单位与 transmit_clock 相同。
 * @return 收发时刻、几何距离、引力时延距离、两端钟差和伪距，以及迭代次数与末次变化量。
 * @details 先反求接收 TDB，再迭代发射 TDB；单次时刻变化不超过 1e-12 s 时收敛。
 * 两端 E-BCRS 位置分别加各自事件的地球 BCRS 位置。引力天体取接收 TDB，跳过地球，
 * 太阳项的对数分子、分母均加 2GM/c²。伪距等于几何距离、引力时延距离和接收减发射钟差距离之和。
 * 设备时延、天线偏移和测量噪声取零。十轮内未收敛时报告错误。
 */
OneWayRange ComputeCislunarOneWayRange(const std::array<double, 2>& reference_mjd, double clock_tag_day,
    const CislunarInterpolant& transmit_orbit, const CislunarInterpolant& receive_orbit,
    const CislunarInterpolant& transmit_clock, const CislunarInterpolant& receive_clock);

}
#endif
