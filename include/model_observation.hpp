/**
 * @file model_observation.hpp
 * @brief 原始观测与测站资料，以及定轨使用的理论观测值和几何偏导。
 * @details observation 中的记录保持文件原单位；determination 中的记录保存计算所用 rad、m 及 sigma。
 * 资料选用、坐标转换和算法适用性检查由使用阶段完成。
 *
 * @par 输入输出约定
 * 测站表每行为站号、坐标类型和三个坐标分量；观测行以 UTC MJD、类型、目标号、观测者号开头。
 * 原始资料角度为 deg，距离及三维位置为 m，距离率为 m/s；目标位置类型 7 为 GCRS，类型 8 为 WGS84 地固 XYZ。
 * 理论观测计算接收同历元的归一化 GCRS 轨道及计算观测，返回理论观测量和瞬时几何偏导。
 * 归一化长度尺度为 ae，时间尺度为 UNITTIME=sqrt(ae³/GE)，角度为 rad。
 */
#ifndef CDC_MODEL_OBSERVATION_HPP
#define CDC_MODEL_OBSERVATION_HPP

#include "utils_overwatch.hpp"
#include "orbit_statedef.hpp"
#include <Eigen/Core>
#include <array>
#include <ostream>
#include <string>
#include <vector>

namespace observation
{

/**
 * @brief 取得观测类型的统一英文名称。
 * @param[in] type OBS 文件中的类型编号，1—8 对应地基测角、测距、测速、天基测角或目标位置观测。
 * @return 类型名称的常量字符串；编号不在 1—8 内时返回 Unknown。
 * @details 供原始资料概览、定轨逐轮统计及最终残差统计共用，不读取或修改观测数据。
 */
const char* ObservationTypeName(int type);

// 一条地面测站记录，保持文件中的原始单位。
struct ObserverRecord
{
    std::string id;

    // 1：大地坐标；2：ITRS 直角坐标。
    int coordinate_type = 0;

    // 类型 1 保留原文件顺序：大地纬度 phi、经度 lambda、高度，单位 deg、deg、m。
    // 读取时不做转换；使用阶段采用项目 constants 中现有参数定义的参考椭球。
    double geodetic_phi_deg = 0.0;
    double geodetic_lambda_deg = 0.0;
    double geodetic_height_m = 0.0;

    // 类型 2 使用 ITRS 下的 X、Y、Z，单位为 m。
    std::array<double, 3> position_itrs_m = {};
};

// 管理整份地面测站资料。
class Observers
{
private:
    std::vector<ObserverRecord> m_records;

public:
    /**
     * @brief 读取地面测站表，保存站号、坐标类型和原始坐标。
     * @param[in,out] tracer 记录文件打开、列格式、非有限坐标及重复站号错误。
     * @param[in] file 测站表路径。每行依次为站号、坐标类型和三个坐标分量。
     * 类型 1 为大地纬度 phi(deg)、经度 lambda(deg)、高度(m)；类型 2 为 ITRS XYZ(m)。
     * @details 跳过空行和 # 注释，检查列数、坐标有限性及站号唯一性。
     * 全部读取成功后替换 m_records；失败时保留已有测站资料。读取阶段保存原单位。
     */
    void ReadFile(bastion::Overwatch& tracer, const std::string& file);

    /**
     * @brief 取得当前保存的原始测站记录。
     * @return m_records 的只读引用；记录顺序与站表一致，大地坐标为 deg、deg、m，ITRS 坐标为 m。
     */
    const std::vector<ObserverRecord>& GetRecords() const;
};

// 一条原始观测记录。
// 观测类型决定哪些字段有效，其余字段不参与计算。
struct ObservationRecord
{
    double mjd_utc = 0.0;
    int type = 0;

    std::string object_id;
    std::string observer_id;

    // 类型 1、6：赤经和赤纬，单位为 deg。
    double ra_deg = 0.0;
    double dec_deg = 0.0;

    // 类型 2、3：方位角和高度角，单位为 deg。
    double az_deg = 0.0;
    double el_deg = 0.0;

    // 类型 3、4：距离；类型 5：距离率。
    double range_m = 0.0;
    double range_rate_mps = 0.0;

    // 类型 7、8：目标自身的位置，顺序为 X、Y、Z。
    std::array<double, 3> position_gcrs_m = {};
    // 类型 8 保留原始 WGS84 地固 XYZ，不在读取和写出时执行基准转换。
    // 此资料约定与 Geodetic 类采用的项目参考椭球分别说明。
    std::array<double, 3> position_wgs84_ecef_m = {};

    // 类型 6：观测历元的平台 GCRS 位置，顺序为 X、Y、Z。
    std::array<double, 3> observer_position_gcrs_m = {};

    // 类型 1—5：关联得到的原始地面测站记录。
    ObserverRecord ground_observer;

    // 类型 1—5 表示已匹配地面站；类型 6 表示平台位置已提供。
    // 类型 7、8 不使用此标志。
    // 此标志不表示已转换坐标或已通过定轨适用性检查。
    bool is_observer_attached = false;
};

// 读写和管理整份原始观测资料，并按需建立观测与测站的关联。
class Observations
{
private:
    std::vector<ObservationRecord> m_records;

public:
    /**
     * @brief 读取类型 1—8 的原始观测资料，保存文件顺序、编号和数值。
     * @param[in,out] tracer 记录文件打开、未知类型、列格式及非有限数值错误。
     * @param[in] file 观测文件路径。每行前四列为 UTC MJD、类型、目标号和观测者号，
     * 后续列按类型解释为角度(deg)、距离(m)、距离率(m/s)或三维位置(m)。
     * @details 跳过空行和 # 注释，按类型检查列数及数值有限性。
     * 全部读取成功后替换 m_records；失败时保留已有观测。读取后可另行关联测站，
     * 资料选用、物理范围检查及坐标和单位转换由使用阶段完成。
     */
    void ReadFile(bastion::Overwatch& tracer, const std::string& file);

    /**
     * @brief 向资料末尾添加一条已经准备好的原始观测。
     * @param[in] record 待复制的观测记录，保留其中的 UTC MJD、类型、编号、原单位数值及测站关联状态。
     * @details 将 record 的副本追加到 m_records；本函数只保存记录，不执行格式或物理范围检查。
     */
    void AddRecord(const ObservationRecord& record);

    /**
     * @brief 将当前观测按类型 1—8 的原始资料格式写入文件。
     * @param[in,out] tracer 记录编号、类型、非有限数值以及文件打开和写入错误。
     * @param[in] file 输出文件路径；打开成功后覆盖同名文件。
     * @details 按 m_records 顺序写出 UTC MJD、类型、目标号、观测者号及对应类型的数据列。
     * 角度为 deg，距离和位置为 m，距离率为 m/s；数值采用 17 位有效数字，字段以空格分隔。
     * 写出前检查编号、类型及数值有限性，观测内容和测站关联状态保持不变。
     */
    void WriteFile(bastion::Overwatch& tracer, const std::string& file) const;

    /**
     * @brief 取得当前保存的原始观测记录。
     * @return m_records 的只读引用，保留资料顺序、原单位数值及当前测站关联状态。
     */
    const std::vector<ObservationRecord>& GetRecords() const;

    /**
     * @brief 按观测者编号，为地基观测关联原始测站记录。
     * @param[in,out] tracer 读取调用链的已有错误状态，并记录本次测站关联的调用过程。
     * @param[in] observers 地面测站集合；按站号读取记录，并复制其坐标类型和原单位坐标。
     * @details 逐条清除 m_records 中原有的 ground_observer 和 is_observer_attached，再重新关联。
     * 类型 1—5 按站号查找并复制测站，未匹配时标记为 false；类型 6 标记为 true，使用记录自带的平台位置。
     * 类型 7、8 保持 false。匹配状态供使用阶段检查，本函数不转换坐标或单位。
     */
    void AttachObservers(bastion::Overwatch& tracer, const Observers& observers);

    /**
     * @brief 输出已保存观测的总条数、各类型条数及历元范围。
     * @param[in,out] output 接收文字概览的输出流；本函数临时调整名称对齐和历元精度，结束时恢复流的格式。
     * @details 统计整份 m_records，类型编号应在 1—8 内；遍历全部记录确定最早、最晚 UTC MJD。
     * 历元按 12 位小数显示，资料为空时输出提示。
     */
    void PrintSummary(std::ostream& output) const;
};
} // namespace observation

// 原始记录在 observation 中保存文件单位；以下计算记录供定轨使用 rad 和 m。
namespace determination
{
// 进入定轨后才使用这份计算表示：测角 rad，位置 m；原始文件仍由 Observations 保存。
enum class ObservationType { GroundRaDec = 1, SpaceRaDec = 6, InertialPosition = 7 };
struct Observation
{
    double mjd_utc = 0;
    ObservationType type = ObservationType::GroundRaDec;
    std::string object_id, observer_id;
    // 类型 1/6 为 [RA,Dec,0](rad)，类型 7 为目标 GCRS [X,Y,Z](m)。
    Eigen::Vector3d value = Eigen::Vector3d::Zero();

    // 测角采用东西、南北方向的标准差，单位 rad；惯性位置为 m。
    Eigen::Vector3d sigma = Eigen::Vector3d::Ones();

    // 类型 1/6 的观测者 GCRS 位置，单位 m；位置观测不使用此量。
    Eigen::Vector3d observer_position_gcrs_m = Eigen::Vector3d::Zero();
};

/**
 * @brief 计算一条观测的理论观测量和瞬时几何偏导。
 * @param[in,out] tracer 记录观测类型或几何定义域错误，调用方通过 Code() 检查。
 * @param[in] obs 计算观测记录。类型 1/6 读取观测者 GCRS 位置(m)及观测 Dec(rad)，
 * 类型 7 只使用观测类型；sigma 由后续加权阶段使用。
 * @param[in] state 观测历元的归一化 GCRS Cartesian 状态，分量为 [x,y,z,vx,vy,vz]，
 * 长度尺度为 ae，时间尺度为 UNITTIME。
 * @param[out] calc_obs 理论观测量。类型 1/6 为 [RA,Dec,0](rad)，类型 7 为目标 GCRS [X,Y,Z](m)。
 * @param[out] partials 3×6 的瞬时几何偏导，列对应归一化 Cartesian 六分量，速度列为零。
 * 测角前两行分别为 cos(观测Dec)*dRA/dx、dDec/dx，供东西、南北方向残差使用；
 * 位置观测的三行为 d[X,Y,Z]/dx，单位为 m。
 * @details 采用同历元几何模型，计算理论值及其几何偏导；O-C 由批处理残差循环计算。
 * 测角第一行偏导含观测 Dec 的余弦因子，calc_obs 中仍保存 RA 本身。
 * 模型不含光行时、光行差和折射，零距离或极点附近的测角几何通过 tracer 报错。
 */
void Compute_ObsAndPartials(bastion::Overwatch& tracer, const Observation& obs, const orbit::Cartesian& state,
    Eigen::Vector3d& calc_obs, Eigen::Matrix<double, 3, 6>& partials);
} // namespace determination

#endif
