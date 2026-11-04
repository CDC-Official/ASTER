# 公共接口参考

本页面向要在 C++ 中复用 ASTER 的调用者。先说明如何链接和检查错误，再按轨道计算流程介绍主要入口。参数说明中的“输入”由调用方准备，“输出”由函数写入，“输入输出”表示函数读取后还会修改。完整声明与公式注释保留在各头文件，文件协议见[数据格式](05_data-formats.md)。

## 按任务找接口

| 要完成的计算 | 主要入口 | 说明位置 |
| --- | --- | --- |
| 状态表示与单位转换 | `orbit::Keplerian`、`Cartesian`、`NonSingular2` | [状态对象](#states) |
| 历元、EOP、测站坐标 | `UTC2TT`、`EOP::Get_EOP`、`Geodetic::ToGCRS` | [时间与坐标](#frames) |
| 地月状态插值与单程测距 | `Cislunar::InterpolateState`、`ComputeOneWayRange` | [地月测距](#cislunar) |
| 读写原始观测和站表 | `observation::Observations`、`Observers` | [原始资料](#observations) |
| 求理论观测和几何偏导 | `determination::Compute_ObsAndPartials` | [观测模型](#observation-model) |
| 整段预报及观测时刻取样 | `OrbitPropagation_NSG2`、`Interpolate_CartesianOrbit` | [数值预报](#propagation) |
| GL 初轨确定 | `iod_solver::InitialOrbitDetermination` | [初轨求解](#iod) |
| 批处理定轨和协方差传递 | `BatchOrbitDetermination_NSG2`、`MovePosteriorToEpoch` | [批处理](#batch) |
| 密切根数转均根数 | `orbit::OsculatingToMean` | [瞬转平](#averaging) |
| JPL 日月/行星状态 | `Read_JPLEPH`、`Get_JPLEPH` | [JPL 星历](#ephemeris) |
| 读配置与检查数值 | `io::SimpleInput` | [配置读取](#settings) |
| 动力学右端、单步积分和 STM | `Acceleration_*`、`RKF78_*`、`Adams_*` 等 | [较低层计算](#lower-level) |

## 链接公共库

工程定义 `ASTER::core`，包含科学实现和 Eigen、SOFA、数学后端的传递依赖。自己的 CMake 目标链接它即可使用对应头文件及编译条件。

例如，在 `build/interface-demo/` 放一份 `CMakeLists.txt` 和下面的 `main.cpp`：

```cmake
cmake_minimum_required(VERSION 3.26)
project(orbit_demo LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(BUILD_TESTING OFF CACHE BOOL "Build ASTER tests" FORCE)
add_subdirectory("${ASTER_REPO}" aster-build)
add_executable(orbit_demo main.cpp)
target_link_libraries(orbit_demo PRIVATE ASTER::core)
set_target_properties(orbit_demo PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/aster-build/bin")
```

从 ASTER 根目录构建：

```sh
cmake -S build/interface-demo -B build/interface-demo/build \
    -DASTER_REPO="$PWD" -DCMAKE_BUILD_TYPE=Release
cmake --build build/interface-demo/build --target orbit_demo --parallel 3
build/interface-demo/build/aster-build/bin/orbit_demo
```

Boost 路径按[构建指南](02_quick-start.md)追加。示例将程序放在 OpenLibm 所在的 `aster-build/bin/`；自己的输出目录若不同，需把对应私有共享库放到程序旁边。常用头文件可通过 `#include "ASTER.h"` 引入；GL 求解器额外包含 `core_iod_solver.hpp`。

下文示例通过 `using namespace constants::Scale;` 引入归一化单位，通过 `using namespace constants::Conversion;` 引入换算系数，物理值和来源见[常数说明](11_constants.md)。

## 错误怎样返回

多数科学函数接收 `bastion::Overwatch& tracer`，通过 `tracer.Code()` 传递错误。正常值为 `bastion::ERR_NOERROR`；每次调用后立即检查，出错就停止使用本次结果。库函数通常返回调用者，主程序再决定如何保存日志或结束进程。

返回值还可能承载科学状态：BODN 的 `false` 可以表示迭代上限，GL 的负 flag 可以表示未收敛。先看 tracer 是否报错，再判断求解状态。JPL 和日历函数则使用自身的整数状态码，见对应小节。

<a id="states"></a>
## 状态对象与单位转换

声明：[orbit_statedef.hpp](../include/orbit_statedef.hpp)；实现：[orbit_statedef.cpp](../src/orbit_statedef.cpp)。

三个对象均保存六分量和归一化标志，不保存历元或参考系；调用方另外维护这些信息。构造或 `Set` 的最后一个参数 `true` 表示已经归一化，`false` 表示尚未归一化。单位由随后使用的长度、时间尺度配套决定。

| 接口 | 参数方向与单位 | 结果用途 |
| --- | --- | --- |
| `Keplerian(a,e,i,OMEGA,omega,M,is_normalized)` | 输入六根数；未归一化角度为 deg，归一化角度为 rad | 保存同历元根数，M 为平近点角 |
| `Cartesian(x,y,z,vx,vy,vz,is_normalized)` | 输入同系位置和速度 | 保存同历元 Cartesian |
| `NonSingular2(tracer,a,xi,eta,h,k,lambda,is_normalized)` | 输入 NSG2，只有 lambda 是角度 | 适用于当前 NSG2 计算链 |
| `kep.Normalize(tracer,uLength)` | 输入输出 kep；a 除以 uLength，角度 deg→rad | 长度原为 km 时用 UNIT_LENGTH_M/1000，原为 m 时用 UNIT_LENGTH_M |
| `cart.Normalize(tracer,uLength,uTime)` | 输入输出 cart；位置除以 uLength，速度乘 uTime/uLength | 时间尺度通常为 UNIT_TIME_S |
| `Restore(...)` | 输入输出对象，按相同尺度恢复长度、速度与角度 | 在输出副本上恢复物理单位 |
| `ToCartesian` / `ToKeplerian` / `ToNonSingular2` | 输入归一化对象，输出数组或对象 | 状态表示转换，输出仍归一化，参考系和历元由调用方保持 |
| `SyncAnomaly(tracer)` | 输入输出 Keplerian/NSG2 | 同步近点角等派生量；NSG2 转换和读取派生量前使用 |
| `ToArray(double*)` | 输出至少六个连续 double | 按[分量顺序](05_data-formats.md)导出，保留对象当前单位 |

下面是一份可编译运行的 `main.cpp`，使用预报算例中的初始根数。后面的小节给出可放在同一文件 `main` 之前的示例函数，按需在 `main` 中调用。

```cpp
#include "ASTER.h"
#include "core_iod_solver.hpp"
#include <iomanip>
#include <iostream>

using namespace constants::Scale;
using namespace constants::Conversion;

int main()
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    orbit::Keplerian kep(6800.0, 0.01, 45.0, 10.0, 10.0, 10.0, false);
    kep.Normalize(tracer, UNIT_LENGTH_M / 1000.0);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    orbit::Cartesian cart;
    kep.ToCartesian(tracer, cart);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    // 保留计算用归一化状态，在副本中恢复米和米每秒。
    orbit::Cartesian cart_si = cart;
    cart_si.Restore(tracer, UNIT_LENGTH_M, UNIT_TIME_S);
    if (tracer.Code() != bastion::ERR_NOERROR) return 1;

    double values[6];
    cart_si.ToArray(values);
    std::cout << std::setprecision(17);

    for (double value : values) std::cout << value << ' ';

    std::cout << '\n';
    return 0;
}
```

NSG2 定义为 xi=e·cos(omega+OMEGA)、eta=e·sin(omega+OMEGA)、h=sin(i/2)·cos(OMEGA)、k=sin(i/2)·sin(OMEGA)、lambda=M+omega+OMEGA。当前常用转换链使用椭圆轨道，瞬转平另有 a、i 等范围要求。低层近点角函数位于 [orbit_basic.hpp](../include/orbit_basic.hpp)，输入输出角度为 rad，`MeanAnomaly2EccentricAnomaly(M,e)` 等函数用于同一轨道的角度转换。

<a id="frames"></a>
## 时间、EOP 与测站坐标

声明：[time_transform.hpp](../include/time_transform.hpp)、[earth_eop.hpp](../include/earth_eop.hpp)、[frame_definition.hpp](../include/frame_definition.hpp)、[frame_earth.hpp](../include/frame_earth.hpp)。实现位于对应同名 `src/*.cpp`。

| 接口 | 输入 | 输出与状态 |
| --- | --- | --- |
| `Calendar2MJD(year,month,day,hour,minute,second,mjd)` | 同一时间尺度的日历各分量，秒可含小数 | 写入 MJD；返回 SOFA 日历状态，0 正常，负值表示日历错误；它仅改变日期表示 |
| `UTC2TT(mjd_utc)` / `TT2UTC(mjd_tt)` | 指定尺度的 MJD | 返回另一尺度的 MJD |
| `TT2TCG` / `TCG2TT`、`TDB2TCB` / `TCB2TDB` | 函数名所指尺度的 MJD | 返回目标尺度的 MJD |
| `Utc2Else(dut1,mjd_utc)` | UT1−UTC（s）、UTC MJD | 返回当前实现填写的 UTC、UT1、TAI、TT 字段；TCG/TDB/TCB 另行计算 |
| `EOP::Read_EOP(tracer,path)` | 输入路径为可修改字符串引用，读取现有 EOP 文本 | 保存原始 s、arcsec 记录，通过 tracer 报错 |
| `EOP::Get_EOP(tracer,mjd_utc,eop,fit)` | 输入 UTC MJD，通常 fit=1 插值 | 输出 aEOP；dut1/dat 为 s，xp/yp/dx00/dy00 为 rad；失败时保留传入输出 |
| `Geodetic(height,lambda,phi,false)` | 高度 m、经度 deg、大地纬度 deg | 构造站址；顺序与测站表不同 |
| `site.Normalize(tracer,UNIT_LENGTH_M)` / `site.ToGCRS(...)` | 输入输出站址归一化；输入 UTC MJD 和同历元 EOP | 输出归一化 GCRS 位置，无速度 |
| `GCRS_to_ITRS_Matrix(utc,eop,matrix)` | UTC MJD、s/rad 的 aEOP | 输出 3×3 位置旋转矩阵，`r_itrs=matrix*r_gcrs`；反向接口使用转置 |

EOP 插值要求连续日记录，非整日通常还要有下一日；超出覆盖时失败。站址类的 ITRS/GCRS 相互转换保持长度单位，椭球与归一化约定见[数据格式](05_data-formats.md)。

```cpp
void ExampleStation(bastion::Overwatch& tracer, double mjd_utc, frame::GCRS& station_norm)
{
    EOP table;
    std::string path = "data/eop/EOP-All_250516.txt";
    table.Read_EOP(tracer, path);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    aEOP eop;
    table.Get_EOP(tracer, mjd_utc, eop, 1);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    frame::Geodetic site(943.0, 86.174633, 41.725892, false);
    site.Normalize(tracer, UNIT_LENGTH_M);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    site.ToGCRS(tracer, mjd_utc, eop, station_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;
}
```

从仓库根目录调用，例如 `ExampleStation(tracer, 58406.8859375, station)`，输出位置除以 UNIT_LENGTH_M。GCRS/ITRS 矩阵仅旋转位置；完整速度变换要加入旋转矩阵的时间导数。站心地平矩阵 `Horizontal_to_GCRS_Matrix` 采用教材北—西—天轴，不含原点平移。

<a id="observations"></a>
## 原始观测与测站资料

声明：[model_observation.hpp](../include/model_observation.hpp)、[io_observation_arc.hpp](../include/io_observation_arc.hpp)；实现：[model_observation.cpp](../src/model_observation.cpp)、[io_observation_arc.cpp](../src/io_observation_arc.cpp)。

`observation::Observers::ReadFile` 和 `Observations::ReadFile` 接收 tracer 与文件路径。前者保存原始站表，后者保存类型 1—8 观测。角度 deg、位置/距离 m、距离率 m/s；完整读入成功后替换对象原有记录，失败时保留旧资料。

`GetRecords()` 返回只读 vector 引用。`AttachObservers(tracer,stations)` 按站号复制原始站址并设置关联标记，类型 6 使用自带平台位置。调用方在选定资料后检查匹配、目标号、时间顺序和算法支持的观测类型，再转换成计算单位。

```cpp
void ExampleReadObservations(bastion::Overwatch& tracer, observation::Observations& data)
{
    observation::Observers stations;
    stations.ReadFile(tracer, "tests/iod_gl/input/observers.txt");
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    data.ReadFile(tracer, "tests/iod_gl/input/observations.obs");
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    data.AttachObservers(tracer, stations);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    data.PrintSummary(std::cout);
}
```

`Observations::AddRecord(record)` 追加一条原单位记录，`WriteFile(tracer,path)` 按原始协议覆盖写出，保留资料顺序。`observation::WriteObservationArc` 接收非空 `vector<SimulatedObservationRecord>`、已存在结果目录、目标/站号、从 1 开始的弧号和两个真值输出开关；成功返回 OBS 路径，错误返回空字符串并设置 tracer。每条记录内的 OBS 为 UTC MJD、deg，kep 为 km、deg，cart 为 km、km/s；三类文件逐行对应。

可见性筛选、分弧和生成器专有输出位于[观测主程序](../apps/observation/main_generate_observations.cpp)，用法见[观测生成](06_observations.md)。

<a id="observation-model"></a>
## 理论观测和几何偏导

`determination::Compute_ObsAndPartials(tracer,obs,state,calc_obs,partials)` 接收同历元的计算观测与归一化 GCRS Cartesian，输出理论量和 3×6 瞬时几何偏导。它支持 `GroundRaDec`（1）、`SpaceRaDec`（6）和 `InertialPosition`（7）。

| 对象 | 方向与单位 |
| --- | --- |
| `obs` | 输入；历元 UTC MJD。类型 1/6 的 `value=[RA,Dec,0]` 为 rad，观测者 GCRS 位置为 m；类型 7 的 value 为目标 GCRS XYZ（m） |
| `obs.sigma` | 输入；东西/南北标准差 rad 或 XYZ 标准差 m，留给后续加权阶段使用 |
| `state` | 输入；归一化 GCRS `[x,y,z,vx,vy,vz]`，长度 UNIT_LENGTH_M、时间 UNIT_TIME_S |
| `calc_obs` | 输出；1/6 为 `[RA,Dec,0]`（rad），7 为 XYZ（m） |
| `partials` | 输出；列对应六个归一化 Cartesian 分量，速度列为 0；测角前两行是 `cos(Dec_obs)*dRA/dx`、`dDec/dx`，位置为 `dXYZ_m/dx_norm` |

```cpp
void ExampleObservationModel(bastion::Overwatch& tracer,
    const determination::Observation& obs, const orbit::Cartesian& state_norm)
{
    Eigen::Vector3d calculated;
    Eigen::Matrix<double, 3, 6> partials;
    determination::Compute_ObsAndPartials(tracer, obs, state_norm, calculated, partials);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    std::cout << calculated.transpose() << '\n';
}
```

理论 RA 本身保存在 `calculated[0]`，余弦因子已经进入第一行偏导。几何模型采用同历元目标减观测者的位置，不含光行时、光行差和折射；零距离或极点附近的测角偏导报错。单纯把视线转成距离/赤经/赤纬可用 `LOS_to_Equatorial`，其距离沿用输入单位，角度为 rad。

<a id="propagation"></a>
## 整段数值预报与轨迹插值

声明：[core_propagation.hpp](../include/core_propagation.hpp)、[model_propagation.hpp](../include/model_propagation.hpp)；实现：[core_propagation.cpp](../src/core_propagation.cpp)。

`OrbitPropagationJob` 分别保存 `opn_model`（力模型、卫星参数）、`opn_setting`（起点、时长、积分设置）、EOP/大气资料和两种归一化初态。GCRS 任务的 UTC、TT 起点须对应同一时刻；duration_norm 和 h0_norm 都以 UNIT_TIME_S 为尺度。力模型默认关闭摄动，保留中心二体。

| 高层接口 | 输入 | 输出及约定 |
| --- | --- | --- |
| `OrbitPropagation_NSG2(tracer,job,id,epochs,nsg2,cart)` | 准备好的任务、101/102/201/202 积分器编号 | 覆盖三个数组；UTC MJD、归一化 NSG2 与 Cartesian，含初态和完整接受节点 |
| `OrbitPropagation_RKF78_Cart/NSG2` | 任务，模式由 job 的 rkf 设置决定 | 向历元/状态 vector 追加记录；检查错误后使用有效结果 |
| `OrbitPropagation_Adams_Cart/NSG2` | 任务，Adams 节点及模式设置 | 同样追加完整接受节点；不足整步尾段由 RKF78 处理 |
| `Interpolate_CartesianOrbit(tracer,nodes,states,requests)` | 严格递增 UTC 节点与逐一配对状态；请求非递减且在覆盖内，允许重复 | 返回与请求一一对应的归一化 Cartesian；错误时空数组并设置 tracer |
| `determination::PropagateObservationEpochs(tracer,job,id,requests)` | 正向传播任务，覆盖请求历元 | 先整段积分，再调用公共插值，返回观测时刻状态 |
| `OrbitPropagation_RKF78_Cart_EBCRS` / `OrbitPropagation_Adams_Cart_EBCRS` | `OrbitPropagationJob_EBCRS`，JPL 已加载 | 追加 E-BCRS Cartesian 和 TDB MJD，状态仍采用 UNIT_LENGTH_M、UNIT_TIME_S |

下面从现成预报初态做一次 600 s 二体传播。示例函数输出完整节点；需要摄动力时，按[预报主程序](../apps/opn/earth/main_opn_nonsingular2.cpp)准备对应资料和卫星归一化参数。

```cpp
void ExamplePropagation(bastion::Overwatch& tracer, std::vector<double>& epochs_utc,
    std::vector<orbit::Cartesian>& cart_norm)
{
    OrbitPropagationJob job{};
    std::string path = "data/eop/EOP-All_250516.txt";
    job.opn_EOP.Read_EOP(tracer, path);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    job.initKep_norm.Set(6800.0, 0.01, 45.0, 10.0, 10.0, 10.0, false);
    job.initKep_norm.Normalize(tracer, UNIT_LENGTH_M / 1000.0);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    job.initKep_norm.SyncAnomaly(tracer);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    job.initKep_norm.ToCartesian(tracer, job.initCart_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    job.opn_setting.t0_MJD_UTC = 55255.99923398;
    job.opn_setting.t0_MJD_TT = UTC2TT(job.opn_setting.t0_MJD_UTC);
    job.opn_setting.duration_norm = 600.0 / UNIT_TIME_S;
    job.opn_setting.h0_norm = 10.0 / UNIT_TIME_S;
    job.opn_setting.rkf.maxTE_norm = 1e-12;
    job.opn_setting.display_progress = 0;
    std::vector<orbit::NonSingular2> nsg2_norm;
    OrbitPropagation_NSG2(tracer, job, 101, epochs_utc, nsg2_norm, cart_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 在已有轨迹上取中间时刻，无需再次积分。
    std::vector<double> request{epochs_utc.front() + 300.0 / 86400.0};
    auto sampled = Interpolate_CartesianOrbit(tracer, epochs_utc, cart_norm, request);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    std::cout << "nodes=" << epochs_utc.size() << ", samples=" << sampled.size() << '\n';
}
```

传播方向由 duration_norm 决定，零时长仅保留初态。高层预报允许反向传播，但公共插值要求输入节点严格递增。插值使用附近最多四个节点的位置与速度构造最高七次 Hermite 多项式，速度取同一位置多项式的导数，命中节点时直接返回已有状态。

`OrbitPropagationKep_2Body(tracer,t0,t1,kep0,kep1)` 另提供归一化 Kepler 二体推进，直接用两个 MJD 的差更新平近点角；时间尺度由调用方保持一致，函数内部不做 UTC→TT 转换。

<a id="iod"></a>
## GL 初轨求解

声明：[core_iod_solver.hpp](../include/core_iod_solver.hpp)、[core_iod.hpp](../include/core_iod.hpp)；实现：[core_iod_solver.cpp](../src/core_iod_solver.cpp)。

```text
iod_solver::InitialOrbitDetermination(tracer, job, iteration_output = nullptr)
```

tracer 和 job 为输入输出，第三项是可选 `std::ostream*`。函数写入 `job.solution`，包括参考历元轨道、O−C、RMS、J、flag 和逐轮历史。轨道/测站在 GCRS 中，位置除以 UNIT_LENGTH_M，速度以 UNIT_LENGTH_M/UNIT_TIME_S 归一化，测角与残差为 rad。

| 调用前的 job 字段 | 必须准备的内容 |
| --- | --- |
| `epochMJDUTC` | 待求状态的 UTC MJD |
| `obs.number`、`obs.type` | 观测条数 N，当前 GL 为 `OBS_TYPE_RADEC` |
| `obs.time_MJDUTC`、`ra_rad`、`dec_rad` | 长度均为 N，UTC MJD、rad、rad，属于同一目标 |
| `obs.tau` | `(time_MJDUTC−epochMJDUTC)*DAY_TO_NORMALIZED_TIME` |
| `obs.sigma_RaDec` | 东西、南北共同标准差，有限正数，rad |
| `observer.positions_gcrs_norm_emat` | 3×N，每列对应一条观测同历元测站位置/UNIT_LENGTH_M |
| `observer.position_gcrs_epoch_norm` | 第一条观测所属测站在参考历元的位置/UNIT_LENGTH_M，用于几何诊断 |
| `dynamicsModel`、`iterMethod`、`constraintType`、`constraintEcc` | 模型、迭代和约束，取值见[配置参考](04_configuration.md) |
| `initValue.type` | 默认线性初值取 `IOD_INITIAL_VALUE_DEFAULT`；指定初值取 `IOD_INITIAL_VALUE_SPECIFIED` 并准备 initCart |
| `iterCountMax`、`iterTOL`、`printIterInfo` | 最大轮数、归一化最大状态改正门限、屏幕开关 |

以下函数演示将已准备的观测数组交给 GL。调用前保证数组配对、站址参考系及单位正确；原始文件读取和站址准备可参照[GL 主程序](../apps/iod/main_iod_gl.cpp)。

```cpp
void ExampleGl(bastion::Overwatch& tracer, const Eigen::ArrayXd& epochs_utc,
    const Eigen::ArrayXd& ra_rad, const Eigen::ArrayXd& dec_rad,
    const Eigen::MatrixXd& stations_norm, const Eigen::Vector3d& station_at_epoch_norm,
    iod::IodJob& job)
{
    job = iod::IodJob{};
    job.epochMJDUTC = epochs_utc(0);
    job.obs.type = OBS_TYPE_RADEC;
    job.obs.number = static_cast<int>(epochs_utc.size());
    job.obs.time_MJDUTC = epochs_utc;
    job.obs.tau = (epochs_utc - job.epochMJDUTC) * DAY_TO_NORMALIZED_TIME;
    job.obs.ra_rad = ra_rad;
    job.obs.dec_rad = dec_rad;
    job.obs.sigma_RaDec = 3.0 * ARCSEC_TO_RAD;

    job.observer.positions_gcrs_norm_emat = stations_norm;
    job.observer.position_gcrs_epoch_norm = station_at_epoch_norm;

    job.dynamicsModel = IOD_DYNAMICS_J2_SOLAR_LUNAR;
    job.iterMethod = IOD_ITERATION_NEWTON;
    job.constraintType = IOD_CONSTRAINT_ORIENTATION;
    job.initValue.type = IOD_INITIAL_VALUE_DEFAULT;
    job.iterCountMax = 40;
    job.iterTOL = 1e-10;
    job.printIterInfo = false;

    iod_solver::InitialOrbitDetermination(tracer, job);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    std::cout << "flag=" << job.solution.flag << ", RMS(rad)=" << job.solution.RMS_angle << '\n';
}
```

正 flag 表示达到状态改正门限的轮数，−1 表示达到迭代上限，−2 表示数值失败。求解会重置本次结果和历史；Ljung–Box 的遗留 Q/p 当前为 −1。

GLML 的候选生成、边界处理、最小 J 筛选和专有输出保留在[主程序](../apps/iod/main_iod_glml.cpp)。需要完整 GLML 行为时使用该命令行入口；公共 GL 接口用于其中每个来源轨道的求解。

指定偏心率的来源用带约束的 GL 方程生成，迭代检查归一化状态改正量；最终 J 则在边界处理后按测角残差和共同 sigma 重算。生成阶段的方程残差、来源收敛状态与最终似然分别反映不同信息。外层继续允许有效的迭代上限候选按有限 J 参选，实验中同时保存所选轨道、J 和来源 flag。

<a id="batch"></a>
## 批处理定轨与后验历元转换

声明：[core_bodn.hpp](../include/core_bodn.hpp)；实现：[core_bodn.cpp](../src/core_bodn.cpp)。

```text
determination::BatchOrbitDetermination_NSG2(tracer, input, observations, settings, result)
determination::MovePosteriorToEpoch(tracer, input, settings, reference_epoch_utc, result)
```

| 参数 | 方向、单位和要求 |
| --- | --- |
| `input` | 输入 OrbitPropagationJob；initKep_norm/initCart_norm 已在首观测 UTC 历元，UTC/TT 起点一致，动力学及外部资料已准备 |
| `observations` | 输入 `vector<determination::Observation>`；单目标、UTC 非递减，类型 1/6/7，允许同历元多条；单位见[观测模型](#observation-model) |
| `settings` | 输入 BatchSetting；积分器、STM、差分步长、RMS/DRMS、剔除、协方差及可选 AMR 估计 |
| `result` | 输出 BatchResult，每次拟合先重置；保存状态、逐条残差和采用标记、完整积分节点、观测时刻轨道、历史与协方差 |

返回 `true` 表示满足停止条件；`false` 且 tracer 无错误表示迭代上限，保留最近结果；计算错误时返回 `false` 并设置 tracer。六状态或七参数问题均需采用分量数多于参数数且设计矩阵满列秩。

```cpp
bool ExampleBatch(bastion::Overwatch& tracer, const OrbitPropagationJob& prior_at_first_obs,
    const std::vector<determination::Observation>& observations,
    double output_epoch_utc, determination::BatchResult& result)
{
    determination::BatchSetting settings;
    settings.integrator_id = 101;
    settings.stm_method = 2;

    bool converged = determination::BatchOrbitDetermination_NSG2(
        tracer, prior_at_first_obs, observations, settings, result);
    if (tracer.Code() != bastion::ERR_NOERROR) return false;

    determination::MovePosteriorToEpoch(
        tracer, prior_at_first_obs, settings, output_epoch_utc, result);
    if (tracer.Code() != bastion::ERR_NOERROR) return false;

    return converged;
}
```

库返回的拟合历元最初是首观测历元。`MovePosteriorToEpoch` 把后验状态与协方差传到指定 UTC 历元，更新 `fit_epoch_utc`，保留观测弧段的节点、残差和采用集合。BODN 主程序使用它把结果传回输入先验历元。

`BatchResult::residuals` 的测角为东西/南北 rad，位置为 m；`covariance_nsg2` 采用归一化 NSG2，`covariance_cart_si` 采用 m、m/s 的乘积。联合协方差前六项为归一化 NSG2，可选末项为 AMR（m²/kg），含交叉项。协方差按已知观测 sigma 计算，输出中 `converged`、停止原因与形式不确定度分别读取。

<a id="averaging"></a>
## 瞬转平

声明：[analyt_averaging.hpp](../include/analyt_averaging.hpp)；实现：[analyt_averaging.cpp](../src/analyt_averaging.cpp)。

```text
orbit::OsculatingToMean(tracer, mjd_utc, osculating, options, mean)
```

输入为 UTC MJD、GCRS 密切 NSG2 和 `AveragingOptions`，输出同历元、同系、同尺度的均 NSG2。状态的 a 以 UNIT_LENGTH_M 归一化、lambda 为 rad，其余无量纲；输入要求椭圆、i<pi、a≤90000 km。`options` 的物理 AMR 使用 m²/kg，Cr 无量纲，和状态的归一化单位分开。

函数固定次数求解 `mean=osculating−short(mean)`，模型开关、展开阶数和月球模式决定均根数含义。失败时保持传入 mean 并设置 tracer。默认模型及 e≥0.5 时 JLM 的处理见[配置参考](04_configuration.md)。

```cpp
void ExampleAveraging(bastion::Overwatch& tracer, double mjd_utc,
    const orbit::Keplerian& kep_norm, orbit::NonSingular2& mean_norm)
{
    orbit::NonSingular2 osculating;
    kep_norm.ToNonSingular2(tracer, osculating);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    orbit::AveragingOptions options;
    options.srp_area_mass_m2_kg = 0.01;
    options.radiation_coefficient = 1.44;

    orbit::OsculatingToMean(tracer, mjd_utc, osculating, options, mean_norm);
    if (tracer.Code() != bastion::ERR_NOERROR) return;
}
```

<a id="ephemeris"></a>
## JPL 星历

声明：[model_ephemeris_jpl.hpp](../include/model_ephemeris_jpl.hpp)；实现：[model_ephemeris_jpl.cpp](../src/model_ephemeris_jpl.cpp)。

| 接口 | 输入 | 输出与状态 |
| --- | --- | --- |
| `Read_JPLEPH(path)` | 当前读取器兼容的 JPL 二进制文件 | 1 新加载，2 已加载同一路径，−1 失败或请求切换文件 |
| `Get_JPLEPH(mjd_tdb,target,center,state)` | TDB MJD、目标号、中心号、六个 double 的输出空间 | 返回 1 并写入目标相对中心的归一化位置/速度；未加载返回 −1 且保留输出 |
| `Get_Libration(mjd_tdb,Omega_p,i_s,Lambda)` | TDB MJD，星历已加载 | 输出月球天平动三个欧拉角，rad；未加载返回 −1 |

常用编号为地球 3、月球 10、太阳 11、太阳系质心 12。状态轴向沿用 JPL 惯性轴；`center=3` 取得地心相对状态，`center=12` 取得太阳系质心状态。在文件覆盖内调用，覆盖和校验见[科学数据](../data/README.md)。

```cpp
int ExampleSun(double mjd_tdb, double state_norm[6])
{
    int status = Read_JPLEPH("data/ephemeris/JPLEPH_C.430");
    if (status < 0) return status;

    return Get_JPLEPH(mjd_tdb, 11, 3, state_norm);
}
```

读取器使用进程内全局句柄，同一进程继续使用同一路径。长度尺度为 UNIT_LENGTH_M，速度尺度为 UNIT_LENGTH_M/UNIT_TIME_S；上例输出为太阳相对地球的六状态。

<a id="settings"></a>
## 配置读取

声明：[io_settings.hpp](../include/io_settings.hpp)；实现：[io_settings.cpp](../src/io_settings.cpp)。`SimpleInput` 保存关键词及文本，不自动解释文件路径或物理单位；路径基准由主程序设置。

| 接口 | 参数与结果 |
| --- | --- |
| `ReadSimpleInput(tracer,path,strict=false)` | 输入配置路径字符串，成功后替换记录；严格模式检查行格式、空键和重复键，失败时保留旧记录 |
| `GetString` / `GetDouble` / `GetInteger` | 返回值和输出引用两种形式；缺项为空字符串或 0，普通数值解析可读取文本开头的数值 |
| `GetDoubleStrict(tracer,key,value)` | 输出引用，要求完整、有限 double；失败保持原值并设置 tracer |
| `GetIntegerStrict(...)` | 进一步检查整数性和 int 范围；2、2.0、2e0 均可，2.5 报错 |
| `GetBoolean(...)` | 严格布尔读取，支持 true/1/on、false/0/off，区分大小写 |
| `HasKeyword(key)` | 查询可选项是否存在，不产生缺项警告；值为空也算存在 |
| `Count()` / `GetRecord(tracer,index,record)` | 记录数及单条键值，GetRecord 下标从 1 开始 |

```cpp
bool ExampleReadStep(bastion::Overwatch& tracer, std::string setting_file, double& step_s)
{
    io::SimpleInput input;
    input.ReadSimpleInput(tracer, setting_file, true);
    if (tracer.Code() != bastion::ERR_NOERROR) return false;

    input.GetDoubleStrict(tracer, "integrator_stepsize", step_s);
    if (tracer.Code() != bastion::ERR_NOERROR) return false;

    return true;
}
```

这里的 step_s 单位由所读关键词约定为 s。主程序还要检查正负、组合条件和资料覆盖，具体入口规则见[配置参考](04_configuration.md)。

<a id="lower-level"></a>
## 动力学、单步积分和矩阵接口

需要改变积分循环或研究偏导时，可以继续使用以下入口。常规整弧计算先使用上面的高层接口，统一准备任务和外部资料。

| 入口与声明 | 科学关系、单位和返回方式 |
| --- | --- |
| [Acceleration_Cart / Acceleration_NSG2](../include/accel_complete.hpp) | 输入 tracer、GCRS 任务、实际阶段 MJD(TT)、n=6、归一化状态 x；输出相对归一化时间的六维导数 f，Cartesian 为 `[vx,vy,vz,ax,ay,az]` |
| [Acceleration_Cart_EBCRS](../include/accel_complete.hpp) | 相同结构，时刻为 TDB、状态为 E-BCRS；预先加载 JPL 与 EOP，按开关加入 EIH 修正 |
| [RKF78_Cart / RKF78_NSG2 / RKF78_Cart_EBCRS](../include/integrator_rkf78.hpp) | 输入当前 TT/TDB MJD、任务、归一化非零试探步长，输入输出六状态；返回实际接受步长，可选返回归一化局部误差。失败时步长为 0、误差 NaN、状态保持原值并设置 tracer |
| [Adams_Cart / Adams_NSG2 / Adams_Cart_EBCRS](../include/integrator_adams.hpp) | 输入固定 TT/TDB 起点、归一化经过时间和固定步长，输入输出六状态与 AdamsHistory，输出实际步长及起步误差；每条轨道使用自己的历史，切换模型/资料后 Reset |
| [GeometryPartialsNSG2](../include/orbit_geometry.hpp) | 输入归一化 NSG2，返回 6×6 `dCartesian_norm/dNSG2_norm`，只含瞬时几何关系 |
| [AnalyticStmNSG2](../include/orbit_stm.hpp) | 输入归一化 NSG2 六向量、归一化经过时间和 AnalyticStmSetting，返回同尺度 6×6 STM；二体＋可选 J2 长期项/日月平均项，模型项独立选择 |
| [RTN_Matrix_from_Cart / RTN_Matrix_from_Kep](../include/orbit_geometry.hpp) | 输入同系 Cartesian 或 rad 根数，输出惯性分量到 R/T/N 的旋转矩阵；三行分别为 R、T、N 单位向量，T=N×R |
| [aster::math_backend::BuildDescription](../include/math_backend.hpp) | 无参数，返回构建后端与浮点配置的静态字符串，供版本和实验记录使用 |

GL 还在 `core_iod_solver.hpp` 暴露首轮、方程构造和残差函数，供研究者检查求解过程；完整任务从 `InitialOrbitDetermination` 进入。摄动力的分项公式位于 `accel_*.hpp`，解析短周期项位于 `analyt_*.hpp`，其调用组合由当前力模型或瞬转平接口负责。

<a id="cislunar"></a>
## 地月状态与单程测距

声明：[model_cislunar_observation.hpp](../include/model_cislunar_observation.hpp)，也由 `ASTER.h` 引入。状态使用绝对 TDB 日历历元和 E-BCRS 下的米制位置速度。调用时显式使用 `Cislunar::`。

| 对象或接口 | 输入与输出 |
| --- | --- |
| `Cislunar::StateRecord` | 年月日时分、带小数的秒、三维位置 m、三维速度 m/s、总钟差 s |
| `Cislunar::OneWayRangeRecord` | 真实收发状态、接收钟日历读数、几何距离、引力时延、钟差贡献和伪距；主程序另填编号和误差 |
| `InterpolateState(tracer,records,state)` | 输入递增状态表及 state 的 TDB 历元，写回该事件的位置、速度和钟差；覆盖范围包含首末端点 |
| `ComputeGravityDelay(tracer,transmit,receive,bodies,ids,delay)` | 输入两个事件的 BCRS 位置 m、接收时刻各天体位置及 JPL 编号；输出引力时延 s |
| `ComputeOneWayRange(tracer,receive,query_transmit,ids,record)` | 输入真实接收状态和按历元填写发射状态的函数；输出光行时解、发射事件、接收钟读数及理论伪距 |

`ComputeOneWayRange` 调用前先载入覆盖计算时段的 JPL 二进制历表。空引力天体列表关闭引力时延；命令行程序默认选择除地球外的十个主要天体。光行时按秒迭代，末次变化不超过 1e−12 s，最多二十轮。噪声与系统差由调用方在理论伪距上加入。

下面的示例用离散状态构成发射端查询函数。接收状态由调用方在真实接收时刻准备，包含位置和钟差；JPL 也由调用方预先加载。

```cpp
/**
 * @brief 根据接收真值和发射端状态序列计算一条理论单程测距。
 * @param[in,out] tracer 记录插值及光行时求解的错误。
 * @param[in] receive 真实 TDB 接收历元、E-BCRS 位置 m 和钟差 s。
 * @param[in] transmit_states 发射端按 TDB 递增的状态，位置 m、速度 m/s、钟差 s。
 * @param[out] record 输出收发事件、接收钟读数和各项理论伪距贡献。
 */
void ExampleOneWayRange(bastion::Overwatch& tracer, const Cislunar::StateRecord& receive,
    const std::vector<Cislunar::StateRecord>& transmit_states, Cislunar::OneWayRangeRecord& record)
{
    const auto query_transmit = [&transmit_states](bastion::Overwatch& trace, Cislunar::StateRecord& state)
    {
        Cislunar::InterpolateState(trace, transmit_states, state);
    };

    const std::vector<int> gravity_bodies{1, 2, 4, 5, 6, 7, 8, 9, 10, 11};
    Cislunar::ComputeOneWayRange(tracer, receive, query_transmit, gravity_bodies, record);
    if (tracer.Code() != bastion::ERR_NOERROR) return;
}
```
