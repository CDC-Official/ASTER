# ASTER 与 Orekit 的常数使用及轨道预报对照

比较双方为 **ASTER 当前 `fix/earth-gm-conventions` 工作区**与 **Orekit 13.1、Hipparchus 4.0.1**。ASTER 的提交起点为 `4062942`，包含之后的常数规范化、命名空间简化及注释修改；逐文件源码 SHA-256 见验证目录 `baseline.json`。本次没有改动 ASTER 计算源码。

## 常数及模型如何使用

| 计算 | Orekit 13.1 | 当前 ASTER |
| --- | --- | --- |
| 大地经纬高与直角坐标转换 | `OneAxisEllipsoid(a, f, bodyFrame)`，由调用者选择椭球；`transform()` 使用该对象的半径、扁率 | `Geodetic::ToITRS()`、`ITRS::ToGeodetic()` 使用 WGS84 半径和扁率 |
| WGS84 椭球参数 | 6378137 m，1/298.257223563 | 数值完全相同 |
| 球谐非球形加速度 | `HolmesFeatherstoneAttractionModel` 从重力场 provider 取得 GM、参考半径及 C/S；由地固 XYZ 计算球坐标 | `Earth_Nonspherical_Gravitational_Acceleration()` 使用 GGM03C 的配套参数和地固 XYZ |
| 地球姿态 | Frame 独立选择 IERS 约定及 EOP 处理 | `GCRS_to_ITRS_Matrix()` 调用 SOFA 并读取 EOP |

Orekit 的 `Constants` 也提供 IERS2010 椭球参数：6378136.6 m、1/298.25642；`OneAxisEllipsoid` 构造说明明确列出 WGS84、GRS80、IERS 等选择。因此，使用参考椭球需要先选定椭球，而非必然采用 WGS84。这里确认的是旧扁率与 Orekit 的 IERS 常数数值一致，不能据此确定用户旧常数的历史出处。

Orekit 的 `IERS2010_EARTH_MU` 为 3.986004418e14；本次预报明确传入与 ASTER 相同的 TT 地球 GM 3.986004415e14，没有直接采用该常数。地球形状的选择也没有改变球谐模型的 GM、半径和系数。

源码依据：[Constants](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/utils/Constants.java)、[OneAxisEllipsoid](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/bodies/OneAxisEllipsoid.java)、[HolmesFeatherstoneAttractionModel](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/gravity/HolmesFeatherstoneAttractionModel.java)、[ICGEMFormatReader](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/gravity/potential/ICGEMFormatReader.java)。

## 相同条件

- 初始历元：MJD(UTC) 55255.99923398；Orekit 以 ASTER `UTC2TT()` 得到的同一 TT 历元为积分起点。
- 初始位置（km）：5937.3763146517094，2716.356956252755，1644.074813304293。
- 初始速度（km/s）：−3.5084830876880702，4.594869895007859，5.1343051912384086。
- 地球 GM：3.986004415e14 m³/s²；球谐参考半径：6378136.3 m。
- 球谐系数直接导出 ASTER 的 GGM03C 表，保留 80 阶完全归一化 C/S；本次不作潮汐转换，Orekit provider 的潮汐元数据标为 UNKNOWN。四个固定 ECEF 点另查 10 阶和 80 阶。
- 预报 86400 s，每 10 s 对照一次，共 8641 点。纯二体和仅地球 80 阶重力分别运行。
- 没有日月、行星摄动、大气阻力、光压或其他附加力。两边读取同一份 EOP；本次没有使用 JPL 天体位置参与力计算。
- ASTER 主对照采用 RKF78、10 s 试探步长、归一化局部误差门限 1e-12。另以 5 s、1e-14 检查步长敏感性。
- Orekit 采用 Cartesian 状态和 Dormand–Prince 8(5,3)；收敛检查将步长上限由 60 s 依次减为 30、15、7.5 s。最终共享姿态结果使用位置绝对门限 1e-7 m、速度门限 1e-10 m/s、相对门限 1e-14、7.5 s 上限。两库的误差门限定义不同，以收敛检查判断误差量级。

主对照的时间轴按已确认的输出步序对应经过秒数，避免再次用输出的完整 MJD 相减来对齐历元。ASTER 5 s 算例最后一步另有累计时间舍入，见下节。

## 固定 ECEF 点的非球形加速度

双方输入相同 XYZ、GM、半径及球谐系数；都只输出非球形部分。坐标点为 (7000000,2000000,1000000)、(−4000000,5000000,−3000000)、(1000000,2000000,7000000)、(42000000,10000000,−2000000) m。

| 比较双方 | 阶次 | 最大分量差（m/s²） |
| --- | ---: | ---: |
| 当前 ASTER vs Orekit Holmes–Featherstone | 10 | 6.94e-18 |
| 当前 ASTER vs Orekit Holmes–Featherstone | 80 | 1.21e-17 |

这组结果核对的是两套球谐加速度实现，未经过天球/地固旋转。

## 一天轨道预报

| 比较双方及条件 | 全天最大位置差 | 一天末位置差 |
| --- | ---: | ---: |
| ASTER 10 s 二体 vs Orekit 数值二体（60 s 步长上限） | 7.32 μm | 7.32 μm |
| ASTER 10 s 二体 vs Orekit 二体解析解 | 2.41 μm | 2.32 μm |
| Orekit 数值二体（60 s 上限）vs Orekit 二体解析解 | 9.72 μm | 9.64 μm |
| ASTER 10 s 的 80 阶预报 vs Orekit 7.5 s 上限，共用 ASTER 姿态 | 6.74 μm | 6.74 μm |
| ASTER 10 s 的 80 阶预报 vs Orekit 15 s 上限，各用自己的姿态 | 5.997 mm | 5.997 mm |
| Orekit 15 s 上限，共用 ASTER 姿态 vs 改用自己的 IERS2010 姿态 | 5.992 mm | 5.992 mm |

共用姿态时，验证程序通过一个本地调用直接取得 ASTER 的 `GCRS_to_ITRS_Matrix()` 结果；没有采样插值旋转矩阵。轨道积分和球谐求和由 Orekit 独立完成。矩阵进入 Orekit Frame 后的方向核对通过，旋转夹角差为零。

独立姿态一侧为 Orekit `IERS_2010`、`simpleEOP=true`，使用同源逐日 EOP，关闭短周期 EOP 潮汐；ASTER 本算例调用 `iauXys00b` 并线性插值 EOP。两侧姿态角差在每 600 s 采样中最大为 0.522 mas。这里的约 6 mm 是姿态算法及 EOP 处理的合并影响，没有进一步拆成某一个章动项或插值项。

二体解析解用于评估积分误差；80 阶算例核对两库的数值一致性。当前证据表明，引力核与积分在统一模型下相符，独立姿态对照的毫米差主要由姿态处理差异引起。

## 收敛与末步时间舍入

Orekit 的 80 阶共享姿态预报在 60→30 s 上限时相差最多 2.56 mm，30→15 s 为 12.43 μm，15→7.5 s 为 1.89 μm。因此正式引力对照采用 7.5 s 上限结果，没有把初始设置下的毫米差归因于物理常数。

ASTER 5 s、1e-14 算例的二体末点相对解析解为 0.154 mm；80 阶末点相对 ASTER 10 s 算例为 0.166 mm。核对 `OrbitPropagation_RKF78_Cart()` 的时间累计和末步截短后，5 s 算例最后一步实际为 4.9999999782081685 s，累计步长对应的物理时长比 86400 s 短约 21.78 ns。以末速度乘这段时间，预测位置差为 0.165464 mm，与实际差向量相减后剩余约 0.085 μm。

该项属于累计时间与末步截短的浮点舍入，原始输出和诊断保存在 `step_time_roundoff.json`。本次只记录，没有修改积分调度或时间表示。

## 结果与复现

完整目录：

```text
/Users/chendingcun/Projects/五院地月空间滤波项目/04_代码工作/实现观测量仿真/Orekit对照_010
```

- `baseline.json`：ASTER 工作区源码校验值；`dependencies.json`：Orekit、Hipparchus 下载地址和 SHA-256。
- `aster_setting.txt`、`ggm03c80.txt`：共用初值、开关和球谐系数；EOP 路径在配置中。
- `runs/aster_*/setting.txt`、`cart.txt`：ASTER 的四组实际运行配置和输出。
- `runs/*.txt`：Orekit 数值及解析输出；所有状态采用 m、m/s，相对时间采用 s。
- `fixed_ecef.txt`、`orientation.txt`：加速度及姿态控制对照。
- `comparison.json`、各项 CSV：每组比较双方、逐时刻差异及汇总。
- `AsterNative.cpp`、`OrekitComparison.java`：验证源码；`math_backend_java.cpp` 仅把数学库查找路径改为当前 ASTER 的 OpenLibm 绝对路径，以适应 Java 进程。
- `native_build.log`、`orekit_run.log`、`orekit_refine.log`、`aster_run.log`：构建和运行记录，均成功完成。

Java 使用 Amazon Corretto 8；`run_aster.py` 运行 ASTER，Java 主程序不带第二个参数运行初始组，增加 `refine` 参数运行收敛组，`compare.py` 汇总所有结果。`build_aster_native.py` 保存本机原生验证库的完整构建参数。
