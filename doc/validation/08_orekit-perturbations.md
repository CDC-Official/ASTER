# ASTER 与 Orekit 的近地摄动逐项对照

比较双方为 **ASTER 当前 `fix/earth-gm-conventions` 工作区**与 **Orekit 13.1 / Hipparchus 4.0.1**。ASTER 的提交起点为 `4062942`，包含随后完成的常数规范化修改；具体源码由 `baseline.json` 中的逐文件 SHA-256 确定。验证期间计算源码未变。

本轮在相同初态、常数及环境输入下，分别打开现有近地预报程序的各项摄动，再运行全开组合。八组一天预报的全天最大位置差为 0.71～12.72 μm；全开组为 6.62 μm，Orekit 再减半步长后为 10.55 μm。这里比较的是同模型下的计算结果，实际预报精度还取决于力模型、环境资料和初始轨道的准确度。

## 输入与开关

沿用上一轮 Orekit 对照的近地轨道：初始 MJD(UTC) 55255.99923398，位置为 (5937.3763146517094, 2716.356956252755, 1644.074813304293) km，速度为 (−3.5084830876880702, 4.594869895007859, 5.1343051912384086) km/s。两边从 ASTER 转换出的同一 TT 历元出发，预报 86400 s，每 10 s 比较一次，共 8641 点。

| 参数或资料 | 本轮取值与使用方式 |
| --- | --- |
| 地球中心 GM、球谐 GM | GGM03C，3.986004415e14 m³/s² |
| 地球球谐参考半径和系数 | 6378136.3 m，GGM03C 完全归一化系数，80 阶；直接从同一张系数表导出 |
| 太阳、月球点质量 GM | DE440 常数转为 TT 相容值，分别为 1.3271244200652094e20、4.902800190601983e12 m³/s² |
| 日月位置 | ASTER 当前 `JPLEPH_C.430` 数值历表，共享同一查询结果；光压单项组使用当前程序的解析太阳 |
| 月球球谐 | `grail.txt`，20 阶；模型 GM=4.90279980693169e12 m³/s²，参考半径=1738000 m |
| 地球姿态、月球姿态 | 共用 ASTER 的地球旋转矩阵和 JPL 月球天平动矩阵 |
| 大气 | 同一张密度表、同一初始轨道所确定的参考密度与标高；WGS84 扁率修正参考密度面，加入相同的周日变化 |
| 卫星与光压 | 质量 1 kg，阻力及光压面积均 0.01 m²，Cd=2.2，Cr=1.44；1 AU 处光压 4.56e-6 N/m²，全程按照明计算 |

每个单项组都保留地球中心引力。月球非球形开关要求同时打开 JPL 月球点质量；“全开”指地球 80 阶、太阳第三体、月球第三体、月球 20 阶、光压及大气阻力这六项。

ASTER 使用 RKF78、10 s 试探步长、归一化误差门限 1e-12。Orekit 使用 Cartesian 状态、Dormand–Prince 8(5,3)，最大步长 7.5 s，位置和速度绝对门限分别为 1e-7 m、1e-10 m/s，相对门限为 1e-14。全开组另将 Orekit 最大步长减为 3.75 s。

ASTER 的完整 MJD 输出标签与 10 s 步序的偏差最多为 0.93 μs。本轮核对输出点数和步序后，用经过秒数对应状态，避免用完整 MJD 再相减造成时刻对齐误差。

## 一天轨道差

下表每一行均为 **ASTER 与 Orekit 7.5 s 最大步长结果的比较**。位置差是三维差向量的模长；最大值取全天 8641 个输出点。

| 开启的摄动 | 全天最大位置差（μm） | 一天末位置差（μm） | 全天最大速度差（m/s） |
| --- | ---: | ---: | ---: |
| 二体基线 | 12.72 | 12.72 | 1.41e-8 |
| 地球非球形 80 阶 | 6.74 | 6.74 | 7.54e-9 |
| 太阳第三体 | 9.71 | 9.61 | 1.08e-8 |
| 月球第三体 | 9.24 | 8.77 | 1.04e-8 |
| 光压 | 10.94 | 10.82 | 1.23e-8 |
| 大气阻力 | 0.71 | 0.45 | 7.71e-10 |
| 月球第三体＋月球非球形 20 阶 | 5.08 | 2.67 | 5.67e-9 |
| 六项全开 | 6.62 | 6.62 | 7.26e-9 |

全开组的 **Orekit 7.5 s 与 Orekit 3.75 s** 两次预报，全天最大位置差为 4.00 μm；**ASTER 与 Orekit 3.75 s** 的最大位置差为 10.55 μm。细化后仍在十微米量级，当前积分和浮点舍入对比较结果的影响也处于这一量级。

本轮统一了地球姿态。上一轮双方各用自己的地球姿态时，地球 80 阶预报相差约 6 mm；那组姿态算法和 EOP 处理差异另见 `07_orekit.md`，没有混入本表。

## 同一状态下的摄动力

取 ASTER 全开轨道在 0、21600、43200、86400 s 的位置和速度，传给两边的力计算函数。天体位置、姿态、常数和卫星参数均相同，直接比较输出加速度，借此检查差异是否来自力计算。

| 摄动项 | ASTER 与 Orekit 最大加速度分量差（m/s²） |
| --- | ---: |
| 太阳第三体 | 2.54e-18 |
| 月球第三体 | 7.94e-21 |
| 光压，解析太阳 | 5.29e-23 |
| 光压，JPL 太阳 | 2.65e-23 |
| 大气阻力 | 1.14e-19 |
| 月球非球形 20 阶 | 7.57e-28 |

两边按相同公式分别计算的大气密度，最大差为 2.84e-25 kg/m³。地球非球形的固定 ECEF 点验证已在上一轮完成，80 阶最大分量差为 1.21e-17 m/s²。

日月第三体调用 Orekit 的 `ThirdBodyAttraction`；地月球谐调用 `HolmesFeatherstoneAttractionModel`。阻力调用 `DragForce` 与 `IsotropicDrag`，密度函数由验证程序按 ASTER 当前公式配置。光压使用相同的无地影通量，再调用 `IsotropicRadiationSingleCoefficient` 计算卫星响应。这两项分别验证当前大气公式配合阻力计算，以及无地影光压计算。

## 当前模型中需要区分的设置

源程序有三处直接影响这次比较的具体行为，验证程序均按实际调用配置：

1. JPL 接口参数以 TDB 定义，近地右函数按调用处注释，以 TT 近似传入；月球天平动也使用 TT。本轮共享这些查询结果，并保留这一已知近似。
2. 仅开光压时，太阳位置来自解析历表；全开时来自 JPL。光压函数当前关闭地影，两组均按持续照明计算。
3. 月球非球形项计算卫星处的球谐加速度，尚未扣除地球处对应加速度。地心方程中这一间接项属于需要核定的模型问题，本轮没有改动。双方按当前同一表达式比较，验证了球谐求值及单位转换。

本轮没有发现所测力计算之间的数值不一致；上面列出的时间、地影和地心方程设置另行记录，不能与浮点积分差混为一谈。

## 文件与复现

完整验证目录：

```text
/Users/chendingcun/Projects/五院地月空间滤波项目/04_代码工作/实现观测量仿真/Orekit摄动开关对照_011
```

`cases.json` 保存八组实际开关；`runs/aster_*/setting.txt`、`cart.txt`、`console.txt` 分别为 ASTER 配置、轨道和日志。`runs/orekit_*.txt` 保存 Orekit 轨道，其中 `all_fine` 为 3.75 s 最大步长。`forces.txt` 保存逐分量加速度，`comparison.json` 与各组 CSV 保存汇总和逐时刻差异。

`AsterNative.cpp` 传递环境输入并调用 ASTER 力函数；`OrekitPerturbations.java` 执行 Orekit 对照。`math_backend_java.cpp` 将 OpenLibm 查找路径指定到本次构建目录，以便 Java 进程加载同一数学库。`baseline.json`、`dependencies.json` 和 `inputs.json` 保存源码、依赖及输入文件校验值。JNI、Java 编译和八组 ASTER、九组 Orekit 预报全部成功完成。

在项目工作目录运行 `python3 04_代码工作/实现观测量仿真/Orekit摄动开关对照_011/build_aster_native.py` 构建验证库；进入验证目录后执行：

```sh
python3 run_aster.py
/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/javac -encoding UTF-8 -cp 'lib/*' OrekitPerturbations.java
/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/java -cp '.:lib/*' OrekitPerturbations "$PWD"
/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/java -cp '.:lib/*' OrekitPerturbations "$PWD" 3.75
python3 compare.py
```

Orekit 源码依据：[第三体引力](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/gravity/ThirdBodyAttraction.java)、[球谐引力](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/gravity/HolmesFeatherstoneAttractionModel.java)、[大气阻力](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/drag/DragForce.java)、[球形卫星光压](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/radiation/IsotropicRadiationSingleCoefficient.java)。
