# ASTER 与 Orekit 独立环境模型的一天预报

比较双方为当前 ASTER `fix/earth-gm-conventions` 工作区与 Orekit 13.1。ASTER 的提交起点为 `4062942`，包含此后的常数规范化修改；具体计算源码与上一轮一致，92 个源码文件的 SHA-256 核对通过。

这次保留相同初态、引力常数和重力场系数，让 Orekit 使用自己的地球姿态、历表读取、大气及光压模型。一天最大位置差为 **3021.952 m**，主要来自大气模型的不同。上一轮在相同模型下得到的十微米级差异，核对了计算实现；本轮量化了模型选择对这条近地轨道的影响。

## 共同输入与不同模型

两边从 MJD(UTC) 55255.99923398 对应的同一物理初态出发，预报 86400 s，每 10 s 输出一次，共 8641 点。初始位置和速度见上一轮 `08_orekit-perturbations.md`。地球 GM、80 阶 GGM03C 系数及半径、日月 TT 相容 GM、卫星质量、面积、Cd 和 Cr 均相同。

| 项目 | ASTER 当前近地程序 | 本轮 Orekit |
| --- | --- | --- |
| 地球姿态 | SOFA `iauXys00b`，线性插值 EOP | IERS 2010，独立插值同源逐日 EOP；`simpleEOP=false`，启用短周期 EOP 修正 |
| 日月位置 | 自有读取器查询 DE430，以 TT 近似传入 TDB 参数 | `CelestialBodyFactory` 直接读取同一个 `JPLEPH_C.430`，自行按 TDB 查询 |
| 大气密度 | 表格确定初始轨道的参考密度、标高，加入高度和周日变化 | `HarrisPriester` 内置平均太阳活动密度表，默认余弦指数 n=4 |
| 大气速度 | 按固定角速度绕输入惯性系 z 轴共转 | 随 Orekit 地固系共转，由坐标变换取得惯性速度 |
| 光压 | 球形卫星模型，全程照明 | `SolarRadiationPressure`＋`IsotropicRadiationSingleCoefficient`，计入地球本影和半影 |
| 月球非球形 | 20 阶，卫星处直接项 | 同一系数、月球姿态与直接项表达式，由 Orekit 球谐函数计算；地心间接项另做单独试验 |

Harris–Priester 采用 Orekit 随附的密度表，适用高度为 100～1000 km。它是本轮选定的独立大气模型，输入不包含实测空间天气。光压按两边相同的质量 1 kg、面积 0.01 m²、Cr=1.44 计算，光压基准也相同。模型差主要进入密度、共转速度及遮挡因子。

接口 `Get_JPLEPH(const double MJD_TDB, ...)` 的参数定义确为 TDB；近地调用处也明确注释了以 TT 近似传入。这是当前近地程序采用的有意近似，本轮保留，不修改时间转换或接口。

## 完整预报与模型替换对照

完整比较：**ASTER 当前六项摄动全开 vs Orekit 采用上述独立环境、大气和光压模型**。

| 比较双方 | 全天最大位置差 | 一天末位置差 |
| --- | ---: | ---: |
| ASTER vs Orekit，Orekit 最大步长 7.5 s | 3021.952166 m | 3021.952166 m |
| ASTER vs Orekit，Orekit 最大步长 3.75 s | 3021.951515 m | 3021.951515 m |
| Orekit 7.5 s vs Orekit 3.75 s | 0.651 mm | 0.651 mm |

继续减半 Orekit 步长只改变约 0.65 mm，远小于此次模型替换产生的千米差。完整跨库比较的全天最大速度差为 3.348 m/s。

为定位来源，下面每行均比较 **上一轮与 ASTER 对齐的 Orekit 全开模型 vs 仅作指定替换的 Orekit 模型**。这几组都使用 Orekit 积分器，最大步长为 7.5 s。

| 单独替换的内容 | 全天最大位置变化 | 一天末位置变化 |
| --- | ---: | ---: |
| 地球姿态及独立 DE430 查询 | 8.538 mm | 8.537 mm |
| 带地影的光压 | 13.433 m | 5.476 m |
| Harris–Priester 大气及地固系共转速度 | 3016.824 m | 3016.824 m |
| 上述三种替换同时采用 | 3021.952 m | 3021.952 m |
| 仅补月球非球形的地心间接项 | 0.176 mm | 0.025 mm |

在同一条 ASTER 轨道上逐点计算，Harris–Priester 密度为 ASTER 密度的 0.691～1.250 倍；初始密度分别为 1.333e-11 和 1.110e-11 kg/m³。8641 个输出时刻中，Orekit 判断 3401 点处于本影、28 点处于半影。这些密度和遮挡差异在轨道积分中积累，解释了大气和光压替换后的主要变化。

独立查询 DE430 得到的太阳位置相对 ASTER 相差约 41.17～42.85 m，月球约 1.35～1.40 m；这包含两边查询时间处理的差别。它们与姿态处理合在上述“环境”替换组中，没有把该组的 8.538 mm 全部归因于 TT/TDB 近似。

这组比较没有真实轨道作为裁判。结论是当前实现的同模型计算能够对齐，而改用另一套大气及遮挡模型会显著改变预报；哪套模型更贴近实测轨道，需要结合密度资料、卫星参数或观测残差判断。

## 月球非球形为何要减去地心项

在牛顿阶近似下，地心相对位置等于卫星位置减地球位置，因此其二阶导数也要相减。月球非球形引力对地心相对运动的贡献为：

\[
\boldsymbol a_{\text{相对地心}}
=\boldsymbol a_{\text{月球非球形}}(\text{卫星处})
-\boldsymbol a_{\text{月球非球形}}(\text{地心处})。
\]

两项都由同一月球重力场求得，并转换到同一惯性轴向。第二项扣掉地心自身受到的加速度；它不是额外作用在卫星上的力。这一关系对 GCRS 的牛顿阶地心相对方程同样成立，GCRS 的名称不会使减法自动消失。ASTER 现有点质量第三体函数本来就采用这种相对加速度。

源码核对结果如下：

- ASTER `src/accel_complete.cpp:1079`：EBCRS 函数明确计算地心相对月心的位置，再以卫星处直接项减地心处间接项。
- ASTER `src/accel_complete.cpp:794`：近地 GCRS Cartesian 函数只计算卫星处直接项，随后旋回惯性系并加入右函数；近地 NonSingular2 路径也采用直接项。
- LHH `local/src/program1/EOM_EBCRS.FOR:789`：计算 `-XM` 对应的地心位置及月球非球形引力；第 800 行执行 `FNS_M=FNS_M-VEC2`。第二个右函数在第 1770 行也执行同样的减法。
- ASTER `src/accel_thirdbody.cpp:145` 和 Orekit `ThirdBodyAttraction`：点质量第三体已经同时计算卫星处作用与地心处作用。

本近地算例的五个采样时刻中，地心项模长为 4.29e-13～4.46e-13 m/s²；卫星处直接项为 4.13e-13～4.60e-13 m/s²；相减后为 1.62e-14～3.07e-14 m/s²。共同加速度大部分被扣掉，剩下空间变化造成的相对作用。

只在 Orekit 验证程序中补上这个减法，一天最大位置变化约 0.176 mm。数值很小，但若近地模型保留月球非球形摄动，方程中应采用其直接项与间接项之差。本轮没有修改 ASTER 动力学代码。

## 近地间接项修正后的复核

在上述对照完成后，已在 `Acceleration_Cart()` 和 `Acceleration_NSG2()` 中补齐地心间接项。两处均按现有 EBCRS 的表达式，计算卫星处与地心处的月球非球形加速度之差，再旋回惯性系并归一化。TT 近似查询保持原样。

本轮生产代码修改为 `src/accel_complete.cpp` 中的两段计算，以及 `include/accel_complete.hpp` 中对应的中文参数说明。声明与实现的说明完全一致；三个局部风格整理文件为 `src/frame_definition.cpp`、`src/model_atmosphere.cpp` 和 `src/accel_solar_radiation.cpp`，仅调整常数用途注释、空行和长公式换行，去除注释与空白后的代码标记未变。未增加生产函数或命名空间，继续使用现有 Tracer 入口和错误处理。

`tests/core/gravity_conventions.cpp` 增加了独立 J2 解析式检查：给定仅含 C20 的月球场，在近地和高轨状态分别核对卫星处减地心处的加速度，并检查 NonSingular2 的导数。Cartesian 最大加速度分量差为 8.31e-17 m/s²，NonSingular2 归一化导数差为 0，CTest 通过。

一天预报使用此前同一初态及六项摄动全开配置，共 8641 个输出点：

| 比较双方 | 全天最大位置差 | 一天末位置差 |
| --- | ---: | ---: |
| 修正后 ASTER Cartesian vs 修正前 ASTER Cartesian | 0.171 mm | 0.025 mm |
| 修正后 ASTER Cartesian vs 包含地心间接项的 Orekit 同模型结果 | 1.14 μm | 1.11 μm |
| 修正后 ASTER NonSingular2 vs 包含地心间接项的 Orekit 同模型结果 | 4.85 μm | 3.39 μm |
| 修正后 ASTER NonSingular2 vs 修正后 ASTER Cartesian | 5.01 μm | 4.48 μm |

另运行仅地球 80 阶、关闭月球非球形的控制组，输出与修改前逐字节一致。三个受影响目标编译成功，`git diff --check` 通过。修改保留在当前分支工作区，供用户验收。

本轮配置、输出、日志及单独 diff 位于 `月球间接项收尾_013`，与本页下方的验证目录同级。[修正验证记录](evidence/12_lunar-indirect-fix.json) 保存比较双方、结果、源码前后校验值和文件校验值；`changes.diff` 只包含此次收尾的 6 个代码文件。

## 结果文件与复现

完整目录：

```text
/Users/chendingcun/Projects/五院地月空间滤波项目/04_代码工作/实现观测量仿真/Orekit独立模型对照_012
```

`runs/independent.txt` 和 `independent_fine.txt` 是完整独立模型两档步长的轨道。其余四个文件分别保存环境、光压、大气、月球间接项替换结果。`environment.txt` 给出同一轨道处两种密度、照明比例、日月位置差；`lunar_indirect.txt` 给出卫星处、地心处及相减后的加速度模长。`comparison.json` 和 CSV 记录比较双方与逐点差异。

本轮直接读取上一轮 `aster_reference` 指向目录的 ASTER 全开轨道及共用资料，未重复运行未变化的 ASTER。新建的 `OrekitModelComparison.java` 调用 Orekit 自身模型，复用上一轮验证程序的数据准备及球谐读取。JNI 模块用于单项控制组，以及双方继续共用的月球姿态。六组 Orekit 预报均完成，Java 构建和运行日志保存在本目录。

进入本目录后可运行：

```sh
/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/javac -encoding UTF-8 -cp 'aster_reference:lib/*' OrekitModelComparison.java
/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/java -cp '.:aster_reference:lib/*' OrekitModelComparison "$PWD"
/Library/Java/JavaVirtualMachines/amazon-corretto-8.jdk/Contents/Home/bin/java -cp '.:aster_reference:lib/*' OrekitModelComparison "$PWD" refine
python3 compare.py
```

模型依据：[Orekit Harris–Priester](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/models/earth/atmosphere/HarrisPriester.java)、[太阳光压与地球遮挡](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/radiation/SolarRadiationPressure.java)、[地心相对第三体加速度](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/gravity/ThirdBodyAttraction.java)。
