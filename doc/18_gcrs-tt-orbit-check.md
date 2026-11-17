# GCRS/TT 轨道预报核对分支

本分支为双方程序核对提供一个独立入口 `aster_opn_gcrs_tt_check`。主程序按“读取 Cartesian 初态—准备 DE440 和 EOP—列出有效 GM—RKF78 积分—保存 Cartesian 结果”展开。分支名为 `codex/gcrs-tt-orbit-check`，基于 ASTER `c87a531ecdef2bba14e1cf8980dca08c8d024dd6` 建立。主程序源码见 [main_opn_gcrs_tt_check.cpp](../apps/opn/earth/main_opn_gcrs_tt_check.cpp)。

## Windows 下载、编译与运行

下载时选择 **codex/gcrs-tt-orbit-check** 分支，再使用 GitHub 的 Code → Download ZIP，或执行：

```powershell
git clone --branch codex/gcrs-tt-orbit-check --single-branch https://github.com/CDC-Official/ASTER.git
cd ASTER
```

准备 Visual Studio 2022 的“使用 C++ 的桌面开发”组件、CMake 3.26 或以上，以及 Boost 头文件。Eigen 3.4.0、SOFA、DE440 二进制历表、EOP 和月球重力场均已随库提供。以下命令在仓库根目录执行，Boost 路径改为本机含有 `boost/` 子目录的目录：

```powershell
cmake --preset gcrs-tt-check-windows -DASTER_BOOST_ROOT="C:/local/boost_1_87_0"
cmake --build --preset gcrs-tt-check-windows
.\build\client\bin\Release\aster_opn_gcrs_tt_check.exe -I tests/gcrs_tt_orbit_check/setting.txt
```

**本次核对统一使用系统数学后端 `SYSTEM`。** 上述专用 preset 已设置 `ASTER_MATH_BACKEND=SYSTEM`，使用 Windows 工具链的数学函数。请使用这个 preset；仓库通用的 `release` preset 默认选择 OpenLibm。测试入口若发现 OpenLibm 构建会报错并指出所需选项。

手动配置时的等价命令为：

```powershell
cmake -S . -B build/client -G "Visual Studio 17 2022" -A x64 -DASTER_MATH_BACKEND=SYSTEM -DBUILD_TESTING=OFF -DASTER_BOOST_ROOT="C:/local/boost_1_87_0"
cmake --build build/client --config Release --target aster_opn_gcrs_tt_check
```

macOS/Linux 使用 `cmake --preset gcrs-tt-check` 和 `cmake --build --preset gcrs-tt-check`，随后运行 `build/client/bin/aster_opn_gcrs_tt_check -I tests/gcrs_tt_orbit_check/setting.txt`。SYSTEM 的浮点实现随工具链变化，双方比较时保存各自 `--version` 或摘要中的构建信息。

## 默认算例与输出

| 项目 | 取值 |
| --- | --- |
| 初始历元 | UTC 2028-05-01 00:00:00，MJD 61892.0 |
| GCRS 位置，km | −172078.764589166，316878.877578190，136905.589039608 |
| GCRS 速度，km/s | −0.343284732289，−0.175709556631，−0.114628368645 |
| 积分时间 | TT，时长 864000 s |
| 积分器 | RKF78，固定步长 60 s；末段不足一步时补齐 |
| 地球引力 | 地球质点 + GGM03C 12×12 阶非球形 |
| 第三体 | DE440 太阳、月球质点 |
| 月球非球形、其他行星 | 默认关闭，可分别开启 |
| 光压、阻力、固体潮、1PN 项 | 关闭 |
| EOP | 资料覆盖内插值；缺测时定向参数置零 |

配置见 [setting.txt](../tests/gcrs_tt_orbit_check/setting.txt)，所有关键词均须给出，数值单位写在关键词或注释中。`earth_degree` 和 `moon_degree` 支持 0 或 2—80，同时指定带谐和田谐截断；0 表示关闭对应非球形项。开启 `moon_degree` 时须保留 `use_moon=true`。`use_planets=true` 加入水星、金星、火星、木星、土星、天王星、海王星和冥王星；火星至冥王星采用行星系统质心和系统 GM。

文件路径相对于配置文件所在目录。默认输出到 `build/gcrs-tt-check-output/`：

- `cartesian.csv`：含初态共 14401 行；列为 `elapsed_TT_s,MJD_UTC,MJD_TT,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s`。
- `summary.txt`：构建信息、模型开关、GM 表、地月质量参数比、初末状态与零 EOP 查询次数。
- `program.log`：原始配置、事件和失败调用链。

积分器的嵌入误差阈值为归一化状态的 `1e-12`，固定步长模式将其用于监测。`max_local_TE_normalized` 是局部误差估计；位置差应使用两套程序在同一历元输出的 Cartesian 坐标计算。运行使用的 RKF78 是现有 13 级实现，输出采用其八阶更新。

## 本分支如何使用 GM

以下数值单位均为 **km³/s²**。GM 是参与坐标运动方程的质量参数；TT、TDB 的“相容”表示配套的时间与空间尺度约定。

### 地球质点与非球形项

地球有效 GM 取 **398600.4415**，来自本库的 IERS 2010 TT 约定，与 GGM03C 的模型 GM 数值一致；球谐参考半径为 **6378.1363 km**。地球中心引力与非球形项共同采用这一 GM。

本程序以该地球 GM 和参考半径定义归一化单位。因此，内部地球 GM 写成 1，时间单位仍包含 `sqrt(R³/GM_E)`；恢复到 km、s 后，运动仍由上述有量纲 GM 决定。第三体使用 `GM_body/GM_E` 的比值。

### 日月及其他行星

日月和行星的基础数值来自 DE440/DE441 论文表 2 [3]，集中定义在 `constants::DE440::GM_TDB`。本入口使用：

\[
GM^{\mathrm{TT}}=\frac{1-L_G}{1-L_B}GM^{\mathrm{TDB}},\qquad
L_G=6.969290134\times10^{-10},\quad L_B=1.550519768\times10^{-8}.
\]

这里的上标表示质量参数的时空尺度约定。比例约为 `1.000000014808269`；对应的尺度讨论见 [1] 的 II.2、II.4 节及 [2]。

| 天体 | DE440 TDB 相容 GM | 本入口使用的 TT 相容 GM |
| --- | ---: | ---: |
| 月球 | 4902.800118 | 4902.800190601984 |
| 太阳 | 132712440041.279419 | 132712442006.52094 |

其他行星也统一使用该比例。运行摘要逐项打印 TDB 原值与 TT 换算值，便于直接与对方的常数表逐行核对；本分支直接使用源码中的 DE440 常数表，JPL 二进制文件提供位置、速度和月球天平动。

**实际地月 GM 比由实际采用的两项常数相除：**

- `GM_E / GM_M = 81.30056824752191`；
- `GM_M / GM_E = 0.012300037029943891`。

地球采用 GGM03C/IERS 的模型值，月球采用 DE440 值，因此实际比值与 DE440 文件头 `EMRAT=81.300568221497215` 略有差异。DE440 地球 GM 单独换算为 TT 后约为 `398600.4414095825`；摘要把它列为对照值，地球力模型实际使用前述 `398600.4415`。核对时须区分数据模型选择与 TT/TDB 尺度换算这两个步骤。

### 月球非球形项

打开 `moon_degree` 后，程序读取 `data/gravity/grail.txt` 的参考半径和无量纲归一化 C/S；文件头 GM 读取后保存在摘要中供核对，**有效 GM 随即覆盖为 DE440 月球的 TT 值 `4902.800190601984`**。这保证本入口的月球质点项与非球形项使用同一个有效月球 GM。

配套文件为 GRGM660PRIM 的 80 阶截取，文件头 GM 为 `4902.79980693169`，参考半径为 `1738 km`。本次约定保留原参考半径及 C/S，所以相对于文件原模型，球谐加速度幅度乘以“有效 GM / 文件 GM”。摘要在开启后打印两个 GM 和该幅度比。主库其他入口仍使用其已有的文件 GM 逻辑。

月球非球形摄动计算卫星处与地心处的球谐加速度之差，再转回积分轴系。GRGM660PRIM 发布时使用 DE421 月球主轴系，本算例采用 DE440 月球天平动旋转，并保留 C/S；双方开启该项时应共同采用这一姿态约定。默认算例将该项关闭。

## 时间、历表与 EOP 的具体约定

DE440 为仓库随附的原生二进制文件 `data/ephemeris/JPLEPH_C.440`，覆盖 JED 2433264.5—2506352.5；读取器支持这种原生格式。输入时检查编号必须为 440，且覆盖完整预报时段。

本入口沿用当前 Cartesian 右函数的约定：**以 TT 近似 JPL 查询参数 TDB，直接使用 JPL 的地心位置差矢量**。GM 采用上述 TT 尺度，动力学包含牛顿引力和所选球谐项。完整的相对论 BCRS—GCRS 坐标变换及 1PN 动力学见 [1]；本交付程序只实现这里列出的比较模型。双方核对时须对齐这两处历表近似以及力模型开关。

随库 EOP 文件截止到 2025-11-12，故本次 2028 年算例会将 `UT1−UTC`、`xp`、`yp`、`dX`、`dY` 全部置零；`UT1=UTC`，岁差、章动和地球自转模型照常计算。`file_eop` 留空也可显式采用整段零定向参数；给出的路径不存在或文件格式损坏时报告输入错误。

**TAI−UTC（DAT）按 SOFA 闰秒表处理。** 随库 SOFA 表对本算例给出 DAT=37 s，因此 TT−UTC=69.184 s。2028 年算例据此假定起始时刻前仍为该闰秒差；双方比较应使用同一假定。EOP 置零只改变地球定向参数。

## 参考文献

1. Turyshev, S. G., Bar-Sever, Y. E., & Bertiger, W. I. (2026). *Relativistic Framework for High-precision GNSS Processing in Geocentric and Barycentric Celestial Reference Systems with Extension to Cislunar Space*. **The Astronomical Journal, 172**(1), 58. [DOI: 10.3847/1538-3881/ae7436](https://doi.org/10.3847/1538-3881/ae7436)。[作者公开全文，arXiv:2511.12058v3](https://arxiv.org/html/2511.12058v3)。用于查阅 TT/TDB 尺度、GCRS/BCRS 动力学和地月空间扩展。
2. Klioner, S. A. (2008). *Relativistic scaling of astronomical quantities and the system of astronomical units*. **Astronomy & Astrophysics, 478**(3), 951–958. [DOI: 10.1051/0004-6361:20077786](https://doi.org/10.1051/0004-6361:20077786)。[正式期刊全文](https://www.aanda.org/articles/aa/pdf/2008/06/aa7786-07.pdf)。对应此前的 arXiv:astro-ph/0508292，本文引用正式发表版本；用于说明时间、空间坐标和质量参数的配套尺度变换。
3. Park, R. S., Folkner, W. M., Williams, J. G., & Boggs, D. H. (2021). *The JPL Planetary and Lunar Ephemerides DE440 and DE441*. **The Astronomical Journal, 161**, 105. [DOI: 10.3847/1538-3881/abd414](https://doi.org/10.3847/1538-3881/abd414)。表 2 提供本程序采用的 DE440 GM。
4. Petit, G., & Luzum, B. (eds.) (2010). *IERS Conventions (2010)*. IERS Technical Note 36. [IERS 官方全文](https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote36/tn36.pdf)。地球常数和重力场模型约定见第 1、6 章。

实际验证条件和数值结果见[本分支验证记录](validation/18_gcrs-tt-orbit-check.md)。
