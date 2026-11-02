# 观测生成

`aster_generate_observations` 从目标星历和一个地面测站生成理论 RA/Dec，按可见性和输入 ID 连续性划分弧段，并输出同历元轨道真值。随后用 `tools/add_angular_noise.py` 加入随机误差。两部分来自 CDC_lib 实际工作区的 `apps/og/main_arcs.cpp` 和配套脚本，源码及原验收输入的 SHA-256 记录在 [source.json](../tests/generate_observations/source.json)。

## 运行原 LEO 算例

从 ASTER 根目录执行：

```sh
cmake -S . -B build/v1.2.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.2.0 --parallel 3
mkdir -p build/observation-demo
build/v1.2.0/bin/aster_generate_observations -I tests/generate_observations/setting.txt
python3 tools/add_angular_noise.py build/observation-demo/obs_*.txt \
    --sigma-arcsec 3 --samples 2 --seed 42 --output build/observation-noisy
```

配置相对路径从配置所在目录解释；`dir_output` 必须事先创建，生成器沿用覆盖同名弧段文件的行为。加噪的 `--output` 必须是新目录，样本写入 `sample_000001/` 等子目录。程序支持 `-I`、`-input`、`--input`、`--help` 和 `--version`；省略参数时显示用法。

原配置位于 [tests/generate_observations/setting.txt](../tests/generate_observations/setting.txt)。它保留原二体 LEO 星历的 8,641 条采样、测站、最低高度角 5°，默认关闭昼夜和地影筛选。该配置输出 6 段、364 条可见观测；星历是历史仿真算例，供功能和回归检查。

## 输入、筛选与输出

星历两种格式二选一，每个非空行恰好八列：

```text
ID MJD_UTC x y z vx vy vz
ID MJD_UTC a e i OMEGA omega M
```

目标状态采用地心 J2000 惯性轴作为 GCRS 近似。`is_input_normalized=0` 时使用 km、km/s、deg；为 1 时长度和时间尺度为 `ae`、`UNITTIME`，角度为 rad，历元始终为 UTC MJD。根数中的 M 是平近点角，要求 a>0、0≤e<1。输入时刻严格递增，程序逐采样使用星历，不重新积分或插值轨道。

测站表为 `站号 坐标类型 三个分量`。类型 1 为大地纬度(deg)、经度(deg)、高度(m)，类型 2 为 ITRS X、Y、Z(m)。`observer_id` 原样进入观测正文，数字站号仅在文件名补足四位。站址使用公共坐标转换；每个观测历元独立插值 EOP，无覆盖时失败。

当前只生成 `observation_type=1` 的地基几何 RA/Dec。最低高度角使用教材球形可见性公式；`use_station_night` 和 `use_earth_shadow` 分别启用测站黑夜和柱形地影筛选。启用任一项时，`solar_ephemeris=1` 使用随库兼容的 JPL 二进制与 TDB，`=2` 使用现有解析太阳历表与 UT1。夜间筛选还需 `max_sun_elevation_deg`（−90°≤值<0°）。太阳计算和现有坐标转换均通过公共数学后端构建。

数值关键词保留普通 SimpleInput 的读取规则：缺项或无法读取时警告并返回零，再由主程序检查业务范围。因此归一化和输出/筛选开关缺省为 0，最低高度角缺省为 0°；观测类型、启用筛选时的太阳模型和夜间阈值应显式填写。完整取值与条件项见配置注释。

文件共用 `目标号_站号_弧号_起始UTC_MJD.txt` 后缀：

| 文件 | 每行内容 |
| --- | --- |
| `obs_*.txt` | UTC MJD、类型、目标号、站号、RA(deg)、Dec(deg) |
| `kep_*.txt` | UTC MJD、a(km)、e、i、OMEGA、omega、M(deg) |
| `cart_*.txt` | UTC MJD、x、y、z(km)、vx、vy、vz(km/s) |
| `all_*.txt` | 与 OBS 相同的前六列，再加可见标记和输入 ID；包含不可见点 |

`need_output_Keplerian`、`need_output_Cartesian`、`need_output_all` 分别控制可选输出。弧段中的真值和观测逐行对应；可见 ID 中断时分弧，单点和双点弧也会保留。标准输出仅给出 OBS 文件数和路径，运行摘要写入标准错误。`all_` 是检查资料，应将 `obs_` 或加噪样本交给 GL/GLML。

## 随机误差与复现

加噪脚本仅用 Python 标准库。对每条理论观测，先独立抽取 ηE、ηN∼N(0,σ²)，再计算：

```text
RA = wrap(RA0 + ηE / cos(Dec0))
Dec = Dec0 + ηN
```

σ 由 `--sigma-arcsec` 指定，内部为 rad，输出为 deg；分母使用理论赤纬。`--seed` 必须为非负整数；脚本保留 `random.Random`/MT19937/`gauss` 及 SHA-256 派生种子的原实现，输入文件内容摘要和样本号决定各弧随机流。`noise.json` 记录模型、seed、源资料与脚本摘要和 Python 版本。移动文件、改变输入顺序、追加样本数都不改变已有样本。零噪声原样复制文件，真值保留在理论资料目录。

极点处该小角度模型无定义；正噪声在理论极点或加噪赤纬越极时明确失败。当前不含光行时、光行差和折射。2026-10-30 论文批次采用相同东西/南北误差公式，但使用 NumPy `SeedSequence/default_rng`、先整批抽取东西再抽取南北分量；本脚本保留 CDC_lib 的标准库随机流，相同 seed 不对应论文批次的相同样本。

## 检查

```sh
ctest --test-dir build/v1.2.0 -R '^observation_' --output-on-failure
python3 tests/check_observation_handoff.py --bin-dir build/v1.2.0/bin \
    --output build/v1.2.0/acceptance/observations
```

交接检查的输出目录须为新目录。它复用原 LEO 配置，两次生成理论值、以 seed 42 生成两组 3″ 样本，并将理论及首组带噪观测分别交给现有 GL、GLML 的十条观测配置。完整科学输出、日志和 `summary.json` 保存到指定构建目录。几何、输入约束和加噪检查复用 CDC_lib 原验收脚本。

## 本轮 Mac 核对结果

- 原运行检查 239 项、独立几何检查 62 项和原加噪检查通过；OPENLIBM、SYSTEM 两种构建均通过。原算例的可见记录数在无筛选、黑夜、地影和联合筛选下分别为 364、106、306、48，解析和 JPL 太阳历表在该算例中一致。
- SYSTEM 与 CDC_lib 原入口的 19 份观测、完整观测和真值文件逐字节一致。相同理论输入、3″、seed 42、两组样本下，原/新加噪脚本的 12 份文件逐字节一致。
- 默认 OPENLIBM 与旧 SYSTEM 相比，全部历元最大角差为 3.98×10⁻¹¹ deg，Cartesian 真值位置分量差不超过 1.01×10⁻¹¹ km、速度分量差不超过 1.07×10⁻¹⁴ km/s；历元、ID、可见标记、分弧和 Keplerian 真值相同。独立站心矢量反算与生成角度的最大差为 5.69×10⁻¹⁴ deg。
- 理论值及 seed 42 的带噪资料分别通过现有 GL、GLML 的十条观测入口，四次运行均写出有限轨道与十条残差；原算例单次生成约 0.1 s。结果用于功能和迁移核对。

加噪派生种子包含源文件内容摘要。因此切换数学后端造成的理论文件末位变化，也会改变派生随机流；逐样本复现需要相同输入字节、脚本、Python 环境及 seed。详细本地结果位于 `build/v1.2.0/acceptance/observations-final/summary.json` 和 `build/v1.2.0/acceptance/migration-vs-cdc-lib/migration-summary.json`。

## Linux 与跨平台核对

Linux x86_64 / GCC 13.3 的原观测运行、几何、加噪和坐标/SOFA 五项 CTest 全部通过，四次 GL/GLML 资料交接也通过。Mac arm64 / AppleClang 21 与 Linux 共用 OpenLibm 0.8.8、Eigen 3.4.0 和同一原 LEO 输入，输出比较如下：

| 输出 | 文件数 | Mac/Linux 逐字节一致 |
| --- | --- | --- |
| 理论 OBS、完整观测及两类真值 | 19 | 19 |
| 3″、seed 42、两组带噪样本 | 12 | 12 |
| 理论/带噪资料的 GL、GLML 轨道解和残差 | 8 | 8 |

本轮加噪实际使用 Mac Python 3.9.6、Linux Python 3.12.3；结论限于这组已比较的输入与环境。两端内部重复生成、反转加噪输入顺序均保持文件一致，四次求解的收敛标志相同；生成原算例的实际单次耗时约 0.12—0.15 s。

Linux 还以字节不变的原 `main_arcs.cpp`，通过四行头文件适配链接当前 ASTER 公共 core/math_backend，进行同后端的旧新入口受控比较，19 份输出逐字节一致。该对照只在构建目录内编译原主程序。Mac 使用旧 CDC_lib 原二进制的独立 SYSTEM 对照见上一节。

两端日常仓库均使用 `feature/observation-generation`，开发版本为 1.2.0；main 与 v1.1.0 仍在 `c18178a`。跨平台哈希、Linux 测试日志及旧入口摘要保存在 `build/v1.2.0/acceptance/cross-platform/`，其中 `comparison.json` 为逐类比较结果。
