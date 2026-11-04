# 观测生成

`aster_generate_observations` 从目标星历和一个地面测站生成理论 RA/Dec，按可见性和输入 ID 连续性划分弧段，并输出同历元轨道真值。随后用 `tools/add_angular_noise.py` 加入随机误差，供 GL/GLML 使用；真值用于在相同历元评价定轨误差。来源及原验收输入的 SHA-256 见 [source.json](../tests/generate_observations/source.json)。

## 运行原 LEO 算例

从 ASTER 根目录执行：

```sh
cmake -S . -B build/v1.3.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.3.0 --parallel 3
mkdir -p build/observation-demo
build/v1.3.0/bin/aster_generate_observations -I tests/generate_observations/setting.txt
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

目标状态采用地心 J2000 惯性轴作为 GCRS 近似。`is_input_normalized=0` 时使用 km、km/s、deg；为 1 时长度和时间尺度为 `UNIT_LENGTH_M`、`UNIT_TIME_S`，角度为 rad，历元始终为 UTC MJD。根数中的 M 是平近点角，要求 a>0、0≤e<1。输入时刻严格递增，程序逐采样使用星历，不重新积分或插值轨道。

测站表为 `站号 坐标类型 三个分量`。类型 1 为大地纬度(deg)、经度(deg)、高度(m)，类型 2 为 ITRS X、Y、Z(m)。`observer_id` 原样进入观测正文，数字站号仅在文件名补足四位。站址使用公共坐标转换；每个观测历元独立插值 EOP，无覆盖时失败。

当前只生成 `observation_type=1` 的地基几何 RA/Dec。最低高度角使用教材球形可见性公式；`use_station_night` 和 `use_earth_shadow` 分别启用测站黑夜和柱形地影筛选。启用任一项时，`solar_ephemeris=1` 使用随库兼容的 JPL 二进制与 TDB，`=2` 使用现有解析太阳历表与 UT1。夜间筛选还需 `max_sun_elevation_deg`（−90°≤值<0°）。太阳计算和现有坐标转换均通过公共数学后端构建。

数值关键词保留普通 SimpleInput 的读取规则：缺项或无法读取时警告并返回零，再由主程序检查业务范围。因此归一化和输出/筛选开关缺省为 0，最低高度角缺省为 0°；观测类型、启用筛选时的太阳模型和夜间阈值应显式填写。完整取值与条件项见配置注释。

`obs_`、`kep_`、`cart_` 共用 `目标号_站号_弧号_起始UTC_MJD.txt` 后缀，弧号为五位，MJD 保留十位小数。`all_` 保存整份输入，其文件名为 `all_目标号_站号_首条UTC_MJD.txt`。

| 文件 | 每行内容 |
| --- | --- |
| `obs_*.txt` | UTC MJD、类型、目标号、站号、RA(deg)、Dec(deg) |
| `kep_*.txt` | UTC MJD、a(km)、e、i、OMEGA、omega、M(deg) |
| `cart_*.txt` | UTC MJD、x、y、z(km)、vx、vy、vz(km/s) |
| `all_*.txt` | 与 OBS 相同的前六列，再加可见标记和输入 ID；包含不可见点 |

`need_output_Keplerian`、`need_output_Cartesian`、`need_output_all` 分别控制可选输出。弧段中的真值和观测逐行对应；可见 ID 中断时分弧，单点和双点弧也会保留。GL/GLML 使用前应选择至少三条观测的弧段。标准输出仅给出 OBS 文件数和路径，运行摘要写入标准错误。`all_` 是检查资料，定轨读取 `obs_` 或加噪样本。

## 随机误差与复现

加噪脚本仅用 Python 标准库。对每条理论观测，先独立抽取 ηE、ηN∼N(0,σ²)，再计算：

```text
RA = wrap(RA0 + ηE / cos(Dec0))
Dec = Dec0 + ηN
```

σ 由 `--sigma-arcsec` 指定，内部为 rad，输出为 deg；分母使用理论赤纬。`--seed` 必须为非负整数；脚本保留 `random.Random`/MT19937/`gauss` 及 SHA-256 派生种子的原实现，输入文件内容摘要和样本号决定各弧随机流。`noise.json` 记录模型、seed、源资料与脚本摘要和 Python 版本。移动文件、改变输入顺序、追加样本数都不改变已有样本。零噪声原样复制文件，真值保留在理论资料目录。

从旧配置迁移时，删除理论生成 setting 中的 `noise_sigma_as` 和 `random_seed`，把噪声大小与种子填写到加噪命令。当前生成器会对这两个遗留字段报错。加噪种子控制观测误差，GLML 配置的 `random_seed` 控制求解回退，实验中分别记录。

极点处该小角度模型无定义；正噪声在理论极点或加噪赤纬越极时明确失败。当前不含光行时、光行差和折射。2026-10-30 论文批次采用相同东西/南北误差公式，但使用 NumPy `SeedSequence/default_rng`、先整批抽取东西再抽取南北分量；本脚本保留 CDC_lib 的标准库随机流，相同 seed 不对应论文批次的相同样本。

## 检查

```sh
ctest --test-dir build/v1.3.0 -R '^observation_' --output-on-failure
python3 tests/check_observation_handoff.py --bin-dir build/v1.3.0/bin \
    --output build/v1.3.0/acceptance/observations
```

交接检查的输出目录须为新目录。它复用原 LEO 配置，两次生成理论值、以 seed 42 生成两组 3″ 样本，并将理论及首组带噪观测分别交给现有 GL、GLML 的十条观测配置。完整科学输出、日志和 `summary.json` 保存到指定构建目录。几何、输入约束和加噪检查复用 CDC_lib 原验收脚本。

## 发布验收记录

1.2.0 的迁移、独立几何和跨平台结果集中保存在[观测生成验收记录](validation/04_observations.md)。日常检查入口见[测试说明](../tests/README.md)。

<a id="cislunar"></a>
## 地月单程测距仿真

地月测距使用独立程序 `aster_generate_cislunar_observations`。最小计算单元是一条真实接收事件：已知接收端的时刻、位置和总钟差，并提供能够按时刻查询发射端状态的序列或函数。程序迭代求发射时刻和光行时，计算传播距离、引力时延及两端钟差贡献，输出接收钟读数和伪距。

导航星发射与用户星发射使用同一套计算，只需交换两端的职责。命令行程序读取两份同格式状态表和一份接收事件表，核心科学函数接收结构体与状态查询函数。文件读取、误差生成和输出写出由主程序组织。

| 资料或结果 | 含义 |
| --- | --- |
| 两颗星的状态序列 | TDB 年月日时分秒，E-BCRS 位置和速度，以及卫星钟读数减 TDB 的总钟差 |
| 接收事件表 | 真实 TDB 接收时刻和收发对象编号 |
| JPL 历表 | 在各自事件时刻查询地球 BCRS 位置，并准备引力时延所需的天体位置 |
| `observations.txt` | 接收钟读数、收发编号和最终伪距 |
| `components.csv` | 真实收发事件及几何距离、引力时延、钟差和误差等计算分项 |

两端的 E-BCRS 位置分别加上各自时刻的地球 BCRS 位置，从而在同一 BCRS 轴向下计算两个事件之间的距离。状态表用附近最多十点作 Neville 插值，位置、速度和总钟差分别求值。主程序可按开关计入除地球外十个主要天体的引力时延，并对最终伪距加入系统差和高斯白噪声。

[程序使用](03_programs.md#地月单程测距两颗星两个方向)提供独立运行副本的准备命令；[配置参考](04_configuration.md#地月单程测距仿真)列出全部关键词；[数据格式](05_data-formats.md#cislunar)解释输入输出历元。核心接口见 [Cislunar 科学计算](07_interfaces.md#cislunar)。

```sh
ctest --test-dir build/v1.3.0 \
    -R '^(cislunar_interpolation|core_cislunar_observation|cislunar_observation_runtime)$' \
    --output-on-failure
```

这些检查覆盖解析插值、跨日事件、钟差符号、引力时延、正反两向最小观测、系统差、噪声复现和 Tracer 失败追踪。1.3.0 的整体结果见[发布验收](validation/11_release-1.3.0.md)。
