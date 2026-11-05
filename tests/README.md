# 算例与检查

`tests/` 保存小输入、配置、科学参考和检查程序。自动检查的运行结果写入构建目录；手动地月最小案例按使用要求写入自身的 `output/`。只想使用程序时，先按[程序指南](../doc/03_programs.md)准备配置副本；这里说明检查各自比较什么。

## 日常选择检查

构建当前版本后，从仓库根目录执行对应命令。以下一组覆盖核心计算、SOFA、输入约束、数据校验、初轨到批处理接续，以及观测几何、加噪和地月测距：

```sh
ctest --test-dir build/v1.3.0 \
    -R '^(core_|math_backend_sofa|invalid_inputs|data_integrity|iod_batch_handoff|observation_|cislunar_|gravity_conventions)' \
    --output-on-failure
```

`ctest --preset release` 使用同一组当前回归。按改动范围选择正则即可，例如 `-R '^observation_'` 运行三项光学观测检查，`-R '^core_frame_coordinates$'` 只检查坐标。

| 检查 | 关注内容 |
| --- | --- |
| `core_frame_coordinates` | 测站、地固/天球位置、单位和观测几何 |
| `core_orbit_averaging` | 瞬转平公式、模型选项与状态转换 |
| `core_batch_determination` / `core_batch_parameters` | 批处理轨道、权重、剔除、协方差及附加 AMR |
| `core_rkf78_safety` / `core_adams_safety` | 实际阶段历元、反向/零时长、接受步长、起步与错误状态 |
| `core_iod_eccentricity` | 偏心率约束与独立参考表 |
| `core_earth_solid_tide` | 固定位置与参数的 NAO 二阶潮汐参考、独立势梯度、两极位置、永久部分及两套动力学接入 |
| `math_backend_sofa` | SOFA 本身的验证入口，与当前后端一起构建 |
| `invalid_inputs` | 九个原有主程序的格式、有限值、版本和写出失败路径；新增生成器的输入检查由 observation_runtime 覆盖 |
| `data_integrity` | `data/manifest.json` 中的大小与 SHA-256 |
| `iod_batch_handoff` | 已保存 GL/GLML 轨道文件与显式相同初态进入 BODN 后，残差一致 |
| `observation_runtime` / `observation_geometry` | 观测、真值和分弧，筛选、单位、EOP、输入约束和独立几何 |
| `observation_angular_noise` | 误差模型、固定种子、输入顺序、样本扩展、零噪声与失败情况 |

## 1.1.0 数学后端验收追溯

以下入口保存 1.1.0 数学后端切换时的参考与门限，使用相应版本构建执行。当前发布检查采用本文末尾的十一程序入口。四个单 GL 的历史检查使用已保存的观测和科学输出哈希，同时核对逐轮状态与停止分支：

```sh
python3 tests/math_backend_cases/run_checks.py --bin-dir build/v1.1.0/bin \
    --output build/v1.1.0/acceptance/gl --only-gl
```

需要完整候选验收时，去掉 `--only-gl`，执行四个单 GL 和 e010/e070 两个各含 10,001 候选的 GLML；`math_backend_iod` 这一 CTest 调用的就是该完整入口，只在 OPENLIBM 构建中注册。

七个普通程序使用 `cases.json` 中已有输入，按物理量门限比较：

```sh
python3 tests/check_math_backend_apps.py --bin-dir build/v1.1.0/bin \
    --output build/v1.1.0/acceptance/apps
```

门限为位置/半长轴 1 mm、速度 1 μm/s、偏心率与倾角参数 1e−10、角度 1e−7 deg、角残差 1e−10 rad；历元、编号和文本对应，角度按 360° 周期处理。可用 `--native-summary` 指定已有同平台原生运行清单；默认与仓库参考比较。范围、构建条件和历史结果见[数学后端指南](../doc/08_math-backend.md)及[1.1.0 验收](../doc/validation/03_math-backend.md)。

## 观测生成与定轨接续

`generate_observations/` 保存 CDC_lib 原 LEO 配置、8,641 条星历与测站，`source.json` 记录来源与 SHA-256。[观测指南](../doc/06_observations.md)给出直接生成和加噪命令。

```sh
ctest --test-dir build/v1.3.0 -R '^observation_' --output-on-failure
python3 tests/check_observation_handoff.py --bin-dir build/v1.3.0/bin \
    --output build/v1.3.0/acceptance/observation-handoff-new
```

接续检查的输出目录须为新目录。它两次生成理论值，以 seed 42 生成两组 3″ 样本，再把理论及首组带噪观测分别交给现有 GL、GLML 十条观测入口。科学结果、日志和 `summary.json` 一起保存。发布时的原/新入口与 Mac/Linux 对照见[1.2.0 验收记录](../doc/validation/04_observations.md)。

## 1.0.0 原始字节参考

九个原有程序各有一套迁移输入与参考，`cases.json` 记录入口、输入路径关键词和参考文件名。七套为配置型，另两项转换直接接收输入表，合计比较 16 份科学文件。参考包含迁移前已验证的 GLML 角度与随机种子修复。

`run_case.py` 在独立目录运行并逐字节比较原始参考，目录内有准备后的 setting 和 console。它用于 1.0.0 迁移追溯，例如：

```sh
python3 tests/run_case.py --case opn_nonsingular2 \
    --binary build/v1.3.0/bin/aster_opn_nonsingular2 \
    --output build/v1.3.0/legacy-reference/nonsingular2
```

命令执行后仍会比较旧字节；默认 OPENLIBM 与原生数学结果可能不同。无过滤的 `ctest` 包含这九个旧参考检查，当前后端的科学结果按上面的哈希与物理量比较入口判断。

`core/` 中的偏心率公式使用已保存独立参考，Adams 保留导师参考实现及来源说明。仅用于旧、新 IOD 求解器互比的测试留在旧工程。1.0.0 的迁移范围与当时的 19 项 CTest 结果见[迁移记录](../doc/validation/01_index.md#migration)。

## 地月单程观测最小案例

[`generate_cislunar_observations/`](generate_cislunar_observations/README.md) 集中保存 12、901 两颗星的仿真状态、双向两条真实接收事件和相对路径配置。运行所需 DE440 实体历表也保存在案例的 `input/` 中；输出保存到案例的 `output/`，生成文件不纳入 Git。构建、命令行和 CLion 运行方法见案例说明。

1.3.0 的当前回归还包括 `cislunar_interpolation`、`core_cislunar_observation`、`cislunar_observation_runtime` 和 `gravity_conventions`。光压参数估计采用连续受光的冬季低倾角高轨，按原门限检查轨道和面质比；圆柱地影的受光、阴影内外状态由 `gravity_conventions` 单独检查。

`observation_runtime` 的固定站址按 WGS84 配套，保留大地坐标与 ECEF 两种输入所得科学输出的逐字节比较。发布条件见[1.3.0 发布验收](../doc/validation/11_release-1.3.0.md)。

## 测站位置转换算例

[`convert_coordinates/`](convert_coordinates/README.md)保存九条经纬高输入，前三条来自 GL 的测站与观测时刻。`core_station_coordinates` 对照 SOFA 的 WGS84 与 IAU 2000B 组装结果，检查六向接口、往返和精确极点；`station_coordinates_runtime` 检查六向文件入口与内存结果一致，并覆盖坏输入、EOP 范围和错误追踪。两项检查已纳入当前回归预设。

```sh
ctest --test-dir build/release -R '^(core_station_coordinates|station_coordinates_runtime)$' --output-on-failure
```

## 坐标与轨道状态的最小往返

[`convert_states/`](convert_states/README.md)提供两条椭圆轨道。`conversion_roundtrips` 仅运行当前坐标和状态入口的往返转换，覆盖每个入口的三对表示及度米、度千米两套单位，结果保存在 `build/release/test-output/conversion_roundtrips/`。

## 地球固体潮算例

[`earth_solid_tide/`](earth_solid_tide/README.md)提供一颗低轨卫星的一天预报配置，复用公共 EOP、JPL 和原始 GGM03C。自动核心检查已纳入当前回归预设；`invalid_inputs` 同时检查两套地球入口的开关与资料错误。模型设置、与导师 NAO 的加速度及轨道比较见[固体潮验证](../doc/validation/12_earth-solid-tide.md)。

## 发布代表算例与跨平台核对

`release_smoke.py` 复用九个既有程序的最小输入，并运行地基测角、双向单程测距及六向位置坐标、轨道状态转换。它检查各程序的版本与帮助，保存科学文件的 SHA-256、行数和实际运行配置。

```sh
python3 tests/release_smoke.py --bin-dir build/v1.3.0/bin \
    --output build/v1.3.0/acceptance/local
```

输出目录须为新目录。两端使用同一源码、输入资料和 OpenLibm 构建，比较同名科学文件的 SHA-256；日志和含机器路径的报告随各次运行保存。

## 论文实验线 CPU 对照

`core_iod_reuse` 检查五种动力学的日月量复用、两种迭代方式的完整数值历史及报告重建。
`tests/thesis/run_comparison.py` 比较已发布 v1.4.0 与论文分支的三个完整 GLML 算例，保存实际输入、
逐次耗时和全部候选/轨道/残差校验；使用及验证范围见[论文实验记录](../doc/validation/15_thesis-cpu.md)。
