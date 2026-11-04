# 程序使用

ASTER 的十一个主程序按任务选择。配置关键词见[配置参考](04_configuration.md)，文件列序和单位见[数据格式](05_data-formats.md)。以下命令从仓库根目录执行，程序已按[构建指南](02_quick-start.md)生成。

## 选择入口

| 任务 | 可执行程序 | 仓库内的配置或输入 | 主要结果 |
| --- | --- | --- | --- |
| 星历生成地基光学观测 | `aster_generate_observations` | [generate_observations/setting.txt](../tests/generate_observations/setting.txt) | 理论 OBS 弧段、可选轨道真值 |
| 两星状态生成单程测距 | `aster_generate_cislunar_observations` | [generate_cislunar_observations/aster.settings](../tests/generate_cislunar_observations/aster.settings) | 接收钟读数、伪距、计算分项与日志 |
| GL 初轨确定 | `aster_iod_gl` | [iod_gl/setting.txt](../tests/iod_gl/setting.txt) | 参考历元轨道、测角残差、报告 |
| GLML 初轨确定 | `aster_iod_glml` | [iod_glml/setting.txt](../tests/iod_glml/setting.txt) | 最小 J 轨道、候选表、残差、报告 |
| 地球 Cartesian 数值预报 | `aster_opn_cartesian` | [opn_cartesian/setting.txt](../tests/opn_cartesian/setting.txt) | Cartesian、Keplerian、NSG2 星历 |
| 地球 NSG2 数值预报 | `aster_opn_nonsingular2` | [opn_nonsingular2/setting.txt](../tests/opn_nonsingular2/setting.txt) | 同上，以 NSG2 积分 |
| 地月空间数值预报 | `aster_opn_cislunar` | [opn_cislunar/setting.txt](../tests/opn_cislunar/setting.txt) | E-BCRS/TDB Cartesian 星历 |
| 批处理定轨 | `aster_odn_batch` | [odn_batch/setting.txt](../tests/odn_batch/setting.txt) | 后验轨道、残差、形式协方差、可选 AMR |
| Kepler 根数转直角坐标 | `aster_convert_keplerian_to_cartesian` | [七列根数](../tests/convert_keplerian_to_cartesian/input.txt) | 七列 Cartesian，m、m/s |
| 直角坐标转 Kepler 根数 | `aster_convert_cartesian_to_keplerian` | [七列 Cartesian](../tests/convert_cartesian_to_keplerian/input.txt) | 七列根数，km、deg |
| 密切根数转均根数 | `aster_convert_osculating_to_mean` | [convert_osculating_to_mean/setting.txt](../tests/convert_osculating_to_mean/setting.txt) | 同历元七列均根数 |

每个程序支持 `--version` 和 `--help`。九个配置型程序通过 `-I setting.txt` 启动；两项直接转换通过位置参数接收文件名。配置内的输入和输出相对路径从 setting 所在目录解释，普通程序的父目录由运行前创建；单程测距程序会创建配置指定的输出目录。

## 先准备七套独立配置

以下命令为 GL、GLML、三套预报、BODN 和瞬转平准备运行目录。它使用 `tests/cases.json` 中已有的输入路径清单，把输入改为仓库资料的绝对路径，输出仍写在各自配置旁。运行结果会保存在 `build/` 中。

整段复制到终端执行一次：

```sh
python3 - <<'PYTHON'
import json
from pathlib import Path

root = Path.cwd()
example_root = root / 'build/examples/trial-001'
example_root.mkdir(parents=True, exist_ok=False)
cases = json.loads((root / 'tests/cases.json').read_text())
for name, info in cases.items():
    if info.get('conversion'):
        continue
    fixture = root / 'tests' / name
    run = example_root / name
    run.mkdir()
    lines = []
    for line in (fixture / 'setting.txt').read_text().splitlines():
        if '=' not in line or line.lstrip().startswith('#'):
            continue
        key, value = map(str.strip, line.split('=', 1))
        if value and key in info['input_keys']:
            value = str((fixture / value).resolve())
        lines.append(f'{key} = {value}')
    (run / 'setting.txt').write_text('\n'.join(lines) + '\n')
    print(run / 'setting.txt')
PYTHON
```

再次准备一批时，将 `trial-001` 改为 `trial-002`，并在后续命令中使用相同名称。直接运行 `tests/.../setting.txt` 会按原配置把相对输出写到测试目录，因此日常运行使用这里生成的副本。

## GL：从测角资料求初轨

```sh
build/v1.3.0/bin/aster_iod_gl -I build/examples/trial-001/iod_gl/setting.txt
```

输入是一份同目标的地基 RA/Dec 观测、一份测站表和 EOP。可选 `file_initial_orbit` 提供同参考历元的 GCRS 初值；未填写时使用默认线性初值。至少选用三条观测，`observations_to_use` 可指定前 N 条或 `All`。

算例目录中会生成 `gl_solution.txt`、`gl_residuals.txt` 和 `gl_report.txt`。先读轨道文件中的 `reference_epoch`、`state_units` 和 `solution_flag`，再看报告里的迭代与 RMS：正 flag 表示收敛轮次，−1 表示迭代上限，−2 表示数值失败。配置容差约束六个归一化状态改正量，轨道相对真值的精度另行比较。

## GLML：生成候选并按 J 选轨

```sh
OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE build/v1.3.0/bin/aster_iod_glml \
    -I build/examples/trial-001/iod_glml/setting.txt
```

输入与 GL 相同。程序先求基准 GL，再在 0≤e<0.8 内均匀指定 10,000 个偏心率生成候选，做边界处理并重新计算测角似然 J，从含基准解的 10,001 个候选中选择最小的有限 J。`random_seed` 控制失败回退时的随机抽样，省略时为 0，范围为 0—2147483647。

结果包括 `glml_solution.txt`、`glml_candidates.csv`、`glml_residuals.txt` 和报告。轨道文件的 `selected_candidate` 对应候选表的 `id`；候选 0 是基准解。`source_solution_flag` / `source_flag` 保存边界处理前的来源求解状态，当前规则允许达到迭代上限但满足其他有效性条件的有限 J 候选参与比较。应结合 J、残差和来源状态判断结果。

`iteration_maxround` 控制基准 GL 的轮数，指定偏心率候选在当前主程序中固定为 40 轮。两类设置与最终评分的关系见[配置参考](04_configuration.md)；历史诊断见[记录索引](validation/01_index.md)。

## 地球预报：Cartesian 或 NSG2

```sh
build/v1.3.0/bin/aster_opn_cartesian \
    -I build/examples/trial-001/opn_cartesian/setting.txt
build/v1.3.0/bin/aster_opn_nonsingular2 \
    -I build/examples/trial-001/opn_nonsingular2/setting.txt
```

两者使用 GCRS 状态，输入输出历元为 UTC，内部在 TT 中积分；积分状态分别采用 Cartesian 和第二类无奇点根数 NSG2。已有算例从同一初轨预报 600 s，并开启配置中的地球非球形、日月、月球非球形、阻力和光压项。

`reference_Keplerian` 与 `reference_Cartesian` 二选一，另一项留空。`propagation_length` 是带符号的秒数，正向、反向和零时长分别表示向前预报、向后预报和仅保留初态。输出文件为 `propagation_ephemeris_cartesian.txt`、`propagation_ephemeris_keplerian.txt`、`propagation_ephemeris_nsg2.txt`；可留空关闭某一种表示，但至少开启一个星历输出。

输出保存积分接受节点，步长随积分器而定。`integrator_maxTE` 比较归一化局部截断误差估计，其含义见[配置参考](04_configuration.md)。

## 地月预报：E-BCRS / TDB

```sh
build/v1.3.0/bin/aster_opn_cislunar \
    -I build/examples/trial-001/opn_cislunar/setting.txt
```

该程序以地球为原点、保留 BCRS 轴向，使用 TDB 时间。配置 `reference_epoch` 和输出星历首列都为 TDB，日历输入也按 TDB 解释。输入 Cartesian 为 km、km/s，输出单位由 `output_state_units` 决定。与地球入口交换状态前，应先处理参考系和时间尺度。

它使用 JPL 星历，支持配置中的行星第三体和后牛顿修正。当前地球非球形项沿用既有 GCRS 近似处理。结果为 `propagation_ephemeris_cartesian.txt`，任务报告记录实际力模型与时段。

## BODN：用观测改正已有轨道

```sh
build/v1.3.0/bin/aster_odn_batch -I build/examples/trial-001/odn_batch/setting.txt
```

BODN 接收同一目标的地基测角（类型 1）、天基测角（类型 6）和 GCRS 位置（类型 7），允许多测站及同历元多条记录。资料按 UTC 非递减排列。程序准备首观测历元的先验轨道，通过加权批处理估计六状态，按配置增加一个阻力或光压面积质量比 AMR。

已有算例使用类型 6 的两平台测角，输出 `bodn_residuals.txt` 和 `bodn_report.txt`。报告中的停止原因、WRMS、逐轮采用记录和形式标准差帮助判断拟合；RMS 或可比较的 DRMS 达标，并满足最少轮次等条件时停止。形式协方差使用提供的观测 sigma，不按拟合后的 WRMS 再放大。

### 接收 GL / GLML 轨道

在自己的 BODN 配置里填写 `file_apriori_orbit`，指向前一步的 `gl_solution.txt` 或 `glml_solution.txt`；同时将 `apriori_epoch`、`apriori_Keplerian`、`apriori_Cartesian` 留空。目标号必须一致。程序从轨道文件的 `state_units` 识别单位，并将初轨传到首观测历元；最终状态传回先验历元。

可通过 `file_orbit_solution_BODN` 保存后验轨道；三个 `file_aposteriori_eph_..._BODN` 字段分别控制星历。这里的轨道文件是 SimpleInput 文本，源码注释有时简称“SI”；它的单位由 `state_units` 指明。BODN 固定写 `state_units=2`，即 km、km/s、deg。

## 两项直接状态转换

```sh
mkdir -p build/conversion-demo
build/v1.3.0/bin/aster_convert_keplerian_to_cartesian \
    tests/convert_keplerian_to_cartesian/input.txt build/conversion-demo/cartesian.txt
build/v1.3.0/bin/aster_convert_cartesian_to_keplerian \
    build/conversion-demo/cartesian.txt build/conversion-demo/keplerian.txt m
```

根数输入为 UTC MJD、a（km）、e、i、OMEGA、omega、M（deg），适用于椭圆轨道；转出的 Cartesian 为 m、m/s。反向程序最后一项为输入长度单位 `m` 或 `km`，速度同时使用 m/s 或 km/s，输出根数为 km、deg。显式填写单位即可批量运行；省略时会进入交互选择。

每个非空行须恰好七个有限数值，允许尾随空白，额外列和非数值内容报错。转换保持原历元和坐标轴。

## 瞬转平：密切根数转均根数

```sh
build/v1.3.0/bin/aster_convert_osculating_to_mean \
    -I build/examples/trial-001/convert_osculating_to_mean/setting.txt
```

输入输出均为 UTC MJD 加六个 Kepler 根数（km、deg），末项为平近点角 M。程序逐历元删除所选模型的解析短周期项；均根数因此与 J2/J22、JL/JLM、日月和光压开关及展开阶数相关。配置、范围和条件见[配置参考](04_configuration.md)及[瞬转平接口](07_interfaces.md#averaging)。

算例输出 `o2m_ephemeris_mean_keplerian.txt`。输出和 `log_file` 已存在时会拒绝写入；再次运行请使用新的实验目录或输出文件名。

## 输出、报告和日志怎样保存

控制台用于进度和摘要，`file_report` 保存科学设置、计算过程和结果，`file_programlog` 保存原始配置、事件及失败调用链。瞬转平使用 `log_file`。按用途选择输出，路径相互分开。

| 程序 | 同名输出的处理 |
| --- | --- |
| GL / GLML、三套预报 | 科学结果和报告覆盖写入；`file_programlog` 已存在时拒绝 |
| BODN | 启用的输出按覆盖方式写入，包括事件日志 |
| 两项直接转换 | 覆盖输出文件 |
| 瞬转平 | 结果和日志拒绝覆盖 |
| 观测生成 / 加噪 | 生成器覆盖同名弧段；加噪要求新的结果目录 |

准备自己的资料时，以对应配置副本为起点，依次修改输入路径、目标与历元、模型参数、输出路径。保留实际运行的 setting，并记录 `--version`、Git 提交、线程数和随机种子。检查输出格式后再增加数据量。

<a id="cislunar"></a>
## 地月单程测距：两颗星、两个方向

`aster_generate_cislunar_observations` 读取导航星和用户星的离散状态，并按给定的真实接收事件生成单程观测。每条事件指定发射与接收编号，因此同一份配置可以处理导航星发射、用户星接收，也可以处理反方向。最小算例中每个方向各有一条观测。

先把仓库配置准备成独立运行副本：

```sh
python3 - <<'PYTHON'
from pathlib import Path
case = Path('tests/generate_cislunar_observations').resolve()
run = Path('build/examples/trial-001/cislunar').resolve()
run.mkdir(parents=True, exist_ok=False)
settings = dict(line.split('=', 1) for line in (case / 'aster.settings').read_text().splitlines() if '=' in line)
settings = {key.strip(): value.strip() for key, value in settings.items()}
for key in ['file_navigation_state', 'file_user_state', 'file_receive_events', 'file_JPLEPH']:
    settings[key] = str(case / settings[key])
settings['directory_output'] = str(run / 'output')
(run / 'setting.txt').write_text('\n'.join(f'{key} = {value}' for key, value in settings.items()) + '\n')
PYTHON
build/v1.3.0/bin/aster_generate_cislunar_observations \
    -I build/examples/trial-001/cislunar/setting.txt
```

状态历元采用 TDB 年月日时分秒，位置、速度采用 E-BCRS 下的 m、m/s，钟差为卫星钟读数减 TDB，单位 s。输出的 `observations.txt` 保存接收钟读数和伪距；`components.csv` 保存真实发收时刻、几何距离、引力时延、两端钟差、误差和迭代结果；`run.log` 保存运行与错误追踪。同名输出已存在时，重新指定输出目录再运行。

完整输入列序见[数据格式](05_data-formats.md#cislunar)，模型与计算流程见[观测生成](06_observations.md#cislunar)，库接口见[公共接口](07_interfaces.md#cislunar)。
