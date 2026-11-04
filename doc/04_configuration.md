# 配置参考

本页按主程序解释 `setting.txt`。先从[程序指南](03_programs.md)生成一份已有算例的配置，再修改与当前任务有关的字段。关键词区分大小写，例如 `use_SRP_perturb` 和 `use_srp_perturb` 是不同文本。

## 文件写法与路径

通常一行填写一个 `关键词 = 值`，向量分量以空格分开。注释单独成行，字符串路径直接填写，不加引号。例如：

```text
# 输入是一份地基赤经赤纬观测
file_observation = input/observations.obs
noise_sigma_arcsecond = 3.0
reference_epoch = First
```

命令行 `-I` 后的相对路径从启动工作目录解释；配置内部的输入、输出路径从配置所在目录解释。可选文件路径留空表示关闭相应输入或输出。程序按指定路径打开文件，运行前创建父目录。

除瞬转平外，配置使用普通 `SimpleInput` 规则：忽略空行、以 `#`、`!`、`*` 开头的注释行和没有等号的行；行尾文字保留在值中。普通数值读取在缺项或无法读取时警告并返回 0，字符串缺项返回空。因此，表中的“必填”表示成功使用该程序需要给出该值，现有算例是完整起点。

瞬转平启用严格配置：可用行尾 `#` 注释，拒绝重复键、空键、没有等号的非空行及未知关键词，数值须完整且有限。布尔值支持 `true/false`、`1/0`、`on/off`。其他程序的数值开关按各自表格使用 `0/1`。读取接口细节见[公共接口](07_interfaces.md#settings)。

## GL 与 GLML

依据：[GL 主程序](../apps/iod/main_iod_gl.cpp)、[GLML 主程序](../apps/iod/main_iod_glml.cpp)。两者均以 GCRS、UTC MJD 组织同一目标的地基 RA/Dec，文件角度为 deg，内部为 rad。

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `file_observation` | 必填，类型 1 的 OBS 文件 |
| `file_observer` | 必填，所选观测涉及的测站表 |
| `file_EOP` | 必填，覆盖参考历元的 EOP；当前 IOD 在各观测时刻转换站址，但共用参考历元 EOP |
| `observations_to_use` | 必填，`All` 或整数 N，取文件中前 N 条，N≥3 |
| `noise_sigma_arcsecond` | 必填，有限正数；东西、南北两个方向共同的标准差，arcsec |
| `reference_epoch` | 必填，`First`、`Mid`、`Last` 或 UTC MJD；`Mid` 取下标 N/2 的观测时刻（从 0 计），不是首尾平均时刻 |
| `dynamics_model_GL` | 必填，1 二体；2 J2；3 J2＋太阳；4 J2＋月球；5 J2＋日月，采用当前 f/g 模型 |
| `iteration_method_GL` | 必填，1 简单迭代；2 微分改正；带约束时使用 2 |
| `constraint_type_GL` | 0 无约束；1 方向；2 圆轨道；3 方向＋圆轨道；4 方向＋偏心率；缺项按 0；方法 1 只支持类型 0 |
| `constraint_eccentricity` | 类型 4 时读取，范围 0≤e<1；缺项为 0，其他约束下不使用 |
| `iteration_maxround` | 必填，正整数；GL 的迭代上限，GLML 中控制基准 GL |
| `iteration_state_tolerance` | 必填，有限正数；六个归一化 Cartesian 改正分量的最大绝对值门限 |
| `display_iterations` | 0 关闭、1 显示屏幕逐轮信息；缺项为 0，任务结果仍保留迭代记录 |
| `output_state_units` | 必填，1 归一化/rad；2 km、km/s、deg；3 m、m/s、deg |
| `file_orbit_solution_GL` | GL 必填，参考历元轨道解的输出路径 |
| `file_orbit_solution_GLML` | GLML 必填，最终选轨的输出路径 |
| `file_residual` | 可选，逐条 O-C 残差（rad） |
| `file_report` | 可选，设置、迭代和科学结果报告 |
| `file_programlog` | 可选，原始配置、事件和失败调用链；拒绝覆盖已有日志 |
| `file_initial_orbit` | 仅 GL 可选，SimpleInput 初轨；目标号和 UTC 历元须与本次任务相符，历元允许差 1e−10 日；省略或留空使用默认线性初值 |
| `file_candidates_GLML` | 仅 GLML 可选，候选 CSV；其状态固定使用 km、km/s、deg，不随输出单位码改变 |
| `random_seed` | 仅 GLML 可选，整数 0—2147483647，默认 0；控制失败回退抽样 |

GLML 的指定偏心率网格在主程序中固定为 `[0,0.8)` 的 10,000 点，另加一个基准解。每个指定偏心率来源使用方向＋偏心率约束、40 轮上限和当前状态容差。最终 J 在边界处理后统一重算，来源达到上限的有效候选继续参与最小 J 选优。来源 flag 和最终轨道的字段解释见[数据格式](05_data-formats.md)。

## 三套数值预报

依据：[Cartesian](../apps/opn/earth/main_opn_cartesian.cpp)、[NSG2](../apps/opn/earth/main_opn_nonsingular2.cpp)、[地月](../apps/opn/cislunar/main_opn_cislunar.cpp)主程序。

### 初始状态、时段与积分器

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `reference_epoch` | 必填，单个 JD/MJD 或 `年 月 日 时 分 秒`；地球两入口按 UTC，地月按 TDB；单数值大于 2400000.5 时按 JD 转 MJD |
| `reference_Keplerian` | 六项 `a e i OMEGA omega M`，km、deg；与 Cartesian 二选一，另一项留空 |
| `reference_Cartesian` | 六项 `x y z vx vy vz`，km、km/s；地球入口为 GCRS，地月为 E-BCRS |
| `propagation_length` | 必填，传播时长 s；正向、反向和零时长均支持 |
| `integrator_id` | 必填，101 自适应 RKF78；102 定步长 RKF78；201 Adams PE；202 Adams PECE |
| `integrator_stepsize` | 必填，s；自适应 RKF78 的初始试探步长，或定步长方法的步长大小；非零时长时不能为 0 |
| `integrator_maxTE` | 必填，有限正数；RKF78 六分量归一化局部误差绝对值之和的门限，Adams 起步与尾段也使用 |
| `integrator_adams_steps` | Adams 模式读取，整数 2—9；已有算例填 9 |
| `output_state_units` | 必填，1 归一化/rad；2 km、km/s、deg；3 m、m/s、deg；历元始终为 MJD |

自适应 RKF78 超限后减半步长重算；定步长 RKF78 保留步长并报告超限。Adams 使用 RKF78 起步，末尾不足整步时由 RKF78 处理。输出保存接受节点；要在指定历元取状态，可调用[轨迹插值接口](07_interfaces.md#propagation)。

### 资料、力模型与卫星参数

下列地球引力、日月、阻力、光压和卫星字段也用于 BODN。BODN 读取自己的配置，支持的月球项为质点引力；月球非球形三项属于预报入口。

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `file_EOP` | 必填，地球定向参数；应覆盖实际计算中的 UTC 时段 |
| `file_JPLEPH` | 地球入口在使用 JPL 日月、SRP 或月球非球形时需要；地月入口始终加载；兼容 JPL 二进制文件，非 SPICE BSP |
| `file_atmosphere_density` | 开启阻力时需要，大气密度表 |
| `earth_zonal_order` | 地球带谐最高 n，0—80；0 关闭这一部分 |
| `earth_tesseral_degree` | 地球田谐最高 n，0—80，含 m=1…n；0 关闭这一部分 |
| `solar_gravity` / `lunar_gravity` | 0 关闭，1 JPL，2 解析日月星历 |
| `use_lunar_NSG` | 预报入口的月球非球形开关，0/1；开启需 JPL 月球状态、天平动和月球重力场 |
| `lunar_zonal_order` / `lunar_tesseral_degree` | 开启月球非球形时使用的截断设置，按现有重力场表范围填写 |
| `file_lunar_gravity_field` | 月球非球形系数文件 |
| `use_drag_perturb` / `use_SRP_perturb` | 0 关闭，1 开启大气阻力/太阳光压 |
| `satellite_id` | 预报入口的卫星标识，用于模型显示 |
| `satellite_mass` | 质量，kg；启用阻力或光压时为正 |
| `satellite_drag_cross_area` / `satellite_srp_cross_area` | 迎风/受光面积，m²；AMR 由面积除以质量得到 |
| `satellite_drag_Cd` | 阻力系数 Cd，无量纲 |
| `satellite_srp_eta` | 光压反射参数 eta，无量纲；此字段沿用当前动力学实现，和瞬转平 `reflectivity` 分别设置 |

数值预报的光压采用地球圆柱地影：进入阴影后光压为零，离开后恢复。测角观测生成中的 `use_earth_shadow` 则用于筛选受光目标。月球非球形加速度采用卫星处与地心处加速度之差。

地月入口另读 `planetary_gravity`（0/1）、`planetary_body`（空格分隔的 JPL 行星编号，1、2、4、5、6、7、8、9）和 `use_post_newtonian`（0/1）。前者指定行星第三体，后者加入当前 EIH 后牛顿修正。

### 星历和显示

| 关键词 | 用途 |
| --- | --- |
| `file_eph_Cartesian` | 三入口均支持；地月入口必须填写 |
| `file_eph_Keplerian` / `file_eph_NonSingular2` | 两套地球入口支持；与 Cartesian 中至少开启一种 |
| `file_report` | 可选科学报告 |
| `file_programlog` | 可选事件日志，拒绝覆盖 |
| `display_model` / `display_setting` | 地球入口的模型/积分设置显示，0/1 |
| `display_summary` | 地月入口的摘要显示，0/1 |
| `display_progress` | 0 静默；1 约每 10% 显示进度；2 显示各步及试算信息；缺项按普通读取规则为 0 |

## BODN 批处理定轨

依据：[BODN 主程序](../apps/odn/main_odn_batch.cpp)、[BatchSetting 与 BatchResult](../include/core_bodn.hpp)。积分器、EOP、JPL、地球引力、日月、阻力、光压和卫星物理参数使用上节对应字段。

### 观测与先验

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `file_observation` | 必填，类型 1、6、7 可混合，单目标、UTC 非递减 |
| `file_observer` | 含类型 1 时必填；仅类型 6、7 时可留空 |
| `observations_to_use` | 必填，`All` 或前 N 条，N≥3；有效观测分量数还须多于待估参数数 |
| `noise_sigma_arcsecond` | 必填，测角东西/南北标准差，正数，arcsec |
| `noise_sigma_position_m` | 必填，GCRS XYZ 各分量标准差，正数，m |
| `apriori_epoch` | 内联先验的 UTC MJD；使用轨道文件时留空 |
| `apriori_Keplerian` / `apriori_Cartesian` | 内联六分量先验，二选一，km、deg 或 km、km/s；使用轨道文件时均留空 |
| `file_apriori_orbit` | 可选，GL、GLML 或 BODN 的 SimpleInput 轨道解；填写后内联先验三项均须留空 |
| `use_apriori_orbit_amr` | 可选，默认 0；1 从先验文件导入有效 AMR，同时更新面积；质量、Cd、eta 仍来自配置 |

输入先验可以早于或晚于观测弧段。主程序先把先验传到首观测历元用于拟合，再把后验及协方差传回输入先验历元；观测残差和星历保留在各自观测时刻。

### 估计、剔除与停止

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `stm_method` | 可选，默认 2；1 为 NSG2 解析近似，2 为完整传播中心差分 |
| `stm_use_J2` / `stm_use_lunisolar` | 方法 1 可选，默认均为 0；控制解析 STM 的 J2 长期项/日月平均项，独立于实际传播力模型 |
| `stm_fd_step` | 方法 2 可选，默认 1e−6；da=step·max(1,abs(a))，其他归一化 NSG2 分量扰动为 step |
| `estimate_parameter` | 可选，默认 0；0 六状态，1 阻力 AMR，2 光压 AMR，3 自动；同估 AMR 需相应力项开启且使用数值 STM |
| `parameter_fd_step` | AMR 数值偏导可选，默认 1e−3；对 AMR/先验 AMR 的比值施加中心差分 |
| `iteration_maxround` / `iteration_minround` | 必填，最大/最少轮数；已有算例 12/2 |
| `iteration_rms_tolerance` | 必填，标准化残差 WRMS 的门限，无量纲 |
| `iteration_drms_tolerance` | 必填，相邻可比较轮次 WRMS 的变化门限，无量纲 |
| `iteration_drms_relative` | 必填，0 为绝对变化；1 再除以本轮 WRMS |
| `rejection_scale` | 可选，默认 3；按同类型原始残差 RMS 的倍数剔除，0 关闭 |
| `rejection_start` | 剔除开启时可选，默认第 1 轮 |
| `rejection_angle_window_rad` | 剔除开启时可选，默认 1 rad；测角分量硬窗口 |
| `rejection_position_window_m` | 剔除开启时可选，默认一个归一化长度单位（6378136.3 m）；位置分量硬窗口 |
| `compute_covariance` | 可选，默认 1；0 关闭、1 计算已知 sigma 下的形式协方差 |
| `display_iterations` | 必填，0/1；仅控制屏幕，报告仍保留逐轮记录 |

自动 AMR 模式在不足三圈时只估计六状态。满三圈后，先验半长轴对应高度达到 2000 km 时选光压，低于该高度且已开启阻力时选阻力，否则仍只估计六状态；实际选中的 AMR 仍需相应力项开启。显式估计一种 AMR 时，质量与力系数固定。

WRMS 由采用分量的 `(O−C)/sigma` 计算。只有前后两轮采用集合相同时才比较 DRMS；RMS 或 DRMS 任一达标，并满足最少轮次和剔除启动要求后停止。剔除按观测类型跨站分组，使用剔除前原始残差 RMS，任一有效分量超限就剔除整条记录，每轮允许重新纳入。

### 输出

| 关键词 | 用途 |
| --- | --- |
| `file_orbit_solution_BODN` | 可选，先验历元的后验轨道，固定 `state_units=2`；可交给下一次 BODN |
| `file_aposteriori_eph_Cartesian_BODN` | 可选，观测弧段的 Cartesian 后验星历，km、km/s |
| `file_aposteriori_eph_Keplerian_BODN` | 可选，同历元 Kepler 后验星历，km、deg |
| `file_aposteriori_eph_NonSingular2_BODN` | 可选，同历元 NSG2 后验星历，km、deg |
| `ephemeris_interval` | 可选，默认 0，输出完整积分节点；正值须≥1e−5 s，按 UTC 秒等间隔取样并补末历元 |
| `file_residual` | 可选，含采用标记的逐条残差；测角为 rad，类型 7 位置残差除以 `UNIT_LENGTH_M` |
| `file_report` / `file_programlog` | 可选，科学报告/事件日志；现有入口均覆盖写入 |

## 瞬转平

依据：[主程序](../apps/conversion/main_convert_osculating_to_mean.cpp)、[AveragingOptions](../include/analyt_averaging.hpp)。使用 GCRS、UTC MJD 的七列密切 Kepler 根数，a>0、0≤e<1、0≤i<180°、a≤90000 km。输出为同历元均根数。

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `input_file` / `output_file` | 必填，七列输入/输出；结果拒绝覆盖 |
| `log_file` | 可选，模型、警告和任务状态；拒绝覆盖 |
| `gravity_degree` | 可选，整数 2—10，默认 10；JL/JLM 的最高次数 |
| `iteration_count` | 可选，整数 1—100，默认 3；固定点更新次数 |
| `hansen_expansion_order` | 可选，整数 0—4，默认 4；Hansen 偏心率展开及 q 截断 |
| `use_J2` / `use_J22` | 布尔值，默认均开启；分别扣除相应短周期项 |
| `use_frame` | 布尔值，默认开启；岁差章动与 J2 耦合项，需 J2 同时开启 |
| `use_JL` / `use_JLM` | 布尔值，默认均开启；l≥3 的带谐/田谐短周期项 |
| `use_lunar` / `use_solar` / `use_SRP` | 布尔值，默认均开启；月球、太阳、光压短周期项 |
| `lunar_short_period` | 默认 `single_fastarg`；`double_fastarg` 使用双快变量公式替代月球单快变量结果 |
| `SRP_AMR_m2_kg` / `SRP_AMR_normalized` | 可选且只能出现一个，非负；物理 m²/kg 或归一化 AMR；都省略时使用归一化默认 1468568181.9309881，约 0.01 m²/kg |
| `drag_AMR_normalized` | 可选，默认同上；仅写日志，不参与瞬转平 |
| `reflectivity` | 可选，Cr，范围 0—2，默认 1.44 |
| `display_records` | 布尔值，默认关闭；开启后在控制台逐条显示结果 |

启用 JLM 且输入 e≥0.5 时，当前调用关闭 JLM 并将最高次数设为 4。日月位置采用现有解析平位置，UT1 取 UTC；岁差章动使用 TT。有地影时省略光压短周期项。均根数的定义与上述模型选择一起保存，方便后续实验比较。

## 观测生成与加噪

[观测生成指南](06_observations.md)给出完整输入、筛选公式、输出和命令；[原配置](../tests/generate_observations/setting.txt)逐项注明条件。

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `object_id` | 必填，目标标识，用于观测和文件名 |
| `file_eph_Cartesian` / `file_eph_Keplerian` | 二选一，另一项留空；每行整数 ID、UTC MJD、六分量状态 |
| `is_input_normalized` | 0 为 km、km/s、deg；1 为归一化长度、时间及 rad；缺项为 0 |
| `file_EOP` | 必填，逐观测历元插值使用 |
| `file_observer` / `observer_id` | 必填，测站表与选用的原始站号 |
| `observation_type` | 必填，当前取 1，地基几何 RA/Dec |
| `min_elevation_deg` | 0≤值<90°，默认 0；LEO 算例填 5° |
| `use_station_night` / `use_earth_shadow` | 0/1，默认均为 0；测站黑夜/目标受光筛选 |
| `max_sun_elevation_deg` | 黑夜筛选开启时填写，−90°≤值<0° |
| `solar_ephemeris` | 任一光学筛选开启时填写；1 JPL/TDB，2 解析太阳/UT1 |
| `file_JPLEPH` | 筛选开启且太阳模型为 1 时填写 |
| `dir_output` | 必填，事先创建的结果目录 |
| `need_output_Keplerian` / `need_output_Cartesian` | 0/1，默认 0；逐行对应的轨道真值 |
| `need_output_all` | 0/1，默认 0；包含不可见点的完整检查表 |
| `file_programlog` | 可选，事件日志，拒绝覆盖 |

加噪通过 `tools/add_angular_noise.py` 单独运行，使用 `--sigma-arcsec`、`--samples`、`--seed` 和 `--output`；输出目录须为新目录。从旧生成配置迁移时，删除 `noise_sigma_as` 和 `random_seed`，改在加噪命令填写。加噪 seed 与 GLML 回退 seed 控制不同阶段，实验记录中分别保存。

## 地月单程测距仿真

依据：[单程测距主程序](../apps/observation/main_generate_cislunar_observations.cpp)。配置使用普通 SimpleInput 文件规则，数值字段通过严格读取接口检查。相对路径以配置目录为基准。

| 关键词 | 取值、单位与用途 |
| --- | --- |
| `file_navigation_state` | 必填，导航星的 TDB/E-BCRS 状态与总钟差文件 |
| `file_user_state` | 必填，用户星的同格式状态与总钟差文件 |
| `file_receive_events` | 必填，真实 TDB 接收时刻、发射编号和接收编号 |
| `file_JPLEPH` | 必填，覆盖事件及发射迭代时刻的 JPL 二进制历表；用于地球位置平移和引力时延 |
| `directory_output` | 必填，程序创建目录，并写出观测、分项与日志；已有同名结果时结束并报告路径 |
| `navigation_id` / `user_id` | 必填，两个不同的正整数，事件中的发收编号须为这两个编号 |
| `gravity_delay` | 可选，0 关闭、1 开启，默认 1；计入除地球外的十个主要天体，GM 采用 DE440 |
| `systematic_bias_m` | 可选，统一加入各条伪距的有限系统差，单位 m，默认 0 |
| `noise_sigma_m` | 可选，高斯白噪声标准差，单位 m，有限且非负，默认 0 |
| `noise_seed` | 可选，非负整数，默认 1；同一平台与构建下重复运行可复现噪声序列 |

两份状态表使用相同结构，每颗星至少两个严格递增的采样点。采样时段应覆盖接收事件及光行时求解访问的发射时刻。每条事件的方向由编号决定。输入列序和输出时刻的物理含义见[地月测距数据](05_data-formats.md#cislunar)。
