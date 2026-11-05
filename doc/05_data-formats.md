# 数据格式与计算约定

本页解释文件各列与计算中的单位，供准备输入、接续程序和读取结果时查阅。字段以当前主程序和[公共头文件](../include/ASTER.h)为依据；配置关键词另见[配置参考](04_configuration.md)。

## 时间、参考系和归一化

### 历元属于哪一种时间

MJD 是以日为单位的历元数值，`MJD = JD − 2400000.5`。UTC、TT、TDB 是时间尺度，使用相同的 MJD 写法并不表示它们数值相同。

| 使用位置 | 文件或接口历元 | 内部计算 |
| --- | --- | --- |
| GL / GLML | UTC MJD | 观测相对参考历元的 UTC 日差，再除以 `UNIT_TIME_S` |
| 地球 Cartesian / NSG2 预报 | 输入输出 UTC MJD | TT 积分 |
| BODN | 输入输出 UTC MJD | TT 积分，在观测时刻插值 |
| 地月预报 | 输入输出 TDB MJD，日历初始时刻也按 TDB | E-BCRS/TDB 积分 |
| 地基光学观测生成 | UTC MJD | 逐点读轨道，逐历元准备 EOP；JPL 太阳采用 TDB |
| 地月单程测距 | 输入真实 TDB 年月日时分秒；输出接收钟年月日时分秒 | 光行时迭代及地球 BCRS 位置平移 |
| 瞬转平 | UTC MJD | 岁差章动用 TT，日月解析平位置与格林尼治角取 UT1=UTC |
| 测站位置转换 | UTC MJD | 涉及 GCRS 时按同历元 EOP 计算地球定向；其余转换保留历元 |
| 两项直接状态转换 | 入口约定为 UTC MJD | 保留历元数值，逐行转换状态 |
| `Get_JPLEPH` / `Get_Libration` | TDB MJD | JPL 二进制星历插值 |

普通配置中的时间、文件中的 MJD 都保持“日”，不因状态单位码为 1 而归一化。传播时长和积分步长在配置中用 s，进入任务结构时再除以 `UNIT_TIME_S`。

### 状态属于哪一个参考系

地球轨道、测站惯性位置、初轨与 BODN 使用 GCRS。观测生成的目标星历沿用地心 J2000 惯性轴作为 GCRS 近似。ITRS 是地固直角坐标；测站表应随资料说明所用 ITRF 实现。地月入口的 E-BCRS 是把 BCRS 原点平移至地心而保留轴向，时间为 TDB。

Cartesian/Keplerian/NSG2 是状态的表示形式。转换这些表示时保留参考系，GCRS 与地月 E-BCRS 之间的转换需单独处理。GCRS↔ITRS 的矩阵接口用于位置，速度还需考虑矩阵随时间的变化。

测站 `Geodetic` 类采用项目常数中的参考椭球：WGS84 赤道半径 `6378137 m`，扁率 `1/298.257223563`。原始观测类型 8 标记为 WGS84 地固 XYZ，它的读写保持资料原值；当前 BODN 接收类型 1、6、7。

### 归一化的具体含义

常数定义见 [utils_constants.hpp](../include/utils_constants.hpp)，计算用途和来源见[物理常数与模型参数](11_constants.md)。下列名称位于 `constants::Scale` 中：

```text
长度单位：UNIT_LENGTH_M = 6378136.3 m
GM 单位：UNIT_GM_M3_S2 = 3.986004415e14 m^3/s^2
时间单位：UNIT_TIME_S = sqrt(UNIT_LENGTH_M^3 / UNIT_GM_M3_S2)
质量单位：UNIT_MASS_KG = UNIT_GM_M3_S2 / G
位置_norm = 位置_m / UNIT_LENGTH_M
速度_norm = 速度_m_s * UNIT_TIME_S / UNIT_LENGTH_M
时间差_norm = 时间差_s / UNIT_TIME_S
AMR_norm = AMR_m2_kg * UNIT_MASS_KG / UNIT_LENGTH_M^2
```

| 状态单位码 | Cartesian | Keplerian | NSG2 |
| --- | --- | --- | --- |
| 1 | 位置/UNIT_LENGTH_M，速度/(UNIT_LENGTH_M/UNIT_TIME_S) | a/UNIT_LENGTH_M，角度 rad | a/UNIT_LENGTH_M，lambda 为 rad，其余无量纲 |
| 2 | km、km/s | a 为 km，角度 deg | a 为 km，lambda 为 deg，其余无量纲 |
| 3 | m、m/s | a 为 m，角度 deg | a 为 m，lambda 为 deg，其余无量纲 |

这张表用于 GL/GLML 的 `state_units`、`output_state_units` 以及三套预报的输出。其他入口有固定文件单位：两项直接转换、观测生成真值、BODN 和瞬转平分别按下文读取。

## 统一轨道状态转换表

`aster_convert_states` 每行保存 UTC MJD 和[下表](#轨道状态的分量顺序)所列的六个状态分量。`--units deg-m`（默认）或 `--units deg-km` 同时决定输入输出的长度、速度单位；角度均为度。Cartesian 的位置、速度分别为 m、m/s 或 km、km/s；Keplerian 和 NSG2 的 a 为 m 或 km。NSG2 的 lambda 为度，xi、eta、h、k 无量纲。

## 位置坐标转换表

`aster_convert_coordinates` 的输入、输出都为四列。大地坐标排列为 `MJD_UTC longitude_deg latitude_deg height_m`，ITRS、GCRS 排列为 `MJD_UTC X_m Y_m Z_m`。上述列名采用默认的 `--units deg-m`。可选 `--units deg-km` 时，经纬度仍为度，输入输出的高度及 XYZ 均改为千米，输出列名相应使用 `_km`。东经为正，纬度是 WGS84 椭球法线定义的大地纬度，高度是椭球高。反算经度为 −180° 至 180°，精确自转轴上取 0。

历元列始终按 UTC MJD 原值保存；同一行的 EOP 由现有读取器插值取得。文件采用空白分隔，允许空行、整行或行末 `#` 注释。输出带有说明坐标、单位和转换模型的 `#` 表头，可直接用 `numpy.loadtxt` 读取。

此文件的经纬高顺序是经度、纬度、高度；`frame::Geodetic` 构造函数及 `ToArray` 的顺序为高度、经度、纬度，批量主程序在入口和出口完成排列转换。

## 轨道状态的分量顺序

| 表示 | 六分量 | 含义 |
| --- | --- | --- |
| Cartesian | `x y z vx vy vz` | 同历元的位置和速度 |
| Keplerian | `a e i OMEGA omega M` | 半长轴、偏心率、倾角、升交点赤经、近地点幅角、平近点角 |
| NSG2 | `a xi eta h k lambda` | 第二类无奇点根数，定义如下 |

```text
xi = e*cos(OMEGA + omega)
eta = e*sin(OMEGA + omega)
h = sin(i/2)*cos(OMEGA)
k = sin(i/2)*sin(OMEGA)
lambda = M + omega + OMEGA
```

NSG2 的 xi、eta、h、k 都是无量纲数，只有 lambda 是角度。M 为平近点角，f 为真近点角，E 为偏近点角；输入表末项填写 M。

## 测站表

每个非空行五列，可用空行和 `#` 注释：

```text
站号 坐标类型 分量1 分量2 分量3
```

| 坐标类型 | 后三列 |
| --- | --- |
| 1 | 大地纬度 phi（deg）、经度 lambda（deg）、高度（m） |
| 2 | ITRS X、Y、Z（m） |

站号是字符串，观测文件须按原样匹配，例如 `1` 与 `0001` 是不同站号。站号须唯一，列数与坐标有限性在读入时检查。

C++ 的 `frame::Geodetic(height, lambda, phi, false)` 参数顺序是“高度、经度、纬度”，与表格顺序不同。读入原始记录后按字段名赋值，示例见[接口参考](07_interfaces.md#frames)。

## OBS 原始观测

每行前四列固定为：

```text
MJD_UTC type object_id observer_id
```

后续列按 type 决定。编号保留为字符串，数值使用空格分隔，读取时跳过空行和 `#` 注释，检查列数及有限值。

| type | 后续列，按文件中的顺序 | 单位 |
| --- | --- | --- |
| 1 地基赤经赤纬 | `RA Dec` | deg、deg |
| 2 地基方位高度 | `Az El` | deg、deg |
| 3 距离＋方位高度 | `range Az El` | m、deg、deg |
| 4 距离 | `range` | m |
| 5 距离率 | `range_rate` | m/s |
| 6 天基赤经赤纬 | `RA Dec observer_X observer_Y observer_Z` | deg、deg，平台 GCRS 位置 m |
| 7 目标惯性位置 | `X Y Z` | 目标自身的 GCRS 位置 m |
| 8 目标地固位置 | `X Y Z` | 目标自身的 WGS84 地固位置 m |

原始资料类支持八种类型，实际计算入口按任务选用：GL/GLML 使用类型 1，BODN 使用 1、6、7，当前观测生成器输出类型 1。含地基观测时，站号须能在测站表中找到。测角记录的赤纬严格位于 −90° 与 90° 之间。

GL/GLML 使用同一目标、至少三条观测，历元非递减；BODN 同样处理一个目标，允许混合支持的类型和多站，历元非递减。同历元多站观测可以相邻保存。

## 星历、转换输入和真值

| 入口或文件 | 一行的列序 | 单位及时间 |
| --- | --- | --- |
| 根数转 Cartesian 的输入 | `MJD_UTC a e i OMEGA omega M` | km、deg |
| 根数转 Cartesian 的输出 | `MJD_UTC x y z vx vy vz` | m、m/s |
| Cartesian 转根数的输入 | `MJD_UTC x y z vx vy vz` | 命令行选 `m` 或 `km`，速度相应为 m/s 或 km/s |
| Cartesian 转根数的输出 | `MJD_UTC a e i OMEGA omega M` | km、deg |
| 瞬转平输入/输出 | `MJD_UTC a e i OMEGA omega M` | 密切/均根数，km、deg |
| 地球预报的三种星历 | `MJD_UTC` 加对应六分量 | 按 `output_state_units` |
| 地月预报星历 | `MJD_TDB x y z vx vy vz` | E-BCRS，按 `output_state_units` |
| 观测生成器输入 | `ID MJD_UTC` 加 Cartesian 或 Keplerian 六分量 | 8 列；0 模式为 km、km/s、deg，1 模式为归一化/rad |
| 观测生成器 `kep_*.txt` | `MJD_UTC a e i OMEGA omega M` | km、deg |
| 观测生成器 `cart_*.txt` | `MJD_UTC x y z vx vy vz` | km、km/s |

两项直接转换的每个非空行恰好七个有限数值，文件中不加标题行。观测生成器的 ID 为整数，历元严格递增，逐采样生成观测；可见点的 ID 连续性用于分弧，弧段中的观测与真值逐行对应。

把地球预报星历交给生成器时，先选择单位码 2，并在七列星历前添加从 1 递增的整数 ID；`file_eph_Cartesian` 或 `file_eph_Keplerian` 与状态表示一致，`is_input_normalized=0`。转换程序输出的 Cartesian 是米，需先转成 km、km/s 再作为同模式生成器输入。

观测生成的 `all_*.txt` 含所有采样：`MJD_UTC type object observer RA_deg Dec_deg visible input_id`。它用于检查可见性；定轨输入使用六列 `obs_*.txt` 或加噪后的六列观测。

BODN 后验星历的列数与普通预报不同：Cartesian 是 `index MJD_UTC x_km y_km z_km vx_km_s vy_km_s vz_km_s`；Keplerian 是 `object MJD_UTC a_km e i_deg OMEGA_deg omega_deg M_deg`；NSG2 是 `object MJD_UTC a_km xi eta h k lambda_deg`。三者使用同一组采样历元。

## 初轨与后验轨道解

GL、GLML、BODN 以 SimpleInput 文本写出轨道。通用字段如下：

```text
object_id = 目标号
state_units = 2
reference_epoch = UTC_MJD
reference_Keplerian = a e i OMEGA omega M
reference_Cartesian = x y z vx vy vz
```

这里是字段示意，数值和状态须来自实际轨道。`state_units` 使用前文的 1/2/3 约定，BODN 固定输出 2。GL 的 `file_initial_orbit` 读取同参考历元初值；BODN 的 `file_apriori_orbit` 可读不同历元先验并传播。文件同时有两种状态时，现有读取器优先使用 Cartesian。

其他结果字段按入口解释：

| 字段 | 含义 |
| --- | --- |
| `rms_angle_rad` | 两个测角方向全部分量的合成 RMS，rad |
| `solution_flag` | GL：正数为收敛轮次；−1 为达到迭代上限；−2 为数值失败 |
| `negative_log_likelihood` | GLML 最终选轨的 J，越小越优 |
| `selected_candidate` | GLML 候选编号，与 CSV 的 id 一致 |
| `random_seed` | GLML 回退使用的基础种子 |
| `source_solution_flag` | GLML 边界处理前的来源求解状态；编码同 GL flag |
| `rms_position_m` | BODN 有位置观测时写出的物理位置分量 RMS，m |
| `drag_amr_m2_kg` / `srp_amr_m2_kg` | BODN 实际启用力项的后验 AMR，m²/kg；也可作为下一次定轨的先验 |

## 残差、似然和求解状态

测角使用观测减计算（O−C），东西方向取赤经差主值并乘观测赤纬的余弦：

```text
r_E = wrap(RA_obs − RA_calc)*cos(Dec_obs)
r_N = Dec_obs − Dec_calc
RMS_angle = sqrt((sum(r_E^2) + sum(r_N^2)) / (2*N))
J = N*ln(2*pi*sigma^2) + (sum(r_E^2) + sum(r_N^2))/(2*sigma^2)
```

角度和 sigma 在这些式子中均为 rad，N 是观测次数，每次有两个分量。J 的常数项可为负，负的 J 是正常值；候选排序在同一份观测、同一 sigma 下进行。当前东西/南北等标准差的加噪方式见[观测生成](06_observations.md)。GL 结束后的遗留 Ljung–Box 检验已停用，Q/p 字段为 −1，GLML 按上述 J 选优。

GL/GLML 残差文件为六列：

```text
MJD_UTC type object observer r_E_rad r_N_rad
```

BODN 残差文件为九列：

```text
index MJD_UTC type accepted r1 r2 r3 object observer
```

`accepted=1` 表示采用，0 表示剔除，文件保留原始记录顺序。类型 1/6 的三个残差分量为 `r_E(rad) r_N(rad) NaN`；类型 7 为 `dX/UNIT_LENGTH_M dY/UNIT_LENGTH_M dZ/UNIT_LENGTH_M`，要恢复米制残差需乘 UNIT_LENGTH_M。库内 `BatchResult::residuals` 的类型 7 分量直接为 m。

GLML 候选 CSV 的表头为：

```text
id,a(km),e,i(deg),RAAN(deg),omega(deg),M(deg),RMS,J,x(km),y(km),z(km),vx(km/s),vy(km/s),vz(km/s),source_flag
```

RMS 为 rad，RAAN 即 OMEGA；状态使用边界处理后的轨道，source_flag 属于处理前的来源。无法完成边界处理或残差评价的候选以 NaN 标记，J 仅在有限时参与选优。

## 科学数据文件

EOP、JPL、月球重力场和大气密度表保存在 `data/`。[数据说明](../data/README.md)记录覆盖范围，[manifest.json](../data/manifest.json)记录来源、大小和 SHA-256。JPL 文件是当前读取器支持的二进制格式；EOP 插值需要连续日记录和必要的下一日数据，超出覆盖时通过 tracer 报错。

<a id="cislunar"></a>
## 地月单程测距状态、事件与观测

状态文件每行十三列：

```text
year month day hour minute second x y z vx vy vz clock_offset_s
```

年月日时分由整数保存，秒保留小数；历元为真实 TDB，位置和速度采用 E-BCRS，单位为 m 和 m/s。`clock_offset_s` 为该事件的卫星钟读数减 TDB，单位 s。位置、速度和钟差对应同一事件，各行按历元严格递增。

接收事件文件每行八列：

```text
year month day hour minute second transmitter_id receiver_id
```

这里的时刻是真实 TDB 接收时刻。发射端状态由光行时迭代按需查询。状态和事件文件均可使用空行及 `#` 注释。

输出 `observations.txt` 每行九列：

```text
year month day hour minute second transmitter_id receiver_id pseudorange_m
```

输出年月日时分秒是接收钟读数，即真实接收 TDB 加接收钟差；末列为加入系统差和白噪声后的伪距，单位 m。钟差跨越午夜时，输出日历会相应进位或借位。

`components.csv` 用列名保存真实发收历元、光行时、两个发收事件的 BCRS 几何距离、引力时延、发射钟差、接收钟差、钟差距离贡献、理论伪距、系统差、噪声、最终伪距和迭代量。其 `rx_*`、`tx_*` 历元字段均为真实 TDB。文件开头的注释记录所用历表版本、GM 来源和引力天体。

理论伪距等于几何距离，加上光速乘引力时延，再加上光速乘“接收钟差减发射钟差”。最小算例和状态来源说明见 [generate_cislunar_observations](../tests/generate_cislunar_observations/README.md)。
