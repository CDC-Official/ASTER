# 物理常数与模型参数

ASTER 将全程序共用的常数直接放在 `constants` 中，来源写在定义旁；需要区分用途或配套使用的模型参数保留子命名空间。重力场的 GM、参考半径和系数共同决定引力势；参考椭球的半径和扁率共同决定大地坐标；程序归一化尺度只决定数值单位。调用处可在无歧义的源文件作用域内用 `using` 引入常量；不同来源的 GM 和参考半径保留明确的模型名称。

## 开源实现中采用的组织方式

- [Orekit 13.1 Constants](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/utils/Constants.java) 用 WGS84、GRS80、IERS 等来源区分常数；[OneAxisEllipsoid](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/bodies/OneAxisEllipsoid.java) 同时接收半径和扁率。
- [GMAT GravityField](https://github.com/nasa/GMAT/blob/main/src/base/forcemodel/GravityField.cpp) 从重力场模型取得配套 GM 和参考半径。
- [Tudat 球谐重力场](https://py.api.tudat.space/en/latest/dynamics/environment_setup/gravity_field.html) 将 GM、参考半径、系数和固连坐标系放在同一模型中；[形状模型](https://py.api.tudat.space/en/latest/dynamics/environment_setup/shape.html) 单独定义椭球。

ASTER 用命名空间管理固定参数，继续用现有 `LunarGravityField` 保存从文件读取的模型参数。来源确定的固定参数放在 `utils_constants.hpp`；大型系数表和特定近似模型的参数留在所属模块。

## 常数存放与调用

π、光速、万有引力常数、天文单位及通用时间常数直接通过 `constants::名称` 调用。下表其余分组是 `constants` 的子命名空间。

| 分组 | 参数及单位 | 来源与用途 |
| --- | --- | --- |
| `constants`（通用） | 真空光速 299792458 m/s | SI 定义值，用于传播时延、测距和后牛顿项 |
| `constants`（通用） | G=6.67430×10⁻¹¹ m³/(kg·s²) | [CODATA 2022](https://arxiv.org/abs/2409.03787)，用于质量尺度及其单位换算 |
| `constants`（通用） | AU=149597870700 m | [IAU 2012 B2](https://www.iau.org/static/resolutions/IAU2012_English.pdf)，用于解析日月模型中的长度换算 |
| `IERS2010` | TT 地球 GM=3.986004415×10¹⁴ m³/s²；名义地球角速度 7.292115×10⁻⁵ rad/s | [IERS 2010](https://iers-conventions.obspm.fr/content/tn36.pdf) 第 1 章表 1.1 及其后的 TT 换算说明、第 6.1 节；分别用于 GCRS/TT 动力学和共转大气速度近似 |
| `DE440` | 11 个 TDB 相容 GM，m³/s² | [Park 等，2021，表 2](https://doi.org/10.3847/1538-3881/abd414)；数组按 JPL 编号 1—11 排列，火星至冥王星采用行星系统 GM |
| `GGM03C` | GM、参考半径、完全归一化 C/S，以及由 C20 得到的 J2 | [CSR GGM03 发布资料](https://download.csr.utexas.edu/pub/grace/GGM03/GGM03_Archive.zip)；GM=3.986004415×10¹⁴ m³/s²，参考半径 6378136.3 m；系数位于 `data_gravityfield.hpp` |
| `Wgs84` | 赤道半径 6378137 m，扁率 1/298.257223563，第一偏心率平方 | WGS84 参考椭球；与固定版本 SOFA 的 `eform.c`、Orekit WGS84 参数一致，用于测站大地坐标、大气参考面与地影几何 |
| `SolarRadiation` | 参考距离 149597870000 m，参考压力 4.56×10⁻⁶ N/m² | [Orekit 13.1 SolarRadiationPressure](https://github.com/CS-SI/Orekit/blob/13.1/src/main/java/org/orekit/forces/radiation/SolarRadiationPressure.java) 默认参数；数值预报和解析光压共同使用 |
| `Scale` | 长度、GM、时间、质量尺度及归一化后的天体 GM | 长度取 GGM03C 参考半径，GM 取 IERS TT 地球 GM；时间为 `sqrt(UNIT_LENGTH_M³/UNIT_GM_M3_S2)`，质量为 `UNIT_GM_M3_S2/G` |
| `Conversion` | 速度、长度、时间、角度、面质比换算系数 | 从 `Scale` 和单位定义推导 |
| `constants`（通用） | 日长、世纪长度、日期原点、TT/TDB 比例及 π | 日期常数与 L_G、L_B 使用固定版本 SOFA 定义；π 为数学常数 |

`Scale::NORMALIZED_BODY_GM_TDB` 用 DE440 GM 除以程序 GM 尺度。`Scale::NORMALIZED_MOON_GM_TT` 和 `Scale::NORMALIZED_SUN_GM_TT` 则先做 TT 相容换算，再除以尺度。两者对应不同的动力学时间变量。

GGM03C 和 IERS 的地球 GM 数值相同，但分别表示重力场模型的配套参数和 TT 动力学约定。修改一项时，应根据模型关系检查另一项，而非仅因数值相等而合并来源。

## 各计算场景使用哪组参数

| 场景 | 实际调用 | 代码位置 |
| --- | --- | --- |
| GCRS/TT 地球质点引力、轨道根数转换 | IERS TT 地球 GM 所定义的归一化单位，地球质点引力参数为 1 | `accel_complete.cpp`、`orbit_basic.cpp` |
| E-BCRS/TDB 地球质点引力、行星摄动、后牛顿项 | `Scale::NORMALIZED_BODY_GM_TDB`；地月程序的根数输入还换算中心 GM 对速度的影响 | `accel_complete.cpp`、`main_opn_cislunar.cpp` |
| GCRS/TT 日月摄动 | `Scale::NORMALIZED_MOON_GM_TT`、`Scale::NORMALIZED_SUN_GM_TT` | `accel_thirdbody.cpp`、`accel_complete.cpp` |
| 地球球谐引力 | `GGM03C::GM`、`REFERENCE_RADIUS`、`COEFFICIENTS`；TDB 入口同时缩放 GM 与半径 | `accel_gravity.cpp`、`accel_complete.cpp` |
| 初轨 J2、解析 STM、瞬转平 | J2 从同一 C20 推导；GM 比与参考半径幂显式参与归一化 | `core_iod.cpp`、`orbit_stm.hpp`、`analyt_averaging.cpp` |
| 月球球谐引力 | `LunarGravityField` 中从文件读取的 GM、参考半径和 C/S；用于二阶及以上项 | `accel_gravity.hpp/.cpp`、`accel_complete.cpp` |
| 测站大地坐标与地固坐标 | `Wgs84` 半径、扁率及偏心率 | `frame_definition.cpp` |
| 天球与地固坐标旋转 | EOP 和 SOFA 的 IAU 2000/2006 实现 | `frame_earth.cpp`；椭球常数不参与旋转矩阵 |
| 大气阻力 | `Wgs84` 参考面、IERS 名义角速度、输入密度表、卫星 Cd 与面质比 | `model_atmosphere.cpp`、`accel_drag.cpp` |
| 数值和解析光压 | 同一组 `SolarRadiation` 参数及卫星 Cr、面质比；由物理量换算加速度 | `accel_solar_radiation.cpp`、`analyt_averaging.cpp` |
| 地影 | `Wgs84` 半径除以长度尺度 | `in_shadow`、`Check_OrbitEclipse`；数值预报采用地球圆柱地影，进入阴影后光压为零 |
| 单程测距及引力传播时延 | 真空光速、DE440 TDB GM | `model_cislunar_observation.cpp` |
| 时间转换 | SOFA 的定义常数和 HF2002 的模型系数 | `time_transform.cpp` |

## 专用模型参数与数据

`model_ephemeris_analytic.hpp` 的 `MeanOrbit` 保存解析模型共用的日月半长轴和 J2000 平黄赤交角；根数及其变化率见汤靖师《人造卫星精密定轨(II)》2025 讲义第 104 页（项目 `part2.pdf` 第 95 页）。同一组距离和黄赤交角也用于解析 STM 的日月平均近似。该近似把日月平均轨道面取为平黄道面。各根数系数留在解析历表实现中，并与其模型说明对应。

解析田谐项使用 GMST 多项式在 J2000 的斜率，和共转大气的名义角速度分开。多项式对应 IAU 1982 平恒星时表达，系数可对照固定版本 SOFA `gmst82.c`；这项解析近似不改变预报坐标旋转所用的 IAU 2000/2006。

HF2002 的频率、幅值、起始历元和线性修正是拟合模型的一部分，集中保存在 `time_transform.cpp`，对应 IERS 第 10 章的 HF2002 实现。它们与通用的日期原点、日长和 IAU 定义速率分别管理。

月球重力数据已经与 NASA PDS 的 GRGM660PRIM 逐行核对。随库文件保留 80 阶，GM、参考半径及标准差的表头使用 SI 单位，转换步骤和 SHA-256 见 [科学数据说明](../data/README.md)与 `data/manifest.json`。其发布坐标系为 DE421 月球主轴系，当前计算使用所加载历表的天平动。GM 和参考半径按文件发布值使用，TT/TDB 相容约定需结合该模型的动力学定义确认。

大气表头明确注明数据来自中科院国家天文台 H. W. Cheng、由南京大学 J. S. Tang 生成，条件为 Ap=10、F10.7=150。密度和标高属于该输入模型，来源、单位及校验值见数据目录。EOP 和 JPL 文件也由各自读取器处理；JPL 文件头的 AU、EMRAT 用于解析该历表，独立于动力学选择的 DE440 GM。

## 使用与维护

新增科学计算先确定时间尺度、参考系和模型，再显式引用对应分组。例如测站转换使用 `constants::Wgs84::ECCENTRICITY_SQUARED`；TDB 测距使用 `constants::DE440::GM_TDB[body - 1]`；地球球谐场同时传入 `constants::GGM03C::GM` 与 `constants::GGM03C::REFERENCE_RADIUS`。

采用其他重力场时，连同 GM、参考半径、归一化方式和潮汐约定一起核对。采用其他参考椭球时，同时替换半径和扁率。修改数值尺度时，从有量纲公式检查面质比、密度、光压和各阶重力项的单位换算。

此次修改的控制对照和数值影响见 [常数规范化验证](validation/06_constants.md)。

## 从旧常数名称迁移

直接调用库的代码按下表更新名称。配置文件继续填写物理单位，归一化换算集中使用当前单位尺度。

| 旧名称 | 当前调用 |
| --- | --- |
| `pi` | `constants::PI` |
| `ae` | `constants::Scale::UNIT_LENGTH_M`；若表示测站椭球半径，使用 `constants::Wgs84::EQUATORIAL_RADIUS` |
| `UNITTIME` | `constants::Scale::UNIT_TIME_S` |
| `GE`，作为积分 GM 单位 | `constants::Scale::UNIT_GM_M3_S2` |
| `MASS_OF_EARTH`，作为质量单位 | `constants::Scale::UNIT_MASS_KG` |
| `Units` 命名空间 | `constants::Conversion` |

速度、长度、时间和面质比的换算系数由单位尺度计算。需要模型物理 GM 时，按上文选用 `IERS2010`、`DE440` 或重力场自身的参数。
