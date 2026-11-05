# 测站位置转换与初轨输入一致性

六个方向的测站位置转换已具备内存接口和批量入口。以现有 GL 算例的站址及观测时刻生成 GCRS 位置时，一步接口与初轨主程序使用的分步调用逐分量相同；批量输出与内存结果也保持在纳米量级的一致性。以同一份 EOP 和测角观测比较 ASTER 与 Gooding 时，可以直接共用这组测站位置。

## 本轮实现及输入约定

本轮在 `feature/earth-solid-tide` 工作区完成，基线提交为 `3e88b79cca030c65405bbd495234483a480ed832`。原有 `frame::Geodetic`、`frame::ITRS`、`frame::GCRS` 已包含六向方法，本轮补齐了组合调用的错误传播和完整接口说明，并使用库内 SOFA `iauGc2gde` 完成地固位置到大地经纬高的反算。独立主程序为 `aster_convert_coordinates`，文件读取、批量循环和输出直接在主程序中组织。

大地坐标采用 WGS84 的 `a=6378137 m`、`f=1/298.257223563`。文件按 `UTC MJD、经度(deg)、大地纬度(deg)、椭球高(m)` 或 `UTC MJD、X(m)、Y(m)、Z(m)` 排列。四个涉及 GCRS 的方向使用 `data/eop/EOP-All_250516.txt`，由原有读取器插值出同历元 EOP，保持原有 IAU 2000B 地球定向矩阵。内存接口的椭球计算仍使用 `UNIT_LENGTH_M` 归一化高度和位置、rad 制角度。

GL、GLML 的求解流程与正向 WGS84、ITRS→GCRS 计算公式未改动。位置类、EOP 插值和转换矩阵均与批量入口共用；本轮另行补齐了错误输入的 Tracer。

## 参考与检查结果

九条小输入包括 GL 算例的前三次观测、一般站址、赤道、负高程、两极和近极点。SOFA `iauGd2gc(WGS84,...)` 以米制直接计算地固参考；地球定向参考由 SOFA `iauXys00b`、`iauC2ixys`、`iauPom00`、`iauC2tcio` 组合，再取逆旋转得到 GCRS。它们与归一化类接口分别计算。

| 比较双方或检查项目 | 最大差异 |
| --- | --- |
| ASTER 大地→ITRS 与 SOFA WGS84，位置分量 | 9.3132×10⁻¹⁰ m |
| ASTER 大地→GCRS 与 SOFA 合成矩阵参考，位置分量 | 9.3132×10⁻¹⁰ m |
| 内存六向接口的地固位置往返，位置分量 | 2.1243×10⁻⁹ m |
| 内存反算椭球高与输入高度 | 1.7395×10⁻⁹ m |
| 大地→GCRS 一步接口与初轨分步接口 | 全部分量相同 |
| 文件入口六向转换与内存参考，直角坐标分量 | 1.8626×10⁻⁹ m |
| 文件反算椭球高与输入高度 | 2.7137×10⁻⁹ m |

精确北极、南极另外用 `X=Y=0` 检查，经度取 0、纬度为 ±90°，恢复的椭球高通过 1 μm 门限。文件检查同时覆盖十八种错误输入与路径情况，包括非法纬度、缺列、非有限值、EOP 缺失或越界；失败时保留已有输出，并显示调用链。

GL 算例首条观测为 `58406.8859375000 UTC MJD`，测站经度为 `86.174633°`、大地纬度为 `41.725892°`、椭球高为 `943 m`。转换所得 GCRS 位置为：

```text
X = 1698803.440193041 m
Y = 4457900.520865699 m
Z = 4220621.665928426 m
```

两份参考输出与输入一起保存在[最小算例目录](../../tests/convert_coordinates/README.md)：[ITRS](../../tests/convert_coordinates/reference_itrs.txt)、[GCRS](../../tests/convert_coordinates/reference_gcrs.txt)。该目录同时提供六向命令与 Python 调用示例。

## 构建及现有程序对照

Mac arm64 / AppleClang 21、默认 OpenLibm Release 全目标构建通过。直接相关的六项检查通过：`core_frame_coordinates`、`core_station_coordinates`、`station_coordinates_runtime`、`observation_geometry`、`gravity_conventions` 和 `iod_batch_handoff`。

更新后的发布代表算例检查通过十二程序，共保存四十三份科学输出。将其中原有三十七份与 **v1.3.0 发布记录中的 Mac 科学输出 SHA-256** 比较，全部相同，覆盖 GL、GLML、批处理、三套预报及原有转换与观测入口。新增六份输出对应测站转换的六个方向。

本轮的源码快照及实际运行保存在 `build/release/station-coordinate-smoke/`，六向输出和数值检查保存在 `build/release/test-output/station_coordinates/`。输入校验、各分项数值、源码哈希和比较双方见[原始证据](evidence/19_station-coordinates.json)。

## 两套单位与状态转换入口

位置坐标入口统一命名为 `aster_convert_coordinates`，默认采用 `deg-m`，可选 `deg-km`；输入输出共用所选单位，经纬度始终为度。米制与千米制的输入分别位于 `tests/convert_coordinates/geodetic.txt`、`geodetic_km.txt`。

新增 `aster_convert_states` 复用三种现有轨道类，在 Cartesian、Keplerian、NSG2 之间转换，同样只提供度米、度千米两套单位。NSG2 的 a 为长度、lambda 为角度，其余四项无量纲。最小验证仅用当前两套入口做三对表示的往返，两个单位模式均通过；输入、命令和各分量见[状态算例](../../tests/convert_states/README.md)，结果保存在 `build/release/test-output/conversion_roundtrips/summary.json`。
