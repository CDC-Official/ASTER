# 位置坐标转换最小算例

从 ASTER 根目录运行。`geodetic.txt` 保存九条 UTC 时刻与 WGS84 经纬高，前三条使用 GL 初轨算例的测站和前三次观测时刻，适合准备 Gooding 所需的测站位置。其余记录覆盖南半球、赤道、负高程、两极和近极点。

## 单位选项

命令末尾可加 `--units deg-m`（度、米，默认）或 `--units deg-km`（度、千米）。所选单位同时用于输入和输出；经纬度始终为度，椭球高及 XYZ 采用所选长度单位。下表和六向示例采用默认的度、米。

千米制输入见 [geodetic_km.txt](geodetic_km.txt)，与米制输入逐行对应。使用方式为：

```sh
mkdir -p build/station-demo
build/release/bin/aster_convert_coordinates geodetic gcrs \
    tests/convert_coordinates/geodetic_km.txt build/station-demo/gcrs_km.txt \
    data/eop/EOP-All_250516.txt --units deg-km
```

## 输入与输出

| 坐标名称 | 每行四列 | 单位与约定 |
| --- | --- | --- |
| `geodetic` | `MJD_UTC longitude latitude height` | 日、度、度、米；东经为正，纬度为椭球大地纬度，高度为 WGS84 椭球高 |
| `itrs` | `MJD_UTC X Y Z` | 日、米、米、米；地心地固直角坐标 |
| `gcrs` | `MJD_UTC X Y Z` | 日、米、米、米；地心天球直角坐标 |

空行及 `#` 注释均可使用，输出保持输入的历元与行序。大地坐标结果的经度位于 −180° 至 180°，精确自转轴上经度取 0。不同 ITRF 实现的站址应在实验记录中注明实际来源。

椭球采用 `constants::Wgs84` 的 `a=6378137 m`、`f=1/298.257223563`。涉及 GCRS 的四个方向还需提供覆盖各时刻的 EOP，读取和插值方式与 GL、GLML 相同；地球定向使用现有 IAU 2000B 矩阵，包括极移与 CIP 改正。EOP 插值通常需要下一日记录。

[ITRS 参考输出](reference_itrs.txt)和[GCRS 参考输出](reference_gcrs.txt)与输入逐行对应，均保留四列表头及单位。

## 运行六个方向

```sh
mkdir -p build/station-demo
build/release/bin/aster_convert_coordinates geodetic itrs \
    tests/convert_coordinates/geodetic.txt build/station-demo/itrs.txt
build/release/bin/aster_convert_coordinates geodetic gcrs \
    tests/convert_coordinates/geodetic.txt build/station-demo/gcrs.txt \
    data/eop/EOP-All_250516.txt
build/release/bin/aster_convert_coordinates itrs geodetic \
    build/station-demo/itrs.txt build/station-demo/geodetic_from_itrs.txt
build/release/bin/aster_convert_coordinates itrs gcrs \
    build/station-demo/itrs.txt build/station-demo/gcrs_from_itrs.txt \
    data/eop/EOP-All_250516.txt
build/release/bin/aster_convert_coordinates gcrs itrs \
    build/station-demo/gcrs.txt build/station-demo/itrs_from_gcrs.txt \
    data/eop/EOP-All_250516.txt
build/release/bin/aster_convert_coordinates gcrs geodetic \
    build/station-demo/gcrs.txt build/station-demo/geodetic_from_gcrs.txt \
    data/eop/EOP-All_250516.txt
```

## 从 Python 准备 Gooding 输入

`epochs_utc` 填入本次测角观测的 UTC MJD，每个时刻对应一行站址。下面使用 GL 算例的三次观测，程序执行完后 `station_gcrs_m` 就是对应的三行 GCRS 位置，单位为米。

```python
from pathlib import Path
import subprocess
import numpy as np

run = Path('build/station-demo')
run.mkdir(parents=True, exist_ok=True)
epochs_utc = [58406.8859375000, 58406.8859722222, 58406.8860069444]
rows = [[utc, 86.174633, 41.725892, 943.0] for utc in epochs_utc]
np.savetxt(run / 'station_geodetic.txt', rows, fmt='%.17g')
subprocess.run([
    'build/release/bin/aster_convert_coordinates', 'geodetic', 'gcrs',
    str(run / 'station_geodetic.txt'), str(run / 'station_gcrs.txt'),
    'data/eop/EOP-All_250516.txt'
], check=True)
station_gcrs_m = np.loadtxt(run / 'station_gcrs.txt')[:, 1:4]
```

这些站址与原有测角观测一起交给采用相同 GCRS 轴向的 Gooding 调用。内存接口和归一化规则见[公共接口](../../doc/07_interfaces.md#frames)，自动检查与数值结果见[测站转换验证](../../doc/validation/13_station-coordinates.md)。
