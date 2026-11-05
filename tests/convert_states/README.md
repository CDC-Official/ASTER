# 轨道状态转换最小算例

`aster_convert_states` 用两个参数选择输入、输出表示：`cart`、`kep`、`nsg2`。三种表示可以两两转换，每行七列，第一列为 UTC MJD，后六列按下表排列。转换保留历元与参考系。

| 表示 | 六分量顺序 | 量纲 |
| --- | --- | --- |
| `cart` | `x y z vx vy vz` | 位置为长度，速度为长度/秒 |
| `kep` | `a e i OMEGA omega M` | a 为长度，e 无量纲，其余为角度；M 为平近点角 |
| `nsg2` | `a xi eta h k lambda` | a 为长度，lambda 为角度，xi、eta、h、k 无量纲 |

NSG2 采用现有定义：`xi=e*cos(omega+OMEGA)`、`eta=e*sin(omega+OMEGA)`、`h=sin(i/2)*cos(OMEGA)`、`k=sin(i/2)*sin(OMEGA)`、`lambda=OMEGA+omega+M`。其中 OMEGA 为升交点赤经，omega 为近地点幅角。

单位只选两种，输入输出共用一套：

- `--units deg-m`：度、米，默认；Cartesian 速度为 m/s。
- `--units deg-km`：度、千米；Cartesian 速度为 km/s。

程序使用库内的地球 GM，即 `constants::Scale::UNIT_GM_M3_S2`，配套现有归一化尺度，处理地球椭圆轨道。`keplerian_km.txt` 保存两条普通椭圆轨道，a 为千米。

从 ASTER 根目录运行：

```sh
build/release/bin/aster_convert_states kep nsg2 \
    tests/convert_states/keplerian_km.txt build/states_nsg2.txt --units deg-km
build/release/bin/aster_convert_states nsg2 cart \
    build/states_nsg2.txt build/states_cartesian.txt --units deg-km
build/release/bin/aster_convert_states cart kep \
    build/states_cartesian.txt build/states_keplerian.txt --units deg-km
```

文件允许空行和 `#` 注释，输出附有列名、单位和 GM。同名输出在全部输入转换成功后覆盖。

最小检查只运行本程序的三对表示往返，并在两套单位下重复；位置坐标入口也按相同方式检查：

```sh
ctest --test-dir build/release -R '^conversion_roundtrips$' --output-on-failure
```

输出及摘要保存在 `build/release/test-output/conversion_roundtrips/`。
