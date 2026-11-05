# 地球固体潮最小算例

本算例在原始 GGM03C 静态场上加入日月二阶固体潮修正，预报同一颗低轨卫星一天。初态采用 GCRS 位置、速度，历元采用 UTC；程序内部使用 GCRS/TT。配置中的位置单位为 km，速度单位为 km/s。

模型取 `k2=0.3`、零滞后角。新增修正扣除了 GGM03C 已包含的永久潮贡献，原始球谐系数表保持不变。把 `use_solidtide_perturb` 改为 `0`，或者删去这一项，即只保留静态场中已有的永久贡献。

## 运行

在仓库根目录执行：

```sh
cmake --preset release
cmake --build --preset release
mkdir -p build/release/solid-tide-example
build/release/bin/aster_opn_cartesian -I tests/earth_solid_tide/setting.txt
```

`setting.txt` 引用仓库中的 `data/eop/EOP-All_250516.txt` 与 `data/ephemeris/JPLEPH_C.430`。GGM03C 系数由库提供。算例开启地球 80×80 非球形引力、JPL 日月质点引力与固体潮修正，步长为 10 s，使用定步长 RKF78。

输出集中在 `build/release/solid-tide-example/`：`cart.txt` 保存 8,641 个 UTC 历元的 GCRS 状态，`report.txt` 保存模型和积分设置，`program.log` 保存原始配置和运行事件。重复运行时，将配置中的三个输出路径改到新的目录。使用 NSG2 入口时也采用这份配置，并为它设置单独的输出路径。

## 科学检查

```sh
ctest --test-dir build/release -R '^core_earth_solid_tide$' --output-on-failure
```

`nao_acceleration.tsv` 保存固定卫星、日月位置后，由 NAO 原始 `AccelSolidTideGeneral` 计算的完整二阶潮汐加速度。文件采用 SI 单位，列序及生成参数见 `source.json`。生成时将 NAO 的长度、GM 单位和日月 GM 与本轮 ASTER 固定点检查统一，Love 数取 0.3，滞后角取零。

检查程序将 ASTER 的附加修正加回永久部分，再与这份完整潮汐参考比较；另用独立的潮汐势函数数值求导，核对加速度、两极位置和两套动力学方程的接入。一天预报的完整比较见[验证记录](../../doc/validation/12_earth-solid-tide.md)。
