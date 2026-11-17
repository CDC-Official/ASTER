# GCRS/TT 核对分支验证

## 源码与计算条件

测试入口实现提交为 **598e1813837c3f037d4a9b7b3dfb1d8c9ae9119b**，基于主分支 `c87a531ecdef2bba14e1cf8980dca08c8d024dd6`，分支为 `codex/gcrs-tt-orbit-check`。完整参考轨迹在实现提交的干净工作区中重新构建、运行生成；随后提交的参考资料与本文记录保持同一套计算源码。

运行平台为 macOS 26.6.2 / arm64、Apple clang 21.0.0、Release、**SYSTEM 数学后端**。Windows 的 Visual Studio 2022 preset、UTF-8 源码选项和构建说明已提供；本轮实际编译与运行在 macOS 上完成。

初值为 UTC MJD 61892.0，GCRS 位置 `[-172078.764589166,316878.877578190,136905.589039608] km`，速度 `[-0.343284732289,-0.175709556631,-0.114628368645] km/s`。模型为地球质点、GGM03C 12×12 阶、DE440 日月质点；月球非球形、其他行星、光压、阻力、潮汐、1PN 均关闭。GM、历表查询和缺测 EOP 约定见[交付说明](../18_gcrs-tt-orbit-check.md)。

## 已完成检查

SYSTEM 全目标构建通过。以下 5 项 CTest 全部通过：

| 检查 | 主要覆盖内容 |
| --- | --- |
| `core_frame_coordinates` | 现有坐标转换回归 |
| `core_rkf78_safety` | RKF78 解析参考、拒绝与失败行为等现有检查 |
| `gravity_conventions` | 现有地球和月球引力尺度约定检查 |
| `core_gcrs_tt_check` | 默认严格 EOP 行为、零定向回退、DAT 保留、插值和末日边界、非法历元；八颗行星加速度与独立牛顿直接／间接项公式比较 |
| `gcrs_tt_check_runtime` | 125 s 时长的 60+60+5 s 末步、显式零 EOP、圆轨道解析解、月球文件头 GM 改变后的输出一致性、月球和行星开关、零时长与非法输入 |

月球文件头 GM 覆盖检查把文件头 GM 人为加倍，保留参考半径和所有 C/S，开启 12 阶月球球谐。两次轨迹逐值相同，说明球谐计算实际使用入口指定的 DE440 TT 月球 GM。该检查的修改文件只存在于构建目录；仓库的月球场数据保持原样。

独立圆轨道检查采用地球 GM=398600.4415 km³/s²、半径 7000 km、解析圆轨道速度，预报 600 s、步长 5 s，位置分量与解析解的差异均小于 `1e-7 km`，速度分量均小于 `1e-10 km/s`。

## 十天参考轨迹

RKF78 固定 60 s，预报 864000 TT s，共输出 **14401** 行；末点 UTC MJD **61902.0**，TT MJD **61902.000800740738**。全时段 EOP 均超出文件覆盖范围，187200 次右函数查询按约定使用零定向参数，DAT=37 s。

末状态为：

```text
r_GCRS (km)   = [-80006.714747643215, -345511.81276058516, -172415.17830684679]
v_GCRS (km/s) = [0.34752202371829721, -0.18836278607179291, -0.058359701674503849]
```

完整输出见 [reference_cartesian.csv](../../tests/gcrs_tt_orbit_check/reference_cartesian.csv)，运行参数与实际 GM 见 [reference_summary.txt](../../tests/gcrs_tt_orbit_check/reference_summary.txt)，输入和输出校验值见 [reference_evidence.json](../../tests/gcrs_tt_orbit_check/reference_evidence.json)。

## 60 s 与 30 s 步长复算

双方均为上述同一 ASTER 实现，保持初值、DE440、模型和 GM 相同；仅将 RKF78 固定步长从 60 s 改为 30 s。按每 60 s 的共同历元比较 GCRS 位置的欧氏距离：

| 统计量 | 位置差 |
| --- | ---: |
| 全时段最大值 | **0.001940942696 m** |
| 最大值所在 TT 累计秒 | **858300 s** |
| 十天末点 | **0.001895079656 m** |

这一复算包含步长变化、分级历元取整和浮点舍入的综合影响，用于检查本算例对步长的敏感性。甲方程序的轨迹尚未纳入本次分支验证。

## 复现命令

从仓库根目录配置 SYSTEM，并保留测试目标：

```sh
cmake -S . -B build/client -DASTER_MATH_BACKEND=SYSTEM -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/client --parallel 8
ctest --test-dir build/client -C Release -R '^(core_gcrs_tt_check|gcrs_tt_check_runtime|core_rkf78_safety|gravity_conventions|core_frame_coordinates)$' --output-on-failure
build/client/bin/aster_opn_gcrs_tt_check -I tests/gcrs_tt_orbit_check/setting.txt
```

30 s 复算使用同目录中复制的配置，将 `integrator_stepsize` 改为 30，并另设输出目录；取其每两行一个节点与 60 s 结果对齐。Windows 使用专用 preset，执行路径见交付说明。


## 恢复原 setting 关键词后的核对

测试入口已恢复原 Cartesian 配置中的 `reference_Cartesian`、`propagation_length`、`integrator_stepsize`、`file_JPLEPH` 等字段，地月带谐与田谐继续分别设置。配置旁补充单位、GM 和模型说明；主程序、运行检查和使用文档同步更新。

恢复后的十天计算与上述完整参考 CSV 逐字节一致，14401 行的各历元和六维状态全部相同。运行检查还覆盖原 `reference_epoch` 的 MJD、JD、年月日时分秒形式，以及分别设置地球带谐、田谐的行为。
