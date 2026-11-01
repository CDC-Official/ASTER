# ASTER 1.0.0

ASTER 是用于轨道研究的个人 C++ 程序库。本版包含 GL、GLML 初轨确定，三套轨道预报、批处理定轨，以及三项轨道状态转换程序。科学计算采用现有工程中已整理的实现，最小算例集中在 `tests/`，配套科学数据保存在 `data/`。

## 构建和检查

需要 CMake 3.26 或更新版本、支持 C++17 的 C/C++ 编译器、Boost 头文件和 OpenMP。Eigen 3.4.0 与 SOFA 2018-01-30 源码已随库保存。

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

若 Boost 不在默认搜索路径中，配置时添加 `-DASTER_BOOST_ROOT=/path/to/boost`。Mac 上的 OpenMP 由现有 libomp 提供；具体工具链要求见 [构建与迁移](doc/build.md)。

先运行一个 GL 算例：

```sh
python3 tests/run_case.py --case iod_gl --binary build/release/bin/aster_iod_gl --output build/manual/iod_gl
```

脚本会创建本次运行目录，准备配置并比较科学结果。输出中的目录路径可直接用于查看配置、控制台和轨道解。

## 程序入口

| 可执行程序 | 主程序 |
| --- | --- |
| `aster_iod_gl` | `apps/iod/main_iod_gl.cpp` |
| `aster_iod_glml` | `apps/iod/main_iod_glml.cpp` |
| `aster_opn_cartesian` | `apps/opn/earth/main_opn_cartesian.cpp` |
| `aster_opn_nonsingular2` | `apps/opn/earth/main_opn_nonsingular2.cpp` |
| `aster_opn_cislunar` | `apps/opn/cislunar/main_opn_cislunar.cpp` |
| `aster_odn_batch` | `apps/odn/main_odn_batch.cpp` |
| `aster_convert_keplerian_to_cartesian` | `apps/conversion/main_convert_keplerian_to_cartesian.cpp` |
| `aster_convert_cartesian_to_keplerian` | `apps/conversion/main_convert_cartesian_to_keplerian.cpp` |
| `aster_convert_osculating_to_mean` | `apps/conversion/main_convert_osculating_to_mean.cpp` |

每个程序支持 `--version`。配置型程序使用 `-I setting.txt`；两项 Keplerian/Cartesian 转换直接接收输入和输出文件路径，详见 [程序使用](doc/programs.md)。

## 目录

`apps/` 保存九个主程序；`include/`、`src/` 保存公共接口与实现；`cmake/` 接入依赖；`third_party/` 保存固定版本源码及许可。`tests/` 包含小算例、科学结果参考和检查程序；`data/` 保存 EOP、JPL、重力场及大气数据；`doc/` 记录使用方式和版本来源。

- [科学数据及校验](data/README.md)
- [测试和验证范围](tests/README.md)
- [迁移来源与验收](doc/migration.md)
- [版本记录](CHANGELOG.md)
- [后续维护规则](AGENTS.md)
