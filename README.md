# ASTER 1.2.0（开发中）

ASTER 是用于轨道研究的个人 C++ 程序库。本版包含 GL、GLML 初轨确定，三套轨道预报、批处理定轨，三项轨道状态转换程序，以及地基光学观测生成和配套加噪工具。科学计算采用现有工程中已整理的实现，最小算例集中在 `tests/`，配套科学数据保存在 `data/`。

## 构建和检查

需要 CMake 3.26 或更新版本、支持 C++17 的 C/C++ 编译器、Boost 头文件和 OpenMP。Eigen 3.4.0、SOFA 2018-01-30 和官方 OpenLibm 0.8.8 源码及许可已随库保存。默认数学后端为 OPENLIBM，覆盖各 C++ 程序与公共/SOFA 路径；已验证平台为 Mac arm64 / AppleClang 21 与 Linux x86_64 / GCC 13.3。

```sh
cmake --preset release
cmake --build --preset release
ctest --test-dir build/release -R '^(core_|math_backend_|invalid_inputs|data_integrity)' --output-on-failure
python3 tests/check_math_backend_apps.py --bin-dir build/release/bin --output build/release/acceptance/apps
```

若 Boost 不在默认搜索路径中，配置时添加 `-DASTER_BOOST_ROOT=/path/to/boost`。Mac 上的 OpenMP 由现有 libomp 提供；具体工具链要求见 [构建与迁移](doc/build.md)。

先检查四个已有单 GL 样本：

```sh
python3 tests/math_backend_cases/run_checks.py --bin-dir build/release/bin --output build/release/acceptance/gl --only-gl
```

脚本使用仓库内的原样观测，创建独立运行目录并核对既有科学结果。完整旧 CTest 中的 1.0.0 字节参考继续保留；1.1.0 的后端验收入口和物理门限见 [测试说明](tests/README.md)。

使用系统数学时另建构建目录，配置 `-DASTER_MATH_BACKEND=SYSTEM`；默认值为 `OPENLIBM`。OPENLIBM 的程序与同目录 `libopenlibm.dylib` / `libopenlibm.so` 一起移动。构建与开关详见 [math_backend 说明与验收](doc/math-backend.md)。

## 程序入口

| 可执行程序 | 主程序 |
| --- | --- |
| `aster_generate_observations` | `apps/observation/main_generate_observations.cpp` |
| `aster_iod_gl` | `apps/iod/main_iod_gl.cpp` |
| `aster_iod_glml` | `apps/iod/main_iod_glml.cpp` |
| `aster_opn_cartesian` | `apps/opn/earth/main_opn_cartesian.cpp` |
| `aster_opn_nonsingular2` | `apps/opn/earth/main_opn_nonsingular2.cpp` |
| `aster_opn_cislunar` | `apps/opn/cislunar/main_opn_cislunar.cpp` |
| `aster_odn_batch` | `apps/odn/main_odn_batch.cpp` |
| `aster_convert_keplerian_to_cartesian` | `apps/conversion/main_convert_keplerian_to_cartesian.cpp` |
| `aster_convert_cartesian_to_keplerian` | `apps/conversion/main_convert_cartesian_to_keplerian.cpp` |
| `aster_convert_osculating_to_mean` | `apps/conversion/main_convert_osculating_to_mean.cpp` |

每个程序支持 `--version`。配置型程序使用 `-I setting.txt`；观测生成配置位于 `tests/generate_observations/setting.txt`，用法见 [观测生成](doc/observations.md)。两项 Keplerian/Cartesian 转换直接接收输入和输出文件路径，详见 [程序使用](doc/programs.md)。

## 目录

`apps/` 保存十个主程序；`include/`、`src/` 保存公共接口与实现；`cmake/` 接入依赖；`third_party/` 保存固定版本源码及许可。`tests/` 包含小算例、科学结果参考和检查程序；`data/` 保存 EOP、JPL、重力场及大气数据；`doc/` 记录使用方式和版本来源。

- [观测生成、加噪和最小验证](doc/observations.md)
- [科学数据及校验](data/README.md)
- [测试和验证范围](tests/README.md)
- [迁移来源与验收](doc/migration.md)
- [服务器使用与数值对照](doc/server.md)
- [1.1.0 更新公告](https://github.com/CDC-Official/ASTER/releases/tag/v1.1.0)
- [版本记录](CHANGELOG.md)
- [后续维护规则](AGENTS.md)
