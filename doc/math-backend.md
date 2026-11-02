# math_backend 开发说明

`dev/1.1.0` 从 `963e44c` 建立。开发目录为 Mac 的 `/Users/chendingcun/CLionProjects/ASTER-dev-1.1.0` 和服务器的 `/home/cdc/projects/ASTER-dev-1.1.0`。原 `main` 工作区及其未提交试验保留，现有 `v1.0.0` 标签不动。

## 接入与构建

`include/math_backend.hpp`、`src/math_backend.cpp` 和 `cmake/MathBackend.cmake` 集中管理公共入口、版本信息与浮点编译条件。九个正式程序共用 `aster_core` 和 `aster_sofa`，通过 `ASTER::math_backend` 传播对象和规则；不再构建 IOD 专属 core/SOFA。新增调用者链接 `ASTER::core` 或 `ASTER::math_backend` 即可取得相同条件。所有可执行文件放在构建目录的 `bin` 中。

```sh
cmake -S . -B build/openlibm -DCMAKE_BUILD_TYPE=Release
cmake --build build/openlibm -j 8
build/openlibm/bin/aster_iod_gl --version

cmake -S . -B build/system -DCMAKE_BUILD_TYPE=Release -DASTER_MATH_BACKEND=SYSTEM
cmake --build build/system -j 8
```

服务器追加 `-DASTER_BOOST_ROOT=/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0`，复用既有 Boost。Eigen 实际版本固定为 3.4.0，SOFA 保持 2018-01-30。

默认 `OPENLIBM` 使用已验证的官方 OpenLibm 0.8.8。来源、提交和归档 SHA-256 见 `third_party/openlibm-source.json`；许可见 `third_party/openlibm/LICENSE.md`，上游源文件及其 CMake 保持原样。OpenLibm 按其上游规则构建为邻接私有共享库，Linux 增加 `-Wl,-Bsymbolic-functions`，防止库内调用回到主程序入口。

公共 core、SOFA、主程序和数学入口使用 `-fno-builtin-sin -fno-builtin-cos -fno-fast-math -ffp-contract=off`，C++ 同时定义 `EIGEN_DONT_VECTORIZE`。只限制 sin/cos 的内建替换及合并，保留硬件 sqrt、常数幂与小内存复制优化。上游 OpenLibm 自身仍有其 `-fno-builtin`，该选项不传播给 ASTER。`SYSTEM` 使用原生数学及工具链默认浮点/Eigen 优化。两种后端都将 `__FILE__` 的源码路径映射为相对路径。各程序 `--version` 显示实际后端、构建类型和主要浮点条件。

标准 double 函数入口涵盖源码实际使用的 sin/cos、tan、asin/acos/atan/atan2、sqrt、pow、exp/log、hypot、floor/ceil、fmod/remainder、copysign、frexp/ldexp/scalbn、modf，以及常用其他指数、双曲和特殊函数。sin/cos 使用此前已验收的完整 binary64 键缓存，每函数每线程 65536 项，碰撞重新求值。矩阵四则运算仍为 double。float/long double 重载不属于本次统一入口，当前正式源码数值路径使用 double。

程序通过自身位置加载同目录的 `libopenlibm.dylib` 或 `libopenlibm.so`。移动程序时一起移动该库；不安装或替换系统 libm。链接 SYSTEM 的程序不加载 OpenLibm。

Ljung–Box 调用按既定决定注释保留，任务入口把四个遗留 Q/p 字段置为 −1。科学公式、收敛条件、迭代上限及按最终 J 选优的规则不变。

## 有限验收命令

先运行四个原有单 GL，再按需要运行完整 GLML；输入来自已保存的观测弧段。

```sh
python3 tests/math_backend_cases/run_checks.py --bin-dir build/openlibm/bin --output build/openlibm/acceptance/gl --only-gl
python3 tests/math_backend_cases/run_checks.py --bin-dir build/openlibm/bin --output build/openlibm/acceptance/iod
python3 tests/check_math_backend_apps.py --bin-dir build/openlibm/bin --output build/openlibm/acceptance/apps
ctest --test-dir build/openlibm -R '^(core_|math_backend_sofa)' --output-on-failure
```

四个 GL 和两个各含 10001 候选的 GLML 复用已验收科学输出的 SHA-256，要求精确相同。七个普通程序各用既有代表算例，比较历元、编号和文本，并按 `doc/linux-validation.json` 的原有物理量阈值验收：位置/半长轴 1 mm、速度 1 μm/s、偏心率/倾角参数 1e−10、角度 1e−7 度、角残差 1e−10 rad。脚本可用 `--native-summary` 指定此前同平台原生运行清单，读取其基线输出；默认与仓库 fixture 比较。NaN 必须对应 NaN，角度按 360 度周期比较。

旧 `iod_gl`、`iod_glml` 字节 fixture 保留，属于此前系统数学路径。`math_backend_iod` 仅为 OPENLIBM 注册。完整旧 CTest 不作为新后端通过标准；本次实际验收结果与性能记录另见 `doc/math-backend-validation.json`。
