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


## 九程序验收结果（2026-11-02）

全库公共后端的跨平台一致性目标通过：九个程序在相同输入、参数和线程数下，OPENLIBM 的科学输出在 Mac arm64 / AppleClang 21 与 Linux x86_64 / GCC 13.3 间全部逐字节一致，数值差均为 0。四个单 GL 的逐轮状态、停止分支及 Q/p=−1 检查通过；两个完整 GLML 样本 e010、e070 的各 10001 行候选、最终解与残差均复用既有验收哈希通过。两端各 8 项核心/SOFA 检查通过，SYSTEM 的坐标/SOFA 各 2 项通过。全套旧 CTest 未作为通过标准。

下表的同平台差异使用既有代表算例，Δr、Δv 为位置、速度的最大绝对分量，Δres 为逐条角残差最大差，Δa 为半长轴差，Δ角按 360° 周期计算。位置/半长轴 1 mm、速度 1 μm/s、角残差 1e−10 rad、角度 1e−7°、偏心率/倾角参数 1e−10、J 1e−6 的原门限保持不变。完整分量及文件检查见 [机器可读结果](math-backend-validation.json)。

| 程序 | OPENLIBM Mac↔Linux | Mac SYSTEM→OPENLIBM 科学差异 | Linux SYSTEM→OPENLIBM 科学差异 | 同平台近似门限 |
|---|---|---|---|---|
| GL (`iod_gl`) | 0，字节一致 | Δr 21208.6 m；Δv 10.6849 m/s；Δres 2.25e-08 rad | Δr 3198.92 m；Δv 1.61193 m/s；Δres 3.4e-09 rad | 两端均超限 |
| GLML (`iod_glml`) | 0，字节一致 | Δr 81887.5 m；Δv 41.275 m/s；Δres 8.77e-08 rad；ΔJ 0.000348 | Δr 28030.7 m；Δv 14.1414 m/s；Δres 3.06e-08 rad；ΔJ 0.000158 | 两端均超限 |
| BODN (`odn_batch`) | 0，字节一致 | 角残差 3.48e-15 rad | 角残差 5.69e-15 rad | 两端均通过 |
| Cartesian 预报 (`opn_cartesian`) | 0，字节一致 | Δr 1.82e-09 m；Δv 2.22e-12 m/s；Δa 6.37e-09 m；Δ角 1.92e-12 ° | Δr 0 m；Δv 0 m/s；Δa 0 m；Δ角 4.97e-14 ° | 两端均通过 |
| NSG2 预报 (`opn_nonsingular2`) | 0，字节一致 | Δr 5e-09 m；Δv 4.88e-12 m/s；Δa 0 m；Δ角 6.11e-13 ° | Δr 1.82e-09 m；Δv 2.66e-12 m/s；Δa 0 m；Δ角 5.15e-14 ° | 两端均通过 |
| 地月预报 (`opn_cislunar`) | 0，字节一致 | Δr 0 m；Δv 0 m/s | Δr 0 m；Δv 0 m/s | 两端均通过 |
| 瞬时根数→平根数 (`convert_osculating_to_mean`) | 0，字节一致 | Δa 1.82e-09 m；Δ角 9.35e-12 ° | Δa 0 m；Δ角 1.02e-13 ° | 两端均通过 |
| 根数→Cartesian (`convert_keplerian_to_cartesian`) | 0，字节一致 | Δr 2.33e-10 m；Δv 0 m/s | Δr 0 m；Δv 0 m/s | 两端均通过 |
| Cartesian→根数 (`convert_cartesian_to_keplerian`) | 0，字节一致 | Δa 0 m；Δ角 3.41e-13 ° | Δa 0 m；Δ角 0 ° | 两端均通过 |

七个普通程序均在原门限内，偏心率/倾角参数最大差也低于 6.11e−16。e070 GL 与 GLML 的原生后端对照有明显轨道变化，不能按“只有末位变化”验收；两种后端的来源求解均为迭代上限 flag=−1。完整 GLML 在两平台、两后端都选中候选 5190，已核对最终 J 确为全表有限 J 的最小值。Mac 的最终 GLML 位置向量差为 97.301 km，Linux 为 33.307 km；单 GL 分别为 25.201 km、3.801 km。两种原生配置切换到统一浮点配置后，敏感迭代沿不同数值轨迹推进，既有残差/J、迭代上限候选参与选优的规则均保持。当前通过的是统一后端后的平台复现目标；未开展新的全仿真真值精度统计，尚不能据此认定绝对定轨精度提高。

GLML 的原生→OPENLIBM 全候选表也有显著变化：Mac/Linux 分别有 297/200 行出现有限值与非有限值差异，来源 flag 分别改变 2/1 行；共同有限状态的最大位置分量差分别为 25533.7/20890.4 km。OPENLIBM 两平台的全表中这些差异均为 0。

## 同平台耗时对照

两端均为 Release、`-O3 -DNDEBUG`，Eigen 3.4.0、SOFA 2018-01-30、Boost 1.87.0。SYSTEM 保留工具链默认 FMA/Eigen/builtin；OPENLIBM 使用上文统一条件。因此耗时是整套后端配置切换的成本。GL 为 1 线程，GLML 为 4 线程，其余为 2 线程，`OMP_DYNAMIC=FALSE`。三次交错运行取中位数，包含进程启动、动态库加载、计算及输出。输入与参数指纹相同；每次科学输出哈希都与已有结果相同，两平台各 54 次全部通过。

| 程序 | Mac SYSTEM / s | Mac OPENLIBM / s | 变化 | Linux SYSTEM / s | Linux OPENLIBM / s | 变化 |
|---|---:|---:|---:|---:|---:|---:|
| GL | 0.061289 | 0.064232 | +4.8% | 0.044677 | 0.050667 | +13.4% |
| GLML | 5.063107 | 6.371624 | +25.8% | 4.108581 | 3.995077 | -2.8% |
| BODN | 0.410261 | 0.460386 | +12.2% | 0.757422 | 0.791226 | +4.5% |
| Cartesian 预报 | 0.093699 | 0.110721 | +18.2% | 0.159972 | 0.351709 | +119.9% |
| NSG2 预报 | 0.095759 | 0.108558 | +13.4% | 0.160245 | 0.348806 | +117.7% |
| 地月预报 | 0.058542 | 0.061681 | +5.4% | 0.051283 | 0.057889 | +12.9% |
| 瞬时根数→平根数 | 0.018115 | 0.048909 | +170.0% | 0.032419 | 0.075266 | +132.2% |
| 根数→Cartesian | 0.002840 | 0.003024 | +6.5% | 0.003024 | 0.004396 | +45.4% |
| Cartesian→根数 | 0.002803 | 0.003081 | +9.9% | 0.001197 | 0.002672 | +123.3% |

本次完整 e070 GLML 在 Mac 增加 25.8%，Linux 减少 2.8%；Linux 两个地球预报程序增加约 118%—120%，瞬时根数转平根数在两端也有明显成本。短程序受启动、缓存与系统负载影响，尤其毫秒级转换程序的百分比波动较大；这些数字描述本次代表输入，未测长时预报或全库总体吞吐量。旧测时含进程超时轮询造成的量化误差，本表采用补测的管道计时结果。

配对测时驱动为本任务工作目录中的 `measure_backends.py`，复用已有输入与输出证据；原始三次秒数、可执行文件 SHA-256、实际编译命令、链接路由审计及差异明细均保存在机器可读结果中。源码与输入已核对，不改变科学公式和验收门限。开发分支保留供用户验收，未合并 main、未 push、未创建 v1.1.0 标签。
