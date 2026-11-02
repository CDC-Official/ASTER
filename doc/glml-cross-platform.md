# Linux Clang 对照与 GLML 短弧差异定位

2026-11-02 在 `cdc-114` 完成 Clang 18 对照。Mac 与 Linux 的源码提交均为 `4e0e1f4cb72b7ada6c659eb6fd9a405fd2b53948`，Eigen 为 3.4.0，SOFA 为 2018-01-30。本轮保留科研源码、候选接受规则与已有数值门限。

Linux 的 Clang 和 GCC 在九个代表算例中生成的 16 份科学输出逐字节一致。Clang 的七项核心检查、数据完整性、输入契约及初轨到批处理接口检查也全部通过。与保存的 Mac 参考结果比较，八个程序通过数值检查，10 条观测的 GLML 短弧算例仍未通过。

后续的停止判据、两个目标函数及 40／80／160 轮控制试验见 [GLML 收敛诊断](glml-convergence.md)。该检查保留了最终选择规则，并进一步说明为什么仅筛已收敛候选会改变科学结果。

## 环境与 OpenMP

| 项目 | Mac | Linux |
| --- | --- | --- |
| 架构 | arm64 | x86_64 |
| 编译器 | AppleClang 21.0.0 | GCC 13.3.0、Clang 18.1.3 |
| 构建 | `build/release` | `build/linux-gcc-release`、`build/linux-clang-release` |
| 优化 | Release，`-O3 -DNDEBUG` | Release，`-O3 -DNDEBUG` |
| OpenMP | Homebrew libomp | GCC 使用系统 libgomp；Clang 使用用户目录 libomp 18 |
| 候选运行设置 | 种子 0；分别使用 1、2 线程 | 种子 0；分别使用 1、2 线程 |

Clang 所缺的 OpenMP 从 Ubuntu 正式发行源取得，包名为 `libomp-18-dev_18.1.3-1ubuntu1_amd64.deb` 与 `libomp5-18_18.1.3-1ubuntu1_amd64.deb`，版本均为 `1:18.1.3-1ubuntu1`，与服务器 Clang 匹配。下载总计 1,576,106 字节；Mac 和服务器分别按服务器已缓存的 APT 元数据核对 SHA256，均通过。

包来自 [Ubuntu 开发包](https://archive.ubuntu.com/ubuntu/pool/universe/l/llvm-toolchain-18/libomp-18-dev_18.1.3-1ubuntu1_amd64.deb)和 [Ubuntu 运行库包](https://archive.ubuntu.com/ubuntu/pool/universe/l/llvm-toolchain-18/libomp5-18_18.1.3-1ubuntu1_amd64.deb)。包及来源清单位于服务器实验目录的 `dependencies/packages/`，通过 `dpkg-deb -x` 解包到 `dependencies/openmp18/`。Clang 程序实际加载该目录的 `libomp.so.5`，CMake 的 OpenMP 动态调度链接检查通过。整个过程使用用户目录。

## 编译器、采样与线程对照

测试复用 `tests/run_case.py` 的算例准备函数，保持原始 10 条观测、动力学模型 5、迭代方法 2、基准约束类型 1、40 轮上限与 `1e-10` 状态容差。固定偏心率候选数仍为 10000，另有一个基准候选。运行时显式设置 `random_seed=0`、`OMP_DYNAMIC=FALSE`，并限制实际线程数。

为取得边界处理前的来源轨道和逐轮状态，在实验目录生成独立的主程序诊断副本，链接各自现有的核心库。诊断副本仅增加原始状态输出和可选的初值回读。其候选表、残差与最终解在 Mac、Linux GCC、Linux Clang 上均与各自正式程序逐字节一致。

| 检查 | 结果 |
| --- | --- |
| 同机 GCC / Clang 的九个算例 | 16 份科学输出逐字节一致 |
| 同机 GCC / Clang 的 GLML 来源轨道和抽查迭代记录 | 逐字节一致 |
| Mac 的 1 / 2 线程 GLML 输出 | 逐字节一致 |
| Linux GCC 的 1 / 2 线程 GLML 输出 | 逐字节一致 |
| Linux Clang 的 1 / 2 线程 GLML 输出 | 逐字节一致 |
| 三个平台的基准 GL | 均在第 16 轮收敛，未触发随机初值回退 |
| 三个平台的全部 10001 个候选 | 边界处理中的随机抽样次数均为 0 |

因此，本算例的差异不来自 GCC 与 Clang 的切换、随机数分布实现或 OpenMP 调度。原始 RA、Dec、观测时刻和归一化时间差完全相同；准备后的测站位置分量最大相差 `1.99e-14` 个归一化长度单位，基准轨道的位置相差约 `1.65e-5 m`、速度相差约 `7.41e-8 m/s`。

## 未收敛迭代如何放大差异

两端均有 8593 个来源解收敛、1408 个达到迭代上限，数值失败数为 0；这些计数包含基准解。当前规则按边界处理后的最小有限 J 选择候选，因此达到迭代上限的轨道仍可胜出。

| 项目 | Mac | Linux GCC / Clang |
| --- | --- | --- |
| 所选候选 | 592 | 335 |
| 来源状态 | -1，达到迭代上限 | -1，达到迭代上限 |
| 来源轮数 | 40 | 40 |
| 指定偏心率 | 0.04728 | 0.02672 |
| 来源实际偏心率 | 0.1134011362 | 0.1130575476 |
| 末轮最大归一化状态改正 | 0.0651161 | 0.0725682 |
| J | -198.0177118202 | -198.0172445332 |

最终所选轨道的位置相差 `15387.3352 m`，速度相差 `31.0198138 m/s`。两端所选来源的末轮改正量远大于 `1e-10` 的停止门限。

对同一个候选 592，微小初值差异在反复迭代中逐渐放大：

| 轮数 | Mac / Linux 的同一候选位置差 |
| --- | --- |
| 0 | 0.0000165 m |
| 10 | 0.00476 m |
| 20 | 0.260 m |
| 30 | 444.2 m |
| 39 | 260.27 km |
| 40 | 43.29 km |

候选 335 也未收敛，其第 40 轮在两端相差约 339.75 m。J 排序随后在这些不同的末轮轨道中选取最优项，导致两端最终采用不同候选。

同机控制进一步验证了这种敏感性。保持 Linux GCC、观测、种子、线程与全部设置不变，只把初值 x 分量向正无穷方向移动一个 double 的末位，输入位置仅改变 `1.77e-10 m`，最终仍选候选 335，但位置改变 `1793.3555 m`、速度改变 `3.61147 m/s`。原样回读同一初值文件时，三份科学输出逐字节不变。将候选初值改为 Mac 保存的初值后，Linux 改选候选 596，最终位置相差 `10938.8454 m`。

这些结果将问题定位到本短弧中未收敛候选对微小浮点差异的放大，以及这些候选参与最小 J 选择的规则；未进一步归因到某一条浮点指令或数学库函数。

## 已收敛候选的诊断对照

从同一批已保存候选表中，仅在来源状态为正数的候选内选择最小 J，Mac、Linux GCC、Linux Clang 及两组初值扰动均选中候选 1447。Mac 与 Linux 的位置差为 `2.61e-5 m`，速度差为 `7.64e-8 m/s`，J 差为 `1.95e-9`，均在原有门限内；两端来源解都在第 15 轮收敛。

这项筛选仅保留为诊断对照，没有写入正式程序。其最小 J 约为 -197.9857569，比原规则所选轨道的 J 略高。作者已确认允许保留迭代上限处的轨道，在原有其他有效性要求下按最终 J 选优；不采用来源必须收敛的筛选。完整决定见 [GLML 收敛诊断](glml-convergence.md)。

Linux GCC / Clang 同机对照已通过；Mac/Linux 的 10 条观测 GLML 短弧在原门限下仍未通过，数值敏感性尚未解决。数值门限继续采用位置 1 mm、速度 `1e-6 m/s`、角残差 `1e-10 rad` 和 J `1e-6`。

## 原始记录与复现

机器可读结果见 [clang-glml-validation.json](clang-glml-validation.json)。源码及输入校验值、OpenMP 包校验值、16 份输出的逐项比较、原始来源轨道、逐轮状态和扰动对照均已保存。服务器实验目录为 `/home/cdc/projects/aster-clang-20261102`；Mac 实验目录为 `/Users/chendingcun/Documents/Codex/2026-11-02/task-2/diagnostics/20261102`。

服务器从固定源码重建并复现：

```sh
ssh -o ClearAllForwardings=yes cdc-114
cd /home/cdc/projects/aster-clang-20261102
bash reproduce_server.sh
```

其中 Clang 使用独立构建目录，编译并行度为 2。脚本复用已校验并解包的用户目录 OpenMP。单独重跑九个代表算例及同机比较：

```sh
python3 run_cases.py --project /home/cdc/projects/ASTER --bin /home/cdc/projects/ASTER/build/linux-gcc-release/bin --output repeat-gcc
python3 run_cases.py --project /home/cdc/projects/ASTER --bin /home/cdc/projects/ASTER/build/linux-clang-release/bin --output repeat-clang
python3 compare_cases.py --reference repeat-gcc --actual repeat-clang --output repeat-gcc-vs-clang.json
```

与 Mac 参考结果比较时，GLML 仍应报告失败，并返回状态 1：

```sh
python3 compare_cases.py --reference /home/cdc/projects/ASTER/tests --actual repeat-clang --output repeat-mac-vs-clang.json
```

Mac 上重新汇总已保存的完整诊断记录：

```sh
cd /Users/chendingcun/Documents/Codex/2026-11-02/task-2
python3 diagnostics/20261102/analyze_glml.py --root diagnostics/20261102 --output diagnostics/20261102/glml-diagnostics.json
```
