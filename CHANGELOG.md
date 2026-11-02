# 版本记录

## 1.2.0 — 2026-11-02

- 从 CDC_lib 当前实际入口迁入 `aster_generate_observations`：读取 Cartesian/Keplerian 星历和测站，逐历元使用 EOP，执行高度角、可选黑夜/地影筛选，生成理论 RA/Dec 弧段、同历元真值及可选完整观测。输入单位、ID 分弧和输出格式保持原定义。
- 迁入 `tools/add_angular_noise.py`，保留东西/南北等标准差小角度模型、标准库随机流、原种子派生与来源记录。原 NumPy 论文批次使用另一随机流，seed 不作样本对应。
- 新入口接入公共 OPENLIBM/SYSTEM 数学后端及版本查询；CMake 正式版本为 1.2.0。沿用 ASTER 的显式 `-I` 配置入口，替代旧程序的本机默认路径；弧段写出进入公共 `include/`、`src/`。
- 复用原 LEO 算例与运行、几何、加噪验收，补充理论/带噪观测至 GL、GLML 的实际交接检查；用法见 [观测生成](doc/observations.md)。
- Mac/Linux 原观测检查和 GL/GLML 资料交接通过；OPENLIBM 下 19 份理论/真值、12 份固定 seed 带噪文件及 8 份定轨科学输出逐字节一致。Mac SYSTEM 与旧库原二进制的 19 份输出一致，Linux 同后端旧/新入口受控比较的 19 份输出一致；范围和末位差异见 [观测生成](doc/observations.md)。
- 发布与维护：两端 main 使用同一套源码和同一发布提交；v1.0.0、v1.1.0 标签原指向保留。main 为唯一长期分支，feature/、fix/ 分支在验收合并后清理；未合并的独有修改保留。

## 1.1.0 — 2026-11-02

- **公共数学后端**：九程序、公共 core 与 SOFA 共用 `math_backend`，通过 `ASTER_MATH_BACKEND=OPENLIBM/SYSTEM` 选择，默认官方 OpenLibm 0.8.8；`--version` 显示版本、后端和构建条件。
- **计算规则**：Ljung–Box 调用注释保留，未计算 Q/p 为 −1。科学公式、残差/J、原门限、迭代上限及按最小 J 选优的规则保持不变。
- **平台与依赖**：已验证 Mac arm64 / AppleClang 21、Linux x86_64 / GCC 13.3，Boost 1.87.0 与现有 OpenMP。需要 CMake ≥3.26、C++17；Eigen 3.4.0、SOFA 2018-01-30、OpenLibm 0.8.8 固定源码及许可随库保存。
- **验证范围**：代表输入下九程序 OPENLIBM 的 Mac/Linux 科学输出逐字节一致；七个非 IOD 程序同平台前后差异满足原物理门限。两端各 8 项核心/SOFA、4 个单 GL 逐轮检查、2 个各 10001 候选的完整 GLML 通过，SYSTEM 坐标/SOFA 各 2 项通过。敏感 GL/GLML 的轨迹与候选可能改变。e010 四样本中普通 GL 和三组 GLML 基本不变，3″、27秒短弧 GLML 换候选后误差下降，但新来源仍达 40 轮上限，位置误差约 93.720 km；这些结果不推广为所有场景更准确。
- **效率备忘**：九程序已有计时保留在[后端验收](doc/math-backend.md)。完整 e070 GLML 的 Mac 中位耗时 5.063→6.372 s（+25.8%），Linux 4.109→3.995 s（−2.8%）；Linux 两个地球预报增加约 118%—120%，瞬时根数转平根数的 Mac/Linux 增加约 170%/132%。毫秒级程序百分比噪声较大。后续优先分析这些热点，并继续保持现有数值验收门限；本版不追加性能试验。
- **版本管理**：每台机器日常只用原 ASTER 仓库，通过标签和对应构建目录切换版本。`v1.0.0` 保持原提交，`v1.1.0` 为正式 annotated 标签。历史试验保存在本地备份/stash，不混入 main；GitHub 自动提供源码归档，本版不附未经验证的安装器或二进制包。

默认 OPENLIBM 构建：

```sh
cmake -S . -B build/v1.1.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.1.0 --parallel 3
build/v1.1.0/bin/aster_iod_gl --version
```

SYSTEM 使用单独构建目录，并添加 `-DASTER_MATH_BACKEND=SYSTEM`。移动 OPENLIBM 程序时，一起移动同目录私有 `libopenlibm.dylib` / `libopenlibm.so`。测试脚本现从 CMake 读取项目版本；旧 1.0.0 字节 fixture 保留，完整旧 CTest 不作为数学路径切换后的通过标准，当前命令见[测试说明](tests/README.md)。

## 1.0.0 — 2026-11-01

建立独立 ASTER 工程，收录九个轨道研究程序。GL、GLML 采用已整理的实验版实现，预报和批处理保留当前科研流程与输出规则。状态转换的短输入检查在主程序中展开。

统一构建入口、程序名称与三段版本号，最小算例集中到 `tests/`。JPL、EOP、月球重力场和大气数据随 `data/` 保存，记录来源和校验值。Eigen 固定 3.4.0，SOFA 固定 2018-01-30。

GLML 包含迁移前已验证的弧度回退修复，以及可记录、可复现的 `random_seed`。源码来源和实际验收结果见 `doc/migration.md`。
