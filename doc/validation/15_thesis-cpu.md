# 论文实验线的 GLML CPU 优化

本实验线为 `experiment/thesis`，从正式 v1.4.0 的
`f2c2a8228040f3b5d9b86c3c02b4ec6b7e13435c` 创建，不修改 main、版本号和发布标签。

GLML 候选保留完整数值迭代记录，选出最优候选后生成其来源报告；日月距离按本次任务的观测历元一次计算，
以只读数组交给候选求解、方程构造和边界后的残差评价。GL 的默认调用方式与即时报告继续保留。

专项检查 `core_iod_reuse` 覆盖五种动力学、简单迭代与微分改正、负/零/正相对时间、重复历元，
逐位比较传播系数与数值历史，核对报告重建、显式流与延后留存，并检查无效模型、非有限时间和数组长度错误。
完整候选与耗时对照使用仓库保存的 10 条敏感短弧、e070 的 20 条短弧以及 e010 的 100 条较长弧。

## 复现入口

在两个已提交的源码工作树分别构建相同的 Release / OpenLibm 配置，再运行：

```sh
python3 tests/thesis/run_comparison.py \
    --baseline-source /home/cdc/projects/ASTER \
    --baseline-bin /home/cdc/projects/ASTER/build/v1.4.0/bin \
    --candidate-source /home/cdc/projects/ASTER/build/experiment-thesis-src \
    --candidate-bin /home/cdc/projects/ASTER/build/experiment-thesis/bin \
    --output /home/cdc/projects/ASTER/build/thesis-validation-20261105 \
    --threads 1 32 64 --repeats 5
```

输出目录须为新目录。每组先各运行一次预热，计时阶段交错运行两个独立进程；预热单独保存并排除。
精确墙钟包括进程启动、解析、GL、10,000 个候选、边界、评价、选优、科学输出及退出。
每次保留 setting、输入/输出 SHA-256、控制台、GNU time、所选候选及不含机器路径和耗时的迭代报告。
比较三份科学文件（全候选表、最终轨道、逐条残差）和基准/所选来源的逐轮报告；完整任务报告另含路径、
源码行号、耗时和内存，这些运行信息不做整文件字节比较。专项检查验证其数值轨迹等价。

## 实际构建与验收

实测基线为 v1.4.0 `f2c2a8228040f3b5d9b86c3c02b4ec6b7e13435c`，优化代码提交为
`19bab2408ee54519eb78f7646079cb292340e91d`。后续仅补充本页及验证证据的提交不改变所测计算源码；
实际使用时仍以 `git rev-parse HEAD` 保存分支完整提交，不能仅凭 `--version` 的 1.4.0 区分。
分支保留 v1.4.0 版本号，不建立新的发布标签。

服务器为双路 Xeon Gold 6430、64 个物理核 / 128 个逻辑线程、Linux x86_64；
GCC 13.3、CMake 3.28.3、Release（`-O3 -DNDEBUG`）、Boost 1.87.0、Eigen 3.4.0、SOFA 2018-01-30、
OpenLibm 0.8.8。两版均使用 `-fno-builtin-sin -fno-builtin-cos -fno-fast-math -ffp-contract=off`
和 `EIGEN_DONT_VECTORIZE`。线程绑定为 `OMP_DYNAMIC=FALSE OMP_PROC_BIND=spread OMP_PLACES=cores`，
并按组固定 `OMP_NUM_THREADS`。完整环境、二进制及库校验见[原始 JSON](evidence/21_thesis-cpu.json)。

服务器全目标构建及 23 项当前回归通过，含本次 `core_iod_reuse`；Mac 的 GL、GLML、专项检查也已编译，
专项检查通过。这里不将历史版本的字节参考用作当前版本的物理验证。

## 完整任务耗时

每格为 5 次独立完整进程墙钟的中位数（s）；每组各一次预热不计入。相同线程数逐组比较，
前后顺序交替，不与回归或其他本次算例并行运行。逐次耗时及范围均保存在 JSON 中。

| 算例与跨度 | OpenMP 线程 | v1.4.0 / s | 论文实验线 / s | 加速比 | 两版所选候选 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 敏感 10 条 / 27 s | 1 | 4.059 | 2.042 | 1.99× | 379 |
| 敏感 10 条 / 27 s | 32 | 0.339 | 0.278 | 1.22× | 379 |
| 敏感 10 条 / 27 s | 64 | 0.347 | 0.321 | 1.08× | 379 |
| e070 20 条 / 57 s | 1 | 15.304 | 8.127 | 1.88× | 7156 |
| e070 20 条 / 57 s | 32 | 0.714 | 0.467 | 1.53× | 7156 |
| e070 20 条 / 57 s | 64 | 0.557 | 0.408 | 1.36× | 7156 |
| e010 100 条 / 297 s | 1 | 31.713 | 18.547 | 1.71× | 0 |
| e010 100 条 / 297 s | 32 | 1.341 | 0.834 | 1.61× | 0 |
| e010 100 条 / 297 s | 64 | 0.918 | 0.637 | 1.44× | 0 |

## 数值与报告对照

三个输入、三种线程数、两版各一次预热及五次计时，共 108 次完整 GLML 运行：
每次的 10,001 行候选表、最终轨道及残差均与对应基线逐字节一致，所选候选编号及来源 flag 一致。
候选表包含全部根数、Cartesian、RMS、J 和来源状态，因此未收敛候选也在本次比较中。
基准和所选来源的逐轮报告文字均逐字节一致；保留的整份报告可追溯各次路径、事件、时间与内存。

`core_iod_reuse` 对五模型、两种迭代方式逐条核对完整数值轨迹，另比较 80 组传播系数。
完整 GLML 对照核对全部候选最终输出与实际展示的迭代文字，未另导出全部 10,001 条候选的完整逐轮二进制历史。

验证范围为已保存的论文代表输入及核心专项检查，未遍历全部论文蒙特卡洛批次；
“较长弧”在此指 100 条、297 s 的相对较长输入。速度是本机本次构建下的完整任务结果。

## 服务器保存与后续使用

- 原 main 源码：`/home/cdc/projects/ASTER`，保持 v1.4.0；原发布构建继续位于 `build/v1.4.0`。
- 分支工作树：`/home/cdc/projects/ASTER/build/experiment-thesis-src`。
- 论文构建：`/home/cdc/projects/ASTER/build/experiment-thesis`，程序在其中的 `bin/`。
- 全部运行证据：`/home/cdc/projects/ASTER/build/thesis-validation-20261105`。
- 构建/回归/计时日志：原仓库 `build/experiment-thesis-*.log`；构建目录 `environment.json` 保存实际环境。
- GitHub：[experiment/thesis](https://github.com/CDC-Official/ASTER/tree/experiment/thesis)，与 main 并列保存。

在服务器重新构建：

```sh
cmake -S /home/cdc/projects/ASTER/build/experiment-thesis-src \
    -B /home/cdc/projects/ASTER/build/experiment-thesis \
    -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DASTER_MATH_BACKEND=OPENLIBM \
    -DASTER_BOOST_ROOT=/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0
cmake --build /home/cdc/projects/ASTER/build/experiment-thesis --parallel 3
```

后续实验把已有配置复制到新的实验目录，再调用此构建：

```sh
git -C /home/cdc/projects/ASTER/build/experiment-thesis-src rev-parse HEAD
OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OMP_PROC_BIND=spread OMP_PLACES=cores \
    /home/cdc/projects/ASTER/build/experiment-thesis/bin/aster_iod_glml -I /absolute/path/to/experiment/setting.txt
```

保存实际 setting、输入 SHA-256、求解/加噪种子、线程数和完整 Git 提交；既有配置无需新增优化开关。
