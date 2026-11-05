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

服务器构建和实测结果将在验收后记录于本页及对应 JSON。
