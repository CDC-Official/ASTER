# 1.2.0 观测生成验收记录

本篇保存 1.2.0 发布时的迁移、几何、加噪及 Mac/Linux 对照结果。操作和输入格式见[观测生成指南](../06_observations.md)，来源及输入 SHA-256 见[source.json](../../tests/generate_observations/source.json)。文中构建目录是原始运行证据的位置。

## Mac 核对结果

- 原运行检查 239 项、独立几何检查 62 项和原加噪检查通过；OPENLIBM、SYSTEM 两种构建均通过。原算例的可见记录数在无筛选、黑夜、地影和联合筛选下分别为 364、106、306、48，解析和 JPL 太阳历表在该算例中一致。
- SYSTEM 与 CDC_lib 原入口的 19 份观测、完整观测和真值文件逐字节一致。相同理论输入、3″、seed 42、两组样本下，原/新加噪脚本的 12 份文件逐字节一致。
- 默认 OPENLIBM 与旧 SYSTEM 相比，全部历元最大角差为 3.98×10⁻¹¹ deg，Cartesian 真值位置分量差不超过 1.01×10⁻¹¹ km、速度分量差不超过 1.07×10⁻¹⁴ km/s；历元、ID、可见标记、分弧和 Keplerian 真值相同。独立站心矢量反算与生成角度的最大差为 5.69×10⁻¹⁴ deg。
- 理论值及 seed 42 的带噪资料分别通过现有 GL、GLML 的十条观测入口，四次运行均写出有限轨道与十条残差；原算例单次生成约 0.1 s。结果用于功能和迁移核对。

加噪派生种子包含源文件内容摘要。因此切换数学后端造成的理论文件末位变化，也会改变派生随机流；逐样本复现需要相同输入字节、脚本、Python 环境及 seed。详细本地结果位于 `build/v1.2.0/acceptance/observations-final/summary.json` 和 `build/v1.2.0/acceptance/migration-vs-cdc-lib/migration-summary.json`。

## Linux 与跨平台核对

Linux x86_64 / GCC 13.3 的原观测运行、几何、加噪和坐标/SOFA 五项 CTest 全部通过，四次 GL/GLML 资料交接也通过。Mac arm64 / AppleClang 21 与 Linux 共用 OpenLibm 0.8.8、Eigen 3.4.0 和同一原 LEO 输入，输出比较如下：

| 输出 | 文件数 | Mac/Linux 逐字节一致 |
| --- | --- | --- |
| 理论 OBS、完整观测及两类真值 | 19 | 19 |
| 3″、seed 42、两组带噪样本 | 12 | 12 |
| 理论/带噪资料的 GL、GLML 轨道解和残差 | 8 | 8 |

本轮加噪实际使用 Mac Python 3.9.6、Linux Python 3.12.3；结论限于这组已比较的输入与环境。两端内部重复生成、反转加噪输入顺序均保持文件一致，四次求解的收敛标志相同；生成原算例的实际单次耗时约 0.12—0.15 s。

Linux 还以字节不变的原 `main_arcs.cpp`，通过四行头文件适配链接当前 ASTER 公共 core/math_backend，进行同后端的旧新入口受控比较，19 份输出逐字节一致。该对照只在构建目录内编译原主程序。Mac 使用旧 CDC_lib 原二进制的独立 SYSTEM 对照见上一节。

跨平台哈希、Linux 测试日志及旧入口摘要保存在 `build/v1.2.0/acceptance/cross-platform/`，其中 `comparison.json` 为逐类比较结果。
