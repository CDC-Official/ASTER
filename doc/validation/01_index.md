# 科学验证资料索引

本目录保存后续实验仍会用到的独有科学验证资料：GLML 迭代敏感性、数学后端的数值与耗时比较、观测生成的迁移和几何检查，以及对应原始 JSON。每份记录标明所用版本和计算条件。日常操作从[文档首页](../01_overview.md)进入。

## 按问题查证据

| 要核对什么 | 说明 | 原始证据 |
| --- | --- | --- |
| 1.0.0 源码从哪里来，迁移时怎样验证 | [迁移来源与核对](#migration) | [文件映射和 SHA-256](evidence/01_source-map.json)、[迁移验收](evidence/02_migration-validation.json) |
| 首次 Linux 数值对照及既有门限 | [早期 Linux 对照](#linux) | [Linux 原始结果](evidence/03_linux-validation.json) |
| 原生数学下的 GLML 短弧为何选出不同候选 | [GLML 敏感性诊断](02_glml-sensitivity.md) | [编译器、线程与扰动](evidence/04_clang-glml-validation.json)、[收敛与控制试验](evidence/05_glml-convergence.json) |
| 1.1.0 统一后端的结果与成本 | [后端数值和计时](03_math-backend.md) | [构建、科学差异和计时](evidence/06_math-backend-validation.json) |
| 1.2.0 观测生成是否保持原行为 | [观测生成验收](04_observations.md) | [迁移来源与输入校验](../../tests/generate_observations/source.json)；完整运行证据位置见说明 |
| 地球 GM 如何按 TT/TDB 取值，以及对轨道的影响 | [GM 约定与轨道对照](05_earth-gm.md) | [参数、输入校验和数值结果](evidence/07_earth-gm-validation.json) |

<a id="migration"></a>
## 1.0.0 迁移来源与核对

ASTER 从 CDC_lib 实际工作区整理而来，来源提交起点为 `96dbe9d`，并包含该提交之后作者与 Agent 已完成的工作区修改，以及迁移前验证的 GLML 角度和随机种子修复。文件对应和源码校验值保存在 `01_source-map.json`。

迁移时建立九个主程序与公共核心，保留输入关键词、科学单位、计算流程和输出职责。独立构建核对提交 `de55d5150f96af2f13ca0b4aa9a00475c789135f` 的 19 项 CTest 通过，九程序的 16 份科学结果逐字节一致。具体工具链、输入错误检查、数据大小等保存在 `02_migration-validation.json`。观测生成在 1.2.0 迁入，来源另见对应记录。

<a id="linux"></a>
## 早期 Linux 对照

`03_linux-validation.json` 对应 1.0.0 后的 `4b58175` 阶段，记录 GCC 13.3 的既有数值门限和逐文件结果。该阶段的标准指针声明、单字符分隔符、同平台空白输入对照修复分别由 `d2cb63d`、`b4713fc`、`4b58175` 追溯。后续 Clang/OpenMP 对照与短弧候选控制归入 GLML 诊断。

旧原生数学下，GLML 的 10 条短弧没有通过 Mac/Linux 原门限；1.1.0 统一 OPENLIBM 及浮点设置后，已检查代表算例的两端科学输出逐字节一致。两阶段分别读取各自记录，当前构建配置见[数学后端指南](../08_math-backend.md)。

## 旧指南和重复报告放在哪里

旧构建、部署操作和迁移过程保存在 Git 历史。当前有效的使用说明已进入编号指南，GLML 两份关联诊断合为一篇；版本变更统一在根 [CHANGELOG](../../CHANGELOG.md)。需要查看重组前的完整说明时，可在仓库根目录执行：

```sh
git show d138c98:doc/migration.md
git show d138c98:doc/server.md
```

`evidence/` 的六份 JSON 只调整目录和文件名，内容及 SHA-256 保持不变。它们内部的提交和机器路径记录当时的证据位置；本索引给出现在的对应入口。
