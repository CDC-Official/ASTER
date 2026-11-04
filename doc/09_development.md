# 开发与 Git 维护

日常只保留一份 ASTER 源码，版本交给 Git，编译与实验结果交给 `build/`。准备改动时先找到对应主程序、配置和测试，了解当前科学行为后再修改。

## 源码各放在哪里

| 目录 | 职责 |
| --- | --- |
| `apps/` | 十一个主程序，按读资料、准备状态、调用计算、写结果展开 |
| `include/` / `src/` | 公共声明和科学实现；同名文件通常对应同一模块 |
| `cmake/` | 固定依赖、OpenMP、SOFA 和数学后端接入 |
| `third_party/` | Eigen 3.4.0、SOFA 2018-01-30、OpenLibm 0.8.8 源码及许可 |
| `tests/` | 小输入、参考科学结果和检查程序 |
| `data/` | 共享科学数据及来源/校验清单 |
| `tools/` | 当前测角加噪工具 |
| `doc/` | 按阅读顺序编号的使用与接口指南；独有科学验证资料放 `doc/validation/` |
| `build/` | 本机编译、实际运行、检查与实验输出，Git 忽略 |

公共模块的入口从[接口参考](07_interfaces.md)查。`core_iod_solver` 负责 GL 求解，GLML 的候选生成、边界处理和筛选留在主程序；`core_propagation` 负责整弧积分与插值；`core_bodn` 负责批处理；`analyt_averaging` 负责瞬转平。

新增共享科学接口时，同时说明参数方向、单位、参考系、时间尺度、失败方式和结果用途。涉及公式、单位、选择规则或模型变化时，单独说明依据，并执行相应数值检查。具体维护规则见 [AGENTS.md](../AGENTS.md)。

## 选哪些检查

[tests/README.md](../tests/README.md)列出实际命令。先做当前改动需要的检查：

| 改动范围 | 对应检查 |
| --- | --- |
| 文档、配置示例、路径 | 本地链接、输入路径和示例命令；涉及 C++ 示例时编译小示例 |
| 轨道状态、坐标、瞬转平 | 对应 `core_...` 检查及相关应用算例 |
| 积分、观测模型、BODN | 积分安全、批处理核心、输入与资料接续检查 |
| 数学后端或依赖 | 统一后端的单 GL、完整 GLML 与普通应用比较，按变动范围执行 |
| 观测生成或加噪 | `observation_...`、固定种子复现及 GL/GLML 接续检查 |

测试输入与参考保存在 `tests/`，运行结果保存在构建目录。增加新的科学数据时，同时维护 `data/manifest.json` 的来源、大小和 SHA-256。

## 先理解几个 Git 名称

- **提交（commit）**：某次修改完成后的源码记录，有唯一编号，可以回看。
- **分支（branch）**：指向某条开发进度的名称。切换分支是在同一个目录选择另一条历史。
- **main**：唯一长期主分支，保存验收后的工作。
- **feature/...、fix/...**：临时功能或修复分支，完成后合回 main。
- **标签（tag）**：标记正式版本的固定提交，例如 v1.2.0；后续修改不会自动移动它。
- **origin/main**：本地最近一次获取到的 GitHub main 位置；`fetch` 后才反映远端更新。

Mac、Linux 和 GitHub 的 main 同步到同一提交后，三处代码一致。每台机器保留现有一个源码目录，各系统自己的构建目录保存在本机。同步时核对 `git rev-parse HEAD`，源码一致后再在目标平台编译。

## 做一项修改的顺序

先保存当前工作区，再从最新 main 建临时分支。下面的分支名是示例，使用时改成此次工作的名称：

```sh
git status --short --branch
git switch main
git pull --ff-only origin main
git switch -c feature/example-change
```

完成改动与相关检查后，先审阅差异，再按本次实际文件名提交：

```sh
git diff --stat
git diff
git add doc/03_programs.md
git diff --cached
git commit -m "docs: explain program inputs and outputs"
```

上例 `git add` 只暂存一份文档，实际提交时逐项填写本次要保存的文件。`git status --short` 中未暂存的修改不会随这次提交保存。

验收完成后再合并并推送：

```sh
git switch main
git merge --ff-only feature/example-change
git push origin main
git branch -d feature/example-change
```

`--ff-only` 要求历史可直接前进；若因两端有独立提交而失败，先查看并整合双方改动。`git branch -d` 只删除已经合并的分支名称，已进入 main 的提交和代码仍在。已有 stash 独立保存，合并或删分支不会自动处理它。

## 看版本和历史

```sh
git log --graph --decorate --oneline --all -20
git branch -avv
git tag --list
git show --stat v1.2.0
```

`log --graph` 左边的线表示提交之间的关系，`main`、`origin/main`、`tag: v1.2.0` 是相应提交旁的名称。正式记录为 v1.0.0、v1.1.0、v1.2.0、v1.3.0，变更内容统一查 [CHANGELOG](../CHANGELOG.md)。

需要临时运行旧版时，先确保工作区改动已保存，切到标签并使用对应构建目录：

```sh
git switch --detach v1.1.0
cmake -S . -B build/v1.1.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.1.0 --parallel 3
git switch main
```

服务器构建仍按需加 Boost 路径。`--detach` 用于查看或构建固定版本，旧版命令和验收按该标签中的文档使用；准备在其基础上改代码时建立临时分支。

## 版本发布与两端同步

### 提交与发布是两件事

提交用于保存一次完成的修改，发布用于给一份验收后的源码确定正式版本。main 可以继续积累已验收的修复和文档，准备对外或供后续实验固定使用时，再统一发布。每一次提交都不必新建版本号、标签或 GitHub Release。

例如，本次 OpenLibm 单目标构建修复与文档整理分别保存一笔提交，可一起验收后集成到 main；正式 `v1.2.0` 仍指向 `d138c98`。两项变化记入根 [CHANGELOG](../CHANGELOG.md) 的“未发布”，本次继续使用 1.2.0 的程序版本号。

### 正式版本号怎样递增

ASTER 使用语义化版本 `MAJOR.MINOR.PATCH`，三段都是非负整数：主版本、次版本、修订版本。确定发布时，根据实际变化选择：

| 下次发布包含的变化 | 增加哪一段 | 从 1.2.0 举例 |
| --- | --- | --- |
| 保持现有接口和输入约定的修错 | PATCH | 1.2.1 |
| 保持兼容的新增功能 | MINOR，并将 PATCH 归零 | 1.3.0 |
| 不兼容的接口或输入约定变化 | MAJOR，并将 MINOR、PATCH 归零 | 2.0.0 |

版本号不是小数，也不逢十进一：`1.2.9` 后可为 `1.2.10`，`1.9.0` 后可为 `1.10.0`。文档和日常修复可先积累在“未发布”中，正式发布时再统一确定版本及变更范围。

程序三段版本号来自根 `CMakeLists.txt` 的 `project(... VERSION ...)`。正式发布先完成相关检查、更新版本号和 CHANGELOG 并形成验收提交，再为该提交建立 `vMAJOR.MINOR.PATCH` 标签。标签固定指向该发布提交；main 的后续提交继续前进。

### 实验记录要能找到实际源码

版本标签说明正式起点，提交编号定位实验实际使用的源码。在标签之后的 main 上运行实验，必须记录实际提交编号。先提交实验所用源码、重新编译，再记录：

```sh
git rev-parse HEAD
git describe --tags --always --dirty
build/v1.3.0/bin/aster_iod_gl --version
```

上例程序路径随实际构建目录填写。实验同时保存实际 setting、输入数据来源与校验值、数学后端、线程数，以及加噪和求解各自的随机种子。这样即使程序版本行仍显示 1.2.0，也能明确区分发布提交与后续修复。

### 同步已验收提交

Mac 以 GitHub 为 origin。服务器能正常连接时使用普通 Git 快进；现有离线中转方式见[服务器指南](10_server.md)。同步源码之后，在目标平台的构建目录重新编译。实验记录保留使用的提交和 `--version`，方便以后还原相同计算条件。
