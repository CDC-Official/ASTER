# 服务器使用与同步

本页记录现有 `cdc-114` 工作目录中的日常操作。两端共用跨平台源码，分别编译；通用构建和程序用法见[快速上手](02_quick-start.md)、[程序使用](03_programs.md)。

## 登录与构建

```sh
ssh -o ClearAllForwardings=yes cdc-114
cd /home/cdc/projects/ASTER
git status --short --branch
```

服务器源码目录固定为 `/home/cdc/projects/ASTER`，日常使用 main。1.2.0 发布时，两端 main 和 v1.2.0 均指向 `d138c98`。程序在各自的 `build/v1.3.0/bin/`，旧版源码通过标签查看，旧版构建可继续保留在对应 `build/` 子目录。

在服务器执行：

```sh
cmake -S . -B build/v1.3.0 -DCMAKE_BUILD_TYPE=Release \
    -DASTER_BOOST_ROOT=/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0
cmake --build build/v1.3.0 --parallel 3
build/v1.3.0/bin/aster_iod_gl --version
```

发布构建使用 Linux x86_64 / GCC 13.3、既有 libgomp 和 Boost 1.87.0。Eigen、SOFA、OpenLibm 随仓库源码构建。早期 Linux 数值对照见[验证索引](validation/01_index.md#linux)，编译器与短弧候选对照见[GLML 诊断](validation/02_glml-sensitivity.md)，统一后端的后续结果见[1.1.0 验收](validation/03_math-backend.md)。

## 在哪里运行实验

小算例按[程序指南](03_programs.md)复制配置到 `build/examples/`；每次实际实验使用独立目录，保留 setting、输入来源、程序版本、提交编号和随机种子。结果路径按配置决定，构建或切换分支不会自动搬移实验资料。

GLML 的 OpenMP 线程数可在启动时填写，例如：

```sh
OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE build/v1.3.0/bin/aster_iod_glml \
    -I build/examples/trial-001/iod_glml/setting.txt
```

先准备配置再执行。这个示例使用两个线程；实验批次的并发进程和每进程线程数应在实验记录里分别注明。

## GitHub 可用时同步

仓库为 [CDC-Official/ASTER](https://github.com/CDC-Official/ASTER)，两端 origin 指向同一地址。开始工作前检查当前分支和工作区，再快进 main：

```sh
git status --short --branch
git switch main
git pull --ff-only origin main
```

临时功能分支、提交和合并步骤见[Git 维护](09_development.md)。两端 main 相同时已完成源码同步，程序在本机重新构建。

## 服务器无法直连 GitHub 时中转

使用现有 SSH 传输 Git bundle。先在 Mac 将已验收 main 推送到 GitHub，确认 Mac main 包含服务器 main，再生成增量包。以下在 Mac 的 ASTER 根目录执行：

```sh
aster_server_head=$(ssh -o ClearAllForwardings=yes cdc-114 'git -C /home/cdc/projects/ASTER rev-parse main')
git log --oneline "${aster_server_head}..main"
git merge-base --is-ancestor "$aster_server_head" main
```

若没有新提交，两端已相同。若祖先检查成功且确有新提交，继续：

```sh
git bundle create /tmp/ASTER-main.bundle main "^${aster_server_head}"
scp -o ClearAllForwardings=yes /tmp/ASTER-main.bundle cdc-114:/home/cdc/projects/aster-sync/ASTER-main.bundle
```

随后在服务器的干净 main 上执行：

```sh
cd /home/cdc/projects/ASTER
git status --short --branch
git switch main
git bundle verify /home/cdc/projects/aster-sync/ASTER-main.bundle
git fetch --no-tags /home/cdc/projects/aster-sync/ASTER-main.bundle main:refs/remotes/origin/main
git merge --ff-only origin/main
git log -1 --oneline
```

bundle 传送提交历史，origin 的地址继续保持 GitHub。服务器有独立新提交时，先从服务器生成该分支的 bundle 传回 Mac，在 Mac 获取并整合，再推送和同步；保留双方的新提交。祖先检查或快进失败时，从 `git log --graph --decorate --oneline --all` 查看分叉位置。

这里的增量包只传 main 新提交。发布标签随发布流程单独同步，日常文档或功能同步不移动已有版本标签。
