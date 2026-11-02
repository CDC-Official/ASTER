# 实验室服务器

## 1.1.0 当前使用

日常源码统一使用原 `/home/cdc/projects/ASTER` 的 main，Mac 使用原 `/Users/chendingcun/CLionProjects/ASTER`；正式程序位于各自 `build/v1.1.0/bin/`。版本以标签回切，临时 dev worktree 在发布后收束，必要验证材料保存在原仓库的构建目录。当前数学后端、构建与验证范围见 [构建说明](build.md)、[后端验收](math-backend.md)和 [版本更新记录](../CHANGELOG.md)。

## 1.0.0 部署与兼容性历史

ASTER 已部署到 `cdc-114` 的 `/home/cdc/projects/ASTER`。主机为 lenovo-ThinkStation-PX，Ubuntu x86_64，用户为 cdc。GCC 构建目录为 `build/linux-gcc-release`，Clang 构建目录为 `build/linux-clang-release`，各自的九个程序位于对应目录的 `bin/`。

## 版本与依赖

原始发布标签 `v1.0.0` 对应 `4460983`。本次部署在该基线上增加了三个独立提交：`d2cb63d` 将输出数组参数改为标准指针声明；`b4713fc` 将单字符分隔符的字符串复制改为字符赋值；`4b58175` 让空白输入检查使用同一平台的普通行作为对照。科研算法和原标签均未改变。

服务器使用 GCC 13.3.0、CMake 3.28.3、Python 3.12.3，Eigen 3.4.0 和 SOFA 2018-01-30 随库构建。Boost 1.87.0 头文件和原许可从 Mac 离线传入，位于 `/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0`；15,899 个头文件逐项校验通过。GCC 的 OpenMP 使用服务器已有 libgomp。Clang 18.1.3 使用匹配的 Ubuntu 官方 OpenMP 18 包，解包在 `/home/cdc/projects/aster-clang-20261102/dependencies/openmp18`；包来源、SHA256 及用户目录构建方式见 [Clang 与 GLML 对照记录](glml-cross-platform.md)。

部署材料、构建日志和原始结果保存在 `/home/cdc/projects/aster-deployment-20261101`。JPL、EOP、重力场和大气数据已随仓库到位，四份文件的校验值与 Mac 一致。

## 使用

```sh
ssh cdc-114
cd /home/cdc/projects/ASTER
build/linux-gcc-release/bin/aster_iod_gl --version
python3 ../aster-deployment-20261101/run_server_case.py --case opn_cartesian
```

算例脚本为每次运行创建独立目录，并在输出中给出路径。查看该目录的 `setting.txt` 和 `console.txt` 后，可以复制配置用于自己的资料。

重新构建：

```sh
cmake -S . -B build/linux-gcc-release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DASTER_BOOST_ROOT=/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0
cmake --build build/linux-gcc-release --parallel 2
ctest --test-dir build/linux-gcc-release --output-on-failure -R '^(core_|data_integrity|invalid_inputs|iod_batch_handoff)'
```

## 验收结果

九个程序均已编译并运行成功。七项核心数值检查、数据完整性、初轨到批处理的文件接口以及输入契约检查通过。另用冒号配置验证了单字符分隔符修复。Mac 上重新编译后，12 项受影响入口和接口检查也全部通过。

原测试中的应用参考文件来自 Mac，采用逐字节比较。首次 Linux 全套 CTest 为 10 项通过、9 项失败；其中输入契约的比较方式已修正并重跑通过，其余八项字节比较仍报告差异。跨平台数值对照保存在 `linux-validation.json`，同时保留了原始输出。

本轮状态对照门限为位置 1 mm、速度 1e-6 m/s，角残差为 1e-10 rad。八个程序的代表结果在相应门限内：GL 位置差 1.65e-5 m、速度差 7.41e-8 m/s；BODN 角残差最大差 2.89e-15 rad；两套地球预报位置最大差小于 6e-9 m；地月星历逐字节一致。状态转换和瞬转平的差异也在门限内。

2026-11-02 的同机 Linux GCC / Clang 对照已通过：九个代表算例的 16 份科学输出逐字节一致，Clang 的十项核心、数据和接口检查全部通过。

Mac/Linux 的 10 条观测短弧 GLML 仍未通过原有数值一致性门限：Mac 选候选 592，Linux 选候选 335，位置差 15,387.34 m、速度差 31.0198 m/s，J 相差 0.0004673。两端来源均达到 40 轮上限。候选迭代对微小浮点差异的放大仍存在；本算例未发生随机回退，单、双线程结果一致。证据与复现命令见 [专项记录](glml-cross-platform.md)，逐项结果见 [clang-glml-validation.json](clang-glml-validation.json)。

作者已确认保留现有选择规则：达到迭代上限的轨道在满足原有其他有效性要求后，仍按最终 J 参与选优，来源收敛标志仅作诊断。决定及停止判据、目标函数和有限延长试验见 [GLML 收敛诊断](glml-convergence.md)。

服务器尚未启动正式批量实验，运行资料与输出放在各自实验目录。

## GitHub 与日常同步

私有仓库为 [CDC-Official/ASTER](https://github.com/CDC-Official/ASTER)。Mac 与服务器均以该仓库作为 `origin`，日常使用 `main`；原 `v1.0.0` 标签保留。只推送 `main` 和需要发布的明确标签，临时部署分支留在本地。

2026-11-02，Mac 的既有认证可访问该私有仓库，服务器无法解析 GitHub 域名，正常 HTTPS 访问失败。本次服务器通过 Mac 中转的增量 Git bundle 同步；这不表示服务器已经能直接从 GitHub pull。

Mac 开始修改前执行 `git pull --ff-only`，完成修改后按明确文件名 `git add`、`git commit`，再执行 `git push origin main`。同步前先检查 `git status --short`，工作区中的修改应先提交或单独保存。

服务器当前使用以下中转方式。在 Mac 的 ASTER 根目录执行；两端 main 已相同时无需传送：

```sh
aster_server_head=$(ssh -o ClearAllForwardings=yes cdc-114 'git -C /home/cdc/projects/ASTER rev-parse main')
git bundle create /tmp/ASTER-main.bundle main "^${aster_server_head}"
scp -o ClearAllForwardings=yes /tmp/ASTER-main.bundle cdc-114:/home/cdc/projects/aster-sync/ASTER-main.bundle
```

随后在服务器执行：

```sh
cd /home/cdc/projects/ASTER
git bundle verify /home/cdc/projects/aster-sync/ASTER-main.bundle
git fetch --no-tags /home/cdc/projects/aster-sync/ASTER-main.bundle main:refs/remotes/origin/main
git merge --ff-only origin/main
```

bundle 中的 main 来自 Mac 已推送的提交，服务器的 `origin` 始终指向 GitHub。若服务器产生了独立的新提交，应先在服务器创建 main 的 bundle 传回 Mac，再在 Mac 快进整合并推送；不要覆盖任一端的新提交。服务器以后恢复正常 GitHub 网络和既有认证后，可直接使用 `git pull --ff-only`。
