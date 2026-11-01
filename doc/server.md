# 实验室服务器

ASTER 已部署到 `cdc-114` 的 `/home/cdc/projects/ASTER`。主机为 lenovo-ThinkStation-PX，Ubuntu x86_64，用户为 cdc。构建目录是 `build/linux-gcc-release`，九个程序位于其中的 `bin/`。

## 版本与依赖

原始发布标签 `v1.0.0` 对应 `4460983`。本次部署在该基线上增加了三个独立提交：`d2cb63d` 将输出数组参数改为标准指针声明；`b4713fc` 将单字符分隔符的字符串复制改为字符赋值；`4b58175` 让空白输入检查使用同一平台的普通行作为对照。科研算法和原标签均未改变。

服务器使用 GCC 13.3.0、CMake 3.28.3、Python 3.12.3，Eigen 3.4.0 和 SOFA 2018-01-30 随库构建。Boost 1.87.0 头文件和原许可从 Mac 离线传入，位于 `/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0`；15,899 个头文件逐项校验通过。OpenMP 使用服务器已有 libgomp。Clang 18.1.3 已安装，但缺少配套 OpenMP 头文件及运行库，本次未安装额外工具链。

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

10 条观测的短弧 GLML 未通过数值一致性检查：Mac 选候选 592，Linux 选候选 335，位置差 15,387.34 m、速度差 31.0198 m/s，J 相差 0.0004673。两端所选候选的来源状态均为 -1，即达到迭代上限。当前规则允许这类有限 J 的候选参与选择，本次保持该规则。进行短弧 GLML 批量实验前，需要单独处理这项数值敏感性。

服务器尚未启动正式批量实验。后续同步继续传送 Git bundle 或通过已授权的 SSH 路径更新固定提交，运行资料与输出放在各自实验目录。
