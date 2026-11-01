# 构建与迁移

## 本地构建

本版使用 C++17、CMake 3.26 及以上版本。Eigen 固定为随库的 3.4.0，SOFA 固定为 2018-01-30。Boost 使用头文件；GLML 的并行候选生成使用 OpenMP。构建过程读取本地依赖。

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

Mac 配置会先尝试 CMake 的 OpenMP 查找，再查询已安装 Homebrew 的 libomp 路径。当前验证使用 Apple Clang、Boost 1.87.0 和本机 libomp；更换编译器后请使用新的构建目录。仅需应用时可设置 `-DBUILD_TESTING=OFF`；启用测试需要 Python 3。

编译结果位于构建目录的 `bin/`。可显式构建单个目标，例如 `cmake --build build/release --target aster_iod_glml`。

## 实验室 Linux

将源码与 `data/` 一起传到服务器，在服务器上重新配置、编译。提前准备 CMake、C/C++ 编译器、Boost 头文件和对应的 OpenMP 运行库；Python 3 用于检查算例。校园网内可通过 SSH 访问 Git 仓库，也可以在 Mac 上生成 bundle 后传输：

```sh
git bundle create ASTER-1.0.0.bundle --all
# 将 bundle 文件传到服务器后执行：
git clone ASTER-1.0.0.bundle ASTER
cd ASTER
git checkout v1.0.0
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build/release --parallel 3
ctest --test-dir build/release --output-on-failure
```

数据校验清单在 `data/manifest.json`。正式批量实验前，应完成 Linux 短算例检查，记录编译器、提交、配置、数据校验值和线程数。首次服务器构建是下一阶段的工作。

## 版本管理

`CMakeLists.txt` 的项目版本是程序版本的唯一来源。`v1.0.0` 标签对应通过验收的提交，后续修正使用新的提交与版本号。大型实验从固定标签建立独立运行目录，结果和日志保存到实验目录，日常开发继续在主仓库进行。
