# 快速上手与构建

本页从现有 ASTER 源码目录开始，完成编译、版本查询和一次根数转换。日常实验继续使用这一份源码，程序和结果分别放在 `build/` 的子目录中。

## 1. 确认自己在仓库根目录

Mac 的现有目录是 `/Users/chendingcun/CLionProjects/ASTER`，服务器是 `/home/cdc/projects/ASTER`。进入相应目录后执行：

```sh
pwd
git status --short --branch
```

根目录应能看到 `CMakeLists.txt`、`apps/`、`include/`、`src/` 和 `tests/`。`git status` 会显示当前分支及尚未提交的修改。

## 2. 编译当前版本

所需环境为 CMake 3.26 及以上、支持 C++17 的 C/C++ 编译器、Boost 头文件、OpenMP，以及运行检查和加噪工具所需的 Python 3。仓库自带 Eigen 3.4.0、SOFA 2018-01-30 和 OpenLibm 0.8.8，来源和许可见[第三方说明](../third_party/README.md)。

1.3.0 已使用 Mac arm64 / AppleClang 21 和 Linux x86_64 / GCC 13.3 完成全目标构建与代表算例检查。两端共 37 份科学输出逐字节一致，具体范围见[发布验收](validation/11_release-1.3.0.md)。

### Mac

在现有 Mac 环境中执行：

```sh
cmake -S . -B build/v1.3.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.3.0 --parallel 3
```

本机 Boost 头文件位于 `/usr/local/include/boost_1_87_0`，CMake 已包含这个搜索位置。OpenMP 先通过 CMake 查找；Apple 平台查找失败时，再查询现有 Homebrew libomp。

### Linux 服务器

服务器复用已保存的 Boost，配置时明确给出路径：

```sh
cmake -S . -B build/v1.3.0 -DCMAKE_BUILD_TYPE=Release \
    -DASTER_BOOST_ROOT=/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0
cmake --build build/v1.3.0 --parallel 3
```

`-S .` 表示读取当前目录的源码，`-B` 指定构建目录，`Release` 启用发布构建，`--parallel 3` 同时使用三个编译任务。配置成功后再构建，生成的程序位于 `build/v1.3.0/bin/`。

已有 `release` preset 也可使用：`cmake --preset release` 和 `cmake --build --preset release`。它将程序放在 `build/release/bin/`；后续命令应使用与所选构建目录对应的路径。

## 3. 查询程序版本

```sh
build/v1.3.0/bin/aster_iod_gl --version
```

首行应包含程序名和 `1.3.0`。多数程序的后续文字显示实际数学后端、构建类型和浮点条件。默认后端为 `OPENLIBM-0.8.8`。

OPENLIBM 程序需要同目录中的 `libopenlibm.dylib`（Mac）或 `libopenlibm.so`（Linux）。移动程序时一起移动该库。切换后端和链接公共库的方法见[数学后端](08_math-backend.md)。

## 4. 跑出第一份结果

先用不需要配置文件的根数转换熟悉调用。以下命令读取仓库内的一行 Kepler 根数，写出 Cartesian 状态，再转回根数：

```sh
mkdir -p build/quick-start
build/v1.3.0/bin/aster_convert_keplerian_to_cartesian \
    tests/convert_keplerian_to_cartesian/input.txt build/quick-start/cartesian.txt
build/v1.3.0/bin/aster_convert_cartesian_to_keplerian \
    build/quick-start/cartesian.txt build/quick-start/keplerian.txt m
```

成功后，`cartesian.txt` 是七列：UTC MJD、三个位置分量（m）和三个速度分量（m/s）；`keplerian.txt` 也是七列：UTC MJD、a（km）、e、i、OMEGA、omega、M（deg）。M 是平近点角。两次转换保持原历元和参考系，数值格式见[数据约定](05_data-formats.md)。这两项转换会覆盖同名输出。

下一步可按[程序使用](03_programs.md)准备 GL、预报、BODN 或瞬转平算例；要生成仿真观测，按[观测生成](06_observations.md)运行 LEO 算例。

## 5. 在 CLion 中运行

用 CLion 打开包含根 `CMakeLists.txt` 的 ASTER 目录，选择现有本机 C/C++ 工具链。沿用已有 Release 配置的 `cmake-build-release`，Debug 使用 `cmake-build-debug`；程序分别在这些目录的 `bin/` 中。Boost 不在默认搜索位置时，将上述 `-DASTER_BOOST_ROOT=...` 填入 CMake 选项。

CLion 的 Ninja 配置与终端的 Unix Makefiles 配置分别使用自己的构建目录。修改 CMake 后重新载入项目；当前公共数学接口会传播 OpenLibm 的构建依赖，单独构建 GLML 等运行目标时也会生成邻接共享库。

重新载入 CMake 后，选择要运行的 `aster_...` 目标。把工作目录设为仓库根目录，命令行参数按程序分别填写，例如：

```text
-I build/examples/trial-001/iod_gl/setting.txt
```

这份配置须先按[程序使用](03_programs.md)中的准备命令生成。两项直接转换填写输入、输出路径，Cartesian 转根数还要填写 `m` 或 `km`。工作目录决定命令行相对路径的位置；配置内部的相对文件路径从配置所在目录解释。

## 6. 做一次小范围检查

```sh
ctest --test-dir build/v1.3.0 \
    -R '^(core_frame_coordinates|math_backend_sofa|data_integrity)$' --output-on-failure
```

这三项分别检查坐标转换、SOFA/数学后端和科学数据校验值。其他检查入口、原始字节参考与当前后端验收的区别见[测试说明](../tests/README.md)。

## 常见错误

| 看到的现象 | 检查和处理 |
| --- | --- |
| CMake 找不到 Boost | `ASTER_BOOST_ROOT` 应指向含 `boost/` 的目录或其安装根目录；核对实际头文件位置后重新配置 |
| OpenMP 或动态调度链接失败 | 核对编译器与 OpenMP 运行库是否配套。Mac 用现有 libomp，服务器 GCC 用系统 libgomp；GLML 需要动态调度链接检查通过 |
| 启动时报找不到 OpenLibm | 检查可执行文件旁的私有共享库是否存在，并确认使用了预期构建目录中的程序 |
| 找不到 setting 或资料文件 | 先看程序打印的 setting 绝对路径，再按该配置所在目录展开内部相对路径 |
| 无法写出结果 | 先创建父目录，再核对同名文件的处理规则；瞬转平和加噪的规则见[程序指南](03_programs.md) |
| EOP 无覆盖或历元越界 | 核对时段与数据覆盖范围；需要插值的时刻还要有下一日 EOP，详见[科学数据](../data/README.md) |
