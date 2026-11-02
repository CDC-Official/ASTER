# 构建与迁移

## 当前版本构建

1.1.0 使用 C++17、CMake 3.26 及以上版本。已验证 Mac arm64 / AppleClang 21、Linux x86_64 / GCC 13.3，Boost 1.87.0 头文件及现有 OpenMP。Eigen 固定 3.4.0，SOFA 固定 2018-01-30，OpenLibm 固定官方 0.8.8，源码与许可均随库保存。

```sh
cmake -S . -B build/v1.1.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.1.0 --parallel 3
build/v1.1.0/bin/aster_iod_gl --version
ctest --test-dir build/v1.1.0 -R '^(core_frame_coordinates|math_backend_sofa|invalid_inputs|data_integrity)$' --output-on-failure
```

完整有限验收见 [测试说明](../tests/README.md)。默认 `ASTER_MATH_BACKEND=OPENLIBM`；使用系统数学时另建 `build/v1.1.0-system`，配置时添加 `-DASTER_MATH_BACKEND=SYSTEM`。两种配置共用相同源码。OPENLIBM 程序加载同目录私有共享库，移动程序时一起移动该库；不会替换系统 libm。

Boost 不在默认路径时添加 `-DASTER_BOOST_ROOT=/path/to/boost`。Mac 先使用 CMake 的 OpenMP 查找，再查询现有 Homebrew libomp；Linux 使用现有 libgomp。服务器继续复用 `/home/cdc/projects/aster-deployment-20261101/dependencies/boost_1_87_0`。仅构建应用可加 `-DBUILD_TESTING=OFF`；启用测试需要 Python 3。

## 同一仓库管理版本

Mac 日常源码目录为 `/Users/chendingcun/CLionProjects/ASTER`，Linux 为 `/home/cdc/projects/ASTER`。版本由分支和标签保留，不为每个版本建立长期源码副本。切换前保存工作区改动，再使用对应构建目录，例如：

```sh
git switch --detach v1.0.0
cmake -S . -B build/v1.0.0 -DCMAKE_BUILD_TYPE=Release
cmake --build build/v1.0.0 --parallel 3
git switch main
```

`CMakeLists.txt` 是程序三段版本号的唯一来源。`v1.0.0` 保持原提交，`v1.1.0` 指向本次正式发布。历史试验使用本地备份/stash 保留，不混入正式 main。

## 离线同步服务器

继续使用现有 SSH 和 Git bundle，不需要更改网络或新增凭据：

```sh
git bundle create /tmp/ASTER-main.bundle main v1.1.0
# 传输后，在服务器现有 ASTER 仓库内核对 bundle，再正常快进 main。
git bundle verify /path/to/ASTER-main.bundle
git fetch /path/to/ASTER-main.bundle main:refs/remotes/origin/main
git merge --ff-only origin/main
```

数据校验见 `data/manifest.json`。1.1.0 的统一后端已在代表例上实现 Mac/Linux 科学输出一致；[旧平台差异记录](glml-cross-platform.md)保留 SYSTEM 路径的历史问题，当前范围见 [后端验收](math-backend.md)。
