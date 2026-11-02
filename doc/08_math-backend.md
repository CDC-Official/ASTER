# 数学后端

ASTER 默认使用 `OPENLIBM`，把十个主程序、公共核心和 SOFA 中的 double 数学路径纳入同一后端与浮点构建设置。`SYSTEM` 使用本机工具链的原生数学实现。运行前用程序的 `--version` 确认实际配置，并把它和实验输入一起记录。

## 选择构建

正常使用按[构建指南](02_quick-start.md)生成默认 OPENLIBM。需要和系统数学比较时，为 SYSTEM 使用另一个构建目录：

```sh
cmake -S . -B build/system -DCMAKE_BUILD_TYPE=Release -DASTER_MATH_BACKEND=SYSTEM
cmake --build build/system --parallel 3
build/system/bin/aster_iod_gl --version
```

`ASTER_MATH_BACKEND` 只接受 `OPENLIBM` 或 `SYSTEM`。两种配置共用源码，各自保存 CMake 缓存、程序和检查结果。服务器配置时仍需填写既有 Boost 路径。

OPENLIBM 使用固定的官方 OpenLibm 0.8.8，来源提交、归档 SHA-256 和许可见 [openlibm-source.json](../third_party/openlibm-source.json)及 [LICENSE.md](../third_party/openlibm/LICENSE.md)。程序从自身所在目录加载 `libopenlibm.dylib` 或 `libopenlibm.so`，移动时一起移动；SYSTEM 程序使用系统数学库。

## C++ 目标怎样接入

[math_backend.hpp](../include/math_backend.hpp)、[math_backend.cpp](../src/math_backend.cpp)和 [MathBackend.cmake](../cmake/MathBackend.cmake)集中管理公共入口、构建信息和浮点条件。科学程序链接 `ASTER::core`；只需要数学入口的目标可链接 `ASTER::math_backend`。对象和规则通过 CMake 目标传递，示例见[公共接口参考](07_interfaces.md)。

OPENLIBM 的共享库依赖通过公共接口传播。Ninja/CLion 单独构建某个应用目标时，也会生成所需的 OpenLibm；实际加载仍从可执行文件旁查找。

OPENLIBM 的构建条件为：

- 公共核心、主程序、SOFA 和数学入口使用 `-fno-builtin-sin -fno-builtin-cos -fno-fast-math -ffp-contract=off`，C++ 定义 `EIGEN_DONT_VECTORIZE`。
- double 入口覆盖当前源码使用的三角、反三角、平方根、幂、指数、对数、取整、余数和相关函数；float/long double 保留自身路径，当前正式科学路径使用 double。
- sin/cos 使用完整 binary64 键缓存，每函数每线程 65,536 项，碰撞时重新求值；矩阵四则运算仍为 double。
- 上游 OpenLibm 使用自己的构建选项，源码保持原样。Linux 的私有库增加 `-Wl,-Bsymbolic-functions`，将库内调用固定在 OpenLibm 内。

SYSTEM 使用工具链默认浮点/Eigen 条件。两种后端都把日志 `__FILE__` 中的源码根路径映射为相对路径。Eigen 固定 3.4.0，SOFA 固定 2018-01-30。

## 为什么实验记录需要后端信息

部分短弧 GL/GLML 候选迭代会放大很小的浮点差异，改变达到迭代上限时的轨道和候选评分。1.1.0 统一数学库与编译设置后，已检查的 Mac/Linux 代表算例得到逐字节相同的科学输出；1.2.0 的观测生成、加噪与定轨接续样本也完成了两端比较。

后端变化的影响应按实际结果比较。已有 e070 短弧的 SYSTEM→OPENLIBM 对照出现明显轨道变化，七个普通程序的代表结果仍在既有物理量门限内。具体数值、计时和四场景真值误差见[1.1.0 后端验收](validation/03_math-backend.md)，观测功能的文件比较见[1.2.0 验收](validation/04_observations.md)。

GL/GLML 的科学公式、状态容差和有限 J 选优规则由当前求解器和主程序定义。Ljung–Box 旧调用已停用，遗留 Q/p 字段为 −1。

## 检查入口

日常可先运行坐标/SOFA 和四个已有单 GL：

```sh
ctest --test-dir build/v1.2.0 -R '^(core_frame_coordinates|math_backend_sofa)$' --output-on-failure
python3 tests/math_backend_cases/run_checks.py --bin-dir build/v1.2.0/bin \
    --output build/v1.2.0/acceptance/gl --only-gl
```

四个单 GL 和两个完整 GLML 使用保存的 OPENLIBM 科学输出哈希。七个普通程序由 `tests/check_math_backend_apps.py` 按既有物理量门限比较；完整命令与旧 fixture 的分工见[测试说明](../tests/README.md)。1.0.0 原生数学的跨平台门限保存在[原始 JSON](validation/evidence/03_linux-validation.json)，1.1.0 的后端结果保存在[验收 JSON](validation/evidence/06_math-backend-validation.json)。
