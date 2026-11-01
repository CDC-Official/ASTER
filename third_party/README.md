# 第三方源码

Eigen 3.4.0 与 SOFA 2018-01-30 沿用原工程的固定快照，算法源码保持原样，原始许可和来源记录随目录保存。Eigen 的主要许可见 `eigen/COPYING.README`、`eigen/COPYING.MPL2`，其他组成部分的许可保留在对应文件中。SOFA 的使用说明及许可见 `sofa/00READ.ME` 和源码文件头。

CMake 接入逻辑放在 `cmake/`。Boost 与 OpenMP 使用平台已有依赖；Mac 的动态库和构建缓存不作为源码依赖迁入。依赖升级在后续版本中单独验证。
