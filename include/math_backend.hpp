/**
 * @file math_backend.hpp
 * @brief 提供 ASTER 公共数学后端的构建信息。
 * @details 通过 ASTER::math_backend 链接标准 double 数学入口及浮点编译规则。
 * OPENLIBM 使用邻接的固定版本库；SYSTEM 使用工具链原生数学实现。
 */
#pragma once
#include "math_backend_config.hpp"

namespace aster::math_backend
{
/**
 * @brief 返回本次构建的数学后端和浮点配置。
 * @return 供版本输出和诊断使用的静态字符串。
 */
const char* BuildDescription() noexcept;
}
