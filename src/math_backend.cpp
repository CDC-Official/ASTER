/**
 * @file math_backend.cpp
 * @brief 将 ASTER 的标准双精度数学入口交给 OpenLibm 0.8.8。
 * @details 公共 core、SOFA 和应用限制 sin/cos 内建替换及
 * FMA 收缩，sqrt 和常数幂可用编译器优化；其余函数调用交给独立 OpenLibm。
 * 矩阵计算和基本四则运算仍使用 double。
 * sin/cos 以完整 binary64 参数缓存，每函数每线程最多 65536 项，冲突时重新求值。
 */
#include "math_backend.hpp"

#if ASTER_MATH_OPENLIBM
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dlfcn.h>
#include <filesystem>

#ifdef __APPLE__
#include <mach-o/dyld.h>
#include <vector>
#endif

#ifdef __APPLE__
#define ASTER_MATH_NOEXCEPT
#else
#define ASTER_MATH_NOEXCEPT noexcept
#endif

namespace
{
/**
 * @brief 取得与当前可执行文件放在同一目录的私有数学库路径。
 * @return 当前 ASTER 构建附带的 OpenLibm 路径，可随 bin 目录一起移动。
 */
std::filesystem::path OpenLibmPath()
{
    std::filesystem::path executable;
    std::error_code error;

#ifdef __APPLE__
    std::vector<char> path(1024);
    std::uint32_t size = static_cast<std::uint32_t>(path.size());

    if (_NSGetExecutablePath(path.data(), &size) != 0)
    {
        path.resize(size);

        if (_NSGetExecutablePath(path.data(), &size) != 0)
        {
            std::fprintf(stderr, "Cannot locate the ASTER executable.\n");
            std::abort();
        }
    }

    executable = path.data();
#else
    executable = std::filesystem::read_symlink("/proc/self/exe", error);
#endif

    executable = std::filesystem::canonical(executable, error);

    if (error)
    {
        std::fprintf(stderr, "Cannot locate the ASTER executable: %s\n", error.message().c_str());
        std::abort();
    }

    return executable.parent_path() / ASTER_OPENLIBM_NAME;
}

/**
 * @brief 从本次构建的 OpenLibm 取得指定数学入口。
 * @param[in] name OpenLibm 导出的 C 函数名。
 * @return 对应函数地址；库或函数不可用时报告具体原因并终止计算。
 */
void* OpenLibmFunction(const char* name)
{
    static void* library = dlopen(OpenLibmPath().c_str(), RTLD_NOW | RTLD_LOCAL);

    if (!library)
    {
        std::fprintf(stderr, "Cannot load OpenLibm: %s\n", dlerror());
        std::abort();
    }

    void* function = dlsym(library, name);

    if (!function)
    {
        std::fprintf(stderr, "OpenLibm function is unavailable: %s\n", name);
        std::abort();
    }

    return function;
}

struct TrigCacheEntry
{
    std::uint64_t argument = 0;
    double result = 0;
    bool valid = false;
};

/**
 * @brief 计算一个正弦或余弦并缓存同一 binary64 输入的结果。
 * @param[in] angle 弧度制参数，按完整64位表示区分输入。
 * @return OpenLibm 的 double 结果，缓存不截断输入或输出。
 */
template<bool Sine> double CachedTrig(double angle)
{
    using Function = double(*)(double);
    static const Function function = reinterpret_cast<Function>(OpenLibmFunction(Sine ? "sin" : "cos"));
    thread_local std::array<TrigCacheEntry, 65536> cache;
    std::uint64_t bits;
    std::memcpy(&bits, &angle, sizeof(bits));
    const std::uint64_t hash = (bits ^ (bits >> 32)) * 0x9e3779b97f4a7c15ULL;
    TrigCacheEntry& entry = cache[(hash ^ (hash >> 32)) & (cache.size() - 1)];

    if (entry.valid && entry.argument == bits) return entry.result;

    const double result = function(angle);
    entry = {bits, result, true};
    return result;
}
}

// 这些入口只绑定本目标的数学调用。参数、定义域和返回值由 OpenLibm 决定。
#define ASTER_UNARY(name) \
extern "C" double name(double x) ASTER_MATH_NOEXCEPT \
{ \
    using Function = double(*)(double); \
    static const Function function = reinterpret_cast<Function>(OpenLibmFunction(#name)); \
    return function(x); \
}

#define ASTER_BINARY(name) \
extern "C" double name(double x, double y) ASTER_MATH_NOEXCEPT \
{ \
    using Function = double(*)(double, double); \
    static const Function function = reinterpret_cast<Function>(OpenLibmFunction(#name)); \
    return function(x, y); \
}

ASTER_UNARY(tan) ASTER_UNARY(asin) ASTER_UNARY(acos) ASTER_UNARY(atan)
ASTER_UNARY(sqrt) ASTER_UNARY(cbrt)
ASTER_UNARY(exp) ASTER_UNARY(exp2) ASTER_UNARY(expm1)
ASTER_UNARY(log) ASTER_UNARY(log2) ASTER_UNARY(log10) ASTER_UNARY(log1p)
ASTER_UNARY(sinh) ASTER_UNARY(cosh) ASTER_UNARY(tanh)
ASTER_UNARY(asinh) ASTER_UNARY(acosh) ASTER_UNARY(atanh)
ASTER_UNARY(floor) ASTER_UNARY(ceil) ASTER_UNARY(round)
ASTER_UNARY(trunc) ASTER_UNARY(rint) ASTER_UNARY(nearbyint)
ASTER_UNARY(erf) ASTER_UNARY(erfc) ASTER_UNARY(lgamma) ASTER_UNARY(tgamma)
ASTER_BINARY(pow) ASTER_BINARY(atan2) ASTER_BINARY(hypot)
ASTER_BINARY(fmod) ASTER_BINARY(remainder) ASTER_BINARY(copysign)

extern "C" double sin(double angle) ASTER_MATH_NOEXCEPT { return CachedTrig<true>(angle); }
extern "C" double cos(double angle) ASTER_MATH_NOEXCEPT { return CachedTrig<false>(angle); }

/**
 * @brief 将数分解为二进制尾数和整数指数。
 * @param[in] value 待分解的 double 数值。
 * @param[out] exponent 接收二进制指数。
 * @return OpenLibm 计算的尾数。
 */
extern "C" double frexp(double value, int* exponent) ASTER_MATH_NOEXCEPT
{
    using Function = double(*)(double, int*);
    static const Function function = reinterpret_cast<Function>(OpenLibmFunction("frexp"));
    return function(value, exponent);
}

/**
 * @brief 按二进制整数幂缩放一个数。
 * @param[in] value 待缩放的 double 数值。
 * @param[in] exponent 以2为底的整数指数。
 * @return OpenLibm 计算的 value*2^exponent。
 */
extern "C" double scalbn(double value, int exponent) ASTER_MATH_NOEXCEPT
{
    using Function = double(*)(double, int);
    static const Function function = reinterpret_cast<Function>(OpenLibmFunction("scalbn"));
    return function(value, exponent);
}

/**
 * @brief 将数分解为整数部分和小数部分。
 * @param[in] value 待分解的 double 数值。
 * @param[out] integer 接收与输入同号的整数部分。
 * @return OpenLibm 计算的小数部分。
 */
extern "C" double modf(double value, double* integer) ASTER_MATH_NOEXCEPT
{
    using Function = double(*)(double, double*);
    static const Function function = reinterpret_cast<Function>(OpenLibmFunction("modf"));
    return function(value, integer);
}

// binary64 的基数为2，ldexp与scalbn都执行精确的整数二进制幂缩放。
extern "C" double ldexp(double value, int exponent) ASTER_MATH_NOEXCEPT { return scalbn(value, exponent); }

#endif // ASTER_MATH_OPENLIBM

/**
 * @brief 返回本次构建的数学后端和浮点配置。
 * @return 供版本输出和诊断使用的静态字符串。
 */
const char* aster::math_backend::BuildDescription() noexcept
{
    return ASTER_MATH_BUILD_DESCRIPTION "; build=" ASTER_MATH_BUILD_TYPE;
}
