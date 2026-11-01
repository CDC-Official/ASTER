/**
 * @file orbit_expansion.hpp
 * @brief Kaula 倾角函数、导数及 Hansen 偏心率展开。
 */
#pragma once

double Flmp(int l, int m, int p, double i);

// 未归一化 Kaula 倾角函数及对 i(rad) 的一、二阶导数；半角多项式在 i=0 也可求值。
void FlmpDerivatives(int l, int m, int p, double i, double& F, double& dF, double& ddF);

// Hansen X^{l,p}_q 的 0..4 阶偏心率展开。dX 保留参考解析定轨算法的截断阶次：
// q-p=0、±2 仅保留 O(e)，±1 保留 O(e²)，±3、±4 置零，并非 X 展开的完整导数。
void HansenSeries(int l, int p, int q, double e, int order, double& X, double& dX);
