//
// Created by 陈酊存 on 2024/10/16.
//

#ifndef CDC_LIB_INCLUDE_ORBIT_EQUATION_HPP_
#define CDC_LIB_INCLUDE_ORBIT_EQUATION_HPP_

#include "orbit_statedef.hpp"
#include "utils_overwatch.hpp"
using namespace orbit;
using namespace Identifier;

void GaussPlanetaryEquation_RTNtoNSG2(bastion::Overwatch& tracer,
	const NonSingular2& nsg2, const double* RTN, double* f);
#endif //CDC_LIB_INCLUDE_ORBIT_EQUATION_HPP_
