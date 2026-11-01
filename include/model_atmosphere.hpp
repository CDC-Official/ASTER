//
// Created by 陈酊存 on 2024/1/23.
//

#ifndef CDC_MODEL_ATMOSPHERE_HPP
#define CDC_MODEL_ATMOSPHERE_HPP
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <Eigen/Dense>
#include "utils_constants.hpp"
#include "orbit_statedef.hpp"
#include "math_basic.hpp"
#include "model_ephemeris_analytic.hpp"
#include "orbit_basic.hpp"

using namespace Eigen;

struct AtmosphereData
{
	double h_km;
	double rho_kgm3;
	double H_km;
	double dh_km;//h的间隔，之前命名为step，更名为dh
};

struct ReferenceDensity
{
	//处初始近地点的参考密度、标高等，归一化
	double cos_i0, sin_i0, sin_w0, rp0 = -1.0, hp0, sigma0;
	double rho, H = -1.0, Hrate;
	bool isSet = false;
};

/**
 * @brief 读取分层大气密度表及表头注明的空间天气条件。
 * @param[in] atm_path 密度表文件路径；数值记录依次为参考高度 km、密度 kg/m³、标高 km。
 * @param[in,out] atm_list 在数组末尾追加读到的分层记录；每条记录同时保存高度间距 dh_km。
 * @param[out] dh 从表头读取的高度间距，单位 km；读取失败时不保证有有效值。
 * @param[out] ap 可选接收表头的 Ap 指数，无量纲；传入 nullptr 时忽略，缺失或无效时写入 NaN。
 * @param[out] f107 可选接收表头的 F10.7 数值，按文件原标度保存；传入 nullptr 时忽略，缺失或无效时写入 NaN。
 * @return 成功返回 1；无法打开文件或数值记录读取失败时返回 -1。
 * @details Ap、F10.7 仅说明本表采用的固定条件，不改变密度值、插值方式或参考密度计算。
 */
int Read_AtmosphereDensity(const std::string& atm_path, std::vector<AtmosphereData>& atm_list, double& dh, double* ap = nullptr, double* f107 = nullptr);

//设置大气参考密度
void Set_ReferenceDensity(const std::vector<AtmosphereData>& atm_list,
	const orbit::Keplerian& initKep_norm,
	ReferenceDensity& refDnesity);

//计算实时位置的大气密度，基于参考密度
double Get_AtmosphereDensity(
	const double* cart_norm,//实时的位置速度，需要计算成根数
	const ReferenceDensity& refDensity,
	double MJD_UT1,//用于计算太阳位置，计算密度周日峰
	const std::vector<AtmosphereData>& atm_list);

//旧重载函数，直接在内部计算参考密度，繁琐过时，暂且保留，不推荐使用
double Get_AtmosphereDensity(
	const double* rv_Norm,//实时的位置速度，需要计算成根数
	const double* initKepArray_norm,//初始时刻的开普勒根数，归一化
	double MJD_UT1,//用于计算太阳位置，计算密度周日峰
	const std::vector<AtmosphereData>& atm_list);

//旧重载函数，直接在内部计算参考密度，繁琐过时，暂且保留，不推荐使用
double Get_AtmosphereDensity(const double* rv_Norm,
	const orbit::Keplerian& initKep_norm,
	double MJD_UT1,
	const std::vector<AtmosphereData>& atm_list);



#endif //CDC_MODEL_ATMOSPHERE_HPP
