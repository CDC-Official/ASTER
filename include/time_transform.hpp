//
// Created by 陈酊存 on 2024/2/2.
//

#ifndef CDC_LIB_INCLUDE_TIME_TRANSFORM_HPP_
#define CDC_LIB_INCLUDE_TIME_TRANSFORM_HPP_

#include <iostream>
#include <vector>
#include <sofa.h>
#include "utils_constants.hpp"

struct TimesMJD
{
	double UTC;
	double UT1;
	double TT;
	double TAI;
	// newly added, for EBCRS time transformations
	double TCG;
	double TDB;
	double TCB;
};

int Calendar2MJD(const int& year, const int& month, const int& day, const double& day_fraction, double& MJD);

int Calendar2MJD(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& second, double& MJD);

int Calendar2JD(const int& year, const int& month, const int& day, const double& day_fraction, double& JD);

int Calendar2JD(const int& year, const int& month, const int& day, const int& hour, const int& minute, const double& second, double& JD);

int MJD2Calendar(const double& MJD, int& year, int& month, int& day, int& hour, int& minute, double& second);

double MJD2JD(double MJD);

double JD2MJD(double JD);

TimesMJD Utc2Else(double dut1, double MJD_UTC);

double UTC2TT(double MJD_UTC);

double UTC2TAI(double MJD_UTC);

double TT2TAI(double MJD_TT);

double TT2UTC(double MJD_TT);

/**
 * Convert TT to TCG using the IAU SOFA linear time-scale transformation.
 *
 * @param MJD_TT Input Terrestrial Time as Modified Julian Date.
 * @return Geocentric Coordinate Time as Modified Julian Date.
 *
 * @note This is a convenience wrapper around iauTttcg. It accepts and returns
 *       a single MJD value instead of SOFA's two-part Julian Date format.
 */
double TT2TCG(double MJD_TT);

/**
 * Convert TCG to TT using the IAU SOFA linear time-scale transformation.
 *
 * @param MJD_TCG Input Geocentric Coordinate Time as Modified Julian Date.
 * @return Terrestrial Time as Modified Julian Date.
 *
 * @note This is a convenience wrapper around iauTcgtt. It accepts and returns
 *       a single MJD value instead of SOFA's two-part Julian Date format.
 */
double TCG2TT(double MJD_TCG);

/**
 * Convert TDB to TCB using the IAU SOFA linear time-scale transformation.
 *
 * @param MJD_TDB Input Barycentric Dynamical Time as Modified Julian Date.
 * @return Barycentric Coordinate Time as Modified Julian Date.
 *
 * @note This is a convenience wrapper around iauTdbtcb. It accepts and returns
 *       a single MJD value instead of SOFA's two-part Julian Date format.
 */
double TDB2TCB(double MJD_TDB);

/**
 * Convert TCB to TDB using the IAU SOFA linear time-scale transformation.
 *
 * @param MJD_TCB Input Barycentric Coordinate Time as Modified Julian Date.
 * @return Barycentric Dynamical Time as Modified Julian Date.
 *
 * @note This is a convenience wrapper around iauTcbtdb. It accepts and returns
 *       a single MJD value instead of SOFA's two-part Julian Date format.
 */
double TCB2TDB(double MJD_TCB);

/**
 * Compute the geocentric TCB-TCG offset using the HF2002 approximation.
 *
 * @param JD Julian Date of the input epoch. The model is formally expressed
 *           as a function of TDB, while the IERS routine notes that TT may be
 *           used with no practical effect on the prediction accuracy.
 * @return Geocentric TCB-TCG in seconds.
 *
 * @note This function is a CDC_lib C++ adaptation of the IERS Conventions
 *       HF2002_IERS.F routine. It approximates the TE405 time ephemeris and
 *       includes the IERS linear c^-4 correction term.
 * @note This function only evaluates the geocentric term. The position term
 *       v_E dot (r-r_E) / c^2 is added by TCG2TCB and TCB2TCG.
 */
double GeocenterTCBMinusTCG_HF2002(double JD);

/**
 * Convert TCG to TCB at an arbitrary position near the Earth.
 *
 * @param MJD_TCG Input Geocentric Coordinate Time as Modified Julian Date.
 * @param earth_vx_mps Earth BCRS velocity x-component in m/s.
 * @param earth_vy_mps Earth BCRS velocity y-component in m/s.
 * @param earth_vz_mps Earth BCRS velocity z-component in m/s.
 * @param dx_m x-component of target position relative to the geocenter,
 *             r-r_E, in meters.
 * @param dy_m y-component of target position relative to the geocenter,
 *             r-r_E, in meters.
 * @param dz_m z-component of target position relative to the geocenter,
 *             r-r_E, in meters.
 * @return Barycentric Coordinate Time as Modified Julian Date.
 *
 * @note The caller is responsible for providing the Earth velocity and target
 *       geocentric position in SI units. This function does not read any
 *       ephemeris file and performs no unit normalization.
 * @note The geocentric TCB-TCG term is evaluated with
 *       GeocenterTCBMinusTCG_HF2002. The local position correction is
 *       v_E dot (r-r_E) / c^2.
 */
double TCG2TCB(double MJD_TCG,
               const double earth_vx_mps, const double earth_vy_mps, const double earth_vz_mps,
               const double dx_m, const double dy_m, const double dz_m);

/**
 * Convert TCB to TCG at an arbitrary position near the Earth.
 *
 * @param MJD_TCB Input Barycentric Coordinate Time as Modified Julian Date.
 * @param earth_vx_mps Earth BCRS velocity x-component in m/s.
 * @param earth_vy_mps Earth BCRS velocity y-component in m/s.
 * @param earth_vz_mps Earth BCRS velocity z-component in m/s.
 * @param dx_m x-component of target position relative to the geocenter,
 *             r-r_E, in meters.
 * @param dy_m y-component of target position relative to the geocenter,
 *             r-r_E, in meters.
 * @param dz_m z-component of target position relative to the geocenter,
 *             r-r_E, in meters.
 * @return Geocentric Coordinate Time as Modified Julian Date.
 *
 * @note The caller is responsible for providing the Earth velocity and target
 *       geocentric position in SI units. This function does not read any
 *       ephemeris file and performs no unit normalization.
 * @note The geocentric TCB-TCG term is evaluated with
 *       GeocenterTCBMinusTCG_HF2002. The local position correction is
 *       v_E dot (r-r_E) / c^2.
 */
double TCB2TCG(double MJD_TCB,
               const double earth_vx_mps, const double earth_vy_mps, const double earth_vz_mps,
               const double dx_m, const double dy_m, const double dz_m);

/**
 * Convert TDB to TT at a target position expressed in EBCRS.
 *
 * @param MJD_TDB Input Barycentric Dynamical Time as Modified Julian Date.
 * @param dx_m Target x-position in EBCRS, i.e. the x-component of r-r_E,
 *             in meters.
 * @param dy_m Target y-position in EBCRS, i.e. the y-component of r-r_E,
 *             in meters.
 * @param dz_m Target z-position in EBCRS, i.e. the z-component of r-r_E,
 *             in meters.
 * @return Terrestrial Time as Modified Julian Date at the target position.
 *
 * @note This function uses the already implemented chain
 *       TDB -> TCB -> TCG -> TT.
 * @note The Earth BCRS velocity is read internally from the loaded JPL
 *       ephemeris with Get_JPLEPH(MJD_TDB, 3, 12, ...), where target 3 is
 *       the Earth and center 12 is the solar-system barycenter in this
 *       wrapper.
 * @note The caller must load the JPL ephemeris with Read_JPLEPH before
 *       calling this function.
 */
double TDB2TT(double MJD_TDB,const double dx_m, const double dy_m, const double dz_m);


#endif //CDC_LIB_INCLUDE_TIME_TRANSFORM_HPP_
