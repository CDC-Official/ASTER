#ifndef CDC_MODEL_EPHEMERIS_JPL_HPP
#define CDC_MODEL_EPHEMERIS_JPL_HPP

#include "utils_constants.hpp"
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#define DENUM 440

#if   DENUM == 200
#define KSIZE 1652
#elif DENUM == 403 || DENUM == 405 || DENUM == 421 || DENUM == 422 || DENUM == 430 || DENUM == 431 || DENUM == 440 || DENUM == 441
#define KSIZE 2036
#elif DENUM == 404 || DENUM == 406
#define KSIZE 1456
#endif

#define NRECL 4
#define RECSIZE (NRECL * KSIZE)
#define NCOEFF (KSIZE / 2)
#define TRUE 1
#define FALSE 0

#define JPLEPH_NMAX 1000
#define OLDMAX 400

#pragma pack(push, 1)

struct rec1 {
	char ttl[3][84];
	char cnam[OLDMAX][6];
	double ss[3];
	int ncon;
	double au;
	double emrat;
	int ipt[12][3];
	int numde;
	int lpt[3];
	int rpt[3];
	int tpt[3];
	char cnam2[(JPLEPH_NMAX - OLDMAX)][6];
};

struct JPLEPH_R1_Struct {
	rec1 r1;
	char spare[RECSIZE - sizeof(rec1)];
};

struct rec2 {
	double cval[JPLEPH_NMAX];
};

struct JPLEPH_R2_Struct {
	rec2 r2;
	char spare[RECSIZE - sizeof(rec2)];
};

#pragma pack(pop)

// 全局变量声明
extern JPLEPH_R1_Struct JPLEPH_R1;
extern JPLEPH_R2_Struct JPLEPH_R2;

void constan(char nam[][6], double val[], double sss[], int* n);
void state(double et2[], int list[], double pv[][6], double nut[]);
void split(double tt, double fr[]);
void interp(double buf[], double t[], int ncf, int ncm, int na, int ifl, double pv[]);
void pleph(double et, int ntarg, int ncent, double rrd[]);

// Load: 1=new, 2=same path already loaded, -1=failure/different path.
// File switching is unsupported; load failure may be retried.
// Get wrappers return -1 before loading and leave output unchanged.
// 自己封装的代码
int Read_JPLEPH(std::string JPL_ephfile);
int Get_JPLEPH(const double MJD_TDB, int n_target, int n_center, double* r_v_Norm);
int Get_Libration(const double MJD_TDB, double& Omega_p, double& i_s, double& Lambda);

#endif // CDC_MODEL_EPHEMERIS_JPL_HPP
