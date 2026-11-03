/*****************************************************************************
*        *****    jpl planetary and lunar ephemerides    *****     C ver.1.6 *
******************************************************************************
* program testeph                                                            *
*                                                                            *
*                                                                            *
* Testeph tests the jpl ephemeris reading and interpolating routine using    *
* examples computed from the original ephemeris.                             *
*                                                                            *
* Testeph contains the reading and interpolating subroutines that are of     *
* eventual interest to the user.  Once testeph is working correctly, the     *
* user can extract those subroutines and the installation process is         *
* complete.                                                                  *
*                                                                            *
* You must allow access to "testpo.xxx" to testeph.                          *
* "testpo.xxx" is the specially formatted text file that contains the test   *
* cases for the ephemeris, dexxx.                                            *
*                                                                            *
* After the initial identifying text which is concluded by an "EOT" in       *
* columns 1-3, the test file contains the following quantities:              *
*                                                                            *
*     JPL ephemeris number                                                   *
*     calendar date                                                          *
*     julian ephemeris date                                                  *
*     target number (1-mercury, ...,3-earth, ,,,9-pluto, 10-moon, 11-sun,    *
*                    12-solar system barycenter, 13-earth-moon barycenter    *
*                    14-nutations, 15-librations)                            *
*     center number (same codes as target number)                            *
*     coordinate number (1-x, 2-y, ... 6-zdot)                               *
*     coordinate  [au, au/day].                                              *
*                                                                            *
* For each test case input, testeph                                          *
*                                                                            *
*     - computes the corresponding state from data contained                 *
*       in dexxx,                                                            *
*                                                                            *
*     - compares the two sets,                                               *
*                                                                            *
*     - writes an error message if the difference between                    *
*       any of the state components is greater than 10**(-13).               *
*                                                                            *
*     - writes state and difference information for every npt'th             *
*       test case processed.                                                 *
*                                                                            *
*                                                                            *
*  This program was written in standard fortran-77 and it was manually       *
*  translated to the C language by Piotr A. Dybczynski (dybol@amu.edu.pl).   *
*                                                                            *
*  This is version 1.6 of this C translation, use jplbin.h version 1.6       *
*
******************************************************************************
*                 Last modified: April 18, 2021 by PAD                       *
******************************************************************************
**   PAD: dr. Piotr A. Dybczynski,               e-mail: dybol@amu.edu.pl   **
**  Astronomical Observatory Institute, A.Mickiewicz Univ., Poznan, Poland  **
**       https://apollo.astro.amu.edu.pl/PAD/pmwiki.php?n=Dybol.JPLEph      **
*****************************************************************************/

#include "model_ephemeris_jpl.hpp"

// 240924 使用TZB时期的cpp版本覆盖掉原本的c版本，并且增加了对多种编译器的支持，对比之前的c版本没问题。
// 定义全局变量
JPLEPH_R1_Struct JPLEPH_R1;
JPLEPH_R2_Struct JPLEPH_R2;

int JPLEPH_BARY = 0; // 函数中需要在外部定义的变量
FILE* JPLEPH_F1 = nullptr; // JPL eph文件指针

double JPLEPH_PVSUN[6]; // 太阳位置速度
int JPLEPH_KM = 1; // KM=1时，单位为KM KM/s；为0时，单位为AU AU/day 默认为1

int Read_JPLEPH(std::string JPL_ephfile)
{
    // The underlying interpolator caches record indices: switching files in one
    // process is unsupported. A failed first load remains retryable.
    static std::string loaded_file;
    if (JPLEPH_F1 != nullptr) return loaded_file == JPL_ephfile ? 2 : -1;
    JPLEPH_F1 = std::fopen(JPL_ephfile.c_str(), "rb");
    if (!JPLEPH_F1) {
        std::cerr << "Cannot open JPL ephemeris: " << JPL_ephfile << std::endl;
        return -1;
    }
    char nams[JPLEPH_NMAX][6];
    double ss[3], vals[JPLEPH_NMAX];
    int nvs = -1;
    constan(nams, vals, ss, &nvs);
    if (nvs < 0) {
        std::cerr << "Invalid/truncated JPL header: " << JPL_ephfile << std::endl;
        std::fclose(JPLEPH_F1);
        JPLEPH_F1 = nullptr;
        return -1;
    }
    loaded_file = JPL_ephfile;
    return 1;
}

int Get_JPLEPH(const double MJD_TDB, int n_target, int n_center, double* r_v_Norm)
{
	if (!JPLEPH_F1) return -1;
	// 260209 严格来说，JPL是太阳系的历表，因此应该用太阳系的质心动力学时，TDB
	// 之前传入的MJD_TT 严格来说是错误的，但是在实际使用中，TT和TDB的差异非常小（通常在1.6毫秒范围内），对于大多数应用来说，这种差异是可以忽略的，因此之前传入MJD_TT虽然不完全正确，但在实际应用中不会导致显著的误差。
	// 现已经修改为传入MJD_TDB，确保使用正确的时间标准。
	pleph(MJD_TDB + START_OF_MJD, n_target, n_center, r_v_Norm);
	r_v_Norm[0] *= KM_TO_NORMALIZED_LENGTH;
	r_v_Norm[1] *= KM_TO_NORMALIZED_LENGTH;
	r_v_Norm[2] *= KM_TO_NORMALIZED_LENGTH;
	r_v_Norm[3] *= KMPS_TO_NORMALIZED_VELOCITY;
	r_v_Norm[4] *= KMPS_TO_NORMALIZED_VELOCITY;
	r_v_Norm[5] *= KMPS_TO_NORMALIZED_VELOCITY;

	return 1;
}

int Get_JPLEPH_ronly(double MJD_TDB, int n_target, int n_center, double* r_Norm)
{
	// 260209 严格来说，JPL是太阳系的历表，因此应该用太阳系的质心动力学时，TDB
	// 之前传入的MJD_TT 严格来说是错误的，但是在实际使用中，TT和TDB的差异非常小（通常在1.6毫秒范围内），对于大多数应用来说，这种差异是可以忽略的，因此之前传入MJD_TT虽然不完全正确，但在实际应用中不会导致显著的误差。
	// 现已经修改为传入MJD_TDB，确保使用正确的时间标准。
	double temp_r_v_Norm[6];
	pleph(MJD_TDB + START_OF_MJD, n_target, n_center, temp_r_v_Norm);
	r_Norm[0] = temp_r_v_Norm[0] * KM_TO_NORMALIZED_LENGTH;
	r_Norm[1] = temp_r_v_Norm[1] * KM_TO_NORMALIZED_LENGTH;
	r_Norm[2] = temp_r_v_Norm[2] * KM_TO_NORMALIZED_LENGTH;

	return 1;
}

int Get_Libration(const double MJD_TDB, double& Omega_p, double& i_s, double& Lambda)
{
	if (!JPLEPH_F1) return -1;
	// 260210 参考汤老师的代码 获取月球天平动的三个欧拉角，单位弧度
	// euler angles are interpreted in the order of 0~2: Omega', i_s, Lambda
	double rrd[6];
	// 设定ntar=15，ncent=0，获取月球天平动的三个欧拉角
	pleph(MJD_TDB + START_OF_MJD, 15, 0, rrd);

	Omega_p = rrd[0];
	i_s = rrd[1];
	Lambda = rrd[2];
	return 1;
}

/****************************************************************************/
/*****************************************************************************
**                         pleph(et,ntar,ncent,rrd)                         **
**                 给出JD(TDB)=et时刻，ntar相对ncent的位置速度                   **
******************************************************************************
**                                                                          **
**    This subroutine reads the jpl planetary ephemeris                     **
**    and gives the position and velocity of the point 'ntarg'              **
**    with respect to 'ncent'.                                              **
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**      et = (double) julian ephemeris date at which interpolation          **
**           is wanted.                                                     **
**                                                                          **
**    ntarg = integer number of 'target' point.                             **
**                                                                          **
**    ncent = integer number of center point.                               **
**                                                                          **
**    The numbering convention for 'ntarg' and 'ncent' is:                  **
**                                                                          **
**            1 = mercury           8 = neptune                             **
**            2 = venus             9 = pluto                               **
**            3 = earth            10 = moon                                **
**            4 = mars             11 = sun                                 **
**            5 = jupiter          12 = solar-system barycenter             **
**            6 = saturn           13 = earth-moon barycenter               **
**            7 = uranus           14 = nutations (longitude and obliq)     **
**                                 15 = librations, if on eph. file         **
**                                                                          **
**            (If nutations are wanted, set ntarg = 14.                     **
**             For librations, set ntarg = 15. set ncent= 0)                **
**                                                                          **
**     rrd = output 6-element, double array of position and velocity        **
**           of point 'ntarg' relative to 'ncent'. The units are au and     **
**           au/day. For librations the units are radians and radians       **
**           per day. In the case of nutations the first four words of      **
**           rrd will be set to nutations and rates, having units of        **
**           radians and radians/day.                                       **
**                                                                          **
**           The option is available to have the units in km and km/sec.    **
**			单位说明
**           for this, set km=TRUE at the beginning of the program.         **
*****************************************************************************/
void pleph(double et, int ntarg, int ncent, double rrd[])
{
    pleph(et, 0.0, ntarg, ncent, rrd);
}

/**
 * @brief 用拆分儒略日查询 JPL 天体相对状态。
 * @param[in] jd1 TDB 儒略日的整数日或参考日，与 jd2 相加为完整 JD。
 * @param[in] jd2 相对 jd1 的 TDB 小数日，独立保存亚秒时间精度。
 * @param[in] ntarg JPL 目标天体编号，1—13；14、15 分别查询章动和天平动。
 * @param[in] ncent JPL 中心天体编号，12 为太阳系质心。
 * @param[out] rrd 六维相对位置和速度；默认 JPLEPH_KM=1 时为 km、km/s，角运动为 rad、rad/s。
 * @details 使用已载入的历表，将两部分日期直接传给 state 的插值过程。
 */
void pleph(double jd1, double jd2, int ntarg, int ncent, double rrd[])
{
	double et2[2], pv[13][6];
	/* pv is the position/velocity array
                             NUMBERED FROM ZERO: 0=Mercury,1=Venus,...
                             8=Pluto,9=Moon,10=Sun,11=SSBary,12=EMBary
                             First 10 elements (0-9) are affected by state(),
                             all are adjusted here.                         */
	int bsave, i, k;
	int list[12];
	/* list is a vector denoting, for which "body"
						   ephemeris values should be calculated by state():
						   0=Mercury,1=Venus,2=EMBary,...,8=Pluto,
						   9=geocentric Moon, 10=nutations in long. & obliq.
						   11= lunar librations  */
	/*  initialize et2 for 'state' and set up component count   */
	et2[0] = jd1;
	et2[1] = jd2;

	for (i = 0; i < 6; ++i) rrd[i] = 0.0;

	if (ntarg == ncent) return;

	for (i = 0; i < 12; ++i) list[i] = 0;


	/*   check for nutation call    */
	if (ntarg == 14)
	{
		if (JPLEPH_R1.r1.ipt[11][1] > 0)
		{
			list[10] = 2;
			state(et2, list, pv, rrd);
		}
		else puts("***** no nutations on the ephemeris file  ******\n");
		return;
	}

	/*  check for librations   */
	if (ntarg == 15)
	{
		if (JPLEPH_R1.r1.lpt[1] > 0)
		{
			list[11] = 2;
			state(et2, list, pv, rrd);
			for (i = 0; i < 6; ++i) rrd[i] = pv[10][i];
		}
		else puts("*****  no librations on the ephemeris file  *****\n");
		return;
	}

	/*  force barycentric output by 'state'     */
	bsave = JPLEPH_BARY;
	JPLEPH_BARY = 1;

	/*  set up proper entries in 'list' array for state call     */
	for (i = 0; i < 2; ++i)
	{
		k = ntarg - 1;
		if (i == 1) k = ncent - 1;
		if (k <= 9) list[k] = 2;
		if (k == 9) list[2] = 2;
		if (k == 2) list[9] = 2;
		if (k == 12) list[2] = 2;
	}

	/*   make call to state   */
	state(et2, list, pv, rrd);

	/* Solar System barycentric Sun state goes to pv[10][] */
	if (ntarg == 11 || ncent == 11)
		for (i = 0; i < 6; ++i)
			pv[10][i] = JPLEPH_PVSUN[i];

	/* Solar System barycentric EMBary state:  */
	if (ntarg == 12 || ncent == 12)
		for (i = 0; i < 6; ++i)
			pv[11][i] = 0.0;

	/* if moon from earth or earth from moon ..... */
	if (ntarg == 13 || ncent == 13)
		for (i = 0; i < 6; ++i)
			pv[12][i] = pv[2][i];

	if ((ntarg * ncent) == 30 && (ntarg + ncent) == 13)
		for (i = 0; i < 6; ++i)
			pv[2][i] = 0.0;
	else
	{
		if (list[2] == 2)/* calculate earth state from EMBary */
			for (i = 0; i < 6; ++i)
				pv[2][i] -= pv[9][i] / (1.0 + JPLEPH_R1.r1.emrat);

		if (list[9] == 2)/* calculate Solar System barycentric moon state */
			for (i = 0; i < 6; ++i)
				pv[9][i] += pv[2][i];
	}

	for (i = 0; i < 6; ++i)
		rrd[i] = pv[ntarg - 1][i] - pv[ncent - 1][i];
	JPLEPH_BARY = bsave;
}


/*****************************************************************************
**                     interp(buf,t,ncf,ncm,na,ifl,pv)                      **
******************************************************************************
**                                                                          **
**    this subroutine differentiates and interpolates a                     **
**    set of chebyshev coefficients to give position and velocity           **
**                                                                          **
**    calling sequence parameters:                                          **
**                                                                          **
**      input:                                                              **
**                                                                          **
**        buf   1st location of array of d.p. chebyshev coefficients        **
**              of position                                                 **
**                                                                          **
**          t   t[0] is double fractional time in interval covered by       **
**              coefficients at which interpolation is wanted               **
**              (0 <= t[0] <= 1).  t[1] is dp length of whole               **
**              interval in input time units.                               **
**                                                                          **
**        ncf   # of coefficients per component                             **
**                                                                          **
**        ncm   # of components per set of coefficients                     **
**                                                                          **
**         na   # of sets of coefficients in full array                     **
**              (i.e., # of sub-intervals in full interval)                 **
**                                                                          **
**         ifl  integer flag: =1 for positions only                         **
**                            =2 for pos and vel                            **
**                                                                          **
**                                                                          **
**      output:                                                             **
**                                                                          **
**        pv   interpolated quantities requested.  dimension                **
**              expected is pv(ncm,ifl), dp.                                **
**                                                                          **
*****************************************************************************/
void interp(double coef[], double t[2], int ncf, int ncm, int na, int ifl, double posvel[6])
{
	static double pc[18], vc[18];
	static int np = 2, nv = 3, first = 1;
	static double twot = 0.0;
	double dna, dt1, temp, tc, vfac, temp1;
	int l, i, j;

	if (first)
	{
		pc[0] = 1.0;
		pc[1] = 0.0;
		vc[1] = 1.0;
		first = 0;
	}

	dna = (double) na;
	modf(t[0], &dt1);
	temp = dna * t[0];
	l = (int) (temp - dt1);

	tc = 2.0 * (modf(temp, &temp1) + dt1) - 1.0;

	if (tc != pc[1])
	{
		np = 2;
		nv = 3;
		pc[1] = tc;
		twot = tc + tc;
	}

	if (np < ncf)
	{
		for (i = np; i < ncf; ++i) pc[i] = twot * pc[i - 1] - pc[i - 2];
		np = ncf;
	}

	for (i = 0; i < ncm; ++i)
	{
		posvel[i] = 0.0;
		for (j = ncf - 1; j >= 0; --j)
			posvel[i] = posvel[i] + pc[j] * coef[j + i * ncf + l * ncf * ncm];
	}

	if (ifl <= 1) return;

	vfac = (dna + dna) / t[1];
	vc[2] = twot + twot;
	if (nv < ncf)
	{
		for (i = nv; i < ncf; ++i) vc[i] = twot * vc[i - 1] + pc[i - 1] + pc[i - 1] - vc[i - 2];
		nv = ncf;
	}

	for (i = 0; i < ncm; ++i)
	{
		posvel[i + ncm] = 0.0;
		for (j = ncf - 1; j > 0; --j)
			posvel[i + ncm] = posvel[i + ncm] + vc[j] * coef[j + i * ncf + l * ncf * ncm];
		posvel[i + ncm] = posvel[i + ncm] * vfac;
	}
}


/****************************************************************************
****                       split(tt,fr)                                  ****
*****************************************************************************
****  this subroutine breaks a d.p. number into a d.p. integer           ****
****  and a d.p. fractional part.                                        ****
****                                                                     ****
****  calling sequence parameters:                                       ****
****                                                                     ****
****    tt = d.p. input number                                           ****
****                                                                     ****
****    fr = d.p. 2-word output array.                                   ****
****         fr(1) contains integer part                                 ****
****         fr(2) contains fractional part                              ****
****                                                                     ****
****         for negative input numbers, fr(1) contains the next         ****
****         more negative integer; fr(2) contains a positive fraction.  ****
****************************************************************************/
void split(double tt, double fr[2])
{
	fr[1] = modf(tt, &fr[0]);

	if (tt >= 0.0 || fr[1] == 0.0) return;

	fr[0] = fr[0] - 1.0;
	fr[1] = fr[1] + 1.0;
}

/*****************************************************************************
**                        state(et2,list,pv,nut)                            **
******************************************************************************
** This subroutine reads and interpolates the jpl planetary ephemeris file  **
 * 读取并且插值！
**                                                                          **
**    Calling sequence parameters:                                          **
**                                                                          **
**    Input:                                                                **
**                                                                          **
**        et2[] double, 2-element JED epoch at which interpolation          **
**              is wanted.  Any combination of et2[0]+et2[1] which falls    **
**              within the time span on the file is a permissible epoch.    **
**                                                                          **
**               a. for ease in programming, the user may put the           **
**                  entire epoch in et2[0] and set et2[1]=0.0               **
**                                                                          **
**               b. for maximum interpolation accuracy, set et2[0] =        **
**                  the most recent midnight at or before interpolation     **
**                  epoch and set et2[1] = fractional part of a day         **
**                  elapsed between et2[0] and epoch.                       **
**                                                                          **
**               c. as an alternative, it may prove convenient to set       **
**                  et2[0] = some fixed epoch, such as start of integration,**
**                  and et2[1] = elapsed interval between then and epoch.   **
**                                                                          **
**       list   12-element integer array specifying what interpolation      **
**              is wanted for each of the "bodies" on the file.             **
**                                                                          **
**                        list[i]=0, no interpolation for body i            **
**                               =1, position only                          **
**                               =2, position and velocity                  **
**                                                                          **
**              the designation of the astronomical bodies by i is:         **
**                                                                          **
**                        i = 0: mercury                                    **
**                          = 1: venus                                      **
**                          = 2: earth-moon barycenter                      **
**                          = 3: mars                                       **
**                          = 4: jupiter                                    **
**                          = 5: saturn                                     **
**                          = 6: uranus                                     **
**                          = 7: neptune                                    **
**                          = 8: pluto                                      **
**                          = 9: geocentric moon                            **
**                          =10: nutations in longitude and obliquity       **
**                          =11: lunar librations (if on file)              **
**                                                                          **
**                                                                          **
**    output:                                                               **
**                                                                          **
**    pv[][6]   double array that will contain requested interpolated       **
**              quantities.  The body specified by list[i] will have its    **
**              state in the array starting at pv[i][0]  (on any given      **
**              call, only those words in 'pv' which are affected by the    **
**              first 10 'list' entries (and by list(11) if librations are  **
**              on the file) are set.  The rest of the 'pv' array           **
**              is untouched.)  The order of components in pv[][] is:       **
**              pv[][0]=x,....pv[][5]=dz.                                   **
**                                                                          **
**              All output vectors are referenced to the earth mean         **
**              equator and equinox of epoch. The moon state is always      **
**              geocentric; the other nine states are either heliocentric   **
**              or solar-system barycentric, depending on the setting of    **
**              global variables (see below).                               **
**                                                                          **
**              Lunar librations, if on file, are put into pv[10][k] if     **
**              list[11] is 1 or 2.                                         **
**                                                                          **
**        nut   dp 4-word array that will contain nutations and rates,      **
**              depending on the setting of list[10].  the order of         **
**              quantities in nut is:                                       **
**                                                                          **
**                       d psi  (nutation in longitude)                     **
**                       d epsilon (nutation in obliquity)                  **
**                       d psi dot                                          **
**                       d epsilon dot                                      **
**                                                                          **
**    Global variables:                                                     **
**                                                                          **
**         KM   logical flag defining physical units of the output          **
**              states. KM = TRUE, km and km/sec                            **
**                         = FALSE, au and au/day                           **
**              default value = FALSE.  (Note that KM determines time unit  **
**              for nutations and librations. Angle unit is always radians.)**
**                                                                          **
**       BARY   logical flag defining output center.                        **
**              only the 9 planets (entries 0-8) are affected.              **
**                       BARY = TRUE =\ center is solar-system barycenter   **
**                            = FALSE =\ center is sun                      **
**              default value = FALSE                                       **
**                                                                          **
**      PVSUN   double, 6-element array containing                          **
**              the barycentric position and velocity of the sun.           **
**                                                                          **
*****************************************************************************/
void state(double et2[2], int list[12], double pv[][6], double nut[4])
{
	int i, j;
	static int ipt[13][3], first = 1;
	int nr;
	static int nrl = 0;
	double pjd[4];
	static double buf[NCOEFF];
	double s, t[2], aufac;
	double pefau[6];

	if (first)
	{
		first = 0;
		for (i = 0; i < 3; ++i)
		{
			for (j = 0; j < 12; ++j) ipt[j][i] = (int) JPLEPH_R1.r1.ipt[j][i];
			ipt[12][i] = (int) JPLEPH_R1.r1.lpt[i];
		}
	}

	s = et2[0] - 0.5;
	split(s, &pjd[0]);
	split(et2[1], &pjd[2]);
	pjd[0] = pjd[0] + pjd[2] + 0.5;
	pjd[1] = pjd[1] + pjd[3];
	split(pjd[1], &pjd[2]);
	pjd[0] = pjd[0] + pjd[2];

	if ((pjd[0] + pjd[3]) < JPLEPH_R1.r1.ss[0] || (pjd[0] + pjd[3]) > JPLEPH_R1.r1.ss[1])
	{
		puts("Requested JED not within ephemeris limits.\n");
		return;
	}

	nr = (int) ((pjd[0] - JPLEPH_R1.r1.ss[0]) / JPLEPH_R1.r1.ss[2]) + 2;
	if (pjd[0] == JPLEPH_R1.r1.ss[1]) nr = nr - 1;
	t[0] = (pjd[0] - ((1.0 * nr - 2.0) * JPLEPH_R1.r1.ss[2] + JPLEPH_R1.r1.ss[0]) + pjd[3]) / JPLEPH_R1.r1.ss[2];

	if (nr != nrl)
	{
		nrl = nr;
		fseek(JPLEPH_F1, nr * RECSIZE, SEEK_SET);
		fread(buf, sizeof(buf), 1, JPLEPH_F1);
	}

	if (JPLEPH_KM)
	{
		t[1] = JPLEPH_R1.r1.ss[2] * 86400.0;
		aufac = 1.0;
	}
	else
	{
		t[1] = JPLEPH_R1.r1.ss[2];
		aufac = 1.0 / JPLEPH_R1.r1.au;
	}

	interp(&buf[ipt[10][0] - 1], t, ipt[10][1], 3, ipt[10][2], 2, pefau);

	for (i = 0; i < 6; ++i) JPLEPH_PVSUN[i] = pefau[i] * aufac;

	for (i = 0; i < 10; ++i)
	{
		if (list[i] == 0) continue;

		interp(&buf[ipt[i][0] - 1], t, ipt[i][1], 3, ipt[i][2], list[i], pefau);

		for (j = 0; j < 6; ++j)
		{
			if (i < 9 && !JPLEPH_BARY)
				pv[i][j] = pefau[j] * aufac - JPLEPH_PVSUN[j];
			else
				pv[i][j] = pefau[j] * aufac;
		}
	}

	if (list[10] > 0 && ipt[11][1] > 0)
		interp(&buf[ipt[11][0] - 1], t, ipt[11][1], 2, ipt[11][2], list[10], nut);

	if (list[11] > 0 && ipt[12][1] > 0)
	{
		interp(&buf[ipt[12][0] - 1], t, ipt[12][1], 3, ipt[12][2], list[11], pefau);
		for (j = 0; j < 6; ++j) pv[10][j] = pefau[j];
	}
}


/****************************************************************************
**                        constan(nam,val,sss,n)                           **
*****************************************************************************
**                                                                         **
**    this function obtains the constants from the ephemeris file          **
**    calling sequence parameters (all output):                            **
**      char nam[][6] = array of constant names (max 6 characters each)    **
**      double val[] = array of values of constants                        **
**      double sss[3] = JED start, JED stop, step of ephemeris             **
**      int n = number of entries in 'nam' and 'val' arrays                **
** Note: we changed name of this routine because const is a reserved word  **
**       in the C language.                                                **
*****************************************************************************
**    external variables:                                                  **
**         struct R1 and R2 (first two records of ephemeris)               **
**         defined in file:     jplbin.h                                   **
**         F1 = ephemeris binary file pointer (obtained from fopen() )     **
****************************************************************************/
void constan(char nam[][6], double val[], double sss[], int* n)
{
	int i, j;

    *n = -1;
    if (!JPLEPH_F1 ||
        std::fread(&JPLEPH_R1, sizeof(JPLEPH_R1), 1, JPLEPH_F1) != 1 ||
        std::fread(&JPLEPH_R2, sizeof(JPLEPH_R2), 1, JPLEPH_F1) != 1) return;
    const auto& header = JPLEPH_R1.r1;
    if (header.ncon < 0 || header.ncon > JPLEPH_NMAX ||
        !std::isfinite(header.ss[0]) || !std::isfinite(header.ss[1]) || !std::isfinite(header.ss[2]) ||
        header.ss[1] <= header.ss[0] || header.ss[2] <= 0 ||
        !std::isfinite(header.au) || header.au <= 0 ||
        !std::isfinite(header.emrat) || header.emrat <= 0) return;
    // DENUM selects record size, not a single allowable DE release.
    // Preserve DE430/440 etc. files using the same layout.
    int expected_ksize = 0;
    switch (header.numde) {
        case 200: expected_ksize = 1652; break;
        case 404: case 406: expected_ksize = 1456; break;
        case 403: case 405: case 421: case 422: case 430: case 431: case 440: case 441:
            expected_ksize = 2036; break;
        default: return;
    }
    if (expected_ksize != KSIZE) return;
    *n = header.ncon;

	for (i = 0; i < 3; ++i) sss[i] = JPLEPH_R1.r1.ss[i];

	for (i = 0; i < OLDMAX; ++i)
		for (j = 0; j < 6; ++j) nam[i][j] = JPLEPH_R1.r1.cnam[i][j];

	if (*n > OLDMAX)
		for (i = OLDMAX; i < *n; ++i)
			for (j = 0; j < 6; ++j) nam[i][j] = JPLEPH_R1.r1.cnam2[i - OLDMAX][j];

	for (i = 0; i < *n; ++i) val[i] = JPLEPH_R2.r2.cval[i];
}
