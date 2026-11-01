#ifndef CDC_EARTH_EOP_HPP
#define CDC_EARTH_EOP_HPP

#include <string>
#include <vector>
#include "utils_overwatch.hpp"

using namespace bastion;

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 关于读取struct EOP的声明 @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

struct aEOP
{
	//全部初始化为0！注意单位不做硬性要求！可以为原始的s、as，也可以为s、rad
	double MJD_UTC = 0;
	double dut1 = 0;
	double dat = 0;
	double xp = 0;
	double yp = 0;
	double dx00 = 0;
	double dy00 = 0;
};

class EOP
{
 private:
	std::vector<aEOP> mEOP;//默认读进来的单位是 s as！
	bool m_isSet = false;
 public:
	EOP();
	~EOP();
	//读进来的单位是 s as，不做调整！！
	void Read_EOP(Overwatch& tracer, std::string& fn_EOP);
	// Output: s and rad. Records must be consecutive UTC days.
	// Interpolation needs the next day, except at exact midnight.
	// Error: tracer is set; output is unchanged. No extrapolation.
	void Get_EOP(Overwatch& tracer, double MJD_UTC, aEOP& eopNeed, int fit_or_not = 1) const;

};


#endif // CDC_EARTH_EOP_HPP
