#include "earth_eop.hpp"
#include "io_common.hpp"
#include "utils_constants.hpp"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sofa.h>

EOP::EOP()
{
	m_isSet = false;
}

EOP::~EOP()
{
	// do nothing
}

void EOP::Read_EOP(Overwatch& tracer, std::string& fn_EOP)
{
	// 240204 读EOP 到vector struct中
	// 240813 改为读取celestrak整合发布的EOP文件，格式有区别，但内容还是一样的
	// 241024 封装为类函数

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	if (!m_isSet)
	{
		// 没读过，读
		std::ifstream ifs_EOP(fn_EOP);
		if (!ifs_EOP.is_open())
		{
			std::string msg = "Failed to open EOP ifs_EOP: " + fn_EOP + " !";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, EventType::ERROR);
			return;
		}

		std::vector<aEOP> records;
		std::string line;
		bool data_section_found = false;

		while (std::getline(ifs_EOP, line))
		{
			// 检查是否到达数据部分
			if (line.find("BEGIN OBSERVED") != std::string::npos || line.find("BEGIN PREDICTED") != std::string::npos)
			{
				data_section_found = true;
				continue;
			}

			if (line.find("END OBSERVED") != std::string::npos || line.find("END PREDICTED") != std::string::npos)
			{
				data_section_found = false;
				continue;
			}

			if (!data_section_found)
			{
				continue; // 跳过非数据部分
			}

			line = strip(line);
			if (line.empty() || line[0] == '#')
			{
				continue; // 跳过空行和注释行
			}

			std::istringstream iss(line);
			int year, month, day;
			aEOP eop;
			double lod, dPsi, dEpsilon;

			// 解析每一行的EOP数据
			if (!(iss >> year >> month >> day >> eop.MJD_UTC >> eop.xp >> eop.yp >> eop.dut1 >> lod >> dPsi >> dEpsilon
					  >> eop.dx00 >> eop.dy00 >> eop.dat))
			{
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Malformed EOP record: " + line, EventType::ERROR);
                return;
            }
            if (!std::isfinite(eop.MJD_UTC) || eop.MJD_UTC != std::floor(eop.MJD_UTC) ||
                !std::isfinite(eop.xp) || !std::isfinite(eop.yp) || !std::isfinite(eop.dut1) ||
                !std::isfinite(eop.dx00) || !std::isfinite(eop.dy00) || !std::isfinite(eop.dat) ||
                (!records.empty() && eop.MJD_UTC != records.back().MJD_UTC + 1))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "EOP records must be finite, ordered and daily without gaps.", EventType::ERROR);
                return;
            }
            records.push_back(eop);
		}
        if (ifs_EOP.bad() || records.empty())
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "EOP input is empty or unreadable.", EventType::ERROR);
            return;
        }
        mEOP.swap(records);
		ifs_EOP.close();
		m_isSet = true;
	}
	else
	{
		// do nothing
	}

	tracer.CheckOut(__func__);
}

void EOP::Get_EOP(Overwatch& tracer, double MJD_UTC, aEOP& eopNeed, int fit_or_not) const
{
	/*
	 * fit_or_not 默认为1;为0时，不拟合，直接取0h的值，否则做线性拟合
	 * MJD_UTC 需要读取EOP的UTC
	 * eopNeed 需要的EOP  单位 s rad
	 * 其中会令eopNeed.MJD_UTC=MJD_UTC,
	 * 因此传入参数的时候需要注意最好不要循环，除非事先定义了eopNeed中的UTC
	 */

	tracer.CheckIn(__FILE__, __func__, __LINE__);

    const double day = std::floor(MJD_UTC);
    if (!std::isfinite(MJD_UTC) || mEOP.empty() || day < mEOP.front().MJD_UTC || day > mEOP.back().MJD_UTC) {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid EOP coverage.", EventType::ERROR); return;
    }
    const double offset = day - mEOP.front().MJD_UTC;
    if (offset < 0 || offset >= static_cast<double>(mEOP.size())) {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid EOP coverage.", EventType::ERROR); return;
    }
    const size_t index = static_cast<size_t>(offset);
    const double MJD_int_part = day;
    const bool interpolate = fit_or_not != 0 && MJD_UTC != day;
    if (mEOP[index].MJD_UTC != day || (interpolate && (index+1 >= mEOP.size() || mEOP[index+1].MJD_UTC != day+1))) {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid EOP coverage.", EventType::ERROR); return;
    }

	if (!interpolate)//不需要拟合
	{
		eopNeed = mEOP[index];

		eopNeed.MJD_UTC = MJD_UTC;//修改一下，还是包含小数部分好
		// 已经从celestrak版本里读了，DAT是直接读进list的，不需要再额外修改
		//EOPs_needed.dat = DAT_2024;

		eopNeed.xp *= DAS2R;
		eopNeed.yp *= DAS2R;
		eopNeed.dx00 *= DAS2R;
		eopNeed.dy00 *= DAS2R;
	}
	else//需要拟合
	{
		double xp, yp, dut1, dat, dx00, dy00;  // utc为数组
		double MJD_frac_part = MJD_UTC - MJD_int_part;//MJD的小数部分

		//x=x0+kt
		eopNeed.MJD_UTC = MJD_UTC;
		eopNeed.xp =
			mEOP[index].xp + MJD_frac_part * (mEOP[index + 1].xp - mEOP[index].xp);
		eopNeed.yp =
			mEOP[index].yp + MJD_frac_part * (mEOP[index + 1].yp - mEOP[index].yp);
		eopNeed.dut1 =
			mEOP[index].dut1 + MJD_frac_part * (mEOP[index + 1].dut1 - mEOP[index].dut1);
		eopNeed.dx00 =
			mEOP[index].dx00 + MJD_frac_part * (mEOP[index + 1].dx00 - mEOP[index].dx00);
		eopNeed.dy00 =
			mEOP[index].dy00 + MJD_frac_part * (mEOP[index + 1].dy00 - mEOP[index].dy00);

		// 修改成读取celestrak后，DAT没有赋值，加上！
		// FIXME:240928 DAT就先不插值了，很少会遇到碰到闰秒的时候吧？？真是这样的话以后得改成使用TT进行积分了！留意一下。
		// EOPs_needed.dat = DAT_2024;
		eopNeed.dat = mEOP[index].dat;
		eopNeed.xp *= DAS2R;
		eopNeed.yp *= DAS2R;
		eopNeed.dx00 *= DAS2R;
		eopNeed.dy00 *= DAS2R;
	}
	tracer.CheckOut(__func__);
}
