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
    bool m_zero_fallback = false; // 显式启用后，覆盖范围外采用零定向参数。
    mutable size_t m_zero_fallback_count = 0;
 public:
	EOP();
	~EOP();
	//读进来的单位是 s as，不做调整！！
	void Read_EOP(Overwatch& tracer, std::string& fn_EOP);
    /**
     * @brief 设置覆盖范围外的零 EOP 策略，并重置调用计数。
     * @param[in] enabled true 时令缺测时刻的 dut1、xp、yp、dx00、dy00 为零，DAT 仍取 SOFA 闰秒表；false 恢复严格覆盖检查。
     */
    void SetZeroFallback(bool enabled);

    /**
     * @brief 返回自上次设置策略以来采用零 EOP 的查询次数。
     * @return 采用零定向参数的查询次数，包含积分器内部评估时刻。
     */
    size_t ZeroFallbackCount() const { return m_zero_fallback_count; }

    /**
     * @brief 查询 UTC 时刻的地球定向参数，按显式策略处理缺测。
     * @param[in,out] tracer 记录覆盖、时间转换错误及首次零 EOP 事件。
     * @param[in] MJD_UTC 查询时刻，单位 UTC MJD。
     * @param[out] eopNeed 输出 dut1、dat（s）及 xp、yp、dx00、dy00（rad），并记录查询时刻；失败时保持原值。
     * @param[in] fit_or_not 0 取当日零时记录，其余值做日间线性插值；整日零时无需下一日记录。
     * @details 默认要求连续每日记录。启用 SetZeroFallback 后，覆盖外或未载入资料时定向参数置零，DAT 从 SOFA 获取。
     */
	void Get_EOP(Overwatch& tracer, double MJD_UTC, aEOP& eopNeed, int fit_or_not = 1) const;

};


#endif // CDC_EARTH_EOP_HPP
