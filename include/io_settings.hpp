#ifndef CDC_IO_SETTINGS_HPP
#define CDC_IO_SETTINGS_HPP

#define DEFAULT_SIMPLEINPUT_DELIMITER               '='
#define DEFAULT_SIMPLEINPUT_COMMENTCHAR             "#!*"
#define DEFAULT_SIMPLEINPUT_BLANKCHAR               "\t\n\r "
#define DEFAULT_SIMPLEINPUT_WARN_KEYWORD_NOTFOUND   true
#define WARNING_SIMPLEINPUT_KEYWORD_NOT_FOUND   1
#include <string>
#include <vector>
#include "utils_overwatch.hpp"

namespace io
{

	struct SimpleRecord
	{
		std::string kw, va;//keyword = value
	};

	class SimpleInput
	{
		// This class defines a simple input file, with various member variables and functions.
		// Users can use ReadSimpleInput to fill the variables with the simple input file.
		// Input values can be accessed using various Get** functions.
		// Warnings are associated with flags which determine whether warnings are raised at all.
	 private:
		std::vector<SimpleRecord> m_rec;//用于记录SimpleRecord kw=va 的vector
		std::string m_fn;//文件路径
		char m_delimiter = DEFAULT_SIMPLEINPUT_DELIMITER;//界定符号，默认为等号
		std::string m_commentchar = DEFAULT_SIMPLEINPUT_COMMENTCHAR;//注释符号，默认为#!*
		std::string m_blankchar = DEFAULT_SIMPLEINPUT_BLANKCHAR;//空白符号，\t \n \r
		bool m_warnkeywordnotfound = DEFAULT_SIMPLEINPUT_WARN_KEYWORD_NOTFOUND;

	 public:
        /**
         * @brief 读取关键词与文本值，并在成功后替换当前记录。
         * @param[in,out] tracer 记录文件读取错误，以及严格模式的格式、空键和重复键错误。
         * @param[in] fn_input 输入配置文件路径；本函数保持该字符串不变。
         * @param[in] strict false 按普通 SimpleInput 规则读取；true 启用严格配置格式检查。
         * @details 普通模式忽略空行、以 #、! 或星号开头的注释行和无分隔符行，保留同名记录及行尾文字。
         * 严格模式只将 # 及其后文字视为注释，要求非空行有分隔符且键非空，拒绝重复键；值允许为空。
         * 两种模式都去除键值边缘空白。读取失败时保留此前的记录；具体关键词的类型和取值由读取接口及调用方检查。
         */
        void ReadSimpleInput(bastion::Overwatch& tracer, std::string& fn_input, bool strict = false);

		/* -------- functions below return various values by given keyword --------- */

        // 三种类型均保留返回值和接收变量两类接口；缺项返回零或空字符串。
        // 数值解析沿用字符串流规则，缺项警告由原有开关控制。
        double GetDouble(bastion::Overwatch& tracer, const std::string& kw) const;
		void GetDouble(bastion::Overwatch& tracer, const std::string& kw, double& x) const;

		int GetInteger(bastion::Overwatch& tracer, const std::string& kw) const;
		void GetInteger(bastion::Overwatch& tracer, const std::string& kw, int& x) const;

		std::string GetString(bastion::Overwatch& tracer, const std::string& kw) const;
		void GetString(bastion::Overwatch& tracer, const std::string& kw, std::string& x) const;

        /**
         * @brief 读取完整且有限的浮点配置值。
         * @param[in,out] tracer 记录关键词缺失、文本格式及取值错误，失败时保留调用链。
         * @param[in] kw 待读取的关键词，区分大小写。
         * @param[out] x 成功时写入解析值；失败时保持传入值。
         * @details 仅当整段文本都可解析为有限 double 时写回；空值、缺项、尾随文字、溢出及 NaN/Inf 均报错。数值的物理单位由关键词定义。
         */
        void GetDoubleStrict(bastion::Overwatch& tracer, const std::string& kw, double& x) const;

        /**
         * @brief 读取无小数损失且在 int 范围内的整数配置值。
         * @param[in,out] tracer 记录关键词缺失、文本格式及取值错误，失败时保留调用链。
         * @param[in] kw 待读取的关键词，区分大小写。
         * @param[out] x 成功时写入解析值；失败时保持传入值。
         * @details 先按严格浮点规则解析整段文本，再检查整数性和 int 范围；2、2.0、2e0 均可表示整数，2.5 报错。
         */
        void GetIntegerStrict(bastion::Overwatch& tracer, const std::string& kw, int& x) const;

        /**
         * @brief 读取明确的布尔开关配置值。
         * @param[in,out] tracer 记录关键词缺失、文本格式及取值错误，失败时保留调用链。
         * @param[in] kw 待读取的关键词，区分大小写。
         * @param[out] x 成功时写入解析值；失败时保持传入值。
         * @details 接受 true/1/on 和 false/0/off，区分大小写；空值、缺项及其他文本均报错。
         */
        void GetBoolean(bastion::Overwatch& tracer, const std::string& kw, bool& x) const;

		void GetRecord(bastion::Overwatch& tracer, const size_t idx, SimpleRecord& rec) const;

		/* -------- functions below provide auxillary functionality --------- */

		void SetDelimiter(bastion::Overwatch& tracer, char delimiter);
		void SetWarning(bastion::Overwatch& tracer, int warningId, bool flag = true);
		size_t Count() const;

        // 查询可选关键词是否存在，不产生“缺少关键词”警告；值为空仍算存在。
        bool HasKeyword(const std::string& keyword) const;

	};

}

#endif // CDC_IO_SETTINGS_HPP
