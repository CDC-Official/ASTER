/**
 * @file io_settings.cpp
 * @brief 读取 SimpleInput 关键词配置并提供普通及严格的类型读取接口。
 * @details 普通数值接口接受前缀数值，缺项警告并返回零；严格接口要求完整有效值，失败时保留接收变量。
 * 配置文件读取先建立临时记录，成功后整体替换。关键词的物理单位和业务取值范围由调用方定义。
 */
#include "io_settings.hpp"
#include "io_common.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

/**
 * @brief 读取关键词与文本值，并在成功后替换当前记录。
 * @param[in,out] tracer 记录文件读取错误，以及严格模式的格式、空键和重复键错误。
 * @param[in] fn_input 输入配置文件路径；本函数保持该字符串不变。
 * @param[in] strict false 按普通 SimpleInput 规则读取；true 启用严格配置格式检查。
 * @details 普通模式忽略空行、以 #、! 或星号开头的注释行和无分隔符行，保留同名记录及行尾文字。
 * 严格模式只将 # 及其后文字视为注释，要求非空行有分隔符且键非空，拒绝重复键；值允许为空。
 * 两种模式都去除键值边缘空白。读取失败时保留此前的记录；具体关键词的类型和取值由读取接口及调用方检查。
 */
void io::SimpleInput::ReadSimpleInput(bastion::Overwatch& tracer, std::string& fn_input, bool strict)
{
    tracer.CheckIn(__FILE__, __func__, __LINE__);

    std::ifstream file(fn_input);

    if (!file.is_open())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open file: " + fn_input + " !", bastion::ERROR);
        return;
    }

    // 先收集本次记录，完整读取成功后才替换对象中的配置。
    std::vector<SimpleRecord> records;
    std::string line;
    std::size_t line_number = 0;

    while (std::getline(file, line))
    {
        ++line_number;

        // 严格模式允许行尾 # 注释，路径和数值都在去掉注释后解释。
        if (strict)
        {
            const std::size_t comment = line.find('#');

            if (comment != std::string::npos) line.erase(comment);
        }

        line = strip(line);

        if (line.empty()) continue;

        if (!strict && m_commentchar.find(line[0]) != std::string::npos) continue;

        const std::size_t delimiter_position = line.find(m_delimiter);

        if (delimiter_position == std::string::npos)
        {
            if (strict)
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Configuration line " + std::to_string(line_number) + " has no delimiter: " + fn_input, bastion::ERROR);
                return;
            }

            continue;
        }

        // SimpleInput 保存文本值，类型检查在 GetDoubleStrict 等读取接口中完成。
        const std::string key = strip(line.substr(0, delimiter_position));
        const std::string value = strip(line.substr(delimiter_position + 1));

        if (strict)
        {
            if (key.empty())
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__,
                    "Configuration key is empty at line " + std::to_string(line_number) + ": " + fn_input, bastion::ERROR);
                return;
            }

            for (const SimpleRecord& record : records)
            {
                if (record.kw == key)
                {
                    tracer.LogEvent(__FILE__, __func__, __LINE__, "Duplicate configuration key: " + key, bastion::ERROR);
                    return;
                }
            }
        }

        records.push_back(SimpleRecord{key, value});
    }

    if (file.bad())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed to read configuration: " + fn_input, bastion::ERROR);
        return;
    }

    m_rec.swap(records);
    m_fn = fn_input;
    tracer.CheckOut(__func__);
}

double io::SimpleInput::GetDouble(bastion::Overwatch& tracer, const std::string& kw) const
{
	// 两类接口使用相同的查找与解析过程，确保缺项返回值和警告行为一致。
	double value = 0.0;
	GetDouble(tracer, kw, value);
	return value;
}

void io::SimpleInput::GetDouble(bastion::Overwatch& tracer, const std::string& kw, double& x) const
{
	// 按记录顺序查找关键词，用字符串流读取开头的浮点数，允许数值后附带文本。
	// 同名记录解析失败后警告并归零，继续尝试后续记录，首次成功便停止。
	// 没有记录解析成功时按原有开关补充缺项警告；关闭警告仍返回零。

	tracer.CheckIn(__FILE__, __func__, __LINE__);

	double value = 0.0;
	bool isFound = false;

	for (const auto& record : m_rec)
	{    // 按原记录顺序查找，失败后继续尝试后续同名记录。
		// 直接比较关键字

		if (record.kw == kw)
		{    // 找到同名记录后尝试读取数值，失败时仍继续遍历。
			std::stringstream ss(record.va);  // 用 stringstream 解析字符串

			// 将字符串转换为 double
			if (ss >> value)
			{
				isFound = true;
				break;
			}
			else
			{
				// 当前记录不能读成浮点数，警告并归零，再尝试后续同名记录。
				std::string msg = "Failed to parse value for keyword: " + kw + ", zero is returned.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
				value = 0.0;
			}
		}
	}

	// 未读取到有效值时，沿用原版的缺项警告开关。
	if (!isFound && m_warnkeywordnotfound)
	{
		std::string msg = "Keyword not found: " + kw + ", zero is returned.";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
		value = 0.0;
	}

	tracer.CheckOut(__func__);
	x = value;
}

int io::SimpleInput::GetInteger(bastion::Overwatch& tracer, const std::string& kw) const
{
	// 两类接口使用相同的查找与解析过程，确保缺项返回值和警告行为一致。
	int value = 0;
	GetInteger(tracer, kw, value);
	return value;
}

void io::SimpleInput::GetInteger(bastion::Overwatch& tracer, const std::string& kw, int& x) const
{
	// 按记录顺序查找关键词，直接用整数流提取，允许数值后附带文本。
	// 值中含小数点时先警告，随后读取开头的整数部分，不做额外舍入。
	// 解析失败时警告并归零，再尝试后续同名记录；全部失败则按开关补充缺项警告。

	tracer.CheckIn(__FILE__, __func__, __LINE__);
	int value = 0;
	bool isFound = false;
	for (const auto& record : m_rec)
	{
		// 直接比较关键字
		if (record.kw == kw)
		{
			// 检查是否包含小数点，发出警告
			if (record.va.find('.') != std::string::npos)
			{
				std::string msg;
				msg = "Value for keyword '" + kw + "' contains a decimal, returning integer part only.";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
			}

			std::stringstream ss(record.va);  // 用 stringstream 解析字符串

			// 将字符串转换为整数
			if (ss >> value)
			{
				isFound = true;
				break;
			}
			else
			{
				std::string msg;
				msg = "Failed to parse value for keyword: " + kw + ", zero is returned!";
				tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
				value = 0;
			}
		}
	}

	// 未读取到有效值时，沿用原版的缺项警告开关。
	if (!isFound && m_warnkeywordnotfound)
	{

		std::string msg;
		msg = "Keyword not found: " + kw + ", zero is returned!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
		value = 0;

	}

	tracer.CheckOut(__func__);
	x = value;
}

std::string io::SimpleInput::GetString(bastion::Overwatch& tracer, const std::string& kw) const
{
	// 两类接口使用相同的查找与解析过程，确保缺项返回值和警告行为一致。
	std::string value = "";
	GetString(tracer, kw, value);
	return value;
}

void io::SimpleInput::GetString(bastion::Overwatch& tracer, const std::string& kw, std::string& x) const
{
	// 直接读取第一条同名记录的字符串，空值也表示已经找到关键词。
	// 没有匹配记录时按原有开关发出警告，并返回空字符串。

	tracer.CheckIn(__FILE__, __func__, __LINE__);
	std::string value = "";
	bool isFound = false;

	for (const auto& record : m_rec)
	{
		// 直接比较关键字
		if (record.kw == kw)
		{
			value = record.va;//va原本就是string类，因此不需要额外的类型判断
			isFound = true;
			break;
		}
	}

	// 未找到关键词时按开关提示，空字符串记录仍视为找到。
	if (!isFound && m_warnkeywordnotfound)
	{
		std::string msg = "Keyword not found: " + kw + ", returning empty string!";
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
	}

	x = value;
	tracer.CheckOut(__func__);
}


bool io::SimpleInput::HasKeyword(const std::string& keyword) const
{
    for (const auto& record : m_rec)
        if (record.kw == keyword) return true;
    return false;
}

void io::SimpleInput::GetRecord(bastion::Overwatch& tracer, const size_t idx, SimpleRecord& rec) const
{
	// ************************************************************************
	//    Function: GetRecord
	// Description: Get the SimpleRecord with index idx.
	//    Last Rev: Jingshi @ 2020/08/24
	//    Modified: CDC @ 2024/10/21
	//
	//       Calls: <vector>.size
	//   Called By: Called when needed.
	//
	//       Input: idx - int, index of requested object, this index should start from ONE (1), not 0.
	//      Output: rec - SimpleRecord object, by reference
	//                    holds requested SimpleRecord object if found, untouched if out of bound.
	//      Return: Void
	//
	//        Note: 1. If requested index exceeds collection size, a warning is popped.
	// ************************************************************************

	tracer.CheckIn(__FILE__, __func__, __LINE__);
	// GetRecord的索引从1开始而不是0
	if (idx > m_rec.size())//超过最大值
	{
		std::string msg;
		msg = "Requested index exceeds number of records, nothing returned!";
		//std::cerr << msg;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
		//return; warning应当直接到checkout
	}
	else if (idx <= 0)//小于等于0
	{
		std::string msg;
		msg = "Index should be positive, and it is counted from 1 instead of 0!";
		//std::cerr << msg;
		tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::ERROR);
		return;//error直接return
	}
	else//合法的
	{
		rec = m_rec[idx - 1];
	}

	tracer.CheckOut(__func__);
}

void io::SimpleInput::SetDelimiter(bastion::Overwatch& tracer, char delimiter)
{
	// 此函数基本没报错的可能
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	strcpy(&(this->m_delimiter), &delimiter);

	tracer.CheckOut(__func__);
}

void io::SimpleInput::SetWarning(bastion::Overwatch& tracer, int warningId, bool flag)
{
	// 不知道干啥用的，先写着
	tracer.CheckIn(__FILE__, __func__, __LINE__);

	switch (warningId)
	{
		case WARNING_SIMPLEINPUT_KEYWORD_NOT_FOUND:
			m_warnkeywordnotfound = flag;
			break;
		default:
			std::string msg;
			msg = "Unknown warning Id for SimpleInput. NOOP.";
			tracer.LogEvent(__FILE__, __func__, __LINE__, msg, bastion::EventType::WARNING);
			break;
	}

	tracer.CheckOut(__func__);
}

size_t io::SimpleInput::Count() const
{
	return m_rec.size();
}

/**
 * @brief 读取完整且有限的浮点配置值。
 * @param[in,out] tracer 记录关键词缺失、文本格式及取值错误，失败时保留调用链。
 * @param[in] kw 待读取的关键词，区分大小写。
 * @param[out] x 成功时写入解析值；失败时保持传入值。
 * @details 仅当整段文本都可解析为有限 double 时写回；空值、缺项、尾随文字、溢出及 NaN/Inf 均报错。数值的物理单位由关键词定义。
 */
void io::SimpleInput::GetDoubleStrict(bastion::Overwatch& tracer, const std::string& kw, double& x) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (!HasKeyword(kw))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Missing numeric keyword: " + kw, bastion::ERROR);
        return;
    }

    std::string text;
    GetString(tracer, kw, text);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    // 先解析到局部值，失败时保留接收变量；stod 的越界及格式异常转为 tracer 错误。
    double value;
    std::size_t used = 0;

    try
    {
        value = std::stod(text, &used);
    }
    catch (const std::invalid_argument&)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Invalid numeric value for keyword: " + kw, bastion::ERROR);
        return;
    }
    catch (const std::out_of_range&)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Numeric value out of range for keyword: " + kw, bastion::ERROR);
        return;
    }

    if (used != text.size() || !std::isfinite(value))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected one finite number for keyword: " + kw, bastion::ERROR);
        return;
    }

    x = value;

    tracer.CheckOut(__func__);
}

/**
 * @brief 读取无小数损失且在 int 范围内的整数配置值。
 * @param[in,out] tracer 记录关键词缺失、文本格式及取值错误，失败时保留调用链。
 * @param[in] kw 待读取的关键词，区分大小写。
 * @param[out] x 成功时写入解析值；失败时保持传入值。
 * @details 先按严格浮点规则解析整段文本，再检查整数性和 int 范围；2、2.0、2e0 均可表示整数，2.5 报错。
 */
void io::SimpleInput::GetIntegerStrict(bastion::Overwatch& tracer, const std::string& kw, int& x) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    double value;
    GetDoubleStrict(tracer, kw, value);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    if (std::floor(value) != value || value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max())
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected an integer in int range for keyword: " + kw, bastion::ERROR);
        return;
    }

    x = static_cast<int>(value);

    tracer.CheckOut(__func__);
}

/**
 * @brief 读取明确的布尔开关配置值。
 * @param[in,out] tracer 记录关键词缺失、文本格式及取值错误，失败时保留调用链。
 * @param[in] kw 待读取的关键词，区分大小写。
 * @param[out] x 成功时写入解析值；失败时保持传入值。
 * @details 接受 true/1/on 和 false/0/off，区分大小写；空值、缺项及其他文本均报错。
 */
void io::SimpleInput::GetBoolean(bastion::Overwatch& tracer, const std::string& kw, bool& x) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    if (!HasKeyword(kw))
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Missing boolean keyword: " + kw, bastion::ERROR);
        return;
    }

    std::string text;
    GetString(tracer, kw, text);
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    if (text == "true" || text == "1" || text == "on")
    {
        x = true;
    }
    else if (text == "false" || text == "0" || text == "off")
    {
        x = false;
    }
    else
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected true/false, 1/0 or on/off for keyword: " + kw, bastion::ERROR);
        return;
    }

    tracer.CheckOut(__func__);
}
