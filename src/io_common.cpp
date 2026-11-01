//
//  io_common.cpp
//  Test_dylib
//  控制输入输出、文件读写的函数
//  Created by 陈酊存 on 2023/8/30.
//

#include "io_common.hpp"
#include <cmath>

//底层函数

void writeDouble(std::ostream& out, double number, const int width, const int precision, bool ENDL)
{
	//输出一个值，默认输出宽度为22，默认输出精度为16，默认不换行
	out << std::setiosflags(std::ios::left | std::ios::showpoint) << std::setw(width) << std::setprecision(precision)
		<< number;
	if (ENDL)
	{
		out << std::endl;
	}
}

void writeInteger(std::ostream& out, int number, const int width, bool ENDL)
{
	// 输出一个整数，默认输出宽度为6，默认不换行
	out << std::setiosflags(std::ios::left) << std::setw(width) << number;
	if (ENDL)
	{
		out << std::endl;
	}
}

void WriteString(std::ostream& out, std::string str, const int width, bool ENDL)
{
	//输出str
	out << std::setiosflags(std::ios::left | std::ios::showpoint) << std::setw(width) << str;
	if (ENDL == true)
	{
		out << std::endl;
	}
}

//int line_count(std::ifstream& infile)
//{
//	std::string line;
//	int lineCnt = 0;  // 统计非空行数
//
//	// 逐行读取文件内容
//	while (std::getline(infile, line))
//	{
//		// 去除行首尾的空白字符，判断是否为空行
//		if (!line.empty())  // 如果该行非空
//		{
//			lineCnt++;
//		}
//	}
//	return lineCnt;
//}

int line_count(std::string file_path)
{
	//统计文件内数据行数（文件末尾如果有空行，则不统计在内）
	std::ifstream infile(file_path.c_str());

	int lineCnt = 0;  // 统计非空行数
	std::string line;
	// 逐行读取文件内容
	while (std::getline(infile, line))
	{
		// 去除行首尾的空白字符，判断是否为空行
		if (!line.empty())  // 如果该行非空
		{
			lineCnt++;
		}
	}
	return lineCnt;
}


// =============================================================================
// C++字符串分割函数 - 与Python行为保持一致
// =============================================================================

/**
 * 按指定分隔符分割字符串（模拟Python str.split(separator)）
 */
void split(const std::string& str, const std::string& delimiter, std::vector<std::string>& tokens)
{
	/*
	 * 功能：使用指定的分隔符将字符串分割成多个子字符串
	 *
	 * 行为特点（与Python str.split(separator) 一致）：
	 * - 如果找不到分隔符，返回包含原字符串的数组
	 * - 连续的分隔符会产生空字符串元素
	 * - 分隔符前后的空字符串也会被保留
	 *
	 * Python对比示例：
	 *   "a,b,c".split(",")     → ['a', 'b', 'c']     正常分割
	 *   "hello".split(",")     → ['hello']           没找到分隔符返回原字符串
	 *   "a,,b".split(",")      → ['a', '', 'b']      连续分隔符产生空字符串
	 *   ",a,".split(",")       → ['', 'a', '']       首尾分隔符保留空字符串
	 *   "".split(",")          → ['']                空字符串
	 */

	tokens.clear();  // 清空输出容器

	// 特殊情况：分隔符为空，返回原字符串（Python会抛出ValueError）
	if (delimiter.empty()) {
		tokens.push_back(str);
		return;
	}

	std::size_t start = 0;  // 当前搜索的起始位置
	std::size_t end = str.find(delimiter);  // 第一个分隔符的位置

	// 循环查找所有分隔符位置
	while (end != std::string::npos) {
		// 提取从start到end之间的子字符串（不包含分隔符）
		// 注意：即使子字符串为空，也要添加到结果中（与Python一致）
		tokens.push_back(str.substr(start, end - start));

		// 移动到下一个分隔符之后的位置
		start = end + delimiter.length();

		// 从新位置开始查找下一个分隔符
		end = str.find(delimiter, start);
	}

	// 添加最后一段（最后一个分隔符之后的部分，或整个字符串如果没有分隔符）
	// 这确保了与Python行为一致：
	// - 如果没有找到分隔符，整个字符串作为唯一元素
	// - 如果最后一个字符是分隔符，会添加一个空字符串
	tokens.push_back(str.substr(start));
}

/**
 * 按空白字符分割字符串（模拟Python str.split()）
 */
void split(const std::string& str, std::vector<std::string>& tokens)
{
	/*
	 * 功能：按任意空白字符（空格、tab、换行等）分割字符串
	 *
	 * 行为特点（与Python str.split() 一致）：
	 * - 自动跳过字符串前后的空白字符
	 * - 将连续的空白字符视为单一分隔符
	 * - 不会产生空字符串元素
	 * - 如果字符串为空或只包含空白字符，返回空数组
	 *
	 * 空白字符包括：空格(' ')、制表符('\t')、换行符('\n')、回车符('\r')、
	 *              换页符('\f')、垂直制表符('\v')
	 *
	 * Python对比示例：
	 *   "hello world".split()      → ['hello', 'world']     正常分割
	 *   "  hello   world  ".split() → ['hello', 'world']    自动去除前后空白
	 *   "hello\t\nworld".split()   → ['hello', 'world']     处理各种空白字符
	 *   "   ".split()              → []                     纯空白返回空数组
	 *   "".split()                 → []                     空字符串返回空数组
	 *   "hello".split()            → ['hello']              单词返回单元素数组
	 */

	tokens.clear();  // 清空输出容器

	std::size_t start = 0;  // 当前处理位置

	// 遍历整个字符串寻找单词
	while (start < str.size()) {
		// 第一阶段：跳过所有空白字符，找到下一个单词的开始位置
		// std::isspace() 识别所有标准空白字符
		while (start < str.size() && std::isspace(static_cast<unsigned char>(str[start]))) {
			start++;
		}

		// 如果已经到达字符串末尾，说明后面都是空白字符
		// 这种情况下不添加任何元素，与Python行为一致
		if (start >= str.size()) {
			break;
		}

		// 第二阶段：从单词开始位置找到单词的结束位置
		std::size_t end = start;
		while (end < str.size() && !std::isspace(static_cast<unsigned char>(str[end]))) {
			end++;  // end指向第一个空白字符的位置或字符串末尾
		}

		// 第三阶段：提取单词并加入结果数组
		// substr(start, end-start) 提取从start开始，长度为(end-start)的子字符串
		// 由于我们跳过了空白字符，这里提取的一定是非空字符串
		tokens.push_back(str.substr(start, end - start));

		// 第四阶段：移动到下一个处理位置
		start = end;  // 从当前单词结束位置开始下一轮搜索
	}
}


std::string strip(const std::string& str)
{
	// 简单的去除字符串前后空白字符的函数，功能类似python的stripe()
	size_t start = 0;
	// 找到第一个非空白字符的位置
	while (start < str.size() && isspace(static_cast<unsigned char>(str[start])))
	{
		start++;
	}
	size_t end = str.size();
	// 找到最后一个非空白字符的位置
	while (end > start && isspace(static_cast<unsigned char>(str[end - 1])))
	{
		end--;
	}
	// 返回去除空白后的字符串
	return str.substr(start, end - start);
}
