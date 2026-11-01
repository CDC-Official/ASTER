#include "utils_overwatch.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdlib>  // for std::exit

namespace bastion
{

	Overwatch::Overwatch()
	{
		// 默认初始化，如果这样的话需要手动再Set
		//do nothing
	}
	// 构造函数，接受文件名、函数名、行号，以及可选的 ID
	Overwatch::Overwatch(const std::string& file, const std::string& func, int line, int id)
	{
		// 初始化的时候直接CheckIn，记录创建对象时的调用，id默认参数为0
		ID = id;
		CheckIn(file, func, line);
	}

	Overwatch::~Overwatch()
	{
		// 自动清理调用链和事件日志，无需手动删除
	}

	void Overwatch::Set(const std::string& file, const std::string& func, int line, int id)
	{
		//Set功能完全类似于手动赋值的初始化，再无参数初始化后调用，id默认参数为0
		ID = id;
		CheckIn(file, func, line);
	}

	int Overwatch::Code() const
	{
		return this->ErrorCode;
	}
	int Overwatch::GetId() const
	{
		return ID;
	}

	// 记录调用链
	void Overwatch::CheckIn(const std::string& file, const std::string& func, int line)
	{
		// C++11引入的列表初始化（list initialization）功能，允许使用花括号{}的形式为类或结构体的成员赋值并直接创建对象。
		callStack.push_back({ file, func, line });
	}
	// 退出调用链
	void Overwatch::CheckOut(const std::string& func)
	{
		// 检查调用栈是否为空
		if (callStack.empty())
		{
			// 打印错误信息到标准错误输出
			std::cerr << "Emergency Exit: Call stack is empty during CheckOut!" << std::endl;
			// 直接紧急退出
			exit(EXIT_FAILURE);
		}

		// 获取调用栈顶的函数信息
		const auto& topCall = callStack.back();

		// 验证当前传入的函数名是否与调用栈顶记录的函数名匹配
		if (topCall.function != func)
		{
			// 打印错误信息
			std::cerr << "Emergency Exit: Function mismatch in CheckOut! Expected: "
					  << topCall.function << ", Got: " << func << std::endl;
			// 直接紧急退出
			exit(EXIT_FAILURE);
		}
		// 没问题的话移除调用栈栈顶函数，不用移除事件
		callStack.pop_back();
	}

	/**
	 * @brief 保存一条运行事件，并按事件类型和屏显开关即时提示。
	 * @param[in] file 产生事件的源文件路径。
	 * @param[in] func 产生事件的函数名。
	 * @param[in] line 产生事件的源文件行号。
	 * @param[in] msg 事件说明，包含定位问题所需的参数、记录号或文件路径。
	 * @param[in] type MESSAGE、WARNING、ERROR 或对应的用户事件编号。
	 * @param[in] force_alert FORCE_ALERT_DEFAULT 遵守默认开关，ON 强制提示，OFF 仅保存事件。
	 * @details 所有事件均留在 eventLog；ERROR 设置错误码并记录是否已在控制台提示。
	 * 失败出口据此补充尚未提示的错误。设置展示和逐轮计算文字由各程序直接写入控制台及 report。
	 */
	void Overwatch::LogEvent(const std::string& file, const std::string& func, int line,
		const std::string& msg, int type, int force_alert)
	{
		// 先push_back一个event，先不分类type
		Event newEvent = { file, func, line, msg, type };
		eventLog.push_back(newEvent);

		//细分事件的类型，以及是否需要警告、是否用户定义
		int baseType = EventType::UNDEF;//注意，此处没有对UNDEF的事件特别处理！谨慎生成这些事件！
		bool noAlert = false;
		bool isUserDefined = false;

		if (type == EventType::MESSAGE || type == EventType::WARNING || type == EventType::ERROR)
		{
			// 三种默认类型之一
			baseType = type;
		}
		else if (type >= EventType::MESSAGE_USER_TYPE_MIN && type <= EventType::MESSAGE_USER_TYPE_MAX)
		{
			// 用户定义的message类型
			baseType = EventType::MESSAGE;
			isUserDefined = true;
			nUserEvent[type]++;
			noAlert = (ALERT_SUPPRESS_REPEATED_USER_MESSAGE) ? (nUserEvent[type] > 1) : false;
		}
		else if (type >= EventType::WARNING_USER_TYPE_MIN && type <= EventType::WARNING_USER_TYPE_MAX)
		{
			baseType = EventType::WARNING;
			isUserDefined = true;
			nUserEvent[type]++;
			noAlert = (ALERT_SUPPRESS_REPEATED_USER_WARNING) ? (nUserEvent[type] > 1) : false;
		}

		// 记录事件个数
		nEvent++;
		switch (baseType)
		{
			case EventType::MESSAGE:
			{
				nMessage++;
				break;
			}
			case EventType::WARNING:
			{
				nWarning++;
				break;
			}
			case EventType::ERROR:
			{
				ErrorCode = ERR_UNSPECIFIC;  // 未指定的错误
				nError++;
				break;
			}
			default:
				break;
		}

		//是否输出
		if (force_alert != FORCE_ALERT_OFF && (!noAlert))
		{
			//如果没有强制关掉输出的话
			switch (baseType)
			{
				case EventType::MESSAGE:
				{
					if ((force_alert == FORCE_ALERT_ON) || ALERT_MESSAGE)  // 条件判断
					{
						// 强制打开输出且需要输出信息，缺一不可
						std::cerr << std::setw(5) << std::setfill('0') << this->ID << "> " << msg;
						if (isUserDefined && ALERT_SUPPRESS_REPEATED_USER_MESSAGE)
						{
							std::cerr << " (further repeated warnings will be suppressed)";
						}
						std::cerr << std::endl;  // 换行
					}
					break;
				}
				case EventType::WARNING:
				{
					if ((force_alert == FORCE_ALERT_ON) || ALERT_WARNING)  // 条件判断
					{
						// 使用 C++17 的 std::filesystem 处理文件路径
						std::filesystem::path filePath(file);
						std::string f = filePath.filename().string();  // 获取文件名，包括扩展名

						// 输出警告信息
						std::cerr << std::setw(5) << std::setfill('0') << this->ID
								  << "> [WARNING] " << msg
								  << " (" << func << " in " << f << " @ #" << line << ")";

						if (isUserDefined && ALERT_SUPPRESS_REPEATED_USER_WARNING)
						{
							std::cerr << " (further repeated warnings will be suppressed)";
						}
						std::cerr << std::endl;  // 换行
					}

					break;
				}
				case EventType::ERROR:
				{
					ErrorCode = ERR_UNSPECIFIC;  // 未指定的错误，之前已经赋值过，这里其实可有可无
					if ((force_alert == FORCE_ALERT_ON) || ALERT_ERROR)  // 条件判断
					{
						// 使用 C++17 的 std::filesystem 处理文件路径
						std::filesystem::path filePath(file);
						std::string f = filePath.filename().string();  // 获取文件名，包括扩展名
						// 输出错误信息
						std::cerr << std::setw(5) << std::setfill('0') << this->ID
								  << "> [ERROR] " << msg
								  << " (" << func << " in " << f << " @ #" << line << ")\n";
						eventLog.back().alerted = true;
					}
					break;
				}
				default:
					break;
			}
		}
	}

	void Overwatch::ReportEvents(std::ostream& os) const
	{
		os << "------------------------ Event Log ------------------------" << std::endl;

		// 遍历所有事件
		for (size_t i = 0; i < eventLog.size(); ++i)
		{
			const auto& event = eventLog[i];
			std::string eventType;

			// 根据事件类型设置前缀和标签
			switch (event.type)
			{
				case EventType::MESSAGE:
					eventType = "[MESSAGE]";
					break;
				case EventType::WARNING:
					eventType = "[WARNING]";
					break;
				case EventType::ERROR:
					eventType = "[ERROR]  ";
					break;
			}

            if (event.type >= MESSAGE_USER_TYPE_MIN && event.type <= MESSAGE_USER_TYPE_MAX)
                eventType = "[MESSAGE]";
            else if (event.type >= WARNING_USER_TYPE_MIN && event.type <= WARNING_USER_TYPE_MAX)
                eventType = "[WARNING]";
            else if (eventType.empty()) eventType = "[UNKNOWN]";

            const auto fileBasename = std::filesystem::path(event.file).filename().string();
            // Avoid truncating long paths/messages, and leave the caller's stream formatting unchanged.
            std::ostringstream line;
            line << '#' << std::setw(5) << std::setfill('0') << i + 1 << ' '
                 << eventType << ' ' << event.message << " (" << event.function
                 << " in " << fileBasename << " @ #" << event.line << ')';
            os << line.str() << std::endl;
		}
	}
	// 打印所有日志事件到指定文件
	void Overwatch::ReportEvents(const std::string& fn_report) const
	{
		std::ofstream file(fn_report);
		if (file.is_open())
		{
			ReportEvents(file);  // 使用文件流输出事件日志
			file.close();
		}
		else
		{
			std::cerr << "Could not open report file: " << fn_report << std::endl;
		}
	}

	// 打印当前调用链到控制台
	void Overwatch::ReportTrace(std::ostream& os) const
	{
		os << "\"------------------------ Call Stack ------------------------ (ID: " << ID << ")" << std::endl;  // 打印 ID
		for (const auto& call : callStack)
		{
			os << "In function: " << call.function
			   << " (File: " << call.file << ", Line: " << call.line << ")" << std::endl;
		}
	}
	// 打印当前调用链到指定文件
	void Overwatch::ReportTrace(const std::string& fn_report) const
	{
		std::ofstream file(fn_report);
		if (file.is_open())
		{
			ReportTrace(file);  // 使用文件流输出调用链
			file.close();
		}
		else
		{
			std::cerr << "Error: Could not open file: " << fn_report << std::endl;
		}
	}

	/**
	 * @brief 显示尚未提示的错误和当前调用链，不重放普通事件。
	 * @param[in,out] os 接收错误文字和调用位置的输出流，默认使用标准错误流。
	 * @details 仅补充 eventLog 中尚未屏显的 ERROR，随后输出当前调用链。
	 * 写入标准输出或标准错误流后标记该错误已经提示；写入其他流时保留其屏显状态。
	 * 本函数不退出程序，不输出 MESSAGE、WARNING 历史，也不改动错误码。
	 */
	void Overwatch::ReportFailure(std::ostream& os)
	{
		// 子函数可以暂时关闭错误提示，失败出口补上这些错误，已显示的错误不再重复。
		for (Event& event : eventLog)
		{
			if (event.type != EventType::ERROR || event.alerted) continue;

			std::ostringstream message;
			message << std::setfill('0') << std::setw(5) << ID << "> [ERROR] " << event.message
					<< " (" << event.function << " in " << std::filesystem::path(event.file).filename().string()
					<< " @ #" << event.line << ")\n";
			os << message.str();

			if (&os == &std::cerr || &os == &std::cout) event.alerted = true;
		}

		ReportTrace(os);
	}

	/**
	 * @brief 在控制台报告失败，并向已启用的程序日志保存诊断后退出。
	 * @param[in,out] log 可选的程序日志文件流；已打开时追加完整事件和调用链并关闭文件。
	 * @details 控制台只补充尚未提示的错误及调用链；日志未打开时仍正常报告失败。
	 * 设置展示和迭代记录由各程序写入 report，不由失败出口重新显示。
	 */
	void Overwatch::EmergencyExit(std::ofstream& log)
	{
		ReportFailure();

		if (log.is_open())
		{
			ReportEvents(log);
			ReportTrace(log);
			log.close();
		}

		std::exit(EXIT_FAILURE);
	}

	/**
	 * @brief 报告任务失败并以失败状态退出程序。
	 * @param[in,out] os 接收失败信息的输出流，默认使用标准错误流。
	 * @details 标准输出和标准错误流只接收尚未提示的错误及调用链。
	 * 其他输出流保存完整事件和调用链；不自动向控制台重放普通事件。
	 */
	void Overwatch::EmergencyExit(std::ostream& os)
	{
		if (&os == &std::cerr || &os == &std::cout)
		{
			ReportFailure(os);
		}
		else
		{
			ReportEvents(os);
			ReportTrace(os);
		}

		std::exit(EXIT_FAILURE);
	}
	void Overwatch::EmergencyExit(const int exit_flag, std::ostream& os)
	{
		// 如果是需要输出事件
		if (exit_flag == EXIT_FLAG_REPORT_EVENT)
		{
			ReportEvents(os);
		}
			// 如果是需要输出调用追踪
		else if (exit_flag == EXIT_FLAG_REPORT_TRACE)
		{
			ReportTrace(os);
		}

		// 检查 os 是否为文件流，如果是并关闭文件
		if (&os != &std::cerr && &os != &std::cout)  // 确保 os 不是标准输出或标准错误流
		{
			if (std::ofstream* fileStream = dynamic_cast<std::ofstream*>(&os))
			{
				if (fileStream->is_open())  // 如果文件流已打开
				{
					fileStream->close();    // 关闭文件
				}
			}
		}

		// 退出程序
		std::exit(EXIT_FAILURE);

	}
	void Overwatch::EmergencyExit(const int exit_flag1, const int exit_flag2, std::ostream& os)
	{
		if (exit_flag1 == exit_flag2)
		{
			EmergencyExit(exit_flag1, os);
		}
		else
		{
			int flags[2] = { exit_flag1, exit_flag2 };
			for (int j = 0; j < 2; j++)
			{
				switch (flags[j])
				{
					case EXIT_FLAG_REPORT_EVENT:
						ReportEvents(os);
						break;
					case EXIT_FLAG_REPORT_TRACE:
						ReportTrace(os);
						break;
				}
			}
		}
		// 检查 os 是否为文件流，如果是并关闭文件
		if (&os != &std::cerr && &os != &std::cout)  // 确保 os 不是标准输出或标准错误流
		{
			if (std::ofstream* fileStream = dynamic_cast<std::ofstream*>(&os))
			{
				if (fileStream->is_open())  // 如果文件流已打开
				{
					fileStream->close();    // 关闭文件
				}
			}
		}
		// 退出程序
		std::exit(EXIT_FAILURE);
	}

}// namespace bastion
