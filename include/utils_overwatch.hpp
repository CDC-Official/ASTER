#ifndef CDC_UTILS_OVERWATCH_HPP
#define CDC_UTILS_OVERWATCH_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <filesystem>

// 退出时输出事件还是调用追踪
#define EXIT_FLAG_REPORT_EVENT 1
#define EXIT_FLAG_REPORT_TRACE 2

// 定义控制台输出开关
#define ALERT_MESSAGE false //默认是关闭信息的输出的，改为true则输出信息
#define ALERT_WARNING true //默认是输出警告的
#define ALERT_ERROR true     //默认是输出错误的

//定义强制开关
#define FORCE_ALERT_DEFAULT  0
#define FORCE_ALERT_ON 1
#define FORCE_ALERT_OFF 2

// user-defined events
// alert is instantly prompted to stderr; report is output upon request of report
// message is message; warning is warning; error causes almost immediate stop, no need to worry about repeated msg
#define ALERT_SUPPRESS_REPEATED_USER_MESSAGE      true
#define ALERT_SUPPRESS_REPEATED_USER_WARNING      true
#define REPORT_SUPPRESS_REPEATED_USER_MESSAGE     true
#define REPORT_SUPPRESS_REPEATED_USER_WARNING     true

namespace bastion
{

	// 定义日志事件类型
	enum EventType
	{   // 必须通过名空间bastion访问
		// 默认类型
		UNDEF = 0,
		MESSAGE,
		WARNING,
		ERROR,

		// 自定义message类型
		MESSAGE_USER_TYPE_MIN = 10000,
		MESSAGE_ORBITAL_RESONANCE_DETECTED,
		MESSAGE_ADAMS_STEPSIZE_CHANGE,
		MESSAGE_SPACEWEATHER_BEYOND_SCOPE,
		MESSAGE_USER_TYPE_MAX = 19999,

		// 自定义warning类型
		WARNING_USER_TYPE_MIN = 20000,
		WARNING_MEANPROP_OVER_MILD_ECCENTRICITY,
		WARNING_MEANPROP_OVER_CRITICAL_ECCENTRICITY,
		WARNING_AVERAGING_OVER_MILD_ECCENTRICITY,
		WARNING_AVERAGING_OVER_CRITICAL_ECCENTRICITY,
		WARNING_AVERAGING_OVER_CRITICAL_ECCENTRICITY_JLSHORT,
		WARNING_RKF78_TRUNCATIONERROR_EXCEEDING_LIMIT,
		WARNING_AVERAGING_ZEROINC_MITIGATION_ACTIVATED,
		WARNING_MODEL_NO_EXTERNAL_EOP,
		WARNING_ORBIT_ZERO_ECCENTRICITY,
		WARNING_ORBIT_ZERO_INCLINATION,
		WARNING_ORBIT_ZERO_ECCENTRICITY_AND_INCLINATION,
		WARNING_FRAME_EARTH_GST_FOR_ERA,
		WARNING_FRAME_EARTH_IAU_CARTOGRAPHIC_MODEL,
		WARNING_FRAME_MOON_IAU_CARTOGRAPHIC_MODEL,
		WARNING_UNHANDLED_TESSERAL_ORBITAL_RESONANCE = 21000,
		WARNING_USER_TYPE_MAX = 29999,

		//最大的用户定义数量
		USER_TYPE_MAX = 65535
	};

	enum ErrCode
	{
		// 必须通过名空间bastion访问
		ERR_NOERROR = 0,
		ERR_UNSPECIFIC = 1
	};

	// 事件记录结构体
	struct Event
	{
		std::string file;
		std::string function;
		int line;
		std::string message;
		int type;
		bool alerted = false; // 错误是否已在控制台提示，供失败出口避免重复打印。
	};

// 调用链记录结构体
	struct Call
	{
		std::string file;
		std::string function;
		int line;
	};

// Overwatch 类，用于监控和日志管理
	class Overwatch
	{
	 public:
		Overwatch();

		// 构造函数，初始化时传入文件名、函数名、行号和可选的 ID，默认 ID 为 0
		Overwatch(const std::string& file, const std::string& func, int line, int id = 0);
		~Overwatch();

		void Set(const std::string& file, const std::string& func, int line, int id = 0);

		int Code() const;
		int GetId() const;        // 获取实例的 ID

		// 已记录的警告总数，供逐条资料处理时统计新增警告。
		int WarningCount() const { return nWarning; }

		void CheckIn(const std::string& file, const std::string& func, int line);        // 记录调用链
		void CheckOut(const std::string& func);        // 退出调用链

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
		void LogEvent(const std::string& file, const std::string& func, int line,
			const std::string& msg, int type, int force_alert = FORCE_ALERT_DEFAULT);

		// 打印所有日志事件，重载支持控制台输出和文件输出
		void ReportEvents(std::ostream& os = std::cerr) const;
		void ReportEvents(const std::string& fn_report) const;

		// 打印当前调用链，重载支持控制台输出和文件输出
		void ReportTrace(std::ostream& os = std::cerr) const;
		void ReportTrace(const std::string& fn_report) const;

		/**
		 * @brief 显示尚未提示的错误和当前调用链，不重放普通事件。
		 * @param[in,out] os 接收错误文字和调用位置的输出流，默认使用标准错误流。
		 * @details 仅补充 eventLog 中尚未屏显的 ERROR，随后输出当前调用链。
		 * 写入标准输出或标准错误流后标记该错误已经提示；写入其他流时保留其屏显状态。
		 * 本函数不退出程序，不输出 MESSAGE、WARNING 历史，也不改动错误码。
		 */
		void ReportFailure(std::ostream& os = std::cerr);

		/**
		 * @brief 报告任务失败并以失败状态退出程序。
		 * @param[in,out] os 接收失败信息的输出流，默认使用标准错误流。
		 * @details 标准输出和标准错误流只接收尚未提示的错误及调用链。
		 * 其他输出流保存完整事件和调用链；不自动向控制台重放普通事件。
		 */
		void EmergencyExit(std::ostream& os = std::cerr);

		/**
		 * @brief 在控制台报告失败，并向已启用的程序日志保存诊断后退出。
		 * @param[in,out] log 可选的程序日志文件流；已打开时追加完整事件和调用链并关闭文件。
		 * @details 控制台只补充尚未提示的错误及调用链；日志未打开时仍正常报告失败。
		 * 设置展示和迭代记录由各程序写入 report，不由失败出口重新显示。
		 */
		void EmergencyExit(std::ofstream& log);
		void EmergencyExit(const int exit_flag, std::ostream& os = std::cerr);
		void EmergencyExit(const int exit_flag1, const int exit_flag2, std::ostream& os = std::cerr);

	 private:
		// 线程安全性：如果在多线程环境下使用 Overwatch，需要注意线程安全问题。std::vector 不是线程安全的，如果多个线程共享一个 Overwatch 实例，需要添加互斥锁等机制。
		// ID的使用：在并行计算中，为每个线程或任务分配不同的 ID，有助于在日志中区分不同的执行路径。
		// 错误处理：EmergencyExit() 会直接调用 std::exit(EXIT_FAILURE); 终止程序。在实际使用中，可以根据需要修改为抛出异常或其他错误处理机制。
		int ID = 0;  // 用于区分不同的 Overwatch 实例的ID，默认值为0注意事项：

		int ErrorCode = ErrCode::ERR_NOERROR;
		std::vector<Call> callStack;  // 调用链记录
		std::vector<Event> eventLog;  // 事件日志

		int nEvent = 0;    // 事件数
		int nMessage = 0;    // 信息数
		int nWarning = 0;    // 警告数
		int nError = 0;    // 错误数，应当不超过1

		// Store only event types actually encountered; one tracer is created per candidate.
		std::unordered_map<int, std::size_t> nUserEvent;
	};

} // namespace bastion

#endif // CDC_UTILS_OVERWATCH_HPP
