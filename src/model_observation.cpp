/**
 * @file model_observation.cpp
 * @brief 读取、关联原始观测与测站资料，并计算理论观测值及几何偏导。
 * @details observation 保存文件原始单位；determination 的计算记录使用 rad 和 m。
 * 使用阶段负责坐标和单位转换；理论观测计算接收观测历元的归一化 GCRS 状态。
 *
 * @par 输入输出约定
 * 测站表每行为站号、坐标类型和三个坐标分量；观测行以 UTC MJD、类型、目标号、观测者号开头。
 * 原始资料角度为 deg，距离及三维位置为 m，距离率为 m/s；目标位置类型 7 为 GCRS，类型 8 为 WGS84 地固 XYZ。
 * 理论观测计算接收同历元的归一化 GCRS 轨道及计算观测，返回理论观测量和瞬时几何偏导。
 * 归一化长度尺度为 ae，时间尺度为 UNITTIME=sqrt(ae³/GE)，角度为 rad。
 */
#include "model_observation.hpp"
#include "orbit_geometry.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

using namespace constants::Scale;

namespace observation
{

	/**
	 * @brief 取得观测类型的统一英文名称。
	 * @param[in] type OBS 文件中的类型编号，1—8 对应地基测角、测距、测速、天基测角或目标位置观测。
	 * @return 类型名称的常量字符串；编号不在 1—8 内时返回 Unknown。
	 * @details 供原始资料概览、定轨逐轮统计及最终残差统计共用，不读取或修改观测数据。
	 */
	const char* ObservationTypeName(int type)
	{
		const std::array<const char*, 9> names = {
			"Unknown", "Ground-based RA/Dec", "Ground-based Az/El", "Ground-based range + Az/El",
			"Ground-based range", "Ground-based range rate", "Space-based RA/Dec",
			"GCRS XYZ position", "WGS84 XYZ position"
		};

		if (type < 1 || type > 8) return names[0];

		return names[type];
	}

	/**
	 * @brief 读取地面测站表，保存站号、坐标类型和原始坐标。
	 * @param[in,out] tracer 记录文件打开、列格式、非有限坐标及重复站号错误。
	 * @param[in] file 测站表路径。每行依次为站号、坐标类型和三个坐标分量。
	 * 类型 1 为大地纬度 phi(deg)、经度 lambda(deg)、高度(m)；类型 2 为 ITRS XYZ(m)。
	 * @details 跳过空行和 # 注释，检查列数、坐标有限性及站号唯一性。
	 * 全部读取成功后替换 m_records；失败时保留已有测站资料。读取阶段保存原单位。
	 */
	void Observers::ReadFile(bastion::Overwatch& tracer, const std::string& file)
	{
		if (tracer.Code() != bastion::ERR_NOERROR) return;

		tracer.CheckIn(__FILE__, __func__, __LINE__);

		// 打开本次站表。打不开时报告错误，已有测站资料保持不变。
		std::ifstream input(file);
		if (!input.is_open())
		{
			tracer.LogEvent( __FILE__, __func__, __LINE__, "Cannot open observer file: " + file, bastion::EventType::ERROR);
			return;
		}

		// 先将站表读入临时列表，全部成功后再替换当前资料。
		std::vector<ObserverRecord> records;
		std::string line;
		std::size_t line_number = 0;

		while (std::getline(input, line))
		{
			++line_number;

			// 删除 # 开始的说明，允许整行注释和行末注释。
			const std::size_t comment = line.find('#');
			if (comment != std::string::npos)
			{
				line.erase(comment);
			}

			// 跳过空行或只含空白的行。
			std::istringstream row(line);
			row >> std::ws;
			if (row.eof()) continue;

			ObserverRecord record;
			double coord1 = 0.0;
			double coord2 = 0.0;
			double coord3 = 0.0;
			std::string extra;

			// 站表每行固定五列：站号、坐标类型、三个坐标分量。
			// 暂时按列号接收，后面再按类型存入有明确名称的字段。
			const std::string location = file + ": line " + std::to_string(line_number) + ": ";

			if (!(row >> record.id >> record.coordinate_type
					  >> coord1 >> coord2 >> coord3) ||
				(row >> extra))
			{
				tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Expected: id coordinate_type coord1 coord2 coord3.", bastion::EventType::ERROR);
				return;
			}

			// 这里只检查数值是否有限，不做单位或坐标转换。
			if (!std::isfinite(coord1) || !std::isfinite(coord2) || !std::isfinite(coord3))
			{
				tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Observer coordinates must be finite.", bastion::EventType::ERROR);
				return;
			}

			// 类型 1：大地坐标，纬度、经度、高度，保留 deg、deg、m。
			if (record.coordinate_type == 1)
			{
				record.geodetic_phi_deg = coord1;
				record.geodetic_lambda_deg = coord2;
				record.geodetic_height_m = coord3;
			}
			else if (record.coordinate_type == 2)
			{
				// 类型 2：地固坐标，ITRS 的 X、Y、Z，保留 m。
				record.position_itrs_m[0] = coord1;
				record.position_itrs_m[1] = coord2;
				record.position_itrs_m[2] = coord3;
			}
			else
			{
				tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Observer coordinate type must be 1 or 2.", bastion::EventType::ERROR);
				return;
			}

			// 检查站号唯一性，使每条地基观测能够关联唯一的测站记录。
			for (const ObserverRecord& previous : records)
			{
				if (previous.id == record.id)
				{
					tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Duplicate observer ID: " + record.id, bastion::EventType::ERROR);
					return;
				}
			}

			records.push_back(record);
		}

		// 正常读到文件末尾后才更新资料；底层读取失败则保留旧资料。
		if (input.bad())
		{
			tracer.LogEvent( __FILE__, __func__, __LINE__, "Failed while reading observer file: " + file, bastion::EventType::ERROR);
			return;
		}

		m_records = records;
		tracer.CheckOut(__func__);
	}

	/**
	 * @brief 取得当前保存的原始测站记录。
	 * @return m_records 的只读引用；记录顺序与站表一致，大地坐标为 deg、deg、m，ITRS 坐标为 m。
	 */
	const std::vector<ObserverRecord>& Observers::GetRecords() const
	{
		return m_records;
	}

	/**
	 * @brief 读取类型 1—8 的原始观测资料，保存文件顺序、编号和数值。
	 * @param[in,out] tracer 记录文件打开、未知类型、列格式及非有限数值错误。
	 * @param[in] file 观测文件路径。每行前四列为 UTC MJD、类型、目标号和观测者号，
	 * 后续列按类型解释为角度(deg)、距离(m)、距离率(m/s)或三维位置(m)。
	 * @details 跳过空行和 # 注释，按类型检查列数及数值有限性。
	 * 全部读取成功后替换 m_records；失败时保留已有观测。读取后可另行关联测站，
	 * 资料选用、物理范围检查及坐标和单位转换由使用阶段完成。
	 */
	void Observations::ReadFile(bastion::Overwatch& tracer, const std::string& file)
	{
		if (tracer.Code() != bastion::ERR_NOERROR) return;

		tracer.CheckIn(__FILE__, __func__, __LINE__);

		// 打开观测文件。失败时保留当前对象已有的资料。
		std::ifstream input(file);
		if (!input.is_open())
		{
			tracer.LogEvent( __FILE__, __func__, __LINE__, "Cannot open observation file: " + file, bastion::EventType::ERROR);
			return;
		}

		// 先完整读取到临时列表，全部成功后再替换已有资料。
		std::vector<ObservationRecord> records;
		std::string line;
		std::size_t line_number = 0;

		while (std::getline(input, line))
		{
			++line_number;

			// 删除 # 开始的说明，允许整行注释和行末注释。
			const std::size_t comment = line.find('#');
			if (comment != std::string::npos)
			{
				line.erase(comment);
			}

			// 空行和只含空白的行不产生观测记录。
			std::istringstream row(line);
			row >> std::ws;
			if (row.eof()) continue;

			// 每行重新建立记录，避免上一条观测的字段混入当前记录。
			ObservationRecord record;
			const std::string location = file + ": line " + std::to_string(line_number) + ": ";

			// 所有类型的前四列相同：UTC MJD、类型、目标号、观测者号。
			if (!(row >> record.mjd_utc >> record.type
					  >> record.object_id >> record.observer_id))
			{
				tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Expected: MJD_UTC type object_id observer_id.", bastion::EventType::ERROR);
				return;
			}

			// 按文件定义解释后面的列，只填写原始数值。
			// 这里识别的是文件格式，不判断某个定轨程序能否使用它。
			switch (record.type)
			{
				case 1:
					// 地基赤经、赤纬，deg。
					row >> record.ra_deg >> record.dec_deg;
					break;

				case 2:
					// 地基方位角、高度角，deg。
					row >> record.az_deg >> record.el_deg;
					break;

				case 3:
					// 距离(m)、方位角(deg)、高度角(deg)。
					row >> record.range_m >> record.az_deg >> record.el_deg;
					break;

				case 4:
					// 距离，m。
					row >> record.range_m;
					break;

				case 5:
					// 距离率，m/s。
					row >> record.range_rate_mps;
					break;

				case 6:
					// 天基赤经、赤纬(deg)，随后是平台 GCRS XYZ(m)。
					row >> record.ra_deg >> record.dec_deg
						>> record.observer_position_gcrs_m[0]
						>> record.observer_position_gcrs_m[1]
						>> record.observer_position_gcrs_m[2];
					break;

				case 7:
					// 目标自身的 GCRS XYZ，m。
					row >> record.position_gcrs_m[0]
						>> record.position_gcrs_m[1]
						>> record.position_gcrs_m[2];
					break;

				case 8:
					// 目标自身的 WGS84 地固 XYZ，m。
					row >> record.position_wgs84_ecef_m[0]
						>> record.position_wgs84_ecef_m[1]
						>> record.position_wgs84_ecef_m[2];
					break;

				default:
					// 未定义的编号无法确定各列含义，作为文件格式错误报告。
					tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Undefined observation type: " + std::to_string(record.type), bastion::EventType::ERROR);
					return;
			}

			// 检查本类型所需的列是否读全，且没有多余字段。
			std::string extra;
			if (!row || (row >> extra))
			{
				tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Invalid values or column count for observation type " + std::to_string(record.type), bastion::EventType::ERROR);
				return;
			}

			// 未使用的字段默认是零，因此可以统一检查各数值是否有限。
			// 赤纬、距离等物理取值范围留到读取后的使用检查阶段。
			bool values_finite =
				std::isfinite(record.mjd_utc) &&
				std::isfinite(record.ra_deg) &&
				std::isfinite(record.dec_deg) &&
				std::isfinite(record.az_deg) &&
				std::isfinite(record.el_deg) &&
				std::isfinite(record.range_m) &&
				std::isfinite(record.range_rate_mps);

			for (int i = 0; i < 3; ++i)
			{
				if (!std::isfinite(record.position_gcrs_m[i]) ||
					!std::isfinite(record.position_wgs84_ecef_m[i]) ||
					!std::isfinite(record.observer_position_gcrs_m[i]))
				{
					values_finite = false;
				}
			}

			if (!values_finite)
			{
				tracer.LogEvent( __FILE__, __func__, __LINE__, location + "Observation values must be finite.", bastion::EventType::ERROR);
				return;
			}

			// 保留原单位和文件顺序，此时尚未关联地面测站。
			records.push_back(record);
		}

		// 底层读取失败不能当作正常文件结束。
		if (input.bad())
		{
			tracer.LogEvent( __FILE__, __func__, __LINE__, "Failed while reading observation file: " + file, bastion::EventType::ERROR);
			return;
		}

		m_records = records;
		tracer.CheckOut(__func__);
	}

	/**
	 * @brief 向资料末尾添加一条已经准备好的原始观测。
	 * @param[in] record 待复制的观测记录，保留其中的 UTC MJD、类型、编号、原单位数值及测站关联状态。
	 * @details 将 record 的副本追加到 m_records；本函数只保存记录，不执行格式或物理范围检查。
	 */
	void Observations::AddRecord(const ObservationRecord& record)
	{
		m_records.push_back(record);
	}

	/**
	 * @brief 将当前观测按类型 1—8 的原始资料格式写入文件。
	 * @param[in,out] tracer 记录编号、类型、非有限数值以及文件打开和写入错误。
	 * @param[in] file 输出文件路径；打开成功后覆盖同名文件。
	 * @details 按 m_records 顺序写出 UTC MJD、类型、目标号、观测者号及对应类型的数据列。
	 * 角度为 deg，距离和位置为 m，距离率为 m/s；数值采用 17 位有效数字，字段以空格分隔。
	 * 写出前检查编号、类型及数值有限性，观测内容和测站关联状态保持不变。
	 */
	void Observations::WriteFile(bastion::Overwatch& tracer, const std::string& file) const
	{
		if (tracer.Code() != bastion::ERR_NOERROR) return;

		tracer.CheckIn(__FILE__, __func__, __LINE__);

		// 打开输出文件前，先检查各条记录的编号、类型和数值。
		for (const ObservationRecord& record : m_records)
		{
			if (record.object_id.empty() || record.observer_id.empty() ||
				record.object_id.find_first_of(" \t\r\n#") != std::string::npos ||
				record.observer_id.find_first_of(" \t\r\n#") != std::string::npos ||
				record.type < 1 || record.type > 8)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot write observation: require two non-empty IDs without whitespace/# and type 1..8.", bastion::EventType::ERROR);
				return;
			}

			// 检查记录内的各个数值，物理范围由使用这些资料的程序判断。
			bool values_finite = std::isfinite(record.mjd_utc) &&
								 std::isfinite(record.ra_deg) && std::isfinite(record.dec_deg) &&
								 std::isfinite(record.az_deg) && std::isfinite(record.el_deg) &&
								 std::isfinite(record.range_m) && std::isfinite(record.range_rate_mps);

			for (int i = 0; i < 3; ++i)
			{
				values_finite = values_finite &&
								std::isfinite(record.position_gcrs_m[i]) &&
								std::isfinite(record.position_wgs84_ecef_m[i]) &&
								std::isfinite(record.observer_position_gcrs_m[i]);
			}

			if (!values_finite)
			{
				tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot write observation: values must be finite.", bastion::EventType::ERROR);
				return;
			}
		}

		std::ofstream output(file);
		if (!output.is_open())
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot open observation file to write: " + file, bastion::EventType::ERROR);
			return;
		}

		// 17 位有效数字供 double 读写往返；用空格分列，编号较长也不会和下一列粘连。
		output << std::setprecision(17);

		for (const ObservationRecord& record : m_records)
		{
			output << record.mjd_utc << ' ' << record.type << ' '
				   << record.object_id << ' ' << record.observer_id;

			switch (record.type)
			{
				case 1: // 地基赤经、赤纬，deg。
					output << ' ' << record.ra_deg << ' ' << record.dec_deg;
					break;
				case 2: // 地基方位、高度角，deg。
					output << ' ' << record.az_deg << ' ' << record.el_deg;
					break;
				case 3: // 距离 m，随后为方位、高度角 deg。
					output << ' ' << record.range_m << ' ' << record.az_deg << ' ' << record.el_deg;
					break;
				case 4: // 距离，m。
					output << ' ' << record.range_m;
					break;
				case 5: // 距离率，m/s。
					output << ' ' << record.range_rate_mps;
					break;
				case 6: // 天基赤经赤纬 deg，随后为平台 GCRS XYZ，m。
					output << ' ' << record.ra_deg << ' ' << record.dec_deg;

					for (int i = 0; i < 3; ++i) output << ' ' << record.observer_position_gcrs_m[i];
					break;
				case 7: // 目标 GCRS XYZ，m。
					for (int i = 0; i < 3; ++i) output << ' ' << record.position_gcrs_m[i];
					break;
				case 8: // 目标 WGS84 地固 XYZ，m。
					for (int i = 0; i < 3; ++i) output << ' ' << record.position_wgs84_ecef_m[i];
					break;
			}

			output << '\n';
		}

		output.close();
		if (!output)
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing observation file: " + file, bastion::EventType::ERROR);
			return;
		}

		tracer.CheckOut(__func__);
	}

	/**
	 * @brief 取得当前保存的原始观测记录。
	 * @return m_records 的只读引用，保留资料顺序、原单位数值及当前测站关联状态。
	 */
	const std::vector<ObservationRecord>& Observations::GetRecords() const
	{
		return m_records;
	}

	/**
	 * @brief 按观测者编号，为地基观测关联原始测站记录。
	 * @param[in,out] tracer 读取调用链的已有错误状态，并记录本次测站关联的调用过程。
	 * @param[in] observers 地面测站集合；按站号读取记录，并复制其坐标类型和原单位坐标。
	 * @details 逐条清除 m_records 中原有的 ground_observer 和 is_observer_attached，再重新关联。
	 * 类型 1—5 按站号查找并复制测站，未匹配时标记为 false；类型 6 标记为 true，使用记录自带的平台位置。
	 * 类型 7、8 保持 false。匹配状态供使用阶段检查，本函数不转换坐标或单位。
	 */
	void Observations::AttachObservers(bastion::Overwatch& tracer, const Observers& observers)
	{
		if (tracer.Code() != bastion::ERR_NOERROR) return;

		tracer.CheckIn(__FILE__, __func__, __LINE__);

		const std::vector<ObserverRecord>& observer_records = observers.GetRecords();

		// 范围 for：依次访问 m_records 中的每条观测，不需要手动维护下标。
		// & 表示引用原记录，不产生副本；下面的修改会直接保存到 m_records 中。
		for (ObservationRecord& record : m_records)
		{
			// 清除上次关联的站址和状态，本次重新匹配。
			record.ground_observer = ObserverRecord{};
			record.is_observer_attached = false;

			// 类型 1—5 是地基观测，需要从站表取得测站资料。
			if (record.type >= 1 && record.type <= 5)
			{
				// const& 表示只读地引用原测站记录，不复制，也不允许修改它。
				for (const ObserverRecord& observer : observer_records)
				{
					if (observer.id == record.observer_id)
					{
						// 复制匹配的原始测站记录，保留其坐标表达和原单位。
						// 这里保存的是站表的副本；站表更新后需再次调用本函数。
						record.ground_observer = observer;
						record.is_observer_attached = true;
						break;
					}
				}
			}
			else if (record.type == 6)
			{
				// ReadFile 已完整读入本条的平台 GCRS 位置，不查地面站表。
				// 此标志只表示位置已提供，不表示已通过物理合理性检查。
				record.is_observer_attached = true;
			}

			// 类型 7、8 是目标自身的位置观测，不需要关联外部测站。
			// 地基观测未找到站号时保留 false，由主程序检查实际选用的记录。
		}

		tracer.CheckOut(__func__);
	}

	/**
	 * @brief 输出已保存观测的总条数、各类型条数及历元范围。
	 * @param[in,out] output 接收文字概览的输出流；本函数临时调整名称对齐和历元精度，结束时恢复流的格式。
	 * @details 统计整份 m_records，类型编号应在 1—8 内；遍历全部记录确定最早、最晚 UTC MJD。
	 * 历元按 12 位小数显示，资料为空时输出提示。
	 */
	void Observations::PrintSummary(std::ostream& output) const
	{
		output << "\n----------------------------- Observation data -----------------------------\n";
		output << "Loaded records     : " << m_records.size() << '\n';

		if (m_records.empty())
		{
			output << "No observation records.\n";
			return;
		}

		// ReadFile 只接受类型 1—8，因此可以直接用类型编号作为数组下标。
		// 下标 0 不使用，数组中保存每种类型的观测条数。
		std::array<std::size_t, 9> type_counts = {};
		double min_mjd_utc = m_records[0].mjd_utc;
		double max_mjd_utc = m_records[0].mjd_utc;

		// 逐条读取原记录；const& 表示不复制记录，也不允许在循环中修改它。
		for (const ObservationRecord& record : m_records)
		{
			++type_counts[record.type];

			// 文件不一定按时间排序，因此比较全部历元，而不只取首末两条。
			if (record.mjd_utc < min_mjd_utc)
			{
				min_mjd_utc = record.mjd_utc;
			}

			if (record.mjd_utc > max_mjd_utc)
			{
				max_mjd_utc = record.mjd_utc;
			}
		}

		// 暂时调整名称对齐和历元格式；结束时恢复调用者的输出格式和精度。
		const std::ios_base::fmtflags old_flags = output.flags();
		const std::streamsize old_precision = output.precision();
		const char old_fill = output.fill();
		output.fill(' ');

		// 只列出文件中实际出现的类型，统计覆盖整份已读入的观测资料。
		output << "Observation types  :\n";

		for (int type = 1; type <= 8; ++type)
		{
			if (type_counts[type] == 0) continue;

			// setw 只控制紧接着的名称字段；left 让名称左对齐，使后面的冒号对齐。
			output << "  Type " << type << " - " << std::left << std::setw(28)
				   << ObservationTypeName(type) << ": " << type_counts[type] << " records\n";
		}

		// 最早、最晚观测历元均为 UTC MJD，显示时保留 12 位小数。
		output << "\nEpochs (MJD UTC)\n"
			   << "  Earliest observation : " << std::fixed << std::setprecision(12)
			   << min_mjd_utc << '\n'
			   << "  Latest observation   : " << max_mjd_utc << '\n';

		output.flags(old_flags);
		output.precision(old_precision);
		output.fill(old_fill);
	}
} // namespace observation

namespace determination
{
	/**
	 * @brief 计算一条观测的理论观测量和瞬时几何偏导。
	 * @param[in,out] tracer 记录观测类型或几何定义域错误，调用方通过 Code() 检查。
	 * @param[in] obs 计算观测记录。类型 1/6 读取观测者 GCRS 位置(m)及观测 Dec(rad)，
	 * 类型 7 只使用观测类型；sigma 由后续加权阶段使用。
	 * @param[in] state 观测历元的归一化 GCRS Cartesian 状态，分量为 [x,y,z,vx,vy,vz]，
	 * 长度尺度为 ae，时间尺度为 UNITTIME。
	 * @param[out] calc_obs 理论观测量。类型 1/6 为 [RA,Dec,0](rad)，类型 7 为目标 GCRS [X,Y,Z](m)。
	 * @param[out] partials 3×6 的瞬时几何偏导，列对应归一化 Cartesian 六分量，速度列为零。
	 * 测角前两行分别为 cos(观测Dec)*dRA/dx、dDec/dx，供东西、南北方向残差使用；
	 * 位置观测的三行为 d[X,Y,Z]/dx，单位为 m。
	 * @details 采用同历元几何模型，计算理论值及其几何偏导；O-C 由批处理残差循环计算。
	 * 测角第一行偏导含观测 Dec 的余弦因子，calc_obs 中仍保存 RA 本身。
	 * 模型不含光行时、光行差和折射，零距离或极点附近的测角几何通过 tracer 报错。
	 */
	void Compute_ObsAndPartials(bastion::Overwatch& tracer, const Observation& obs, const orbit::Cartesian& state,
		Eigen::Vector3d& calc_obs, Eigen::Matrix<double, 3, 6>& partials)
	{
		if (tracer.Code() != bastion::ERR_NOERROR) return;

		tracer.CheckIn(__FILE__, __func__, __LINE__);

		// 将归一化位置恢复为 GCRS 米制位置，与观测者位置使用相同单位。
		std::array<double, 6> state_values;
		state.ToArray(state_values.data());
		Eigen::Vector3d position_gcrs_m(state_values[0] * UNIT_LENGTH_M, state_values[1] * UNIT_LENGTH_M, state_values[2] * UNIT_LENGTH_M);

		// calc_obs 和偏导矩阵按三分量准备，测角只使用前两行；瞬时观测的速度偏导为零。
		calc_obs.setZero();
		partials.setZero();

		if (obs.type == ObservationType::InertialPosition)
		{
			// 位置观测直接取目标 XYZ，ae 将偏导恢复为对归一化位置的米制偏导。
			calc_obs = position_gcrs_m;
			partials.block<3, 3>(0, 0) = UNIT_LENGTH_M * Eigen::Matrix3d::Identity();
			tracer.CheckOut(__func__);
			return;
		}

		if (!(obs.type == ObservationType::GroundRaDec || obs.type == ObservationType::SpaceRaDec))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "Unsupported observation type.", bastion::EventType::ERROR);
			return;
		}

		// 目标位置减去观测者位置，得到 GCRS 下的相对位置向量，单位 m。
		const Eigen::Vector3d dr_gcrs_m = position_gcrs_m - obs.observer_position_gcrs_m;
		// 计算得到赤道坐标，长度单位为m（同输入参数的长度单位），角度单位为弧度，范围为 RA [0,2π)，Dec [-π/2,π/2]。
		const EquatorialCoordinates calc_equatorial = LOS_to_Equatorial(dr_gcrs_m);

		// rho 为斜距，xy 为相对位置在赤道面内的长度，均为 m。
		const double rho = calc_equatorial.rho;
		const double xy = std::hypot(dr_gcrs_m[0], dr_gcrs_m[1]);

		if (!(rho > 0 && xy > rho * 1e-10 && std::abs(std::cos(obs.value[1])) > 1e-10))
		{
			tracer.LogEvent(__FILE__, __func__, __LINE__, "RA/Dec geometry too close to pole or zero range.", bastion::EventType::ERROR);
			return;
		}

		// 保存理论赤经、赤纬，供残差循环与实测值作差。
		calc_obs[0] = calc_equatorial.right_ascension;
		calc_obs[1] = calc_equatorial.declination;

		// 第一行偏导乘 cos(观测Dec)，对应东西方向；第二行对应南北方向。
		// ae 将对米制位置的偏导转为对归一化位置的偏导。
		const double cos_dec = std::cos(obs.value[1]);
		partials.block<1, 3>(0, 0) = cos_dec * UNIT_LENGTH_M * Eigen::RowVector3d(-dr_gcrs_m[1] / (xy * xy), dr_gcrs_m[0] / (xy * xy), 0);
		partials.block<1, 3>(1, 0) = UNIT_LENGTH_M * Eigen::RowVector3d(-dr_gcrs_m[0] * dr_gcrs_m[2] / xy, -dr_gcrs_m[1] * dr_gcrs_m[2] / xy, xy) / (rho * rho);

		tracer.CheckOut(__func__);
	}

} // namespace determination
