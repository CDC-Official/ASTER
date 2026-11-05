/**
 * @file frame_definition.hpp
 * @brief 地面站址的大地坐标、ITRS 和 GCRS 位置对象及六个方向的转换。
 * @details 大地坐标采用 WGS84。对象未归一化时高度、XYZ 为 m，角度为 deg；
 * 椭球计算使用按 constants::Scale::UNIT_LENGTH_M 归一化的高度、XYZ 和 rad 制角度。
 * 地固与天球位置转换接收 UTC MJD 和同历元 EOP，保持输入的长度单位。
 */
#ifndef CDC_FRAME_DEFINITION_HPP
#define CDC_FRAME_DEFINITION_HPP

#include "earth_eop.hpp"
#include "utils_overwatch.hpp"

using namespace bastion;

// 地面站址的位置表示：大地坐标、ITRS 地固直角坐标和 GCRS 地心天球直角坐标。
// 这些类不存储速度；直角坐标可以是米制或归一化值，由 m_is_normalized 区分。
// ITRS 类型不指定 ITRF 实现版本；测站数据实际采用的版本应在数据来源中说明。

namespace frame
{
	//只讨论地球坐标系
	class Geodetic;
	class GCRS;
	class ITRS;

	// 大地坐标采用 WGS84 参考椭球。
	// 分量顺序为 height、lambda（经度）、phi（大地纬度）。
	// 未归一化时为 m、deg、deg；归一化后为 height/UNIT_LENGTH_M、rad、rad。
	// phi 由椭球法线定义，不是地心纬度或天文纬度。
	class Geodetic
	{
	 private:
		double m_coord[3];// height、lambda、phi，单位由归一化标志决定。
		bool m_is_normalized;

	 public:
		Geodetic();
		// 构造顺序保持高度、经度、纬度，与站表的纬度、经度、高度顺序不同。
		Geodetic(double height, double lambda, double phi, bool is_normalized);
		//析构函数
		~Geodetic();
		//修改设置成员值
		void Set(double height, double lambda, double phi, bool is_normalized);
		//转为一般的数组
		void ToArray(double* geodetic_array) const;
		//获取成员值
		double GetValue(Overwatch& tracer, int idx) const;
		// 高度除以 uLength，角度由 deg 转为 rad；进行椭球转换时 uLength 应取 constants::Scale::UNIT_LENGTH_M。
		void Normalize(Overwatch& tracer, double uLength);
		// 高度乘回 uLength，角度由 rad 恢复为 deg。
		void Restore(Overwatch& tracer, double uLength);

		/**
		 * @brief 将 WGS84 大地经纬高转换为 ITRS 地固位置。
		 * @param[in,out] tracer 记录本次转换的调用链和输入错误。
		 * @param[out] itrs_position 地固位置 [X,Y,Z]，以 constants::Scale::UNIT_LENGTH_M 归一化。
		 * @details 当前对象须已归一化：高度除以 UNIT_LENGTH_M，经度和大地纬度为 rad。
		 * 使用 constants::Wgs84 的赤道半径与扁率；失败时保留传入的输出对象。
		 */
		void ToITRS(Overwatch& tracer, ITRS& itrs_position) const;

		/**
		 * @brief 将 WGS84 大地经纬高转换为指定时刻的 GCRS 地心天球位置。
		 * @param[in,out] tracer 记录椭球转换和地球定向旋转的调用链及错误。
		 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
		 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
		 * @param[out] gcrs_position GCRS 位置 [X,Y,Z]，以 constants::Scale::UNIT_LENGTH_M 归一化。
		 * @details 当前对象须已归一化：高度除以 UNIT_LENGTH_M，角度为 rad。
		 * 依次调用 ToITRS 和 ITRS::ToGCRS；采用与初轨程序相同的 IAU 2000B 地球定向旋转。
		 * 失败时保留传入的输出对象。
		 */
		void ToGCRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, GCRS& gcrs_position) const;

	};

	class ITRS
	{
	 private:
		double m_coord[3];//X Y Z 暂不考虑速度
		bool m_is_normalized;

	 public:
		ITRS();
		//自定义的构造函数 X Y Z
		ITRS(double X, double Y, double Z, bool is_normalized);//参数中删除了时间
		//析构函数
		~ITRS();

		//修改设置成员值
		void Set(double X, double Y, double Z, bool is_normalized);
		//转为一般的数组
		void ToArray(double* itrs_array) const;
		//获取成员值
		double GetValue(Overwatch& tracer, int idx) const;
		//归一化
		void Normalize(Overwatch& tracer, double uLength);
		//去归一化
		void Restore(Overwatch& tracer, double uLength);

		/**
		 * @brief 将 ITRS 地固位置转换为 WGS84 大地经纬高。
		 * @param[in,out] tracer 记录本次反算的调用链和输入错误。
		 * @param[out] geodetic_position 大地坐标 [高度,经度,大地纬度]；高度以
		 * constants::Scale::UNIT_LENGTH_M 归一化，角度为 rad。
		 * @details 当前对象须为以 UNIT_LENGTH_M 归一化的有限、非零地固位置。
		 * 使用 SOFA iauGc2gde 和 WGS84 椭球参数；自转轴上经度取 0，失败时保留输出。
		 */
		void ToGeodetic(Overwatch& tracer, Geodetic& geodetic_position) const;

		/**
		 * @brief 将 ITRS 地固位置旋转为指定时刻的 GCRS 地心天球位置。
		 * @param[in,out] tracer 记录本次旋转的调用链和输入错误。
		 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
		 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
		 * @param[out] gcrs_position GCRS 位置 [X,Y,Z]，长度单位及归一化状态与当前对象相同。
		 * @details 复用 ITRS_to_GCRS_Matrix 的 IAU 2000B 位置旋转；失败时保留输出对象。
		 */
		void ToGCRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, GCRS& gcrs_position) const;
	};

	class GCRS
	{
	 private:
		double m_coord[3];//X Y Z 暂不考虑速度
		bool m_is_normalized;

	 public:
		GCRS();
		//自定义的构造函数 X Y Z
		GCRS(double X, double Y, double Z, bool is_normalized);//参数中删除了时间
		//析构函数
		~GCRS();
		//修改设置成员值
		void Set(double X, double Y, double Z, bool is_normalized);
		//转为一般的数组
		void ToArray(double* gcrs_array) const;

		//获取成员值
		double GetValue(Overwatch& tracer, int idx) const;
		//归一化
		void Normalize(Overwatch& tracer, double uLength);
		//去归一化
		void Restore(Overwatch& tracer, double uLength);

		/**
		 * @brief 将指定时刻的 GCRS 地心天球位置转换为 WGS84 大地经纬高。
		 * @param[in,out] tracer 记录地球定向旋转和椭球反算的调用链及错误。
		 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
		 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
		 * @param[out] geodetic_position 大地坐标 [高度,经度,大地纬度]；高度以
		 * constants::Scale::UNIT_LENGTH_M 归一化，角度为 rad。
		 * @details 当前对象须已按 UNIT_LENGTH_M 归一化。依次调用 ToITRS 和 ITRS::ToGeodetic，
		 * 采用 IAU 2000B 位置旋转及 WGS84 椭球；失败时保留输出对象。
		 */
		void ToGeodetic(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, Geodetic& geodetic_position) const;

		/**
		 * @brief 将 GCRS 地心天球位置旋转为指定时刻的 ITRS 地固位置。
		 * @param[in,out] tracer 记录本次旋转的调用链和输入错误。
		 * @param[in] MJD_UTC 转换时刻，以 UTC MJD 表示，单位 d。
		 * @param[in] EOP_needed 同历元 EOP；dut1、dat 为 s，xp、yp、dx00、dy00 为 rad。
		 * @param[out] itrs_position 地固位置 [X,Y,Z]，长度单位及归一化状态与当前对象相同。
		 * @details 复用 GCRS_to_ITRS_Matrix 的 IAU 2000B 位置旋转；失败时保留输出对象。
		 */
		void ToITRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, ITRS& itrs_position) const;
	};
}

#endif // CDC_FRAME_DEFINITION_HPP
