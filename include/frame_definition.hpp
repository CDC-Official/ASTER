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

	// 大地坐标使用 constants 中现有常量所定义的参考椭球。
	// 分量顺序为 height、lambda（经度）、phi（大地纬度）。
	// 未归一化时为 m、deg、deg；归一化后为 height/ae、rad、rad。
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
		// 高度除以 uLength，角度由 deg 转为 rad；进行椭球转换时 uLength 应取 ae。
		void Normalize(Overwatch& tracer, double uLength);
		// 高度乘回 uLength，角度由 rad 恢复为 deg。
		void Restore(Overwatch& tracer, double uLength);
		// 按项目参考椭球将大地经纬高转为 ITRS XYZ；长度须先按 ae 归一化。
		void ToITRS(Overwatch& tracer, ITRS& itrs_position) const;
		//转为GCRS
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
		//转为大地坐标
		void ToGeodetic(Overwatch& tracer, Geodetic& geodetic_position) const;
		//转为GCRS，储存到GCRS类
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
		//转为大地坐标
		void ToGeodetic(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, Geodetic& geodetic_position) const;
		//转为ITRS，储存到ITRS类
		void ToITRS(Overwatch& tracer, double MJD_UTC, const aEOP& EOP_needed, ITRS& itrs_position) const;
	};
}

#endif // CDC_FRAME_DEFINITION_HPP
