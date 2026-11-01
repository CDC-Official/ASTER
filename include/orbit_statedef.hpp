//
//  orbit_statedef.hpp
//  CDC_lib
//
//  Created by 陈酊存 on 2023/9/5.
//

#ifndef orbit_statedef_hpp
#define orbit_statedef_hpp

#include <cmath>
#include <iostream>
#include "utils_constants.hpp"
#include "orbit_basic.hpp"
#include "utils_overwatch.hpp"
#include "algorithm"
using namespace bastion;

namespace orbit
{
	class Cartesian;
	class Keplerian;
	//class NonSingular1;
	class NonSingular2;

	class Cartesian
	{
	 private:
		//为了避免定义构造函数时，构造函数的参数和成员名一样，因此在成员名中加m_的前缀，或者加后缀_
		//double m_t_MJD_UTC;//MJD UTC!! 暂时删除时间，因为感觉不太方便创建和修改
		double m_x, m_y, m_z, m_vx, m_vy, m_vz;
		double m_r, m_v;
		bool m_is_normalized;  //是否归一化
	 public:
		//默认构造函数
		Cartesian();
		//自定义的构造函数
		Cartesian(double x, double y, double z, double vx, double vy, double vz, bool is_normalized);//参数中删除了时间
		//析构函数
		~Cartesian();
		// show
		void Show(std::ostream& os = std::cout, int precision = 16) const;         // const成员函数，保证该函数不会修改调用的对象
		//修改设置成员值
		void Set(double x, double y, double z, double vx, double vy, double vz, bool is_normalized);
		void ToArray(double* cart_array) const;

		double GetValue(Overwatch& tracer, int varId) const;//获取成员值
		void Normalize(Overwatch& tracer, double uLength, double uTime);//归一化
		void Normalized(Overwatch& tracer, double uLength, double uTime,
			Cartesian& new_cart) const;//归一化生成新的Cartesian，不改变本身

		void Restore(Overwatch& tracer, double uLength, double uTime);//去归一化

		void ToKeplerian(Overwatch& tracer, double* kep_array) const;     //转为Kepler根数，储存到数组
		void ToKeplerian(Overwatch& tracer, Keplerian& kep_class) const;     //转为Kepler根数，储存到Keplerian类

		void ToNonSingular2(Overwatch& tracer, double* nsg2_array) const;     //转为nsg2根数，储存到数组
		void ToNonSingular2(Overwatch& tracer, NonSingular2& nsg2_class) const;     //转为nsg2根数，储存到Keplerian类



	};

	class Keplerian
	{
	 private:
		//double m_t_MJD_UTC;
		double m_a, m_e, m_i, m_O, m_o, m_M;
		double m_f, m_E, m_u;//几个近点角，以及u=omega+f
		double m_r;//地心距离
		bool m_is_normalized;
		bool m_is_fullySet = false;

		// 避免数组和对象重载维护两套转换公式。
		bool ToNonSingular2Array(Overwatch& tracer, double* nsg2_array, bool strict_check) const;

	 public:
		Keplerian();
		//初始化还要考虑mod 2pi
		Keplerian(double a, double e, double i, double OMEGA, double omega, double M, bool is_normalized);

		~Keplerian();

		void Set(double a, double e, double i, double OMEGA, double omega, double M, bool is_normalized);

		void Show(std::ostream& os = std::cout, int precision = 16) const;         // const成员函数，保证该函数不会修改调用的对象

		void SyncAnomaly(Overwatch& tracer);//同步近点角,必须要手动执行！

		double GetValue(Overwatch& tracer, int varId) const;

		void Normalize(Overwatch& tracer, double uLength);//归一化 自身
		void Restore(Overwatch& tracer, double uLength);//去归一化

		void Normalized(Overwatch& tracer, double uLength, Keplerian& new_kep) const;//归一化生成新的Keplerian，不改变本身

		void ToCartesian(Overwatch& tracer, double* cart_array) const;     //转为直角坐标，储存到数组
		void ToCartesian(Overwatch& tracer, Cartesian& cart_class) const;     //转为直角坐标，储存到Cartesian类

		void ToNonSingular2(Overwatch& tracer, double* nsg2_array, bool strict_check = false) const;
		void ToNonSingular2(Overwatch& tracer, NonSingular2& nsg2_class, bool strict_check = false) const;

		void ToArray(double* kep_array) const;

		bool IsNormalized() const;

	};

	class NonSingular2
	{
	 private:
		// definition elements
		double m_a, m_xi, m_eta, m_h, m_k, m_lambda;
		// derived elements
		double m_r, m_uTilde, m_u;// uTilde = m_Eow = E + omega + w		u = fow = f + omega + w
		// whether the elements are fully set
		bool m_is_normalized;
		bool m_is_fullySet = false;
		bool m_strict_check = false;

		// 避免数组和对象重载维护两套转换公式。
		bool ToCartesianArray(Overwatch& tracer, double* cart_array) const;
		bool ToKeplerianArray(Overwatch& tracer, double* kep_array,
			double epsZeroEcc, double epsZeroInc) const;

	 public:
		NonSingular2();
		~NonSingular2();

		void Show(std::ostream& os = std::cout, int precision = 16) const;

		void ToArray(double* nsg2_array) const;

		NonSingular2(Overwatch& tracer, double a, double xi, double eta,
			double h, double k, double lambda, bool is_normalized);

		void Set(Overwatch& tracer, double a, double xi, double eta,
			double h, double k, double lambda, bool is_normalized);

		void SyncAnomaly(Overwatch& tracer);
		double GetValue(Overwatch& tracer, int varId) const;//获取成员值
		void Normalize(Overwatch& tracer, double uLength);//归一化
		void Restore(Overwatch& tracer, double uLength);//去归一化

		//转化为cart和kep
		void ToCartesian(Overwatch& tracer, Cartesian& cart_class) const;

		void ToCartesian(Overwatch& tracer, double* cart_array) const;

		void ToKeplerian(Overwatch& tracer, Keplerian& kep_class,
			double epsZeroEcc = 1E-6, double epsZeroInc = 1E-6) const;
		void ToKeplerian(Overwatch& tracer, double* kep_array,
			double epsZeroEcc = 1E-6, double epsZeroInc = 1E-6) const;

	};

}  // namespace orbit
#endif /* orbit_statedef_hpp */
