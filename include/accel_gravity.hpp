//
// Created by 陈酊存 on 2024/1/23.
//

#ifndef CDC_ACCEL_GRAVITY_HPP
#define CDC_ACCEL_GRAVITY_HPP

#include <cmath>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "utils_constants.hpp"
#include "orbit_statedef.hpp"
#include "data_gravityfield.hpp"

using namespace Eigen;


// [补丁修改] 结构体改名为 LunarGravityField (月球引力场)
struct LunarGravityField
{
	double GM;          // 系数文件的 GM 尺度，单位 m³/s²，用于二阶及以上球谐项；质点引力使用公共 DE440 GM。
	double R_ref;       // 参考半径 m
	int max_degree;     // 最大阶数
	bool is_ready = false; // 是否加载
	// 系数表 (C_nm, S_nm)
	std::vector<std::vector<double>> C;
	std::vector<std::vector<double>> S;

	LunarGravityField() : GM(0.0), R_ref(0.0), max_degree(0)
	{}

	/**
	 * @brief 读取月球重力场常数和归一化球谐系数。
	 * @param[in] filepath 月球重力场文件路径，支持 ICGEM 或 NASA PDS/SHADR 文本格式。
	 * @return 文件可读、GM 非零且 max_degree 大于零时返回 true，否则返回 false。
	 * @details 写入 GM（m³/s²）、R_ref（m）、max_degree 及无量纲系数 C、S，成功后置 is_ready。
	 * 加载成功事件由调用方通过 tracer 记录；本函数只在读取失败时输出诊断信息。
	 */
	bool Load(const std::string& filepath)
	{
		std::ifstream file(filepath);
		if (!file.is_open())
		{
			std::cerr << "[LunarField] Error: Cannot open " << filepath << std::endl;
			return false;
		}

		// 预分配内存
		int alloc_size = 720;
		C.assign(alloc_size, std::vector<double>(alloc_size, 0.0));
		S.assign(alloc_size, std::vector<double>(alloc_size, 0.0));

		std::string line;
		bool is_icgem = false; // 标记是否为 ICGEM 格式

		while (std::getline(file, line))
		{
			// 1. 预处理：把 Fortran 的 'D' 换成 'E'，把逗号换成空格
			for (char& ch : line)
			{
				if (ch == 'D' || ch == 'd') ch = 'E';
				if (ch == ',') ch = ' ';
			}

			// 跳过空行或注释行
			if (line.empty() || line[0] == '#' || line.substr(0, 3) == "***") continue;

			std::stringstream ss(line);
			std::string tag;

			// 尝试读取第一个单词看看是什么
			// 这里使用 peek 检查是否是数字
			ss >> std::ws; // 跳过前导空格
			if (ss.eof()) continue;

			char first_char = ss.peek();

			// === 分支 A: ICGEM 格式 (以单词开头) ===
			if (isalpha(first_char))
			{
				ss >> tag;
				if (tag == "earth_gravity_constant" || tag == "gravity_constant")
				{
					ss >> GM;
				}
				else if (tag == "radius")
				{
					ss >> R_ref;
				}
				else if (tag == "gfc")
				{
					int n, m;
					double c_val, s_val;
					if (ss >> n >> m >> c_val >> s_val)
					{
						if (n < alloc_size && m < alloc_size)
						{
							C[n][m] = c_val;
							S[n][m] = s_val;
							if (n > max_degree) max_degree = n;
						}
					}
				}
			}
				// === 分支 B: NASA PDS/SHADR 格式 (纯数字) ===
			else
			{
				// 如果还没读取到 Header 常数，且这一行只有几个数，认为是 Header
				// PDS Header: R_ref, GM, 1.0, MaxN, ...
				if (GM == 0.0 && R_ref == 0.0)
				{
					// 尝试读取 Header
					double v1, v2;
					if (ss >> v1 >> v2)
					{
						// 简单的启发式判断：半径通常在 1e6 量级，GM 在 1e12 量级
						// 或者是 PDS 格式: R, GM
						if (v1 < v2)
						{ // 也就是 1.7e6 < 4.9e12
							R_ref = v1;
							GM = v2;
						}
					}
				}
					// 否则认为是系数行: n, m, C, S, ...
				else
				{
					int n, m;
					double c_val, s_val; // c_err, s_err
					if (ss >> n >> m >> c_val >> s_val)
					{
						if (n < alloc_size && m < alloc_size)
						{
							C[n][m] = c_val;
							S[n][m] = s_val;
							if (n > max_degree) max_degree = n;
						}
					}
				}
			}
		}

		// 补全中心项
		if (C[0][0] == 0.0) C[0][0] = 1.0;

		// 如果还没读到，说明彻底失败了
		if (max_degree == 0 || GM == 0.0)
		{
			std::cerr << "[LunarField] Critical Warning: Failed to parse file content correctly!" << std::endl;
			return false;
		}

		is_ready = true;
		return true;
	}
};

// [补丁修改] 全局变量改名为 moon_gravity_field
extern LunarGravityField moon_gravity_field;


//递推计算数组的，高效
int Norm_LP_Array(int l, double mu, double* LP);

int Norm_LPD_Array(int l, double mu, const double* LP, double* LPD);

int Norm_LPA_Array(const int l, double mu, std::vector<std::vector<double> >& LPA);//安全，不会越界
int Norm_LPAD_Array(int l,
	double mu,
	const std::vector<std::vector<double>>& LPA,
	std::vector<std::vector<double>>& LPAD);//安全，不会越界

//递归计算单个值的，低效
double Norm_LP_Single(int l, double mu);

double Norm_LPD_Single(int l, double mu);

double Norm_LPA_Single(int l, int m, double mu);

double Norm_LPAD_Single(int l, int m, double mu);

//计算非球形引力加速度，归一化
int Earth_Nonspherical_Gravitational_Acceleration(double* a_celestial,
	const double* R_terrestrial,
	const Matrix3d& HG_transpose,
	int field_id, int zonal_order, int tesseral_degree);

/**
 * @brief 计算地固轴系中的 GGM03C 非球形引力加速度。
 * @param[out] a_terrestrial 三维摄动加速度 [ax, ay, az]，以 GE/ae² 归一化，不含中心引力。
 * @param[in] R_terrestrial 三维地心位置 [x, y, z]，按调用方的时空约定取值，并以 ae 归一化。
 * @param[in] field_id 引力场编号，当前实现采用 GGM03C。
 * @param[in] zonal_order 带谐项 C_n0 的最高次数 n，支持 0—80；0 或 1 关闭带谐项。
 * @param[in] tesseral_degree 田谐项和扇谐项的最高次数 n，保留 m=1—n，支持 0—80。
 * @param[in] field_gm 与球谐系数配套的 GM，单位 m³/s²，默认取 GGM03C 文件头的值。
 * @param[in] reference_radius 与球谐系数配套的参考半径，单位 m，默认取 GGM03C 文件头的值。
 * @return 成功返回 1，截断次数为负时返回 -1。
 * @details 先将位置转为参考半径归一化，再计算二阶及以上球谐项，最后转为 GE/ae² 尺度。
 *          field_gm、reference_radius 和位置须采用同一尺度约定；TT 调用使用默认参数，
 *          TDB 调用传入按 (1-L_B)/(1-L_G) 缩放的参数。2/0 的截断设置仅保留 J2 项。
 */
int Earth_Nonspherical_Gravitational_Acceleration(double* a_terrestrial, const double* R_terrestrial,
	int field_id, int zonal_order, int tesseral_degree,
	double field_gm = constants::GGM03C::GM, double reference_radius = constants::GGM03C::REFERENCE_RADIUS);


/**
 * @brief 计算月固坐标系下，月球非球形引力摄动加速度
 * * @param[out] a_moon_fixed_moon_norm
 * 输出参数：月固系下的**月球归一化**摄动加速度向量 [ax, ay, az]。
 * - 物理含义：实际加速度 / 月球表面重力加速度(GM_moon / R_moon^2)。
 * - 归一化基准：1.0 代表月球表面的重力加速度 (约 1.62 m/s^2)。
 * - 注意：在用于地球归一化积分器前，需乘以缩放因子 (GM_m/GM_e)*(R_e/R_m)^2。
 * * @param[in] rsat_moon_fixed_moon_norm
 * 输入参数：月固系下的**月球归一化**位置向量 [x, y, z]。
 * - 物理含义：卫星在月固系的位置坐标(米) / 月球参考半径 R_moon (1738000.0米)。
 * - 归一化基准：1.0 代表距离月心 1 个月球半径。
 * - 注意：必须是 pos_body_EarthNorm * (R_e/R_m) 后的值。
 * * @param[in] zonal_order     带谐项 (Zonal) 计算阶数 (如 20)
 * @param[in] tesseral_degree 田谐项 (Tesseral) 计算阶数 (如 20)
 * * @return int 返回 1 表示计算成功，返回 -1 表示模型未加载或参数错误
 */
int Lunar_Nonspherical_Gravitational_Acceleration(
	double* a_moon_fixed_moon_norm,
	const double* rsat_moon_fixed_moon_norm,
	int zonal_order,
	int tesseral_degree
);

#endif //CDC_LIB_INCLUDE_GRAVITY_HPP_
