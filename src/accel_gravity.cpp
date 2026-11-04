//
// Created by 陈酊存 on 2024/1/23.
//

#include "accel_gravity.hpp"


// [补丁修改] 全局变量实体定义
LunarGravityField moon_gravity_field;


int Norm_LP_Array(int l, double mu, double* LP)
{
	/*
	 * 计算归一化的勒让德多项式数组
	 * int l:阶次，l>=0
	 * double x:自变量
	 * double *LP: 归一化勒让德多项式的数组，从P_0(mu)一直到P_l(mu) ,维度必须为l+1，即LP[l+1]
	 * 返回1 正常工作
	 * 返回-1 l<0
	 * 如果数组A在函数外定义，并且作参数传入函数，怎么在函数里求其长度？
	 * ——不好求，一般我们会再多传入一个参数记录长度。因此LP必须要事先定义为 double LP[L+1].
	 * 摘自blog：如果在函数中定义数组A，sizeof(A)是所有数组元素所占用的内存空间字节数，
	 * 因为编译器在将 C 代码转换成汇编代码时，自动将其替换成了实际的数值，
	 * 因此可以通过sizeof(A)/sizeof(A[0])计算数组大小。
	 * 但是A实际是一个指针，表示的是这个数组首个元素的地址，
	 * 不信你可以试着尝试cout<<A<<endl，你会发现输出的是一个地址。
	 * 所以如果此数组在函数外定义，然后作为函数的输入参数传入，
	 * 则传入的A只能表示指针，所以此时sizeof(A)为A[0]的地址的长度(x86地址长度为32位，即4个字节)
	 * 所以sizeof(A)/sizeof(A[0])为1（x86）
	 * 所以在我们传入这种数组做参数时，通常同时会再多传一个参数记录数组的长度。
	 * 经过测试，P_30(0.5)=LP[30]=1.170356658012492,正确
	 */
	//int size = sizeof(LP) / sizeof(double); 这样不行
	if (l < 0)//先判断l的值是否合法
	{
		std::cout << "Invalid value of l < 0!" << std::endl;
		return -1;
	}
	else if (l == 0)
	{
		LP[0] = 1;
		return 1;
	}
	else if (l == 1)
	{
		LP[0] = 1;
		LP[1] = sqrt(3.0) * mu;
		return 1;
	}
	else//对应l>=2
	{
		LP[0] = 1;
		LP[1] = sqrt(3.0) * mu;
		double coeff, first_term, second_term;
		for (int i = 2; i <= l; i++)
		{
			coeff = sqrt((2.0 * i + 1.0) / (2.0 * i - 1.0));
			first_term = (2.0 - 1.0 / i) * mu * LP[i - 1];
			second_term =
				-sqrt((2.0 * i - 1.0) / (2.0 * i - 3.0)) * (1.0 - 1.0 / i) * LP[i - 2];
			LP[i] = coeff * (first_term + second_term);
		}
		return 1;
	}
}

int Norm_LPD_Array(int l, double mu, const double* LP, double* LPD)
{
	/*
	 * 已知勒让德多形式数组，计算归一化的勒让德多项式导数数组
	 * 	注意程序暂时无法判断两个函数的自变量mu是否一致，以及维度，因此需要得到保证
	 * 	P'_30(0.5) = LPD[30] = -10.76087468041399 通过验证
	 */

	if (l < 0)//先判断l的值是否合法
	{
		std::cout << "Invalid value of l < 0!" << std::endl;
		return -1;
	}

	else if (l == 0)
	{
		LPD[0] = 0;
		return 1;
	}
	else if (l == 1)
	{
		LPD[0] = 0;
		LPD[1] = sqrt(3.0);
		return 1;
	}
	else
	{
		LPD[0] = 0;
		LPD[1] = sqrt(3.0);
		double coeff, first_term, second_term;
		for (int i = 2; i <= l; i++)
		{
			coeff = i / (1.0 - mu * mu);
			first_term = sqrt((2.0 * i + 1.0) / (2.0 * i - 1.0)) * LP[i - 1];
			second_term = -mu * LP[i];
			LPD[i] = coeff * (first_term + second_term);
		}
		return 1;
	}
}

int Norm_LPA_Array(const int l, double mu, std::vector<std::vector<double> >& LPA)
{
	/*
	 * 计算连带勒让德多项式矩阵LPA 必须初始化为全0，因为没有刻意填充为0的上三角
	 * LPA维度为l+1 * l+1，从LPA[0][0]一直到LPA[l][l],
	 * 要求l必须大于1，P_l0即第一列应当为勒让德函数的值，
	 * 使用动态数组，调用了Norm_LP_Array用于填充。
	 * 正常需要的部分为从P_11一直到P_ll的下三角（其余为0）
	 * 注意下标！
	 */
	if (l <= 0)
	{
		std::cout << "Invalid value of l <= 0!" << std::endl;
		return -1;
	}
	else if (l == 1)
	{  // LPA 2*2
		LPA[0][0] = 1;//普通勒让德
		LPA[0][1] = 0;
		LPA[1][0] = sqrt(3.0) * mu;//普通勒让德
		LPA[1][1] = sqrt(3.0) * sqrt(1.0 - mu * mu);
		return 1;
	}
	else
	{  // l >= 2
		LPA[1][1] = sqrt(3.0) * sqrt(1.0 - mu * mu);
		for (int i = 2; i <= l; i++)
		{  // 先算主对角线元素 从P_11 -> P_ll
			LPA[i][i] = sqrt((2.0 * i + 1.0) / (2.0 * i)) * sqrt(1.0 - mu * mu) * LPA[i - 1][i - 1];
		}

		double first_term, second_term;
		for (int col = 1; col < l; col++)//从第一列开始
		{  // 每一列往下递推 按照列 不包括最后一列，因为仅有Pll
			for (int row = col + 1; row <= l; row++)//从主对角线下面的元素开始
			{
				first_term = sqrt(((2.0 * row + 1.0) * (2.0 * row - 1.0)) / double((row + col) * (row - col))) * mu *
							 LPA[row - 1][col];
				second_term = -sqrt(((2.0 * row + 1.0) * (row - 1.0 + col) * (row - 1.0 - col)) /
									double((2.0 * row - 3.0) * (row + col) * (row - col))) * LPA[row - 2][col];
				LPA[row][col] = first_term + second_term;
			}
		}

		//new一个数组出来，填充P[0][0]一直到P[l][0]，即普通的勒让德
		double* dynamic_LP_Array = new double[l + 1];
		int a = Norm_LP_Array(l, mu, dynamic_LP_Array);//此时不会有隐患，不必处理报错
		for (int i = 0; i <= l; i++)
		{
			LPA[i][0] = dynamic_LP_Array[i];//填充P[0][0]一直到P[l][0]
		}
		delete[] dynamic_LP_Array;//释放内存
		return 1;
	}
}

int Norm_LPAD_Array(int l,
	double mu,
	const std::vector<std::vector<double>>& LPA,
	std::vector<std::vector<double>>& LPAD)
{
	/*
	 * 计算连带勒让德多项式导数矩阵LPAD 初始化为全0，下面会说明
	 * LPAD维度为l+1 * l+1，从LPAD[0][0]一直到LPAD[l][l],
	 * 要求l必须大于1，P'_l0即第''0''列应当为勒让德函数导数的值，
	 * 使用Norm_LPD_Array()进行了填充
	 * 需要的部分为从P_11一直到P_ll的下三角（其余为0）
	 * 注意下标！
	 */
	if (l < 0)
	{
		std::cout << "Invalid value of l < 0!" << std::endl;
		return -1;
	}
	else if (l == 0)
	{
		LPAD[0][0] = 0;
		return 1;
	}
	else if (l == 1)
	{
		LPAD[0][0] = 0;
		LPAD[0][1] = 0;
		LPAD[1][0] = sqrt(3.0);//即LPD[1]
		LPAD[1][1] = -mu / (1 - mu * mu) * LPA[1][1];

		return 1;
	}
	else
	{
		double coeff = 1.0 / sqrt(1.0 - mu * mu);
		double coeff_squared = 1.0 / (1.0 - mu * mu);
		double first_term, second_term;
		// 先把主对角线的计算
		for (int i = 1; i <= l; i++)
		{
			//first_term = 0;  // P_i,i+1=0
			//second_term = -i * mu * coeff * LPA[i][i];
			//LPAD[i][i] = coeff * (first_term + second_term);
			//写法等价，效率稍高
			LPAD[i][i] = -i * mu * coeff_squared * LPA[i][i];
		}
		// 再处理除了对角线之外的下三角
		for (int col = 1; col < l; col++)
		{  // 每一列往下递推 按照列 不包括最后一列
			for (int row = col + 1; row <= l; row++)
			{
				first_term = sqrt((row + col + 1.0) * (row - col)) * LPA[row][col + 1];
				second_term = -col * mu * coeff_squared * LPA[row][col];
				LPAD[row][col] = coeff * first_term + second_term;
			}
		}

		//填充P'[0][0]到P'[l][0]
		double* dynamic_LP_Array = new double[l + 1];
		double* dynamic_LPD_Array = new double[l + 1];

		// 使用 std::fill 初始化数组为全0
		std::fill(dynamic_LP_Array, dynamic_LP_Array + l + 1, 0.0);
		std::fill(dynamic_LPD_Array, dynamic_LPD_Array + l + 1, 0.0);

		int a = Norm_LP_Array(l, mu, dynamic_LP_Array);//此时不会有隐患，不必处理报错
		int b = Norm_LPD_Array(l, mu, dynamic_LP_Array, dynamic_LPD_Array);
		for (int i = 0; i <= l; i++)
		{
			LPAD[i][0] = dynamic_LPD_Array[i];
		}
		delete[] dynamic_LP_Array;
		delete[] dynamic_LPD_Array;
		return 1;
	}
}

double Norm_LP_Single(int l, double mu)
{
	//Norm_Legendre_Polynomials_Single 返回归一化的勒让德多形式值\bar{P}_{l}(\mu)
	if (l < 0)
	{
		std::cout << "Invalid value of l < 0!" << std::endl;
		exit(-1);
	}
	else if (l == 0)
	{
		return 1;
	}
	else if (l == 1)
	{
		return sqrt(3.0) * mu;
	}
	else
	{
		double coeff = sqrt((2.0 * l + 1.0) / (2.0 * l - 1.0));
		double first_term = (2.0 - 1.0 / l) * mu * Norm_LP_Single(l - 1, mu);
		double second_term =
			-sqrt((2.0 * l - 1.0) / (2.0 * l - 3.0)) * (1.0 - 1.0 / l) * Norm_LP_Single(l - 2, mu);
		return coeff * (first_term + second_term);
	}
}

double Norm_LPD_Single(int l, double mu)
{
	//Norm_Legendre_Polynomials_Derivatives_Single
	if (l < 0)
	{
		std::cout << "Invalid value of l < 0!" << std::endl;
		exit(-1);
	}
	else if (l == 0)
	{
		return 0;//导数=0
	}
	else
	{
		double coeff = l / (1.0 - mu * mu);
		double first_term = sqrt((2.0 * l + 1.0) / (2.0 * l - 1.0)) * Norm_LP_Single(l - 1, mu);
		double second_term = -mu * Norm_LP_Single(l, mu);
		return coeff * (first_term + second_term);
	}

}

double Norm_LPA_Single(int l, int m, double mu)
{
	//Norm_Legendre_Polynomials_Associated_Single
	if (abs(mu) > 1)//取值非法
	{
		std::cout << "Invalid value of mu!" << std::endl;
		exit(-1);
	}
	else if (l < 0 || m < 0)//取值非法
	{
		std::cout << "Invalid value of l < 0 or m < 0!" << std::endl;
		//return -1;
		exit(-1);
	}
	else if (l >= 0 && m == 0)//m=0
	{
		return Norm_LP_Single(l, mu);
	}
	else if (l == 1 && m == 1)
	{
		return sqrt(3.0) * sqrt(1.0 - mu * mu);
	}
	else if (l < m)
	{
		return 0;
	}
	else if (l == m)
	{
		return sqrt((2.0 * l + 1.0) / (2.0 * l)) * sqrt(1.0 - mu * mu)
			   * Norm_LPA_Single(l - 1, m - 1, mu);
	}
	else
	{
		double first_term = sqrt(((2.0 * l + 1.0) * (2.0 * l - 1.0)) / double((l + m) * (l - m))) * mu
							* Norm_LPA_Single(l - 1, m, mu);
		double second_term =
			-sqrt(((2.0 * l + 1.0) * (l - 1.0 + m) * (l - 1.0 - m)) / double((2.0 * l - 3.0) * (l + m) * (l - m)))
			* Norm_LPA_Single(l - 2, m, mu);
		return first_term + second_term;
	}
}

double Norm_LPAD_Single(int l, int m, double mu)
{
	if (abs(mu) > 1)//取值非法
	{
		std::cout << "Invalid value of mu = " << mu << " !" << std::endl;
		exit(-1);
	}
	else if (l < 0 || m < 0)//取值非法
	{
		std::cout << "Invalid value of l = " << l << ", " << "m = " << m << " !" << std::endl;
		//return -1;
		exit(-1);
	}
	else if (l == 0 && m == 0)//m=0
	{
		std::cout << "Special value of l = " << l << ", " << "m = " << m << " ! Return 0!" << std::endl;
		return 0;
	}
	else if (l == 0 && m > 0)
	{
		std::cout << "Dangerous value of l = " << l << ", " << "m = " << m << " ! Return 0!" << std::endl;
		return 0;
	}
	else if (l > 0 && m == 0)//m=0 退化为勒让德多项式
	{
		return Norm_LPD_Single(l, mu);
	}
	else //l>0,m>0
	{
		double coeff = 1.0 / sqrt(1.0 - mu * mu);
		double first_term = sqrt((l + m + 1.0) * (l - m)) * Norm_LPA_Single(l, m + 1, mu);
		double second_term = -m * mu / sqrt(1 - mu * mu) * Norm_LPA_Single(l, m, mu);
		return coeff * (first_term + second_term);
	}
}


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
	int field_id, int zonal_order, int tesseral_degree, double field_gm, double reference_radius)
{
	// 球谐展开中的距离以模型参考半径为单位，外部积分状态仍以 ae 为单位。
	const double length_scale = ae / reference_radius;
	double X = R_terrestrial[0] * length_scale;
	double Y = R_terrestrial[1] * length_scale;
	double Z = R_terrestrial[2] * length_scale;

	//R R'
	double R = sqrt(X * X + Y * Y + Z * Z);
	double R_prime = sqrt(X * X + Y * Y);

	// 这里的 phi 是地心纬度：sin(phi)=Z/R；lambda 为地固球坐标经度。
	// 直接用于球谐函数，不需要先转换成大地经纬高。
	double sin_phi = Z / R;
	double cos_lambda = X / R_prime;
	double sin_lambda = Y / R_prime;

	// R\vec  k\hat G\vec
	Vector3d k_hat(0, 0, 1);
	Vector3d R_vec(X, Y, Z);
	Vector3d G_vec(-sin_lambda, cos_lambda, 0);

	// l m
	int l = 0, m = 0;

	//先处理带谐项
	Vector3d a_zonal(0, 0, 0);//带谐项加速度
	if (zonal_order < 0)//非法的阶次
	{
		std::cout << "Invalid value of zonal order < 0!" << std::endl;
		return -1;
	}
	else if (zonal_order == 0 || zonal_order == 1)//zonal_order=0 or 1 无摄动加速度，a_zonal=(0,0,0);
	{
		std::cout << "Caution: zonal order = " << zonal_order
				  << ", perturbation acceleration due to non-spherical gravity is equal to zero!" << std::endl;
	}
	else
	{
		if (zonal_order > 80)//超出上线，重新设置为80
		{
			std::cout << "Zonal order = " << zonal_order
					  << ", exceeded the maximum value of 80, zonal order is reset to 80 here." << std::endl;
			zonal_order = 80;
		}
		//计算P与P'
		double* LP = new double[zonal_order + 1];
		double* LPD = new double[zonal_order + 1];

		// 使用 std::fill 初始化数组为全0
		std::fill(LP, LP + zonal_order + 1, 0.0);
		std::fill(LPD, LPD + zonal_order + 1, 0.0);

		//	不会有隐患，没事，不用判断ab的值
		int a = Norm_LP_Array(zonal_order, sin_phi, LP);
		int b = Norm_LPD_Array(zonal_order, sin_phi, LP, LPD);

		double temp_Cl0;
		for (l = 2; l <= zonal_order; l++)
		{
			temp_Cl0 = GGM03C_arr[(3 + l) * (l - 2) / 2][0];
			//注意是减号
			a_zonal -= temp_Cl0 / pow(R, l + 3)
					   * (((l + 1) * LP[l] + sin_phi * LPD[l]) * R_vec - R * LPD[l] * k_hat);
		}

		delete[] LP;
		delete[] LPD;

	}

	//再处理田谐项
	Vector3d a_tesseral(0, 0, 0);
	if (tesseral_degree < 0)
	{
		std::cout << "Invalid value of tesseral order < 0!" << std::endl;
		return -1;
	}
	else if (tesseral_degree == 0 || tesseral_degree == 1)//tesseral_order=0 or 1 无摄动加速度，a_tesseral=(0,0,0);
	{
		std::cout << "Caution: tesseral order = " << tesseral_degree
				  << ", perturbation acceleration due to non-spherical gravity is equal to zero!" << std::endl;
	}
	else
	{
		if (tesseral_degree > 80)//超出上线，重新设置为80
		{
			std::cout << "Tesseral order = " << tesseral_degree
					  << ", exceeded the maximum value of 80, tesseral order is reset to 80 here." << std::endl;
			tesseral_degree = 80;
		}
		//LPA LPAD
		std::vector<std::vector<double> > LPA(tesseral_degree + 1, std::vector<double>(tesseral_degree + 1, 0.0));
		std::vector<std::vector<double> > LPAD(tesseral_degree + 1, std::vector<double>(tesseral_degree + 1, 0.0));
		// 240927 发现田谐项摄动加速度为0！问题在这里，连带勒让德多项式LPA和其导数LPAD只初始化，并没有赋值！！！
		// 240927 添加田谐项的连带勒让德计算，未经验证
		int a = Norm_LPA_Array(tesseral_degree, sin_phi, LPA);
		int b = Norm_LPAD_Array(tesseral_degree, sin_phi, LPA, LPAD);

		double* sin_lambda_array = new double[tesseral_degree + 1];
		double* cos_lambda_array = new double[tesseral_degree + 1];

		sin_lambda_array[0] = 0;
		sin_lambda_array[1] = sin_lambda;

		cos_lambda_array[0] = 1;
		cos_lambda_array[1] = cos_lambda;

		//sin_lambda_array[m]=sin(m*Lambda_G) 1.138
		for (m = 2; m <= tesseral_degree; m++)
		{
			sin_lambda_array[m] = 2 * cos_lambda * sin_lambda_array[m - 1] - sin_lambda_array[m - 2];
			cos_lambda_array[m] = 2 * cos_lambda * cos_lambda_array[m - 1] - cos_lambda_array[m - 2];
		}

		//开始(1.153)
		double temp_C_lm, temp_S_lm;
		for (l = 2; l <= tesseral_degree; l++)
		{

			for (m = 1; m <= l; m++)
			{
				temp_C_lm = GGM03C_arr[(3 + l) * (l - 2) / 2 + m][0];
				temp_S_lm = GGM03C_arr[(3 + l) * (l - 2) / 2 + m][1];
				//注意是减号
				a_tesseral -= 1.0 / (pow(R, l + 3))
							  * (((l + 1) * LPA[l][m] + sin_phi * LPAD[l][m]) * R_vec - R * LPAD[l][m] * k_hat)
							  * (temp_C_lm * cos_lambda_array[m] + temp_S_lm * sin_lambda_array[m])
							  + m / R_prime / pow(R, l + 1) * LPA[l][m]
								* (temp_C_lm * sin_lambda_array[m] - temp_S_lm * cos_lambda_array[m]) * G_vec;
				// FIXME: DEBUG用的输出
				// std::cout << a_tesseral[0] << '\t' << a_tesseral[1] << '\t' << a_tesseral[2] << std::endl;
			}
		}

		delete[] sin_lambda_array;
		delete[] cos_lambda_array;
	}


	// 将模型的 GM/reference_radius² 加速度尺度换成积分器使用的 GE/ae²。
	const double acceleration_scale = (field_gm / GE) * length_scale * length_scale;
	Vector3d a_total_terrestrial = (a_zonal + a_tesseral) * acceleration_scale;
	a_terrestrial[0] = a_total_terrestrial[0];
	a_terrestrial[1] = a_total_terrestrial[1];
	a_terrestrial[2] = a_total_terrestrial[2];

	return 1;
}

// [补丁] 月球非球形引力计算函数 (参数命名已强化月球归一化语义)
int Lunar_Nonspherical_Gravitational_Acceleration(double* a_moon_fixed_moon_norm, const double* rsat_moon_fixed_moon_norm, int zonal_order, int tesseral_degree)
{
	// 0. 安全检查
	if (moon_gravity_field.max_degree == 0)
	{
		std::cerr << "Error: Lunar gravity field not loaded! Call moon_gravity_field.Load() first." << std::endl;
		return -1;
	}

	// 1. 提取坐标 (明确这是以月球半径为单位的)
	double X = rsat_moon_fixed_moon_norm[0];
	double Y = rsat_moon_fixed_moon_norm[1];
	double Z = rsat_moon_fixed_moon_norm[2];

	// 2. 计算辅助量
	double R = sqrt(X * X + Y * Y + Z * Z); // 此时 R=1.0 代表在月球表面
	double R_prime = sqrt(X * X + Y * Y);

	// 这里以月心为原点，phi、lambda 为月固轴系中的月心球坐标纬度和经度。
	// 保留与地球球谐公式相同的简短数学符号，其原点与轴向由本计算块说明。
	double sin_phi = Z / R;
	if (sin_phi > 1.0) sin_phi = 1.0;
	if (sin_phi < -1.0) sin_phi = -1.0;

	double cos_lambda = (R_prime == 0) ? 1.0 : X / R_prime;
	double sin_lambda = (R_prime == 0) ? 0.0 : Y / R_prime;

	Vector3d k_hat(0, 0, 1);
	Vector3d R_vec(X, Y, Z);
	Vector3d G_vec(-sin_lambda, cos_lambda, 0);

	int l = 0, m = 0;

	// ==========================================
	// 4. 处理带谐项 (Zonal Terms)
	// ==========================================
	Vector3d a_zonal(0, 0, 0);

	if (zonal_order > moon_gravity_field.max_degree)
	{
		zonal_order = moon_gravity_field.max_degree;
	}

	if (zonal_order >= 2)
	{
		double* LP = new double[zonal_order + 1];
		double* LPD = new double[zonal_order + 1];
		std::fill(LP, LP + zonal_order + 1, 0.0);
		std::fill(LPD, LPD + zonal_order + 1, 0.0);

		Norm_LP_Array(zonal_order, sin_phi, LP);
		Norm_LPD_Array(zonal_order, sin_phi, LP, LPD);

		for (l = 2; l <= zonal_order; l++)
		{
			double C_l0 = moon_gravity_field.C[l][0];

			// 注意：这里的 pow(R, l+3) 中的 R 是月球归一化距离
			a_zonal -= C_l0 / pow(R, l + 3)
					   * (((l + 1) * LP[l] + sin_phi * LPD[l]) * R_vec - R * LPD[l] * k_hat);
		}
		delete[] LP;
		delete[] LPD;
	}

	// ==========================================
	// 5. 处理田谐项 (Tesseral Terms)
	// ==========================================
	Vector3d a_tesseral(0, 0, 0);

	if (tesseral_degree > moon_gravity_field.max_degree)
	{
		tesseral_degree = moon_gravity_field.max_degree;
	}

	if (tesseral_degree >= 2)
	{
		std::vector<std::vector<double> > LPA(tesseral_degree + 1, std::vector<double>(tesseral_degree + 1, 0.0));
		std::vector<std::vector<double> > LPAD(tesseral_degree + 1, std::vector<double>(tesseral_degree + 1, 0.0));

		Norm_LPA_Array(tesseral_degree, sin_phi, LPA);
		Norm_LPAD_Array(tesseral_degree, sin_phi, LPA, LPAD);

		double* sin_mL = new double[tesseral_degree + 1];
		double* cos_mL = new double[tesseral_degree + 1];
		sin_mL[0] = 0;
		cos_mL[0] = 1;
		sin_mL[1] = sin_lambda;
		cos_mL[1] = cos_lambda;

		for (m = 2; m <= tesseral_degree; m++)
		{
			sin_mL[m] = 2 * cos_lambda * sin_mL[m - 1] - sin_mL[m - 2];
			cos_mL[m] = 2 * cos_lambda * cos_mL[m - 1] - cos_mL[m - 2];
		}

		for (l = 2; l <= tesseral_degree; l++)
		{
			for (m = 1; m <= l; m++)
			{
				double C_lm = moon_gravity_field.C[l][m];
				double S_lm = moon_gravity_field.S[l][m];

				double term_CS = (C_lm * cos_mL[m] + S_lm * sin_mL[m]);
				double term_CS_deriv = (C_lm * sin_mL[m] - S_lm * cos_mL[m]);

				a_tesseral -= 1.0 / pow(R, l + 3)
							  * (((l + 1) * LPA[l][m] + sin_phi * LPAD[l][m]) * R_vec - R * LPAD[l][m] * k_hat)
							  * term_CS
							  + (double) m / R_prime / pow(R, l + 1) * LPA[l][m]
								* term_CS_deriv * G_vec;
			}
		}
		delete[] sin_mL;
		delete[] cos_mL;
	}

	// 6. 输出赋值
	Vector3d a_total = a_zonal + a_tesseral;

	// 输出的是：以月球表面重力加速度为 1.0 的归一化值 后续在GCRS积分时需要转换为地球归一化的单位
	a_moon_fixed_moon_norm[0] = a_total[0];
	a_moon_fixed_moon_norm[1] = a_total[1];
	a_moon_fixed_moon_norm[2] = a_total[2];

	return 1;
}
