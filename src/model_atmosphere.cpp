//
// Created by 陈酊存 on 2024/1/23.
//

#include "model_atmosphere.hpp"
#include <limits>

using namespace constants::Scale;
using namespace constants::Conversion;

/**
 * @brief 读取分层大气密度表及表头注明的空间天气条件。
 * @param[in] atm_path 密度表文件路径；数值记录依次为参考高度 km、密度 kg/m³、标高 km。
 * @param[in,out] atm_list 在数组末尾追加读到的分层记录；每条记录同时保存高度间距 dh_km。
 * @param[out] dh 从表头读取的高度间距，单位 km；读取失败时不保证有有效值。
 * @param[out] ap 可选接收表头的 Ap 指数，无量纲；传入 nullptr 时忽略，缺失或无效时写入 NaN。
 * @param[out] f107 可选接收表头的 F10.7 数值，按文件原标度保存；传入 nullptr 时忽略，缺失或无效时写入 NaN。
 * @return 成功返回 1；无法打开文件或数值记录读取失败时返回 -1。
 * @details Ap、F10.7 仅说明本表采用的固定条件，不改变密度值、插值方式或参考密度计算。
 */
int Read_AtmosphereDensity(const std::string& atm_path, std::vector<AtmosphereData>& atm_list, double& dh, double* ap, double* f107)
{
    // 可选元数据先置为空值；密度数值读取不依赖表头是否含有这两项说明。
    if (ap) *ap = std::numeric_limits<double>::quiet_NaN();

    if (f107) *f107 = std::numeric_limits<double>::quiet_NaN();

    std::ifstream file(atm_path);

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << atm_path << std::endl;
        return -1;
    }

    // 跳过说明行时，同时取得形如 Ap = 10, F10.7 = 150 的文件元数据。
    std::string line;

    while (std::getline(file, line))
    {
        const std::size_t ap_pos = line.find("Ap");
        const std::size_t f107_pos = line.find("F10.7");

        if (ap && ap_pos != std::string::npos)
        {
            const std::size_t equal = line.find('=', ap_pos + 2);
            double value;

            if (equal != std::string::npos && std::istringstream(line.substr(equal + 1)) >> value && std::isfinite(value)) *ap = value;
        }

        if (f107 && f107_pos != std::string::npos)
        {
            const std::size_t equal = line.find('=', f107_pos + 5);
            double value;

            if (equal != std::string::npos && std::istringstream(line.substr(equal + 1)) >> value && std::isfinite(value)) *f107 = value;
        }

        // 高度间距仍从既有标记之后的一行读取。
        if (line.find("# The step (km) is") != std::string::npos)
        {
            std::getline(file, line);
            std::istringstream input(line);

            if (!(input >> dh))
            {
                std::cerr << "Error reading dh_km parameter from file: " << atm_path << std::endl;
                return -1;
            }

            break;
        }
    }

    // 跳过列标题，按原顺序追加参考高度、密度和标高。
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::istringstream input(line);
        AtmosphereData data;
        data.dh_km = dh;

        if (!(input >> data.h_km >> data.rho_kgm3 >> data.H_km))
        {
            std::cerr << "Error reading data from file: " << atm_path << std::endl;
            return -1;
        }

        atm_list.push_back(data);
    }

    return 1;
}

void Set_ReferenceDensity(const std::vector<AtmosphereData>& atm_list,
	const orbit::Keplerian& initKep_norm,
	ReferenceDensity& refDnesity)
{
/*!
 * 由初始的kepler根数计算大气密度参考值
 * @param[in] atm_list			读取完成的大气参数vector,单位km kg/m^3
 * @param[in] initKep_norm		初始的kepler根数，用于预报，归一化单位
 * @param[out] refDnesity		待设定的参考密度结构体，均为归一化单位
 */

	double kep[6];
	initKep_norm.ToArray(kep);
	double sin_i0 = sin(kep[2]);
	double sin_w0 = sin(kep[4]);
	double rp0_norm = kep[0] * (1.0 - kep[1]);    //a(1-e)

	// 参考椭球的物理半径固定，转成积分器使用的长度单位后计算近地点高度。
	const double ellipsoid_shape = 1.0 - constants::Wgs84::FLATTENING * pow(sin_i0 * sin_w0, 2);
	double rRefEllipsoid_norm = constants::Wgs84::EQUATORIAL_RADIUS / UNIT_LENGTH_M * ellipsoid_shape;
	double hp0_km = (rp0_norm - rRefEllipsoid_norm) * NORMALIZED_LENGTH_TO_KM;//当前计算的初始近地点距离参考椭球面的高度km 1.225

	//refdensity赋值，都是归一化
	refDnesity.sin_i0 = sin_i0;
	refDnesity.sin_w0 = sin_w0;
	refDnesity.cos_i0 = cos(kep[2]);
	refDnesity.rp0 = rp0_norm;//参考点的地心距r
	refDnesity.hp0 = hp0_km * KM_TO_NORMALIZED_LENGTH;

	// 密度参考面的赤道半径；乘当前纬度的形状因子后得到该处的参考地心距。
	refDnesity.sigma0 = refDnesity.rp0 / ellipsoid_shape;

	if (hp0_km < atm_list[0].h_km)//初始近地点高度已经小于表格中的最小值
	{
		std::cerr << "Warning: The initial perigee altitude " << hp0_km
				  << " km is smaller than the minimum value " << atm_list[0].h_km
				  << "km in the table! ";
		std::cerr << "Return values at the lowest altitude in the table!" << std::endl;
		refDnesity.rho = atm_list[0].rho_kgm3;//有待归一化
		refDnesity.H = atm_list[0].H_km;//有待归一化
	}
	else if (hp0_km > atm_list[atm_list.size() - 1].h_km)//初始近地点高度已经大于表格中的最小值
	{
		std::cerr << "Warning: The initial perigee altitude " << hp0_km
				  << " km is larger than the maximum value " << atm_list[atm_list.size() - 1].h_km
				  << "km in the table! ";
		std::cerr << "Return values at the highest altitude in the table!!" << std::endl;
		refDnesity.rho = atm_list[atm_list.size() - 1].rho_kgm3;//有待归一化
		refDnesity.H = atm_list[atm_list.size() - 1].H_km;//有待归一化
	}
	else//在表格范围中
	{
		//准备插值
		// interpolate (1/2): find the best insertion point
		double dh_km = atm_list[0].dh_km;
		// instead of j = floor(), this makes sure that j+1 <= nRec
		int j = (int)ceil((hp0_km - atm_list[0].h_km) / dh_km) - 1;
		double k = (hp0_km - atm_list[j].h_km) / dh_km;
		refDnesity.H = atm_list[j].H_km + k * (atm_list[j + 1].H_km - atm_list[j].H_km);//线性插值得到标高，有待归一化
		//使用小区间计算更为精确的标高Hc，再用指数模型计算参考点处的大气密度
		double Hc = -dh_km / log(atm_list[j + 1].rho_kgm3 / atm_list[j].rho_kgm3);
		refDnesity.rho = atm_list[j].rho_kgm3 * exp(-(hp0_km - atm_list[j].h_km) / Hc);//有待归一化
	}

	// 密度和标高的归一化
	refDnesity.H *= KM_TO_NORMALIZED_LENGTH;
	refDnesity.rho *= UNIT_LENGTH_M * UNIT_LENGTH_M * UNIT_LENGTH_M / UNIT_MASS_KG;//kg/m^3的归一化！

}

double Get_AtmosphereDensity(const double* cart_norm,
	const ReferenceDensity& refDensity,
	double MJD_UT1,
	const std::vector<AtmosphereData>& atm_list)
{
	double rho_norm;
	//先考虑简单的周日变化(1.222)
	double r;
	Vector3d r_hat, rm_hat;
	r_hat << cart_norm[0], cart_norm[1], cart_norm[2];//量纲是归一化，向量还没归一化
	r = r_hat.norm();//位置的模长
	r_hat.normalize();//向量归一化

	double R_Sun[3] = { 0 };//平均太阳的位置
	Cart_MeanSun(MJD_UT1, R_Sun);//用分析历表获取
	double lambda_m = 30 * DD2R;//30度的弧度

	rm_hat << R_Sun[0] * cos(lambda_m) - R_Sun[1] * sin(lambda_m),
		R_Sun[1] * cos(lambda_m) + R_Sun[0] * sin(lambda_m),
		R_Sun[2];// 1.222
	rm_hat.normalize(); //再归一化

	double cos_psi_star = r_hat.dot(rm_hat);
	//double F_star = 0.3;//松鼠建议取值(m_rhoratio-1.0)/(m_rhoratio+1.0);
	//以下两个变量不用经常重复计算，可以设置为static
	static double rhoratio = 1.6;//rhomax/rhomin，取1.6，老师的取值
	static double F_star = (rhoratio - 1.0) / (rhoratio + 1.0);//老师modelatm.hpp line154取值，试试

	//再考虑高度变化
	//1.223
	double kep[6];//当前实时位置的根数
	CartesianToKeplerian(cart_norm[0], cart_norm[1], cart_norm[2],
		cart_norm[3], cart_norm[4], cart_norm[5], kep);
	double f = MeanAnomaly2TrueAnomaly(kep[5], kep[1]);
	// sin(phi)=sin(i)sin(omega+f)，是相对于输入赤道平面的地心纬度正弦。
	// 此处用它修正参考密度面的半径，没有求大地纬度或另转瞬时地固轴系。
	double sin_phi = sin(kep[2]) * sin(f + kep[4]);
	double sigma_p0 = refDensity.sigma0 * (1.0 - constants::Wgs84::FLATTENING * sin_phi * sin_phi);

	//1.220,参考密度*周日变化系数*高度衰减系数
	rho_norm = refDensity.rho * (1.0 + F_star * cos_psi_star) * exp(-(r - sigma_p0) / refDensity.H);
	return rho_norm;
}

//计算大气密度
double Get_AtmosphereDensity(const double* rv_Norm,
	const double* initKepArray_norm,
	double MJD_UT1,
	const std::vector<AtmosphereData>& atm_list)
{
	/*
	 * 计算大气密度，归一化
	 * rho_Norm 归一化的大气密度
	 * r_Norm 归一化的卫星J2000位置
	 * MJD_UT1 计算太阳平轨道用，严格来说（似乎）要用UT1，但是大气的不确定因素太大，UTC也完全没问题
	 * r_p0_list			参考高度		km
	 * density_p0_list		参考密度		kg/m^3
	 * H_p0_list			参考标高		km
	 *
	 */

	double rho_Norm;//最终返回的密度值

	//static变量仅会被初始化一次，用于存储查表得到的值 km kg/m^3
	static double h_p0_static = 0, rho_p0_static = 0, H_p0_static = 0;
	static double h_p0_static_Norm = 0, rho_p0_static_Norm = 0, H_p0_static_Norm = 0;//对应的归一化值

	double h_p0_km_temp;//当前计算的初始近地点高度km 1.225
	h_p0_km_temp = initKepArray_norm[0] * (1.0 - initKepArray_norm[1])
		- constants::Wgs84::EQUATORIAL_RADIUS / UNIT_LENGTH_M * (1.0 - constants::Wgs84::FLATTENING * pow(sin(initKepArray_norm[2]) * sin(initKepArray_norm[4]), 2));
	h_p0_km_temp *= NORMALIZED_LENGTH_TO_KM;//转化为km

	//初始时刻static高度为0，按理说会进入下面的块
	if (abs(h_p0_static - h_p0_km_temp) > 1E-10)//如果h_p0_km_static不等于当前的h_p0_km_temp，需要开始/重新查表
	{
		//查表，找参考密度，当然需要在密度表的step为整10倍的前提下
		int integerPart = static_cast<int>(h_p0_km_temp);  // 将double转换为整数部分
		int tens = integerPart / 10;  // 除以10得到十位数
		int reference_value = tens * 10;  // h_p0_km_temp取整为10的倍数后的参考值

		//判断当前初始近地点高度是否在表中
		/*
		 * std::find_if：这是一个算法函数，用于在给定范围内查找满足特定条件的元素。
		atm_list.begin() 和 atm_list.end()：这两个函数分别返回指向atm_list向量的开头和结尾的迭代器。
		 它们定义了搜索的范围，即从atm_list的第一个元素开始，直到最后一个元素之后的位置。
		[h0](const AtmosphereData& data) { return data.h == h0; }：这是一个lambda函数，用于定义搜索的条件。
		 在这个lambda函数中，我们使用[h0]来捕获变量h0，以便在lambda函数体内使用。
		 const AtmosphereData& data指定了lambda函数的参数类型和参数名，这里表示迭代器指向的元素类型为AtmosphereData。
		 函数体return data.h == h0;表示我们将检查每个元素的h成员是否等于h0。
		auto it =：这是一个自动类型推导，auto关键字告诉编译器根据右侧表达式的类型自动推导出it的类型。
		 在这里，它将推导为指向满足条件的元素的迭代器类型。
		综上所述，这行代码的作用是在atm_list向量中查找满足条件data.h == h0的元素。如果找到匹配的元素，它将返回指向第一个匹配元素的迭代器，并将其保存在it中。如果未找到匹配的元素，则it将等于atm_list.end()，表示搜索失败。
		 */
		auto it = std::find_if(atm_list.begin(), atm_list.end(), [reference_value](const AtmosphereData& data)
		{
		  return data.h_km == reference_value;
		});
		int index;
		if (it == atm_list.end())//表格中找不到对应的高度
		{
			if (reference_value < atm_list[0].h_km)//初始近地点小于表格中的最小值
			{
				std::cerr << "Warning: The initial perigee altitude " << h_p0_km_temp
						  << " km is smaller than the minimum value " << atm_list[0].h_km
						  << "km in the table!" << std::endl;
				return -1;
			}
			else//初始近地点大于表格中的最小值
			{
				std::cout << "Warning: The initial perigee altitude " << h_p0_km_temp
						  << " km is larger than the maximum value " << atm_list[atm_list.size() - 1].h_km
						  << "km in the table! Atmosphere drag can be ignored!" << std::endl;
				return 0;//直接返回密度0，即没有大气
			}
		}
		else//能找到，返回索引
			index = std::distance(atm_list.begin(), it);

		//给参考值赋值，以及相应的归一化的值！
		h_p0_static = atm_list[index].h_km;
		rho_p0_static = atm_list[index].rho_kgm3;
		H_p0_static = atm_list[index].H_km;

		h_p0_static_Norm = h_p0_static * KM_TO_NORMALIZED_LENGTH;
		rho_p0_static_Norm = rho_p0_static * UNIT_LENGTH_M * UNIT_LENGTH_M * UNIT_LENGTH_M / UNIT_MASS_KG;//kg/m^3的归一化！
		H_p0_static_Norm = H_p0_static * KM_TO_NORMALIZED_LENGTH;

	}

	//先考虑简单的周日变化(1.222)
	double r;
	Vector3d r_hat, rm_hat;
	r_hat << rv_Norm[0], rv_Norm[1], rv_Norm[2];//是量纲归一化，暂时还没向量归一化
	r = r_hat.norm();//位置的归一化量纲模长
	r_hat.normalize();//向量归一化

	double R_Sun[3];//平均太阳的位置
	Cart_MeanSun(MJD_UT1, R_Sun);
	double lambda_m = 30 * DD2R;//30度的弧度

	rm_hat << R_Sun[0] * cos(lambda_m) - R_Sun[1] * sin(lambda_m),
		R_Sun[1] * cos(lambda_m) + R_Sun[0] * sin(lambda_m),
		R_Sun[2];// 1.222
	rm_hat.normalize(); //再归一化

	double cos_phi_star = r_hat.dot(rm_hat);
	double F_star = 0.3;//松鼠建议取值
	double coeff1 = 1.0 + F_star * cos_phi_star;//周日变化系数(1.220)

	//再考虑高度变化
	//1.223
	double keplerian_Norm[6];//当前实时位置的根数
	CartesianToKeplerian(rv_Norm[0], rv_Norm[1], rv_Norm[2], rv_Norm[3], rv_Norm[4], rv_Norm[5], keplerian_Norm);
	double f = MeanAnomaly2TrueAnomaly(keplerian_Norm[5], keplerian_Norm[1]);
	// phi 是相对于输入赤道平面的地心纬度；phi0 对应初始轨道近地点方向。
	// 两者通过轨道根数直接计算，不调用大地坐标转换。
	double sin_phi = sin(keplerian_Norm[2]) * sin(f + keplerian_Norm[4]);
	double sin_phi0 = sin(initKepArray_norm[2]) * sin(initKepArray_norm[4]);
	double r_p0 = initKepArray_norm[0] * (1 - initKepArray_norm[1]);//a0(1-e0)
	double
		sigma_p0 = r_p0 * (1.0 - constants::Wgs84::FLATTENING * sin_phi * sin_phi) / (1.0 - constants::Wgs84::FLATTENING * sin_phi0 * sin_phi0);
	double mu = 0.1;
	double H_r = H_p0_static_Norm + mu / 2.0 * (r - r_p0);//1.224
	double coeff2 = exp(-(r - sigma_p0) / H_r);//高度变化系数

	//1.220
	rho_Norm = rho_p0_static_Norm * coeff1 * coeff2;//参考密度*周日变化系数*高度变化系数
	return rho_Norm;
}

//重载函数，使用kep类
double Get_AtmosphereDensity(const double* rv_Norm,
	const orbit::Keplerian& initKep_norm,
	double MJD_UT1,
	const std::vector<AtmosphereData>& atm_list)
{
	/*
	 * 计算大气密度，归一化
	 * rho_Norm 归一化的大气密度
	 * r_Norm 归一化的卫星J2000位置
	 * MJD_UT1 计算太阳平轨道用，严格来说（似乎）要用UT1，但是大气的不确定因素太大，UTC也完全没问题
	 * r_p0_list			参考高度		km
	 * density_p0_list		参考密度		kg/m^3
	 * H_p0_list			参考标高		km
	 *
	 */

	double rho_Norm;//最终返回的密度值
	double initKepArray_norm[6];
	initKep_norm.ToArray(initKepArray_norm);
	//static变量仅会被初始化一次，用于存储查表得到的值 km kg/m^3
	static double h_p0_static = 0, rho_p0_static = 0, H_p0_static = 0;
	static double h_p0_static_Norm = 0, rho_p0_static_Norm = 0, H_p0_static_Norm = 0;//对应的归一化值

	double h_p0_km_temp;//当前计算的初始近地点高度km 1.225
	h_p0_km_temp = initKepArray_norm[0] * (1.0 - initKepArray_norm[1])
		- constants::Wgs84::EQUATORIAL_RADIUS / UNIT_LENGTH_M * (1.0 - constants::Wgs84::FLATTENING * pow(sin(initKepArray_norm[2]) * sin(initKepArray_norm[4]), 2));
	h_p0_km_temp *= NORMALIZED_LENGTH_TO_KM;//转化为km

	//初始时刻static高度为0，按理说会进入下面的块
	if (abs(h_p0_static - h_p0_km_temp) > 1E-10)//如果h_p0_km_static不等于当前的h_p0_km_temp，需要开始/重新查表
	{
		//查表，找参考密度，当然需要在密度表的step为整10倍的前提下
		int integerPart = static_cast<int>(h_p0_km_temp);  // 将double转换为整数部分
		int tens = integerPart / 10;  // 除以10得到十位数
		int reference_value = tens * 10;  // h_p0_km_temp取整为10的倍数后的参考值

		//判断当前初始近地点高度是否在表中
		/*
		 * std::find_if：这是一个算法函数，用于在给定范围内查找满足特定条件的元素。
		atm_list.begin() 和 atm_list.end()：这两个函数分别返回指向atm_list向量的开头和结尾的迭代器。
		 它们定义了搜索的范围，即从atm_list的第一个元素开始，直到最后一个元素之后的位置。
		[h0](const AtmosphereData& data) { return data.h == h0; }：这是一个lambda函数，用于定义搜索的条件。
		 在这个lambda函数中，我们使用[h0]来捕获变量h0，以便在lambda函数体内使用。
		 const AtmosphereData& data指定了lambda函数的参数类型和参数名，这里表示迭代器指向的元素类型为AtmosphereData。
		 函数体return data.h == h0;表示我们将检查每个元素的h成员是否等于h0。
		auto it =：这是一个自动类型推导，auto关键字告诉编译器根据右侧表达式的类型自动推导出it的类型。
		 在这里，它将推导为指向满足条件的元素的迭代器类型。
		综上所述，这行代码的作用是在atm_list向量中查找满足条件data.h == h0的元素。如果找到匹配的元素，它将返回指向第一个匹配元素的迭代器，并将其保存在it中。如果未找到匹配的元素，则it将等于atm_list.end()，表示搜索失败。
		 */
		auto it = std::find_if(atm_list.begin(), atm_list.end(), [reference_value](const AtmosphereData& data)
		{
		  return data.h_km == reference_value;
		});
		int index;
		if (it == atm_list.end())//表格中找不到对应的高度
		{
			if (reference_value < atm_list[0].h_km)//初始近地点小于表格中的最小值
			{
				std::cerr << "Warning: The initial perigee altitude " << h_p0_km_temp
						  << " km is smaller than the minimum value " << atm_list[0].h_km
						  << "km in the table!" << std::endl;
				return -1;
			}
			else//初始近地点大于表格中的最小值
			{
				std::cout << "Warning: The initial perigee altitude " << h_p0_km_temp
						  << " km is larger than the maximum value " << atm_list[atm_list.size() - 1].h_km
						  << "km in the table! Atmosphere drag can be ignored!" << std::endl;
				return 0;//直接返回密度0，即没有大气
			}
		}
		else//能找到，返回索引
			index = std::distance(atm_list.begin(), it);

		//给参考值赋值，以及相应的归一化的值！
		h_p0_static = atm_list[index].h_km;
		rho_p0_static = atm_list[index].rho_kgm3;
		H_p0_static = atm_list[index].H_km;

		h_p0_static_Norm = h_p0_static * KM_TO_NORMALIZED_LENGTH;
		rho_p0_static_Norm = rho_p0_static * UNIT_LENGTH_M * UNIT_LENGTH_M * UNIT_LENGTH_M / UNIT_MASS_KG;//kg/m^3的归一化！
		H_p0_static_Norm = H_p0_static * KM_TO_NORMALIZED_LENGTH;

	}

	//FIXME:先作弊直接用老师插值算好的标高差错，记得改掉！！！
	// FIXME:251103发现这里还没改回去，然后加了注释了 有空检查下是否有问题
	//rho_p0_static_Norm = 8.2760764121592779E-17;
	//h_p0_static_Norm = 0.070321166578502958;
	//H_p0_static_Norm = 0.0097621054545199894;
	//先考虑简单的周日变化(1.222)
	double r;
	Vector3d r_hat, rm_hat;
	r_hat << rv_Norm[0], rv_Norm[1], rv_Norm[2];//是量纲归一化，暂时还没向量归一化
	r = r_hat.norm();//位置的归一化量纲模长
	r_hat.normalize();//向量归一化

	double R_Sun[3];//平均太阳的位置
	Cart_MeanSun(MJD_UT1, R_Sun);
	double lambda_m = 30 * DD2R;//30度的弧度

	rm_hat << R_Sun[0] * cos(lambda_m) - R_Sun[1] * sin(lambda_m),
		R_Sun[1] * cos(lambda_m) + R_Sun[0] * sin(lambda_m),
		R_Sun[2];// 1.222
	rm_hat.normalize(); //再归一化

	double cos_phi_star = r_hat.dot(rm_hat);
	//double F_star = 0.3;//松鼠建议取值(m_rhoratio-1.0)/(m_rhoratio+1.0);
	double rhoratio = 1.6;
	double F_star = (rhoratio - 1.0) / (rhoratio + 1.0);//老师modelatm.hpp line154取值，试试
	double coeff1 = 1.0 + F_star * cos_phi_star;//周日变化系数(1.220)

	//再考虑高度变化
	//1.223
	double keplerian_Norm[6];//当前实时位置的根数
	CartesianToKeplerian(rv_Norm[0], rv_Norm[1], rv_Norm[2], rv_Norm[3], rv_Norm[4], rv_Norm[5], keplerian_Norm);
	double f = MeanAnomaly2TrueAnomaly(keplerian_Norm[5], keplerian_Norm[1]);
	// phi 是相对于输入赤道平面的地心纬度；phi0 对应初始轨道近地点方向。
	// 两者通过轨道根数直接计算，不调用大地坐标转换。
	double sin_phi = sin(keplerian_Norm[2]) * sin(f + keplerian_Norm[4]);
	double sin_phi0 = sin(initKepArray_norm[2]) * sin(initKepArray_norm[4]);
	double r_p0 = initKepArray_norm[0] * (1 - initKepArray_norm[1]);//a0(1-e0)
	double
		sigma_p0 = r_p0 * (1.0 - constants::Wgs84::FLATTENING * sin_phi * sin_phi) / (1.0 - constants::Wgs84::FLATTENING * sin_phi0 * sin_phi0);
	double mu = 0.1;
	double H_r = H_p0_static_Norm + mu / 2.0 * (r - r_p0);//1.224
	double coeff2 = exp(-(r - sigma_p0) / H_r);//高度变化系数

	//1.220
	rho_Norm = rho_p0_static_Norm * coeff1 * coeff2;//参考密度*周日变化系数*高度变化系数
	return rho_Norm;
}
