/**
 * @file model_propagation.cpp
 * @brief 显示地球与地月轨道预报的动力学模型和积分设置。
 * @details GCRS 与 E-BCRS 使用相同的分区、层级和物理单位显示方式，显示函数不执行轨道传播。
 * 任务历元、初末状态及运行摘要由各自主程序组织。
 */
#include "model_propagation.hpp"
#include <sstream>

/**
 * @brief 显示当前积分器的步长、误差控制与起步设置。
 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
 * @param[in] integrator_id 积分器编号：101 为自适应 RKF78，102 为定步长 RKF78，201/202 为 Adams PE/PECE。
 * @param[in,out] output 接收设置文本的输出流；函数保持该流原有的格式状态。
 * @param[in] precision 普通设置值的有效数字位数；步长和 TOL 固定保留两位小数与两位科学记数法尾数。
 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
 * @details 读取 h0_norm、rkf 和 adams，不修改设置或执行积分。步长乘 UNITTIME 恢复为 s。
 * TOL 比较内部归一化状态的局部截断误差估计；定步长 RKF78 仅将其作为告警阈值。
 * Adams 显示 RKF78 起步、尾段所用的容差及试算次数上限。
 */
void OrbitPropagationSetting::Display(bastion::Overwatch& tracer, int integrator_id, std::ostream& output, int precision, int label_width) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 只展开本次生效的积分设置；列宽由调用程序与其余任务信息统一。
    label_width = std::max(label_width, 31);
    const char* name = integrator_id == 101 || integrator_id == 102 ? "RKF78"
                     : integrator_id == 201 ? "Adams PE" : "Adams PECE";
    std::ostringstream text;
    text << std::defaultfloat << std::setprecision(precision) << "\n================================= Integrator settings ==================================\n\n"
         << std::left << std::setw(label_width) << "- Integrator" << ": " << integrator_id << " — " << name
         << '\n' << std::setw(label_width) << "  > Step mode" << ": " << (integrator_id == 101 ? "Adaptive" : "Fixed")
         << '\n' << std::setw(label_width) << (integrator_id == 101 ? "  > Initial step [s]" : "  > Step [s]")
         << ": " << std::fixed << std::setprecision(2) << std::fabs(h0_norm * UNITTIME);

    // Adams 的 TOL 用于 RKF78 起步及尾段；普通多步公式没有相应的单步误差估计。
    if (integrator_id == 201 || integrator_id == 202)
    {
        text << '\n' << std::setw(label_width) << "  > History nodes" << ": " << adams.steps
             << '\n' << std::setw(label_width) << "  > Correction" << ": " << (adams.use_correction ? "One correction per step" : "None")
             << '\n' << std::setw(label_width) << "  > Startup / Tail method" << ": RKF78"
             << '\n' << std::setw(label_width) << "    > TOL" << ": " << std::scientific << std::setprecision(2) << rkf.maxTE_norm
             << '\n' << std::setw(label_width) << "    > Maximum trials" << ": " << rkf.iterCount_max;
    }
    else if (integrator_id == 102)
    {
        text << '\n' << std::setw(label_width) << "  > TOL (warning threshold)" << ": " << std::scientific << std::setprecision(2) << rkf.maxTE_norm;
    }
    else
    {
        text << '\n' << std::setw(label_width) << "  > TOL" << ": " << std::scientific << std::setprecision(2) << rkf.maxTE_norm
             << '\n' << std::setw(label_width) << "  > Maximum trials" << ": " << rkf.iterCount_max;
    }

    output << text.str() << '\n';

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing integrator settings.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 显示实际启用的引力、阻力和光压模型及其物理参数。
 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
 * @param[in,out] output 接收模型文本的输出流；函数保持该流原有的格式状态。
 * @param[in] precision 面积、AMR、系数和大气表头条件等普通数值的有效数字位数。
 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
 * @details 读取地月非球形开关及截断、日月历表模式、阻力和光压开关、sat 和密度表条件。
 * 面积乘 ae² 恢复为 m²，AMR 乘 ae²/MASS_OF_EARTH 恢复为 m²/kg。
 * 月球引力截断为 0 / 0 时注明仅考虑质点引力；缺少密度表头条件时显示 Not provided。
 * 光压使用 Cannonball 模型，当前加速度实现中地影衰减关闭。函数不修改动力学参数。
 */
void OrbitPropagationNModel::Display(bastion::Overwatch& tracer, std::ostream& output, int precision, int label_width) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 逐项显示模型开关，生效时再展开参数；资料路径由主程序的读取事件记录。
    label_width = std::max(label_width, 31);
    std::ostringstream text;
    text << std::defaultfloat << std::setprecision(precision) << std::left
         << "\n==================================== Dynamics model ====================================\n\n"
         << std::setw(label_width) << "- Earth non-spherical gravity" << ": " << (use_earth_NSG ? "ON" : "OFF");

    if (use_earth_NSG)
    {
        text << '\n' << std::setw(label_width) << "  > Gravity model" << ": GGM03C"
             << '\n' << std::setw(label_width) << "  > Zonal / Tesseral" << ": " << earth_zonal_order << " / " << earth_tesseral_degree;
    }

    // 月球质点与非球形项共同归入月球引力；仅启用质点项时截断为 0 / 0。
    text << "\n\n" << std::setw(label_width) << "- Lunar gravity" << ": " << (lunar_gravity_config ? "ON" : "OFF");

    if (lunar_gravity_config)
    {
        text << '\n' << std::setw(label_width) << "  > Zonal / Tesseral" << ": " << (use_lunar_NSG ? lunar_zonal_order : 0)
             << " / " << (use_lunar_NSG ? lunar_tesseral_degree : 0);

        if (!use_lunar_NSG || (lunar_zonal_order == 0 && lunar_tesseral_degree == 0)) text << " (point-mass only)";

        text << '\n' << std::setw(label_width) << "  > Ephemeris" << ": " << (lunar_gravity_config == LUNAR_GRAVITY_JPL ? "JPL" : "Analytical");
    }

    text << "\n\n" << std::setw(label_width) << "- Solar point-mass gravity" << ": " << (solar_gravity_config ? "ON" : "OFF");

    if (solar_gravity_config)
    {
        text << '\n' << std::setw(label_width) << "  > Ephemeris" << ": " << (solar_gravity_config == SOLAR_GRAVITY_JPL ? "JPL" : "Analytical");
    }

    // Ap、F10.7 是密度表提供的条件；未提供时分别注明，不能补成固定默认值。
    text << "\n\n" << std::setw(label_width) << "- Atmospheric drag" << ": " << (use_drag ? "ON" : "OFF");

    if (use_drag)
    {
        text << '\n' << std::setw(label_width) << "  > Ap / F10.7" << ": ";

        if (std::isfinite(atmosphere_Ap))
        {
            text << atmosphere_Ap;
        }
        else
        {
            text << "Not provided";
        }

        text << " / ";

        if (std::isfinite(atmosphere_F107))
        {
            text << atmosphere_F107;
        }
        else
        {
            text << "Not provided";
        }

        // UTF-8 的上标 ² 占两个字节、一个显示列，setw 相应多留一个字节。
        text << '\n' << std::setw(label_width + 1) << "  > Cross-section [m²]" << ": " << sat.area_drag * ae * ae
             << '\n' << std::setw(label_width + 1) << "  > AMR [m²/kg]" << ": " << sat.AMR_drag * ae * ae / MASS_OF_EARTH
             << '\n' << std::setw(label_width) << "  > Cd" << ": " << sat.Cd;
    }

    // 系数直接显示当前反射参数；地影开关对应实际加速度实现。
    text << "\n\n" << std::setw(label_width) << "- SRP (Cannonball)" << ": " << (use_SRP ? "ON" : "OFF");

    if (use_SRP)
    {
        text << '\n' << std::setw(label_width + 1) << "  > Cross-section [m²]" << ": " << sat.area_SRP * ae * ae
             << '\n' << std::setw(label_width + 1) << "  > AMR [m²/kg]" << ": " << sat.AMR_SRP * ae * ae / MASS_OF_EARTH
             << '\n' << std::setw(label_width) << "  > Reflectivity coefficient" << ": " << sat.eta_SRP
             << '\n' << std::setw(label_width) << "  > Shadow attenuation" << ": OFF";
    }

    output << text.str() << '\n';

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing dynamics model.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 显示 E-BCRS/TDB 积分器的步长、误差控制与起步设置。
 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
 * @param[in] integrator_id 积分器编号：101 为自适应 RKF78，102 为定步长 RKF78，201/202 为 Adams PE/PECE。
 * @param[in,out] output 接收设置文本的输出流；函数保持该流原有的格式状态。
 * @param[in] precision 普通设置值的有效数字位数；步长和 TOL 固定保留两位小数与两位科学记数法尾数。
 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
 * @details 读取 h0_norm、rkf 和 adams，不修改设置或执行积分。步长乘 UNITTIME 恢复为 s，积分历元采用 TDB。
 * TOL 比较内部归一化状态的局部截断误差估计；定步长 RKF78 仅将其作为告警阈值。
 * Adams 显示 RKF78 起步、尾段所用的容差及试算次数上限。
 */
void OrbitPropagationSetting_EBCRS::Display(bastion::Overwatch& tracer, int integrator_id, std::ostream& output, int precision, int label_width) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 只展开本次生效的积分设置；列宽由调用程序与其余任务信息统一。
    label_width = std::max(label_width, 31);
    const char* name = integrator_id == 101 || integrator_id == 102 ? "RKF78"
                     : integrator_id == 201 ? "Adams PE" : "Adams PECE";
    std::ostringstream text;
    text << std::defaultfloat << std::setprecision(precision) << "\n================================= Integrator settings ==================================\n\n"
         << std::left << std::setw(label_width) << "- Integrator" << ": " << integrator_id << " — " << name
         << '\n' << std::setw(label_width) << "  > Step mode" << ": " << (integrator_id == 101 ? "Adaptive" : "Fixed")
         << '\n' << std::setw(label_width) << (integrator_id == 101 ? "  > Initial step [s]" : "  > Step [s]")
         << ": " << std::fixed << std::setprecision(2) << std::fabs(h0_norm * UNITTIME);

    // Adams 的 TOL 用于 RKF78 起步及尾段；普通多步公式没有相应的单步误差估计。
    if (integrator_id == 201 || integrator_id == 202)
    {
        text << '\n' << std::setw(label_width) << "  > History nodes" << ": " << adams.steps
             << '\n' << std::setw(label_width) << "  > Correction" << ": " << (adams.use_correction ? "One correction per step" : "None")
             << '\n' << std::setw(label_width) << "  > Startup / Tail method" << ": RKF78"
             << '\n' << std::setw(label_width) << "    > TOL" << ": " << std::scientific << std::setprecision(2) << rkf.maxTE_norm
             << '\n' << std::setw(label_width) << "    > Maximum trials" << ": " << rkf.iterCount_max;
    }
    else if (integrator_id == 102)
    {
        text << '\n' << std::setw(label_width) << "  > TOL (warning threshold)" << ": " << std::scientific << std::setprecision(2) << rkf.maxTE_norm;
    }
    else
    {
        text << '\n' << std::setw(label_width) << "  > TOL" << ": " << std::scientific << std::setprecision(2) << rkf.maxTE_norm
             << '\n' << std::setw(label_width) << "  > Maximum trials" << ": " << rkf.iterCount_max;
    }

    output << text.str() << '\n';

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing integrator settings.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    tracer.CheckOut(__func__);
}

/**
 * @brief 显示 E-BCRS 实际启用的引力、阻力、光压及后牛顿修正。
 * @param[in,out] tracer 记录显示过程和输出错误；已有错误时直接返回。
 * @param[in,out] output 接收模型文本的输出流；函数保持该流原有的格式状态。
 * @param[in] precision 面积、AMR、系数和大气表头条件等普通数值的有效数字位数。
 * @param[in] label_width 冒号前的字段显示列宽，包含层级符号和缩进，最小采用 31。
 * @details 读取地月非球形开关及截断、日月与行星引力开关、行星 JPL 编号、阻力、光压、后牛顿开关、sat 和密度表条件。
 * 面积乘 ae² 恢复为 m²，AMR 乘 ae²/MASS_OF_EARTH 恢复为 m²/kg。
 * 月球引力截断为 0 / 0 时注明仅考虑质点引力；缺少密度表头条件时显示 Not provided。
 * 日月与行星使用 JPL；摄动天体按 JPL 编号及英文名称显示，非法编号显示 Unknown。
 * 后牛顿开关对应当前 EIH 修正，使用 11 个 JPL 大天体，独立于经典行星摄动列表。
 * 光压采用 Cannonball 模型，地影衰减关闭；函数不修改动力学参数或转换轨道坐标。
 */
void OrbitPropagationNModel_EBCRS::Display(bastion::Overwatch& tracer, std::ostream& output, int precision, int label_width) const
{
    if (tracer.Code() != bastion::ERR_NOERROR) return;

    tracer.CheckIn(__FILE__, __func__, __LINE__);

    // 逐项显示模型开关，生效时再展开参数；资料路径由主程序的读取事件记录。
    label_width = std::max(label_width, 31);
    std::ostringstream text;
    text << std::defaultfloat << std::setprecision(precision) << std::left
         << "\n==================================== Dynamics model ====================================\n\n"
         << std::setw(label_width) << "- Earth non-spherical gravity" << ": " << (use_earth_NSG ? "ON" : "OFF");

    if (use_earth_NSG)
    {
        text << '\n' << std::setw(label_width) << "  > Gravity model" << ": GGM03C"
             << '\n' << std::setw(label_width) << "  > Zonal / Tesseral" << ": " << earth_zonal_order << " / " << earth_tesseral_degree;
    }

    // 月球质点与非球形项共同归入月球引力；仅启用质点项时截断为 0 / 0。
    text << "\n\n" << std::setw(label_width) << "- Lunar gravity" << ": " << (lunar_gravity_config ? "ON" : "OFF");

    if (lunar_gravity_config)
    {
        text << '\n' << std::setw(label_width) << "  > Zonal / Tesseral" << ": " << (use_lunar_NSG ? lunar_zonal_order : 0)
             << " / " << (use_lunar_NSG ? lunar_tesseral_degree : 0);

        if (!use_lunar_NSG || (lunar_zonal_order == 0 && lunar_tesseral_degree == 0)) text << " (point-mass only)";

        text << '\n' << std::setw(label_width) << "  > Ephemeris" << ": JPL";
    }

    text << "\n\n" << std::setw(label_width) << "- Solar point-mass gravity" << ": " << (solar_gravity_config ? "ON" : "OFF");

    if (solar_gravity_config)
    {
        text << '\n' << std::setw(label_width) << "  > Ephemeris" << ": JPL";
    }

    // 经典行星摄动按输入的 JPL 编号逐项列名，地球、月球和太阳由各自的引力项处理。
    text << "\n\n" << std::setw(label_width) << "- Planetary point-mass gravity" << ": " << (use_planetary_gravity ? "ON" : "OFF");

    if (use_planetary_gravity)
    {
        const std::array<const char*, 10> body_names = {
            "Unknown", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
        text << '\n' << std::setw(label_width) << "  > Perturbing bodies" << ": ";

        for (std::size_t i = 0; i < planetary_body_id.size(); ++i)
        {
            if (i > 0) text << ' ';

            const int id = planetary_body_id[i];
            const char* name = id > 0 && id < static_cast<int>(body_names.size()) ? body_names[id] : "Unknown";
            text << id << " (" << name << ")";
        }

        text << '\n' << std::setw(label_width) << "  > Ephemeris" << ": JPL";
    }

    // Ap、F10.7 是密度表提供的条件；未提供时分别注明，不能补成固定默认值。
    text << "\n\n" << std::setw(label_width) << "- Atmospheric drag" << ": " << (use_drag ? "ON" : "OFF");

    if (use_drag)
    {
        text << '\n' << std::setw(label_width) << "  > Ap / F10.7" << ": ";

        if (std::isfinite(atmosphere_Ap))
        {
            text << atmosphere_Ap;
        }
        else
        {
            text << "Not provided";
        }

        text << " / ";

        if (std::isfinite(atmosphere_F107))
        {
            text << atmosphere_F107;
        }
        else
        {
            text << "Not provided";
        }

        // UTF-8 的上标 ² 占两个字节、一个显示列，setw 相应多留一个字节。
        text << '\n' << std::setw(label_width + 1) << "  > Cross-section [m²]" << ": " << sat.area_drag * ae * ae
             << '\n' << std::setw(label_width + 1) << "  > AMR [m²/kg]" << ": " << sat.AMR_drag * ae * ae / MASS_OF_EARTH
             << '\n' << std::setw(label_width) << "  > Cd" << ": " << sat.Cd;
    }

    // 系数直接显示当前反射参数；地影开关对应实际加速度实现。
    text << "\n\n" << std::setw(label_width) << "- SRP (Cannonball)" << ": " << (use_SRP ? "ON" : "OFF");

    if (use_SRP)
    {
        text << '\n' << std::setw(label_width + 1) << "  > Cross-section [m²]" << ": " << sat.area_SRP * ae * ae
             << '\n' << std::setw(label_width + 1) << "  > AMR [m²/kg]" << ": " << sat.AMR_SRP * ae * ae / MASS_OF_EARTH
             << '\n' << std::setw(label_width) << "  > Reflectivity coefficient" << ": " << sat.eta_SRP
             << '\n' << std::setw(label_width) << "  > Shadow attenuation" << ": OFF";
    }

    // 当前 EIH 实现只有启用开关，具体天体求和与地心相对加速度在动力学中计算。
    text << "\n\n" << std::setw(label_width) << "- Post-Newtonian correction" << ": " << (use_post_newtonian ? "ON" : "OFF");
    output << text.str() << '\n';

    if (!output)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed while writing dynamics model.", bastion::ERROR);
        tracer.CheckOut(__func__);
        return;
    }

    tracer.CheckOut(__func__);
}
