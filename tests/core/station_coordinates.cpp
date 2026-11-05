/**
 * @file station_coordinates.cpp
 * @brief 用 SOFA 参考、六向转换和往返结果检查测站位置接口。
 * @details 使用初轨算例的测站与观测时刻，并覆盖赤道、两极和近极点。
 * 可选的参考输出供文件入口核对，保存 UTC MJD、ITRS XYZ(m)、GCRS XYZ(m)。
 */
#include "frame_definition.hpp"
#include "utils_constants.hpp"
#include <sofa.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

using namespace bastion;
using constants::Scale::UNIT_LENGTH_M;

int main(int argc, char* argv[])
{
    Overwatch tracer(__FILE__, __func__, __LINE__);

    if (argc != 2 && argc != 3)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Expected EOP path and optional reference output path.", EventType::ERROR);
        tracer.EmergencyExit();
    }

    // 前三行取自 GL 算例的测站和前三次观测，其他位置覆盖椭球反算的典型情况。
    const std::array<std::array<double, 4>, 9> cases = {{
        {58406.8859375000, 86.174633, 41.725892, 943.0},
        {58406.8859722222, 86.174633, 41.725892, 943.0},
        {58406.8860069444, 86.174633, 41.725892, 943.0},
        {60000.25, 120.0, 45.0, 1000.0},
        {60000.25, -75.0, -32.0, 12.0},
        {60000.25, 35.0, 0.0, -100.0},
        {60000.25, 0.0, 90.0, 1000.0},
        {60000.25, 0.0, -90.0, 2000.0},
        {55255.99923398, 170.0, 89.999999999, 500.0}
    }};

    EOP table;
    std::string eop_path = argv[1];
    table.Read_EOP(tracer, eop_path);
    if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

    std::ostringstream reference;
    reference << std::setprecision(17);
    double max_forward_m = 0.0, max_rotation_m = 0.0, max_round_trip_m = 0.0, max_height_m = 0.0;

    for (const std::array<double, 4>& record : cases)
    {
        // 使用现有 EOP 插值，与 GL、GLML 的逐观测站址转换使用相同资料。
        aEOP eop;
        table.Get_EOP(tracer, record[0], eop, 1);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        // SOFA 以米直接计算 WGS84 地固位置，独立检查类接口的椭球和归一化。
        double expected_itrs[3];
        const int status = iauGd2gc(WGS84, record[1] * constants::Conversion::DEGREE_TO_RAD,
            record[2] * constants::Conversion::DEGREE_TO_RAD, record[3], expected_itrs);

        if (status != 0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "SOFA WGS84 reference failed.", EventType::ERROR);
            tracer.EmergencyExit();
        }

        frame::Geodetic geodetic(record[3], record[1], record[2], false);
        geodetic.Normalize(tracer, UNIT_LENGTH_M);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        frame::ITRS itrs;
        geodetic.ToITRS(tracer, itrs);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        double itrs_norm[3];
        itrs.ToArray(itrs_norm);

        for (int axis = 0; axis < 3; ++axis)
        {
            max_forward_m = std::max(max_forward_m, std::fabs(itrs_norm[axis] * UNIT_LENGTH_M - expected_itrs[axis]));
        }

        // 独立用 SOFA 的合成接口组装天球到地固矩阵，检查既有旋转的方向和 EOP 单位。
        const double tt = (record[0] + eop.dat / constants::SECONDS_PER_DAY) + 32.184 / constants::SECONDS_PER_DAY;
        const double ut1 = record[0] + eop.dut1 / constants::SECONDS_PER_DAY;
        double x, y, s, c2i[3][3], polar_motion[3][3], c2t[3][3], expected_gcrs[3];
        iauXys00b(constants::MJD_OFFSET, tt, &x, &y, &s);
        iauC2ixys(x + eop.dx00, y + eop.dy00, s, c2i);
        iauPom00(eop.xp, eop.yp, iauSp00(constants::MJD_OFFSET, tt), polar_motion);
        iauC2tcio(c2i, iauEra00(constants::MJD_OFFSET, ut1), polar_motion, c2t);
        iauTrxp(c2t, expected_itrs, expected_gcrs);

        // 先走初轨主程序的分步调用，再与大地坐标的一步接口逐分量比较。
        frame::GCRS gcrs, direct_gcrs;
        itrs.ToGCRS(tracer, record[0], eop, gcrs);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        geodetic.ToGCRS(tracer, record[0], eop, direct_gcrs);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        double gcrs_norm[3], direct_norm[3];
        gcrs.ToArray(gcrs_norm);
        direct_gcrs.ToArray(direct_norm);

        for (int axis = 0; axis < 3; ++axis)
        {
            max_rotation_m = std::max(max_rotation_m, std::fabs(gcrs_norm[axis] * UNIT_LENGTH_M - expected_gcrs[axis]));

            if (direct_norm[axis] != gcrs_norm[axis])
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Direct and IOD two-step station positions differ.", EventType::ERROR);
                tracer.EmergencyExit();
            }
        }

        // 检查反向旋转，并分别从 ITRS、GCRS 反算经纬高。
        frame::ITRS recovered_itrs;
        gcrs.ToITRS(tracer, record[0], eop, recovered_itrs);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        double recovered_norm[3];
        recovered_itrs.ToArray(recovered_norm);

        for (int axis = 0; axis < 3; ++axis)
        {
            max_round_trip_m = std::max(max_round_trip_m, std::fabs(recovered_norm[axis] - itrs_norm[axis]) * UNIT_LENGTH_M);
        }

        for (bool from_gcrs : {false, true})
        {
            frame::Geodetic recovered;

            if (from_gcrs)
            {
                gcrs.ToGeodetic(tracer, record[0], eop, recovered);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
            }
            else
            {
                itrs.ToGeodetic(tracer, recovered);
                if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
            }

            // 往返检查用地固位置差，极点处的经度自由度不参与位置评价。
            frame::ITRS round_trip;
            recovered.ToITRS(tracer, round_trip);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

            round_trip.ToArray(recovered_norm);

            for (int axis = 0; axis < 3; ++axis)
            {
                max_round_trip_m = std::max(max_round_trip_m, std::fabs(recovered_norm[axis] - itrs_norm[axis]) * UNIT_LENGTH_M);
            }

            recovered.Restore(tracer, UNIT_LENGTH_M);
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

            double values[3];
            recovered.ToArray(values);
            max_height_m = std::max(max_height_m, std::fabs(values[0] - record[3]));

            if (std::fabs(values[2] - record[2]) > 1e-10
                || (std::fabs(record[2]) < 89.999999 && std::fabs(std::remainder(values[1] - record[1], 360.0)) > 1e-10))
            {
                tracer.LogEvent(__FILE__, __func__, __LINE__, "Recovered WGS84 longitude or latitude differs from input.", EventType::ERROR);
                tracer.EmergencyExit();
            }
        }

        // 内存接口的米制结果供批量程序核对，前三行直接对应 Gooding 的三次观测。
        reference << record[0];

        for (double value : itrs_norm) reference << ' ' << value * UNIT_LENGTH_M;

        for (double value : gcrs_norm) reference << ' ' << value * UNIT_LENGTH_M;

        reference << '\n';
    }

    if (max_forward_m > 1e-6 || max_rotation_m > 1e-6 || max_round_trip_m > 1e-6 || max_height_m > 1e-6)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, "Station position or height check exceeds 1 micrometre.", EventType::ERROR);
        tracer.EmergencyExit();
    }

    // 精确自转轴位置单独检查，北极和南极都应恢复给定椭球高。
    for (double sign : {-1.0, 1.0})
    {
        const double height_m = 123.0;
        const double polar_radius_m = constants::Wgs84::EQUATORIAL_RADIUS * (1.0 - constants::Wgs84::FLATTENING);
        frame::ITRS pole(0.0, 0.0, sign * (polar_radius_m + height_m) / UNIT_LENGTH_M, true);
        frame::Geodetic recovered;
        pole.ToGeodetic(tracer, recovered);
        if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();

        double result[3];
        recovered.ToArray(result);

        if (std::fabs(result[0] * UNIT_LENGTH_M - height_m) > 1e-6 || result[1] != 0.0
            || result[2] != sign * constants::PI / 2.0)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Exact pole conversion failed.", EventType::ERROR);
            tracer.EmergencyExit();
        }
    }

    // 组合转换遇到非法经纬度或 EOP 时保留原输出，并保留底层错误的调用链。
    for (bool invalid_eop : {false, true})
    {
        Overwatch invalid_tracer(__FILE__, __func__, __LINE__);
        aEOP eop;
        eop.dat = 37.0;
        frame::Geodetic input(0.0, 0.0, invalid_eop ? 0.0 : 91.0, false);
        input.Normalize(invalid_tracer, UNIT_LENGTH_M);
        if (invalid_tracer.Code() != ERR_NOERROR) invalid_tracer.EmergencyExit();

        if (invalid_eop) eop.xp = std::numeric_limits<double>::quiet_NaN();

        frame::GCRS unchanged(1.0, 2.0, 3.0, true);
        input.ToGCRS(invalid_tracer, 60000.25, eop, unchanged);
        double result[3];
        unchanged.ToArray(result);
        std::ostringstream trace;
        invalid_tracer.ReportTrace(trace);

        if (invalid_tracer.Code() == ERR_NOERROR || result[0] != 1.0 || result[1] != 2.0 || result[2] != 3.0
            || trace.str().find("ToGCRS") == std::string::npos)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Failed station conversion changed output or lost its trace.", EventType::ERROR);
            tracer.EmergencyExit();
        }
    }

    if (argc == 3)
    {
        std::ofstream output(argv[2]);
        output << reference.str();
        output.close();

        if (!output)
        {
            tracer.LogEvent(__FILE__, __func__, __LINE__, "Cannot save station reference output.", EventType::ERROR);
            tracer.EmergencyExit();
        }
    }

    std::cout << std::setprecision(17) << "PASS six station coordinate directions\n"
              << "SOFA WGS84 max component difference [m]: " << max_forward_m << '\n'
              << "SOFA IAU2000B max component difference [m]: " << max_rotation_m << '\n'
              << "Round trip max component difference [m]: " << max_round_trip_m << '\n'
              << "Recovered height max difference [m]: " << max_height_m << '\n';
    tracer.CheckOut(__func__);
    return 0;
}
