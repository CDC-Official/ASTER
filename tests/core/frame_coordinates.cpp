#include "frame_definition.hpp"
#include "frame_earth.hpp"
#include "orbit_geometry.hpp"
#include "utils_constants.hpp"
#include <sofa.h>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main()
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    try
    {
        aEOP eop;
        eop.dat = 37;
        eop.dut1 = 0.15;
        eop.xp = 0.12 * constants::Conversion::ARCSEC_TO_RAD;
        eop.yp = -0.23 * constants::Conversion::ARCSEC_TO_RAD;
        const double utc = 60000.25;
        Matrix3d c2t;
        GCRS_to_ITRS_Matrix(utc, eop, c2t);
        // 非零三分量位置能检查类接口是否确实使用输入向量，而非未赋值的输出。
        for (bool normalized : {false, true})
        {
            const double scale = normalized ? 1.0 : constants::Scale::UNIT_LENGTH_M;
            Vector3d input(1.2 * scale, -0.7 * scale, 0.4 * scale);
            frame::GCRS celestial(input[0], input[1], input[2], normalized);
            frame::ITRS terrestrial;
            celestial.ToITRS(tracer, utc, eop, terrestrial);
            double output[3];
            terrestrial.ToArray(output);
            Vector3d expected = c2t * input;
            if ((Eigen::Map<Vector3d>(output) - expected).norm() > 1e-12 * scale)
                throw std::runtime_error("Celestial-to-terrestrial class conversion mismatch");
            frame::GCRS recovered;
            terrestrial.ToGCRS(tracer, utc, eop, recovered);
            recovered.ToArray(output);
            if ((Eigen::Map<Vector3d>(output) - input).norm() > 1e-12 * scale)
                throw std::runtime_error("Position round trip mismatch");
        }
        // 独立椭球公式构造 XYZ，检查非零经纬度、高度及不同半球，
        // 避免正反接口同时使用错误单位、却在往返时相互抵消。
        const std::array<std::array<double, 3>, 6> geodetic_cases = {{
            {120.0, 45.0, 1000.0}, {-75.0, -32.0, 12.0},
            {170.0, 80.0, 500.0}, {-170.0, -80.0, 2000.0},
            {35.0, 0.0, -100.0}, {0.0, 20.0, 600000.0}
        }}; // 每行依次为 lambda(deg)、phi(deg)、height(m)。
        for (const std::array<double, 3>& coordinates : geodetic_cases)
        {
            const double lambda = coordinates[0] * constants::Conversion::DEGREE_TO_RAD;
            const double phi = coordinates[1] * constants::Conversion::DEGREE_TO_RAD;
            const double height = coordinates[2];
            // SOFA 按 WGS84 独立计算站址，检查半径与扁率是否配套。
            double xyz[3];
            if (iauGd2gc(WGS84, lambda, phi, height, xyz) != 0)
                throw std::runtime_error("SOFA WGS84 reference conversion failed");
            Vector3d expected(xyz[0], xyz[1], xyz[2]);
            frame::ITRS terrestrial(expected[0], expected[1], expected[2], false);
            terrestrial.Normalize(tracer, constants::Scale::UNIT_LENGTH_M);
            frame::Geodetic recovered;
            terrestrial.ToGeodetic(tracer, recovered);
            double result[3];
            recovered.ToArray(result);
            if (fabs(result[0] * constants::Scale::UNIT_LENGTH_M - height) > 1e-5 ||
                fabs(std::remainder(result[1] - lambda, 2 * constants::PI)) > 1e-12 ||
                fabs(result[2] - phi) > 1e-12)
                throw std::runtime_error("ITRS to geodetic units or coordinates mismatch");
            frame::ITRS round_trip;
            recovered.ToITRS(tracer, round_trip);
            round_trip.Restore(tracer, constants::Scale::UNIT_LENGTH_M);
            round_trip.ToArray(result);
            if ((Eigen::Map<Vector3d>(result) - expected).norm() > 1e-5)
                throw std::runtime_error("Geodetic to ITRS round trip mismatch");
            recovered.Restore(tracer, constants::Scale::UNIT_LENGTH_M);
            recovered.ToArray(result);
            if (fabs(result[0] - height) > 1e-5 ||
                fabs(std::remainder(result[1] - coordinates[0], 360.0)) > 1e-9 ||
                fabs(result[2] - coordinates[1]) > 1e-9)
                throw std::runtime_error("Restored geodetic degree/meter units mismatch");
        }
        // 独立按教材式（1.50）生成地平视线，检查四个方位及一般方向的反算。
        Matrix3d h2c;
        Horizontal_to_GCRS_Matrix(1.1, 0.6, utc, eop.dut1, eop.dat, h2c);
        Matrix3d c2h = h2c.transpose();
        Vector3d station(1.0, 2.0, 3.0);
        for (double azimuth : {0.0, 0.5 * constants::PI, constants::PI, 1.5 * constants::PI, 0.73})
        {
            const double elevation = 0.31, distance = 2.7;
            Vector3d horizontal(cos(elevation) * cos(azimuth),
                                -cos(elevation) * sin(azimuth), sin(elevation));
            Vector3d target = station + distance * h2c * horizontal;
            double result[3];
            RelativePosition_to_Horizontal(c2h, target[0], target[1], target[2],
                                  station[0], station[1], station[2], result);
            if (fabs(result[0] - distance) > 1e-12 ||
                fabs(std::remainder(result[1] - azimuth, 2 * constants::PI)) > 1e-12 ||
                fabs(result[2] - elevation) > 1e-12)
                throw std::runtime_error("Textbook horizontal convention mismatch");
        }
        // 用教材式（1.54）独立计算角度，避免正反矩阵同时写错却仍能往返通过。
        Matrix3d tod_to_gcrs;
        TOD_to_GCRS_Matrix(utc, eop.dat, tod_to_gcrs);
        const double phi = 0.6;
        const double sidereal = iauGst00b(DJM0, utc + eop.dut1 / DAYSEC) + 1.1;
        for (double ra : {0.2, 1.7, 3.1, 5.4})
        {
            const double dec = 0.4;
            Vector3d tod(cos(dec) * cos(ra), cos(dec) * sin(ra), sin(dec));
            Vector3d target = station + tod_to_gcrs * tod;
            const double north = cos(phi) * sin(dec)
                               - sin(phi) * cos(dec) * cos(sidereal - ra);
            const double east = -cos(dec) * sin(sidereal - ra);
            const double up = sin(phi) * sin(dec)
                            + cos(phi) * cos(dec) * cos(sidereal - ra);
            double result[3];
            RelativePosition_to_Horizontal(c2h, target[0], target[1], target[2],
                                           station[0], station[1], station[2], result);
            if (fabs(std::remainder(result[1] - atan2(east, north), 2 * constants::PI)) > 1e-12 ||
                fabs(result[2] - asin(up)) > 1e-12)
                throw std::runtime_error("Textbook equation 1.54 mismatch");
        }
        if (tracer.Code() != bastion::ERR_NOERROR)
            throw std::runtime_error("Unexpected tracer error");
        std::cout << "PASS coordinate frames\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
