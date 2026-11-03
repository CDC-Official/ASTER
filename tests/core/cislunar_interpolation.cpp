/**
 * @file cislunar_interpolation.cpp
 * @brief 用已知多项式核对插值和窗口边界，并核对 JPL 拆分日期的亚微秒响应。
 */
#include "model_cislunar_observation.hpp"
#include "model_ephemeris_jpl.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv)
{
    try
    {
        // 三次轨道与线性钟差在不同采样窗口内应恢复解析值。
        for (bool is_clock : {false, true})
        {
            observation::CislunarSeries series;
            series.is_clock = is_clock;
            series.values.resize(50, 2);

            for (int i = 0; i < 50; ++i)
            {
                series.days.push_back(i * 60.0 / 86400.0);
                const double t = i * 60.0 / 1000.0;
                series.values(i, 0) = 3e8 + 1e6 * t - 2e5 * t * t + 1e4 * t * t * t;
                series.values(i, 1) = 1e-4 + i * 60.0 * 1e-10;
            }

            for (double tag_s : {5.0, 540.0, 545.0, 1045.0, 2930.0})
            {
                const auto polynomial = observation::PrepareCislunarInterpolant(series, tag_s / 86400.0);
                const double event_s = tag_s - 1.3;
                const auto value = observation::EvaluateCislunarInterpolant(polynomial, event_s / 86400.0);
                const double t = event_s / 1000.0;
                const double expected = 3e8 + 1e6 * t - 2e5 * t * t + 1e4 * t * t * t;

                if (std::abs(value(0) - expected) > 2e-5 || std::abs(value(1) - (1e-4 + event_s * 1e-10)) > 1e-17)
                    throw std::runtime_error("Polynomial interpolation failed across a window boundary.");
            }
        }

        // 1 微秒内地球约移动 3 厘米；拆分日期应保留这一变化。
        if (argc != 2) throw std::runtime_error("Expected JPL ephemeris path.");

        const int loaded = Read_JPLEPH(argv[1]);
        if (loaded < 0) throw std::runtime_error("Cannot load JPL test ephemeris.");

        double first[6], next[6];
        pleph(2462276.5, 0.25, 3, 12, first);
        pleph(2462276.5, 0.25 + 1e-6 / 86400.0, 3, 12, next);

        for (int i = 0; i < 3; ++i)
        {
            const double difference_m = (next[i] - first[i]) * 1000.0;
            const double expected_m = first[i + 3] * 1e-3;

            if (std::abs(difference_m - expected_m) > 1e-4)
                throw std::runtime_error("JPL split date lost the one-microsecond displacement.");
        }

        std::cout << "Polynomial interpolation and split-date ephemeris checks passed.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
