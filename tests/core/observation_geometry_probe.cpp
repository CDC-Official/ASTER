// 为独立几何测试提供同历元的测站、目标和太阳向量，不调用可见性判据。
#include "frame_definition.hpp"
#include "math_basic.hpp"
#include "orbit_statedef.hpp"
#include "model_ephemeris_analytic.hpp"
#include "model_ephemeris_jpl.hpp"
#include "utils_constants.hpp"
#include <sofa.h>
#include <fstream>
#include <iomanip>
#include <stdexcept>

using orbit::Keplerian;

int main(int argc, char** argv)
{
    bastion::Overwatch tracer(__FILE__, __func__, __LINE__);
    try
    {
        if (argc < 4) throw std::runtime_error("Arguments: EOP, Keplerian fixture, Sun mode, [JPL]");
        std::string eop_path = argv[1];
        const int mode = std::stoi(argv[3]);
        if (mode == 1 && (argc != 5 || Read_JPLEPH(argv[4]) < 0)) throw std::runtime_error("JPL fixture unavailable");
        EOP table;
        table.Read_EOP(tracer, eop_path);
        frame::Geodetic site;
        site.Set(943.0, 86.174633, 41.725892, false);
        site.Normalize(tracer, ae);
        std::ifstream input(argv[2]);
        if (!input) throw std::runtime_error("Cannot read Keplerian fixture");
        std::cout << std::setprecision(17) << ae << '\n';
        int id;
        double utc, elements[6];
        while (input >> id >> utc >> elements[0] >> elements[1] >> elements[2] >> elements[3] >> elements[4] >> elements[5])
        {
            aEOP current;
            table.Get_EOP(tracer, utc, current, 1);
            frame::GCRS station;
            site.ToGCRS(tracer, utc, current, station);
            Keplerian state(elements[0], elements[1], elements[2], elements[3], elements[4], elements[5], false);
            state.Normalize(tracer, ae / 1000.0);
            double cart[6], sun[6] = {};
            state.ToCartesian(tracer, cart);
            const double ut1 = utc + current.dut1 / 86400.0;
            if (mode == 2) Cart_MeanSun(ut1, sun);
            else
            {
                // 独立走 SOFA UTC->TAI->TT 链；算例中 TAI-UTC 与 EOP 表相同。
                double tai1, tai2, tt1, tt2, tdb1, tdb2;
                if (iauUtctai(START_OF_MJD, utc, &tai1, &tai2) < 0) throw std::runtime_error("Invalid UTC");
                iauTaitt(tai1, tai2, &tt1, &tt2);
                iauTttdb(tt1, tt2, iauDtdb(tt1, tt2, ut1 - floor(ut1), 0, 0, 0), &tdb1, &tdb2);
                if (Get_JPLEPH(tdb1 - START_OF_MJD + tdb2, 11, 3, sun) < 0) throw std::runtime_error("JPL query failed");
            }
            if (tracer.Code() != ERR_NOERROR) tracer.EmergencyExit();
            std::cout << id << ' ' << utc << ' '
                      << station.GetValue(tracer, Identifier::VAR_CARTX) << ' ' << station.GetValue(tracer, Identifier::VAR_CARTY) << ' ' << station.GetValue(tracer, Identifier::VAR_CARTZ);
            for (int j = 0; j < 3; ++j) std::cout << ' ' << cart[j];
            for (int j = 0; j < 3; ++j) std::cout << ' ' << sun[j];
            std::cout << '\n';
        }
        tracer.CheckOut(__func__);
        return 0;
    }
    catch (const std::exception& error)
    {
        tracer.LogEvent(__FILE__, __func__, __LINE__, error.what(), bastion::EventType::ERROR, FORCE_ALERT_ON);
        tracer.EmergencyExit();
    }
    return 1;
}
