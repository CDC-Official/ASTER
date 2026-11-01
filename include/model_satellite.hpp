//
// Created by 陈酊存 on 25-5-29.
//

#ifndef CDC_LIB_INCLUDE_MODEL_SATELLITE_HPP_
#define CDC_LIB_INCLUDE_MODEL_SATELLITE_HPP_

#include <string>

namespace model
{
	namespace common
	{
		struct SphereSatellite
		{
			double mass;
			double area_drag;        // cross-sectional area (for drag)
			double AMR_drag = 1e9;   // area-to-mass ratio (for drag), separatedly defined to facilitate frequent use

			double area_SRP;         // cross-sectional area (for SRP)
			double AMR_SRP = 1e9;   // area-to-mass ratio (for SRP)
			double eta_SRP = 0.1;

			//double reflectivity = 1.44;   // general reflectivity + 1, 0~2, specular or diffuse reflection not considered
			double Cd = 2.2;        // drag coef is a physical constant, yet we link it to the satellite to facilitate possible "correction"
			int Code = 0;                           // code and id are both supposed to be the unique ID of a satellite,
			std::string Id = "?";   // they are to be used at the user's disposal wherever integer or string is more convenient
		};

	}
}

#endif //CDC_LIB_INCLUDE_MODEL_SATELLITE_HPP_
