#ifndef CDC_ACCEL_DRAG_HPP
#define CDC_ACCEL_DRAG_HPP

int AtmosphereDrag_Acceleration(double* a_celestial,
	const double* rv_celestial,
	double CD,
	double AMR_norm,
	double rho_norm);

#endif // CDC_ACCEL_DRAG_HPP
