#ifndef __curvature_h
#define __curvature_h

#include "constants.h"
#include "property.h"
#include "approximation.h"
#include "using.h"

namespace J0218
{
	constexpr double C_get_Gamma_e = 1.5 * Rlc * Rlc / e;
	constexpr double C_get_E_cur = 1.5 * h_bar * speed_of_light / Rlc;
	constexpr double C_get_F_single = 1.732 * e * e / (2.0 * PI * h_bar * Rlc);
	const double C_get_F_cur = F0 * (2.0 * PI * Rlc) * Blc * L/
		(speed_of_light * e) * Rlc;

	inline double get_Gamma_e(double theta_x, double z_reduced)
	{
		double res = pow(C_get_Gamma_e * get_E_parallel(theta_x, z_reduced), 0.25);
		return res > 1.0 ? res : 1.0;
	}

	inline double get_E_cur(double theta_x, double z_reduced)
	{
		return C_get_E_cur * pow(get_Gamma_e(theta_x, z_reduced), 3);
	}

	inline double get_F_single(double theta_x, double z_reduced, double E_gamma)
	{
		double E_cur = get_E_cur(theta_x, z_reduced);
		//double temp = get_Gamma_e(theta_x, z_reduced);
		return C_get_F_single * get_Gamma_e(theta_x, z_reduced) *
			approximation(E_gamma / E_cur) / E_gamma;
	}


	struct get_F_single_Z;

	double get_F_cur(double E_gamma) noexcept;

}



/*
 for pulsar B1821
*/
namespace B1821
{
	constexpr double C_get_Gamma_e = 1.5 * Rlc * Rlc / e;
	constexpr double C_get_E_cur = 1.5 * h_bar * speed_of_light / Rlc;
	constexpr double C_get_F_single = 1.732 * e * e / (2.0 * PI * h_bar * Rlc);
	const double C_get_F_cur = F0 * (2.0 * PI * Rlc) * Blc * L/
		(speed_of_light * e) * Rlc;

	inline double get_Gamma_e(double theta_x, double z_reduced)
	{
		double res = pow(C_get_Gamma_e * get_E_parallel(theta_x, z_reduced), 0.25);
		return res > 1.0 ? res : 1.0;
	}

	inline double get_E_cur(double theta_x, double z_reduced)
	{
		return C_get_E_cur * pow(get_Gamma_e(theta_x, z_reduced), 3);
	}

	inline double get_F_single(double theta_x, double z_reduced, double E_gamma)
	{
		double E_cur = get_E_cur(theta_x, z_reduced);
		//double temp = get_Gamma_e(theta_x, z_reduced);
		return C_get_F_single * get_Gamma_e(theta_x, z_reduced) *
			approximation(E_gamma / E_cur) / E_gamma;
	}


	struct get_F_single_Z;

	double get_F_cur(double E_gamma) noexcept;

}


// for pulsar J1939
namespace J1939
{
	constexpr double C_get_Gamma_e = 1.5 * Rlc * Rlc / e;
	constexpr double C_get_E_cur = 1.5 * h_bar * speed_of_light / Rlc;
	constexpr double C_get_F_single = 1.732 * e * e / (2.0 * PI * h_bar * Rlc);
	const double C_get_F_cur = F0 * (2.0 * PI * Rlc) * Blc * L/
		(speed_of_light * e) * Rlc;

	inline double get_Gamma_e(double theta_x, double z_reduced)
	{
		double res = pow(C_get_Gamma_e * get_E_parallel(theta_x, z_reduced), 0.25);
		return res > 1.0 ? res : 1.0;
	}

	inline double get_E_cur(double theta_x, double z_reduced)
	{
		return C_get_E_cur * pow(get_Gamma_e(theta_x, z_reduced), 3);
	}

	inline double get_F_single(double theta_x, double z_reduced, double E_gamma)
	{
		double E_cur = get_E_cur(theta_x, z_reduced);
		//double temp = get_Gamma_e(theta_x, z_reduced);
		return C_get_F_single * get_Gamma_e(theta_x, z_reduced) *
			approximation(E_gamma / E_cur) / E_gamma;
	}


	struct get_F_single_Z;

	double get_F_cur(double E_gamma) noexcept;

}
#endif


