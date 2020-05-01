#pragma once

namespace J0218
{
	double d_P_IC_over_d_Omega(double theta_1, double E_gamma, double gamma);

	double get_IC_single(double x_reduced, double z_reduced, double E_gamma);

	double get_IC(double E_gamma);

}


namespace B1821
{
	double d_P_IC_over_d_Omega(double theta_1, double E_gamma, double gamma);

	double get_IC_single(double x_reduced, double z_reduced, double E_gamma);

	double get_IC(double E_gamma);

}

namespace J1939
{
	double d_P_IC_over_d_Omega(double theta_1, double E_gamma, double gamma);

	double get_IC_single(double x_reduced, double z_reduced, double E_gamma);

	double get_IC(double E_gamma);

}