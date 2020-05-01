#include "curvature.h"
#include "runge.h"
#include "approximation.h"

extern _j0218* j0218;

double get_E_cur(double x_reduced, double z_reduced)
{
	return C_get_E_cur * pow(get_Gamma_e(x_reduced, z_reduced), 3.0);
}

double get_F_single(double x_reduced, double z_reduced, double E_gamma)
{
	double E_cur = get_E_cur(x_reduced, z_reduced);
	return C_get_F_single * get_Gamma_e(x_reduced, z_reduced) * approximation(E_gamma / E_cur) / E_gamma;
}

struct get_F_single_X {
	double z_reduced, E_gamma;

	get_F_single_X() = default;

	void set_Z(double z_reduced)
	{
		this->z_reduced = z_reduced;
	}

	void set_E_gamma(double E_gamma)
	{
		this->E_gamma = E_gamma;
	}

	void set_value(double z_reduced, double E_gamma)
	{
		this->z_reduced = z_reduced;
		this->E_gamma = E_gamma;
	}

	double operator()(double x_reduced) const
	{
		return get_F_single(x_reduced, z_reduced, E_gamma);
	}
} obj_X;

const double border = j0218->f * j0218->ratio;

struct get_F_single_Z {
	double E_gamma;

	get_F_single_Z() {}

	void set_E_gamma(double E_gamma)
	{
		this->E_gamma = E_gamma;
	}

	double operator()(double z_reduced)
	{
		obj_X.set_value(z_reduced, E_gamma);
		if (z_reduced < border) {
			return runge(obj_X, x_m_reduced, 1.0, 1) * (1.0 - j0218->g1);
		}
		else {
			return runge(obj_X, x_m_reduced, 1.0, 1) * (1.0 + j0218->g2);
		}
	}
} obj_Z;

double get_F_cur(double E_gamma)
{
	obj_Z.set_E_gamma(E_gamma);
	return C_get_F_cur * runge(obj_Z, 0.0, j0218->f, 80);
}
