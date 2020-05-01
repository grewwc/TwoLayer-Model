#include "curvature.h"
#include "runge.h"
#include "using.h"

#include <iostream>
extern _j0218* j0218;
extern _b1821* b1821;
extern _j1939* j1939;



namespace J0218
{
	struct get_F_single_X
	{
		double z_reduced, E_gamma;

		get_F_single_X() = default;

		void set_Z(double mZ_reduced) noexcept
		{
			this->z_reduced = mZ_reduced;
		}

		void set_E_gamma(double mE_gamma) noexcept
		{
			this->E_gamma = mE_gamma;
		}

		void set_value(double mZ_reduced, double mE_gamma) noexcept
		{
			this->z_reduced = mZ_reduced;
			this->E_gamma = mE_gamma;
		}

		double operator()(double theta_x) const noexcept
		{
			return get_F_single(theta_x, z_reduced, E_gamma)
				* sqrt(1.0 + 3.0 * pow(cos(theta_x - alpha), 2)) *
				sin(theta_x - alpha);
		}
	} obj_X;

	struct get_F_single_Z
	{
		double E_gamma;

		get_F_single_Z() {}

		void set_E_gamma(double mE_gamma) noexcept
		{
			this->E_gamma = mE_gamma;
		}

		double operator()(double z_reduced) const noexcept
		{
			obj_X.set_value(z_reduced, E_gamma);
			if (z_reduced < j0218->f * j0218->ratio)
			{
				return runge(obj_X, theta_m, theta_lc, 400) *
					(1.0 - j0218->g1);
			}
			else
			{
				return runge(obj_X, theta_m, theta_lc, 400) *
					(1.0 + j0218->g2);
			}
		}
	} obj_Z;

	double get_F_cur(double E_gamma) noexcept
	{
		obj_Z.set_E_gamma(E_gamma);

		return C_get_F_cur * runge(obj_Z, 0.0, j0218->f, 80); // 80 original
	}
}


/*
**	for pulsar B1821
*/
namespace B1821
{
	struct get_F_single_X
	{
		double z_reduced, E_gamma;

		get_F_single_X() = default;

		void set_Z(double mZ_reduced) noexcept
		{
			this->z_reduced = mZ_reduced;
		}

		void set_E_gamma(double mE_gamma) noexcept
		{
			this->E_gamma = mE_gamma;
		}

		void set_value(double mZ_reduced, double mE_gamma) noexcept
		{
			this->z_reduced = mZ_reduced;
			this->E_gamma = mE_gamma;
		}

		double operator()(double theta_x) const noexcept
		{
			return get_F_single(theta_x, z_reduced, E_gamma)
				* sqrt(1.0 + 3.0 * pow(cos(theta_x - alpha), 2)) *
				sin(theta_x - alpha);
		}
	} obj_X;

	struct get_F_single_Z
	{
		double E_gamma;

		get_F_single_Z() {}

		void set_E_gamma(double mE_gamma) noexcept
		{
			this->E_gamma = mE_gamma;
		}

		double operator()(double z_reduced) const noexcept
		{
			;
			obj_X.set_value(z_reduced, E_gamma);
			if (z_reduced < b1821->f * b1821->ratio)
			{
				return runge(obj_X, theta_m, theta_lc, 400) * (1.0 - b1821->g1);
			}
			else
			{
				return runge(obj_X, theta_m, theta_lc, 400) * (1.0 + b1821->g2);
			}
		}
	} obj_Z;

	double get_F_cur(double E_gamma) noexcept
	{
		obj_Z.set_E_gamma(E_gamma);

		return C_get_F_cur * runge(obj_Z, 0.0, b1821->f, 80); // 80 original
	}
}


//psr J1939
namespace J1939
{
	struct get_F_single_X
	{
		double z_reduced, E_gamma;

		get_F_single_X() = default;

		void set_Z(double mZ_reduced) noexcept
		{
			this->z_reduced = mZ_reduced;
		}

		void set_E_gamma(double mE_gamma) noexcept
		{
			this->E_gamma = mE_gamma;
		}

		void set_value(double mZ_reduced, double mE_gamma) noexcept
		{
			this->z_reduced = mZ_reduced;
			this->E_gamma = mE_gamma;
		}

		double operator()(double theta_x) const noexcept
		{
			return get_F_single(theta_x, z_reduced, E_gamma)
				* sqrt(1.0 + 3.0 * pow(cos(theta_x - alpha), 2)) *
				sin(theta_x - alpha);
		}
	} obj_X;

	struct get_F_single_Z
	{
		double E_gamma;

		get_F_single_Z() {}

		void set_E_gamma(double mE_gamma) noexcept
		{
			this->E_gamma = mE_gamma;
		}

		double operator()(double z_reduced) const noexcept
		{
			double border = j1939->f * j1939->ratio;
			obj_X.set_value(z_reduced, E_gamma);
			if (z_reduced < border)
			{
				return runge(obj_X, theta_m, theta_lc, 400) * (1.0 - j1939->g1);
			}
			else
			{
				return runge(obj_X, theta_m, theta_lc, 400) * (1.0 + j1939->g2);
			}
		}
	} obj_Z;

	double get_F_cur(double E_gamma) noexcept
	{
		obj_Z.set_E_gamma(E_gamma);
		return C_get_F_cur * runge(obj_Z, 0.0, j1939->f, 80); // 80 original
	}


}