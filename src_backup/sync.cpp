//
// Created by grewwc on 27/12/2017.
//

#include "sync.h"
#include "curvature.h"
#include "constants.h"
#include "runge.h"
#include "approximation.h"
#include <iostream>


const double sin_theta_s = sqrt(2 * j0218->f);
const double E_sync_coef = 3 * h * e * global_B * sin_theta_s /
(4 * PI * m_e * speed_of_light);
const double E_gamma_lo_coef = log(2 * m_e * pow(speed_of_light, 2));

const double calc_tao = 1 - exp(-tao);
const double log_gamma_s_hi = log(gamma_s_hi);

const double all_sync_coef = all_sync_coef_without_sin / sin_theta_s * calc_tao;
inline double F_sync_global_per_E(double E)
{
	return get_F_cur(E) * calc_tao;
}

class integrate_E_gamma
{
private:
	double log_gamma_s;
	double lo, hi;
public:
	integrate_E_gamma(double aLog_gamma_s) noexcept
		: log_gamma_s{ aLog_gamma_s }
	{
		lo = E_gamma_lo_coef + log_gamma_s;
		hi = log(E_gamma_hi);
	}

	integrate_E_gamma() = default;

	void set_gamma_s(double log_gamma_s)
	{
		this->log_gamma_s = log_gamma_s;
		this->lo = E_gamma_lo_coef + log_gamma_s;
	}

	double operator()(double log_E_sync) const noexcept
	{
		double E_sync = exp(log_E_sync);
		return get_F_cur(E_sync);
	}

	double doIntegrate() const noexcept
	{
		return runge(*this, lo, hi, 40);
	}

};


namespace _inner
{
	inline double get_E_sync(double gamma_s)
	{
		return E_sync_coef * pow(gamma_s, 2);
	}
}


class integrate_gamma_s
{
private:
	double E_sync;
	integrate_E_gamma _temp;
public:

	integrate_gamma_s(double aEsync) noexcept
		: E_sync{ aEsync }
	{}

	void setE_sync(double aEsync) noexcept
	{
		this->E_sync = aEsync;
	}

	double operator()(double log_gamma_s) noexcept
	{
		double gamma_s = exp(log_gamma_s);
		_temp.set_gamma_s(log_gamma_s);
		double x = this->E_sync / _inner::get_E_sync(gamma_s);
		return _temp.doIntegrate() * approximation(x) / gamma_s;
	}

	double doIntegrate() const noexcept
	{
		return runge(*this, 0, log_gamma_s_hi, 40) * all_sync_coef;
	}
};





double get_Sync(double E_sync)
{
	integrate_gamma_s _temp_gamma_s{ E_sync };
	return _temp_gamma_s.doIntegrate();
}






