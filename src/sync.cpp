//
// Created by grewwc on 27/12/2017.
//

#include "sync.h"
#include "curvature.h"
#include "constants.h"
#include "runge.h"
#include "approximation.h"
#include <iostream>
#include "using.h"

namespace J0218
{
	const auto sin_theta_s = sqrt(2.0 * j0218->f);
	const auto E_sync_coef = 3.0 * h * e * global_B * sin_theta_s /
		(4.0 * PI * m_e * speed_of_light);
	const auto E_gamma_lo_coef = log(2.0 * m_e * pow(speed_of_light, 2.0));

	const auto calc_tao = 1.0 - exp(-tao);

	const auto all_sync_coef = all_sync_coef_without_sin / sin_theta_s * calc_tao;
}


namespace B1821
{
	const auto sin_theta_s = sqrt(2.0 * b1821->f);
	const auto E_sync_coef = 3.0 * h * e * global_B * sin_theta_s /
		(4.0 * PI * m_e * speed_of_light);
	const auto E_gamma_lo_coef = log(2.0 * m_e * pow(speed_of_light, 2.0));

	const auto calc_tao = 1.0 - exp(-tao);

	const auto all_sync_coef = all_sync_coef_without_sin / sin_theta_s * calc_tao;
}

namespace J1939
{
	const auto sin_theta_s = sqrt(2.0 * j1939->f);
	const auto E_sync_coef = 3.0 * h * e * global_B * sin_theta_s /
		(4.0 * PI * m_e * speed_of_light);
	const auto E_gamma_lo_coef = log(2.0 * m_e * pow(speed_of_light, 2.0));

	const auto calc_tao = 1.0 - exp(-tao);

	const auto all_sync_coef = all_sync_coef_without_sin / sin_theta_s * calc_tao;
}



//psr b1821
namespace B1821
{
	class integrate_E_gamma
	{
	private:
		double log_gamma_s;
		double lo, hi;
	public:
		integrate_E_gamma(double aLog_gamma_s) noexcept
			: log_gamma_s{ aLog_gamma_s }
		{
			//lo = E_gamma_lo_coef + log_gamma_s;
			lo = E_gamma_lo_coef;
			hi = log(E_gamma_hi);
		}

		integrate_E_gamma() noexcept
		{
			hi = log(E_gamma_hi);
		};

		void set_gamma_s(double alog_gamma_s) noexcept
		{
			this->log_gamma_s = alog_gamma_s;
			this->lo = E_gamma_lo_coef + log_gamma_s;
		}

		decltype(auto) operator()(double log_E_gamma) const noexcept
		{
			double E_gamma = exp(log_E_gamma);
			return get_F_cur(E_gamma) * E_gamma;
		}

		decltype(auto) doIntegrate() const noexcept
		{
			return runge(*this, lo, hi, 200);
		}

	};

	namespace _inner
	{
		inline decltype(auto) get_E_sync(double gamma_s) noexcept
		{
			return E_sync_coef * pow(gamma_s, 2);
		}
	}

	class integrate_gamma_s
	{
	private:
		double E_sync;
		mutable integrate_E_gamma _temp;
	public:

		explicit integrate_gamma_s(double aEsync) noexcept
			: E_sync{ aEsync }
		{}

		void setE_sync(double aEsync) noexcept
		{
			this->E_sync = aEsync;
		}

		decltype(auto) operator()(double log_gamma_s) const noexcept
		{
			auto gamma_s = exp(log_gamma_s);
			_temp.set_gamma_s(log_gamma_s);
			auto x = this->E_sync / _inner::get_E_sync(gamma_s);

			/*std::cout << "res : " << res << std::endl;*/
			// return _temp.doIntegrate() * approximation(x) / gamma_s;
			return approximation(x) / gamma_s;
		}

		decltype(auto) doIntegrate() const noexcept
		{
			return runge(*this, 1.0, 9.0, 200) * all_sync_coef;
		}
	};

	double get_Sync(double E_sync) noexcept
	{
		integrate_gamma_s _temp_gamma_s{ E_sync };
		return _temp_gamma_s.doIntegrate() * 6e35;
	}
}


//psr J0218
namespace J0218
{
	class integrate_E_gamma
	{
	private:
		double log_gamma_s;
		double lo, hi;
	public:
		integrate_E_gamma(double aLog_gamma_s) noexcept
			: log_gamma_s{ aLog_gamma_s }
		{
			//lo = E_gamma_lo_coef + log_gamma_s;
			lo = E_gamma_lo_coef;
			hi = log(E_gamma_hi);
		}

		integrate_E_gamma() noexcept
		{
			hi = log(E_gamma_hi);
		};

		void set_gamma_s(double alog_gamma_s) noexcept
		{
			this->log_gamma_s = alog_gamma_s;
			this->lo = log_gamma_s + E_gamma_lo_coef;
		}

		decltype(auto) operator()(double log_E_gamma) const noexcept
		{
			auto E_gamma = exp(log_E_gamma);
			return get_F_cur(E_gamma) * E_gamma;
		}

		decltype(auto) doIntegrate() const noexcept
		{
			return runge(*this, lo, hi, 200);
		}

	};


	namespace _inner
	{
		inline decltype(auto) get_E_sync(double gamma_s)
		{
			return E_sync_coef * pow(gamma_s, 2);
		}
	}


	class integrate_gamma_s
	{
	private:
		double E_sync;
		mutable integrate_E_gamma _temp;
	public:

		explicit integrate_gamma_s(double aEsync) noexcept
			: E_sync{ aEsync }
		{}

		void setE_sync(double aEsync) noexcept
		{
			this->E_sync = aEsync;
		}

		decltype(auto) operator()(double log_gamma_s) const noexcept
		{
			auto gamma_s = exp(log_gamma_s);
			//_temp.set_gamma_s(log_gamma_s);
			auto x = this->E_sync / _inner::get_E_sync(gamma_s);

			// return _temp.doIntegrate() * approximation(x) / gamma_s;
			return approximation(x) / gamma_s;
		}

		decltype(auto) doIntegrate() const noexcept
		{
			return runge(*this, 1.0, 9.0, 200) * all_sync_coef;
		}
	};

	double get_Sync(double E_sync) noexcept
	{
		integrate_gamma_s _temp_gamma_s{ E_sync };
		return _temp_gamma_s.doIntegrate() * 0.45e35;
		// return _temp_gamma_s.doIntegrate();
	}


}


namespace J1939
{
	class integrate_E_gamma
	{
	private:
		double log_gamma_s;
		double lo, hi;
	public:
		integrate_E_gamma(double aLog_gamma_s) noexcept
			: log_gamma_s{ aLog_gamma_s }
		{
			//lo = E_gamma_lo_coef + log_gamma_s;
			lo = E_gamma_lo_coef;
			hi = log(E_gamma_hi);
		}

		integrate_E_gamma() noexcept
		{
			hi = log(E_gamma_hi);
		};

		void set_gamma_s(double alog_gamma_s) noexcept
		{
			this->log_gamma_s = alog_gamma_s;
			this->lo = E_gamma_lo_coef + log_gamma_s;
		}

		decltype(auto) operator()(double log_E_gamma) const noexcept
		{
			double E_gamma = exp(log_E_gamma);
			return get_F_cur(E_gamma) * E_gamma;
		}

		decltype(auto) doIntegrate() const noexcept
		{
			return runge(*this, lo, hi, 200);
		}

	};

	namespace _inner
	{
		inline decltype(auto) get_E_sync(double gamma_s)
		{
			return E_sync_coef * pow(gamma_s, 2);
		}
	}

	class integrate_gamma_s
	{
	private:
		double E_sync;
		mutable integrate_E_gamma _temp;
	public:

		explicit integrate_gamma_s(double aEsync) noexcept
			: E_sync{ aEsync }
		{}

		void setE_sync(double aEsync) noexcept
		{
			this->E_sync = aEsync;
		}

		decltype(auto) operator()(double log_gamma_s) const noexcept
		{
			auto gamma_s = exp(log_gamma_s);
			_temp.set_gamma_s(log_gamma_s);
			auto x = this->E_sync / _inner::get_E_sync(gamma_s);

			/*std::cout << "res : " << res << std::endl;*/
			// return _temp.doIntegrate() * approximation(x) / gamma_s;
			return approximation(x) / gamma_s;
		}

		decltype(auto) doIntegrate() const noexcept
		{
			return runge(*this, 1.0, 9.0, 200) * all_sync_coef;
		}
	};

	double get_Sync(double E_sync) noexcept
	{
		integrate_gamma_s _temp_gamma_s{ E_sync };
		return _temp_gamma_s.doIntegrate() * 7e35;
	}
}
