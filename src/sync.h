//
// Created by grewwc on 27/12/2017.
//

#ifndef TWOLAYER_SYNC_H
#define TWOLAYER_SYNC_H
#include "property.h"
#include "using.h"

extern _j0218 *j0218;
extern _b1821 *b1821;
extern _j1939 *j1939;
constexpr double tao = 0.02;
constexpr double E_gamma_hi = 1.6e-1;



namespace J0218
{
	// const auto global_B = sqrt(Bs * Blc) / 10.0;
	const auto global_B = Blc;
	const auto all_sync_coef_without_sin = 1.732 / 2.0 * 3.0 *
		m_e * m_e * speed_of_light * speed_of_light * speed_of_light
		/ (e * global_B * h);
	double get_Sync(double E_gamma) noexcept;
}

namespace B1821
{
	// const auto global_B = sqrt(Bs * Blc);
	const auto global_B = Blc;
	const auto all_sync_coef_without_sin = 1.732 / 2.0 * 3.0 *
		m_e * m_e * speed_of_light * speed_of_light * speed_of_light
		/ (e * global_B * h);

	double get_Sync(double E_gamma) noexcept;

}


namespace J1939
{
	// const auto global_B = sqrt(Bs * Blc);
	const auto global_B = Blc;
	const auto all_sync_coef_without_sin = 1.732 / 2.0 * 3.0 *
		m_e * m_e * speed_of_light * speed_of_light * speed_of_light
		/ (e * global_B * h);

	double get_Sync(double E_gamma) noexcept;

}
#endif //TWOLAYER_SYNC_H












