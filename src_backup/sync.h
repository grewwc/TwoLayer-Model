//
// Created by grewwc on 27/12/2017.
//

#ifndef TWOLAYER_SYNC_H
#define TWOLAYER_SYNC_H

#include <cmath>
#include "property.h"
#include "../star.h"


extern _j0218 *j0218;

constexpr double tao = 0.02;
constexpr double global_gamma_s = 1e5;
constexpr double global_B = Blc;
constexpr double E_gamma_hi = 4.8e-1;
constexpr double gamma_s_hi = 1e5;
constexpr double all_sync_coef_without_sin = 1.732 / 2 * 3 * m_e * m_e * speed_of_light * speed_of_light * speed_of_light
/ (e * global_B * h);


//struct
//{
//	const double x_reduced = 1;
//	const double z_reduced = j0218->f;
//} sync_pos_reduced;




double get_Sync(double E_gamma);


#endif //TWOLAYER_SYNC_H












