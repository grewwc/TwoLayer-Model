#ifndef __curvature_h
#define __curvature_h

#include "constants.h"
#include "property.h"

constexpr double C_get_Gamma_e = 1.5 * Rlc * Rlc / e;
constexpr double C_get_E_cur = 1.5 * h_bar * speed_of_light / Rlc;
constexpr double C_get_F_single = 1.732 * e * e / (2 * PI * h_bar * Rlc);
constexpr double C_get_F_cur = F0 * (2 * PI * Rlc) * Blc / (speed_of_light * e) * Rlc * Rlc;


inline double get_Gamma_e(double x_reduced, double z_reduced)
{
    double res = pow(C_get_Gamma_e * get_E_parallel(x_reduced, z_reduced), 0.25);
    return res > 1 ? res : 10;
}

double get_E_cur(double x_reduced, double z_reduced);

double get_F_single(double x_reduced, double z_reduced, double E_gamma);

struct get_F_single_Z;

double get_F_cur(double E_gamma);

#endif


