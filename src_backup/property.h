#ifndef __property_h
#define __property_h

#include "constants.h"
#include <cmath>
#include <vector>

constexpr double F0 = 430.461054545748;
constexpr double Rlc = speed_of_light / (2 * PI * F0);
constexpr double Bs = 4.29e8;
constexpr double Rs = 1e6;
constexpr double Blc = Bs * (Rs * Rs * Rs) / (Rlc * Rlc * Rlc);
constexpr double alpha = 0.2;
constexpr double dist = 3.15 * 3.068e21;
constexpr double C_E_parallel = F0 * 2 * PI / (speed_of_light * Rlc);

const double theta_m = tanh(0.5 * (sqrt(8 + pow(3 * tan(alpha), 2)) + 3 * tan(alpha)));
const double theta_lc = PI - tanh((sqrt(9 + 8 * pow(tan(alpha), 2)) + 3) / (4 * tan(alpha)));
const double L = Rlc / pow(sin(theta_lc - alpha), 2.0);
const double x_m_reduced = pow(sin(theta_m), 2.0) * L / Rlc;

struct _j0218
{
    double f = 0.2;
    double g1 = 0.97;
    double ratio = 0.89;
    double omega = 0.2;

    double g2 = 1.0 / (1.0 / pow(ratio, 2.0) - 1.0);
    double h2 = f * Rlc;
    double h1 = ratio * h2;

    void set(const std::vector<double>& params)
    {
        this->f = params[0];
        this->g1 = params[1];
        this->ratio = params[2];
        this->omega = params[3];
        this->g2 = 1.0 / (1.0 / pow(ratio, 2.0) - 1.0);
        this->h2 = f * Rlc;
        this->h1 = ratio * h2;
    }
};



inline double get_B(double x)
{
    return Bs * pow(Rs / x, 3.0);
}

inline double get_julian_density(double x)
{
    return F0 * get_B(x) / speed_of_light;
}

double get_E_parallel(double x_reduced, double z_reduced);


#endif
