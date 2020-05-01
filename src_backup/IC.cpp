#include "constants.h"
#include "curvature.h"
#include "property.h"
#include "runge.h"
#include <functional>
#include <iostream>
#include <tuple>

extern _j0218* j0218;

namespace usage
{

    const double theta_0 = sqrt(2 * j0218->f);
    constexpr double F_400 = 35 * 1e-26 / h;
    constexpr double beta_1 = -3.0;
    constexpr double beta_2 = 0.5;
    constexpr double A = F_400 * (dist * dist / Rlc / Rlc) * 73.5167; //73.5 = pow(0.25, beta_1)
    inline double F(double nu_radio)
    {
        if (nu_radio >= 1e7)
        {
            return A * pow((nu_radio * 1e-8), beta_1);
        } else
        {
            return A * pow(0.1, beta_1) * sqrt(nu_radio * 1e-7);
        }
    }

    inline double beta(double gamma)
    {
        if (gamma > 1e6)
        {
            return 1.0 - 0.5 * pow(gamma, -2.0);
        } else
        {
            return sqrt(1.0 - pow(gamma, -2.0));
        }
    }
}

//according to WiKi
inline double P(double theta, double E_gamma)
{
    return 1.0 / (1.0 + (1.0 - cos(theta)) * E_gamma / conbination::me_c_square);
}

inline std::tuple<double, double> transfer(double freq, double theta, double gamma)
{
    double freq_1 = gamma * freq * (1.0 - usage::beta(gamma) * cos(theta));
    double theta_1 = sin(theta) / (gamma * (cos(theta) - usage::beta(gamma)));
    return std::make_tuple(freq_1, theta_1);
}

inline double gamma_to_radio(double theta_1, double E_gamma, double gamma)
{
    double beta = usage::beta(gamma);
    return E_gamma * (1.0 - beta * cos(theta_1)) / (1.0 - beta * cos(usage::theta_0));
    // return E_gamma / (4 * pow(gamma, 2));
}

//theta_1 => theta''
//according to WiKi
inline double d_sigma_over_d_Omega_mul_radio(double theta_1, double E_gamma, double gamma)
{

    double beta = usage::beta(gamma);
    double theta_prime_prime = tanh(sin(theta_1 / (gamma * (cos(theta_1 - beta)))));
    double theta_prime = tanh(sin(usage::theta_0) / (gamma * (cos(usage::theta_0) - beta)));
    double E_radio = gamma_to_radio(theta_1, E_gamma, gamma);
    double E_gamma_prime = gamma * (1.0 + beta * cos(usage::theta_0));
    double theta = theta_prime + theta_prime_prime;
    double P_store = P(theta, E_gamma_prime);
    double res =
            0.5 * pow(fine_structure_constant * r_c * P_store, 2.0) * (P_store + 1.0 / P_store - pow(sin(theta), 2));
    // std::cout<<E_radio/h<<std::endl;
    return res * usage::F(E_radio / h);
}

inline double d_P_IC_over_d_Omega(double theta_1_log, double E_gamma, double gamma)
{
    double theta_1 = exp(theta_1_log);
    double beta = usage::beta(gamma);
    double D = 2.0 / (1.0 / gamma + gamma * pow(theta_1, 2.0));
    return pow(D, 2.0) * (1.0 - beta * cos(usage::theta_0))
           * d_sigma_over_d_Omega_mul_radio(theta_1, E_gamma, gamma) * sin(theta_1)
           * 2.0 * PI * theta_1;
}

inline double get_IC_single(double x_reduced, double z_reduced, double E_gamma)
{
    double gamma = get_Gamma_e(x_reduced, z_reduced);
    auto res = std::bind(d_P_IC_over_d_Omega, std::placeholders::_1, E_gamma, gamma);

    return runge(res, -30.0, -2, 400);
}

struct get_IC_single_X
{
    double z_reduced, E_gamma;

    get_IC_single_X() = default;

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
        return get_IC_single(x_reduced, z_reduced, E_gamma);
    }
} obj_IC_X;

const double border = j0218->f * j0218->ratio;

struct get_IC_single_Z
{
    double E_gamma;

    get_IC_single_Z()
    {}

    void set_E_gamma(double E_gamma)
    {
        this->E_gamma = E_gamma;
    }

    double operator()(double z_reduced)
    {
        obj_IC_X.set_value(z_reduced, E_gamma);
        if (z_reduced < border)
        {
            return runge(obj_IC_X, x_m_reduced, 1.0, 1) * (1.0 - j0218->g1);
        } else
        {
            return runge(obj_IC_X, x_m_reduced, 1.0, 1) * (1.0 + j0218->g2);
        }
    }
} obj_IC_Z;

double get_IC(double E_gamma)
{
    obj_IC_Z.set_E_gamma(E_gamma);
    return C_get_F_cur * runge(obj_IC_Z, 0.0, j0218->f, 100);
}