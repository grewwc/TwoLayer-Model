#include "constants.h"
#include "curvature.h"
#include "property.h"
#include "runge.h"
#include "using.h"
#include <functional>
#include <iostream>
#include <tuple>
extern _j0218* j0218;
extern _b1821* b1821;
extern _j1939* j1939;
using std::cout;
using std::endl;

namespace J0218 {
constexpr double theta_0_tweek = 1.5;

namespace usage {
    constexpr double F_400 = 35.0 * 1e-26 / h;
    constexpr double beta_1 = -2.9221;
    constexpr double beta_2 = 0.5;
    constexpr double A = F_400 * (dist * dist / Rlc / Rlc) * 73.5167; //73.5 = pow(0.25, beta_1)
    inline double F(double nu_radio)
    {
        if (nu_radio >= 1e7) {
            return A * pow((nu_radio * 1e-8), beta_1);
        } else {
            return A * pow(0.1, beta_1) * sqrt(nu_radio * 1e-7);
        }
    }
}
}

namespace B1821 {
constexpr double theta_0_tweek = 1.15;
namespace usage {
    constexpr double F_400 = 40.0 * 1e-26 / h;
    constexpr double beta_1 = -2.3913;
    constexpr double beta_2 = 0.5;
    constexpr double A = F_400 * (dist * dist / Rlc / Rlc) * 73.5167; //73.5 = pow(0.25, beta_1)
    inline double F(double nu_radio)
    {
        if (nu_radio >= 1e7) {
            return A * pow((nu_radio * 1e-8), beta_1);
        } else {
            return A * pow(0.1, beta_1) * sqrt(nu_radio * 1e-7);
        }
    }
}
}

namespace J1939 {
constexpr double theta_0_tweek = 1.4;
namespace usage {
    constexpr double F_400 = 40.0 * 1e-26 / h;
    constexpr double beta_1 = -2.31522;
    constexpr double beta_2 = 0.5;
    constexpr auto A = F_400 * (dist * dist / Rlc / Rlc) * 73.5167; //73.5 = pow(0.25, beta_1)
    inline double F(double nu_radio)
    {
        if (nu_radio >= 1e7) {
            return A * pow((nu_radio * 1e-8), beta_1);
        } else {
            return A * pow(0.1, beta_1) * sqrt(nu_radio * 1e-7);
        }
    }
}
}

inline double toBeta(double gamma)
{
    // if (gamma > 1e6)
    // {
    // return 1.0 - 0.5 * std::pow(gamma, -2);
    // }
    // else
    // {
    // return sqrt(1.0 - std::pow(gamma, -2));
    // }
    return sqrt(1.0 - pow(gamma, -2));
    //return 1.0 - 0.5 * pow(gamma, -2);
}

//according to WiKi
inline double P(double theta, double E_gamma)
{
    return 1.0 / (1.0 + (1.0 - cos(theta)) * E_gamma / conbination::me_c_square);
}

inline double gamma_to_radio(double theta_1, double E_gamma, double gamma)
{
    double beta = toBeta(gamma);
    double theta_0 = sqrt(2.0 * j1939->f);
    //cout << "radio / gamma" << (ret / E_gamma) << endl;
    return E_gamma * (1.0 - beta * cos(theta_1)) / (1.0 - beta * cos(theta_0));
    //return E_gamma / (4 * pow(gamma, 2));
}

//theta_1 => theta''
//according to WiKi
namespace J0218 {
inline double d_sigma_over_d_Omega_mul_radio(double theta_1, double E_gamma, double gamma)
{
    double beta = toBeta(gamma);
    double theta_0 = sqrt(2.0 * j0218->f) / theta_0_tweek;
    double theta_prime_prime = atan(sin(theta_1) / (gamma * (cos(theta_1) - beta)));
    double theta_prime = atan(sin(theta_0) / (gamma * (cos(theta_0) - beta)));
    double E_radio = gamma_to_radio(theta_1, E_gamma, gamma);
    double E_gamma_prime = E_radio * gamma * (1.0 - beta * cos(theta_0));
    double theta = theta_prime + theta_prime_prime;
    double P_store = P(theta, E_gamma_prime);
    double res = 0.5 * pow(fine_structure_constant * r_c * P_store, 2) * (P_store + 1.0 / P_store - pow(sin(theta), 2));
    // std::cout<<E_radio/h<<std::endl;
    return res * J0218::usage::F(E_radio / h);
}

inline double d_P_IC_over_d_Omega(double theta_1_log, double E_gamma, double gamma)
{
    double theta_1 = exp(theta_1_log);
    double theta_0 = sqrt(2.0 * j0218->f) / theta_0_tweek;
    double beta = toBeta(gamma);
    //double D = 2.0 / (1.0 / gamma + gamma * pow(theta_1, 2));
    double D = 1.0 / (gamma * (1.0 - beta * cos(theta_1)));
    return pow(D, 2) * (1.0 - beta * cos(theta_0))
        * d_sigma_over_d_Omega_mul_radio(theta_1, E_gamma, gamma) * sin(theta_1) * 2.0 * PI * theta_1;
}

inline double get_IC_single(double theta_x, double z_reduced, double E_gamma)
{
    double gamma = get_Gamma_e(theta_x, z_reduced);
    auto res = std::bind(d_P_IC_over_d_Omega, std::placeholders::_1, E_gamma, gamma);
    return runge(res, -50.0, -0.1, 80);
}

struct get_IC_single_X {
    double z_reduced, E_gamma;

    get_IC_single_X() = default;

    void set_Z(double mZ_reduced)
    {
        this->z_reduced = mZ_reduced;
    }

    void set_E_gamma(double mE_gamma)
    {
        this->E_gamma = mE_gamma;
    }

    void set_value(double mZ_reduced, double mE_gamma)
    {
        this->z_reduced = mZ_reduced;
        this->E_gamma = mE_gamma;
    }

    double operator()(double theta_x) const noexcept
    {
        return get_IC_single(theta_x, z_reduced, E_gamma)
            * sqrt(1.0 + 3.0 * pow(cos(theta_x - J0218::alpha), 2)) * sin(theta_x - J0218::alpha);
    }
} obj_IC_X;

struct get_IC_single_Z {
    double E_gamma;
    mutable double border;

    get_IC_single_Z()
    {
    }

    void set_E_gamma(double mE_gamma)
    {
        this->E_gamma = mE_gamma;
    }

    double operator()(double mZ_reduced) const noexcept
    {
        border = j0218->f * j0218->ratio;
        obj_IC_X.set_value(mZ_reduced, this->E_gamma);
        if (mZ_reduced < border) {
            return runge(obj_IC_X, J0218::theta_m, J0218::theta_lc, 80) * (1.0 - j0218->g1);
        } else {
            return runge(obj_IC_X, J0218::theta_m, J0218::theta_lc, 80) * (1.0 + j0218->g2);
        }
    }
} obj_IC_Z;

double get_IC(double E_gamma)
{
    obj_IC_Z.set_E_gamma(E_gamma);
    return C_get_F_cur * runge(obj_IC_Z, 0.0, j0218->f, 12);
}
}

// for pulsar B1821
namespace B1821 {
inline double d_sigma_over_d_Omega_mul_radio(double theta_1, double E_gamma, double gamma)
{
    double beta = toBeta(gamma);
    double theta_0 = std::sqrt(2.0 * b1821->f) / theta_0_tweek;
    double theta_prime_prime = std::atan(sin(theta_1) / (gamma * (cos(theta_1) - beta)));
    double theta_prime = std::atan(sin(theta_0) / (gamma * (cos(theta_0) - beta)));
    double E_radio = gamma_to_radio(theta_1, E_gamma, gamma);
    double E_gamma_prime = E_radio * gamma * (1.0 - beta * std::cos(theta_0));
    double theta = theta_prime + theta_prime_prime;
    double P_store = P(theta, E_gamma_prime);
    double res = 0.5 * pow(fine_structure_constant * r_c * P_store, 2) * (P_store + 1.0 / P_store - pow(sin(theta), 2));
    // std::cout<<E_radio/h<<std::endl;
    return res * B1821::usage::F(E_radio / h);
}

inline double d_P_IC_over_d_Omega(double theta_1_log, double E_gamma, double gamma)
{
    double theta_1 = exp(theta_1_log);
    double theta_0 = sqrt(2.0 * b1821->f) / theta_0_tweek;
    double beta = toBeta(gamma);
    //double D = 2.0 / (1.0 / gamma + gamma * pow(theta_1, 2));
    double D = 1.0 / (gamma * (1.0 - beta * cos(theta_1)));
    return pow(D, 2) * (1.0 - beta * cos(theta_0))
        * d_sigma_over_d_Omega_mul_radio(theta_1, E_gamma, gamma) * sin(theta_1)
        * 2.0 * PI * theta_1;
}

inline double get_IC_single(double theta_x, double z_reduced, double E_gamma)
{
    double gamma = get_Gamma_e(theta_x, z_reduced);
    auto res = std::bind(d_P_IC_over_d_Omega, std::placeholders::_1, E_gamma, gamma);
    return runge(res, -50.0, -0.1, 160);
}

struct get_IC_single_X {
    double z_reduced, E_gamma;

    get_IC_single_X() = default;

    void set_Z(double mZ_reduced)
    {
        this->z_reduced = mZ_reduced;
    }

    void set_E_gamma(double mE_gamma)
    {
        this->E_gamma = mE_gamma;
    }

    void set_value(double mZ_reduced, double mE_gamma)
    {
        this->z_reduced = mZ_reduced;
        this->E_gamma = mE_gamma;
    }

    double operator()(double theta_x) const noexcept
    {
        return get_IC_single(theta_x, z_reduced, E_gamma)
            * sqrt(1.0 + 3.0 * pow(cos(theta_x - B1821::alpha), 2)) * sin(theta_x - alpha);
    }
} obj_IC_X;

struct get_IC_single_Z {
    double E_gamma;
    mutable double border;

    get_IC_single_Z()
    {
    }

    void set_E_gamma(double mE_gamma)
    {
        this->E_gamma = mE_gamma;
    }

    double operator()(double mZ_reduced) const noexcept
    {
        border = b1821->f * b1821->ratio;
        obj_IC_X.set_value(mZ_reduced, this->E_gamma);
        if (mZ_reduced < border) {
            return runge(obj_IC_X, theta_m, theta_lc, 120) * (1.0 - b1821->g1);
        } else {
            return runge(obj_IC_X, theta_m, theta_lc, 120) * (1.0 + b1821->g2);
        }
    }
} obj_IC_Z;

double get_IC(double E_gamma)
{
    obj_IC_Z.set_E_gamma(E_gamma);
    return C_get_F_cur * runge(obj_IC_Z, 0.0, b1821->f, 20);
}
}

//psr J1939
namespace J1939 {
inline double d_sigma_over_d_Omega_mul_radio(double theta_1, double E_gamma, double gamma)
{
    double beta = toBeta(gamma);
    double theta_0 = std::sqrt(2.0 * j1939->f) / theta_0_tweek;
    double theta_prime_prime = std::atan(sin(theta_1) / (gamma * (cos(theta_1) - beta)));
    double theta_prime = std::atan(sin(theta_0) / (gamma * (cos(theta_0) - beta)));
    double E_radio = gamma_to_radio(theta_1, E_gamma, gamma);
    double E_gamma_prime = E_radio * gamma * (1.0 - beta * std::cos(theta_0));
    double theta = theta_prime + theta_prime_prime;
    double P_store = P(theta, E_gamma_prime);
    double res = 0.5 * pow(fine_structure_constant * r_c * P_store, 2) * (P_store + 1.0 / P_store - pow(sin(theta), 2));
    // std::cout<<E_radio/h<<std::endl;
    return res * usage::F(E_radio / h);
}

inline double d_P_IC_over_d_Omega(double theta_1_log, double E_gamma, double gamma)
{
    double theta_1 = std::exp(theta_1_log);
    double theta_0 = std::sqrt(2.0 * j1939->f) / theta_0_tweek;
    double beta = toBeta(gamma);
    //double D = 2.0 / (1.0 / gamma + gamma * pow(theta_1, 2));
    double D = 1.0 / (gamma * (1.0 - beta * std::cos(theta_1)));
    return std::pow(D, 2) * (1.0 - beta * std::cos(theta_0))
        * d_sigma_over_d_Omega_mul_radio(theta_1, E_gamma, gamma) * std::sin(theta_1)
        * 2.0 * PI * theta_1;
}

inline double get_IC_single(double theta_x, double z_reduced, double E_gamma)
{
    double gamma = get_Gamma_e(theta_x, z_reduced);
    auto res = std::bind(d_P_IC_over_d_Omega, std::placeholders::_1, E_gamma, gamma);
    return runge(res, -50.0, -0.1, 160);
}

struct get_IC_single_X {
    double z_reduced, E_gamma;

    get_IC_single_X() = default;

    void set_Z(double mZ_reduced)
    {
        this->z_reduced = mZ_reduced;
    }

    void set_E_gamma(double mE_gamma)
    {
        this->E_gamma = mE_gamma;
    }

    void set_value(double mZ_reduced, double mE_gamma)
    {
        this->z_reduced = mZ_reduced;
        this->E_gamma = mE_gamma;
    }

    double operator()(double theta_x) const noexcept
    {
        return get_IC_single(theta_x, z_reduced, E_gamma)
            * std::sqrt(1.0 + 3.0 * std::pow(cos(theta_x - J1939::alpha), 2)) * std::sin(theta_x - J1939::alpha);
    }
} obj_IC_X;

struct get_IC_single_Z {
    double E_gamma;
    mutable double border;

    get_IC_single_Z()
    {
    }

    void set_E_gamma(double mE_gamma)
    {
        this->E_gamma = mE_gamma;
    }

    double operator()(double mZ_reduced) const noexcept
    {
        border = j1939->f * j1939->ratio;
        obj_IC_X.set_value(mZ_reduced, this->E_gamma);
        if (mZ_reduced < border) {
            return runge(obj_IC_X, J1939::theta_m, J1939::theta_lc, 120) * (1.0 - j1939->g1);
        } else {
            return runge(obj_IC_X, J1939::theta_m, J1939::theta_lc, 120) * (1.0 + j1939->g2);
        }
    }
} obj_IC_Z;

double get_IC(double E_gamma)
{
    obj_IC_Z.set_E_gamma(E_gamma);
    return C_get_F_cur * runge(obj_IC_Z, 0.0, j1939->f, 24);
}
}