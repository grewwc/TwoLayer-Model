#ifndef __property_h
#define __property_h

#include "constants.h"
#include "using.h"
#include <cmath>
#include <vector>

constexpr double Rs = 1e6;
namespace J0218 {
constexpr double F0 = 430.461054545748;
constexpr double Bs = 4.29e8;
constexpr auto Rlc = speed_of_light / (2.0 * PI * F0);
constexpr double Blc = Bs * (Rs * Rs * Rs) / (Rlc * Rlc * Rlc);
constexpr double C_E_parallel = F0 * 2.0 * PI / (speed_of_light * Rlc);
constexpr double alpha = 1e-3;
constexpr double dist = 3.15 * 3.068e21;
const double theta_m = atan(0.5 * (sqrt(8.0 + pow(3.0 * tan(alpha), 2)) + 3.0 * tan(alpha)));
const double theta_lc = PI - atan((sqrt(9.0 + 8.0 * pow(tan(alpha), 2)) + 3.0) / (4.0 * tan(alpha)));

const double L = Rlc / (pow(sin(theta_lc - alpha), 2) * sin(theta_lc));
}

namespace B1821 {
constexpr double Bs = 2.25e9;
constexpr double F0 = 327.40559048006;
constexpr double Rlc = speed_of_light / (2.0 * PI * F0);
constexpr double Blc = Bs * (Rs * Rs * Rs) / (Rlc * Rlc * Rlc);
constexpr double C_E_parallel = F0 * 2.0 * PI / (speed_of_light * Rlc);
constexpr double alpha = 1e-3;
constexpr double dist = 5.5 * 3.068e21;
const auto theta_m = atan(0.5 * (sqrt(8.0 + pow(3.0 * tan(alpha), 2)) + 3.0 * tan(alpha)));
const auto theta_lc = PI - atan((sqrt(9.0 + 8.0 * pow(tan(alpha), 2)) + 3.0) / (4.0 * tan(alpha)));

const double L = Rlc / (pow(sin(theta_lc - alpha), 2) * sin(theta_lc));
//const double alpha_0 = alpha + std::asin(std::sqrt(Rs / L));
}

namespace J1939 {
constexpr auto Bs = 4.09e8;
constexpr auto F0 = 641.92822458236;
constexpr auto Rlc = speed_of_light / (2.0 * PI * F0);
constexpr auto Blc = Bs * (Rs * Rs * Rs) / (Rlc * Rlc * Rlc);
constexpr auto C_E_parallel = F0 * 2.0 * PI / (speed_of_light * Rlc);
constexpr auto alpha = 1e-3;
constexpr double dist = 3.5 * 3.068e21;
const auto theta_m = atan(0.5 * (sqrt(8.0 + pow(3.0 * tan(alpha), 2)) + 3.0 * tan(alpha)));
const auto theta_lc = PI - atan((sqrt(9.0 + 8.0 * pow(tan(alpha), 2)) + 3.0) / (4.0 * tan(alpha)));

const auto L = Rlc / (pow(sin(theta_lc - alpha), 2) * sin(theta_lc));
//const double alpha_0 = alpha + std::asin(std::sqrt(Rs / L));
}

struct _star {
public:
    void virtual set(const std::vector<double>& params) = 0;
};

struct _j0218 : public _star {
private:
    double _ratio_square = ratio * ratio;

public:
    double f = 0.2;
    double g1 = 0.9;
    double ratio = 0.9;
    double omega = 0.2;

    //double g2 = 1.0 / (1.0 / pow(ratio, 2.0) - 1.0);
    double g2 = g1 * _ratio_square / (1.0 - _ratio_square);
    double h2 = f * J0218::Rlc;
    double h1 = ratio * h2;

    void set(const std::vector<double>& params) override;
};

struct _b1821 : public _star {
private:
    double _ratio_square = ratio * ratio;

public:
    double f = 0.2;
    double g1 = 0.9;
    double ratio = 0.9;
    double omega = 0.2;

    //double g2 = 1.0 / (1.0 / pow(ratio, 2.0) - 1.0);
    double g2 = g1 * _ratio_square / (1.0 - _ratio_square);
    double h2 = f * B1821::Rlc;
    double h1 = ratio * h2;

    void set(const std::vector<double>& params) override;
};

struct _j1939 : public _star {
private:
    double _ratio_square = ratio * ratio;

public:
    double f = 0.2;
    double g1 = 0.9;
    double ratio = 0.9;
    double omega = 0.2;

    //double g2 = 1.0 / (1.0 / pow(ratio, 2.0) - 1.0);
    double g2 = g1 * _ratio_square / (1.0 - _ratio_square);
    double h2 = f * J1939::Rlc;
    double h1 = ratio * h2;

    void set(const std::vector<double>& params) override;
};

namespace J0218 {
inline decltype(auto) get_B(double theta_x)
{

    //this sin_2 is only for speed the program
    auto sin_2 = pow(sin(theta_x - alpha), 2);
    auto _R = L * sin_2;
    return Bs * pow(Rs / _R, 3);
}

//inline double get_julian_density(double theta_x)
//{
//	return F0 * get_B(theta_x) / speed_of_light;
//}

double get_E_parallel(double theta_x, double z_reduced);
}

namespace B1821 {
inline decltype(auto) get_B(double theta_x)
{

    //this sin_2 is only for speed the program
    auto sin_2 = pow(sin(theta_x - alpha), 2);
    auto _R = L * sin_2;
    return Bs * pow(Rs / _R, 3);
}

//inline double get_julian_density(double theta_x)
//{
//	return F0 * get_B(theta_x) / speed_of_light;
//}

double get_E_parallel(double theta_x, double z_reduced);
}

namespace J1939 {
inline decltype(auto) get_B(double theta_x)
{

    //this sin_2 is only for speed the program
    auto sin_2 = pow(sin(theta_x - alpha), 2);
    auto _R = L * sin_2;
    return Bs * pow(Rs / _R, 3);
}

//inline double get_julian_density(double theta_x)
//{
//	return F0 * get_B(theta_x) / speed_of_light;
//}

double get_E_parallel(double theta_x, double z_reduced);
}

#endif