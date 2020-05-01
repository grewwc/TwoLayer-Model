#pragma once
#include "using.h"
#include <cmath>

namespace appro {
constexpr double _55_over_72 = 55.0 / 72.0;
constexpr double _10151_over_10368 = 10151.0 / 10368.0;
constexpr double _1_over_3 = 1.0 / 3.0;
constexpr double a = 5.4158;
constexpr double b = 0.96;
constexpr double c = 0.578;
constexpr double d = 1.2786;
}

inline double approximation(double x)
{
    if (x < 1e4) {
        double _temp = x + appro::d;
        return appro::a * pow(1.0 / x + appro::b, -appro::c * x - appro::_1_over_3) * 1.2533141 * exp(-_temp) * sqrt(_temp)
            * (1.0 + appro::_55_over_72 / _temp - appro::_10151_over_10368 * pow(_temp, -2));
    } else {
        return 0.0;
    }
}
