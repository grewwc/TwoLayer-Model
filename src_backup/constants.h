#pragma once

constexpr double speed_of_light = 3e10;
constexpr double PI = 3.1416;
constexpr double e = 4.8032E-10;

constexpr double h_bar = 1.0545716e-27;
constexpr double h = 6.626e-27;
constexpr double r_e = 2.81794e-13;
constexpr double m_e = 9.11e-28;

constexpr double r_c = h_bar / (m_e * speed_of_light);
constexpr double fine_structure_constant = 1.0 / 137.04;

namespace conbination {

    constexpr double me_c_square = m_e * speed_of_light * speed_of_light;
}
