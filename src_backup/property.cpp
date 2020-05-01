#include "property.h"

extern _j0218* j0218;

struct
{
    const double c1_prime = -j0218->g1 * j0218->h1 * (j0218->ratio - 2) - j0218->g2 * pow(j0218->ratio, 2.0);
    const double d1_prime = -(j0218->g1 + j0218->g2) * j0218->h1 * j0218->ratio - j0218->g2 * j0218->h2;
} Const;

double get_E_parallel(double x_reduced, double z_reduced)
{
    double x = x_reduced * Rlc;
    double z = Rlc * z_reduced;

    if (z <= j0218->h1) {
        return C_E_parallel * get_B(x) * ((-j0218->g1 * z + Const.c1_prime)) * z;
    } else {
        double res = C_E_parallel * get_B(x) * (z - j0218->h2) * ((z + j0218->h2) + Const.d1_prime);
        return res > 0 ? res : 0;
    }
}
