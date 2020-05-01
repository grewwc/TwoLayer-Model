//
// Created by grewwc on 27/12/2017.
//

#include "init_total.h"
#include <cmath>
#include <iostream>
#include "IC.h"
#include "curvature.h"

extern _j0218 *j0218;
extern _b1821 *b1821;

namespace J0218
{
	namespace init_total
	{
		void init(double *x, double *y, int N, double lo, double hi)
		{
			double log_lo = log(lo), log_hi = log(hi);
			double step = (log_hi - log_lo) / (N - 1);
			double temp;
			for (int i = 0; i < N; ++i)
			{
				x[i] = exp(log_lo);
				temp = x[i] * get_IC(x[i]) / pow(J0218::dist, 2) / j0218->omega;
				y[i] = temp + get_F_cur(x[i]) * pow(x[i] / J0218::dist, 2) / j0218->omega;
				log_lo += step;
			}
		}

	}
}

namespace B1821
{
	namespace init_total
	{
		void init(double *x, double *y, int N, double lo, double hi)
		{
			double log_lo = log(lo), log_hi = log(hi);
			double step = (log_hi - log_lo) / (N - 1);
			double temp;
			for (int i = 0; i < N; ++i)
			{
				x[i] = exp(log_lo);
				temp = x[i] * get_IC(x[i]) / pow(B1821::dist, 2) / b1821->omega;
				y[i] = temp + get_F_cur(x[i]) * pow(x[i] / B1821::dist, 2) / b1821->omega;
				log_lo += step;
			}
		}

	}
}