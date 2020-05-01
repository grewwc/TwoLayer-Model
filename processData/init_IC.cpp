//
// Created by grewwc on 27/12/2017.
//

#include <IC.h>
#include "property.h"
#include "init_IC.h"


extern _j0218 *j0218;
extern _b1821 *b1821;

namespace J0218
{
	namespace init_IC
	{
		void init(double *x, double *y, int N, double lo, double hi)
		{
			double log_lo = log(lo), log_hi = log(hi);
			double step = (log_hi - log_lo) / (N - 1);

			for (int i = 0; i < N; ++i)
			{
				x[i] = exp(log_lo);
				y[i] = x[i] * get_IC(x[i]) / pow(J0218::dist, 2) / j0218->omega;
				log_lo += step;

			}
		}

	}
}

namespace B1821
{
	namespace init_IC
	{
		void init(double *x, double *y, int N, double lo, double hi)
		{
			double log_lo = log(lo), log_hi = log(hi);
			double step = (log_hi - log_lo) / (N - 1);

			for (int i = 0; i < N; ++i)
			{
				x[i] = exp(log_lo);
				y[i] = x[i] * get_IC(x[i]) / pow(B1821::dist, 2) / b1821->omega;
				log_lo += step;

			}
		}

	}
}