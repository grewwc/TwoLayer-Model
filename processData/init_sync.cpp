#include "init_sync.h"
#include "../src/sync.h"
#include "../src/property.h"
#include <cmath>
extern _j0218* j0218;
extern _b1821* b1821;


namespace J0218
{
	namespace init_sync
	{
		void init(double *x, double *y, int N, double lo /*= 1.6e-10*/, double hi /*= 1.6e-5*/)
		{
			double log_lo = log(lo), log_hi = log(hi);
			double step = (log_hi - log_lo) / (N - 1);
			double temp;
			for (int i = 0; i < N; ++i)
			{
				x[i] = exp(log_lo);
				y[i] = get_Sync(x[i]) * pow(x[i] / J0218::dist, 2) / j0218->omega;
				log_lo += step;
			}
		}

	}
}

namespace B1821
{
	namespace init_sync
	{
		void init(double *x, double *y, int N, double lo /*= 1.6e-10*/, double hi /*= 1.6e-5*/)
		{
			double log_lo = log(lo), log_hi = log(hi);
			double step = (log_hi - log_lo) / (N - 1);
			double temp;
			for (int i = 0; i < N; ++i)
			{
				x[i] = exp(log_lo);
				y[i] = get_Sync(x[i]) * pow(x[i] / B1821::dist, 2) / b1821->omega;
				log_lo += step;
			}
		}

	}
}