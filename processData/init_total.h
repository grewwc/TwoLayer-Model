//
// Created by grewwc on 27/12/2017.
//

#ifndef TWOLAYER_TOTAL_H
#define TWOLAYER_TOTAL_H

namespace init_total_inner
{
	constexpr double LO = 1.6e-5;
	constexpr double HI = 4.8e-1;
}

namespace J0218
{
	namespace init_total
	{
		using namespace init_total_inner;
		void init(double *x, double *y, int N, double lo = LO, double hi = HI);
	}
}

namespace B1821
{
	namespace init_total
	{
		using namespace init_total_inner;
		void init(double *x, double *y, int N, double lo = LO, double hi = HI);
	}
}



#endif //TWOLAYER_TOTAL_H
