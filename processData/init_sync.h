#pragma once


namespace init_sync_inner
{
	constexpr double LO = 1.6e-5;
	constexpr double HI = 4.8e-1;
}


namespace J0218
{
	namespace init_sync
	{
		using namespace init_sync_inner;
		void init(double *x, double *y, int N, double lo = LO, double hi = HI);
	}
}


namespace B1821
{
	namespace init_sync
	{
		using namespace init_sync_inner;
		void init(double *x, double *y, int N, double lo = LO, double hi = HI);
	}
}