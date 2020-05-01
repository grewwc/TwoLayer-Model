#pragma once

#include <functional>
#include <future>
#include <iostream>
#include <cmath>
//const std::function<double(double)>& f
using std::abs;


template<typename Func>
inline double runge_once_with_n(const Func& f, double lo, double hi, int n)
{
	//double k1 = 0.0, k2, k3, k4;
	double k1 = 0.0, k2, k4; // k2=k3, so ignore k3 for speed.
	double tn = lo;
	double y_n = 0.0;
	//double y = 0.0;
	const double h = (hi - lo) / n;
	const double h_half = h * 0.5;
	//const double h_6 = h / 6.0;
	for (int i = 0; i < n; ++i)
	{
		k1 = f(tn);
		tn += h_half;
		k2 = f(tn);
		//k3 = k2;
		tn += h_half;
		k4 = f(tn);
		y_n += k1 + 4.0 * k2 + k4;
		//y = y_n;
		//tn += h;
	}

	return y_n * h / 6.0;
}

template<typename Func>
inline double runge_once(const Func& f, double lo, double hi,
	int n = 4, double relative_err = 1e-5)
{
	double real_err;
	auto upper = static_cast<int>(1e4f);
	double y_2n, y_n = runge_once_with_n(f, lo, hi, n);
	do
	{
		y_2n = runge_once_with_n(f, lo, hi, n <<= 1);
		real_err = abs(y_2n - y_n) / abs(y_2n);
		y_n = y_2n;

	} while (real_err < relative_err && n < upper);


	return y_2n;
}

inline double runge(const std::function<double(double)>& f,
	double lo, double hi, int n = 4, double relative_err = 1e-4)
{
	if (n >= 4)
	{
		using the_func = std::function<double(double)>;
		const double step = (hi - lo) * 0.25;
		const double mid1 = lo + step, mid2 = lo + 2.0 * step, mid3 = hi - step;
		const int each_n = n >> 2;
		std::future<double> fu1 = std::async(std::launch::async,
			runge_once_with_n<the_func>,
			std::ref(f), lo, mid1, each_n);

		std::future<double> fu2 = std::async(std::launch::async, 
			runge_once_with_n<the_func>,
			std::ref(f), mid1, mid2, each_n);

		std::future<double> fu3 = std::async(std::launch::async, 
			runge_once_with_n<the_func>,
			std::ref(f), mid2, mid3, each_n);

		std::future<double> fu4 = std::async(std::launch::async, 
			runge_once_with_n<the_func>,
			std::ref(f), mid3, hi, each_n);

		return fu1.get() + fu2.get() + fu3.get() + fu4.get();
	}
	else
	{
		return runge_once_with_n(f, lo, hi, n);
	}
}

