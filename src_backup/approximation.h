#pragma once


inline double approximation(double x)
{
	const double a = 5.4158;
	const double b = 0.85;
	const double c = 0.578;
	const double d = 1.2786;
	if (x < 1000.0) {
		return a * pow(1.0 / x + b, -c * x - 0.33333333) * 1.2533141 * exp(-x - d) * sqrt(x + d) * (1.0 + 55.0 / (72.0 * (x + d)) - 10151.0 / (10368.0 * pow(x + d, 2.0)));
	}
	else {
		return 0;
	}
}
