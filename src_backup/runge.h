#pragma once

#include <functional>
#include <thread>
#include <future>

inline double runge_once(const std::function<double(double)> &f, double lo, double hi, int n)
{
    double k1 = .0, k2, k3, k4;
    double tn = lo;
    double y_n = .0;
    double y = .0;
    const double h = (hi - lo) / n;
    for (int i = 0; i < n; ++i)
    {
        k1 = f(tn);
        k2 = f(tn + h / 2);
        k3 = k2;
        k4 = f(tn + h);
        y_n = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        y = y_n;
        tn += h;
    }
    return y_n;
}


inline double runge(const std::function<double(double)> &f, double lo, double hi, int n)
{
    if (n >= 4)
    {
        double step = (lo + hi) / 4;
        double mid1 = lo + step, mid2 = lo + 2 * step, mid3 = hi - step;
        std::future<double> fu1 = std::async(std::launch::async, runge_once,
                                             f, lo, mid1, n >> 2);

        std::future<double> fu2 = std::async(std::launch::async, runge_once,
                                             f, mid1, mid2, n >> 2);

        std::future<double> fu3 = std::async(std::launch::async, runge_once,
                                             f, mid2, mid3, n >> 2);

        std::future<double> fu4 = std::async(std::launch::async, runge_once,
                                             f, mid3, hi, n >> 2);
        return fu1.get() + fu2.get() + fu3.get() + fu4.get();
    } else
    {
        return runge_once(f, lo, hi, n);
    }

}