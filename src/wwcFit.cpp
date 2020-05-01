#include "Print.h"
#include "wwcString.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <future>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>
#include "property.h"

extern _j0218* j0218;
extern _j1939* j1939;
extern _b1821* b1821;

_star* star;

#ifdef __linux__
inline bool mynan(double test)
{
    return isnan(test);
}
inline bool myinf(double test)
{
    return isinf(test);
}
#elif _WIN32
inline bool mynan(double test)
{
    return isnan<double>(test);
}
inline bool myinf(double test)
{
    return isnan<double>(test);
}
#elif __APPLE__
inline bool mynan(double test)
{
    return isnan<double>(test);
}
inline bool myinf(double test)
{
    return isnan<double>(test);
}
#endif

constexpr double MAX_ALLOWED_CHANGE_ONCE = 1e-1;
constexpr double MIN_ALLOWED_CHANGE_ONCE = 1e-3;

double alpha = 1.0;


namespace log_helper {
std::ofstream f_log{ "/home/wwc129/Twolayer/log/fit.log", std::ios::out };
void _init()
{
    if (!f_log.is_open()) {
        std::cout << "cannot write to the log file" << std::endl;
        exit(-1);
    }
}
}
int N;

double inverse_N = 1.0;
const double err = 1e-14;

std::vector<double> velocity;

namespace _helper {
inline double find_max(const std::vector<double>& params)
{
    double ret = std::numeric_limits<double>::min();
    for (const auto& e : params) {
        ret = ret <= std::abs(e) ? std::abs(e) : ret;
    }
    return ret;
}


inline std::tuple<int, double> find_max_percentage(const std::vector<double>& params)
{
    std::tuple<int, double> ret = std::make_tuple(-1, std::numeric_limits<double>::min());
    int first;
    double second;
    double each_percent = 0.0;

    for (first = 0; first < params.size(); ++first) {
        each_percent = std::abs(velocity[first]) / std::abs(params[first]);
        second = std::get<1>(ret);
        if (second < each_percent) {
            std::get<1>(ret) = each_percent;
            std::get<0>(ret) = first;
        }
    }
    return ret;
}

inline double find_min(const std::vector<double>& params)
{
    double min = std::numeric_limits<double>().max();
    for (const auto& e : params) {
        if (std::abs(e) <= min) {
            min = std::abs(e);
        }
    }
    return min;
}
}

inline double sum(const std::vector<double>& x,
    const std::vector<double>& y,
    double f(double, const std::vector<double>&),
    const std::function<double(double, std::vector<double>&)>& f_1,
    std::vector<double>& params)
{
    
    if(params[0] == 0)
    {
        cout << "wrong here" << endl;
    }
    

    double res = 0;
    for (int i = 0; i < N; i++) {
        res += (f(x[i], params) - y[i]) * f_1(x[i], params);
    }
    return res * inverse_N * 0.5;
}

class derivative {
private:
    std::function<double(double, std::vector<double>&)> _f;

public:
    std::vector<std::function<double(double, std::vector<double>&)>> all_func;

    derivative(double f(double, const std::vector<double>&),
        std::vector<double>& params)
        : _f{ f }
    {

        for (int i = 0; i < params.size(); ++i) {
            all_func.emplace_back([i, this](double x0, std::vector<double>& parameter) {
                double delta = 1e-4* (parameter[i] == 0.0 ? 1e-4 : parameter[i]);
                double y0 = _f(x0, parameter);
                parameter[i] += delta;
                star->set(parameter);
                
                double y1 = _f(x0, parameter);
               
               parameter[i] -= delta;
                star->set(parameter);
                
                double ret = (y1 - y0) / delta;
                
               /* 
                if(i == 2) {
                    std::cout << "third: " << ret <<  "y0, y1 : " << y0 << y1 << std::endl;
                }
                if(i == 3) {
                    std::cout << "fourth: " << ret << std::endl;
                }

                */

                if (mynan(ret)) {
                    std::ostringstream oss(std::ostringstream::ate);
                    oss << "In function '" << __FUNCTION__
                        << " ' the derivative is NaN." << std::endl
                        << "y1 - y0: " << (y1 - y0) << "  delta: " << delta << std::endl;
                    oss << "y1: " << y1 << " y0: " << y0 << std::endl;
                    oss << "x0: " << x0 << std::endl;
                    for_each(parameter.begin(), parameter.end(), [&oss](double i) {
                        oss << i << std::endl;
                    });
                    Print::frame(oss.str(), log_helper::f_log);
                    exit(-1);
                }

                if (myinf(ret)) {
                    log_helper::f_log << "In function '" << __FUNCTION__
                                      << " ' the derivative is 'inf'." << std::endl;
                    exit(-1);
                }
                return ret;
            });
        }
    }
};

inline double get_fit_once_err(const std::vector<double>& x,
    const std::vector<double>& y,
    double f(double, const std::vector<double>&),
    const std::vector<double>& params, int N)
{
    double this_err = 0.0;

    for (int i = 0; i < N; ++i) {
        this_err += std::pow(f(x[i], params) - y[i], 2);
    }
    return sqrt(this_err / (double) N);
}

inline void limit_fit_change_at_one_time(const std::vector<double>& params, double& alpha)
{
    double second, ratio;
    std::tuple<int, double> check_maximum_change;
    do {
        check_maximum_change = _helper::find_max_percentage(params);
        second = std::get<1>(check_maximum_change);


        if (second > MAX_ALLOWED_CHANGE_ONCE) {
            ratio = MAX_ALLOWED_CHANGE_ONCE / second;
            alpha *= ratio;
            for_each(velocity.begin(), velocity.end(), [ratio](double& i) {
                i *= ratio;
            });
        } else if (second < MIN_ALLOWED_CHANGE_ONCE) {
            ratio = MIN_ALLOWED_CHANGE_ONCE / second;
            alpha *= ratio;
            for_each(velocity.begin(), velocity.end(), [ratio](double& i) {
                i *= ratio;
            });
        }
    } while (second > MAX_ALLOWED_CHANGE_ONCE || second < MIN_ALLOWED_CHANGE_ONCE);

    std::cout << "now, velocity: " ;
    for(int i=0; i<velocity.size(); ++i){
        std::cout << velocity[i] << "   ";
    }
    std::cout << std::endl;
}


inline double fit_once(const std::vector<double>& x,
    const std::vector<double>& y,
    double f(double, const std::vector<double>&),
    std::vector<double>& params, const derivative& all_derivative_func, bool any)
{
    auto fu = std::async(std::launch::async, std::ref(get_fit_once_err), x, y,
        std::ref(f), std::ref(params), N);
    double record_sum;
    double _max_; //for change alpha if alpha is too small
    for (int i = 0; i < params.size(); ++i) {
        record_sum = sum(x, y, f, all_derivative_func.all_func[i], params);

        if (alpha == 1.0) {
            std::cout << (all_derivative_func.all_func[i])(1, params) << std::endl;
            std::cout << "record_sum: " << record_sum << std::endl;
        }

        if (mynan(record_sum)) {
            std::ostringstream oss{ std::ostringstream::app };
            oss << "record_sum is NaN, \n"
                << "line number: " << __LINE__ << std::endl;
            Print::frame(oss.str(), log_helper::f_log);
            exit(-1);
        }

        velocity[i] = alpha * record_sum + 0.9 * velocity[i];
    }

    //just for bound the parameters
    if (alpha == 1.0) {
        double change = sqrt(std::accumulate(velocity.begin(), velocity.end(), 0.0,
            [](double x, double y) {
                return x + pow(y, 2.0);
            }));
        alpha = 0.1 * _helper::find_min(params) / change;
    }
    
    
    
    limit_fit_change_at_one_time(params, alpha);
    
    
    
    for (int i = 0; i < params.size(); ++i) {
        params[i] -= 0.1* velocity[i];

        if (any == false) {
            switch (i) {
            case 0:
                if (params[i] < 0.1)
                    params[i] = 0.1;
                else if (params[i] > 0.45)
                    params[i] = 0.45;
                break;
            case 3:
                if (params[i] < 0.05)
                    params[i] = 0.05;
                else if (params[i] > 2)
                    params[i] = 2;
                break;
            default:
                if (params[i] < 0.25)
                    params[i] = 0.25;
                else if (params[i] > 0.995)
                    params[i] = 0.995;
                else if (mynan(params[i])) {
                    std::ostringstream oss{ std::ostringstream::ate };
                    oss << "params: " << params[i] << " is NaN" << std::endl;
                    Print::frame(oss.str(), log_helper::f_log);
                }

                break;
            }
        }
    }
    star->set(params);

    return fu.get();
}

using std::cout;
using std::endl;
using std::string;

void wwcFit(const std::vector<double>& x,
    const std::vector<double>& y,
    double f(double, const std::vector<double>&),
    std::vector<double>& params, int max_time, bool any, const char* write_to_file, const char* psr_name)
{

    //determine which pulsar
    if(string(psr_name) == string("j0218"))
    {
        cout << "j0218 is beginning" << endl;
        star = dynamic_cast<_j0218*>(j0218);
    }
    else if(string(psr_name) == string("b1821"))
    {
        cout << "b1821 is beginning" << endl;
        star = dynamic_cast<_b1821*>(b1821);
    }
    else if(string(psr_name) == string("j1939"))
    {
        cout << "j1939 is beginning" << endl;
        star = dynamic_cast<_j1939*>(j1939);
    }
       
    N = x.size();
    inverse_N = 1.0 / N;

    std::ofstream fout{ write_to_file, std::ios_base::out };
    
//Begin
    //frame the initial parameters
    std::ostringstream oss;
    oss << "Initial Pamaters: " << endl;
    for (const auto& e : params)
    {
        oss << e << "  ";
    }
    oss << endl;
    Print::frame(oss.str());
//End


    for (int i = 0; i < params.size(); ++i) {
        velocity.emplace_back(0.0);
    }
    
    std::cout << "size: " << velocity.size() << std::endl;

    derivative all_derivative_func(f, params);
    unsigned fit_times = 0;
    double res_fit_once;
    do {
        
        res_fit_once = fit_once(x, y, f, params, all_derivative_func, any);
        if (fit_times % 100 == 0) {
            std::cout << "process: %" << (float)fit_times / max_time * 100 << std::endl;
        }

        fout << fit_times << " :     ";
        fout.precision(8);
        fout.setf(std::ios_base::fixed, std::ios_base::floatfield);
        for (int i = 0; i < params.size(); ++i) {
            fout << params[i] << " ,  ";
        }
        fout << "err (sqrt / N):  " << std::setprecision(20) << std::scientific << res_fit_once << std::endl;
        fout.flush();
    } while (++fit_times < max_time && res_fit_once > err);

    fout.close();
}
