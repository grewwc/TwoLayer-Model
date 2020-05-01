#include "model.h"
#include "src/curvature.h"
#include "src/IC.h"
#include "processData/init_cur.h"
#include "processData/write.h"
#include "src/sync.h"
#include "star.h"
#include "processData/init_sync.h"
#include "wwcString.h"
#include "wwcFit.h"
#include <iomanip>
#include <future>
//#include <Print.h>
#include <memory>
#include <ios>

_j0218 *j0218 = new _j0218;
_b1821 *b1821 = new _b1821;
_j1939 *j1939 = new _j1939;

using namespace std;

std::vector<double> PARAMS {0.2,0.9,0.9,0.2};
std::vector<double> __j0218 {0.33, 0.92, 0.915, 0.507};
std::vector<double> __j1939 {0.29, 0.99, 0.925, 0.25};
std::vector<double> __b1821 {0.24716656, 0.955, 0.92, 0.9};

class _dummy
{
public:
    _dummy()
    {
        j0218->set(PARAMS);
        j1939->set(PARAMS);
        b1821->set(PARAMS);
    }
};

const _dummy dummy;

template<typename T1, typename T2>
bool is_equal(T1 s1, T2 s2);


double f_j0218(double x, const std::vector<double> &params)
{
//    j0218->set(params);
    double x_prime = std::pow(1e5, x) * 1.6e-6;

	double tmp= J0218::get_F_cur(x_prime) * pow(x_prime / J0218::dist, 2) / j0218->omega;
    return tmp;
}

double f_b1821(double x, const std::vector<double> &params)
{
//	b1821->set(params);
    double x_prime = std::pow(1e5, x) * 1.6e-6;
	return B1821::get_F_cur(x_prime) * pow(x_prime / B1821::dist, 2) / b1821->omega;
}

double f_j1939(double x, const std::vector<double> &params)
{
    double x_prime = std::pow(1e5, x) * 1e-6;
	return J1939::get_F_cur(x_prime) * pow(x_prime / J1939::dist, 2) / j1939->omega;
}


void init_data_j0218(std::vector<double>& x, std::vector<double>& y)
{
	double log_lo = log(100.0), log_hi = log(1e5);
	double step = (log_hi - log_lo) / (x.size() - 1);
	double N0 = 1.20328e-11;
	double scale = 735.432312;
	double cutoff = 1988.14;
	double index2 = 1.0;
	double index1 = -1.74313;

    double temp;

	for (int i = 0; i < x.size(); ++i)
	{
        x[i] = (log_lo + i * step) / log_hi;
		temp = exp(x[i] * log_hi);
		y[i] = pow(temp, 2) * N0 * pow(temp / scale, index1) * exp(-temp / cutoff) * 1.6e-6;
	}
}

void init_data_j0218_compare(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z)
{
	double log_lo = log(100.0), log_hi = log(1e5);
	double step = (log_hi - log_lo) / (x.size() - 1);
	double N0 = 1.20328e-11;
	double scale = 735.432312;
	double cutoff = 1988.14;
	double index2 = 1.0;
	double index1 = -1.74313;

    double temp;

	for (int i = 0; i < x.size(); ++i)
	{
        x[i] = (log_lo + i * step) / log_hi;
		temp = exp(x[i] * log_hi);
		y[i] = pow(temp, 2) * N0 * pow(temp / scale, index1) * exp(-temp / cutoff) * 1.6e-6;
        z[i] = f_j1939(x[i], {1,2,3,4});
	}
}


void init_data_b1821(std::vector<double>& x, std::vector<double>& y)
{
	double log_lo = log(100.0), log_hi = log(1e5);
	double step = (log_hi - log_lo) / (x.size() - 1);
	double N0 = 6.074783298e-12;
	double scale = 895.18689;
	double cutoff = 4240.737204;
	double index2 = 1.0;
	double index1 = -1.891028904;
    double temp;

	for (int i = 0; i < x.size(); ++i)
	{
        x[i] = (log_lo + i * step) / log_hi;
		temp = exp(x[i] * log_hi);
		y[i] = pow(temp, 2) * N0 * pow(temp / scale, index1) * exp(-temp / cutoff) * 1.6e-6;
	}
}

void init_data_j1939(std::vector<double>& x, std::vector<double>& y)
{
	double log_lo = log(100.0), log_hi = log(1e5);
	double step = (log_hi - log_lo) / (x.size() - 1);
	double N0 = 5.60746147e-12;
	double scale = 1000;
	double cutoff = 5070.10892;
	double index2 = 1.0;
	double index1 = -2.397103881;
    double temp;

	for (int i = 0; i < x.size(); ++i)
	{
        x[i] = (log_lo + i * step) / log_hi;
		temp = exp(x[i] * log_hi);
		y[i] = pow(temp, 2) * N0 * pow(temp / scale, index1) * exp(-temp / cutoff) * 1.6e-6;
	}
}

string getLastLine(const char* fpath)
{
    ifstream fin {fpath, ios::in};
    if(!fin.is_open())
    {
        cout << "cannot open the file: " << fpath << endl;
        exit(1);
    }
    string line;
    string previous;
    while(getline(fin, line))
    {
        previous = line;
    }
    return previous;
}

decltype(auto) read_from_fit_result(const char* fpath)
{
    auto last_line = getLastLine(fpath);
    if(is_equal(last_line, ""))
    {
        return vector<double>();
    }
    cout << last_line << endl;
    wwcString str {last_line};
    auto res = str.find_all_double();
    std::vector<double> params;
    cout << fpath << " (filename):  " << res.size() << endl;

    if(res.size() == 4)
    {
        for(int i=0; i<res.size(); ++i)
        {
            params.emplace_back(stod(res[i]));
        }
    }
    else 
    {
        for(int i=1; i<res.size(); ++i)
        {
            params.emplace_back(stod(res[i]));
        }
    }

    return params; 
}

double _test(double x, const std::vector<double>& params)
{
    return  params[0] * pow(x, params[1]) + params[2] * x + params[3];
}

void _test_init_data(std::vector<double>& x, std::vector<double>& y)
{
    for (int i=0; i<x.size(); ++i)
    {
        x[i] = (double)i / x.size(); 
        y[i] = 2.1 * pow(x[i], 2.25) + 2.2 * x[i] + 0.9; 
    }
}

template<typename t1, typename t2>
bool is_equal(t1 s1, t2 s2)
{
    string str1, str2;
    try
    {
        str1 = s1;
        str2 = s2;
    } 
    catch(...)
    {
        cout << "cannot transfer to std::string" << endl;
        exit(1);
    }
    return str1 == str2;
}

int main(int argc, char** argv)
{
	std::vector<double> x(20), y(20);

    if(argc < 2)
    {
        cout << "should have 1 parameter " << endl;
    }

    const char* operation = argv[1];
    if(is_equal(operation, "model"))
    {
        cout << "'cur' || 'ic' || 'sync' || 'total' || 'all'" << endl;
        std::string type;
        cin >> type;
        
        //auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        // auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        auto res_j0218 = __j0218;
        auto res_b1821 = __b1821;
        //auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        if(res_j0218.size() != 0)
        {
            j0218->set(res_j0218);
        }
        else
        {
            cout << "file is empty (j0218)" << endl;
        }
        if(res_b1821.size() != 0)
        {
            b1821->set(res_b1821);
        }
        else
        {
            cout << "file is empty (b1821)" << endl;
        }
        if(__j1939.size() != 0)
        {
            j1939->set(__j1939);
        }
        else
        {
            cout << "file is empty (j1939) " << endl;
        }
        if(is_equal(type, "cur"))
        {
            B1821::Cur();
            J0218::Cur();
            J1939::Cur();
        }
        else if(is_equal(type, "IC"))
        {   
            B1821::IC();
            J0218::IC();
            J1939::IC();
        }
        else if(is_equal(type, "sync"))
        {
            B1821::sync();
            J0218::sync(); 
            J1939::sync();
        }
        else if(is_equal(type, "total"))
        {
            J0218::Total();
            B1821::Total();
            J1939::Total();
        }
        else if(is_equal(type, "all"))
        {
            B1821::model();
            J0218::model();
            J1939::model();
        }
        else
            cout << "no such type, choose from 'cur', 'ic', 'sync', 'total', 'all' "<< endl;
    }
    else if (is_equal(operation, "IC_j0218"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        j0218->set(res_j0218); 
        J0218::IC();
    }
    else if(is_equal(operation, "IC_j1939"))
    {
        
        auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        j1939->set(__j1939);
        J1939::IC();
    }
    else if(is_equal(operation, "IC_b1821"))
    {
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        b1821->set(res_b1821);
        B1821::IC();
    }
    else if(is_equal(operation, "sync"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        
        if(res_j0218.size() != 0)
        {
            j0218->set(res_j0218);
        }
        else
        {
            cout << "file is empty (j0218)" << endl;
        }
        if(res_b1821.size() != 0)
        {
            b1821->set(res_b1821);
        }
        else
        {
            cout << "file is empty (b1821)" << endl;
        }
        if(__j1939.size() != 0)
        {
            j1939->set(__j1939);
        }
        else
        {
            cout << "file is empty (j1939) " << endl;
        }

        B1821::sync(); 
        J0218::sync();
        J1939::sync();
    }
   
    else if(is_equal(operation, "IC"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        
        if(res_j0218.size() !=0 )
        {
            j0218->set(res_j0218);
        }
        else
        {
            cout << "file is empty (j0218)" << endl;
        }
        if(res_b1821.size() != 0)
        {
            b1821->set(res_b1821);
        }
        else
        {
            cout << "file is empty (b1821)" << endl;
        }
        if(__j1939.size() != 0)
        {
            j1939->set(__j1939);
        }
        else
        {
            cout << "file is empty (j1939) " << endl;
        }

        B1821::IC(); 
        J0218::IC();
        J1939::IC();
    }
    else if(is_equal(operation, "cur"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        
        if(res_j0218.size() != 0)
        {
            j0218->set(res_j0218);
        }
        else
        {
            cout << "file is empty (j0218)" << endl;
        }
        if(res_b1821.size() != 0)
        {
            b1821->set(res_b1821);
        }
        else
        {
            cout << "file is empty (b1821)" << endl;
        }
        if(__j1939.size() != 0)
        {
            j1939->set(__j1939);
        }
        else
        {
            cout << "file is empty (j1939) " << endl;
        }

        B1821::Cur(); 
        J0218::Cur();
        J1939::Cur();
    }
 
    else if (is_equal(operation, "sync_j0218"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        j0218->set(res_j0218);
        J0218::sync();
    }
    else if(is_equal(operation, "sync_j1939"))
    {
        j1939->set(__j1939); 
        J1939::sync();
    }
    else if(is_equal(operation, "sync_b1821"))
    {
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        b1821->set(res_b1821);
        B1821::sync();
    }
    else if(is_equal(operation, "total"))
    {

        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
 

        if(res_j0218.size() != 0)
        {
            j0218->set(res_j0218);
        }
        else
        {
            cout << "file is empty (j0218)" << endl;
        }
        if(res_b1821.size() != 0)
        {
            b1821->set(res_b1821);
        }
        else
        {
            cout << "file is empty (b1821)" << endl;
        }
        if(__j1939.size() != 0)
        {
            j1939->set(__j1939);
        }
        else
        {
            cout << "file is empty (j1939) " << endl;
        }

        B1821::Total(); 
        J0218::Total();
        J1939::Total();

    }
    else if (is_equal(operation, "total_j0218"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        j0218->set(res_j0218); 
        J0218::Total();
    }
    
    else if (is_equal(operation, "model_j0218"))
    {
        auto res_j0218 = read_from_fit_result("/home/wwc129/test/j0218_best.dat");
        j0218->set(res_j0218); 
        J0218::model();
    }
    else if(is_equal(operation, "total_j1939"))
    {
        j1939->set(__j1939); 
        J1939::Total();
    }
    else if (is_equal(operation, "model_j1939"))
    {
        //auto __j1939 = read_from_fit_result("/home/wwc129/test/j1939_best.dat");
        cout << "here < what "<<endl;
        j1939->set(__j1939); 
        J1939::model();
    }
    else if(is_equal(operation, "total_b1821"))
    {
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        b1821->set(res_b1821);
        B1821::Total();
    }
    else if (is_equal(operation, "model_b1821"))
    {
        auto res_b1821 = read_from_fit_result("/home/wwc129/test/b1821_best.dat");
        b1821->set(res_b1821); 
        B1821::model();
    }
    else if(is_equal(operation, "fit_j0218"))
    {
        cout << "here " << endl;
        std::vector<double> params = {0.28, 0.7, 0.93, 0.15};
        std::vector<double> vec_x(100);
        std::vector<double> vec_y(100);
        
        init_data_j0218(vec_x, vec_y);

        wwcFit(vec_x, vec_y, f_j0218, params, static_cast<int>(1e7), false, "/home/wwc129/test/fit_result_j0218.dat",
            "j0218");
    }
    else if(is_equal(operation, "fit_b1821"))
    {
        std::vector<double> params = {0.28, 0.7, 0.93, 0.15};
        std::vector<double> vec_x(100);
        std::vector<double> vec_y(100);
        
        init_data_b1821(vec_x, vec_y);

        wwcFit(vec_x, vec_y, f_b1821, params, static_cast<int>(1e7), false, "/home/wwc129/test/fit_result_b1821.dat",
            "b1821");
    }
    
    else if(is_equal(operation, "fit_j1939"))
    {
        std::vector<double> params = {0.28, 0.7, 0.93, 0.15};
        std::vector<double> vec_x(100);
        std::vector<double> vec_y(100);
        
        init_data_j1939(vec_x, vec_y);

        wwcFit(vec_x, vec_y, f_j1939, params, static_cast<int>(1e7), false, "/home/wwc129/test/fit_result_j1939.dat", 
            "j1939");
    }
    else if (is_equal(operation, "test"))
    {
        constexpr int N = 100;
        std::vector<double> x(N), y(N), z(N); 
        init_data_j1939(x, y);
        write_data(x.data(), y.data(), N, "/home/wwc129/test/real.dat");
        
        //fit(x, y, _test, params, (int)1e6, true);
    }
    else if (is_equal(operation, "test_fit"))
    {
        _test_init_data(x, y);
        std::vector<double> p {1,1,1,1};
        std::vector<double> params = {0.28, 0.7, 0.9, 0.15};
        std::cout << std::setprecision(20) << f_j0218(0.936842, params) << std::endl;
        params = {0.28, 0.7, 0.15, 0.15};
        std::cout << std::setprecision(20) << f_j0218(0.936842, params) << std::endl;
    } 
    else {
        cout << "no such operation " << endl;
    }
    cout << "finished " << endl;
	return 0;

}
