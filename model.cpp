#include "./src/IC.h"
#include "./src/curvature.h"
#include "./src/runge.h"
#include "./src/sync.h"
#include <fstream>
#include <functional>
#include <future>
#include <iostream>
#include <cmath>
#include <thread>
using namespace std;
extern _j0218* j0218;
extern _b1821* b1821;
extern _j1939* j1939;

namespace J0218
{
	void IC()
	{
		ofstream fout{ "../data/j0218_IC.dat", ios_base::out };

		double lo = log(1.6e-10), hi = log(1.6e-4);
		int N = 50;
		double step = (hi - lo) / (N - 1);
		while (lo < hi)
		{
			double e = exp(lo);
			fout << (e / 1.6e-6) << " " << e * get_IC(e) / pow(dist, 2) << endl;
			lo += step;
		}
		std::cout << "J0218 IC finished" << std::endl;
		fout.close();
	}

	void sync()
	{
		ofstream fout{ "../data/j0218_sync.dat", ios_base::out };

		double lo = log(1.6e-10), hi = log(1.6e-4);
		int N = 30;
		double step = (hi - lo) / N;
		while (lo < hi)
		{
			double e = exp(lo);
			double tmp = e * get_Sync(e) / pow(J0218::dist, 2) / j0218->omega;
			std::cout << "energy: " << (e / 1.6e-6) << "    " << tmp << std::endl;
			fout << (e / 1.6e-6) << " " << tmp << endl;
			lo += step;
		}
		std::cout << "J0218 sync finished" << std::endl;
		fout.close();
	}

	void Cur()
	{
		ofstream fout{ "../data/j0218_cur.dat", ios_base::out };
		double lo = log(1.6e-5), hi = log(4.8e-1);
		int N = 50;
		double step = (hi - lo) / N;
		while (lo < hi)
		{
			double e = exp(lo);
			fout << (e / 1.6e-6) << " " << e * e * get_F_cur(e) / pow(J0218::dist, 2) / j0218->omega << endl;
			lo += step;
		}
		std::cout << "J0218 Curvature finished" << std::endl;

		fout.close();
	}

	namespace _help_total
	{
		double _get_IC(double e)
		{
			return e * get_IC(e) / pow(J0218::dist, 2.0);
		}

		double _get_Sync(double e)
		{
			return e * get_Sync(e) / pow(J0218::dist, 2.0) / j0218->omega;
		}

	}
	void Total()
	{
		ofstream fout{ "../data/j0218_total.dat", ios_base::out };
		double lo = log(1.6e-10), hi = log(4.8e-1);
		int N = 50;
		double step = (hi - lo) / N;
		double e;
		double IC, cur, _sync, total;
		while (lo < hi)
		{
			e = exp(lo);
			std::future<double> fu_IC = std::async(std::launch::async,
				std::ref(_help_total::_get_IC), e);
			std::future<double> fu_Sync = std::async(std::launch::async,
				std::ref(_help_total::_get_Sync), e);
			cur = e * e * get_F_cur(e) / pow(J0218::dist, 2.0) / j0218->omega;
			IC = fu_IC.get();
			_sync = fu_Sync.get();
			total = IC + cur + _sync;
			fout << (e / 1.6e-6) << " " << total << endl;
			lo += step;
		}
		std::cout << "J0218 total finished" << std::endl;

		fout.close();
	}



	void model()
	{
		thread t_cur(Cur);
		thread t_IC(IC);
		thread t_Total(Total);
		thread t_sync(sync);

		t_cur.join();
		t_IC.join();
		t_Total.join();
		t_sync.join();
	}
}


//psr b1821

namespace B1821
{
	void IC()
	{
		ofstream fout{ "../data/b1821_IC.dat", ios_base::out };

		double lo = log(1.6e-10), hi = log(1.6e-4);
		int N = 80;
		double step = (hi - lo) / (N - 1);
		while (lo < hi)
		{
			double e = exp(lo);
			fout << (e / 1.6e-6) << " " << e * get_IC(e) / pow(dist, 2) << endl;
			lo += step;
		}
		std::cout << "B1821 IC finished" << std::endl;
		fout.close();
	}

	void sync()
	{
		ofstream fout{ "../data/b1821_sync.dat", ios_base::out };

		double lo = log(1.6e-10), hi = log(1.6e-4);
		int N = 30;
		double step = (hi - lo) / N;
		while (lo < hi)
		{
			double e = exp(lo);
			double tmp = e * get_Sync(e) / pow(dist, 2) / b1821->omega;
			std::cout << "energy: " << (e / 1.6e-6) << "    " << tmp << std::endl;
			fout << (e / 1.6e-6) << " " << tmp << endl;
			lo += step;
		}
		std::cout << "B1821 sync finished" << std::endl;
		fout.close();
	}

	void Cur()
	{
		ofstream fout{ "../data/b1821_cur.dat", ios_base::out };
		double lo = log(1.6e-5), hi = log(4.8e-1);
		constexpr int N = 50;
		double step = (hi - lo) / N;
		while (lo < hi)
		{
			double e = exp(lo);
			fout << (e / 1.6e-6) << " " << e * e * get_F_cur(e) / pow(dist, 2) 
				/ b1821->omega << endl;
			lo += step;
		}
		std::cout << "B1821 Curvature finished" << std::endl;

		fout.close();
	}

	namespace _help_total
	{
		double _get_IC(double e)
		{
			return e * get_IC(e) / pow(dist, 2.0) ;
		}

		double _get_Sync(double e)
		{
			return e * get_Sync(e) / pow(dist, 2.0) / b1821->omega;
		}
	}
	void Total()
	{
		ofstream fout{ "../data/b1821_total.dat", ios_base::out };
		double lo = log(1.6e-10), hi = log(4.8e-1);
		int N = 100;
		double step = (hi - lo) / N;
		double e;
		double IC, cur, _sync, total;
		while (lo < hi)
		{
			e = exp(lo);
			std::future<double> fu_IC = std::async(std::launch::async,
				std::ref(_help_total::_get_IC), e);
			std::future<double> fu_Sync = std::async(std::launch::async,
				std::ref(_help_total::_get_Sync), e);
			cur = e * e * get_F_cur(e) / pow(dist, 2.0) / b1821->omega;
			IC = fu_IC.get();
			_sync = fu_Sync.get();
			total = IC + cur + _sync;
			fout << (e / 1.6e-6) << " " << total << endl;
			lo += step;
		}
		std::cout << "B1821 total finished" << std::endl;

		fout.close();
	}



	void model()
	{
		thread t_cur(Cur);
		thread t_IC(IC);
		thread t_Total(Total);
		thread t_sync(sync);

		t_cur.join();
		t_IC.join();
		t_Total.join();
		t_sync.join();
	}

}



//psr j1939
namespace J1939
{
	void IC()
	{
		ofstream fout{ "../data/j1939_IC.dat", ios_base::out };

		double lo = log(1.6e-10), hi = log(1.6e-4);
		int N = 50;
		double step = (hi - lo) / (N - 1);
		while (lo < hi)
		{
			double e = exp(lo);
			fout << (e / 1.6e-6) << " " << e * get_IC(e) / pow(dist, 2) << endl;
			lo += step;
		}
		std::cout << "J1939 IC finished" << std::endl;
		fout.close();
	}

	void sync()
	{
		ofstream fout{ "../data/j1939_sync.dat", ios_base::out };

		double lo = log(1.6e-10), hi = log(1.6e-4);
		int N = 30;
		double step = (hi - lo) / N;
		while (lo < hi)
		{
			double e = exp(lo);
			double tmp = e * get_Sync(e) / pow(dist, 2) / j1939->omega;
			std::cout << "energy: " << (e / 1.6e-6) << "    " << tmp << std::endl;
			fout << (e / 1.6e-6) << " " << tmp << endl;
			lo += step;
		}
		std::cout << "J1939 sync finished" << std::endl;
		fout.close();
	}

	void Cur()
	{
		ofstream fout{ "../data/j1939_cur.dat", ios_base::out };
		double lo = log(1.6e-5), hi = log(4.8e-1);
		constexpr int N = 50;
		double step = (hi - lo) / N;
		while (lo < hi)
		{
			double e = exp(lo);
			fout << (e / 1.6e-6) << " " << e * e * get_F_cur(e) / pow(dist, 2)
				/ j1939->omega << endl;
			lo += step;
		}
		std::cout << "J1939 Curvature finished" << std::endl;

		fout.close();
	}

	namespace _help_total
	{
		double _get_IC(double e)
		{
			return e * get_IC(e) / pow(dist, 2.0);
		}

		double _get_Sync(double e)
		{
			return e * get_Sync(e) / pow(dist, 2.0) / j1939->omega;
		}

	}
	void Total()
	{
		ofstream fout{ "../data/j1939_total.dat", ios_base::out };
		double lo = log(1.6e-10), hi = log(4.8e-1);
		int N = 100;
		double step = (hi - lo) / N;
		double e;
		double IC, cur, _sync, total;
		while (lo < hi)
		{
			e = exp(lo);
			std::future<double> fu_IC = std::async(std::launch::async,
				std::ref(_help_total::_get_IC), e);
			std::future<double> fu_Sync = std::async(std::launch::async,
				std::ref(_help_total::_get_Sync), e);
			cur = e * e * get_F_cur(e) / pow(dist, 2.0) / j1939->omega;
			IC = fu_IC.get();
			_sync = fu_Sync.get();
			total = IC + cur + _sync;
			fout << (e / 1.6e-6) << " " << total << endl;
			lo += step;
		}
		std::cout << "J1939 total finished" << std::endl;

		fout.close();
	}



	void model()
	{
		thread t_cur(Cur);
		thread t_IC(IC);
		thread t_Total(Total);
		thread t_sync(sync);

		t_cur.join();
		t_IC.join();
		t_Total.join();
		t_sync.join();
	}

}