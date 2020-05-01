#include "property.h"
#include "using.h"
extern _j0218* j0218;
extern _b1821* b1821;
extern _j1939* j1939;




namespace J0218
{
	struct
	{
		double c1_prime = -(j0218->g1 * (j0218->ratio - 2.0) * j0218->ratio
			+ j0218->g2 * pow(j0218->ratio - 1.0, 2)) * j0218->h2;
		double d1_prime = -(j0218->g1 + j0218->g2) * j0218->h1 * j0218->ratio -
			j0218->g2 * j0218->h2;
	} Const;


	double get_E_parallel(double theta_x, double z_reduced)
	{
		double z = Rlc * z_reduced;
		double res;
		if (z < j0218->h1)
		{
			res = C_E_parallel *
				get_B(theta_x) * ((-j0218->g1 * z + Const.c1_prime)) * z;
		}
		else
		{
			res = C_E_parallel * get_B(theta_x) * (z - j0218->h2) *
				((z + j0218->h2) * j0218->g2 + Const.d1_prime);
		}
		return res > 0.0 ? res : 0.0;
	}
}


// for pulsar B1821
namespace B1821
{
	struct
	{
		double c1_prime = -(b1821->g1 * (b1821->ratio - 2.0) * b1821->ratio
			+ b1821->g2 * pow(b1821->ratio - 1.0, 2)) * b1821->h2;
		double d1_prime = -(b1821->g1 + b1821->g2) * b1821->h1 * b1821->ratio -
			b1821->g2 * b1821->h2;
	} Const;


	double get_E_parallel(double theta_x, double z_reduced)
	{
		double z = Rlc * z_reduced;
		double res;
		if (z < b1821->h1)
		{
			res = C_E_parallel * get_B(theta_x) *
				((-b1821->g1 * z + Const.c1_prime)) * z;
		}
		else
		{
			res = C_E_parallel * get_B(theta_x) * (z - b1821->h2) *
				((z + b1821->h2) * b1821->g2 + Const.d1_prime);
		}
		return res > 0.0 ? res : 0.0;
	}

}

//pulsar J1939
namespace J1939
{
	struct
	{
		double c1_prime = -(j1939->g1 * (j1939->ratio - 2.0) * j1939->ratio
			+ j1939->g2 * pow(j1939->ratio - 1.0, 2)) * j1939->h2;
		double d1_prime = -(j1939->g1 + j1939->g2) * j1939->h1 * j1939->ratio -
			j1939->g2 * j1939->h2;
	} Const;


	double get_E_parallel(double theta_x, double z_reduced)
	{
		double z = Rlc * z_reduced;
		double res;
		if (z < j1939->h1)
		{
			res = C_E_parallel * get_B(theta_x) *
				((-j1939->g1 * z + Const.c1_prime)) * z;
		}
		else
		{
			res = C_E_parallel * get_B(theta_x) * (z - j1939->h2) *
				((z + j1939->h2) * j1939->g2 + Const.d1_prime);
		}
		return res > 0.0 ? res : 0.0;
	}

}

void _j0218::set(const std::vector<double>& params)
{

	this->f = params[0];
	this->g1 = params[1];
	this->ratio = params[2];
	this->_ratio_square = ratio * ratio;
	this->omega = params[3];
	this->g2 = this->g1 * _ratio_square / (1 - _ratio_square);
	this->h2 = f * J0218::Rlc;
	this->h1 = ratio * h2;

	J0218::Const.c1_prime = -(j0218->g1 * (j0218->ratio - 2.0) * j0218->ratio
		+ j0218->g2 * pow(j0218->ratio - 1, 2)) * j0218->h2;
	J0218::Const.d1_prime = -(j0218->g1 + j0218->g2) * j0218->h1 * j0218->ratio -
		j0218->g2 * j0218->h2;
}

void _b1821::set(const std::vector<double>& params) 
{

	this->f = params[0];
	this->g1 = params[1];
	this->ratio = params[2];
	this->_ratio_square = ratio * ratio;
	this->omega = params[3];
	this->g2 = this->g1 * _ratio_square / (1 - _ratio_square);
	this->h2 = f * B1821::Rlc;
	this->h1 = ratio * h2;

	B1821::Const.c1_prime = -(b1821->g1 * (b1821->ratio - 2.0) * b1821->ratio
		+ b1821->g2 * pow(b1821->ratio - 1, 2)) * b1821->h2;
	B1821::Const.d1_prime = -(b1821->g1 + b1821->g2) * b1821->h1 * b1821->ratio -
		b1821->g2 * b1821->h2;

}

void _j1939::set(const std::vector<double>& params) 
{

	this->f = params[0];
	this->g1 = params[1];
	this->ratio = params[2];
	this->_ratio_square = ratio * ratio;
	this->omega = params[3];
	this->g2 = this->g1 * _ratio_square / (1 - _ratio_square);
	this->h2 = f * J1939::Rlc;
	this->h1 = ratio * h2;

	J1939::Const.c1_prime = -(j1939->g1 * (j1939->ratio - 2.0) * j1939->ratio
		+ j1939->g2 * pow(j1939->ratio - 1, 2)) * j1939->h2;
	J1939::Const.d1_prime = -(j1939->g1 + j1939->g2) * j1939->h1 * j1939->ratio -
		j1939->g2 * j1939->h2;

}
