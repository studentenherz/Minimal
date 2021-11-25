#include <fstream>
#include <iostream>

#include "collisions.hpp"
#include "odeint.hpp"

using namespace std;

const double v0 = 1.84142e7; // m/s  (3.54 MeV of a proton)
// const double pre_eta = 806039; // m^6 s^-4
const double n0 = 1e20; // m
const double Omega = 0.00085; // s^-1
const double tau = 93e-3; // s
const double a = 0.5; // m

// adimensional constants 
const double gam = 3.42504e6;
const double eta = 1200.55;

// campo de temperatura (en realidad da la velocidad media en relación a v0)
double Tf(const vector_type& v, double t){
	return 1.61029;
}

// campo de densidad (en relación a n0)
double nf(const vector_type& v, double t){
	return 1;
}

	/* Probando el frenamiento de un ion (un protón) por electrones en un plasma
	** con densidad 10^20 1/m³ y temperatura 2.5 keV	
	** que es igual a una velocidad 978702 m/s = 0.053 v0
	*/

vector<int> q = {1};
vector<double> m = {1};
vector<double> loglambda = {17.5};
vector<scalar_field_type> T = {Tf};
vector<scalar_field_type> n = {nf};

int q_a = 1;
double m_a = 1836;

NormalRand ran(2LL);

Collisions col(q, m, loglambda, T, n, eta, m_a, q_a, ran);

vector_type f(const state_type& x, double t){
	return col.slow_down(x, 0) + col.dispersion(x, 0);
}

// Observer Example
class ObserverVelocity{
	ostream &os;
	double thr_v;
	bool check = true;
public:
	ObserverVelocity (ostream &_os, double threshold_v): os(_os), thr_v(threshold_v) {}
	void operator()(const state_type &x, const double t) {
			os << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
			if (check && mod(get_velocity(x)) < thr_v){
				cout << "Got under v = " << thr_v << " at time t = " << t << '\n';\
				check = false;
			}
	}
};

int main(int argc, char* argv[]){

	MotionEquation eq(gam, null_vector_field, null_vector_field, f);

	state_type x = {1.0, 0.0, 0.0, 0.0, 0.0, 0.15};

	cout << f(x, 0);

	ofstream fo("sldwn.dat");
	ObserverVelocity obs(fo, 0.037581);

	integrate(eq, x, 0, 1.2, 0.000001, obs);
	// gives nice result for slowing down

	cout << "Finished :)\n";

	return 0;
}