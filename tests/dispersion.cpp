/*
 *  The purpose of this code is to test whether the dispersion part of the collisons
 *	code is working properly or not. For this purpose in the method euler_step of Collision
 *	class the effect of slowing down should be commented to isolate the effects.
 *
*/


#include <fstream>
#include <iostream>
#include <string>

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

int main(int argc, char* argv[]){

	/*
		Dispersion of an Hydrogen ion in an homogeneus plasma by the effect of the electrons
		n = 10^{20} m^{-3}, T = 2.5 keV
	*/

	vector<int> q = {1};
	vector<double> m = {1};
	vector<double> loglambda = {17.5};
	vector<scalar_field_type> T = {Tf};
	vector<scalar_field_type> n = {nf};

	int q_a = 1;
	double m_a = 1836;

	// Motion equation (no fields)
	// MotionEquation eq(gam, null_vector_field, null_vector_field);

	for(unsigned long long seed=1; seed<=50; seed++){

		// Collisions operator
		NormalRand ran(seed); // seed
		Collisions col(q, m, loglambda, T, n, eta, m_a, q_a, ran);

		// initial state
		state_type x = {1.0, 0.0, 0.0, 0.0, 0.0, 0.15};

		string fname = "dispersion_only/onlydv/" + to_string(seed) + ".dat";

		ofstream fo(fname);
		Observer obs(fo);

		cout << "Calculating collisions with seed " << seed << '\n';

		coll_integrate(col, 1, x, 0.0, 1200.0, 0.002, obs, 1);
	}

	cout << "Finished :)\n";

	return 0;
}