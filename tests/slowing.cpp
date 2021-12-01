/*
 *  The purpose of this code is to test whether the slowing down part of the collisons
 *	code is working properly or not. For this purpose in the method euler_step of Collision
 *	class the effect of dispersion down should be commented to isolate the effects.
 *
*/


#include <fstream>
#include <iostream>
#include <string>

#include "collisions.hpp"
#include "odeint.hpp"

using namespace std;

const double v0 = 1.29477e7; // m/s  (3.54 MeV of an alpha)
// const double pre_eta = 806039; // m^6 s^-4
const double n0 = 1e20; // m
const double Omega = 0.00085; // s^-1
const double tau = 93e-3; // s
const double a = 0.5; // m

// adimensional constants 
const double gam = 2.408e6;
const double eta = 3453.5;

// campo de temperatura (en realidad da la velocidad media en relación a v0)
double Tf(const vector_type& v, double t){
	return 2.29035;
}

// campo de densidad (en relación a n0)
double nf(const vector_type& v, double t){
	return 1;
}

int main(int argc, char* argv[]){

	/*
		Slowing down of an Hydrogen ion in an homogeneus plasma by the effect of the electrons
		n = 10^{20} m^{-3}, T = 2.5 keV
	*/

	vector<int> q = {1};
	vector<double> m = {1};
	vector<double> loglambda = {17.5};
	vector<scalar_field_type> T = {Tf};
	vector<scalar_field_type> n = {nf};

	// // alpha
	// int q_a = 2;
	// double m_a = 4 * 1836;

	// // H
	// int q_a = 1;
	// double m_a = 1836;

	// D
	int q_a = 1;
	double m_a = 2 * 1836;


	// Motion equation (no fields)
	MotionEquation eq(gam, null_vector_field, null_vector_field);


	// Collisions operator
	NormalRand ran(1LL); // this is not important in this case
	Collisions col(q, m, loglambda, T, n, eta, m_a, q_a, ran);

	// initial state
	// state_type x = {1.0, 0.0, 0.0, 0.0, 0.0, 1}; // alpha from D+T
	// state_type x = {1.0, 0.0, 0.0, 0.0, 0.0, 0.326015}; // H at 93 keV from NBI
	state_type x = {1.0, 0.0, 0.0, 0.0, 0.0, 0.230527}; // D at 93 keV from NBI

	string fname = "slowing_down_D_93kev.dat";
	ofstream fo(fname);
	Observer obs(fo);

	full_integrate(eq, col, 200, x, 0.0, 3, 0.00001, obs, 20);

	cout << "Finished :)\n";

	return 0;
}