#include <iostream>
#include <fstream>

#include "util.hpp"
#include "odeint.hpp"

using namespace std;

const double gam = 0.06234974;
const double BT0 = 17.815757116271065; // in Gauss
const double a = 0.5; // in m


vector_type B(const vector_type&x, const double t){
	return B_Asdex(x[0] * a, x[2] * a)/BT0;
}

int main(int argc, char* argv[]){

	state_type x = load_initial_state("initialcond.dat"); // initial state
	MotionEquation motion_eq(gam, B);

	ofstream fo("out.dat");
	Observer obs(fo);
	int steps_rk4nl = integrate(motion_eq, x, 0.0, 2570.0, 0.2, obs);

	cout << x << '\n';
	cout << "Integrator rk4nl took " << steps_rk4nl << " steps\n";

	return 0;
}