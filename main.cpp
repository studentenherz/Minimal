#include <iostream>
#include <fstream>
#include "util.hpp"

using namespace std;

const double gam = 0.06234974;

int steps_count;

// The observer
class Observer{
	// ostream os;
public:
	bool verbose = false;

	Observer (bool verb = false): verbose(verb) {steps_count = 0;}
	void operator()(const state_type &x, const double t) {
		if (verbose)
			cout << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
		steps_count++;
	}
};

vector_type B(const vector_type&x, const double t){
	return null_vector;
}

int main(int argc, char* argv[]){

	state_type x = load_initial_state("initialcond.dat"); // initial state
	MotionEquation motion_eq(gam, B);


	int steps_rk4nl = integrate(motion_eq, x, 0.0, 3.2e5, 0.2);
	cout << x << '\n';
	cout << "Integrator rk4nl took " << steps_rk4nl << " steps\n";

	return 0;
}