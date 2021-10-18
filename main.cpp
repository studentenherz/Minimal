#include <iostream>

#include "util.hpp"

using namespace std;

// The force field
vector_type force(const state_type &x, const double /* t */){
	vector_type B = B_Asdex(x[0], x[3]);
	vector_type F;
	F[0] = x[4] * B[2] - x[5] * B[1];
	F[1] = x[5] * B[0] - x[3] * B[2];
	F[2] = x[3] * B[1] - x[4] * B[0];
	return F;
}

// The observer
void obs_cout(const state_type &x, const double t){
	cout << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
}

int main(int argc, char* argv[]){

	state_type x = load_initial_state("initialcond.dat"); // initial state
	MotionEquation motion_eq(force); // motion equation using the 'force'

	integrate(motion_eq, x, 0.0, 10.0, 0.01, obs_cout); // integration

	return 0;
}