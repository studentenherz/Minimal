#include <iostream>

#include "util.hpp"

using namespace std;


const double m = 1;

// The force field
vector_type force(const state_type &x, const double /* t */){
	return {0, 0, 1};
}

// The observer
void obs_cout(const state_type &x, const double t){
	cout << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
}

int main(int argc, char* argv[]){

	state_type x = {1, 0, 0, 0, 0, 0}; // initial state

	MotionEquationCylindrical motion_eq(force); // motion equation using the 'force'

	integrate(motion_eq, x, 0.0, 10.0, 0.01, obs_cout); // integration

	return 0;
}