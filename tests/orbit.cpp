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
	if(argc < 2){
		cout << "usage: orbit <input_filename> <output_file(default: out)> nskip\n";
		return 1;
	}

	char* ifname = argv[1];
	char const *ofname = "out";
	if(argc >= 3) ofname = argv[2];


	state_type x = load_initial_state(ifname); // initial state
	MotionEquation motion_eq(gam, B);

	ofstream fo(ofname);
	Observer obs(fo);
	
	integrate(motion_eq, x, 0.0, 2570.0, 0.2, obs);

	cout << "Done :)\n";

	return 0;
}