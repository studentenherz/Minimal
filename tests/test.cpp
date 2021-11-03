#include <iostream>
#include <fstream>

#include "util.hpp"
#include "odeint.hpp"

using namespace std;

const double gam = 0.06234974;

vector_type B(const vector_type&x, const double t){
	return B_Asdex(x[0], x[2]);
}

int main(int argc, char* argv[]){
	if(argc < 2){
		cout << "usage: test <input_filename> <output_file(default: out)>\n";
		return 1;
	}

	char* ifname = argv[1];
	char const *ofname = "out";
	if(argc >= 3) ofname = argv[2];

	ifstream fi(ifname);
	ofstream fo(ofname);

	if(!fi.is_open()){ 
		cerr << "Unable to open file " << ifname << "\n";
		return 2;
	}

	double t, r, q, z, vr, vq, vz;
	MotionEquation eq(gam, B);

	state_type dxdy, x;

	while(fi >> t >> r >> q >> z >> vr >> vq >> vz){
		x = {r, q, z, vr, vq, vz};
		eq(x, dxdy, t);
		fo << dxdy << '\n';
	}

	return 0;
}