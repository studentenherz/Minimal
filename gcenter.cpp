#include <fstream>
#include <iostream>

#include "util.hpp"

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
	vector_type x, v, b, gc;

	while(fi >> t >> r >> q >> z >> vr >> vq >> vz){
		x = {r, q, z};
		v = {vr, vq, vz};
		b = B(x, t);
		gc = g_center(x, v, b, gam);
		fo << gc << '\n';
	}

	return 0;
}