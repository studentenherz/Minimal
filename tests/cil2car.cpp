#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

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
	double x, y, vx, vy;

	while(fi >> t >> r >> q >> z >> vr >> vq >> vz){
		x = r * cos(q);
		y = r * sin(q);
		vx = vr * cos(q) - vq * sin(q);
		vy = vr * sin(q) + vq * cos(q);
		fo << t << ' ' << x << ' ' << y << ' ' << z << ' ' << vx << ' ' << vy << ' ' << vz << '\n';
	}

	return 0;
}