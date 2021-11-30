#include <fstream>
#include <iostream>
#include <cmath>
#include <string>

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 2){
		cout << "usage: cil2car nfiles (i.dat)\n Ouputs to ./cart/";
		return 1;
	}

	int n = atoi(argv[1]);

	for(int i=1; i<=n; i++){
		string ifname;
		ifname = to_string(i) + ".dat";
		ifstream fi(ifname);

		if(!fi.is_open()){ 
			cerr << "Unable to open file " << ifname << "\n";
			return 2;
		}

		string ofname("cart/");
		ofname += ifname;
		ofstream fo(ofname);

		double t, r, q, z, vr, vq, vz;
		double x, y, vx, vy;

		while(fi >> t >> r >> q >> z >> vr >> vq >> vz){
				x = r * cos(q);
				y = r * sin(q);
				vx = vr * cos(q) - vq * sin(q);
				vy = vr * sin(q) + vq * cos(q);
				fo << t << ' ' << x << ' ' << y << ' ' << z << ' ' << vx << ' ' << vy << ' ' << vz << '\n';
		}
	}

	return 0;
}