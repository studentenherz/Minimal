#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "util.hpp"
#include "odeint.hpp"

using namespace std;

int main(int argc, char* argv[]){
	if(argc < 3){
		cout << "usage: process_dispersion nfiles N \nFormat of input: t x y z vx vy vz\n";
		return 1;
	}

	int nfiles = atoi(argv[1]); 
	int N = atoi(argv[2]); // columns in one file
	// double dv = atof(argv[3]);

	// perp and par
	vector<double> diffs_par(N);
	vector<double> diffs_perp(N);

	for(int ifile=1; ifile<=nfiles; ++ifile){

		string ifname = to_string(ifile) + ".dat";

		ifstream fi(ifname);
		if(!fi.is_open()) {
			cerr << "Can't open file " << ifname << '\n';
			return 2;
		}

		// vector of velocities
		vector<vector_type> vs;
		double t, x, y, z, vx, vy, vz;
		while(fi >> t >> x >> y >> z >> vx >> vy >> vz){
			vs.push_back({vx, vy, vz});
		}

		vector_type v_sum, e_v, diff;
		double diff_par, diff_perp;

		for(int start = 0; start < N; start++){
			v_sum = vs[start];

			int stop = start + 1;

			v_sum = v_sum + vs[stop];

			e_v = v_sum / mod(v_sum);

			diff = vs[stop] - vs[start];
			diff_par = dot(diff, e_v);
			diff_perp = mod(cross(diff, e_v));

			diffs_par.push_back(diff_par);
			diffs_perp.push_back(diff_perp);
		}
	}

	sort(diffs_par.begin(), diffs_par.end());
	sort(diffs_perp.begin(), diffs_perp.end());

	vector<int> bins_par(), bins_perp;

	ofstream fo("differences.dat");
	for(int i=0; i<(int)diffs_par.size(); ++i)
		fo << diffs_par[i] << ' ' << diffs_perp[i] << '\n';

	return 0;
}