#include <iostream>
#include <fstream>
#include <vector>

#include "util.hpp"
#include "odeint.hpp"

using namespace std;

void vector_mean(const vector<vector<double>>& vs, vector<double>& means){
	for(auto v : vs){
		double sum = 0;
		for (double x : v) sum += x;
		means.push_back(sum/v.size());
	}
}

int main(int argc, char* argv[]){
	if(argc < 2){
		cout << "usage: process_dispersion <input_file> \nFormat of input: t x y z vx vy vz\n";
		return 1;
	}

	char* ifname = argv[1];

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
	int N = vs.size();

	// perp and par
	vector<vector<double>> sqr_diffs_par(N);
	vector<vector<double>> sqr_diffs_perp(N);

	vector_type v_sum, e_v, diff;
	double diff_par, diff_perp;

	for(int start = 0; start < N; start++){
		v_sum = vs[start];

		for(int stop = start; stop < N; stop++){
			v_sum = v_sum + vs[stop];

			e_v = v_sum / mod(v_sum);

			diff = vs[stop] - vs[start];
			diff_par = dot(diff, e_v);
			diff_perp = mod(cross(diff, e_v));

			sqr_diffs_par[stop-start].push_back(sqr(diff_par));
			sqr_diffs_perp[stop-start].push_back(sqr(diff_perp));
		}
	}

	vector<double> mean_par, mean_perp;

	vector_mean(sqr_diffs_par, mean_par);
	vector_mean(sqr_diffs_perp, mean_perp);

	ofstream fo("mean_sqr_diff.dat");

	for(int i=0; i<N; i++)
		fo << mean_par[i] << ' ' << mean_perp[i] << '\n';

	return 0;
}