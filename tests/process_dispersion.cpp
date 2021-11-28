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
	double t0 = -1, dt;
	while(fi >> t >> x >> y >> z >> vx >> vy >> vz){
		if(t0 == -1) t0 = t;
		vs.push_back({vx, vy, vz});
	}
	int N = vs.size();
	dt = (t - t0) / N;

	// perp and par
	vector<vector<double>> sqr_diffs_par(N);
	vector<vector<double>> sqr_diffs_perp(N);
	vector<vector<double>> v_mod(N);

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
			v_mod[stop-start].push_back(mod(v_sum) / (stop - start + 1)); // module of mean velocity
		}
	}

	vector<double> mean_par, mean_perp, mean_v_mod;

	vector_mean(sqr_diffs_par, mean_par);
	vector_mean(sqr_diffs_perp, mean_perp);
	vector_mean(v_mod, mean_v_mod);

	ofstream fo("mean_sqr_diff.dat");

	cout << "Output to mean_sqr_diff.dat in the format:\nDt mean_Dv_par mea_Dv_perp mean_v\n";

	for(int i=0; i<N; i++)
		fo << dt * i << ' ' << mean_par[i] << ' ' << mean_perp[i] << ' ' << mean_v_mod[i] << '\n';

	return 0;
}