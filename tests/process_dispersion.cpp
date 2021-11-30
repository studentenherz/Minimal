#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

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
	if(argc < 5){
		cout << "usage: process_dispersion nfiles N dv max_dN\nFormat of input: t x y z vx vy vz\n";
		return 1;
	}

	int nfiles = atoi(argv[1]); 
	int N = atoi(argv[2]); // columns in one file
	double dv = atof(argv[3]);
	int max_dN = atoi(argv[4]);

	// perp and par
	vector<vector<double>> sqr_diffs_par(N);
	vector<vector<double>> sqr_diffs_perp(N);
	vector<vector<double>> v_mod(N);

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

			for(int stop = start; stop < N && stop < start + max_dN; stop++){
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
	}

	// all to one vector for sorting
	vector<pair<double, pair<double, double>>> v(N);

	for(int i=1; i<N; ++i)
		for(int j=0; j<(int)v_mod[i].size(); j++)
			v.push_back({v_mod[i][j], {sqr_diffs_par[i][j] / i, sqr_diffs_perp[i][j] / i}});

	// sort by velocities
	sort(v.begin(), v.end());

	vector<double> mean_sqr_par;
	vector<double> mean_sqr_perp;
	vector<double> mean_v;

	double v_med = dv/2;
	double sum_par = 0, sum_perp = 0;
	int n_sum = 0;

	for(int i=0; i<(int)v.size(); ++i){
		if(v[i].first > v_med + dv/2){
			mean_v.push_back(v_med);
			mean_sqr_par.push_back(n_sum == 0 ? 0 : sum_par / n_sum);
			mean_sqr_perp.push_back(n_sum == 0 ? 0 : sum_perp/ n_sum);

			v_med += dv;
			n_sum = 0;
			sum_par = 0;
			sum_perp = 0;
		}

		sum_par += v[i].second.first;
		sum_perp += v[i].second.second;
		n_sum++;
	}

	mean_v.push_back(v_med);
	mean_sqr_par.push_back(n_sum == 0 ? 0 : sum_par / n_sum);
	mean_sqr_perp.push_back(n_sum == 0 ? 0 : sum_perp/ n_sum);

	ofstream fo("process_disp_out.dat");
	for(int i=0; i<(int)mean_v.size(); ++i)
		fo << mean_v[i] << ' ' << mean_sqr_par[i] << ' ' << mean_sqr_perp[i] << '\n';

	return 0;
}