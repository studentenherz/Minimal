#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>

#include "util.hpp"
#include "odeint.hpp"

using namespace std;

double mean(vector<double> v){
	if( v.size() == 0) return -1;
	double sum = 0;
	for (double x : v) sum += x;
	return sum/v.size();
}

double max(vector<double> v){
	double maxx = 0;
	for(auto x: v) maxx = max(maxx, x);
	return maxx;
}

double standard_deviation(vector<double> v, double mean){
	if( v.size() == 0) return -1;
	double sum = 0;
	for (double x : v) sum += sqr(x - mean);
	return sqrt(sum / v.size());
}

double standard_error(vector<double> v, double mean){
	if( v.size() == 0 || v.size() == 1) return -1;
	double sum = 0;
	for (double x : v) sum += sqr(x - mean);
	return sqrt(sum / (v.size() * (v.size() - 1)));
}

int main(int argc, char* argv[]){
	if(argc < 5){
		cout << "usage: process_dispersion nfiles N dv max_dN max_dV\nFormat of input: t x y z vx vy vz\n";
		return 1;
	}

	int nfiles = atoi(argv[1]); 
	int N = atoi(argv[2]); // columns in one file
	double dv = atof(argv[3]);
	int max_dN = atoi(argv[4]);
	double max_dV = MAXFLOAT;
	if(argc > 5){
		max_dV = atof(argv[5]);
	}

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

			for(int stop = start; stop < N && stop <= start + max_dN; stop++){
				v_sum = v_sum + vs[stop];

				e_v = v_sum / mod(v_sum);

				diff = vs[stop] - vs[start];
				diff_par = dot(diff, e_v);
				diff_perp = mod(cross(diff, e_v));

				if(sqr(diff_par) < max_dV && sqr(diff_perp) < max_dV){
					sqr_diffs_par[stop-start].push_back(sqr(diff_par));
					sqr_diffs_perp[stop-start].push_back(sqr(diff_perp));
					v_mod[stop-start].push_back(mod(v_sum) / (stop - start + 1)); // module of mean velocity
				}
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

	vector<double> std_sqr_par;
	vector<double> std_sqr_perp;
	vector<double> std_v;

	vector<int> ns;
	vector<double> maxv;

	double v_med = dv/2;
	vector<double> par_temp, perp_temp, v_temp;

	for(int i=0; i<(int)v.size(); ++i){
		if(v[i].first > v_med + dv/2 || i == (int)v.size() - 1){
			ns.push_back(v_temp.size());
			maxv.push_back(max(par_temp));


			// means
			double mv = mean(v_temp);
			double mpar = mean(par_temp);
			double mperp = mean(perp_temp);

			// standard deviations
			double stdv = standard_error(v_temp, mv);
			double stdpar = standard_error(par_temp, mpar);
			double stdperp = standard_error(perp_temp, mperp);

			// save them
			mean_v.push_back(mv);
			mean_sqr_par.push_back(mpar);
			mean_sqr_perp.push_back(mperp);

			std_v.push_back(stdv);
			std_sqr_par.push_back(stdpar);
			std_sqr_perp.push_back(stdperp);

			v_med += dv;
			// clear temp vectors
			v_temp.clear();
			par_temp.clear();
			perp_temp.clear();
		}

		v_temp.push_back(v[i].first);
		par_temp.push_back(v[i].second.first);
		perp_temp.push_back(v[i].second.second);
	}

	ofstream fo("process_disp_out.dat");
	for(int i=0; i<(int)mean_v.size(); ++i)
		fo << mean_v[i] << ' ' << mean_sqr_par[i] << ' ' << mean_sqr_perp[i] << 
		' ' <<  std_v[i] << ' ' << std_sqr_par[i] << ' ' << std_sqr_perp[i] << ' ' << ns[i] << ' ' << maxv[i] << '\n';

	return 0;
}