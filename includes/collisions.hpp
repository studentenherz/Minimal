#if !defined(COLLISIONS_HPP)
#define COLLISIONS_HPP

#define sqr(a) a * a

#include "types.hpp"
#include <vector>
#include <cmath>

const double two_over_sqrt_pi = 2 / sqrt(M_PI);

double erf_minus_d_erf(double x){
	return (erf(x) - x * two_over_sqrt_pi * exp(- sqr(x)));
}

class Collisions{
	vector<int> q; // charges
	vector<double> m; // masses
	vector<double> logl; // log of Lambda
	vector<scalar_field_type> T; // temperatures
	vector<scalar_field_type> n; // concentrations
	double eta; // adimensional constant
	double m_a; // mass of test particle
	double q_a; // charge of test particle
public:
	Collisions (vector<int> q_, vector<double> m_, vector<double> logl_, vector<scalar_field_type> T_, vector<scalar_field_type> n_, double eta_, double ma, double qa): q(q_), m(m_), logl(logl_), T(T_), n(n_), eta(eta_), m_a(ma), q_a(qa) {
		if(!(q.size() == m.size() && q.size() == T.size() && q.size() == n.size()))
			throw length_error("Los vectores deben tener la misma cantidad de elementos, uno para cada especie en orden");
	}

	vector_type slow_down(const state_type& x, double t){
		vector_type r = get_position(x);
		vector_type v = get_velocity(x);
		double v_mod = mod(v);

		double nu = 0;
		for(int i=0; i<(int)q.size(); ++i){
			double xb = v_mod / T[i](r, t);
			nu += eta * sqr(q[i]) * sqr(q_a) * n[i](r, t) *(1 + m_a/m[i]) * logl[i] * erf_minus_d_erf(xb) / (pow(v_mod, 3) * sqr(m_a));
		}

		vector_type sd = - nu * v;
		return sd;
	}
};
  
#endif // COLLISIONS_HPP
