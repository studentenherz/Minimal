#if !defined(COLLISIONS_HPP)
#define COLLISIONS_HPP

#define sqr(a) a * a

#include <vector>
#include <cmath>

#include "types.hpp"
#include "util.hpp"

const double two_over_sqrt_pi = 2 / sqrt(M_PI);

double erf_minus_d_erf(double x){
	return (erf(x) - x * two_over_sqrt_pi * exp(- sqr(x)));
}

double G(double x){
	return erf_minus_d_erf(x)/(2 * sqr(x));
}

class Collisions{
	vector<int> q; // charges
	vector<double> m; // masses
	vector<double> logl; // log of Lambda
	vector<scalar_field_type> T; // temperatures (actaually velocity, not temperature)
	vector<scalar_field_type> n; // concentrations
	double eta; // adimensional constant
	double m_a; // mass of test particle
	double q_a; // charge of test particle

	NormalRand ran; // random generator

public:
	Collisions (vector<int> q_, vector<double> m_, vector<double> logl_, vector<scalar_field_type> T_, vector<scalar_field_type> n_, double eta_, double ma, double qa, NormalRand ran_): q(q_), m(m_), logl(logl_), T(T_), n(n_), eta(eta_), m_a(ma), q_a(qa), ran(ran_) {
		if(!(q.size() == m.size() && q.size() == T.size() && q.size() == n.size()))
			throw length_error("Los vectores deben tener la misma cantidad de elementos, uno para cada especie en orden");
	}
	Collisions (vector<int> q_, vector<double> m_, vector<double> logl_, vector<scalar_field_type> T_, vector<scalar_field_type> n_, double eta_, double ma, double qa): q(q_), m(m_), logl(logl_), T(T_), n(n_), eta(eta_), m_a(ma), q_a(qa), ran(1LL){
		if(!(q.size() == m.size() && q.size() == T.size() && q.size() == n.size()))
			throw length_error("Los vectores deben tener la misma cantidad de elementos, uno para cada especie en orden");
	}

	// Deterministic slowing down rate
	vector_type slow_down(const state_type& x, double t){
		vector_type r = get_position(x);
		vector_type v = get_velocity(x);
		double v_mod = mod(v);

		double nu = 0;
		for(int i=0; i<(int)q.size(); ++i){
			double xb = v_mod / T[i](r, t);
			// terms that depend on plasma species
			nu +=  sqr(q[i]) *  n[i](r, t) * (1 + m_a/m[i]) * logl[i] * erf_minus_d_erf(xb) ;
		}

		// other terms
		nu *= eta * sqr(q_a) / (pow(v_mod, 3) * sqr(m_a));

		vector_type sd = - nu * v;
		return sd;
	}

	double parallel_dispersion_coeff(const state_type& x, double t){
		vector_type r = get_position(x);
		vector_type v = get_velocity(x);
		double v_mod = mod(v);

		double nu = 0;
		for(int i=0; i<(int)q.size(); ++i){
			double xb = v_mod / T[i](r, t);
			// terms that depend on plasma species
			nu +=  sqr(q[i]) *  n[i](r, t) * logl[i] * G(xb) ;
		}

		// other terms
		nu *= 2 * eta * sqr(q_a) / (pow(v_mod, 3) * sqr(m_a));
		return nu;
	}

	double perpendicular_dispersion_coeff(const state_type& x, double t){
		vector_type r = get_position(x);
		vector_type v = get_velocity(x);
		double v_mod = mod(v);

		double nu = 0;
		for(int i=0; i<(int)q.size(); ++i){
			double xb = v_mod / T[i](r, t);
			// terms that depend on plasma species
			nu +=  sqr(q[i]) *  n[i](r, t) * logl[i] * (erf(xb) - G(xb)) ;
		}

		// other terms
		nu *= 2 * eta * sqr(q_a) / (pow(v_mod, 3) * sqr(m_a));
		return nu;
	}

	// Stochastic dispersion
	vector_type dispersion(const state_type& x, double t){
		vector_type v = get_velocity(x);
		double v_mod = mod(v);

		// nu coefficients
		double nu_parallel = parallel_dispersion_coeff(x, t);
		double nu_perpendicular = perpendicular_dispersion_coeff(x, t);

		// parallel and perpendicular versors
		vector_type e_par = v / v_mod;
		
		// get first perpendicular vector
		vector_type e_perp_1 = cross(e_par, {1, 0, 0});
		if(mod(e_perp_1) == 0)
		e_perp_1 = cross(e_par, {0, 1, 0});

		e_perp_1 = e_perp_1 / mod(e_perp_1); // normalize
		vector_type e_perp_2 = cross(e_perp_1, e_par); // second perpendicular versor


		vector_type dvdt = sqrt(nu_parallel) * v_mod * ran() * e_par + sqrt(nu_perpendicular/2) * v_mod * (ran() * e_perp_1 + ran() * e_perp_2);

		return dvdt;
	}

	// Overall efect of collisions one step 
	void euler_step(state_type &x, const double t, const double dt ){

		vector_type dv = /* slow_down(x, t) * dt + */ dispersion(x, t) * sqrt(dt);

		x[3] += dv[0];
		x[4] += dv[1];
		x[5] += dv[2];
	}

};
  
#endif // COLLISIONS_HPP
