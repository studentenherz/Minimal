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

vector_type slow_down(const state_type &x){
	double nu, Za, Zb, nb, ma, mb, xb;
	double C1, vsb;

	vector_type v = get_velocity(x);
	double vmod = mod(v);
	xb = vmod / vsb;

	nu = C1 * sqr(Za) * sqr(Zb) * (1 + ma/mb) * erf_minus_d_erf(xb);


	return null_vector;
}
  
#endif // COLLISIONS_HPP
