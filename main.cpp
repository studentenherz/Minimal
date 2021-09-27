#include <iostream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::array<double, 6> state_type;			// state given by {rho, theta, z, v_rho, v_theta, v_z}
typedef boost::array<double, 3> vector_type;		// a 3 coordinates vector

vector_type F(const state_type &x, const double /* t */){ // Force field
	return {0, 0, 1};
}

const double m = 1;

void motion_eq(const state_type &x, state_type &dxdt, const double t ){
	dxdt[0] = x[3];																						// d(rho)/dt = v_rho
	dxdt[1] = x[4];																						// d(theta)/dt = v_theta
	dxdt[2] = x[5];																						// dz/dt = v_z
	dxdt[3] = F(x, t)[1] / m + x[0] * x[4] * x[4];							// v_rho
	dxdt[3] = (- F(x, t)[2] / m +  2 * x[3] * x[4]) / x[0];		// v_theta
	dxdt[5] = F(x, t)[2];																			// v_z
}

void obs_cout(const state_type &x, const double t){
	cout << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
}

int main(int argc, char* argv[]){

	state_type x = {1, 0, 0, 0, 0, 0};

	integrate(motion_eq, x, 0.0, 10.0, 0.01, obs_cout);

	return 0;
}