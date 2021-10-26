#include <iostream>
#include <fstream>
#include "util.hpp"

using namespace std;

const double gam = 0.06234974;

int steps_count;

// The observer
class Observer{
	// ostream os;
public:
	bool verbose = false;

	Observer (bool verb = false): verbose(verb) {steps_count = 0;}
	void operator()(const state_type &x, const double t) {
		if (verbose)
			cout << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
		steps_count++;
	}
};

vector_type B(const vector_type&x, const double t){
	return null_vector;
}

template <typename eq_type>
int integrator(eq_type f, state_type &x, double ti, double tf, double dt /*, void f(const state_type&, const double) */){
	double t = ti;
	int steps = 0;
	while(t < tf){
		steps++;
		// obs(x, t);

		double a[6], b[6], c[6];
		// Valores de las constantes del método integrador
		a[0]=0.0;					      b[0]=0.032918605146;	c[0]=0.0;
		a[1]=-0.737101392796;		b[1]=0.823256998200;	c[1]=0.032918605146;
		a[2]=-1.634740794341;		b[2]=0.381530948900;	c[2]=0.249351723343;
		a[3]=-0.744739003780;		b[3]=0.200092213184;	c[3]=0.466911705055;
		a[4]=-1.469897351522;		b[4]=1.718581042715;	c[4]=0.582030414044;
		a[5]=-2.813971388035;		b[5]=0.27;				    c[5]=0.847252983783;

		state_type dxdt;
		state_type xx; // intermediate states

		dxdt = xx = null_state;

		for(int i=0; i<6; i++){
			double tt = t + c[i] * dt;

			// dxdt = f(x, t)
			f(x, dxdt, tt);

			for(int j=0; j<6; j++){
				xx[j] = a[i] * xx[j] + dt * dxdt[j];
				x[j] = x[j] + b[i] * xx[j];
			}
		}

		t += dt;
	}
	return steps;
}

int main(int argc, char* argv[]){

	state_type x = load_initial_state("initialcond.dat"); // initial state
	MotionEquation motion_eq(gam, B);


	int steps_rk4nl = integrator(motion_eq, x, 0.0, 3.2e5, 0.2);
	cout << x << '\n';
	cout << "Integrator rk4nl took " << steps_rk4nl << " steps\n";


	x = load_initial_state("initialcond.dat"); // initial state
	Observer obs;
	// obs.verbose = true;

	// runge_kutta4<state_type> rk4; //stepper this one gives me the same results as rk4nl
	controlled_runge_kutta<runge_kutta_dopri5<state_type>> rk4;
	integrate_adaptive(rk4, motion_eq, x, 0.0, 3.2e5, 0.2, obs); // integration
	cout << x << '\n';
	cout << "boost odeint took " << steps_count - 1 << " steps\n";

	return 0;
}