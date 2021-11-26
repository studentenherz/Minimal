#if !defined(ODEINT)
#define ODEINT

#include <fstream>
#include <stdlib.h>

#include "types.hpp"

using namespace std;

class MotionEquation{
	const double gam;			// adimensional factor
	vector_field_type B;	// magnetic induction field
	vector_field_type E;	// electric field
	force_type F;					// other forces
public:
	MotionEquation(const double _gam, vector_field_type _B = null_vector_field, vector_field_type _E = null_vector_field, force_type _F = null_force): gam(_gam), B(_B), E(_E), F(_F) {}
	void operator()(const state_type &x, state_type &dxdt, const double t ){

		vector_type r = get_position(x);
		vector_type b = B(r, t);
		vector_type e = E(r, t);
		vector_type f = F(x, t);	// other forces

		dxdt[0] = gam * x[3];															// d(rho)/dt = v_rho
		dxdt[1] = gam * x[4] / x[0];											// d(theta)/dt = v_theta / rho
		dxdt[2] = gam * x[5];															// dz/dt = v_z
		dxdt[3] = f[0] + x[4] * b[2] - x[5] * b[1] + e[0] + gam * x[4] * x[4] / x[0];		// v_rho
		dxdt[4] = f[1] + x[5] * b[0] - x[3] * b[2] + e[1] - gam * x[3] * x[4] / x[0];		// v_theta
		dxdt[5] = f[2] + x[3] * b[1] - x[4] * b[0] + e[2];															// v_z
	}
};


// template<typename motion_eq_type, typename collision_op_type>
// class Evolution{
// 	motion_eq_type eq;
// 	collision_op_type col;
// public:
// 	Evolution(motion_eq_type motion_eq, collision_op_type collision_operator): eq(motion_eq), col(collision_operator) {}
// 	void operator()(const state_type &x, state_type &dxdt, const double t ){
// 		eq(x, dxdt, t);
// 		col(x, dxdt, t);
// 	}
// };

// Integrator RK4_NL
template <typename eq_type, typename obs_type>
int integrate(eq_type f, state_type &x, double ti, double tf, double dt , obs_type obs, int obs_interval = 1){
	cout << "===================================\n"
			 << "Integrating:\nti = " << ti << "\ntf = " << tf << "\ndt = " << dt
			 << "\nrecording every " << obs_interval << " steps\n";
	double t = ti;
	int steps = 0;
	while(t < tf){
		if (steps % obs_interval == 0) obs(x, t);

		double a[6], b[6], c[6];
		// Valores de las constantes del método integrador
		a[0]=0.0;								b[0]=0.032918605146;	c[0]=0.0;
		a[1]=-0.737101392796;		b[1]=0.823256998200;	c[1]=0.032918605146;
		a[2]=-1.634740794341;		b[2]=0.381530948900;	c[2]=0.249351723343;
		a[3]=-0.744739003780;		b[3]=0.200092213184;	c[3]=0.466911705055;
		a[4]=-1.469897351522;		b[4]=1.718581042715;	c[4]=0.582030414044;
		a[5]=-2.813971388035;		b[5]=0.27;						c[5]=0.847252983783;

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
		steps++;
	}
	return steps;
}

// NULL Observer
void null_observer(const state_type &x, const double t){
	return;
}

// Overloading for NULL Observer
template <typename eq_type>
int integrate(eq_type f, state_type &x, double ti, double tf, double dt){
	return integrate(f, x, ti, tf, dt, null_observer);
}

// Observer Example
class Observer{
	ostream &os;
public:
	Observer (ostream &_os): os(_os) {}
	void operator()(const state_type &x, const double t) {
			os << t << '\t' << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << ' ' << x[4] << ' ' << x[5] << '\n';
	}
};


// Integrate particles 
template <typename eq_type, typename obs_type, typename collisions_type>
int full_integrate(eq_type f, collisions_type col, int nskip,  state_type &x, double ti, double tf, double dt , obs_type obs, int obs_interval = 1){
	cout << "==============COLLISIONS TOO=====================\n"
			 << "Integrating:\nti = " << ti << "\ntf = " << tf << "\ndt = " << dt
			 << "\nrecording every " << obs_interval << " steps\n";
	double t = ti;
	int steps = 0;
	while(t < tf){

		if (steps % obs_interval == 0) obs(x, t);

		double a[6], b[6], c[6];
		// Valores de las constantes del método integrador
		a[0]=0.0;								b[0]=0.032918605146;	c[0]=0.0;
		a[1]=-0.737101392796;		b[1]=0.823256998200;	c[1]=0.032918605146;
		a[2]=-1.634740794341;		b[2]=0.381530948900;	c[2]=0.249351723343;
		a[3]=-0.744739003780;		b[3]=0.200092213184;	c[3]=0.466911705055;
		a[4]=-1.469897351522;		b[4]=1.718581042715;	c[4]=0.582030414044;
		a[5]=-2.813971388035;		b[5]=0.27;						c[5]=0.847252983783;

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

		if (steps % nskip == 0)
			col.euler_step(x, t, dt * nskip);

		t += dt;
		steps++;
	}
	return steps;
}

#endif // ODEINT