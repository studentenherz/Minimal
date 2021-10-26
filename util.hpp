#if !defined(UTIL)
#define UTIL

#include <boost/numeric/odeint.hpp>
#include <functional>
#include <string>
#include <fstream>
#include <cmath>
#include <stdlib.h>

using namespace std;
using namespace boost::numeric::odeint;

// types
typedef boost::array<double, 6> state_type;			// state given by {rho, theta, z, v_rho, v_theta, v_z}
typedef boost::array<double, 3> vector_type;		// a 3 coordinates vector
typedef function<vector_type(const state_type&, const double)> force_type; // F({x, v}, t)
typedef function<vector_type(const vector_type&, const double)> vector_field_type; // B(x, t)

// Null
const vector_type null_vector = {0, 0, 0};
const state_type null_state = {0, 0, 0, 0, 0, 0};
const force_type null_force = [] (const state_type&, const double) { return null_vector;};
const vector_field_type null_vector_field = [] (const vector_type&, const double) { return null_vector;};


// Operations over types

ostream& operator<<(ostream& out, const state_type& state){
	for(auto x: state)
		out << x << ' ';
	return out; 
}

vector_type get_position(const state_type& state){
	vector_type x;
	for(int i=0; i<3; i++)
		x[i] = state[i];
	return x;
}

vector_type get_velocity(const state_type& state){
	vector_type v;
	for(int i=3; i<6; i++)
		v[i] = state[i];
	return v;
}

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
		vector_type f = F(x, t);	// dispersion term

		dxdt[0] = gam * x[3];															// d(rho)/dt = v_rho
		dxdt[1] = gam * x[4] / x[0];											// d(theta)/dt = v_theta / rho
		dxdt[2] = gam * x[5];															// dz/dt = v_z
		dxdt[3] = f[0] + x[4] * b[2] - x[5] * b[1] + e[0] + gam * x[4] * x[4] / x[0];		// v_rho
		dxdt[3] = f[1] + x[5] * b[0] - x[3] * b[2] + e[1] - gam * x[3] * x[4] / x[0] ;	// v_theta
		dxdt[5] = f[2] + x[3] * b[1] - x[4] * b[0] + e[2];															// v_z
	}
};

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

state_type load_initial_state(string filename){
	state_type initial_state;
	ifstream fi(filename);
	if (fi.is_open())
		for (int i=0; i<(int)initial_state.size(); ++i){
			fi >> initial_state[i];
		}
	else
		cerr << "Unable to open file " << filename << '\n';
	return initial_state; 
}


vector_type B_Asdex(double rp,double zp) {
	const double BT0 = 17.815757116271065;

  double T_a=15.2329;
  double p_a=sqrt(T_a);
  double q_a=p_a/2.0;
  double nu_a=p_a*sqrt(3.0/4.0);

  double cc1=0.4733, cc2=-0.2164,cc3=0.0, cc4=0.0, cc5=0.0, cc6=0.0,cc7=-0.06830, cc8=0.01220, cc9=0.1687;
  double cc10=0.8635, cc11=-1.0682, cc12=0.02166,cc13=-0.002662, cc14=0.1178, cc15=1.4008, cc16=-0.2656,cc17=1.3770, cc18=0.2468;

	double r=rp*0.5;
	double z=zp*0.5;

	double csp=cos(p_a*z);
	double snp=sin(p_a*z);
	double csq=cos(q_a*z);
	double snq=sin(q_a*z);
	double csnu=cos(nu_a*z);
	double snnu=sin(nu_a*z);
	double jb1p=j1(p_a*r);
	double jb1q=j1(q_a*r);
	double jb1nu=j1(nu_a*r);
	double yb1q=y1(q_a*r);
	double yb1nu=y1(nu_a*r);
	double rho=sqrt(r*r+z*z);
	double Br=0.0;

  Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r))/(BT0);


	double jb0p=j0(p_a*r);
	double jb0q=j0(q_a*r);
	double jb0nu=j0(nu_a*r);
	double yb0q=y0(q_a*r);
	double yb0nu=y0(nu_a*r);
	double Bz=0.0;

  Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) )/(BT0);


	double Bt=0.0;
	double u_a=-(cc1*T_a);
	double F0_a= 30.4;
	double Psi= cc1 + cc2*r*r+ r*jb1p*(cc3+cc3*z) + cc5*csp + cc6*snp + r*r*(cc7*csp + cc8*snp) +cc9*cos(p_a*rho) + cc10*sin(p_a*rho) + r*jb1nu*(cc11*csq +cc12*snq) + r*jb1q*(cc13*csnu +cc14*snnu) + r*yb1nu*(cc15*csq + cc16*snq) + r*yb1q*(cc17*csnu+cc18*snnu);
	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r)/(BT0) ;

	vector_type B;
  B[0]=Br;
	B[1]=Bt;
	B[2]=Bz;
	// s_flux[0]=Psi;

  if (rp >= -0.1 && rp < -0.099396) {
    cout << "yb1q = " << yb1q << "c/ r = " << rp << '\n';
  }

  return B;
}

#endif // UTIL
