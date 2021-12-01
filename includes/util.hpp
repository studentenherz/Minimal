#if !defined(UTIL)
#define UTIL

#define sqr(a) (a) * (a)
#define two_pi 6.283185307179586477

#include <string>
#include <fstream>
#include <cmath>

#include "types.hpp"


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


vector_type B_Asdex(double r,double z) {
	// Gives ASDEX-U magnetic field for discharge #10958
	// this is in cgs and real units

	double T_a=15.2329;
	double p_a=sqrt(T_a);
	double q_a=p_a/2.0;
	double nu_a=p_a*sqrt(3.0/4.0);

	double cc1=0.4733, cc2=-0.2164,cc3=0.0, cc4=0.0, cc5=0.0, cc6=0.0,cc7=-0.06830, cc8=0.01220, cc9=0.1687;
	double cc10=0.8635, cc11=-1.0682, cc12=0.02166,cc13=-0.002662, cc14=0.1178, cc15=1.4008, cc16=-0.2656,cc17=1.3770, cc18=0.2468;

	// double r=rp*0.5;
	// double z=zp*0.5;

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

	Br=(-(r*jb1p*cc4-cc5*p_a*snp+cc6*p_a*csp+r*r*p_a*( -cc7*snp+cc8*csp ) - cc9*p_a*sin(p_a*rho)*(z/rho) + cc10*p_a*cos(p_a*rho)*(z/rho)+ r*jb1nu*(-q_a*cc11*snq+cc12*q_a*csq) +r*jb1q*( -cc13*nu_a*snnu +cc14*nu_a*csnu ) + r*yb1nu*( -cc15*q_a*snq+cc16*q_a*csq) + r*yb1q*(-nu_a*cc17*snnu + cc18*nu_a*csnu))*(1.0/r));


	double jb0p=j0(p_a*r);
	double jb0q=j0(q_a*r);
	double jb0nu=j0(nu_a*r);
	double yb0q=y0(q_a*r);
	double yb0nu=y0(nu_a*r);
	double Bz=0.0;

	Bz=(( 2.0*cc2*r  + jb1p*( cc3 + z*cc4) +r*(cc3+cc4*z)*( p_a*jb0p-(jb1p/r) ) +2.0*r*( cc7*csp+cc8*snp) - cc9*sin(p_a*rho)*((p_a*r)/rho) + cc10*cos(p_a*rho)*((p_a*r)/rho) + jb1nu*(cc11*csq+cc12*snq) +r*(cc11*csq +cc12*snq)*( nu_a*jb0nu-(jb1nu/r) ) + jb1q*(cc13*csnu + cc14*snnu) + r*(cc13*csnu + cc14*snnu)*( q_a*jb0q-(jb1q/r) ) + yb1nu*(cc15*csq+cc16*snq) +r*(cc15*csq+cc16*snq)*( nu_a*yb0nu-(yb1nu/r) ) + yb1q*( cc17*csnu +cc18*snnu) + r*( cc17*csnu +cc18*snnu)*( q_a*yb0q-(yb1q/r) )   )*(1.0/r) );


	double Bt=0.0;
	double u_a=-(cc1*T_a);
	double F0_a= 30.4;
	double Psi= cc1 + cc2*r*r+ r*jb1p*(cc3+cc3*z) + cc5*csp + cc6*snp + r*r*(cc7*csp + cc8*snp) +cc9*cos(p_a*rho) + cc10*sin(p_a*rho) + r*jb1nu*(cc11*csq +cc12*snq) + r*jb1q*(cc13*csnu +cc14*snnu) + r*yb1nu*(cc15*csq + cc16*snq) + r*yb1q*(cc17*csnu+cc18*snnu);
	Bt= ((sqrt(T_a*Psi*Psi+2.0*u_a*Psi+ F0_a*F0_a ))/r) ;

	vector_type B;
	B[0]=Br;
	B[1]=Bt;
	B[2]=Bz;
	// s_flux[0]=Psi;

	// if (rp >= -0.1 && rp < -0.099396) {
	// 	cout << "yb1q = " << yb1q << "c/ r = " << rp << '\n';
	// }

	return B;
}

vector_type g_center(const vector_type& r, const vector_type& v,const vector_type& B, double gam){
	double modB = mod(B);
	vector_type dr = gam * cross(v, B) / (modB * modB);
	return r + dr;
}

// ============== Rand 2 =======================
// implementation from Numerical Recipes Ed3
class Ran2
{
	unsigned long long u, v, w;

public:
	Ran2(unsigned long long j) : v(4101842887655102017LL), w(1){
		u = j ^ v;
		int64();
		v = u;
		int64();
		w = v;
		int64();
	}

	inline unsigned long long int64(){
		u = u * 2862933555777941757LL + 7046029254386353087LL;
		v ^= v >> 17;
		v ^= v << 31;
		v ^= v >> 8;
		w = 4294957665U * (w & 0xffffffff) + (w >> 32);
		unsigned long long x = u ^ (u << 21);
		x ^= x >> 35;
		x ^= x << 4;
		return (x + v) ^ w;
	}
	inline unsigned int int32() { return (unsigned int)int64(); }
	inline double doub() { return 5.42101086242752217E-20 * int64(); }
	double random(double xmin, double xmax){
		return xmin + doub() * (xmax - xmin);
	}
};

class NormalRand{
	double sigma;	// std deviation
	double mu;		// mean
	unsigned long long seed;
	Ran2 ran;			// uniform random generator
public:
	NormalRand(unsigned long long seed_, double sigma_ = 1.0, double mu_ = 0): sigma(sigma_), mu(mu_), seed(seed_), ran(seed) {};
	
	double operator()(){
		// This uses Box-Muller algorithm
		double x1, x2;
		x1 = ran.doub(); // 0 to 1
		x2 = ran.doub(); // 0 to 1
		return sigma * sqrt(-2.0 * log(x1)) * cos(two_pi * x2) + mu;
	}
};

#endif // UTIL
