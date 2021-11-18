#if !defined(TYPES)
#define TYPES

#include <functional>
#include <array>
#include <cmath>

using namespace std;

// types
typedef array<double, 6> state_type;		// state given by {rho, theta, z, v_rho, v_theta, v_z}
typedef array<double, 3> vector_type;		// a 3 coordinates vector
typedef function<vector_type(const state_type&, const double)> force_type; // F({x, v}, t)
typedef function<vector_type(const vector_type&, const double)> vector_field_type; // B(x, t)
typedef function<double(const vector_type&, const double)> scalar_field_type; // n(x, t)
// typedef struct {
// 	int Z; // charge number
// 	double m; // mass
// } species_type;

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

ostream& operator<<(ostream& out, const vector_type& vec){
	for(auto x: vec)
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

// Vector operations

double dot(const vector_type& a, const vector_type& b){
	double s = 0;
	for (int i=0; i<(int)a.size(); ++i)
		s += a[i] * b[i];
	return s;
}

vector_type cross(const vector_type& a, const vector_type& b){
	vector_type c;
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
	return c;
}

vector_type operator*(const vector_type& a, double t){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = a[i] * t;
	return c;
}

vector_type operator*(double t, const vector_type& a){
	return a * t;
}

vector_type operator*(const vector_type& a, const vector_type& b){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = a[i] * b[i];
	return c;
}

vector_type operator/(const vector_type& a, double t){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = a[i] / t;
	return c;
}

vector_type operator/(double t, const vector_type& a){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = t / a[i];
	return c;
}

vector_type operator/(const vector_type& a, const vector_type& b){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = a[i] / b[i];
	return c;
}

vector_type operator+(const vector_type& a, const vector_type& b){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = a[i] + b[i];
	return c;
}

vector_type operator-(const vector_type& a, const vector_type& b){
	vector_type c;
	for (int i=0; i<(int)a.size(); ++i)
		c[i] = a[i] - b[i];
	return c;
}

double mod(const vector_type& v){
	return sqrt(dot(v, v));
}

#endif // TYPES
