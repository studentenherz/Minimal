#if !defined(TYPES)
#define TYPES

#include <functional>
#include <array>

using namespace std;

// types
typedef array<double, 6> state_type;		// state given by {rho, theta, z, v_rho, v_theta, v_z}
typedef array<double, 3> vector_type;		// a 3 coordinates vector
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


#endif // TYPES
