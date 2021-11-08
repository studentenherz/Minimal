#if !defined(COLLISIONS_HPP)
#define COLLISIONS_HPP

#include "types.hpp"
#include <vector>
#include <cmath>

class Collisions{
	vector<species_type> plasma_species;
	species_type test_particle;
public:
	Collisions();
	vector_type slowing_down(const vector_type& v){
		vector_type dvdt = null_vector;
		for (species_type s : plasma_species){
			dvdt += nu(s, test_particle);
		}
	}
private:
	double nu(species_type s, species_type test){
		double nu_0 = 4 * M_PI
	}
};

#endif // COLLISIONS_HPP
