#include <iostream>
#include <fstream>

#include "util.hpp"

using namespace std;

int main(){
	NormalRand gauss_rand(12LL);

	ofstream fo("random.dat");

	// for quick copy and plot in Mathematica
	fo << "{";
	for(int i=1; i<10000; i++)
		fo << gauss_rand() << ", ";
	fo << gauss_rand() << "}";

	return 0;
}