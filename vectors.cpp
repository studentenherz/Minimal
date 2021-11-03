#include <iostream>
#include "types.hpp"

using namespace std;

int main(){
	vector_type a = {1, 1, 1};
	vector_type b = {0, 3, 4};

	vector_type c = a * b;

	cout << c << '\n';
	cout << 1/c/3 << '\n';
	
	return 0;
}