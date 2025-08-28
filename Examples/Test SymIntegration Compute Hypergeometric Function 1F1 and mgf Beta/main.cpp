// For direct computation, the series representation can be used, but care must be taken regarding convergence and numerical stability:
// This direct computation method is suitable for |z| < 1. 
// For other cases, analytic continuation or other representations of the hypergeometric function are needed. 
// Libraries like Boost Math or specialized numerical software are recommended for robust and accurate computation
//  across a wider range of parameters.

#include <cmath>
#include <iostream>
#include <stdexcept>
 #include "symintegrationc++.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	double t = 0.5;
	double alpha = 2.0;
	double beta = 3.0;

	Symbolic z("z");

	int n = 5;
	cout << "\nfor n = " << n << endl;
	cout << "\nHypergeometric_1F1(" << alpha << ", " << beta << ", " << z << ") = " << hypergeometric_1F1(alpha, beta, z, n) << endl;
	cout << "\nHypergeometric_1F1(" << alpha << ", " << beta << ", " << t << ") = " << hypergeometric_1F1(alpha, beta, t, n)  << endl;
	
	/*try 
	{
		Symbolic result = hypergeometric_1F1(alpha, beta, z, n);
		cout << "\nHypergeometric_1F1(" << alpha << ", " << beta << ", " << z << ") = " << result << endl;
	} 
	catch (const exception& e) {
		cerr << "Error: " << e.what() << endl;
	}	*/

	Symbolic result2 = hypergeometric_1F1(alpha, alpha+beta, t, n);
	cout << "\nmgf Beta (" << alpha << ", " << beta << ", " << t << ") = " << result2 << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}