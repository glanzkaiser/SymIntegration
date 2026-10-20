// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

// Driver code
int main(int argc, char** argv)
{
		// Get starting timepoint
	auto start = high_resolution_clock::now();

	 // Define two complex numbers (Real, Imaginary)
	complex<double> z1(1.0, -1.0);  // 1 - i
	complex<double> z2(0.0, 0.0);  // 0 + 0i

	// Calculate the distance
	double distance = complexdistance(z2, z1);

	// Print result
	cout << "The distance is: " << distance << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}