// Merci beaucoup Freya et Sentinel
// g++ main.cpp -o result -lsymintegration
#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include "symintegrationc++.h"

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f
#define pi  3.1415926535897

#include <chrono>
#include <string>
#include <bitset>

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>


using namespace std::chrono;
using namespace std;

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	Symbolic x("x"), f;
	f = pow(x, Symbolic(3)) + x - 1;	

	// Set output precision
	cout << fixed << setprecision(6);
	// Initial guesses and tolerance
	double initial_x0 = 0.0;
	double initial_x1 = 1.0;
	double tolerance = 0.0001;
	int maxIter = 20;

	double root = secantmethod(f, x, initial_x0, initial_x1, tolerance, maxIter);
	cout << "Root of the equation: " << root << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}