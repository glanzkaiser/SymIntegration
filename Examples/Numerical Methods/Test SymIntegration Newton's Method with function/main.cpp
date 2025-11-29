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

using namespace std::chrono;
using namespace std;


int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	Symbolic x("x"), f, p0;
	p0 = (pi/4);
	int N = 5;
	
	f = cos(x) - x; // the input x has to be in Radian
	cout << "Newton method : " << newtonmethod(f,x,(pi/4),5) << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}