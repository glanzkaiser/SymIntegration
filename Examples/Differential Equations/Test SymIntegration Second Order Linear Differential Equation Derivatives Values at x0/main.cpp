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

	Symbolic x("x"), y("y");
	Symbolic f1 = df(y[x],x,2) +x*df(y[x],x) + y[x];

	Symbolic y0 = 1, dy0 = 0;
	secondorderlineardiffeq_derivativesvalueatx0(f1,y, x,0,y0,dy0);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}