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

	Symbolic x("x");
	
	cout << "\nTaylor series for sin(x) , x_{0} = 0\n" <<endl;
	cout << taylorseries(sin(x),x, 0, 10) << endl;		

	cout << "\nTaylor series for exp(x) , x_{0} = 0\n" <<endl;
	cout << taylorseries(exp(x),x, 0, 10) << endl;		

	cout << "\nTaylor series for x, x_{0} = 1 \n" <<endl;
	cout << taylorseries(x,x, 1, 10) << endl;		

	cout << "\nTaylor series for x^{2}, x_{0} = -1 \n" <<endl;
	cout << taylorseries(pow(x,Symbolic(2)),x, -1, 10) << endl;		

	cout << "\nTaylor series for ln(x), x_{0} = 1 \n" <<endl;
	cout << taylorseries(ln(x),x, 1, 10) << endl;		

	cout << "\nTaylor series for 1/(1-x), x_{0} = 0 \n" <<endl;
	cout << taylorseries(1/(1-x),x, 0, 10) << endl;		

	cout << "\nTaylor series for 1/(1+x), x_{0} = 0 \n" <<endl;
	cout << taylorseries(1/(1+x),x, 0, 10) << endl;		

	cout << "\nTaylor series for 1/(1-x), x_{0} = 2 \n" <<endl;
	cout << taylorseries(1/(1-x),x, 2, 10) << endl;		

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}