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

	Symbolic x("x"), z("z");
	
	Symbolic fx = pow((x-1),Symbolic(divisiond(-2,3)));
	cout << "\nf(x) = "<< fx <<endl;

	Symbolic y = integrate(fx,x) ;
	cout << "integral of f(x) = " << y << endl;
	cout << "integral of f(x) from 0 to 1^{-} = " << y[x==1] - y[x==0] << endl;
	cout << "integral of f(x) from 1^{+} to 3 = " << y[x==3] - y[x==1] << endl;
	
	cout << "integral of f(x) from 0 to 3 = " << (y[x==1] - y[x==0] ) + (y[x==3] - y[x==1]) << endl;

	cout << "\npow(-1,4) = " << pow(-1,4) << endl;
	cout << "pow(-1,4.1) = " << pow(x,z)[x==-1, z==4.1] << endl;
	cout << "pow(-1,0.5) = " << pow(x,z)[x==-1,z==0.5] << endl;

	cout << "\nfloor(0.5) = " << floor(0.5) << endl;
	cout << "floor(2.1) = " << floor(2.1) << endl;
	cout << "floor(-2.1) = " << floor(-2.1) << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}