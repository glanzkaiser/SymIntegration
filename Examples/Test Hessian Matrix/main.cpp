// g++ -o result main.cpp -lsymintegration -larmadillo 
// Merci beaucoup Freya..
// C++ program to plot the particular and general solution of linear system Ax=b and Ax=0

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector> // For std::vector (example container)
#include "symintegrationc++.h"
#include <algorithm> // For std::sort

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

#include <chrono>
#include <string>
#include <bitset>

using namespace std::chrono;
using namespace std;

double division(double x, double y)
{
	return x/y;
}
Symbolic divisionsym(Symbolic x, Symbolic y)
{
	return x/y;
}
// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic x("x"), y("y");
	
	Symbolic f = 2*sin(3*x) +  5 + x*x*y + cos(y);
	//Symbolic H = (( df(df(f,x),x), df(df(f,x),y) ), (df(df(f,y),x), df(df(f,y),y)) );
	cout << "\nf(x,y) = " << f << endl;
	cout << "\nH = " <<  hessian(f,x,y) << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}