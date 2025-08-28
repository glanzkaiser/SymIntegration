// g++ -o result main.cpp -lsymintegration 
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include "symintegrationc++.h"

#include <chrono>

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

	cout << "\nuniformpdf(1;1,5) = " << uniformpdf(1,1,5) << endl;
	cout << "\nuniformcdf(3;1,5) = " << uniformcdf(3,1,5) << endl;
	cout << "\nuniformmgf(x;1,5) = " << uniformmgf(1,1,5) << endl;
	cout << "\nuniformmean(x;1,5) = " << uniformmean(1,1,5) << endl;
	cout << "\nuniformvar(x;1,5) = " << uniformvar(1,1,5) << endl;
	
	cout << "\nnormalpdf(362;300,50) = " << normalpdf(362,300,50) << endl;
	cout << "\nnormalcdf(362;300,50) = " << normalcdf(362,300,50) << endl;
	cout << "\nnormalmgf(362;300,50) = " << normalmgf(362,300,50) << endl;
	cout << "\nnormalmean(362;300,50) = " << normalmean(362,300,50) << endl;
	cout << "\nnormalvar(362;300,50) = " << normalvar(362,300,50) << endl;

	cout<<"\ngammapdf(1,2,0.2) = " << gammapdf(1,2,0.2) << endl;
	cout<<"\ngammacdf(1,2,0.2) = " << gammacdf(1,2,0.2) << endl;
	cout<<"\ngammamgf(1,2,0.2) = " << gammamgf(1,2,0.2) << endl;
	cout<<"\ngammamean(1,2,0.2) = " << gammamean(1,2,0.2) << endl;
	cout<<"\ngammavar(1,2,0.2) = " << gammavar(1,2,0.2) << endl;

	cout<<"\nexponentialpdf(8,0.2) = " << exponentialpdf(8,0.2) << endl;
	cout<<"\nexponentialcdf(8,0.2) = " << exponentialcdf(8,0.2) << endl;
	cout<<"\nexponentialmgf(8,0.2) = " << exponentialmgf(8,0.2) << endl;
	cout<<"\nexponentialmean(8,0.2) = " << exponentialmean(8,0.2) << endl;
	cout<<"\nexponentialvar(8,0.2) = " << exponentialvar(8,0.2) << endl;

	cout<<"\nbetapdf(0.8,3,2) = " << betapdf(0.8,3,2) << endl;
	cout<<"\nbetacdf(0.8,3,2) = " << betacdf(0.8,3,2) << endl;
	cout<<"\nbetamgf(0.8,3,2) = " << betamgf(0.8,3,2) << endl;
	cout<<"\nbetamean(0.8,3,2) = " << betamean(0.8,3,2) << endl;
	cout<<"\nbetavar(0.8,3,2) = " << betavar(0.8,3,2) << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}