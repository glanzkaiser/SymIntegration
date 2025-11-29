// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include "symintegrationc++.h"

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

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

	cout << "\nbinomialpmf(2;4,0.75) = " << binomialpmf(2,4,0.75) << endl;
	cout << "\nbinomialcdf(8;15,0.4) = " << binomialcdf(8,15,0.4) << endl;
	cout << "\nbinomialmean(8;15,0.4) = " << binomialmean(8,15,0.4) << endl;
	cout << "\nbinomialvar(8;15,0.4) = " << binomialvar(8,15,0.4) << endl;
	cout << "\nbinomialmgf(8;15,0.4) = " << binomialmgf(8,15,0.4) << endl;
	
	cout << "\nnegativebinomialpmf(6;4,0.55) = " << negativebinomialpmf(6,4,0.55) << endl;
	cout << "\nnegativebinomialmean(6;4,0.55) = " << negativebinomialmean(6,4,0.55) << endl;
	cout << "\nnegativebinomialvar(6;4,0.55) = " << negativebinomialvar(6,4,0.55) << endl;
	cout << "\nnegativebinomialmgf(6;4,0.55) = " << negativebinomialmgf(6,4,0.55) << endl;
	
	cout << "\ngeometricpmf(5;0.01) = " << geometricpmf(5,0.01) << endl;
	cout << "\ngeometricmean(5;0.01) = " << geometricmean(5,0.01) << endl;
	cout << "\ngeometricvar(5;0.01) = " << geometricvar(5,0.01) << endl;
	cout << "\ngeometricmgf(5;0.01) = " << geometricmgf(5,0.01) << endl;

	cout << "\nhypergeometricgeometricpmf(1;40,5,3) = " << hypergeometricpmf(1,40,5,3) << endl;
	cout << "\nhypergeometricgeometricmean(1;40,5,3) = " << hypergeometricmean(1,40,5,3) << endl;
	cout << "\nhypergeometricgeometricvar(1;40,5,3) = " << hypergeometricvar(1,40,5,3) << endl;
	
	cout << "\npoissonpmf(6;4) = " << poissonpmf(6,4) << endl;
	cout << "\npoissoncdf(6;4) = " << poissoncdf(6,4) << endl;
	cout << "\npoissonmean(6;4) = " << poissonmean(6,4) << endl;
	cout << "\npoissonvar(6;4) = " << poissonvar(6,4) << endl;
	cout << "\npoissonmgf(6;4) = " << poissonmgf(6,4) << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}