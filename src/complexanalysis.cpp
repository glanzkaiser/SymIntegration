/*
   Thanks Freya the Goddess, RK, Sentinel, Berlin, Mother Mary.
*/
// Search for  // NOT YET to see which function we haven't finish yet

#include "symintegral/symintegrationc++.h"

using namespace std;

// Type alias for convenience
using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_COMPLEXANALYSIS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_COMPLEXANALYSIS_DEFINE

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f
#define pi 3.14159265358979323846

#include "symintegrationc++.h"
#include <string>
#include<iostream>
#include<vector>
#include <fstream>
#include <cmath>
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()
#include <map>
#include <algorithm> // For std::max_element,  std::sort, std::reverse ,  std::for_each
#include <numeric> // For std::accumulate
#include <random> // For random number generation
#include <complex>

using namespace SymbolicConstant;

void complexnumber_analysis(complex<double> &c)
{
	double a = real(c);
	double b = imag(c);

	double modulus_z =sqrt(a*a + b*b);
	double angle_radians = std::arg(c);
	double angle_degrees = divisiond(angle_radians * 180.0 , pi);
	
	cout << "\nThe complex number c is: " << c << endl;
	cout << "\nModulus (|z|) = " << modulus_z << endl;
	cout << "Polar angle (radians): " << angle_radians << endl;
	cout << "Polar angle (degrees): " << angle_degrees << endl;
	
	Symbolic polar_form = modulus_z*(cos(angle_radians)  + SymbolicConstant::i*sin(angle_radians));
	cout << "Polar form : \nz = |z| (cos Φ + i sin Φ) = " << polar_form << endl;
	
}
#endif
#endif