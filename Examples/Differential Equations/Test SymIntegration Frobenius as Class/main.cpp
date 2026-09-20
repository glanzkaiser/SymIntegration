/*
    SymIntegration is branching from SymbolicC++ 3.35
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
// THANKS SENTINEL!!! and Freya too

// g++ -o result main.cpp -lsymintegration

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

#define π 3.1415926535897f


using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;


int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Real, distinct roots, roots differ by an integer
	// Roots: r1 = 1, r2 = -1
 	PolynomialDouble P({0.0, 0.0, 1.0});   // P(x) = x^2
	PolynomialDouble Q({0.0, 1.0});        // Q(x) = x
	PolynomialDouble R({-1.0, 0.0, 1.0}); // R(x) = x^2-1

	// Example ODE: (x - 1)^2 y'' + (x - 1) y' - 0.25 y = 0
	// Shifted regular singular point: x0 = 1.0
	// Polynomial inputs in standard powers of x: x^2 - 2x + 1, x - 1, and -0.25
	PolynomialDouble P0({1.0, -2.0, 1.0});   
	PolynomialDouble Q0({-1.0, 1.0});       
	PolynomialDouble R0({-0.25}); 

	// Repeated roots
	PolynomialDouble P1({0.0, 0.0, 1.0});   // P(x) = x^2
	PolynomialDouble Q1({0.0, 5.0});        // Q(x) = 5x
	PolynomialDouble R1({4.0, 0.0, 0.0}); // R(x) =  4

	// Complex roots, set x0 = 2
	// Roots: r1 =  -2 + 2i, r2 = -2 - 2i. 
	PolynomialDouble P2({4.0, -4.0, 1.0});   // P(x) = x^2 - 4x + 4 
	PolynomialDouble Q2({-10.0, 5.0});        // Q(x) = 5x - 10
	PolynomialDouble R2({8.0, 0.0, 0.0}); // R(x) = 8

	// Example: Bessel's Equation of order 1 -> 2x^2*y'' + 3x*y' + (x^2 - 1)*y = 0
	// Real, distinct roots, roots differ by non-integer
	// Roots: r1 = 0.5, r2 = -1. Difference = 1.5.
	PolynomialDouble P3({0.0, 0.0, 2.0}); // P(x) = 2x^2
	PolynomialDouble Q3({0.0, 3.0});      // Q(x) = 3x
	PolynomialDouble R3({-1.0, 0.0, 1.0}); // R(x) = -1 + x^2

	// Real, distinct roots, roots differ by integer 
	// Example: Bessel-like equation scaled with a shifted regular singular point at x0 = 1.0
	// Equation of the form: (x-1)^2 y'' + (x-1) y' + ((x-1)^2 - 0.25) y = 0
	// Expanded around x0 = 1, substituting t = x - 1:
	// P(x) = x^2 - 2x + 1
	// Q(x) = x - 1
	// R(x) = x^2 - 2x + 0.75
	// Roots: r1 = 0.5, r2 = -0.5. 
	PolynomialDouble P4({1.0, -2.0, 1.0}); 
	PolynomialDouble  Q4({-1.0, 1.0});      
	PolynomialDouble  R4({0.75, -2.0, 1.0}); 

	// Real, distinct roots, roots differ by integer
	// Roots: r1 = 1, r2 = -1. 
	PolynomialDouble P5({0.0, 0.0, 1.0}); 
	PolynomialDouble Q5({0.0, 1.0});     
	PolynomialDouble R5({-1.0, 0.0, 1.0}); 

	// Repeated roots
	PolynomialDouble P6({0.0, 0.0, 1.0}); 
	PolynomialDouble Q6({0.0, 1.0, 0.0});     
	PolynomialDouble R6({0.0, 0.0, 1.0}); 

	// Repeated roots
	PolynomialDouble P7({0.0, 0.0, 1.0}); 
	PolynomialDouble Q7({0.0, 3.0, 0.0});     
	PolynomialDouble R7({1.0, 1.0, 0.0}); 

	// Complex roots r1 = 2i, r2 = -2i
	PolynomialDouble P8({0.0, 0.0, 1.0}); 
	PolynomialDouble Q8({0.0, 1.0, 0.0});     
	PolynomialDouble R8({4.0, 0.0, 1.0}); 

	// Complex roots 
	PolynomialDouble P9({0.0, 0.0, 1.0}); 
	PolynomialDouble Q9({0.0, 2.5, 0.0});     
	PolynomialDouble R9({1.0, 1.0, 0.0}); 

	complex<double> x0(0.0,0.0); // Centered around regular singular point x0 = 0.
	int terms = 16;

	// Provide initial values close to regular singular boundary (e.g., at x = 0.1)
	double x_init = 1.1; 	// t0
	double y_init = 0.84147;  // Approximating sin(1)/sqrt(1) 		y(t0)
	double dy_init = -0.13938;				// y'(t0)
	vector<complex<double>> evaluation_points = {{1.0,0.0}, {1.2,0.0}, {1.4,0.0}, {1.6,0.0}, {1.8,0.0}, {2.0,0.0}};
	complex<double> c0_first = 1.0;
	SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver solver1(P9, Q9, R9, x0, terms);

	solver1.computeSeriesCoefficients();
	solver1.printSolution();
	solver1.printCoefficients();
	solver1.solve_ivp(x_init, y_init, dy_init, evaluation_points);
	
	
	
	/*auto [r1, r2] = solver1.solve_indicial_equation_inpair();

	cout << "\nIndicial Roots Found:\n";
	cout << "r1 = " << r1 << "\n";
	cout << "r2 = " << r2 << "\n\n";
	*/

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0; 
}
