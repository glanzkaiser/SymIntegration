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

	Polynomialcoeff P({1.0, 0.0, 0.0});
	Polynomialcoeff Q({0.0, -2.0, 0.0});
	Polynomialcoeff R({0.0, 0.0, 0.0});

	double x0 = 0.0; 
	// Initial conditions: y(0) = y0, y'(0) = dy0 
	double y0 = 1.0; // y(0) = 1
	double dy0 = 0.0; // y'(0) = 0
	int num_terms = 14;


	SecondOrderODE_Homogeneous_PowerSeriesSolver solver(P, Q, R, x0, y0, dy0);
	solver.computeSeries(num_terms);
    
	solver.printSolution();
	cout << endl;

	Polynomialcoeff P2({1.0, 0.0, 0.0});
	Polynomialcoeff Q2({0.0, -2.0, 0.0});
	Polynomialcoeff R2({2.0, 0.0, 0.0});
	y0 = 0.0, dy0 = 1.0;
	SecondOrderODE_Homogeneous_PowerSeriesSolver solver1(P2, Q2, R2, x0, y0, dy0);
	solver1.computeSeries(num_terms);
    
	solver1.printSolution();
	cout << endl;

	Polynomialcoeff P4({1.0, 0.0, 0.0});
	Polynomialcoeff Q4({0.0, -2.0, 0.0});
	Polynomialcoeff R4({4.0, 0.0, 0.0});
	y0 = 1.0, dy0 = 0.0;

	SecondOrderODE_Homogeneous_PowerSeriesSolver solver2(P4, Q4, R4, x0, y0, dy0);
	solver2.computeSeries(num_terms);
	solver2.printSolution();
	cout << endl;

	Polynomialcoeff P6({1.0, 0.0, 0.0});
	Polynomialcoeff Q6({0.0, -2.0, 0.0});
	Polynomialcoeff R6({6.0, 0.0, 0.0});
	y0 = 0.0, dy0 = 1.0;
	SecondOrderODE_Homogeneous_PowerSeriesSolver solver3(P6, Q6, R6, x0, y0, dy0);
	solver3.computeSeries(num_terms);
    
	solver3.printSolution();
	cout << endl;

	Polynomialcoeff P8({1.0, 0.0, 0.0});
	Polynomialcoeff Q8({0.0, -2.0, 0.0});
	Polynomialcoeff R8({8.0, 0.0, 0.0});
	y0 = 1.0, dy0 = 0.0;

	SecondOrderODE_Homogeneous_PowerSeriesSolver solver4(P8, Q8, R8, x0, y0, dy0);
	solver4.computeSeries(num_terms);
	solver4.printSolution();
	cout << endl;

	Polynomialcoeff P10({1.0, 0.0, 0.0});
	Polynomialcoeff Q10({0.0, -2.0, 0.0});
	Polynomialcoeff R10({10.0, 0.0, 0.0});
	y0 = 0.0, dy0 = 1.0;
	SecondOrderODE_Homogeneous_PowerSeriesSolver solver5(P10, Q10, R10, x0, y0, dy0);
	solver5.computeSeries(num_terms);
    
	solver5.printSolution();
	cout << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0; 
}
