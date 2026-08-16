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

	// Example: Solve y'' + y = 0 with y(0) = 1, y'(0) = 0
	vector<double> ode_coeffs = {1.0, 0.0, 1.0}; 
	vector<double> initial_conditions = {1.0, 0.0}; // start from lowest degree
	int num_terms = 11; // Number of series terms to generate

	// Instantiate class object
	HigherOrderODE_Homogeneous_PowerSeriesSolver solver(ode_coeffs, initial_conditions);

	double test_x = 0.45;

	solver.computeSeries(num_terms);
    
	//cout << std::fixed << std::setprecision(10);
	solver.printSeries();

	double approximation = solver.evaluateAt(test_x, 10);
	cout << "\nx = " << test_x << std::endl;
	cout << "\nPower Series Approximation at x : " << std::fixed << std::setprecision(6) << approximation << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0; 
}
