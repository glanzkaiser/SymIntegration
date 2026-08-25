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

	Polynomialcoeff P({1.0, 0.0, -1.0});
	Polynomialcoeff Q({0.0, -1.0, 0.0});

	double x0 = 0.0; 
	// Initial conditions: y(0) = y0, y'(0) = dy0 -> c0 = ... , c1 = ...
	double y0 = 0.0; // y(0)
	// we set this to 1 to compute the solution of the second kind 
	double dy0 = 1.0; // y'(0) 
	int num_terms = 14;

	double test_x = 0.45;

	for (int i = 1; i<10;)
	{
		Polynomialcoeff R({double(i), 0.0, 0.0});	
		SecondOrderODE_Homogeneous_PowerSeriesSolver solver(P, Q, R, x0, y0, dy0);
		solver.computeSeries(num_terms);
    
		solver.printSolution();
		i=i+2;
	}
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0; 
}
