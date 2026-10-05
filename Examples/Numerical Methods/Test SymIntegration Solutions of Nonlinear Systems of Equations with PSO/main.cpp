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
#include <functional>

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

using ComplexVector = vector<Complex>;
using ComplexVectorFunction = function<ComplexVector(const ComplexVector&)>;

// ============================================================================
// 4. MAIN CONTROLLER
// ============================================================================

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	 // Define the nonlinear system inside main() as a lambda
	ComplexVectorFunction nonlinear_system = [](const ComplexVector& x) -> ComplexVector 
	{
		Complex x1 = x[0];
		Complex x2 = x[1];
		Complex x3 = x[2];

		// Eq 1: 3*x1 - cos(x1*x2) - 0.5 = 0
		Complex eq1 = Complex(3)*x1 - Complex(cos(x2*x3)) - Complex(0.5);

		// Eq 2: x1^2 - 81*(x2 + 0.1)^2 + sin(x3) + 1.06 = 0
		Complex eq2 = x1*x1 -Complex(81)*(x2+ 0.1)*(x2 + Complex(0.1)) + Complex(sin(x3)) + Complex(1.06);

		// Eq 3: e^(-x1*x2) + 20*x3 + (10*pi - 3)/3 = 0
		Complex eq3 = Complex(exp(-x1*x2)) + ( Complex(20) * x3 ) + Complex(divisiond(10*PI-3,3));

		return {eq1, eq2, eq3};
	};

	// Initialize the Multi-Root Finder (3 complex variables, 100 particles, 1000 iterations, lb -4, ub 4)
	PSOMultiRootFinder finder(3, 100, 1000, -4.0, 4.0);

	int total_runs = 30; // Perform multiple search routines to isolate distinct roots	
	finder.findMultipleRoots_printSolution(nonlinear_system, total_runs);

	//vector<ComplexVector> roots = finder.findMultipleRoots(nonlinear_system, total_runs);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
