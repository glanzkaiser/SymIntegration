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

// ============================================================================
// 4. MAIN METHOD DEFINE THE NONLINEAR SYSTEM
// ============================================================================
vector<DualClass> targetNonlinearSystem(const vector<DualClass>& x) {
    DualClass x1 = x[0];
    DualClass x2 = x[1];
    DualClass x3 = x[2];

    vector<DualClass> eq(3);
    // eq1 = 3*x1 - cos(x1*x2) - 0.5 = 0
    eq[0] = x1 * 3.0 - cos(x1 * x2) - 0.5;

    // eq2 = x1^2 - 81*(x2 + 0.1)^2 + sin(x3) + 1.06 = 0
    eq[1] = pow(x1, 2) - 81.0 * pow(x2 + 0.1, 2) + sin(x3) + 1.06;

    // eq3 = e^(-x1*x2) + 20*x3 + (10*pi - 3)/3 = 0
    eq[2] = exp(-1.0 * x1 * x2) + x3 * 20.0 + (10.0 * PI - 3.0) / 3.0;

    return eq;
}

// ==========================================
// 5. Main Execution Block
// ==========================================
int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Initial Guess vector (Complex Domain initialized around realistic roots)
	vector<complex<double>> initial_guess = {
	complex<double>(0.1, 0.0),
	complex<double>(0.1, 0.0),
	complex<double>(-0.1, 0.0)
	};

	cout << "Initial Guess coordinates:\n";
	for(size_t i = 0; i < initial_guess.size(); ++i) 
	{
		cout << "x" << i+1 << " = " << initial_guess[i] << "\n";
	}
	cout << endl;

	// Solve the system, create the object first from the class
	QuasiNewtonSolver solver1(100, 1e-9);
	// Execute root solver
	solver1.solve(targetNonlinearSystem, initial_guess);
	solver1.printSolution();

	/*
	vector<complex<double>> roots = solver1.solve_getvector(targetNonlinearSystem, initial_guess);
	cout << "\nCalculated Complex Roots:\n";
	for (size_t i = 0; i < roots.size(); ++i) 
	{
		cout << "x" << i + 1 << " = " << roots[i] << "\n";
	}

	// Final Verification Check
	vector<DualClass> final_eval_dual(roots.size());

	for (size_t i = 0; i < roots.size(); ++i) 
	{
		final_eval_dual[i] = DualClass(roots[i]);
	}

	vector<DualClass> residuals = targetNonlinearSystem(final_eval_dual);
	cout << "\nSystem Residuals (Should equal close to 0):\n";
	for (size_t i = 0; i < residuals.size(); ++i) 
	{
		cout << "eq" << i + 1 << " residual = " << residuals[i].val << "\n";	
	}*/

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
