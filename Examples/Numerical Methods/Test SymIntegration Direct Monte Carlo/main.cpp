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

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;

// --- Generalized N-Dimensional Objective Function ---
double objective_function(const vector<double>& v, double s, double t) 
{
	double r2 = 0.0;
	// Calculate square of the radius across all N dimensions dynamically
	for (double coord : v) 
	{
		r2 += coord * coord;
	}
    
	double term1 = exp(-r2 / (2.0 * s * s));
	double term2 = exp(-r2 / (2.0 * t * t));
	double penalty = 0.1 * r2;
    
	return term1 - term2 + penalty;
}

// ============================================================================
// MAIN SYSTEM ORCHESTRATOR
// ============================================================================
int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// --- ONLY INPUT CONFIGURATIONS AND INITIALIZATION DATA INSIDE MAIN ---
	int total_samples = 500000; 
	double s = 0.75;
	double t = 1.0;
    
	const int DIMENSIONS = 5;
    
	vector<double> lower_bounds(DIMENSIONS, -2.0);
	vector<double> upper_bounds(DIMENSIONS,  2.0);

	cout << "=====================================================================" << endl;
	cout << "       N-DIMENSIONAL DIRECT MONTE CARLO OPTIMIZATION ENGINE          " << endl;
	cout << "=====================================================================" << endl;
	cout << "Target Space Resolution Size: " << DIMENSIONS << " Dimensions" << endl;
	cout << "Search Envelope Matrix Bounds:  [-2.0, 2.0] across all dimensions" << endl;
	cout << "Running Execution Loop for " << total_samples << " random trial allocations..." << endl;

	// Fix: Wrap the objective function in a lambda matching std::function signature
	auto bound_target_func = [s, t](const vector<double>& v) 
	{
		return objective_function(v, s, t);
	};

	// Call execution function engine passing the lambda functional hook
	DMCMetrics results = run_direct_monte_carlo(
	total_samples,
	lower_bounds,
	upper_bounds,
	bound_target_func
	);

	// --- Performance Summaries Report ---
	cout << "\nOptimization Search Completed." << endl;
	cout << "Total Random Samples Evaluated: \t" << results.total_samples << endl;

	cout << "\n--- Absolute Best Global Minimum Discovered ---" << endl;
	cout << "Minimum Energy Value Found: \t\t" << setprecision(12) << results.best_energy << endl;
	cout << "Best Coordinate Vector Array: \t[ ";
	for (double val : results.best_x) 
	{
 		cout << fixed << setprecision(6) << val << " ";
	}
	cout << "]" << endl;
	cout << "=====================================================================" << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}
