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

	// --- INPUT CONFIGURATIONS AND INITIALIZATION DATA INSIDE MAIN ---
	int total_cycles = 5000;           // Global replica configuration swap cycles
	int local_steps = 10;              // Exploration steps per replica inside each cycle
    
	// Define replica temperature ladders (coldest chain to hot exploration chains)
	vector<double> pt_temps = {0.05, 0.22, 0.93, 4.0}; 
    
	double s = 0.75;
	double t = 1.0;
    
	const int DIMENSIONS = 5;
    
	// Provide a starting coordinate vector point (anchor node) for replica space mapping
	vector<double> initial_base(DIMENSIONS, 1.0); 

	cout << "=====================================================================" << endl;
	cout << "         N-DIMENSIONAL PARALLEL TEMPERING OPTIMIZATION ENGINE        " << endl;
	cout << "=====================================================================" << endl;
	cout << "Target Space Resolution Size: " << DIMENSIONS << " Dimensions" << endl;
	cout << "Active Replica Chain Count:     " << pt_temps.size() << " system nodes" << endl;
	cout << "Running Execution Loop for " << total_cycles << " global swap configurations..." << endl;

	// Wrap the objective function in a lambda matching std::function signature
	auto bound_target_func = [s, t](const vector<double>& v) 
	{
		return objective_function(v, s, t);
	};

	// Call execution function engine passing the lambda functional hook
	PTMetrics results = run_parallel_tempering(
	total_cycles,
	local_steps,
	initial_base,
	pt_temps,
	bound_target_func
	);

	// --- Performance Summaries Report ---
	cout << "\nOptimization Search Completed." << endl;
	cout << "Total Swap Attempts Checked:   \t" << results.swap_attempts << endl;
	cout << "Successful Swap Transactions:  \t" << results.swap_successes << endl;
	if (results.swap_attempts > 0) 
	{
		cout << "Swap Acceptance Ratio:          \t" << fixed << setprecision(2) 
		<< (static_cast<double>(results.swap_successes) / results.swap_attempts) * 100.0 << "%" << endl;
	}

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
