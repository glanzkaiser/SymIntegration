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



// ==========================================
// 5. MAIN SYSTEM ENTRY
// ==========================================

// Your system definition parsed natively using our DualClassTemplate templates
vector<DualClassTemplate<complex<double>>> targetSystem(const vector<DualClassTemplate<complex<double>>>& x) 
{
	auto x1 = x[0];
	auto x2 = x[1];
	auto x3 = x[2];

	vector<DualClassTemplate<complex<double>>> eq(3);
	// eq1 = 3*x1 - cos(x1*x2) - 0.5 = 0
	eq[0] = x1 * 3.0 - cos(x1 * x2) - 0.5;
	// eq2 = x1^2 - 81*(x2+0.1)^2 + sin(x3) + 1.06 = 0
	eq[1] = x1 * x1 - (x2 + 0.1) * (x2 + 0.1) * 81.0 + sin(x3) + 1.06;
	// eq3 = e^(-x1*x2) + 20*x3 + (10*pi-3)/3 = 0
	eq[2] = exp(x1 * x2 * -1.0) + x3 * 20.0 + ((10.0 * PI - 3.0) / 3.0);

	return eq;
}

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Set of diverse initial complex guesses to capture different roots
	vector<vector<complex<double>>> initial_guesses = {
	{ {0.1, 0.0},  {0.1, 0.0},  {-0.1, 0.0} },  // Close to primary real root
	{ {1.0, 1.0},  {-1.0, 1.0}, {0.0, -1.0} },  // Deep Complex Region A
	{ {-0.5, -0.5}, {0.5, 0.5},  {1.0, 1.0} }   // Deep Complex Region B
	};

	RK4MultiRootSolver solver1(initial_guesses, targetSystem, 150);
	solver1.printSolution();
	

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
