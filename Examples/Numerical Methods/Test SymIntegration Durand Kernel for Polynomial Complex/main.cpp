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

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Example: P(x) = x^3 - 3x^2 + 3x - 1 (Roots should be exactly 1, 1, 1)
	// Coeffs: [-1, 3, -3, 1]
	PolynomialComplex p1({-1.0, 3.0, -3.0, 1.0});
	vector<complex<double>> roots1 = p1.find_roots_durand_kerner();

	cout << "Roots for x^3 - 3x^2 + 3x - 1:\n";
	for (const auto& root : roots1) 
	{
		cout << root << "\n";
	}

	PolynomialComplex p2({1.0, 0.0, 1.0});
	vector<complex<double>> roots2 = p2.find_roots_durand_kerner();

	cout << "\nRoots for x^2 + 1:\n";
	for (const auto& root : roots2) 
	{
		cout << root << "\n";
	}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
