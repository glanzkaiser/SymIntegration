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
using namespace SymbolicConstant;


int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Example System: 
	//  10x1 - x2 + 2x3            = 6
	//  -x1 + 11x2 - x3 + 3x4  = 25
	//   2x1 - x2 + 10x3 - x4    = -11
	//            3x2 - x3 + 8x4    = 15
	vector<vector<double>> A = loadMatrixFromFile("matrixA.txt");
	vector<double> b = loadVectorFromFile("vectorb.txt");

	/*vector<vector<double>> A = {
	{10.0, -1.0, 2.0, 0.0},
	{-1.0, 11.0, -1.0, 3.0},
	{2.0, -1.0, 10.0, -1.0},
	{0.0, 3.0, -1.0, 8.0}
	};
	vector<double> b = {6.0, 25.0, -11.0, 15.0};*/

	// 1. Convert Dense to CRS representation
	CRSMatrix A_crs = denseToCRS(A);
    
	int maxIterations = 20;
	double tolerance = 0.001;

	// 2. Solve Ax = b using Conjugate Gradient
	vector<double> x = ConjugateGradient(A_crs, b, tolerance, maxIterations);
    
	// Print results
	cout << "\nSolution vector x:\n";
	cout << std::fixed << std::setprecision(8);
	for (size_t i = 0; i < x.size(); ++i) 
	{
		cout << "x[" << i << "] = " << x[i] << "\n";
	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
