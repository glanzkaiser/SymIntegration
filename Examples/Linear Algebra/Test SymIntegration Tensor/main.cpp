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
#include "vector.h"
#include "vecnorm.h"

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// 1. Declare the 3D Vector (vector of matrices)
	vector<vector<std::vector<double>>> threeDVector;

	// 2. Define a 2D Matrix (vector of vectors)
	std::vector<std::vector<double>> matrix = {
	{1.0, 2.0, 3.0},                
	{4.0, 5.0, 6.0}
	};

	// 3. Push back the 2D Matrix
	threeDVector.push_back(matrix);

	// Alternative: Push back directly (creates temporary)
	threeDVector.push_back({{7.0, 8.0}, {9.0, 10.0}});

	std::cout << "3D Vector Size: " << threeDVector.size() << std::endl;
	printMatrix(threeDVector[0]);
	printMatrix(threeDVector[1]);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
