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

// g++ -o result main.cpp -lsymintegration -lboost_iostreams
// make
// ./main > lattice_2d.dat


#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>

#include "gnuplot-iostream.h"

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Define the dimensions of your grid
	//int numRows = 9;
	//int numCols = 3;

	// Create a 2D vector (vector of vectors) to represent the grid
	// Initialize all cells with a default value (e.g., 0 for integers)
	//vector<vector<int>> grid(numRows, vector<int>(numCols, 1));

	// Access and modify grid elements
	//grid[0][0] = 10; // Set the value of the cell at row 0, column 0
	//grid[2][3] = 25; // Set the value of the cell at row 2, column 3

	int n = 4;
	// Print the grid to the console	
	for (int k = -3; k < n; ++k)
	{
		for (int j = -3; j < n; ++j) 
		{
			cout << k << " " << j << " " << 1 ;
			cout << endl;
		}	
	}

	Gnuplot gp;

	// Don't forget to put "\n" at the end of each line!
	gp << "set xrange [-3:3]\n";
	gp << "set yrange [-3:3]\n";
	gp << "set xtics 1\n";
	gp << "set title '2D Lattice' \n";
	gp << "set key off\n"; // supress the key / no legend
	gp << "plot 'lattice_2d.dat' with points pt 7 lc rgb 'blue'\n";

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	//cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
