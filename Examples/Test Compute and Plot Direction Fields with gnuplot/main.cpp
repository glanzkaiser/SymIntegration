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

#include "symintegrationc++.h"
#include <vector>
#define π 3.1415926535897f
#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>

using namespace std::chrono;
using namespace std;

// Define your differential equation
double f(double x, double y) {
    return y*(3 - x*y); // Example: dy/dx = 5 - 3*sqrt(y)
}

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Open a file to write data for Gnuplot
	ofstream dataFile("direction_field.dat");
	
	double x_min = 0.0, x_max = 10.0;
	double y_min = -5.0, y_max = 5.0;
	double step_size = 0.3;
	double dx, dy, magnitude;
	for (double y = y_min; y <= y_max; y += step_size) 
	{
		for (double x = x_min; x <= x_max; x += step_size) 
			{
				if (std::abs(y) > 0.001) 
				{ // Avoid division by zero for this example
				double slope = f(x, y);
				// Calculate components of a unit vector in the direction of the slope
				dx = 1.0;
				dy = slope;
				magnitude = sqrt(dx * dx + dy * dy);
				dx = 0.2 * dx / magnitude;
				dy = 0.2 * dy / magnitude;
				
				// Write starting point (x,y) and vector components (dx, dy) to file
				dataFile << x << " \t \t  " << y << " \t  \t " << dx << " \t  \t " << dy << endl;
				}
			}	
	}

	dataFile.close();

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	
    return 0;
}