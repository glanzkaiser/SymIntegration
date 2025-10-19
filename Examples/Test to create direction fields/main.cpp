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
#include <vector>
#define π 3.1415926535897f
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// Define your differential equation
double f(double x, double y) {
    return 5 - 3*sqrt(y); // Example: dy/dx = 5 - 3*sqrt(y)
}

int main() {
	// Open a file to write data for Gnuplot
	ofstream dataFile("direction_field.dat");
	ofstream dataFile2("direction_field2.dat");
	
	static double xs[500], ys[500];
	int i = 1;
	double x_min = 0.0, x_max = 10.0;
	double y_min = 0.0, y_max = 10.0;
	double step_size = 0.5;
	double dx, dy, magnitude;
	for (double y = y_min; y <= y_max; y += step_size) {
	for (double x = x_min; x <= x_max; x += step_size) {
		if (std::abs(y) > 0.001) { // Avoid division by zero for this example
		double slope = f(x, y);
		// Calculate components of a unit vector in the direction of the slope
		dx = 1.0;
		dy = slope;
		magnitude = sqrt(dx * dx + dy * dy);
		dx = dx / magnitude;
		dy = dy / magnitude;
		
                // Write starting point (x,y) and vector components (dx, dy) to file
                dataFile << x << " \t \t  " << y << " \t  \t " << dx << " \t  \t " << dy << endl;
            }
		xs[i] = x;
		ys[i] = y;
		cout << "i = " << i << endl;
		cout << "xs[i] ys[i] = " << xs[i] << " , " << ys[i] << endl;
		dataFile2 << i << "\t \t \t "<< xs[i] << " \t \t \t  " << ys[i]  << " \t  \t \t "<< dy/dx << endl;
		i = i+1;
		xs[i] = x+dx;
		ys[i] = y+dy;
		cout << "xs[i] ys[i] = " << xs[i] << " , "  << ys[i] << endl;
		dataFile2 << i << "\t \t \t " << xs[i] << " \t \t \t  " << ys[i]  << " \t  \t \t "<< dy/dx << endl;

		i = i+1;
        }
		
	}
    dataFile.close();

    return 0;
}