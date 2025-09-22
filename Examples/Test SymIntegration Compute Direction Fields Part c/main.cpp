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


int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Open a file to write data for Gnuplot
	Symbolic ys("ys"), ts("ts");
	Symbolic f = -ts*ys + 0.1*(pow(ys,Symbolic(3)));

	directionfield(f,ts,ys,0,17,-8,8,0.4,0.3);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	
    return 0;
}