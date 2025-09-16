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

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic l1("l1"), l2("l2"), θ1("θ1"), θ2("θ2"), θ̇1("θ1'"),  θ̇2("θ2'"), θ̈1("θ1''"),θ̈2("θ2''");
	Symbolic x1, y1, x2, y2;

	x1 = l1*sin(θ1);
	y1 = -l1*cos(θ1);
	x2 = l1*sin(θ1) + l2*sin(θ2);
	y2 = -l1*cos(θ1) - l2*cos(θ2);

	cout << "x1 = " << x1 << endl;	
	cout << "y1 = " << y1 << endl;	
	cout << "x2 = " << x2 << endl;	
	cout << "y2 = " << y2 << endl;	

	cout << "Here we go..." << lagrangian(x1,y1,x2,y2,θ1, θ̇1, θ̈1, θ2, θ̇2,θ̈2) << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
