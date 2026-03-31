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
#include <chrono>

#define π 3.1415926535897f

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

	Symbolic y("y"),t("t");

	cout << "\nSolution for u'' + u' + 1.25u = 3 cos(t), u(0) = 2, u'(0) = 3 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationsivpsolution(1,1,1.25,3*cos(t),2,3,y,t);		

	cout << "\nSolution for u'' + 0.125u' + u = 3 cos(t), u(0) = 2, u'(0) = 0 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationsivpforcedvibrationssolution(1,0.125,1,3*cos(t),2,0,y,t);

	cout << "\nSolution for u'' + u = 0.5 cos(0.8t), u(0) = 0, u'(0) = 0 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationsivpforcedvibrationssolution(1,0,1,0.5*cos(0.8*t),0,0,y,t);

	cout << "\nSolution for u'' + 0.25u' + 2u = 2 cos(1.403t), u(0) = 0, u'(0) = 2 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationsivpforcedvibrationssolution(1,0.25,2,2*cos(1.403*t),0,2,y,t);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	

	return 0; 
}
