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

using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	Symbolic x = 2627.2, y = 100.9, z = 2961.3;
	Symbolic r = 4, t = 2*π/3, zc = 5;
	Symbolic ρ = 3960;
	Symbolic Φ = 0.7261; // elevation
	Symbolic θ = 0.0384; // azimuth

	cout << "Spherical to cartesian (ρ = 3960, θ = 0.0384, Φ=0.7261) = " << sphericaltocartesian(ρ, θ, Φ) << endl;
	cout << "Cartesian to spherical (x = 2627.2, y=100.9 , z = 2961.3) = " << cartesiantospherical(x, y, z) << endl;
	
	cout << "Cylindrical to cartesian (r = 4, θ = 2*π/3, z = 5) = " << cylindricaltocartesian(r, t, zc) << endl;
	cout << "Cartesian to cylindrical (x = -5, θ = -5, z = 2) = " << cartesiantocylindrical(-5, -5, 2) << endl;
	
	//cout << cartesiantocylindrical(-5, -5, 2)(0) << endl;
	return 0; 
}
