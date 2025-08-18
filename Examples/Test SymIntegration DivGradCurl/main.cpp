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
	Symbolic x("x"), y("y"), z("z"), f, f2, i("i"), j("j"), k("k"), î("î"), ĵ("ĵ"), k̂("k̂");
	
	f = -(0.5*x)*y+2*y + z*z;
	f2 = y*i -x*j +y*z*k;
	Symbolic f2_hat = y*î -x*ĵ +y*z*k̂;
	cout << "F(x,y,z) = " << f << endl;
	cout << "F2 (vector field) = " << f2_hat << endl;
	
	cout << "\ndiv(F2) = " << div(f2_hat,x,y,z) <<endl;
	cout << "\ngrad(F) = " << grad(f,x,y,z) <<endl;
	cout << "\ncurl(F2) = " << curl(f2_hat,x,y,z) <<endl;
	
	//cout << "\ndiv(F2) = " << div(y,-x,y*z,x,y,z) <<endl;
	//cout << "\ndiv(F2) = " << div(f2,x,y,z) <<endl;
	//cout << "\ncurl(F2) = " << curl(y,-x,y*z,x,y,z) <<endl;
	//cout << "\ncurl(F2) = " << curl(f2,x,y,z) <<endl;
	return 0; 
}
