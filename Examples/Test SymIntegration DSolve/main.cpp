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
	Symbolic Q("Q"),u("u"), r("r"),t("t"), x("x"), y("y"),yt, C("C"), S("S"), k("k"), ft, S0("S0"), ivp1("ivp1");
	
	// 2y' + ty = 2
	yt = -(0.5*t)*y+1;
	cout << "f(t,y) = " << yt << endl;
	
	cout << "\nDSolve for 2y' + ty = 2\n" <<endl;
	cout << "y(t) = " << dsolve( yt,y,t) <<endl;
	
	cout << "\nDSolve for ty' + 2y = 4t^2\n" <<endl;
	cout << "y(t) = " << dsolve( (-2/t)*y+4*t,y,t) <<endl;	
	
	cout << "\nDSolve for y' + 0.5y = 0.5*exp(t/3)\n" <<endl;
	cout << "y(t) = " << dsolve( -0.5*y + 0.5*exp(division(1,3)*t),y,t) <<endl;	
	
	cout << "\nDSolve for y' - 2y = 4-t\n" <<endl;
	cout << "y(t) = " << dsolve( 2*y + 4-t,y,t) <<endl;	

	cout << "\nDSolve for Q' = r/4 - rQ/100\n" <<endl;
	cout << "Q(t) = " << dsolve( r/4 - r*y/100,y,t,r) <<endl;	

	cout << "\nDSolve for S' = rS \n" <<endl;
	cout << "S(t) = " << dsolve( r*S,S,t,r) <<endl;	

	cout << "\nDSolve for S' = rS - k\n" <<endl;
	cout << "S(t) = " << dsolve( r*S - k,S,t,r) <<endl;			

	
	return 0; 
}
