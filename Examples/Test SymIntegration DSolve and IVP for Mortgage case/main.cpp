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
	Symbolic r("r"),t("t"), y("y"), C("C"), S("S"), k("k"), ft, S0("S0"), ivp1("ivp1");
	
	cout << "\nDSolve for S' = rS - k\n" <<endl;
	cout << "S(t) = " << dsolve( r*S - k,S,t,r) <<endl;	

	ft = dsolve( r*S - k,S,t,r)[t==0];
	cout << "S(t=0) = " << ft <<endl;

	cout << "For ivp S(0) = S0, \nC = " << solve(ft-S0,C).front().rhs <<endl;

	cout << "\nTest ivp, \nS(t) = " << ivp(dsolve( r*S - k,S,t,r),t,C) <<endl;

	ivp1 = ivp(dsolve( r*S - k,S,t,r),t,C)[r==0.06, k ==6000, t==20];
	r = 0.06;
	t = 20;
	cout << "For ivp k = 6000, r = 6%, t = 20, S(t) = 0, \nS(0) = " << solve(ivp1,S0).front().rhs <<endl;


	return 0; 
}
