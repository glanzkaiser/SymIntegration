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
//#include <cmath>

#define pi 3.1415926535897f

using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	Symbolic t("t"), v("v"), a("a"), g("g"),  ρ("ρ"),  ρ2("ρ2"), π("π"), μ("μ"), sol_vt;
	
	cout << "\nDSolve for v' = ((-9*μ)/(2*a*a*ρ))*gv+(1-(ρ2/ρ))*g\n" <<endl;
	
	sol_vt = dsolve( ((-9*μ)/(2*a*a*ρ))*v+(1-(ρ2/ρ))*g,v,t);
	cout << "v(t) = " << sol_vt <<endl;		
	
	Symbolic	E("E"), e("e"); 
	Symbolic w = division(4,3)*π*(a^(3))*ρ*g;
	Symbolic Fe = e*E;
	Symbolic B = division(4,3)*π*(a^(3))*ρ2*g;
	cout << "w - Fe - B = 0\n" << w - Fe - B << " = 0" << endl;		
	cout << "e = " << solve(w - Fe - B,e).front().rhs <<endl;		
	
	//cout << "lim_{t -> ∞} v(t) = " << evalf(sol_vt,t,INFINITY) <<endl;		
	return 0; 
}
