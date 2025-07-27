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

#define π 3.1415926535897f

using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	Symbolic r("r"),t("t"), C("C"), Q("Q"), ft, Q0("Q0"), ivp1("ivp1");
	
	cout << "\nDSolve for Q' = -rQ \n" <<endl;
	cout << "Q(t) = " << dsolve( -r*Q,Q,t,r) <<endl;	

	ft = dsolve( -r*Q,Q,t,r)[t==0];
	cout << "Q(t=0) = " << ft <<endl;

	cout << "For ivp Q(0) = Q0, \nC = " << solve(ft-Q0,C).front().rhs <<endl;

	cout << "\nSolution for the ivp, \nQ(t) = " << ivp(dsolve( -r*Q ,Q,t,r),t,C, Q0) <<endl;

	// Determining r
	ivp1 = ivp(dsolve( -r*Q ,Q,t,r),t,C, Q0)[t==5730];
	cout << "\nHalf-life problem for carbon-14, the ivp becomes:" << endl;
	cout << ivp1-0.5*Q0 << " = 0" << endl;
	
	//cout << "r for t = 5730 \nr = " << solvenonlinear(exp(-5730*r) - 0.5,r) <<endl;
	cout << "\nDetermining the rate by solving the ivp solution \nr = " << solve(exp(-5730*r) - 0.5,r).front().rhs <<endl;
	
	Symbolic r_c14 = solve(exp(-5730*r) - 0.5,r).front().rhs;	
	cout << "\nDetermining the age of the remains with 20% of carbon-14 \nt = " << solve(exp(-r_c14*t)- 0.2,t).front().rhs <<endl;
	

	// Test
	
	/*Symbolic f = exp(-8*r) - 5;
	
	cout << f.coeff(1,0) << endl;     
	cout << f << endl; 

	cout << f.coeff(SymbolicConstant::e,0) << endl;     // -5
	    
	cout << df(f,r).coeff(r,0) << endl;   // -8 
	cout << ln(0.5) << endl;
	cout << log(2) << endl; // exists
	cout << ln(2) << endl; // does not exist
	
	Symbolic g = ln(-8*r) - 5;
	cout << df(g,r) << endl;   // 1/r
	cout << df(g,r).coeff(r,-1) << endl;   // 1
	cout << df(g,r).coeff(r,0) << endl;
	cout << exp(g) << endl; */

	

	return 0; 
}
