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
	Symbolic r("r"), k("k"),  P("P"), V("V"), t("t"), CV, C("C"), c0("c0"), Q("Q"), T_lake, t_tenth, t_half, C2, ct, ct2, ct_final, sol_Qt, sol_Qt2, Q0("Q0"), ivpc,  ivpc2;
	
	cout << "\nDSolve for Q' = rk +  P - rc(t) , \t c(t) = Q(t)/V\n" <<endl;
	
	sol_Qt = dsolve( r*k+P-(r/V)*Q,Q,t,r);
	ct = dsolve( r*k+P-(r/V)*Q,Q,t,r)/V ;
	cout << "Q(t) = " << sol_Qt <<endl;	
	cout << "c(t) = " << ct <<endl;	
	cout << "c(0) = " << ct[t==0] <<endl;
	
	ivpc = c0 - ct[t==0] ;	
	cout << "c0 - c(0) = " << ivpc <<endl;

	CV = solve(ivpc,C*(V^(-1))).front().rhs;
	cout << "For ivp c(0) = c0, \nC/V = " << CV <<endl;
	
	cout << "\nTherefore, the concentration is, \nc(t) = " << ct[C*(V^(-1)) == CV] <<endl;
	ct_final = ct[C*(V^(-1)) == CV] ;
	cout << "\nThe limit of c(t) as t -> ∞ \nlim_{t -> ∞} c(t) = " << ct_final[t==INFINITY]<<endl;
	
	sol_Qt2 = dsolve(-(r/V)*Q,Q,t,r);	
	ct2 = dsolve(-(r/V)*Q,Q,t,r)/V;	
	ivpc2 = c0 - ct2[t==0];
	C2 = solve(ivpc2,C).front().rhs;
	cout << "\nIf the pollution stops, then \nQ(t) = " << sol_Qt2 <<endl;
	cout << "\nThe concentration when the pollution stops is \nc(t) = " << ct2 <<endl;
	cout << "\nFor ivp c(0) = c0, \n" << ivpc2 << " = 0" << endl;
	cout << "\nTherefore the ivp solution will be, \nc(t) = " << ct2[C==C2] <<endl;

	t_half = (0.5*c0 - ct2[C==C2] )/c0;
	t_tenth = (0.1*c0 - ct2[C==C2] )/c0;
	
	cout << "\nSet c(t) = 0.5 c0 and solve for t=T, \nT = " << solve(t_half,t).front().rhs<<endl;
	cout << "\nSet c(t) = 0.1 c0 and solve for t=T, \nT = " << solve(t_tenth,t).front().rhs<<endl;
	
	T_lake = solve(t_tenth,t).front().rhs;
	cout << "\nTime to reduce the contamination in each lakes to 10 percent of the original value" << endl;
	cout << "\nLake Superior: \tT = " << T_lake[V==12200, r==65.2]<<endl;
	cout << "\nLake Michigan: \tT = " << T_lake[V==4900, r==158]<<endl;
	cout << "\nLake Erie: \tT = " << T_lake[V==460, r==175]<<endl;
	cout << "\nLake Ontario: \tT = " << T_lake[V==1600, r==209]<<endl;
	
	return 0; 
}
