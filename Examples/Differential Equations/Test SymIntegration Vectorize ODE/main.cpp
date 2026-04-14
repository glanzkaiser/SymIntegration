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

	Symbolic u1("u1"), u2("u2");
	Symbolic  eq_motion1, eq_motion2,  u1_solution ;

	eq_motion1 = 1*df(df(u1[t],t),t) - 2*(u2-u1) + 3*u1;
	eq_motion2 = 1*df(df(u2[t],t),t) + 2*(u2-u1);
	cout << "\nThe equation of motion : " << endl;
	cout << eq_motion1 << " = 0 " << endl;
	cout << eq_motion2 << " = 0 " << endl;

	Equations F_u2 = solve(eq_motion1,u2);
	cout << "\n" << F_u2.front() << endl;
	
	Symbolic u2_solve = F_u2.front().rhs;
	u2_solve = u2_solve[u1==u1[t]];

	// The equation of motion for mass 2 in u1 terms, all are function of t
	Symbolic eq_motion2_u1terms = 1*df(df(u2_solve[t],t),t) + 2*(u2_solve[t]-u1[t]);
 	
	vector<complex<double>>  v1 = higherorderlineardiffeq_vectorize( eq_motion2_u1terms, u1, t, 4);
	cout << "\nThe vector coefficient of homogeneous ODE with constant coefficients:" << endl;
	printComplexVector(v1);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	

	return 0; 
}
