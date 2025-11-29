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

bool isPerfectSquare(int n) {
    if (n < 0) {
        return false;
    }
    int int_sqrt = static_cast<int>(sqrt(n));
    return (int_sqrt * int_sqrt == n);
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic t("t"), y("y"), dx("x'"),dy("y'"),yt, θ("θ"), k("k");
	
	yt = (1+dy*dy)*y-k*k;
	cout << "f(t,y) = " << yt << endl;	

	cout << "\nSolve for y'" <<endl;
	cout << "\ny' = " << solve( yt,dy) <<endl;			// this form is too complicated
	cout << "Take the positive root only." <<endl;
	Symbolic r1 = solve( yt,dy).front().rhs;
	cout << "\ny' = " << r1 <<endl;
	
	Symbolic y0 = k*k*sin(t)*sin(t);
	Symbolic dy_dt = df(y0,t) ;

	cout << "\ndy/dt = " << dy_dt <<endl;
	
	Symbolic dy_dx = (((k*k)/(y) - 1 )^(Symbolic(0.5)))[y==y0]; //
	cout << "\ndy/dx = " << simplify(dy_dx,t) <<endl;	 
	dy_dx = simplify(dy_dx,t) ;
	Symbolic dx_dt =  simplify((dy_dt)/(dy_dx),t);
	cout << "\ndx/dt = " << dx_dt <<endl;	 
	
	Symbolic x = integrate(dx_dt,t);
	cout << "\nx(t) = " << x <<endl;
	cout << "y(t) = " << simplify(y0,t) <<endl;
			
	cout << "\nx(t) = " << x[t==0.5*θ] <<endl;
	cout << "y(t) = " << simplify(y0,t)[t==0.5*θ] <<endl;
	
	cout << "\nThe cycloid has to go through the point (x=1, y=2) at some value of θ, say θ0." << endl;
	cout << "\n1 = " << x[t==0.5*θ] <<endl;
	cout << "2 = " << simplify(y0,t)[t==0.5*θ] <<endl;
	
	cout << "\n" << 2*x[t==0.5*θ] << " = " << simplify(y0,t)[t==0.5*θ] << endl;
	Symbolic fθ = 2 * (k^(Symbolic(-2))) * ( 2*x[t==0.5*θ] - simplify(y0,t)[t==0.5*θ] ) ;
	cout << "f(θ) = " << fθ << endl;

	x =  x[t==0.5*1.4];
	Symbolic f_final = x - 1;
	cout << "k = " << solve(f_final,k).front().rhs << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
