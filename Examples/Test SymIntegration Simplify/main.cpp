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

	Symbolic t("t"), y("y"), x("x"), k("k");
	
	
	cout << "\n perfect square 4 = " << isPerfectSquare(9) << endl; 
	
	cout << "\nsimplify  = " << simplify( 1/(cos(x)^(2)),x) <<endl;			
	cout << "\nsimplify  = " << simplify( 1/sin(x),x) <<endl;		
	cout << "\nsimplify  = " << simplify( (1-(sin(y)^2)),y) <<endl;	 	
	cout << "\nsimplify  = " << simplify( (1-(cos(y)^2)),y) <<endl;	 	
	cout << "\nsimplify  = " << simplify( 0.5*(exp(x) - exp(-x)),x) <<endl;			
	cout << "\nsimplify  = " << simplify( 0.5*(exp(x) + exp(-x)),x) <<endl;			
	cout << "\nsimplify  = " << simplify( sinh(x)/cosh(x),x) <<endl;			
	cout << "\nsimplify  = " << simplify( cosh(x)/sinh(x),x) <<endl;			
	cout << "\nsimplify  = " << simplify( 1/cosh(x),x) <<endl;			
	cout << "\nsimplify  = " << simplify( 1/sinh(x),x) <<endl;			
	cout << "\nsimplify  = " << simplify( (cosh(x)^2) - (sinh(x)^2) ,x) <<endl;			
	cout << "\nsimplify  = " << simplify( 1 + (sinh(x)^2) ,x) <<endl;			
	cout << "\nsimplify  = " << simplify(sin(acos(x)) ,x) <<endl;			
	cout << "\nsimplify  = " << simplify(cos(asin(x)) ,x) <<endl;			
	cout << "\nsimplify  = " << simplify( sinh(x) + cosh(x),x) <<endl;			
	
	cout << "\nsimplify  = " << simplify( ((sin(x))^(4)) + ((sin(x))^(2)) ,x) <<endl;	 // Simplify THIS!!!!! use your brain..  correct answer..		
	
	cout << sqrt( (((sin(x))^(4)) + ((sin(x))^(2)))*(sin(x)^(2)))  << endl;
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
