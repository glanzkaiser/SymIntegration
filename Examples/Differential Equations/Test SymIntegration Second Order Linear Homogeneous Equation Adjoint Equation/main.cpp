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
	Symbolic x("x"), ν("ν"), α("α");
	
	cout << "\nAdjoint for x^{2}y'' + xy' + (x^{2} - ν^{2})y = 0 (Bessel' equation)\n" <<endl;
	adjointequation( x*x,x,(x*x) - (ν*ν),x);
	
	cout << "\n*****************************************************************\n" <<endl;
	cout << "\nAdjoint for (1-x^{2})y'' - 2xy' + α(α+1)y = 0 (Legendre's equation)\n" <<endl;
	adjointequation( 1 - (x*x),-2*x, α*(α+1),x);
	
	cout << "\n*****************************************************************\n" <<endl;
	cout << "\nAdjoint for y'' - xy' = 0 (Airy's equation)\n" <<endl;
	adjointequation( 1,0, -x, x);
	
	return 0; 
}
