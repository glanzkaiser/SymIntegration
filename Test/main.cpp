/*
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

// g++ -o result integral.cpp -lsymintegrationc++

#include <iostream>
#include "symintegrationc++.h"

#define π 3.1415926535897f

using namespace std;

int main(void)
{
	Symbolic x("x"), z("z"), y, y2;
	//y = (x*x*x)+ (x*x);
	y = sin(x);	
	y2 = sin(2*x);	
	cout << "y(x) = " << y << endl;
	cout << "y (1) = " << y[x == 1] << endl;
	
	cout << "f(x,z) = " << pow(x,z) << endl;
	cout << "f(x,3) = " << pow(x,z)[z==3] << endl;
	cout << "f(2,3) = " << pow(x,z)[x==2,z==3] << endl;
	
	y = integrate(y,x);
	
	cout << "integral of sin(x) = " << y << endl;
	cout << "d_{x} (sin(2x)) = " << df(sin(2*x),x) << endl;
	
	y2 = integrate(y2,x);
	
	cout << "integral of sin(2x) = " << y2 << endl;
	cout << "integral of sin(3x) = " << integrate(sin(3*x),x) << endl;
	cout << "integral of sin(0.5x) = " << integrate(sin(0.5*x),x) << endl;
	
	cout << "integral of sin(π/2) = " << sin(0.5*π) << endl;
	cout << "integral of sin(π/4) = " << sin(0.25*π) << endl;
	
	return 0; 
}
