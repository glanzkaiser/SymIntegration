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

#define π 3.1415926535897f

using namespace std;

int main(void)
{
	Symbolic x("x"), m("m"), n("n"), y;
	
	cout << "\nint cos(mx) * sin(nx) = " << integrate(cos(m*x)*sin(n*x),x) <<endl;
	cout << "\nint cos(mx) * sin(mx) = " << integrate(cos(m*x)*sin(m*x),x) <<endl;
	cout << "\nint cos(5x) * sin(3x) = " << integrate(cos(5*x)*sin(3*x),x) <<endl;
		
	cout << "\nint cos(mx) * cos(nx) = " << integrate(cos(m*x)*cos(n*x),x) <<endl;
	cout << "\nint cos(mx) * cos(mx) = " << integrate(cos(m*x)*cos(m*x),x) <<endl;
	cout << "\nint cos(5x) * cos(3x) = " << integrate(cos(5*x)*cos(3*x),x) <<endl;
	
	cout << "\nint sin(mx) * sin(nx) = " << integrate(sin(m*x)*sin(n*x),x) <<endl;
	cout << "\nint sin(mx) * sin(mx) = " << integrate(sin(m*x)*sin(m*x),x) <<endl;
	cout << "\nint sin(5x) * sin(3x) = " << integrate(sin(5*x)*sin(3*x),x) <<endl;
	
	cout << "\nint 1 = " << integrate(1,x) <<endl;
	cout << "int 2 = " << integrate(2,x) <<endl;
	cout << "int x = " << integrate(x,x) <<endl;
	cout << "int 2x = " << integrate(2*x,x) <<endl;
	
	return 0; 
}
