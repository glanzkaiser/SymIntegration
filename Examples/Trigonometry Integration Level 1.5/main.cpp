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
	cout << "int cos(x) * cos(x) = " << integrate(cos(x)*cos(x),x) <<endl;
	cout << "int (cos(x) * cos(x))^{-1} = " << integrate((cos(x)*cos(x))^(-1),x) <<endl;
	cout << "int (cos(x))^{-1} = " << integrate((cos(x))^(-1),x) <<endl;
	cout << "int (cos(x))^{3} = " << integrate((cos(x))^(3),x) <<endl;
	//cout << "int (cos(x))^{1} = " << integrate((cos(x))^(1),x) <<endl;
	//cout << "int (sin(x))^{1} = " << integrate((sin(x))^(1),x) <<endl;
	
	cout << "\nint sin(x) * sin(x) = " << integrate(sin(x)*sin(x),x) <<endl;
	cout << "int (sin(x) * sin(x))^{-1} = " << integrate((sin(x)*sin(x))^(-1),x) <<endl;
	cout << "int (sin(x))^{-1} = " << integrate((sin(x))^(-1),x) <<endl;
	
	cout << "\nint sec(x) * sec(x) = " << integrate(sec(x)*sec(x),x) <<endl;
	cout << "\nint csc(x) * csc(x) = " << integrate(csc(x)*csc(x),x) <<endl;
	
	cout << "\nint (sin(mx) * sin(mx))^{-1} = " << integrate((sin(m*x)*sin(m*x))^(-1),x) <<endl;
	cout << "int (cos(mx) * cos(mx))^{-1} = " << integrate((cos(m*x)*cos(m*x))^(-1),x) <<endl;
	
	return 0; 
}
