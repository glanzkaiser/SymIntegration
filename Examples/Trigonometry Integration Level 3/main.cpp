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

	cout << "int (sec(x))^{3} = " << integrate((sec(x))^(3),x) <<endl;
	cout << "int (sin(x))^{5} = " << integrate((sin(x))^(5),x) <<endl;
	cout << "int (csc(x))^{7} = " << integrate((csc(x))^(7),x) <<endl;
	cout << "int (cos(x))^{9} = " << integrate((cos(x))^(9),x) <<endl;
	cout << "int (tan(x))^{11} = " << integrate((tan(x))^(11),x) <<endl;
	cout << "int (cot(x))^{13} = " << integrate((cot(x))^(13),x) <<endl;
		
	/*
	cout << "int (sin(x))^{3} = " << integrate((sin(x))^(3),x) <<endl;
	cout << "int (tan(x))^{3} = " << integrate((tan(x))^(3),x) <<endl;
	cout << "int (cot(x))^{3} = " << integrate((cot(x))^(3),x) <<endl;
	cout << "int (sec(x))^{3} = " << integrate((sec(x))^(3),x) <<endl;
	cout << "int (csc(x))^{3} = " << integrate((csc(x))^(3),x) <<endl;
	
	cout << "int (cos(x))^{4} = " << integrate((cos(x))^(4),x) <<endl;
	cout << "int (sin(x))^{4} = " << integrate((sin(x))^(4),x) <<endl;
	cout << "int (tan(x))^{4} = " << integrate((tan(x))^(4),x) <<endl;
	cout << "int (cot(x))^{4} = " << integrate((cot(x))^(4),x) <<endl;
	cout << "int (sec(x))^{4} = " << integrate((sec(x))^(4),x) <<endl;
	cout << "int (csc(x))^{4} = " << integrate((csc(x))^(4),x) <<endl;
	cout << "int (sin(x))^{n} = " << integrate((sin(x))^(n),x) <<endl;
	*/
	
	return 0; 
}
