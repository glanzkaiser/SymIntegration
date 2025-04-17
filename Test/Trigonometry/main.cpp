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
	Symbolic x("x");
	
	cout << "\nintegral of sin(x) = " << integrate(sin(x),x)  << endl;
	cout << "integral of sin(3x) = " << integrate(sin(3*x),x) << endl;
	cout << "integral of sin(3x+5) = " << integrate(sin(3*x+5),x) << endl;
	cout << "integral of sin(10x-2) = " << integrate(sin(10*x-2),x) << endl;
	
	cout << "\nintegral of cos(x) = " << integrate(cos(x),x)  << endl;
	cout << "integral of cos(3x) = " << integrate(cos(3*x),x) << endl;
	cout << "integral of cos(3x+5) = " << integrate(cos(3*x+5),x) << endl;
	cout << "integral of cos(10x-2) = " << integrate(cos(10*x-2),x) << endl;
	
	cout << "\nintegral of sec(x) = " << integrate(sec(x),x)  << endl;
	cout << "integral of sec(3x) = " << integrate(sec(3*x),x) << endl;
	cout << "integral of sec(3x+5) = " << integrate(sec(3*x+5),x) << endl;
	cout << "integral of sec(10x-2) = " << integrate(sec(10*x-2),x) << endl;
	
	cout << "\nintegral of csc(x) = " << integrate(csc(x),x)  << endl;
	cout << "integral of csc(3x) = " << integrate(csc(3*x),x) << endl;
	cout << "integral of csc(3x+5) = " << integrate(csc(3*x+5),x) << endl;
	cout << "integral of csc(10x-2) = " << integrate(csc(10*x-2),x) << endl;
	
	cout << "\nintegral of tan(x) = " << integrate(tan(x),x) << endl;
	cout << "integral of tan(3x) = " << integrate(tan(3*x),x) << endl;
	cout << "integral of tan(3x+5) = " << integrate(tan(3*x+5),x) << endl;
	cout << "integral of tan(10x-2) = " << integrate(tan(10*x-2),x) << endl;
	
	cout << "\nintegral of cot(x) = " << integrate(cot(x),x) << endl;
	cout << "integral of cot(3x) = " << integrate(cot(3*x),x) << endl;
	cout << "integral of cot(3x+5) = " << integrate(cot(3*x+5),x) << endl;
	cout << "integral of cot(10x-2) = " << integrate(cot(10*x-2),x) << endl;
	
	return 0; 
}
