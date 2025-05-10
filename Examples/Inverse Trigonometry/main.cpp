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
	
	cout << "pi = " << SymbolicConstant::pi<< endl;
	cout << "asin(-0.3) = " << asin(-0.3) << endl;
	cout << "asin(0.3) = " << asin(0.3) << endl;
	cout << "acos(-0.5) = " << acos(-0.5) << endl;
	cout << "acos(0.5) = " << acos(0.5) << endl;
	cout << "asin(1) = " << asin(1) << endl;
	cout << "acos(1) = " << acos(1) << endl;
	cout << "atan(1) = " << atan(1) << endl;
	
	cout << "\nintegral of asin(x) = " << integrate(asin(x),x)  << endl;
	cout << "integral of asin(3x) = " << integrate(asin(3*x),x) << endl;
	cout << "integral of asin(3x+5) = " << integrate(asin(3*x+5),x) << endl;
	cout << "integral of asin(10x-2) = " << integrate(asin(10*x-2),x) << endl;
	
	cout << "\nderivative of asin(x) = " << df(asin(x),x)  << endl;
	cout << "derivative of asin(3x) = " << df(asin(3*x),x)  << endl;
	cout << "derivative of asin(3x+5) = " << df(asin(3*x+5),x)  << endl;
	cout << "derivative of asin(10x-2) = " << df(asin(10*x-2),x)  << endl;
	
	cout << "\nintegral of acos(x) = " << integrate(acos(x),x)  << endl;
	cout << "integral of acos(3x) = " << integrate(acos(3*x),x) << endl;
	cout << "integral of acos(3x+5) = " << integrate(acos(3*x+5),x) << endl;
	cout << "integral of acos(10x-2) = " << integrate(acos(10*x-2),x) << endl;
	
	cout << "\nderivative of acos(x) = " << df(acos(x),x)  << endl;
	cout << "derivative of acos(3x) = " << df(acos(3*x),x)  << endl;
	cout << "derivative of acos(3x+5) = " << df(acos(3*x+5),x)  << endl;
	cout << "derivative of acos(10x-2) = " << df(acos(10*x-2),x)  << endl;

	cout << "\nintegral of atan(x) = " << integrate(atan(x),x)  << endl;
	cout << "integral of atan(3x) = " << integrate(atan(3*x),x) << endl;
	cout << "integral of atan(3x+5) = " << integrate(atan(3*x+5),x) << endl;
	cout << "integral of atan(10x-2) = " << integrate(atan(10*x-2),x) << endl;
	
	cout << "\nderivative of atan(x) = " << df(atan(x),x)  << endl;
	cout << "derivative of atan(3x) = " << df(atan(3*x),x)  << endl;
	cout << "derivative of atan(3x+5) = " << df(atan(3*x+5),x)  << endl;
	cout << "derivative of atan(10x-2) = " << df(atan(10*x-2),x)  << endl;

	cout << "\nintegral of acot(x) = " << integrate(acot(x),x)  << endl;
	cout << "integral of acot(3x) = " << integrate(acot(3*x),x) << endl;
	cout << "integral of acot(3x+5) = " << integrate(acot(3*x+5),x) << endl;
	cout << "integral of acot(10x-2) = " << integrate(acot(10*x-2),x) << endl;
	
	cout << "\nderivative of acot(x) = " << df(acot(x),x)  << endl;
	cout << "derivative of acot(3x) = " << df(acot(3*x),x)  << endl;
	cout << "derivative of acot(3x+5) = " << df(acot(3*x+5),x)  << endl;
	cout << "derivative of acot(10x-2) = " << df(acot(10*x-2),x)  << endl;

	cout << "\nintegral of asec(x) = " << integrate(asec(x),x)  << endl;
	cout << "integral of asec(3x) = " << integrate(asec(3*x),x) << endl;
	cout << "integral of asec(3x+5) = " << integrate(asec(3*x+5),x) << endl;
	cout << "integral of asec(10x-2) = " << integrate(asec(10*x-2),x) << endl;
	
	cout << "\nderivative of asec(x) = " << df(asec(x),x)  << endl;
	cout << "derivative of asec(3x) = " << df(asec(3*x),x)  << endl;
	cout << "derivative of asec(3x+5) = " << df(asec(3*x+5),x)  << endl;
	cout << "derivative of asec(10x-2) = " << df(asec(10*x-2),x)  << endl;

	cout << "\nintegral of acsc(x) = " << integrate(acsc(x),x)  << endl;
	cout << "integral of acsc(3x) = " << integrate(acsc(3*x),x) << endl;
	cout << "integral of acsc(3x+5) = " << integrate(acsc(3*x+5),x) << endl;
	cout << "integral of acsc(10x-2) = " << integrate(acsc(10*x-2),x) << endl;
	
	cout << "\nderivative of acsc(x) = " << df(acsc(x),x)  << endl;
	cout << "derivative of acsc(3x) = " << df(acsc(3*x),x)  << endl;
	cout << "derivative of acsc(3x+5) = " << df(acsc(3*x+5),x)  << endl;
	cout << "derivative of acsc(10x-2) = " << df(acsc(10*x-2),x)  << endl;
	return 0; 
}
