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
	Symbolic x("x"), y;

	cout << "\nsin(2) = " << sin(2) << endl;
	cout << "sin(-2) = " << sin(-2) << endl;
	cout << "cos(2) = " << cos(2) << endl;
	cout << "cos(-2) = " << cos(-2) << endl;
	cout << "tan(2) = " << tan(2) << endl;
	cout << "tan(-2) = " << tan(-2) << endl;
	cout << "cot(2) = " << cot(2) << endl;
	cout << "cot(-2) = " << cot(-2) << endl;
	cout << "sec(2) = " << sec(2) << endl;
	cout << "sec(-2) = " << sec(-2) << endl;
	cout << "csc(2) = " << csc(2) << endl;
	cout << "csc(-2) = " << csc(-2) << endl;

	cout << "\nsinh(2) = " << sinh(2) << endl;
	cout << "sinh(-2) = " << sinh(-2) << endl;
	cout << "cosh(2) = " << cosh(2) << endl;
	cout << "cosh(-2) = " << cosh(-2) << endl;
	cout << "tanh(2) = " << tanh(2) << endl;
	cout << "tanh(-2) = " << tanh(-2) << endl;
	cout << "coth(2) = " << coth(2) << endl;
	cout << "coth(-2) = " << coth(-2) << endl;
	cout << "sech(2) = " << sech(2) << endl;
	cout << "sech(-2) = " << sech(-2) << endl;
	cout << "csch(2) = " << csch(2) << endl;
	cout << "csch(-2) = " << csch(-2) << endl;

	cout << "\nsinh(x) = " << sinh(x) << endl;
	cout << "cosh(x) = " << cosh(x) << endl;
	cout << "tanh(x) = " << tanh(x) << endl;
	cout << "coth(x) = " << coth(x) << endl;

	cout << "\nsinh(-x) = " << sinh(-x) << endl;
	cout << "sinh(-x+4) = " << sinh(-x+4) << endl;
	cout << "cosh(-x) = " << cosh(-x) << endl;
	cout << "tanh(-x) = " << tanh(-x) << endl;
	cout << "exp(-x) = " << exp(-x) << endl;
	//cout << "exp(-x) = " << SymbolicConstant::e << endl;
	
	cout << "\nintegral of sinh(x) = " << integrate(sinh(x),x)  << endl;
	cout << "integral of sinh(3x) = " << integrate(sinh(3*x),x) << endl;
	cout << "integral of sinh(3x+5) = " << integrate(sinh(3*x+5),x) << endl;
	cout << "integral of sinh(10x-2) = " << integrate(sinh(10*x-2),x) << endl;
	
	cout << "\nintegral of cosh(x) = " << integrate(cosh(x),x)  << endl;
	cout << "integral of cosh(3x) = " << integrate(cosh(3*x),x) << endl;
	cout << "integral of cosh(3x+5) = " << integrate(cosh(3*x+5),x) << endl;
	cout << "integral of cosh(10x-2) = " << integrate(cosh(10*x-2),x) << endl;
	
	cout << "\nintegral of sech(x) = " << integrate(sech(x),x)  << endl;
	cout << "integral of sech(3x) = " << integrate(sech(3*x),x) << endl;
	cout << "integral of sech(3x+5) = " << integrate(sech(3*x+5),x) << endl;
	cout << "integral of sech(10x-2) = " << integrate(sech(10*x-2),x) << endl;
	
	cout << "\nintegral of csch(x) = " << integrate(csch(x),x)  << endl;
	cout << "integral of csch(3x) = " << integrate(csch(3*x),x) << endl;
	cout << "integral of csch(3x+5) = " << integrate(csch(3*x+5),x) << endl;
	cout << "integral of csch(10x-2) = " << integrate(csch(10*x-2),x) << endl;
	
	cout << "\nintegral of tanh(x) = " << integrate(tanh(x),x) << endl;
	cout << "integral of tanh(3x) = " << integrate(tanh(3*x),x) << endl;
	cout << "integral of tanh(3x+5) = " << integrate(tanh(3*x+5),x) << endl;
	cout << "integral of tanh(10x-2) = " << integrate(tanh(10*x-2),x) << endl;
	
	cout << "\nintegral of coth(x) = " << integrate(coth(x),x) << endl;
	cout << "integral of coth(3x) = " << integrate(coth(3*x),x) << endl;
	cout << "integral of coth(3x+5) = " << integrate(coth(3*x+5),x) << endl;
	cout << "integral of coth(10x-2) = " << integrate(coth(10*x-2),x) << endl;
	
	return 0; 
}
