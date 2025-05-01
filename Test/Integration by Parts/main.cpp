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
	Symbolic x("x"), λ("λ"), Inf("Inf"), y, Ex, Varx, u, dv, v, du;

	y = λ*exp(-λ*x);
	u = x;
	dv = λ*exp(-λ*x);
	du = df(u,x);
	v = integrate(dv,x);
	
	cout << "f(x) = " << y << endl;
	
	for(int i=1;i<=1;i++)
	{
	y = integrate(y,x);
	cout << i << "-st integral of f(x) = " << y << endl;
	}
	cout << "int_{0}^{∞} f(x) = " << y[x==Inf] - y[x==0] <<endl;

	Ex = x*λ*exp(-λ*x);
	
	cout << "\nx f(x) = " << Ex << endl;
	cout << "\nu = " << u << endl;
	cout << "dv = " << dv << endl;
	cout << "du = " << du << endl;
	cout << "v = " << v << endl;
	
	Ex = integrate(Ex,x);
	cout << "\nE(x) = int x f(x) = " << Ex << endl;
	cout << "int_{0}^{∞} x f(x) = " << Ex[x==Inf] - y[x==0] <<endl;

	Varx = x*x*λ*exp(-λ*x);
	Varx = integrate(Varx,x);

	cout << "\nVar(x) = int x^{2} f(x) = " << Varx << endl;
	cout << "int_{0}^{∞} x^{2} f(x) = " << Varx[x==INFINITY, λ==1] - y[x==0] <<endl;

	cout << "\nint x^{-1} = " << integrate(1/x,x) <<endl;
	cout << "\nint exp(3x+5) = " << integrate(exp(3*x+5),x) <<endl;
	cout << "\nint exp(-3x+5) = " << integrate(exp(-3*x+5),x) <<endl;
	cout << "\nint x*x*exp(x) = " << integrate(x*x*exp(x),x) <<endl;
	cout << "\nint x*x*x*exp(x) = " << integrate(x*x*x*exp(x),x) <<endl;
	cout << "\nint x*x*x*x*exp(x) = " << integrate(x*x*x*x*exp(x),x) <<endl;
	cout << "\nint x*x*x*x*x*exp(x) = " << integrate(x*x*x*x*x*exp(x),x) <<endl;

	return 0; 
}
