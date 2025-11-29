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
	Symbolic x("x"), y("y"), f, xnew, xold;
	f = x*x - x + 3;
	cout << "F(x) = " << f << endl;
	
	xold = 2;	
	
	int N = 80;
	float γ = 0.1;
	for(int i = 1 ; i <= N  ; i = i+1)
	{
		double fxold = df(f,x)[x==xold];
		xnew = xold - γ*fxold;

		double fxnew= df(f,x)[x==xnew];
		double err = fxnew - fxold;
		if (abs(err) < pow(10,-5))
		{
			cout << "The gradient descent converges." << endl;			
			break;
		}
		cout << "i = " << i << endl;
		cout << "xnew  = " << xnew << endl;
		cout << "f(xnew)  = " << f[x==xnew] << endl;
		
		xold = xnew;

	}
	
	return 0; 
}
