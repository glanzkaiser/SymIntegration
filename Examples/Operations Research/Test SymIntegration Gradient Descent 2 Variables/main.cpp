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
	Symbolic x("x"), y("y"), z("z"), f, x_new, y_new, x_old, y_old, f_new, f_old;
	f = (0.5*x*x) + 2*x + y*y + y + 3;
	cout << "F(x,y) = " << f << endl;
	
	x_old = -5;	
	y_old =  21;

	int N = 51;
	float γ = 0.1;
	for(int i = 1 ; i <= N  ; i = i+1)
	{
		//double fxnew = df(f,x)[x==x_old];
		//double fynew = df(f,y)[y==y_old];
		double fxnew = grad(f,x,y)(0)[x==x_old];
		double fynew = grad(f,x,y)(1)[y==y_old];
	
		x_new = x_old - γ*fxnew;
		y_new = y_old - γ*fynew;

		f_new = f[x==x_new,y==y_new] ;
		f_old = f[x==x_old,y==y_old]; 
		double err = f_new-f_old;
		if (abs(err) < pow(10,-5))
		{
			cout << "The gradient descent converges." << endl;			
			break;
		}

		if (i==N)
		{
			cout << "Maximum number of iteration reached." << endl;			
			break;
		}

		x_old = x_new;
		y_old= y_new;

		cout << "i = " << i << endl;
		cout << "xnew  = " << x_new << endl;
		cout << "ynew  = " << y_new << endl;
		cout << "F(xnew,ynew)  = " << f[x==x_new,y==y_new] << endl;
	}
	//cout << "\ngrad(F) = " << grad(f,x,y)<<endl;
	//cout << "\ngrad(F)[x==10,y==10] = " << grad(f,x,y)[x==10,y==10] <<endl;
	//cout << "\ngrad(F)[x==10,y==10] = " << grad(f,x,y)(0)[x==10] <<endl;
	
	return 0; 
}
