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

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	Symbolic x("x"), y("y"), z("z"),α("α"), β("β"), fx, gx;
	// The class Equation represents an equality or substitution
	// It has a two data members of type Symbolic for the lhs and the rhs on an equation respectively
	Equations e1 = (x==3, y==5); // Create a list of equations that can be used in substitution

	fx = 1/pow(x,z)[x==2-3*x,z==division(1,1)];
	gx = 2 + x + pow(x,z)[x==3*x,z==2];
	cout << "f(x) = " << fx << endl;
	cout << "g(x) = " << gx << endl;
		
	cout << "\nfor f(x) = 0. \n" << solve(x+5,x) <<endl;
	cout << "for g(x) = 0. \n" << solve(gx,x) <<endl;
	
	cout << "\ncoefficient x^{0} for g(x) : \n" << gx.coeff(x,0) <<endl;
	cout << "coefficient x for g(x) : \n" << gx.coeff(x,1) <<endl;
	cout << "coefficient x^{2} for g(x) : \n" << gx.coeff(x,2) <<endl;
	
	list A = solve(gx,x);
	cout << "\nOnly the right hand side. \n" <<endl;
	cout << "for g(x) = 0. \n" << A.front().rhs <<endl;
	cout << A.back().rhs <<endl;
	
	Symbolic b = A.front().rhs ;
	Symbolic n = A.back().rhs ;
	//double nd, bd = CastPtr<const Number<double> >(b)->n;
	//nd = CastPtr<const Number<double> >(n)->n;
	//nd = double(CastPtr<const Number<Rational<Number<void> > > >(n)->n);

	//double bd, nd = CastPtr<const Number<double> >(n)->n;
	//bd = CastPtr<const Number<double> >(b)->n;
  	//bd = double(CastPtr<const Number<Rational<Number<void> > > >(b)->n);
	
	//double bd, nd = CastPtr<const Number<double> >(n)->n;
	//cout << nd <<endl;
	
	return 0; 
}
