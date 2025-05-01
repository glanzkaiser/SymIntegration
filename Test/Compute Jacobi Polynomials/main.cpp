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

long long factorial(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	long long result = 1;
	for (int i = 2; i <= n; ++i) 
	{
		result *= i;
	}
	return result;
}

long long combinations(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factorial(n) / (factorial(r) * factorial(n - r));
}

Symbolic jacobipolynomials(int n, double a, double b, const Symbolic &s){
	if (n < 0 ) 
	{
		return 0; // Invalid input
	}
	Symbolic result = 0;
	for (int i = 0; i <= n; ++i) 
	{
		result += combinations(n,i) * (tgamma(a+b+n+i+1)/tgamma(a+i+1)) * ( (0.5*(s-1))^(i));
	}
	return (tgamma(a+n+1) / (factorial(n)*tgamma(a+b+n+1)) )*result;
}

int main(void)
{
	Symbolic x("x"), m("m"), n("n"), y, z;
	y = cos(x)*cos(x);
	
	cout << "\ntgamma(1) = " << tgamma(1) << endl;
	cout << "\ntgamma(2) = " << tgamma(2) << endl;
	cout << "\ntgamma(3) = " << tgamma(3) << endl;
	cout << "\ntgamma(4) = " << tgamma(4) << endl;
	cout << "\ntgamma(5) = " << tgamma(5) << endl;
	cout << "\ntgamma(6) = " << tgamma(6) << endl;
	cout << "\ntgamma(7) = " << tgamma(7) << endl;
	cout << "\ntgamma(8) = " << tgamma(8) << endl;
	cout << "\ntgamma(9) = " << tgamma(9) << endl;
	cout << "\ntgamma(10) = " << tgamma(10) << endl;
	
	cout << "\nJacobi Polynomials( n=1, α =1, β = 1, z = cos(x) ) = " << jacobipolynomials(1,1,1,cos(x)) << endl;
	cout << "\nJacobi Polynomials( n=2, α =2, β = 3, z = (cos(x))^{2} ) = " << jacobipolynomials(2,2,3,(cos(x))^(2)) << endl;
	
	return 0; 
}
