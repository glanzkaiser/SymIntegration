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
#include <vector>
#define π 3.1415926535897f

using namespace std;

// Factorial and combinations
Symbolic factorial(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	Symbolic result = 1;
	for (int i = 2; i <= n; ++i) 
	{
		result *= i;
	}
	return result;
}

Symbolic combinations(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factorial(n) / (factorial(r) * factorial(n - r));
}

int main(void)
{
		Symbolic u("u"), x("x");
		int bpowersec = 3;
		int bpowertan = 5;
		int bpower = bpowersec+bpowertan;
		Symbolic result;
		int j = 1;
		int m = (bpowertan-1)/2;
		Symbolic sgn = 1;
		for(int i = bpowersec-1 ; i <= bpower-2 ; i = i+2)
		{
			result += sgn*combinations(m,j-1)*(u^i);
			j = j+1;
			sgn = -sgn;			
		} 
		cout << "for secant power of " << bpowersec << " and tangent power of " << bpowertan << endl;
		cout << "we will compute the integral of "<< result << endl;
		Symbolic f = integrate(result,u);
		cout << "integral = "<< f << endl;
		cout << "Substitute back, integral = "<< f[u==(1/cos(x))] << endl;
	return 0; 
}
