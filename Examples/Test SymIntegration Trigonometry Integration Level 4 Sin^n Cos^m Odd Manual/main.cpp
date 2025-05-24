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
		Symbolic x("x");
		int bpower = 5; // for cos
		int bpower2 = 6; // for sin
		Symbolic sgn = 1;
		Symbolic integral;
		int m = bpower + bpower2;
		int n = 0.5*(bpower-1);
		int n2 = 0.5*(bpower2-1);
		int j = 0;
		if(bpower > bpower2)
		{
			if(bpower2 % 2 == 0) // even case
			{
				if((bpower+1)/2 % 2 == 0)
				{
					sgn = -1;
				}
				else if((bpower+1)/2 % 2 != 0)
				{
					sgn = 1;
				}
			for(int i = m ; i >= bpower2  ; i = i-2)
			{
				integral += sgn*combinations(n,j)*((sin(x)^i)/(i));
				sgn = -sgn;
				j = j+1;
			} 
			}
			j = 0;
			if(bpower2 % 2 != 0) // odd case
			{
				if((bpower2+1)/2 % 2 == 0)
				{
					sgn = 1;
				}
				else if((bpower2+1)/2 % 2 != 0)
				{
					sgn = -1;
				}
			for(int i = m ; i >= bpower2  ; i = i-2)
			{
				integral += sgn*combinations(n2,j)*((cos(x)^i)/(i));
				sgn = -sgn;
				j = j+1;
			} 	
			}
		}
		else if(bpower <= bpower2)
		{
			if(bpower2 % 2 == 0) // even case
			{
				if((bpower+1)/2 % 2 == 0)
				{
					sgn = -1;
				}
				else if((bpower+1)/2 % 2 != 0)
				{
					sgn = 1;
				}
			for(int i = m ; i >= bpower2  ; i = i-2)
			{
				integral += sgn*combinations(n,j)*((sin(x)^i)/(i));
				sgn = -sgn;
				j = j+1;
			} 
			}
			j = 0;
			if(bpower2 % 2 != 0) // odd case
			{
				if((bpower+1)/2 % 2 == 0)
				{
					sgn = -1;
				}
				else if((bpower+1)/2 % 2 != 0)
				{
					sgn = 1;
				}
			for(int i = m ; i >= bpower2  ; i = i-2)
			{
				integral += sgn*combinations(n,j)*((sin(x)^i)/(i));
				sgn = -sgn;
				j = j+1;
			} 	
			}
		}
		cout << "integral  = "<< integral << endl;	

	return 0; 
}
