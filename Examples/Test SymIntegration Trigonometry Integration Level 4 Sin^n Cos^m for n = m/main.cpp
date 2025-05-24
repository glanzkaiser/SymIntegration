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
		Symbolic integral;
		Symbolic sgn = 1;
		int bpower = 8; // slow when try power of 177, SymPy is fast, we need to fix it ASAP
		int bpower2 = 8;
		int j = 0;
		int n = (bpower-1)/2;
		// For the coefficient at the numerator of sine with odd power
		if(bpower % 2 != 0) // for odd power
		{
			if((bpower+1)/2 % 2 == 0) // if the (power+1)/2 is even 
			{
				sgn = -1;
			}
			else if((bpower+1)/2 % 2 != 0)//  if the (power+1)/2 is odd 
			{
				sgn = 1;
			}
			for(int i = 1 ; i <= (bpower+1)/2  ; i = i+1)
			{
				integral += sgn*combinations(n,i-1)*((sin(x)^((2*bpower) - j))/((2*bpower) - j));
				sgn = -sgn;
				j = j+2;	
			} 	
		}
		int  k =1;
		int c = 1;
		int l = 1;
		j = 1;
		Symbolic d0 = (pow(2,bpower)) * (2*bpower) ;
		if(bpower % 2 == 0) // for even power
		{
			for(int i = 1 ; i <= (bpower)/2  ; i = i+1)
			{
				d0 = d0*l;
				integral -= ( c*cos(2*x)* ((sin(2*x)^((bpower) - k))) )/(d0);				
				l = bpower - j -1;
				c = c*(bpower - k);
				k = k+2;	
				j = j+2;
				cout << "c = "<< c << endl;	
				cout << "k = "<< k << endl;
				cout << "l = "<< l << endl;	
				cout << "d0 = "<< d0 << endl;
			} 	
				integral = integral + (2*c*x)/(d0) ;
				
		}
		cout << "integral = "<< integral << endl;	

	return 0; 
}
