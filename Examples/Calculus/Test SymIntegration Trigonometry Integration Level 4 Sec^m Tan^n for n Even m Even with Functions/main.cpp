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

Symbolic integralsecanteven(int n) {
		Symbolic x("x");
		int v[999];
		v[0] = 1;		 
		Symbolic integral; 
		Symbolic d0 = 1;	
		int k = 1, l = 2, c=1;
		for(int i = 1 ; i < n ; i = i+2) // to compute the denominator
		{
			d0 *= i;					
		}
		for(int i = 1 ; i < n  - 1; i = i+2)
		{
			c = v[0];
			int arrsec[999]; // make the size of the array as big as possible
			
			for(int j = 0 ; j < i-1  ; j = j+1)
			{	
				arrsec[j] = v[j]; 
			}			
			int d;
			for(int j = 1 ; j < i  ; j = j+1)
			{		
				d = arrsec[j-1];	
				v[j] = d*l;		
			}
			
			v[0] = c*k;
			v[1] = c*l;
			
			k= k + 2;
			l = l + 2;
		} 	
			
		int j_d = 1 ;
		for(int i = 1 ; i < (0.5*n)+1; i = i+1)
		{
			integral += ( v[i-1] * sin(x) ) / ( d0*(cos(x)^(n-j_d)) ) ;
			j_d = j_d+2;
		} 
	return integral;
}

int main(void)
{
		Symbolic x("x");
		int bpowersec = 4;
		int bpowertan = 6;
		int bpower = bpowersec+bpowertan;
		Symbolic result;
		Symbolic sgn = 1;
		int j = 1;
		int m = bpowertan/2;
		for(int i = bpower ; i >= bpowersec ; i = i-2)
		{
			result += sgn*combinations(m,j-1)*integralsecanteven(i);
			sgn = -sgn;
			j = j+1;
		} 
		cout << "integral = "<< result << endl;
		//cout << "integral with subs = " << result[x==1] << endl;
	return 0; 
}
