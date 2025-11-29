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

int main(void)
{
		Symbolic x("x");
		int bpower = 10;
		int v[999];
		v[0] = 1;		 
		Symbolic integral; 
		Symbolic d0 = 1;	
		int k = 1, l = 2, c=1;
		for(int i = 1 ; i < bpower  ; i = i+2) // to compute the denominator
		{
			d0 *= i;		
			cout << "\nd0 = "<< d0 << endl;
			
		}
		for(int i = 1 ; i < bpower  - 1; i = i+2)
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
			cout << "\ni = "<< i << endl;
			cout << "\nl = "<< l << endl;
			cout << "\nd0 = "<< d0 << endl;
			cout << "\nk = "<< k << endl;
			cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;
		} 	
			
		int j_d = 1 ;
		for(int i = 1 ; i < (0.5*bpower)+1; i = i+1)
		{
			integral += ( v[i-1] * sin(x) ) / ( d0*(cos(x)^(bpower-j_d)) ) ;
			j_d = j_d+2;
		} 
		cout << "Integral = " << integral << endl;
	return 0; 
}
