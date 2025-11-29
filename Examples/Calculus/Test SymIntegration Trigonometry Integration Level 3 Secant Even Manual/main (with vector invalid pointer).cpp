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
		vector<int> v={1};
		vector<int> w={1};  // to copy vector v to vector w
			 
		Symbolic integral; 
		Symbolic d0 = 1;	
		Symbolic sgn = 1;
		int k = 1, l = 2, c=1;
		for(int i = 1 ; i < bpower  ; i = i+2)
		{
			d0 *= i;		
		}
		cout << "d0 = " << d0 << endl;
		for(int i = 1 ; i < bpower  - 1; i = i+2)
		{
			//v[0] = d0 / (bpower-1); 
			cout << "\ni = "<< i<< endl;
			c = v[0];
			//w.clear();
			for(int j = 1 ; j < i  ; j = j+1)
			{					
				w.push_back(v[j-1]);		
			}			
			int d;
			for(int j = 1 ; j < i  ; j = j+1)
			{		
				d = w[j-1];	
				v[j] = d*l;		
			}
			
			v[1] = c*l;
			v[0] = c*k;

			k= k + 2;
			l = l + 2;
			
			cout << "\nw[0] = "<< w[0]<< endl;
			cout << "w[1] = " << w[1]<< endl;
			cout << "w[2] = " << w[2]<< endl;
			cout << "w[3] = " << w[3]<< endl;
			cout << "w[4] = " << w[4]<< endl;

			cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;
			
		} 	
			
		int j_num = 0 ;
		for(int i = 1 ; i < (0.5*bpower)+1; i = i+1)
		{
			cout << sgn * v[i-1] * sin(x) << endl;
			j_num = j_num+1;
		} 
		
	return 0; 
}
