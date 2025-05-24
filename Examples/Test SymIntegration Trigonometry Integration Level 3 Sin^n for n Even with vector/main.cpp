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
		int bpower = 10; // for sin
		int c = 1;
		int d = bpower;
		Symbolic integral; 
		vector<int> v={1};
		vector<int> w={1};
		
		v[0] = c; 
		v.assign({v[0]});
		w[0] = d; 
		w.assign({w[0]});

		for(int i = bpower ; i >= 2 ; i = i-2)
		{	
			d = d*(i-2);
			c = c*(i-1);
			v.push_back(c);
			w.push_back(d);
			if (i==4)
			{
				integral += ( (c*x)/(d) );
			}
			
			/*cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;
			cout << "v[5] = " << v[5]<< endl;
			cout << "v[6] = " << v[6]<< endl;
			cout << "\nw[0] = "<< w[0]<< endl;
			cout << "w[1] = " << w[1]<< endl;
			cout << "w[2] = " << w[2]<< endl;
			cout << "w[3] = " << w[3]<< endl;
			cout << "w[4] = " << w[4]<< endl;
			cout << "w[5] = " << w[5]<< endl;
			cout << "w[6] = " << w[6]<< endl; */
			
		} 
		int  j = 1;
		for(int i = bpower ; i >= 2 ; i = i-2)
		{
			integral -= v[j-1]*cos(x)*((sin(x))^(i-1))/(w[j-1]) ;
			j = j+1;
		} 
	
		cout << "integral = " << integral << endl;
	return 0; 
}
