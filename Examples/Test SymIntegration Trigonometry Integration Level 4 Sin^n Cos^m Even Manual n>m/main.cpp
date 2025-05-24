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

Symbolic factorialoddup(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	Symbolic result = 1;
	int k = 1;
	for (int i = 1; i <= n; i=i+1) 
	{
		result = result*k;
		k=k+2;
	}
	return result;
}

Symbolic factorialodd(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	Symbolic result = 1;
	for (int i = 1; i <= n; i=i+2) 
	{
		result *= i;
	}
	return result;
}

int main(void)
{
		Symbolic x("x");
		int bpower = 4; // for cos
		int bpower2 = 8; // for sin
		vector<int> v={1}; // for the numerator for the term of function with sin^(n+1) and above,
		vector<int> v_temp={1};
		vector<int> w={1}; // for the denominator for the term of function with sin^(n+1) and above,
		w[0] = bpower+bpower2; 
		int n = bpower+bpower2-2;
		w.assign({w[0]});
		
		int c = 1;
		int d = bpower2;
		Symbolic integral; 
		vector<int> w1={1}; // for the numerator for the term of function with sin^(n-1) and below,
		vector<int> w2={1}; // for the denominator for the term of function with sin^(n-1) and below,
		
		w1[0] = c; 
		w1.assign({w1[0]});
		w2[0] = d; 
		w2.assign({w2[0]});

		int m = 2;
		int k1;
		int k = -4;	// to do the trick to obtain the middle entries from v[2] to v[i-2]
		// For the coefficient at the numerator 
		for(int i = 1 ; i <= (bpower)/2  ; i = i+1)
		{
			int n = (m/2)-1;
			int d = bpower2+m;
			
			v[0] = 1; 
			v.assign({v[0]});
			int l = 1;
			if (i == 1)
			{
				v[0] = 1; 
				v.assign({v[0]});
			}
			if(i>1 && i<4)
			{
			for(int j=1 ; j < (m)/2; j = j+1)
			{
				int odd = factorialodd(l);
				v.push_back(d*(i-j)*(v_temp[j-1]) + odd);
				l = l+2;
									
			}
			}
			k1 = k;		
			if(i>=4)
			{				
				int c;
				for(int j=1 ; j < 2; j = j+1)
				{
					v.push_back(d*(i-j) + 1);
				}
				for(int j=2 ; j < i-1; j = j+1)
				{
					v.push_back(d*v_temp[j-1] + ((k1)*v_temp[j-1]) + v_temp[j]);	
					c = v_temp[j];
					k1 -=2;
					
				}
				// For the last entry
				int odd = factorialoddup(i-1);
				v.push_back(d*c+ odd);
			}
			
			v_temp[0] = v[0]; 
			v_temp.assign({v_temp[0]});
			for(int j = 1 ; j < (m)/2  ; j = j+1)
			{
					v_temp.push_back(v[j]);	
			}

			m = m+2;
			k = k+2;
			
			/*cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;
			cout << "v[5] = " << v[5]<< endl;
			cout << "v[6] = " << v[6]<< endl;*/

		} 	
		// For the term of function with sin^(n+1) and above, these are the computations needed for this term
		int  j = 1;
		for(int i = 1 ; i <= (bpower/2) - 1 ; i = i+1)
		{	
			w[j] = w[j-1]*n;
			n = n-2;
			j = j+1;
		} 
		j = 1;
		int l = 1;
		
		Symbolic sgn = 1;
		if((bpower)/2 % 2 == 0)
		{
			sgn = 1;
		}
		else if((bpower)/2 % 2 != 0)
		{
			sgn = -1;
		}
		for(int i = 1 ; i <= (bpower/2) ; i = i+1)
		{	
			integral += sgn*v[j-1]*cos(x)*((sin(x))^(bpower+bpower2-l))/(w[j-1]) ;
			j = j+1;			
			l = l+2;
			sgn = -sgn;
		}

		// For the term of function with sin^(n-1) and below, this is the first computation for sin(x)^n
		m = 2;
		for(int i = bpower2 ; i >= 2 ; i = i-2)
		{	
			d = d*(i-2);
			c = c*(i-1);
			w1.push_back(c);
			w2.push_back(d);
			
			m = m+2;
		} 
		// We will update the term coefficients of function with sin^(n-1) and below, this is the second computation for sin(x)^n * cos(x)^m
		j = 1;		
		m = 2;
		for(int i = 1 ; i <= (bpower)/2  ; i = i+1)
		{
			for(int j = 1 ; j <= (bpower2)/2  ; j = j+1)
			{
				w1[j-1] *= m-1;
				w2[j-1] *= bpower2+m;
			}
			m=m+2;
		}
		
		j = 1;
		for(int i = bpower2 ; i >= 2 ; i = i-2)
		{
			integral -= w1[j-1]*cos(x)*((sin(x))^(i-1))/(w2[j-1]) ;
			j = j+1;
		} 
		
		cout << "integral = " << (w1[bpower2/2 - 1] *x)/w2[bpower2/2 - 1]  + integral << endl;
		
	return 0; 
}
