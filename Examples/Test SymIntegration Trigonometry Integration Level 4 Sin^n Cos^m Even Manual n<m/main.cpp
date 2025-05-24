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
		int bpower = 14; // for cos
		int bpower2 = 10; // for sin
		vector<int> v={1}; // for the numerator for the term of function with sin^(n+1) and above,
		vector<int> v_temp={1};
		int n = bpower+bpower2;
		
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
			cout << "\ni = "<< i<< endl;
			cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;
			cout << "v[5] = " << v[5]<< endl;
			cout << "v[6] = " << v[6]<< endl;

		} 	
		// For the denominator at the term of function with sin^(n+1) and above, these are the computations needed for this term
		vector<int> denom2={1};
		
		denom2[0] = n; 
		denom2.assign({denom2[0]});
		for(int i = 1 ; i < (bpower/2)  ; i = i+1) 
		{	
			n = n-2;
//			denom2[i] = denom2[i-1]*n;
			denom2.push_back(denom2[i-1]*n);
		} 
		cout << "\ndenom2[0] = "<< denom2[0]<< endl;
		cout << "denom2[1] = " << denom2[1]<< endl;
		cout << "denom2[2] = " << denom2[2]<< endl;
		cout << "denom2[3] = " << denom2[3]<< endl;
		cout << "denom2[4] = " << denom2[4]<< endl;
		cout << "denom2[5] = " << denom2[5]<< endl;
		cout << "denom2[6] = " << denom2[6]<< endl;
		
		// To compute the integral with the term of function with sin^(n+1) and above:
		Symbolic sgn = 1;
		if((bpower)/2 % 2 == 0)
		{
			sgn = -1;
		}
		else if((bpower)/2 % 2 != 0)
		{
			sgn = 1;
		}

		int  j = 1;
		int l = 1;
		for(int i = 1 ; i <= (bpower/2) ; i = i+1) // Computing the integral is really slow. Why is this?
		{	
			//integral += sgn*v[j-1]*sin(x)*((cos(x))^(bpower+bpower2-l))/(denom2[j-1]) ;
			j = j+1;			
			l = l+2;
			sgn = -sgn;
		}

		// For the coefficients for the term of function with sin^(n-1) and below, this is the first computation for sin(x)^n
		// w1 is the numerator and w2 is the denominator		
		m = 2;
		for(int i = bpower2 ; i >= 2 ; i = i-2)
		{	
			d = d*(i-2);
			c = c*(i-1);
			w1.push_back(c);
			w2.push_back(d);
			
			m = m+2;
		} 
		/*cout << "\nw1[0] = "<< w1[0]<< endl;
		cout << "w1[1] = " << w1[1]<< endl;
		cout << "w1[2] = " << w1[2]<< endl;
		cout << "w1[3] = " << w1[3]<< endl;
		cout << "w1[4] = " << w1[4]<< endl;
		cout << "w1[5] = " << w1[5]<< endl;
		cout << "w1[6] = " << w1[6]<< endl;
		cout << "\nw2[0] = "<< w2[0]<< endl;
		cout << "w2[1] = " << w2[1]<< endl;
		cout << "w2[2] = " << w2[2]<< endl;
		cout << "w2[3] = " << w2[3]<< endl;
		cout << "w2[4] = " << w2[4]<< endl;
		cout << "w2[5] = " << w2[5]<< endl;
		cout << "w2[6] = " << w2[6]<< endl;*/
		// We will update the term coefficients of function with sin^(n-1) and below, this is the second computation for sin(x)^n * cos(x)^m
		j = 1;		
		m = 2;
		for(int i = 1 ; i <= (bpower)/2  ; i = i+1)
		{
			for(int j = 1 ; j <= (bpower2)/2  ; j = j+1)
			{
				w1[j-1] = w1[j-1]*(m-1);
				w2[j-1] = w2[j-1]*(bpower2+m);
			}
			m=m+2;
		}
		cout << "\nw1[0] = "<< w1[0]<< endl;
		cout << "w1[1] = " << w1[1]<< endl;
		cout << "w1[2] = " << w1[2]<< endl;
		cout << "w1[3] = " << w1[3]<< endl;
		cout << "w1[4] = " << w1[4]<< endl;
		cout << "w1[5] = " << w1[5]<< endl;
		cout << "w1[6] = " << w1[6]<< endl;
		cout << "\nw2[0] = "<< w2[0]<< endl;
		cout << "w2[1] = " << w2[1]<< endl;
		cout << "w2[2] = " << w2[2]<< endl;
		cout << "w2[3] = " << w2[3]<< endl;
		cout << "w2[4] = " << w2[4]<< endl;
		cout << "w2[5] = " << w2[5]<< endl;
		cout << "w2[6] = " << w2[6]<< endl;
		j = 1;
		for(int i = bpower2 ; i >= 2 ; i = i-2)
		{
			//integral -= w1[j-1]*sin(x)*((cos(x))^(i-1))/(w2[j-1]) ;
			j = j+1;
		} 
		
		//cout << "integral = " << (w1[bpower2/2 - 1] *x)/w2[bpower2/2 - 1]  << endl;
		
	return 0; 
}
