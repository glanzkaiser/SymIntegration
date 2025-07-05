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

Symbolic integraltangentodd(int n) {
		Symbolic x("x");
		Symbolic integral;
		Symbolic integral_front;
		Symbolic sgn = 1;
		vector<int> v={1};
		vector<int> v_temp={1};
		vector<int> mc={1}; // to store the middle coefficient
		vector<int> mc_temp={1}; // to store the temporary / new middle coefficient
		int d0 = 2;
		int k = 1, l=2;
		// For the coefficient at the numerator 
		for(int i = 1 ; i <= (n-1)/2  ; i = i+1)
		{
			if (i == 1)
			{
				v[0] = 1; 
				v.assign({v[0]});
			}	

			if (i ==2)
			{		
				mc[0]=1;
				mc.assign({mc[0]});
				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
							
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}
			}
			if (i == 3)
			{ 
				mc[0] = mc[0] + 1; // 
				mc.assign({mc[0]});	

				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}					
			}
			if (i == 4)
			{ 
				mc[0] = mc[0] + 1;	
				mc.assign({mc[0]});
				mc.push_back(mc[0]);
				
				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}					
			}
			if (i >= 5)
			{ 
				mc_temp[0] = mc[0]+1;
				mc_temp.assign({mc_temp[0]});
				for(int ic = 1 ; ic < l  ; ic = ic+1)
				{
					mc_temp[ic] = mc[ic-1] + mc[ic];
				}
				mc[l] = (mc_temp[0]);
				mc[l+1] = (mc_temp[0]);

				mc[0]=mc_temp[0];	
				mc.assign({mc[0]});

				for(int ic = 1 ; ic < l ; ic = ic+1)
				{
					mc.push_back(mc_temp[ic]);
				}
				mc.push_back(mc_temp[0]);
				mc.push_back(0);

				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}

				l = l+1;
			}

			d0 = d0*k;
			k = k+1;	
		} 	
		for(int i = 1 ; i <= (n-1)/2  ; i = i+1)
		{			
			integral += sgn*v[i-1]*(cos(x)^(2*(i-1)));
			sgn=-sgn;
			
		} 	
		sgn = 1;
		for(int i = 1 ; i <= (n-1)/2  ; i = i+1)
		{	
			integral_front = sgn;			
			sgn = -sgn;
		} 	
		integral = integral_front*ln(cos(x)) + (integral)/(d0*(cos(x)^(n-1))) ;	
	return integral;
}

int main(void)
{
		Symbolic x("x");
		int bpowersec = 4;
		int bpowertan = 5;
		int bpower = bpowersec+bpowertan;
		Symbolic result;
		int j = 1;
		int m = bpowersec/2;
		for(int i = bpowertan ; i <= bpower ; i = i+2)
		{
			result += combinations(m,j-1)*integraltangentodd(i);
			j = j+1;
		} 
		cout << "for secant power of " << bpowersec << " and tangent power of " << bpowertan << endl;
		cout << "integral = "<< result<< endl;
	return 0; 
}
