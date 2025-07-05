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

Symbolic numeratorintegralsecantodd(int n) {
		Symbolic x("x");
		int k = 1, l=2;
		vector<int> v={1};
		vector<int> mc={1}; // to store the middle coefficient
		Symbolic sgn = -1;
		Symbolic integral_numerator; 
		Symbolic d0 = 2, d1 = 2;
		int j =1;
		int m = n-(0.5*(n+1));
		int last_coeff;
		int first_coeff = 3;
		
		// For the coefficient at the numerator of sine with odd power
		for(int i = 1 ; i < (n-1)/2  ; i = i+1)
		{
			if (i >= 2)
			{ 
				for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					mc[ic-1] = v[ic-1] + v[ic];
				}
			}

			k = k*l;
			last_coeff = v[j-1]; 
			v[0] = v[0]*(first_coeff+2*(i-1)); 
			v.assign({v[0]});			
					
			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v.push_back(mc[ic-1]*(first_coeff+2*(i-1)));	
				}
			d1 = d1*(2*i);
			d0 = d0*(2+2*i);
			v.push_back(last_coeff*(first_coeff+2*(i-1))+k);	
			l = l+2;	
			j = j+1;	
		} 	
		int j_num = 0 ;
		for(int i = n-2 ; i >= 1 ; i = i-2)
		{
			integral_numerator += sgn*v[j_num]*((sin(x))^(i));
			sgn = -sgn;
			j_num = j_num+1;
		} 
	return integral_numerator;
}

Symbolic denominatorintegralsecantodd(int n) {
		Symbolic x("x");
		int k = 1, l=2;
		vector<int> v={1};
		vector<int> mc={1}; // to store the middle coefficient
		Symbolic sgn = 1;
		Symbolic integral_denominator; 
		Symbolic d0 = 2, d1 = 2;
		int j =1;
		int m = n-(0.5*(n+1));
		int last_coeff;
		int first_coeff = 3;
		
		// For the coefficient at the numerator of sine with odd power
		for(int i = 1 ; i < (n-1)/2  ; i = i+1)
		{
			if (i >= 2)
			{ 
				for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					mc[ic-1] = v[ic-1] + v[ic];
				}
			}

			k = k*l;
			last_coeff = v[j-1]; 
			v[0] = v[0]*(first_coeff+2*(i-1)); 
			v.assign({v[0]});			
					
			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v.push_back(mc[ic-1]*(first_coeff+2*(i-1)));	
				}
			d1 = d1*(2*i);
			d0 = d0*(2+2*i);
			v.push_back(last_coeff*(first_coeff+2*(i-1))+k);	
			l = l+2;	
			j = j+1;	
		} 	
		j = 1;
		for(int i = n-1 ; i >= 0 ; i = i-2)
		{
			integral_denominator += sgn*d0*combinations(m,j-1)*((sin(x))^(i));
			sgn = -sgn;
			j = j+1;
		} 
	return integral_denominator;
}

Symbolic logtermsintegralsecantodd(int n) {
		Symbolic x("x");
		int k = 1, l=2;
		vector<int> v={1};
		vector<int> mc={1}; // to store the middle coefficient
		Symbolic sgn = -1;
		Symbolic integral_denominator; 
		Symbolic d0 = 2, d1 = 2;
		int j =1;
		int m = n-(0.5*(n+1));
		int last_coeff;
		int first_coeff = 3;
		
		// For the coefficient at the numerator of sine with odd power
		for(int i = 1 ; i < (n-1)/2  ; i = i+1)
		{
			if (i >= 2)
			{ 
				for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					mc[ic-1] = v[ic-1] + v[ic];
				}
			}

			k = k*l;
			last_coeff = v[j-1]; 
			v[0] = v[0]*(first_coeff+2*(i-1)); 
			v.assign({v[0]});			
					

			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v.push_back(mc[ic-1]*(first_coeff+2*(i-1)));	
				}
			d1 = d1*(2*i);
			d0 = d0*(2+2*i);
			v.push_back(last_coeff*(first_coeff+2*(i-1))+k);	
			l = l+2;	
			j = j+1;	
		} 	
		d1 = d1*(n-1);
	return  (-v[0]*ln(sin(x)-1))/(d1) + (v[0]*ln(sin(x)+1))/(d1) ;
}

int main(void)
{
		Symbolic x("x");
		int bpowersec = 3;
		int bpowertan = 6;
		int bpower = bpowersec+bpowertan;
		Symbolic sgn = 1;
		Symbolic result;
		int j = 1;
		int m = bpowertan/2;
		for(int i = bpower ; i >= bpowersec ; i = i-2)
		{
			result += sgn*combinations(m,j-1)*( (numeratorintegralsecantodd(i)/denominatorintegralsecantodd(i)) + logtermsintegralsecantodd(i));
			sgn = -sgn;
			j = j+1;
		} 
		//cout << "integral at numerator= "<< numeratorintegralsecantodd(3)<< endl;	
		//cout << "integral at denominator= "<< denominatorintegralsecantodd(3)<< endl;	
		//cout << "log terms= "<< logtermsintegralsecantodd(3)<< endl;	
		cout << "for secant power of " << bpowersec << " and tangent power of " << bpowertan << endl;
		cout << "integral = "<< result << endl;	
		//cout << "integral with subs = " << result[x==1] << endl;
	return 0; 
}
