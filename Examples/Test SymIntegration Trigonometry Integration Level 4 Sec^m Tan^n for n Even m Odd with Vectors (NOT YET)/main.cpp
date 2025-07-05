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
		int bpowertan = 2;
		int bpowersec = 7;
		int bpower = bpowertan + 1;
		int bpowertotal = bpowertan + bpowersec;
		int k = 1, l=2;
		int k_total = 1, l_total=2;
		vector<int> v={1};
		vector<int> v_total={1};
		vector<int> v_final={1};
		vector<int> mc={1}; // to store the middle coefficient
		vector<int> mc_total={1}; // to store the middle coefficient
		vector<int> mc_final={1}; // to store the middle coefficient that will run for the power of secant
		Symbolic sgn = -1;
		Symbolic integral_numerator, integral_denominator; 
		Symbolic d0 = 2, d1 = 2;
		Symbolic d0_total, d1_total ;
		int j =1;
		int m = (0.5*(bpowertotal+1)) - 1;
		int last_coeff;
		int first_coeff = 3;
		
		// For the coefficient at the numerator 
		for(int i = 1 ; i < (bpower-1)/2  ; i = i+1)
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
					
			cout << "\ni = " << i << endl;

			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v.push_back(mc[ic-1]*(first_coeff+2*(i-1)));	
				}
			d1 = d1*(2*i);
			d0 = d0*(2+2*i);
			v.push_back(last_coeff*(first_coeff+2*(i-1))+k);	
			l = l+2;	
			j = j+1;	
	
			cout << "\nd1 = "<< d1 << endl;
			cout << "\nd0 = "<< d0 << endl;
			cout << "\nk = "<< k << endl;
			/*cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;

			cout << "\nmc[0] = "<< mc[0]<< endl;
			cout << "mc[1] = " << mc[1]<< endl;
			cout << "mc[2] = " << mc[2]<< endl;
			cout << "mc[3] = " << mc[3]<< endl;
			cout << "mc[4] = " << mc[4]<< endl;*/
		} 	
		// Doesn't need to be included, only trial and error
		/*
		j = 1;
		for(int i = 1 ; i < (bpowertotal-1)/2  ; i = i+1)
		{
			if (i >= 2)
			{ 
				for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					mc_total[ic-1] = v_total[ic-1] + v_total[ic];
				}
			}

			k_total = k_total*l_total;
			last_coeff = v_total[j-1]; 
			v_total[0] = v_total[0]*(first_coeff+2*(i-1)); 
			v_total.assign({v_total[0]});			
					
			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v_total.push_back(mc_total[ic-1]*(first_coeff+2*(i-1)));	
				}
			v_total.push_back(last_coeff*(first_coeff+2*(i-1))+k_total);	
			l_total = l_total+2;	
			j = j+1;	
	
		} 	
		cout << "\nFor total power" << endl;
		cout << "\nv[0] = "<< v_total[0]<< endl;
		cout << "v[1] = " << v_total[1]<< endl;
		cout << "v[2] = " << v_total[2]<< endl;
		cout << "v[3] = " << v_total[3]<< endl;
		cout << "v[4] = " << v_total[4]<< endl;
		cout << "v[5] = " << v_total[5]<< endl;
		cout << "v[6] = " << v_total[6]<< endl;
		cout << "v[7] = " << v_total[7]<< endl;

		cout << "\nmc[0] = "<< mc_total[0]<< endl;
		cout << "mc[1] = " << mc_total[1]<< endl;
		cout << "mc[2] = " << mc_total[2]<< endl;
		cout << "mc[3] = " << mc_total[3]<< endl;
		cout << "mc[4] = " << mc_total[4]<< endl;
		cout << "mc[5] = " << mc_total[5]<< endl;
		cout << "mc[6] = " << mc_total[6]<< endl; */
		// We copy the vector v into v_final
		v_final[0] = v[0];
		for(int ic = 1 ; ic < (bpowertotal-1)/2 -2  ; ic = ic+1)
		{
			v_final.push_back(v[ic]);	
			cout << "\nv[ic] = " << v[ic] << endl;
		}
		

		int msec = 3;
		
		// this is the for loop to compute when sec^m x has m > 1
		// for the case where n = 2 it has no mirror pattern
		if (bpowertan == 2) // this is the conditional when int tan(x)^n has n=2
		{
		d0_total = d0;
		d1_total = d1;
		for(int i = 1 ; i <= (bpowersec-1)/2  ; i = i+1)
		{
			for(int ic = 1 ; ic <= ((bpowertotal-1)/2 ) - 2  ; ic = ic+1)
			{
				mc_final[ic-1] = v_final[ic-1] + v_final[ic];
			}
			
			v_final[0] = v_final[0]*(msec-2); 
			v_final.assign({v_final[0]});			
					
			
			for(int ic = 1 ; ic < (msec+bpowertan-1)/2-1 ; ic = ic+1)
			{
				v_final.push_back(mc_total[ic-1]);	
			}
			v_final.push_back(v_final[0]);	
			msec = msec+2;	
			d1_total = d1_total*(2*i);
			d0_total = d0_total*(bpowertan+2*i);
	
			cout << "\nd1 = "<< d1_total << endl;
			cout << "\nd0 = "<< d0_total << endl;
			
			cout << "\nv[0] = "<< v_final[0]<< endl;
			cout << "v[1] = " << v_final[1]<< endl;
			cout << "v[2] = " << v_final[2]<< endl;
			cout << "v[3] = " << v_final[3]<< endl;
			cout << "v[4] = " << v_final[4]<< endl;
			cout << "v[5] = " << v_final[5]<< endl;
			cout << "v[6] = " << v_final[6]<< endl;
			cout << "v[7] = " << v_final[7]<< endl;

			cout << "\nmc[0] = "<< mc_final[0]<< endl;
			cout << "mc[1] = " << mc_final[1]<< endl;
			cout << "mc[2] = " << mc_final[2]<< endl;
			cout << "mc[3] = " << mc_final[3]<< endl;
			cout << "mc[4] = " << mc_final[4]<< endl;
			cout << "mc[5] = " << mc_final[5]<< endl;
			cout << "mc[6] = " << mc_final[6]<< endl;
		} 	
		}
		
		if (bpowertan > 2) // this is the conditional when tan^n x has n>2
		{
		d0_total = d0;
		d1_total = d1;
		for(int i = 1 ; i <= (bpowersec-1)/2  ; i = i+1)
		{
			for(int ic = 1 ; ic <= ((bpowertotal-1)/2 ) - 2  ; ic = ic+1)
			{
				mc_final[ic-1] = v_final[ic-1] + v_final[ic];
			}
			
			v_final[0] = v_final[0]*(msec-2); 
			v_final.assign({v_final[0]});			
					
			
			for(int ic = 1 ; ic < (msec+bpowertan-1)/2-1 ; ic = ic+1)
			{
				v_final.push_back(mc_final[ic-1]*(msec-2));	
			}
			v_final.push_back(v_final[0]);	
			msec = msec+2;	
			d1_total = d1_total*(2*i);
			d0_total = d0_total*(bpowertan+2*i);

			cout << "\nd1 = "<< d1_total << endl;
			cout << "\nd0 = "<< d0_total << endl;
			
			cout << "\nv[0] = "<< v_final[0]<< endl;
			cout << "v[1] = " << v_final[1]<< endl;
			cout << "v[2] = " << v_final[2]<< endl;
			cout << "v[3] = " << v_final[3]<< endl;
			cout << "v[4] = " << v_final[4]<< endl;
			cout << "v[5] = " << v_final[5]<< endl;
			cout << "v[6] = " << v_final[6]<< endl;
			cout << "v[7] = " << v_final[7]<< endl;

			cout << "\nmc[0] = "<< mc_final[0]<< endl;
			cout << "mc[1] = " << mc_final[1]<< endl;
			cout << "mc[2] = " << mc_final[2]<< endl;
			cout << "mc[3] = " << mc_final[3]<< endl;
			cout << "mc[4] = " << mc_final[4]<< endl;
			cout << "mc[5] = " << mc_final[5]<< endl;
			cout << "mc[6] = " << mc_final[6]<< endl;
		} 	
		}

		int j_num = 0 ;
		sgn = 1;
		for(int i = bpowertotal-2 ; i >= 1 ; i = i-2)
		{
			integral_numerator += sgn*v_final[j_num]*((sin(x))^(i));
			sgn = -sgn;
			j_num = j_num+1;
		} 
		sgn = 1;
		j = 1;
		for(int i = bpowertotal-1 ; i >= 0 ; i = i-2)
		{
			integral_denominator += sgn*d0_total*combinations(m,j-1)*((sin(x))^(i));
			sgn = -sgn;
			j = j+1;
		} 
		cout << "d1 total = "<< d1_total << endl;	
		d1_total = d1_total*(bpowersec+1)*(bpowertotal+2);
		cout << "d0 total = "<< d0_total << endl;	
		cout << "d1 total = "<< d1_total << endl;	
		cout << "integral at numerator= "<< integral_numerator << endl;	
		cout << "integral at denominator= "<< integral_denominator << endl;	
		cout << "first two terms = "<< (v_final[0]*ln(sin(x)-1))/(d1_total) - (v_final[0]*ln(sin(x)+1))/(d1_total) << endl;	

	return 0; 
}
