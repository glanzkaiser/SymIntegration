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
		int bpower = 7;
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
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
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
	
			/*cout << "\ni = "<< i << endl;
			cout << "\nd0 = "<< d0 << endl;
			cout << "\nv[0] = "<< v[0]<< endl;
			cout << "v[1] = " << v[1]<< endl;
			cout << "v[2] = " << v[2]<< endl;
			cout << "v[3] = " << v[3]<< endl;
			cout << "v[4] = " << v[4]<< endl;
			cout << "v[5] = " << v[5]<< endl;
			cout << "v[6] = " << v[6]<< endl;

			cout << "\nmc[0] = "<< mc[0]<< endl;
			cout << "mc[1] = " << mc[1]<< endl;
			cout << "mc[2] = " << mc[2]<< endl;
			cout << "mc[3] = " << mc[3]<< endl;
			cout << "mc[4] = " << mc[4]<< endl;*/
		} 	
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{			
			integral += sgn*v[i-1]*(cos(x)^(2*(i-1)));
			sgn=-sgn;
			
		} 	
		sgn = 1;
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{	
			integral_front = sgn;			
			sgn = -sgn;
		} 	
		cout << "\nintegral = "<< integral_front*ln(cos(x)) + (integral)/(d0*(cos(x)^(bpower-1))) << endl;	
		
	return 0; 
}
