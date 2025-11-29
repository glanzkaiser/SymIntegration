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
#include <bits/stdc++.h>
#include <cmath>

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
        
	Symbolic x("x"), θ("θ"), s("s"), λ("λ"), L;

	Symbolic f = divisions(θ,pow(x,θ+1));
	cout <<"f(x;θ) = "<< f << endl;
	/*L = 1;
	for (int i = 1; i < 7 ; ++i)
	{
		L *= divisions(θ,pow(x[i],θ+1)) ;
	}
	cout <<"\nL(x1,x2,...,x10;θ) = "<< L << endl;
	
	Symbolic lnL = ln(L);
	cout <<"\nln L(x1,x2,...,x10;θ) = "<< lnL << endl;*/

	Symbolic ml= mle(f,x,θ,s,6);
	cout << "estimated θ = " << ml<< endl;	
	
	dvec sample = loadVectorFromFile("sample.txt");
	
	// Can we make it automated with for loop for the Equations rules? this is silly
	Equations rules = (   x[1]==sample[0],
					 x[2]==sample[1],
					x[3]==sample[2],
					x[4]==sample[3],
					x[5]==sample[4],
					x[6]==sample[5]);
	
	cout << "estimated θ = " << ml.subst_all(rules) << endl;	

	Symbolic f2 = (exp(-λ))*(λ^(x))/(factorialsym(x));
	Symbolic ml2= mle(f2,x,λ,s,6);
	cout <<"\nf(x;λ) = "<< f2 << endl;
	cout << "estimated λ = " << ml2<< endl;	

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	return 0; 
}
