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

	Symbolic p0 = 0.2; // initital approximation of Φ_{0} (t)
	Symbolic p;
	Symbolic rho = 2.8;

	cout << "u_{0} = "<< p0 << endl;
	for(int i=1; i<=10; i++)
	{
		p = rho*p0*(1-p0);
		p0 = p;
		cout << "u_{" <<i <<"} = "<< p << endl;
	}
		cout << "Equilibrium solution u_{n} = "<< divisiond(rho-1,rho) << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
