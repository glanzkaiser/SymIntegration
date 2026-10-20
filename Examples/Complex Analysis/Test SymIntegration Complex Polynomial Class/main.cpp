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
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;


int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Represents: (1 + 0i) + (2 + 1i)x + (3 + 5i)x^2
	PolynomialComplex poly({{1.0, 0.0}, {2.0, 1.0}, {3.0, 5.0}});

	cout << "Polynomial P(x): ";
	poly.print();

	// Evaluate at z = 2.0 + 3.0i
	complex<double> z(2.0, 3.0);
	complex<double> result = poly.evaluateAt(z);

	cout << "Evaluation at " << z << " is: " << result << endl;
	

	// In Frobenius recurrence relations, you often shift r to (r + n)
	// Let's shift by n = 3: P(x + 3)
	complex<double> s(4.0, 1.0);
	PolynomialComplex P_shifted = poly.shift(s);

	cout << "\nShifted Polynomial P(x + (4+i)): ";
	P_shifted.print(); 
	// Analytical expected: (r+3)^2 - 2(r+3) = r^2 + 6r + 9 - 2r - 6 = r^2 + 4r + 3
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0; 
}
