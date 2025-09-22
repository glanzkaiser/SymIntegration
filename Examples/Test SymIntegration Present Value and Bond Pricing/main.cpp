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

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	double r = 0.055;
	int n = 40;
	double C = 50;
	double M = 1000;

	cout << "\nA 20-year 10% coupon bond with yield 11% and semiannual coupon payments"<< endl;
	cout << "Present Value of Coupon payments = " << presentvalue(C,r,n) << endl;	
	cout << "Present value of par (maturity value) = " << bondpricing(C,M,r,n) - presentvalue(C,r,n) << endl;
	cout << "Bond Price = " << bondpricing(C,M,r,n) << endl;
	
	cout << "\nPrice-Yield Relationship for a 20-year 10% coupon bond"<< endl;
	cout << "\nYield" << " \t  \t " << "Price" << endl;
	for (double i = 0.045; i <= 0.17; i += 0.005  )
	{
		cout << i << " \t  \t " << bondpricing(C,M,i/2,n) << endl; // semiannual coupon bond so we use i/2
	}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
