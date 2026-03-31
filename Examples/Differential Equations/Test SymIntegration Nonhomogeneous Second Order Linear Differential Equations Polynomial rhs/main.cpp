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

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic y("y");
	Polynomial<double> t("t");
	Polynomial<double> p = 4.0*(t^2) - 1 ;
	Polynomial<double> p2 = 5.0*(t^3) + 2.0*(t^2) - 1 ;
	Polynomial<double> p3 = 5.0*(t)  - 10 ;
	Polynomial<double> p4 = 50.0 ;
	
	cout << "\nSolution for y'' - 3y' - 4y = 4t^{2} - 1 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationssolution(1,-3,-4,p);		

	cout << "\nSolution for y'' - 3y' - 4y = 5t^{3} + 2t^{2} - 1 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationssolution(1,-3,-4,p2);		
	
	cout << "\nSolution for y'' - 3y' - 4y = 5t - 10 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationssolution(1,-3,-4,p3);		

	cout << "\nSolution for y'' - 3y' - 4y = 50 \n" <<endl;
	secondorderlineardiffeq_nonhomogeneousequationssolution(1,-3,-4,p4);		

	//cout << "p(x) = " << p << endl;
	//cout << "Diff(p) = " << Diff(p,"x") << endl;
	//cout << "Int(p) = " << Int(p,"x") << endl;
	//cout << "p(x)^2 = " << (p^2) << endl << endl;
		//cout << p.terms.front().first << endl;
	//cout << p.terms.front().second << endl;
	//cout << p.terms.size()<< endl;
	//cout << p.terms.back().second << endl;
	//cout << p.variable << endl;
	
	/*for (int i = 0 ;  i < 4 ; ++i)
	{
		cout << "coefficient = " << p.terms.front().first << endl;
		cout << "degree = " << p.terms.front().second << endl;
	
		p.terms.pop_front();
	}*/ // this is for all my stressful day
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	

	return 0; 
}
