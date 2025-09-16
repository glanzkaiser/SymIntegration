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

	Symbolic t("t"), ẋ("ẋ"), dy("y'"), θ1("θ1"), θ2("θ2"), θ̇1("θ1'"),  θ̇2("θ2'"), θ̈1("θ1''"), θ̈2("θ2''"), l1("l1"), l2("l2"), m1("m1"), m2("m2"), g("g");
	Symbolic x1, y1, x2, y2, V, T, L, dx1, dy1, dx2, dy2;

	x1 = l1*sin(θ1);
	y1 = -l1*cos(θ1);
	x2 = l1*sin(θ1) + l2*sin(θ2);
	y2 = -l1*cos(θ1) - l2*cos(θ2);

	cout << "x1 = " << x1 << endl;	
	cout << "y1 = " << y1 << endl;	
	cout << "x2 = " << x2 << endl;	
	cout << "y2 = " << y2 << endl;	

	dx1 = df(x1,θ1)*θ̇1;
	dy1 = df(y1,θ1)*θ̇1;
	dx2 = df(x2,θ1)*θ̇1 + df(x2,θ2)*θ̇2;
	dy2 = df(y2,θ1)*θ̇1 + df(y2,θ2)*θ̇2;

	cout << "dx1 = " << dx1 << endl;	
	cout << "dy1 = " << dy1 << endl;	
	cout << "dx2 = " << dx2 << endl;	
	cout << "dy2 = " << dy2 << endl;	
	
	V = m1*g*y1 + m2*g*y2;
	T = 0.5*m1*(dx1*dx1 + dy1*dy1) + 0.5*m2*(dx2*dx2 + dy2*dy2);
	T = simplify(T,θ1);
	T = simplify(T,θ2);
	cout <<"\nPotential energy = " << V << endl;
	cout <<"Kinetic energy = " << T << endl;
	
	L = T-V;
	cout << "Lagrangian = " << L << endl;

	cout << "\ndL / dθ1 = " << df(L,θ1) << endl;
	cout << "\ndL / dθ1' = " << df(L,θ̇1) << endl;
	cout << "\nd/dt (dL / dθ1') = " << df(df(L,θ̇1),θ1)*θ̇1 + df(df(L,θ̇1),θ2)*θ̇2  +  df(L,θ̇1).coeff(θ̇1,1)*θ̈1 + df(L,θ̇1).coeff(θ̇2,1)*θ̈2<< endl;
	
	cout << "\n\ndL / dθ2 = " << df(L,θ2) << endl;
	cout << "\ndL / dθ2' = " << df(L,θ̇2) << endl;
	cout << "\nd/dt (dL / dθ2') = " << df(df(L,θ̇2),θ1)*θ̇1 + df(df(L,θ̇2),θ2)*θ̇2  +  df(L,θ̇2).coeff(θ̇2,1)*θ̈2 +  df(L,θ̇2).coeff(θ̇1,1)*θ̈1<< endl;
	
	Symbolic Final_Eq1 = df(df(L,θ̇1),θ1)*θ̇1 + df(df(L,θ̇1),θ2)*θ̇2  +  df(L,θ̇1).coeff(θ̇1,1)*θ̈1 + df(L,θ̇1).coeff(θ̇2,1)*θ̈2 - df(L,θ1) ;
	Symbolic Final_Eq2 =  df(df(L,θ̇2),θ1)*θ̇1 + df(df(L,θ̇2),θ2)*θ̇2  +  df(L,θ̇2).coeff(θ̇2,1)*θ̈2 +  df(L,θ̇2).coeff(θ̇1,1)*θ̈1 - df(L,θ2);
	cout << "\n\nFirst Equation of motion = \n" << Final_Eq1 << endl;
	cout << "\nSecond Equation of motion = " << Final_Eq2 << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
