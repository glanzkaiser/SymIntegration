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
//#include <cmath>

#define π 3.1415926535897f

using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	Symbolic k("k"),t("t"), C("C"), T("T"), T0("T0"),ft, Troom("Troom"), Q0("Q0"), ivp0("ivp0");
	
	cout << "\nDSolve for T' = -k(T - Troom) \n" <<endl;
	cout << "T(t) = " << dsolve( -k*(T-Troom),T,t,k) <<endl;	

	ft = dsolve( -k*(T-Troom),T,t,k);
	//cout << "\nT(t=0) = " << ft[t==0] <<endl;
	cout << "\nT(t) = " << ft[Troom==70] <<endl;
	// IVP Solution

	// Alternative 1
	//cout << "\nSolution for the ivp T(0) = 200, \nT(t) = " << ivp(ft,t,C, 200) - ft <<endl;
	//ivp0= ivp(ft,t,C, 200) ;
	//cout << "C = " << solve(ivp0-ft,C).front().rhs[Troom==70, t==0] <<endl;

	// Alternative 2
	cout << "\nSolution for the ivp T(0) = 200, \nC = " << solve(ft-200,C).front().rhs <<endl;
	Symbolic C_value = solve(ft-200,C).front().rhs[Troom == 70, t==0];
	cout << "C = " << C_value <<endl;

	// Determining k
	Symbolic ivp_sol = dsolve( -k*(T-Troom),T,t,k)[Troom==70,C==C_value] ;
	cout << "\nT(t) = " << ivp_sol <<endl;
	cout <<"\nDetemining rate k:" << endl;
	cout << "T(1) = " << ivp_sol[t==1] <<endl;
	cout << "k = " << solve(ivp_sol[t==1]-190,k).front().rhs <<endl;
	Symbolic k_value = solve(ivp_sol[t==1]-190,k).front().rhs ;

	cout <<"\nThe coffee reaches temperature T(t) = 150 Fahrenheit at" << endl;
	cout << "t = " << solve(ivp_sol[k==k_value]-150,t).front().rhs <<  " minutes" << endl;


	return 0; 
}
