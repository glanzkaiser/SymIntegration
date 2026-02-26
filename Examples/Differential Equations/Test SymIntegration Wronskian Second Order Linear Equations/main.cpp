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

using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	Symbolic y("y"),t("t");
	
	cout << "\nWronskian for y'' + 5y'+ 6y = 0 with t_{0} = 0\n" <<endl;
	wronskian(1,5,6,y,t,0);		

	cout << "\nWronskian for 16y'' - 8y'+ 145y = 0 with t_{0} = 0\n" <<endl;
	wronskian(16,-8,145,y,t,0);		

	cout << "\nWronskian for y'' + 4y'+ 4y = 0 with t_{0} = 0\n" <<endl;
	wronskian(1,4,4,y,t,0);		

	cout << "\nWronskian for y'' -  y = 0 with t_{0} = 0\n" <<endl;
	wronskian_fundamentalsetofsolutions(1,0,-1,y,t);		
	return 0; 
}
