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
	Symbolic x("x"), y("y"), C("C"),u("u"),p("p"), t("t"),α("α");
	
	cout << "\nDSolveseparable for dy/dx = (x^2) / (1-y^2) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(1-(y^2), (x^2),y,x) <<endl;	
		
	cout << "\nDSolveseparable for dy/dx = (y-4x) / (x-y) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(x-y, y-4*x,y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (x^2+xy+y^2) / (x^2) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(x^2,(x^2)+x*y+(y^2),y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (4y-3x) / (2x-y) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(2*x-y, 4*y-3*x,y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (4x+3y) / (2x+y) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(-2*x-y, 4*x+3*y,y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (x+3y) / (x-y) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(x-y, x+3*y,y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (x^2+3xy+y^2) / (x^2) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(x^2,(x^2)+3*x*y+(y^2),y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (x^2 - 3y^2) / (2xy) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(2*x*y,(x^2)-3*(y^2),y,x) <<endl;	
	
	cout << "\nDSolveseparable for dy/dx = (3y^2 - x^2 ) / (2xy) \n" <<endl;
	cout << "f(x,y) = " << dsolveseparable(2*x*y,3*(y^2) - (x^2),y,x) <<endl;	
	
	cout << "\nDSolveseparable for du/dt = -α*u^4\n" <<endl;
	cout << "f(u,t) = " << dsolveseparable(pow(u,p)[p==-4],-α,u,t) <<endl;	

	//cout << "\nf(x,y) = " << fractionintegrate(1,(x*x + 2*x+1),x) <<endl;	
	//cout << "\nf(x,y) = " << fractionintegrate(2*x,(x*x + 2*x+1),x) <<endl;	
	
	return 0; 
}
