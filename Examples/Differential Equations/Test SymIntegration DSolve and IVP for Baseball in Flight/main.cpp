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

#define PI 3.1415926535897f

using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

Symbolic division2(Symbolic x, Symbolic y)
{
	return x/y;
}

int main(void)
{
	Symbolic v("v"), t("t"), u("u"), h("h"), w("w"), r("r"), g("g"), A("A"), C1("C1"), C2("C2"), C3("C3"), C4("C4"), C("C");
	Symbolic sol_vt, sol_wt, ivp_vt, ivp_wt, C1_value, C2_value, C3_value, C4_value, vt_for_ivp, wt_for_ivp;
	Symbolic xt,  yt, xt_final, yt_final;
	sol_vt = dsolve( -r*v,v,t);
	sol_wt = dsolve( -g -r*w ,w,t);
	
	cout << "Baseball in flight with dv/dt = -rv and dw/dt = -g - rw "<< endl;	
	cout << "\nv(t) = " << sol_vt[C==C1] << endl;	
	cout << "\nw(t) = " << sol_wt[C==C2] << endl;	
	
	vt_for_ivp = sol_vt[C==C1] ;
	wt_for_ivp = sol_wt[C==C2] ;
	C1_value = solve(vt_for_ivp[t==0]  - (u*cos(A)) , C1).front().rhs;
	C2_value = solve(wt_for_ivp[t==0]  - (u*sin(A)) , C2).front().rhs;
	ivp_vt = ivp(sol_vt,t,C,u*cos(A));
	ivp_wt = ivp(sol_wt,t,C,u*sin(A));
	cout << "\nC1 = " << C1_value << endl;
	cout << "C2 = " << C2_value << endl;
	cout << "\nv(t) = " << ivp_vt << endl;
	cout << "w(t) = " << ivp_wt << endl;

	xt = integrate(ivp_vt,t) + C3;
	yt = integrate(ivp_wt,t) + C4;
	cout << "\nx(t) = " << xt << endl;
	cout << "y(t) = " << yt << endl;

	C3_value = solve(xt[t==0], C3).front().rhs;
	C4_value = solve(yt[t==0] - h, C4).front().rhs;
	cout << "\nC3 = " << C3_value << endl;
	cout << "C4 = " << C4_value << endl;
	xt_final = xt[C3==C3_value];
	yt_final = yt[C4==C4_value];
	
	cout << "\nTherefore, \nx(t) = " << xt_final << endl;
	cout << "y(t) = " << yt_final << endl;

	// Faster to substitute all variables with Equations rules
	Equations rules = (r == 0.2,
                        		u == 125,
                    			h == 3,
		        		g == 32);
	cout << "\nFor r=1/5 s^(-1), u = 125 ft/s , h = 3 ft, and g = 32 ft/s^2, \nx(t) = " << xt_final.subst_all(rules) << endl;
	cout << "y(t) = " << yt_final.subst_all(rules)  << endl;
	Symbolic L = 350, H=10;
	Symbolic f = (1-exp(-0.2*t))*(5*u)*(cos(A)) - L;
	Symbolic t_value = solve(f,t).front().rhs;
	cout << "\nTo clear the wall of distance L=350 ft and height H = 10 ft, \nt  = " << t_value<<endl;
	
	Equations rules2 = (r == 0.2,
                    			h == 3,
					t ==t_value,
		        		g == 32);

	cout << "\ny(t) = 10  = " << yt_final.subst_all(rules2) <<endl;
	//cout << "\ny(t) = " << yt_final.subst_all(rules2)[A==0.6432,u==144.9]  <<endl;

	Symbolic f2 = yt_final.subst_all(rules2) - 10;
	//cout << "\nu  = " << solve(f2,u).front().rhs <<endl;
	
	return 0; 
}
