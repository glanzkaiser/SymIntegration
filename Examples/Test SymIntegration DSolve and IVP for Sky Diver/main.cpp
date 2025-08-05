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

#define pi 3.1415926535897f

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
	Symbolic v("v"), t("t"), W("W"), m("m"), g("g"), C1("C1"), C3("C3"), C4("C4"), C6("C6"), C("C"), sol_vt_before, sol_vt_after, yt_before, yt_after;
	Symbolic yt_final_after, yt_final_before, vt_final_after, vt_final_before, ivp_vt_before;
	Symbolic vt_before_10, vt_after_10;
	Symbolic C1_value, C3_value, C4_value, C6_value, W_value, m_value, t_value;
	Symbolic C3_value_final, C6_value_final;
	Symbolic yt_generalsolution_before, yt_generalsolution_after ;
	Symbolic yt_before_C4, yt_after_C6;
	Symbolic yt_before_10, yt_after_10;
	cout << "\nDSolve for v' + (0.75/m)v = -W/m \n" <<endl;
	
	sol_vt_before = dsolve( -(0.75/m)*v - (W/m),v,t);
	sol_vt_after = dsolve( -(12/m)*v - (W/m),v,t);
	
	cout << "Before parachute opens, \nv(t) = " << sol_vt_before[C==C1] << endl;	
	cout << "\nAfter parachute opens, \nv(t) = " << sol_vt_after[t==t-10, C==C3] << endl;	
	
	yt_before = integrate(sol_vt_before,t);
	yt_after = integrate(sol_vt_after,t);
	cout << "\nBefore parachute opens, \ny(t) = " << yt_before[C==C1] + C4 << endl;	
	cout << "\nAfter parachute opens, \ny(t) = " << yt_after[t==t-10, C==C3] + C6 << endl;	
	
	vt_final_before = sol_vt_before[C==C1] ;
	vt_final_after = sol_vt_after[t==t-10, C==C3] ;
	yt_final_before = yt_before[C==C1] + C4 ;
	yt_final_after = yt_after[t==t-10, C==C3] + C6 ;

	ivp_vt_before = ivp(vt_final_before,t,C1,0);
	vt_final_before = sol_vt_before[C==C1, t==0] ;
	cout << "\nC1 = " << solve(vt_final_before-0,C1).front().rhs << endl;	
	cout << "\nBefore parachute opens, \nv(t) = " << ivp_vt_before << endl;	
	
	vt_before_10 = ivp_vt_before[t==10];
	vt_after_10 = vt_final_after[t==10];
	cout << "\nWhen the parachute opens, \nv(10)_{before} = " << vt_before_10 << endl;	
	cout << "\nv(10)_{after} = " << vt_after_10 << endl;	
	cout << "\nC3 = " << solve(vt_before_10-vt_after_10,C3).front().rhs << endl;	

	C1_value = solve(vt_final_before-0,C1).front().rhs;
	C3_value = solve(vt_before_10-vt_after_10,C3).front().rhs ;	
	cout << "\nSo then, \nv(t)_{after} = " << vt_final_after[C3==C3_value] << endl;	
	
	yt_before_C4 = yt_final_before[t==0, C1==C1_value]; 
	C4_value = solve(yt_before_C4-5000,C4).front().rhs;
	cout << "\nThus, \ny(t)_{before} = " << yt_final_before[C1==C1_value] << endl;	
	cout << "\ny(t)_{after} = " << yt_final_after[C3==C3_value] << endl;	
	cout << "\nThe diver starts at 5000 feet, thus \nC4 = " << C4_value << endl;	

	yt_before_10 = yt_final_before[t==10, C1==C1_value, C4==C4_value]; 
	yt_after_10 = yt_final_after[t==10, C3==C3_value]; 
	C6_value = solve(yt_before_10 - yt_after_10,C6).front().rhs;
	cout << "\nThe position of the diver after the parachute opens at 10 seconds, \ny(10)_{before} = " << yt_before_10 << endl;	
	cout << "\ny(10)_{after} = " << yt_after_10 << endl;	
	cout << "\nThus, \nC6 = " << C6_value << endl;	

	W_value = 180;
	m_value = division(180,(9.81*3.28)) ;
	cout << "\nWe have W = 180 lb, g = 9.81 m/s^2, then\nm = " << m_value << endl;	
	
	cout << "\nAs a result, the constants are \nC1 = " << C1_value[m==m_value, W==W_value] << endl;	
	cout << "C3 = " << C3_value[m==m_value, W==W_value] << " = " << evalf(C3_value[m==m_value, W==W_value],1,1)<< endl;	
	cout << "C4 = " << C4_value[m==m_value, W==W_value] << endl;	
	cout << "C6 = " << C6_value[m==m_value, W==W_value] << " = " << evalf(C6_value[m==m_value, W==W_value],1,1) << endl;	

	C3_value_final = evalf(C3_value[m==m_value, W==W_value],1,1)  ;
	C6_value_final = evalf(C6_value[m==m_value, W==W_value],1,1) ;

	// Faster to substitute all variables with Equations rules
	 Equations rules = (C1 == C1_value,
                        		C3 == C3_value_final,
                    			C4 == C4_value,
		        		C6 == C6_value_final,
					W == W_value,
					m == m_value);

	yt_generalsolution_before = yt_final_before.subst_all(rules);
	yt_generalsolution_after = yt_final_after.subst_all(rules);

	// Finding the time till the diver reaches the ground
	cout << "\n*** Newton' Method ***\n" <<endl;
	Symbolic x("x"), f, fd, fp, fpd;
	double pn;
	double p0 = 20*(pi/4);
	int N = 20;
	
	f = yt_generalsolution_after[t==x];
	fd = df(f,x);	
	fp = f[x==p0] ;
	fpd = fd[x==p0] ;
	cout << "\nf(x) = " << f <<endl;
	cout << endl;
	cout << "f'(x) = " << fd <<endl;
	
	cout << endl;
	cout << setw(6) << "n" << "\t\t" << "p_{n}"   << "\n";
	cout << setprecision(14) << setw(6) << "0" << "\t\t" << p0  << "\n";	
	for (int i = 1; i <=N; i++)
	{
		double fp = evalf(f,x,p0);
		double fpd = evalf(fd,x,p0);
		pn = p0 - (fp/fpd);

		cout << setprecision(14) << setw(6) << i << "\t\t" << pn << "\n";
		double err = p0-pn;
		if (abs(err) < pow(10,-5))
		{
			cout << "The procedure was successful." << endl;	
			cout << "*****************************" << endl;			
			break;
		}
		p0 = pn;
	}
	t_value = pn;
	cout << "\nThe time till the diver reaches the ground: \nt = " << t_value << endl;	
	
	cout << "\nThe general solution for the position: \nFor 0 <= t <= 10, \ny(t) = " << yt_generalsolution_before  << endl;	
	cout << "\nFor 10 <= t <= 266.4, \ny(t) = " << yt_generalsolution_after << endl;	
	
	vt_final_before = sol_vt_before[C==C1] ;
	vt_final_after = sol_vt_after[t==t-10, C==C3] ;
	cout << "\nThe general solution for the velocity: \nFor 0 <= t <= 10, \nv(t) = " << vt_final_before.subst_all(rules) << endl;	
	cout << "\nFor 10 <= t <= 266.4, \nv(t) = " << vt_final_after.subst_all(rules) << endl;	
	
	//cout << "\nAt time 266.4 seconds, the diver height is: \ny(266.4)= " << yt_generalsolution_after[t==266.4] << endl;
	cout << "\nThe speed of the sky diver at 10 second: \nv(10) = " << evalf(vt_final_before.subst_all(rules)[t==10],1,1) << endl;	
	cout << "\nThe position of the sky diver at 10 second: \ny(10)_{before} = y(10)_{after} = " << evalf(yt_generalsolution_before[t==10],1,1) << endl;	
	
	return 0; 
}
