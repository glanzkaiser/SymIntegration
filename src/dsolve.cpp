/*
   Thank you Freya the Goddess, Sentinel, Berlin, all Nature, and Mother Mary from Catholic Divine
*/

#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE

#define pi  3.1415926535897

// for std::chrono
#include<iostream>
#include<vector>
#include <fstream>
#include <cmath> // Required for round and pow
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()
#include <map>
#include <algorithm> // For std::max_element,  std::sort, std::reverse ,  std::for_each
#include <numeric> // For std::accumulate
#include <random> // For random number generation
#include <complex>

double roundToDecimal(double value, int places) 
{
	double multiplier = std::pow(10.0, places);
	return std::round(value * multiplier) / multiplier;
}


Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, mu, C("C");
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, r;
		// Case 1 : ay' + ty = b
		eq = (a*x*y + d).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = exp(integrate(-ap*x,x));
		dsol = (integrate(mu*dp,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : aty' + by = ct^2	/ y' = -(a/t)y + b*t	
		eq = ((a/x)*y + d*x).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = x^(-ap);
		dsol = (integrate(mu*dp*x,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : ay' + by = c*exp(d*t)	/ y' = -ay + b*exp(c*t) 	
		eq = (a*y + b*exp(c*x)).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(-ap*x);
		dsol = (integrate(mu*bp*exp(cp*x),x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : ay' + by = c*t+d	/ y' = -ay + b*t + c 	
		eq = (a*y + b*x+c).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(-ap*x);
		dsol = (integrate(mu*(bp*x+cp),x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 5 :  ay' + by = c / y' = ay +by +c
		eq = (a*y+b).match(fx, (a,b)); // pow(u,r)[r==4]
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
		mu = exp(integrate(-ap,x));
		dsol = (integrate(bp*mu,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 6 :  y' + ay = 0
		eq = (a*y).match(fx, (a,b)); // pow(u,r)[r==4]
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		mu = exp(integrate(-ap,x));
		dsol = (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}
	}
	return dsol;
}

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x, const Symbolic &z) // for 1st order ODE with rate r or k
{
	Symbolic dsol, mu, gs,  C("C"), k("k");
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d;
		// Will work for this model: y' + ry/b = r/c	/ y' = r/b - ry/c 
		// y' = ry - k	
		// k = (fx.coeff(z,0)).coeff(y,0)	
		mu = exp(-fx.coeff(y*z,1)*x*z);
		gs = z*(fx.coeff(z,1)).coeff(y,0);
		dsol = (integrate(mu*(gs + (fx.coeff(z,0)).coeff(y,0)),x))/(mu) + (C)/(mu);
		
		/*// Case 1 : y' - ry = -k	/ y' = ry - k 	Useless this cannot be matched. Still keep it for now maybe will be of use one day.
		eq = (a*y*z - b*k).match(fx, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i,a), bp = rhs(*i, b);
		mu = exp(ap*z*x);
		dsol = (integrate(mu*(bp*z),x))/(mu) + (C)/(mu);
		
		return 8;
		
		} catch(const SymbolicError &se) {}
		}*/
	}
	return dsol;
}

Symbolic dsolvelogistic(const Symbolic &fx, const Symbolic &y, const Symbolic &y0, const Symbolic &t, const Symbolic &r, const Symbolic &K, const Symbolic &T) // for 1st order ODE logistic growth
{
	Symbolic dsol;
 	
	if(fx.coeff(y,1) == r && fx.coeff(y,2) == -r*(K^(-1)))
 	{
		dsol = (y0*K)/(y0 + (K - y0)*exp(-r*t));
	}
	if(fx.coeff(y,1) == -r && fx.coeff(y,2) == r*(T^(-1)))
 	{
		dsol = (y0*T)/(y0 + (T - y0)*exp(r*t));
	}
	if(fx.coeff(y,1) == -r && fx.coeff(y,2) == r*(T^(-1)) + r*(K^(-1)))
 	{
		cout << "For T < y < K and y > K"<< endl;
		cout << "y(t) = " << (y0*K)/(y0 + (K - y0)*exp(-r*t)) << endl;
		cout << "For 0 < y < T, \ny(t) ="<< endl;
		dsol = (y0*T)/(y0 + (T - y0)*exp(r*t));
	}
	return dsol;
}

Symbolic ivp(const Symbolic &fx, const Symbolic &x, const Symbolic &c, const Symbolic &so)
{
	Symbolic ivpsol, f0,  C("C");
 
	if(fx != 0)
 	{
		// Will work for this model: dS/dt = rS-k	
		f0 = fx[x==0];
		C = solve(f0-so,c).front().rhs ;
		ivpsol = fx[c==C];	
	}
	return ivpsol;
}

Symbolic dsolveseparable(const Symbolic &fdy, const Symbolic &fdx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, dsol_y, dsol_x, C("C"), c_y;
 	
	if(fdx == 0 && fdy !=0)
 	{
		dsol = C;
	}
	else if(df(fdx,x) == 0 && fdy !=0)
 	{
		c_y = integrate(fdy,y).coeff(y);
		dsol = integrate(fdy,y) - integrate(fdx,x)  - C;
	}
	else if(fdx != 0 && fdy != 0 && fdy.coeff(x,1) == 0 && fdx.coeff(y,1) == 0 && df(fdx,x) != 0 )
 	{
		dsol = integrate(fdy,y) - integrate(fdx,x)  - C; 		
		
	}
	else if(fdy.coeff(x,1) !=0 && fdx.coeff(y,1) !=0 && fdy.coeff(x,2) ==0 && fdx.coeff(y,2) ==0 )
	{
		dsol_y = fdy/x;
		dsol_x = fdx/x ;		
		dsol = dsol_x / dsol_y;
		//dsol = dsol[y*(x^-1)==x];
		//dsol = dsol -(x*dsol_y[y*(x^-1)==x] / dsol_y[y*(x^-1)==x] );
		//dsol = 1/dsol;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x] - x*dsol_y[y*(x^-1)==x];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
		
	}
	else if(fdy.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 )
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}
	else if((fdy)/(y^(-1)) == 1 && x.coeff(x,-1) != 0) // 2tv' - v = 0 
	{
		double coeff_x = x.coeff(x,-1);
		dsol = C*(x^(coeff_x));		
	}
	/*else if(fdx.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 && fdy.coeff(x,1) !=0 ) // is this necessary?
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}*/
	return dsol;
}

void secondorderlineardiffeq_dsolve(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			y1 = exp((-b/(2*a))*x) *(cos(D_real*x) + SymbolicConstant::i*sin(D_real*x));
			y2 = exp((-b/(2*a))*x) *(cos(D_real*x) - SymbolicConstant::i*sin(D_real*x));

			cout <<"\nThe general solution is:" << endl;
			cout << "\ny_{1} (t) = " << y1 << endl;
			cout << "\ny_{2} (t) = " << y2 << endl;
		}
	}
}

void secondorderlineardiffeq_ivpsolution(double a, double b, double c, const Symbolic &y, const Symbolic &x, double t0, double y0, double dy0)
{
	Symbolic yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

			double c1_ans = y0;
			double c2_ans = dy0 - (r1*c1_ans);
			cout <<"\nThe solution for the initial value problem is:" << endl;
			yt = yt[c1 == c1_ans, c2 == c2_ans] ;
			cout << yt << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

			double c1_ans = divisiond(dy0-(r2*y0),r1-r2)*exp(-r1*t0);
			double c2_ans = divisiond((y0*r1)-dy0,r1-r2)*exp(-r2*t0);
			cout <<"\nThe solution for the initial value problem is:" << endl;
			yt = yt[c1 == c1_ans, c2 == c2_ans] ;
			cout << yt << endl;

			Symbolic dyt = df(yt,x);
			cout << "\ny' = " << dyt << endl;

			double tm = NewtonRaphson(dyt,x,0);

			Equations rules = (  SymbolicConstant::e == exp(1), SymbolicConstant::i == sqrt(-1));
			yt = yt.subst_all(rules);
			double ym = yt[x==tm];

			cout << "\nCritical value t_{m} = " << tm << endl;
			cout << "\nMaximum value y_{m} = " << ym << endl;
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			
		
			y1 = exp((-b/(2*a))*x) *(cos(D_real*x) + SymbolicConstant::i*sin(D_real*x));
			y2 = exp((-b/(2*a))*x) *(cos(D_real*x) - SymbolicConstant::i*sin(D_real*x));

			cout <<"\nThe complex-valued general solution is:" << endl;
			cout << "\ny_{1} (t) = " << y1 << endl;
			cout << "\ny_{2} (t) = " << y2 << endl;

			Symbolic ut = exp((-b/(2*a))*x) *(cos(D_real*x));
			Symbolic vt = exp((-b/(2*a))*x) *(sin(D_real*x));
			Symbolic d_ut = df(ut,x);
			Symbolic d_vt = df(vt,x);

			cout << "\nu(t) = " << ut << endl;
			cout << "\nv(t) = " << vt << endl;

			double a11 = ut[x==t0];
			double a12 = vt[x==t0];
			double a21 = d_ut[x==t0];
			double a22 = d_vt[x==t0];
			
			vector<vector<double>> A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			A[0][0] = a11;
			A[0][1] = a12;
			A[1][0] = a21;
			A[1][1] = a22;
			vec_b[0][0] = y0;
			vec_b[1][0] = dy0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(A,vec_b,c_solution);
			//printVector(c_solution);

			Symbolic y_solution = exp((-b/(2*a))*x) * (c_solution[0]*(cos(D_real*x)) + c_solution[1]*(sin(D_real*x)));
			cout <<"\nThe real-valued initial value problem solution is:" << endl;
			cout << "\ny (t) = " << y_solution << endl;
			
		}
	}
	 
}

void secondorderlineardiffeq_springmasssystem(double a, double b, double c, const Symbolic &y, const Symbolic &x, double y0, double dy0)
{
	Symbolic yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	double t0 = 0;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

			double c1_ans = y0;
			double c2_ans = dy0 - (r1*c1_ans);
			cout <<"\nThe solution for the initial value problem is:" << endl;
			yt = yt[c1 == c1_ans, c2 == c2_ans] ;
			cout << yt << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

			double c1_ans = divisiond(dy0-(r2*y0),r1-r2)*exp(-r1*t0);
			double c2_ans = divisiond((y0*r1)-dy0,r1-r2)*exp(-r2*t0);
			cout <<"\nThe solution for the initial value problem is:" << endl;
			yt = yt[c1 == c1_ans, c2 == c2_ans] ;
			cout << yt << endl;

			Symbolic dyt = df(yt,x);
			cout << "\ny' = " << dyt << endl;

			double tm = NewtonRaphson(dyt,x,0);

			Equations rules = (  SymbolicConstant::e == exp(1), SymbolicConstant::i == sqrt(-1));
			yt = yt.subst_all(rules);
			double ym = yt[x==tm];

			cout << "\nCritical value t_{m} = " << tm << endl;
			cout << "\nMaximum value y_{m} = " << ym << endl;
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			
		
			y1 = exp((-b/(2*a))*x) *(cos(D_real*x) + SymbolicConstant::i*sin(D_real*x));
			y2 = exp((-b/(2*a))*x) *(cos(D_real*x) - SymbolicConstant::i*sin(D_real*x));

			Symbolic ut = exp((-b/(2*a))*x) *(cos(D_real*x));
			Symbolic vt = exp((-b/(2*a))*x) *(sin(D_real*x));
			Symbolic d_ut = df(ut,x);
			Symbolic d_vt = df(vt,x);

			double a11 = ut[x==t0];
			double a12 = vt[x==t0];
			double a21 = d_ut[x==t0];
			double a22 = d_vt[x==t0];
			
			vector<vector<double>> A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			A[0][0] = a11;
			A[0][1] = a12;
			A[1][0] = a21;
			A[1][1] = a22;
			vec_b[0][0] = y0;
			vec_b[1][0] = dy0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(A,vec_b,c_solution);
			//printVector(c_solution);

			Symbolic y_solution = exp((-b/(2*a))*x) * (c_solution[0]*(cos(D_real*x)) + c_solution[1]*(sin(D_real*x)));
			cout <<"\nThe real-valued initial value problem solution is:" << endl;
			cout << "\nu (t) = " << y_solution << endl;
			
			double mu = D_real;
			double Td = divisiond(2*pi,mu);
			double delta = atan(divisiond(c_solution[1],c_solution[0]));
			double t_eq = divisiond(1,mu)*(0.5*pi + delta);
			cout << "\nδ = " << delta << endl;
			cout << "\nT_{d} = " << Td << endl;
			cout << "\nThe time when the mass passes through its equilibrium position:\nt = " << t_eq << endl;

		}
	}
	 
}

void secondorderlineardiffeq_RLCserieselectriccircuit(double R, double L, double C, double y0, double dy0)
{
	Symbolic qt, q1, q2, c1("c1"), c2("c2"), q("q"), t("t");
 	double r1, r2;
	double t0 = 0;
	double a = 1;
	double b = divisiond(R,L);
	double c = divisiond(1,L*C);
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			qt = c1*exp(r1*t) + c2*t*exp(r2*t);
			cout <<"\nThe general solution is:" << endl;
			cout << qt << endl;

			double c1_ans = y0;
			double c2_ans = dy0 - (r1*c1_ans);
			cout <<"\nThe solution for the initial value problem / the charge Q at any time t is:" << endl;
			qt = qt[c1 == c1_ans, c2 == c2_ans] ;
			cout << "\nq(t) = " << qt << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			qt = c1*exp(r1*t) + c2*exp(r2*t);
			cout <<"\nThe general solution is:" << endl;
			cout << qt << endl;

			double c1_ans = divisiond(dy0-(r2*y0),r1-r2)*exp(-r1*t0);
			double c2_ans = divisiond((y0*r1)-dy0,r1-r2)*exp(-r2*t0);
			cout <<"\nThe solution for the initial value problem / the charge Q at any time t is:" << endl;
			qt = qt[c1 == c1_ans, c2 == c2_ans] ;
			cout << "\nq(t) = " << qt << endl;

			/*Symbolic dqt = df(qt,t);
			cout << "\nq'(t) = " << dqt << endl;

			double tm = NewtonRaphson(dqt,t,0);

			Equations rules = (  SymbolicConstant::e == exp(1), SymbolicConstant::i == sqrt(-1));
			qt = qt.subst_all(rules);
			double ym = qt[t==tm];

			cout << "\nCritical value t_{m} = " << tm << endl;
			cout << "\nMaximum value y_{m} = " << ym << endl;*/
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			
		
			q1 = exp((-b/(2*a))*t) *(cos(D_real*t) + SymbolicConstant::i*sin(D_real*t));
			q2 = exp((-b/(2*a))*t) *(cos(D_real*t) - SymbolicConstant::i*sin(D_real*t));

			Symbolic ut = exp((-b/(2*a))*t) *(cos(D_real*t));
			Symbolic vt = exp((-b/(2*a))*t) *(sin(D_real*t));
			Symbolic d_ut = df(ut,t);
			Symbolic d_vt = df(vt,t);

			double a11 = ut[t==t0];
			double a12 = vt[t==t0];
			double a21 = d_ut[t==t0];
			double a22 = d_vt[t==t0];
			
			vector<vector<double>> A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			A[0][0] = a11;
			A[0][1] = a12;
			A[1][0] = a21;
			A[1][1] = a22;
			vec_b[0][0] = y0;
			vec_b[1][0] = dy0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(A,vec_b,c_solution);
			//printVector(c_solution);

			Symbolic y_solution = exp((-b/(2*a))*t) * (c_solution[0]*(cos(D_real*t)) + c_solution[1]*(sin(D_real*t)));
			cout <<"\nThe solution for the initial value problem / the charge Q at any time t is:" << endl;
			cout << "\nq (t) = " << y_solution << endl;
			
		}
	}
	 
}

Symbolic wronskian_resultonly(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic yt, c1("c1"), c2("c2");
 	double r1, r2;
	Symbolic W;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = x*exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));	

			W = D_real*exp(2*(-b/(2*a))*x);		
		}
	}

	return W;
}

void wronskian(double a, double b, double c, const Symbolic &y, const Symbolic &x, double t0)
{
	Symbolic W, yt, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = x*exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));			
			W = D_real*exp(2*(-b/(2*a))*x);
		}
	}
	cout <<"\nThe general solution is:" << endl;
	cout << yt << endl;

	cout <<"\nThe Wronskian is: " << W << endl;

	cout << "\nW(t_{0}) = " << W[x==t0] << endl;
}

void wronskian_fundamentalsetofsolutions(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic yt, c1("c1"), c2("c2"), k1("k1"), k2("k2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));
		}
	}
	double t0 = 0;
	double y0 = 1;
	double dy0 = 0;
	double c1_ans = divisiond(dy0-(r2*y0),r1-r2)*exp(-r1*t0);
	double c2_ans = divisiond((y0*r1)-dy0,r1-r2)*exp(-r2*t0);
	cout <<"\nThe solution that satisfies the initial value problem y(0)=1 and y'(0)=0:" << endl;
	Symbolic yt_1 = yt[c1 == c1_ans, c2 == c2_ans] ;
	cout << yt_1 << endl;


	double y0_2= 0;
	double dy0_2 = 1;
	double c1_ans_2 = divisiond(dy0_2-(r2*y0_2),r1-r2)*exp(-r1*t0);
	double c2_ans_2 = divisiond((y0_2*r1)-dy0_2,r1-r2)*exp(-r2*t0);

	cout <<"\nThe solution that satisfies the initial value problem y(0)=0 and y'(0)=1:" << endl;
	Symbolic yt_2 = yt[c1 == c1_ans_2, c2 == c2_ans_2] ;
	cout << yt_2 << endl;

	cout <<"\nThe general solution :" << endl;
	Symbolic yt_general = k1*yt_1 + k2*yt_2;
	cout << yt_general << endl;

}

#include "polynomial.h"
#include "rational.h"

void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, Polynomial<double> &rhs_function)
{
	Symbolic A("A"), t("t");
	Symbolic lhs_final, Yt_final;

	Polynomial<double> p = rhs_function;
	cout << "p(t) = " << p << endl;
	int n = rhs_function.terms.front().second; // degree of the polynomial
	if (n >= 0 ) // this is working for Polynomial of degree n, n = positive integer
	{
		Symbolic Yt;
		for (int i = 0; i <= n; ++i)
		{
			Yt += A[i]*(t^(n-i)); // A[i] is symbolic A_{i} so we don't need to declare infinite number of Symbolic A1("A1"),... ,An("An")
		}
		Symbolic dy = df(Yt,t);
		Symbolic ddy = df(dy,t);
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;

		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		
		for (int i = 0; i <= n; ++i)
		{
			Symbolic coeff_ti = lhs_final.coeff(t,n-i); 
			for (int j = 0; j <= n; ++j)
			{
				mat_A[i][j] = coeff_ti.coeff(A[j],1);
			}
		}

		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		for (int i = 0 ;  i < n+1 ; ++i)
		{
			Yt_final += c_solution[i]*(t^(n-i));
		}
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}
	
	Symbolic ut, yt, y1, y2, c1s("c1"), c2s("c2");
	double a = lhs_a;
	double b = lhs_b;
	double c = lhs_c;
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1s*exp(r1*t) + c2s*t*exp(r2*t);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
			
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1s*exp(r1*t) + c2s*exp(r2*t);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
			
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			

			yt = exp((-b/(2*a))*t) * (c1s*(cos(D_real*t)) + c2s*(sin(D_real*t)));
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << "\ny(t) = " << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
		}
	}

	/*if (n == 3 )
	{
		Symbolic Yt = A*t*t*t + B*t*t + C*t + D;
		Symbolic dy = 3*A*t*t + 2*B*t +C;
		Symbolic ddy = 6*A*t + 2*B;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t3 =  lhs_final.coeff(t,3) ;
		Symbolic coeff_t2 =  lhs_final.coeff(t,2) ;
		Symbolic coeff_t1 =  lhs_final.coeff(t,1) ;
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t3.coeff(A,1);
		mat_A[0][1] = coeff_t3.coeff(B,1);
		mat_A[0][2] = coeff_t3.coeff(C,1);
		mat_A[0][3] = coeff_t3.coeff(D,1);
		mat_A[1][0] = coeff_t2.coeff(A,1);
		mat_A[1][1] = coeff_t2.coeff(B,1);
		mat_A[1][2] = coeff_t2.coeff(C,1);
		mat_A[1][3] = coeff_t2.coeff(D,1);
		mat_A[2][0] = coeff_t1.coeff(A,1);
		mat_A[2][1] = coeff_t1.coeff(B,1);
		mat_A[2][2] = coeff_t1.coeff(C,1);
		mat_A[2][3] = coeff_t1.coeff(D,1);
		mat_A[3][0] = coeff_t0.coeff(A,1);
		mat_A[3][1] = coeff_t0.coeff(B,1);
		mat_A[3][2] = coeff_t0.coeff(C,1);
		mat_A[3][3] = coeff_t0.coeff(D,1);

		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0]*(t^3) + c_solution[1]*(t^2)+ c_solution[2]*t + c_solution[3];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}
	if (n == 2 )
	{
		Symbolic Yt = A*t*t + B*t + C;
		Symbolic dy = 2*A*t + B;
		Symbolic ddy = 2*A;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t2 =  lhs_final.coeff(t,2) ;
		Symbolic coeff_t1 =  lhs_final.coeff(t,1) ;
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t2.coeff(A,1);
		mat_A[0][1] = coeff_t2.coeff(B,1);
		mat_A[0][2] = coeff_t2.coeff(C,1);
		mat_A[1][0] = coeff_t1.coeff(A,1);
		mat_A[1][1] = coeff_t1.coeff(B,1);
		mat_A[1][2] = coeff_t1.coeff(C,1);
		mat_A[2][0] = coeff_t0.coeff(A,1);
		mat_A[2][1] = coeff_t0.coeff(B,1);
		mat_A[2][2] = coeff_t0.coeff(C,1);

		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0]*(t^2) + c_solution[1]*t+ c_solution[2];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}
	if (n == 1 )
	{
		Symbolic Yt = A*t + B;
		Symbolic dy = A;
		Symbolic ddy = 0;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t1 =  lhs_final.coeff(t,1) ;
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t1.coeff(A,1);
		mat_A[0][1] = coeff_t1.coeff(B,1);
		mat_A[1][0] = coeff_t0.coeff(A,1);
		mat_A[1][1] = coeff_t0.coeff(B,1);
	
		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0]*(t) + c_solution[1];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}

	if (n == 0 )
	{
		Symbolic Yt = A;
		Symbolic dy = 0;
		Symbolic ddy = 0;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t0.coeff(A,1);
	
		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}*/
}

void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const Symbolic &rhs_function, const Symbolic &y, const Symbolic &x)
{
	Symbolic A("A"), B("B");
	
	double c_final, c1, c2, c3, lhs_final, rhs_final;
	Symbolic Yt_final;

	if(rhs_function != 0 )
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, f;
		// Case 1 : g(t) = a*exp(b*t)
		eq = (a*exp(b*x)).match(rhs_function, (a,b));
		
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
		
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);
		if(c_final != INFINITY) 
		{
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		if(c_final ==  INFINITY) 
		{
			Symbolic Yt = x*exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			Symbolic c1 = ddy.coeff(exp(bp*x),1);
			Symbolic c2 = dy.coeff(exp(bp*x),1);
			Symbolic c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c ;
			Symbolic subtract = lhs_final;
			lhs_final = lhs_final - subtract.coeff(x*exp(bp*x),1)*x*exp(bp*x) ;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
	
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : g(t) = exp(b*t)
		eq = (exp(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic bp = rhs(*i, b);
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;
		if(df(rhs(*i, b), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : g(t) = a*exp(t)
		eq = (a*exp(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = exp(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(x),1);
		c2 = dy.coeff(exp(x),1);
		c3 = Yt.coeff(exp(x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;

		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : g(t) = a*sin(t)
		eq = (a*sin(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 5 : g(t) = a*cos(t)
		eq = (a*cos(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 6 : g(t) = a*sin(b*t)
		eq = (a*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 7 : g(t) = a*cos(b*t)
		eq = (a*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 8 : g(t) = a*exp(t)*sin(b*t)
		eq = (a*exp(x)*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 9 : g(t) = a*exp(t)*cos(b*t)
		eq = (a*exp(x)*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 10 : g(t) = a*exp(c*t)*sin(b*t)
		eq = (a*exp(c*x)*sin(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 11 : g(t) = a*exp(c*t)*cos(b*t)
		eq = (a*exp(c*x)*cos(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp= rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 12 : g(t) = a*exp(c*t)*sin(t)
		eq = (a*exp(c*x)*sin(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 13 : g(t) = a*exp(c*t)*cos(t)
		eq = (a*exp(c*x)*cos(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}
	}
	Symbolic ut, yt, y1, y2, c1s("c1"), c2s("c2");
	double a = lhs_a;
	double b = lhs_b;
	double c = lhs_c;
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1s*exp(r1*x) + c2s*x*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
			
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1s*exp(r1*x) + c2s*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
			
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			

			yt = exp((-b/(2*a))*x) * (c1s*(cos(D_real*x)) + c2s*(sin(D_real*x)));
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << "\ny(t) = " << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
		}
	}
}

void secondorderlineardiffeq_nonhomogeneousequationsivpsolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const Symbolic &rhs_function,  double y0, double dy0, const Symbolic &y, const Symbolic &x)
{
	Symbolic A("A"), B("B");

	double c_final, c1, c2, c3, lhs_final, rhs_final;
	Symbolic Yt_final;

	if(rhs_function != 0 )
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, f;
		// Case 1 : g(t) = a*exp(b*t)
		eq = (a*exp(b*x)).match(rhs_function, (a,b));
		
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
		
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);
		if(c_final != INFINITY) 
		{
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		if(c_final ==  INFINITY) 
		{
			Symbolic Yt = x*exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			Symbolic c1 = ddy.coeff(exp(bp*x),1);
			Symbolic c2 = dy.coeff(exp(bp*x),1);
			Symbolic c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c ;
			Symbolic subtract = lhs_final;
			lhs_final = lhs_final - subtract.coeff(x*exp(bp*x),1)*x*exp(bp*x) ;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
	
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : g(t) = exp(b*t)
		eq = (exp(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic bp = rhs(*i, b);
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;
		if(df(rhs(*i, b), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : g(t) = a*exp(t)
		eq = (a*exp(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = exp(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(x),1);
		c2 = dy.coeff(exp(x),1);
		c3 = Yt.coeff(exp(x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;

		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : g(t) = a*sin(t)
		eq = (a*sin(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 5 : g(t) = a*cos(t)
		eq = (a*cos(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 6 : g(t) = a*sin(b*t)
		eq = (a*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 7 : g(t) = a*cos(b*t)
		eq = (a*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 8 : g(t) = a*exp(t)*sin(b*t)
		eq = (a*exp(x)*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 9 : g(t) = a*exp(t)*cos(b*t)
		eq = (a*exp(x)*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 10 : g(t) = a*exp(c*t)*sin(b*t)
		eq = (a*exp(c*x)*sin(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 11 : g(t) = a*exp(c*t)*cos(b*t)
		eq = (a*exp(c*x)*cos(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp= rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 12 : g(t) = a*exp(c*t)*sin(t)
		eq = (a*exp(c*x)*sin(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 13 : g(t) = a*exp(c*t)*cos(t)
		eq = (a*exp(c*x)*cos(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}
	}

	Symbolic ut, ut0, ut1, yt, y1, y2, c1s("c1"), c2s("c2");
	double a = lhs_a;
	double b = lhs_b;
	double c = lhs_c;
	double t0 = 0;
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1s*exp(r1*x) + c2s*x*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			ut0 = ut[x==t0];
			ut1 = df(ut,x);
			ut1 = ut1[x==t0];
			
			double c1_ans = solve(ut0-y0,c1s).front().rhs;
			ut1 = (ut1 - dy0);
			ut1 = ut1[c1s == c1_ans];
			double c2_ans = solve(ut1,c2s).front().rhs;
			cout <<"\nThe solution for the initial value problem is:" << endl;
			ut = ut[c1s == c1_ans, c2s == c2_ans] ;
			cout << ut << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1s*exp(r1*x) + c2s*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			ut0 = ut[x==t0];
			ut1 = df(ut,x);
			ut1 = ut1[x==t0];
			
			double c1_ans = solve(ut0-y0,c1s).front().rhs;
			ut1 = (ut1 - dy0);
			ut1 = ut1[c1s == c1_ans];
			double c2_ans = solve(ut1,c2s).front().rhs;
			cout <<"\nThe solution for the initial value problem is:" << endl;
			ut = ut[c1s == c1_ans, c2s == c2_ans] ;
			cout << ut << endl;

		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			

			yt = exp((-b/(2*a))*x) * (c1s*(cos(D_real*x)) + c2s*(sin(D_real*x)));
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << "\ny(t) = " << yt << endl;

			ut = yt + Yt_final;
			ut0 = ut[x==t0];
			ut1 = df(ut,x);
			ut1 = ut1[x==t0];
			
			double c1_ans = solve(ut0-y0,c1s).front().rhs;
			ut1 = (ut1 - dy0);
			ut1 = ut1[c1s == c1_ans];
			double c2_ans = solve(ut1,c2s).front().rhs;
			cout <<"\nThe solution for the initial value problem is:" << endl;
			ut = ut[c1s == c1_ans, c2s == c2_ans] ;
			cout << ut << endl;
			
		}
	}

}
void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const SymbolicMatrix &Matrix_A, const Symbolic &y, const Symbolic &x)
{ //Code it in 51 minutes, a silly bug on redeclaring Symbolic Yt_final occurs, it should be done in 30 minutes, on February 14th, 2026
	Symbolic A("A"), B("B");
	
	Symbolic Yt_solution, Yt_final;
 	int n_row = Matrix_A.rows();

	for (int i = 0; i < n_row ; ++i)
	{
		double c_final, c1, c2, c3, lhs_final, rhs_final;
		Symbolic rhs_function = Matrix_A[i][0];
		cout << "\n***********************************************************"<< endl;
		cout << "\nFor rhs = " << rhs_function << endl;
		if(rhs_function != 0)
	 	{
			list<Equations> eq;
			list<Equations>::iterator i;
			UniqueSymbol a, b, c, d, f;
			// Case 1 : g(t) = a*exp(b*t)
			eq = (a*exp(b*x)).match(rhs_function, (a,b));
			
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
			Symbolic Yt = exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			c1 = ddy.coeff(exp(bp*x),1);
			c2 = dy.coeff(exp(bp*x),1);
			c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
			Yt_final = c_final*Yt ;

			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n " << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 2 : g(t) = exp(b*t)
			eq = (exp(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic bp = rhs(*i, b);
			Symbolic Yt = exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			c1 = ddy.coeff(exp(bp*x),1);
			c2 = dy.coeff(exp(bp*x),1);
			c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
			Yt_final = c_final*Yt ;

			if(df(rhs(*i, b), x) == 0) 
			{
				cout << "\n " << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 3 : g(t) = a*exp(t)
			eq = (a*exp(x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a);
			Symbolic Yt = exp(x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			c1 = ddy.coeff(exp(x),1);
			c2 = dy.coeff(exp(x),1);
			c3 = Yt.coeff(exp(x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
			rhs_final = rhs_function.coeff(exp(x),1);
			c_final = divisiond(rhs_final,lhs_final);
			Yt_final = c_final*Yt ;

			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 4 : g(t) = a*sin(t)
			eq = (a*sin(x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a);
			Symbolic Yt = A*sin(x) + B*cos(x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(x),1);
			Symbolic coeff_cos = Ly.coeff(cos(x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(sin(x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 5 : g(t) = a*cos(t)
			eq = (a*cos(x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a);
			Symbolic Yt = A*sin(x) + B*cos(x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(x),1);
			Symbolic coeff_cos = Ly.coeff(cos(x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(cos(x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
			
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n " << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 6 : g(t) = a*sin(b*t)
			eq = (a*sin(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(sin(bp*x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
			cout << Yt_final << endl;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 7 : g(t) = a*cos(b*t)
			eq = (a*cos(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(cos(bp*x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 8 : g(t) = a*exp(t)*sin(b*t)
			eq = (a*exp(x)*sin(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
			Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
		 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(exp(x)*sin(bp*x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 9 : g(t) = a*exp(t)*cos(b*t)
			eq = (a*exp(x)*cos(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
			Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
		 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(exp(x)*cos(bp*x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			} 

			// Case 10 : g(t) = a*exp(c*t)*sin(b*t)
			eq = (a*exp(c*x)*sin(b*x)).match(rhs_function, (a,c,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
			Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(bp*x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 11 : g(t) = a*exp(c*t)*cos(b*t)
			eq = (a*exp(c*x)*cos(b*x)).match(rhs_function, (a,c,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
			Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(bp*x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 12 : g(t) = a*exp(c*t)*sin(t)
			eq = (a*exp(c*x)*sin(x)).match(rhs_function, (a,c));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
			Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 13 : g(t) = a*exp(c*t)*cos(t)
			eq = (a*exp(c*x)*cos(x)).match(rhs_function, (a,c));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
			Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}
		}
			cout << "\nThe current solution is\nY(t) = " << Yt_final ;
			Yt_solution += Yt_final;
		}
	cout << "\n***********************************************************"<< endl;
	cout << "\n***********************************************************"<< endl;

	cout << "\nThe particular solution is\nY(t) = " << Yt_solution ;

	Symbolic ut, yt, y1, y2, c1s("c1"), c2s("c2");
	double a = lhs_a;
	double b = lhs_b;
	double c = lhs_c;
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1s*exp(r1*x) + c2s*x*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
			
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1s*exp(r1*x) + c2s*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
			
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			

			yt = exp((-b/(2*a))*x) * (c1s*(cos(D_real*x)) + c2s*(sin(D_real*x)));
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << "\ny(t) = " << yt << endl;

			ut = yt + Yt_final;
			cout <<"\nThe general solution for the nonhomogeneous equation is:" << endl;
			cout << ut << endl;
		}
	}
}

void secondorderlineardiffeq_nonhomogeneousequationsivpforcedvibrationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const Symbolic &rhs_function,  double y0, double dy0, const Symbolic &y, const Symbolic &x)
{
	Symbolic A("A"), B("B");
	Symbolic Yt_final;

	double c_final, c1, c2, c3, lhs_final, rhs_final;
	double F0, omega;
	double m = lhs_a;
	double gamma = lhs_b;
	double k = lhs_c;
	
	if(rhs_function != 0 )
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, f;
		// Case 1 : g(t) = a*exp(b*t)
		eq = (a*exp(b*x)).match(rhs_function, (a,b));
		
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
		
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);
		if(c_final != INFINITY) 
		{
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		if(c_final ==  INFINITY) 
		{
			Symbolic Yt = x*exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			Symbolic c1 = ddy.coeff(exp(bp*x),1);
			Symbolic c2 = dy.coeff(exp(bp*x),1);
			Symbolic c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c ;
			Symbolic subtract = lhs_final;
			lhs_final = lhs_final - subtract.coeff(x*exp(bp*x),1)*x*exp(bp*x) ;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
	
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : g(t) = exp(b*t)
		eq = (exp(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic bp = rhs(*i, b);
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;
		if(df(rhs(*i, b), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : g(t) = a*exp(t)
		eq = (a*exp(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = exp(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(x),1);
		c2 = dy.coeff(exp(x),1);
		c3 = Yt.coeff(exp(x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;

		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : g(t) = a*sin(t)
		eq = (a*sin(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		F0 = ap;
		omega = 1;
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 5 : g(t) = a*cos(t)
		eq = (a*cos(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		F0 = ap;
		omega = 1;
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 6 : g(t) = a*sin(b*t)
		eq = (a*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		F0 = ap;
		omega = bp;
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 7 : g(t) = a*cos(b*t)
		eq = (a*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		F0 = ap;
		omega = bp;
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 8 : g(t) = a*exp(t)*sin(b*t)
		eq = (a*exp(x)*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 9 : g(t) = a*exp(t)*cos(b*t)
		eq = (a*exp(x)*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 10 : g(t) = a*exp(c*t)*sin(b*t)
		eq = (a*exp(c*x)*sin(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 11 : g(t) = a*exp(c*t)*cos(b*t)
		eq = (a*exp(c*x)*cos(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp= rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 12 : g(t) = a*exp(c*t)*sin(t)
		eq = (a*exp(c*x)*sin(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 13 : g(t) = a*exp(c*t)*cos(t)
		eq = (a*exp(c*x)*cos(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}
	}

	Symbolic ut, ut0, ut1, yt, y1, y2, c1s("c1"), c2s("c2");
	double a = lhs_a;
	double b = lhs_b;
	double c = lhs_c;
	double t0 = 0;
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1s*exp(r1*x) + c2s*x*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			ut0 = ut[x==t0];
			ut1 = df(ut,x);
			ut1 = ut1[x==t0];
			
			double c1_ans = solve(ut0-y0,c1s).front().rhs;
			ut1 = (ut1 - dy0);
			ut1 = ut1[c1s == c1_ans];
			double c2_ans = solve(ut1,c2s).front().rhs;
			cout <<"\nThe solution for the initial value problem is:" << endl;
			ut = ut[c1s == c1_ans, c2s == c2_ans] ;
			cout << ut << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1s*exp(r1*x) + c2s*exp(r2*x);
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << yt << endl;

			ut = yt + Yt_final;
			ut0 = ut[x==t0];
			ut1 = df(ut,x);
			ut1 = ut1[x==t0];
			
			double c1_ans = solve(ut0-y0,c1s).front().rhs;
			ut1 = (ut1 - dy0);
			ut1 = ut1[c1s == c1_ans];
			double c2_ans = solve(ut1,c2s).front().rhs;
			cout <<"\nThe solution for the initial value problem is:" << endl;
			ut = ut[c1s == c1_ans, c2s == c2_ans] ;
			cout << ut << endl;

		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			

			yt = exp((-b/(2*a))*x) * (c1s*(cos(D_real*x)) + c2s*(sin(D_real*x)));
			cout <<"\nThe general solution for the homogeneous equation is:" << endl;
			cout << "\ny(t) = " << yt << endl;

			ut = yt + Yt_final;
			ut0 = ut[x==t0];
			ut1 = df(ut,x);
			ut1 = ut1[x==t0];
			
			double c1_ans = solve(ut0-y0,c1s).front().rhs;
			ut1 = (ut1 - dy0);
			ut1 = ut1[c1s == c1_ans];
			double c2_ans = solve(ut1,c2s).front().rhs;
			cout <<"\nThe solution for the initial value problem is:" << endl;
			ut = ut[c1s == c1_ans, c2s == c2_ans] ;
			cout << ut << endl;
			
		}
	}
	double omega0=sqrt(divisiond(k,m));
	double triangle = sqrt(m*m*pow(omega0*omega0 - omega*omega,2) + gamma*gamma*omega*omega);
	double Gamma= divisiond(gamma*gamma,m*k);
	double R = divisiond(F0,triangle);
	double delta = asin(divisiond(gamma*omega,triangle));
	//double delta2 = acos(divisiond(m*(omega0*omega0 - omega*omega),triangle));
	double omega_max = sqrt(omega0*omega0 -divisiond(gamma*gamma,2*m*m));
	double R_max = divisiond(F0,gamma*omega0*sqrt(1-divisiond(gamma*gamma,4*m*k)));
	cout << "\nThe forced vibrations parameters: " << endl;
	cout << "F_{0} = "<< F0 << endl;
	cout << "ω = "<< omega << endl;
	cout << "ω_{0} = "<< omega0 << endl;
	cout << "ω_{max} = "<< omega_max << endl;
	cout << "Γ = "<< Gamma << endl;
	cout << "R = "<< R << endl;
	cout << "R_{max} = "<< R_max << endl;
	cout << "△ = "<< triangle << endl;
	cout << "𝛿 = "<< delta << endl;
	//cout << "𝛿 from acos= "<< delta2 << endl;

}



void secondorderlineardiffeq_nonhomogeneousequationssolution_variationofparameters(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const Symbolic &rhs_function, const Symbolic &y, const Symbolic &x)
{

	Symbolic Yt, yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	double a  = lhs_a, b = lhs_b, c = lhs_c;
	Symbolic W = wronskian_resultonly(lhs_a, lhs_b, lhs_c,y,x);

	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			y1 = c1*exp(r1*x);
			y2 = c2*x*exp(r2*x);
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			y1 = c1*exp(r1*x);
			y2 = c2*exp(r2*x);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));			
			y1 = exp((-b/(2*a))*x) * (c1*(cos(D_real*x))) ;
			y2 = exp((-b/(2*a))*x) * (c2*(sin(D_real*x))) ;
		}
	}	
	cout << "y_{1} (t) = " << y1/c1 << endl;
	cout << "y_{2} (t) = " << y2/c2 << endl;
	//cout << "W = " << W << endl;
	Symbolic Y1 = (y2/c2)*rhs_function/(W);
	Y1 = simplifybeforeintegrate(Y1,x);
	Y1 = -(y1/c1)*integrate(Y1,x);	
	Symbolic Y2 = (y1/c1)*rhs_function/(W);
	Y2 = simplifybeforeintegrate(Y2,x);
	Y2 = (y2/c2)*integrate(Y2,x);	
	//Symbolic du1 = -((y2/c2)*rhs_function)/W;	
	//Symbolic du2 = ((y1/c1)*rhs_function)/W;	
	//cout << "u_{1}' = " << du1 << endl;
	//cout << "u_{2}' = " << du2 << endl;


	Yt = Y1 + Y2;
	cout << "Y_{t} = " << Yt << endl;
	yt = Yt + y1 + y2;
	cout << "\nThe general solution is\ny_{t} = " << yt << endl;

}

void higherorderlineardiffeq_homogeneousequationsivpsolution(const vector<complex<double>> &P, const vector<complex<double>> &vec_ic, const vector<complex<double>> &vec_x0,  int N)
{
// Done on March 23rd, 2026
	int n_Polynomial = P.size();
	int n_ic = vec_ic.size();
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_imagroots;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);
	Symbolic general_solution, ivp_solution, df_solution;

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	if(n_ic != n_Polynomial-1)	
	{
		cerr << "Error: The number of initial conditions has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}
	cout << "P: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}

	cout << "\nInitial guess: " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		// To show the process of the Abert-Ehrlich
		//cout <<"\niteration: " << k << endl;
		//cout << "\nz_{i} new: " << endl;
		//printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
	// We use lround because there is an occurence if the root is obtained at very small decimal 
	// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

		if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2),roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
	// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
	// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
	vector<std::complex<double>> vec_unique;
	vector<std::complex<double>> vec_duplicate;

	int m = vec_update.size();
	for (int i = 0 ; i < m ; ++i)
	{
		double a  = real(vec_update[i]);
		double b  = imag(vec_update[i]);

		std::complex<double> target(a, b);

		// Get number of occurrences
		long count = std::count(vec_update.begin(), vec_update.end(), target);

		//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
	
		if(count == 1)
		{
			vec_unique.push_back(vec_update[i]);
		}
		else if(count > 1)
		{
			vec_duplicate.push_back(vec_update[i]);
		}
	}
	std::cout << "\nUnique vector:" << std::endl;
	printComplexVector(vec_unique);
	std::cout << "\nDuplicate vector:" << std::endl;
    	printComplexVector(vec_duplicate);

	// End of splitting into duplicate and unique vectors

	int n_unique = vec_unique.size();
	int n_duplicate = vec_duplicate.size();
	
	// This is for the case when the roots are unique, no duplicate / repeated roots.
	if(n_duplicate == 0)
	{
		Symbolic t("t"), c("c");
		for(int i = 0; i < n;++i)
		{
			complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
			vec_imagroots.push_back(root) ;
		}
		//cout << "\nThe abs imag roots" << endl ;
		//printComplexVector(vec_imagroots);

		// 1. Sort the vector using a custom comparator
		/*
		std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
		{
		if (a.real() != b.real()) 
		{
			return a.real() < b.real();
		}
		else
		{	
			return a.imag() < b.imag();
		}
		});	*/

		//    Use std::unique to move all non-duplicate elements to the front
		//    and return an iterator to the new logical end of the unique range.
		auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

		//    Erase the duplicate elements from the end of the vector.
		vec_imagroots.erase(last, vec_imagroots.end());

		//cout << "\nThe abs imag roots after delete duplicate" << endl ;
		//printComplexVector(vec_imagroots);
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			double lambda = real(vec_imagroots[i]);
			if(mu != 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
				i=i+1;
			}
			else if(mu == 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
				
			}
		}
		cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
		df_solution = general_solution;
		/*for(int i = 0; i < n;++i)
		{
			df_solution = df(df_solution,t);
			cout << "\ny^{(" << i << ")} = " << df_solution << endl;
			cout << "\ny^{(" << i << ")} (0)= " << df_solution[t==0] << endl;
		}*/

		vector<vector<double>> mat_A(n, vector<double>(n));
		vector<vector<double>> vec_b(n, vector<double>(1));
		for(int i = 0; i < n;++i)
		{
			Symbolic df_solution_ivp = df_solution[t==0] ;
			for(int j = 0; j < n;++j)
			{
				mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
			}	
			df_solution = df(df_solution,t);		
		}		
		for(int i = 0; i < n;++i)
		{
			vec_b[i][0] = real(vec_ic[i]);
		}

		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
		//printVector(c_solution);

		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			double lambda = real(vec_imagroots[i]);
			if(mu != 0)
			{
				ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
				i=i+1;
			}
			else if(mu == 0)
			{
				ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
				
			}
		}
		cout << "\nThe initial value problem solution is: \ny(t) = "<< ivp_solution << endl;
	}

	else if(n_duplicate != 0 )
	{
		
		Symbolic t("t"), c("c");
		Symbolic general_solution, ivp_solution, df_solution;
		if (n_unique != 0 && n_duplicate != n)
		{
			// We handle for the unique roots first
			for(int i = 0; i < n_unique;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			//cout << "\nThe abs imag unique roots after delete duplicate" << endl ;
			//printComplexVector(vec_imagroots);

			/*

				TEST 

			*/

			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			cout << "\nSorted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				//cout << "i = " << i << endl;
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "vec occurence = " << vec_occurence[i] << endl;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last2, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last3 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last3, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size() - n_unique;

			int index_i_continuing;
			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					
				}
				index_i_continuing = i+1;
			}

			int m_index = index_i_continuing;
			int i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
				}
				i_occurence = i_occurence+1;
			}
		
			/*

				END OF TEST 

			*/

			cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
			df_solution = general_solution;
		
			vector<vector<double>> mat_A(n, vector<double>(n));
			vector<vector<double>> vec_b(n, vector<double>(1));
			for(int i = 0; i < n;++i)
			{
				Symbolic df_solution_ivp = df_solution[t==0] ;
				for(int j = 0; j < n;++j)
				{
					mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
				}	
				df_solution = df(df_solution,t);		
			}		
			for(int i = 0; i < n;++i)
			{
				vec_b[i][0] = real(vec_ic[i]);
			}

			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
			//printMatrix(mat_A);
			//printVector(c_solution);

			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
				}
				index_i_continuing = i+1;
			}

			m_index = index_i_continuing;
			i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
				}
				i_occurence = i_occurence+1;
			}
		}
		else if (n_duplicate == n)
		{
			//cout << "\nn unique == 0 "<< endl;
			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
		
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			//cout << "\nSorted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last2, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size();

			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = 0; j <= n_occurence ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = 0; j < n_occurence + n_unique ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
				}
			}
			cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
			df_solution = general_solution;
		
			vector<vector<double>> mat_A(n, vector<double>(n));
			vector<vector<double>> vec_b(n, vector<double>(1));
			for(int i = 0; i < n;++i)
			{
				Symbolic df_solution_ivp = df_solution[t==0] ;
				for(int j = 0; j < n;++j)
				{
					mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
				}	
				df_solution = df(df_solution,t);		
			}		
			for(int i = 0; i < n;++i)
			{
				vec_b[i][0] = real(vec_ic[i]);
			}

			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
			//printMatrix(mat_A);
			//printVector(c_solution);

			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = abs(imag(vec_duplicate[i]));
				double lambda = real(vec_duplicate[i]);
				int n_occurence = vec_occurence[i];
				if(mu != 0)
				{			
					int k = 0;		
					for(int j = 0; j <= n_occurence ; j++)
					{
						ivp_solution += exp(lambda*t)*c_solution[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j+1] *sin(mu*t) * pow(t,Symbolic(k)) ;
						k= k+1;
						j = j+1;
					}
					i=i+n_occurence;
				}
				else if(mu == 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
					
				}
			}
		}
		cout << "\nThe initial value problem solution is: \ny(t) = "<< ivp_solution << endl;
	}
}

void higherorderlineardiffeq_homogeneousequationsivpsolution(const vector<complex<double>> &P, const vector<complex<double>> &vec_ic)
{
	vector<complex<double>> vec_x0;
	int n_vec = P.size() - 1;
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
	for(int i=0; i<n_vec; i++) // we create vector x0 with the size of the highest order of the differential equation
	{
		double real_part = distribution(generator);
		double imag_part = 0;
		complex<double> random_complex(real_part, imag_part);
		vec_x0.push_back(random_complex); 	
	}

	int N = 100;
	int n_Polynomial = P.size();
	int n_ic = vec_ic.size();
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_imagroots;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	if(n_ic != n_Polynomial-1)	
	{
		cerr << "Error: The number of initial conditions has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}

	cout << "\nP: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}

	cout << "\nInitial conditions : " << endl;
	printComplexVector(vec_ic);
	cout << "\nInitial guess for the roots (generated randomly): " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		// To show the process of the Abert-Ehrlich
		//cout <<"\niteration: " << k << endl;
		//cout << "\nz_{i} new: " << endl;
		//printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
	// We use lround because there is an occurence if the root is obtained at very small decimal 
	// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

		if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2),roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
	// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
	// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
	vector<std::complex<double>> vec_unique;
	vector<std::complex<double>> vec_duplicate;

	int m = vec_update.size();
	for (int i = 0 ; i < m ; ++i)
	{
		double a  = real(vec_update[i]);
		double b  = imag(vec_update[i]);

		std::complex<double> target(a, b);

		// Get number of occurrences
		long count = std::count(vec_update.begin(), vec_update.end(), target);

		//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
	
		if(count == 1)
		{
			vec_unique.push_back(vec_update[i]);
		}
		else if(count > 1)
		{
			vec_duplicate.push_back(vec_update[i]);
		}
	}
	//cout << "\nUnique vector:" << std::endl;
	//printComplexVector(vec_unique);
	//cout << "\nDuplicate vector:" << std::endl;
    	//printComplexVector(vec_duplicate);

	// End of splitting into duplicate and unique vectors

	int n_unique = vec_unique.size();
	int n_duplicate = vec_duplicate.size();
	
	// This is for the case when the roots are unique, no duplicate / repeated roots.
	if(n_duplicate == 0)
	{
		Symbolic t("t"), c("c");
		Symbolic general_solution, ivp_solution, df_solution;
		for(int i = 0; i < n;++i)
		{
			complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
			vec_imagroots.push_back(root) ;
		}
		//cout << "\nThe abs imag roots" << endl ;
		//printComplexVector(vec_imagroots);

		// 1. Sort the vector using a custom comparator
		
		std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
		{
		if (a.real() != b.real()) 
		{
			return a.real() < b.real();
		}
		else
		{	
			return a.imag() < b.imag();
		}
		});	
		setprecision(5);
		//    Use std::unique to move all non-duplicate elements to the front
		//    and return an iterator to the new logical end of the unique range.
		auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

		//    Erase the duplicate elements from the end of the vector.
		vec_imagroots.erase(last, vec_imagroots.end());

		// FInd a way to delete duplicate root / the complex conjugate.
		//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
		//printComplexVector(vec_imagroots);
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			//cout << "mu = " << mu << endl;
			double lambda = real(vec_imagroots[i]);
			if(mu != 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
				i=i+1;
			}
			else if(mu == 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
				
			}
		}
		cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
		df_solution = general_solution;
		/*for(int i = 0; i < n;++i)
		{
			df_solution = df(df_solution,t);
			cout << "\ny^{(" << i << ")} = " << df_solution << endl;
			cout << "\ny^{(" << i << ")} (0)= " << df_solution[t==0] << endl;
		}*/

		vector<vector<double>> mat_A(n, vector<double>(n));
		vector<vector<double>> vec_b(n, vector<double>(1));
		for(int i = 0; i < n;++i)
		{
			Symbolic df_solution_ivp = df_solution[t==0] ;
			for(int j = 0; j < n;++j)
			{
				mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
			}	
			df_solution = df(df_solution,t);		
		}		
		for(int i = 0; i < n;++i)
		{
			vec_b[i][0] = real(vec_ic[i]);
		}

		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
		//printVector(c_solution);

		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			double lambda = real(vec_imagroots[i]);
			if(mu != 0)
			{
				ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
				i=i+1;
			}
			else if(mu == 0)
			{
				ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
				
			}
		}
		cout << "\nThe initial value problem solution is: \ny(t) = "<< ivp_solution << endl;
	}
	else if(n_duplicate != 0 )
	{
		Symbolic t("t"), c("c");
		Symbolic general_solution, ivp_solution, df_solution;
		if (n_unique != 0 && n_duplicate != n)
		{
			// We handle for the unique roots first
			for(int i = 0; i < n_unique;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			//cout << "\nThe abs imag unique roots after delete duplicate" << endl ;
			//printComplexVector(vec_imagroots);

			/*

				TEST 

			*/

			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			cout << "\nSorted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				//cout << "i = " << i << endl;
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "vec occurence = " << vec_occurence[i] << endl;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last2, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last3 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last3, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size() - n_unique;

			int index_i_continuing;
			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					
				}
				index_i_continuing = i+1;
			}
			int m_index = index_i_continuing;
			int i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
				}
				i_occurence = i_occurence+1;
			}
		
			/*

				END OF TEST 

			*/

			cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
			df_solution = general_solution;
		
			vector<vector<double>> mat_A(n, vector<double>(n));
			vector<vector<double>> vec_b(n, vector<double>(1));
			for(int i = 0; i < n;++i)
			{
				Symbolic df_solution_ivp = df_solution[t==0] ;
				for(int j = 0; j < n;++j)
				{
					mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
				}	
				df_solution = df(df_solution,t);		
			}		
			for(int i = 0; i < n;++i)
			{
				vec_b[i][0] = real(vec_ic[i]);
			}

			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
			//printMatrix(mat_A);
			//printVector(c_solution);

			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
				}
				index_i_continuing = i+1;
			}

			m_index = index_i_continuing;
			i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
				}
				i_occurence = i_occurence+1;
			}
			
		}
		else if (n_duplicate == n)
		{
			//cout << "\nn unique == 0 "<< endl;
			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
		
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			//cout << "\nSorted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last2, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size();

			int m_index2 = 0;
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
				}
			}
			cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
			df_solution = general_solution;
		
			vector<vector<double>> mat_A(n, vector<double>(n));
			vector<vector<double>> vec_b(n, vector<double>(1));
			for(int i = 0; i < n;++i)
			{
				Symbolic df_solution_ivp = df_solution[t==0] ;
				for(int j = 0; j < n;++j)
				{
					mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
				}	
				df_solution = df(df_solution,t);		
			}		
			for(int i = 0; i < n;++i)
			{
				vec_b[i][0] = real(vec_ic[i]);
			}

			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
			//printMatrix(mat_A);
			//printVector(c_solution);

			m_index2 = 0;
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
				}
			}

		}
		cout << "\nThe initial value problem solution is: \ny(t) = "<< ivp_solution << endl;
	}
}

vector<complex<double>> higherorderlineardiffeq_vectorize(const Symbolic &diffeq, const Symbolic &y, const Symbolic &t, int n)
{
	cout << "\nThe ODE with constant coefficients:" << endl;
	cout << diffeq << " = 0 " << endl;
	
	vector<double> vec_coeff;
	Symbolic dummy;
	for(int i=1; i<=n+1 ; ++i)
	{
		Symbolic var_coeff = df(y[t],t,i);
		Symbolic var_coeff0 = diffeq.coeff(var_coeff,0);
		//cout << "\ni = " << i << endl;		
		//cout << "var = " << var_coeff << endl;
		//cout << "var_coeff0 = " << var_coeff0 << endl;
		if(i==1)
		{
			vec_coeff.push_back(var_coeff0.coeff(df(y[t],t,0),1));
		}
		
		if(i>1)
		{
			dummy = diffeq.coeff(df(y[t],t,i-1),0);
			//cout << "dummy = " << dummy << endl;
			if(var_coeff0 == dummy)
			{
				vec_coeff.push_back(0);
			}
			
			if(var_coeff0 != dummy)
			{
				Symbolic new_term = var_coeff0 - dummy;
				//cout << "new term = " << new_term << endl;
				double cnt = new_term/df(y[t],t,i-1);
				vec_coeff.push_back(cnt);
			}
		}
	}
	reverse(vec_coeff.begin(),vec_coeff.end());	// reverse the vector because we need the highest order at the first entry

	vector<complex<double>> vec_complex;
	
	// Populate complex vector P from vec_coeff
	for (int i = 0; i < n+1; ++i)
	{
		if(vec_coeff[0] != 1.0)
		{
			complex<double> P_entry(divisiond(vec_coeff[i],vec_coeff[0]),0.0);
			vec_complex.push_back(P_entry);
		}
		else 
		{
			complex<double> P_entry(vec_coeff[i],0.0);
			vec_complex.push_back(P_entry);
		}
	}
	return vec_complex;
}


void higherorderlineardiffeq_twospringtwomasssystem(double k1, double k2, double m1, double m2, Symbolic &t, const vector<complex<double>> &vec_ic)
{
	// Created on April 4th, 2026. This code contains how to get the order of derivative and its coefficient.
	// Modified on April 8th, 2026. Adding m_index and m_index2.
	Symbolic u1("u1"), u2("u2");
	Symbolic  eq_motion1, eq_motion2,  u1_solution ;

	eq_motion1 = m1*df(df(u1[t],t),t) - k2*(u2-u1) + k1*u1;
	eq_motion2 = m2*df(df(u2[t],t),t) + k2*(u2-u1);
	cout << "\nThe equation of motion : " << endl;
	cout << eq_motion1 << " = 0 " << endl;
	cout << eq_motion2 << " = 0 " << endl;

	Equations F_u2 = solve(eq_motion1,u2);
	cout << "\n" << F_u2.front() << endl;
	
	Symbolic u2_solve = F_u2.front().rhs;
	u2_solve = u2_solve[u1==u1[t]];

	// The equation of motion for mass 2 in u1 terms, all are function of t
	Symbolic eq_motion2_u1terms = m2*df(df(u2_solve[t],t),t) + k2*(u2_solve[t]-u1[t]);
 	cout << "\nThe homogeneous ODE with constant coefficients:" << endl;
	cout << eq_motion2_u1terms<< " = 0 " << endl;
	vector<double> vec_coeff, vec_coeff_u2;
	Symbolic dummy;
	for(int i=1; i<=5 ; ++i)
	{
		Symbolic var_coeff = df(u1[t],t,i);
		Symbolic var_coeff0 = eq_motion2_u1terms.coeff(var_coeff,0);
		//cout << "\ni = " << i << endl;		
		//cout << "var = " << var_coeff << endl;
		//cout << "var_coeff0 = " << var_coeff0 << endl;
		if(i==1)
		{
			vec_coeff.push_back(var_coeff0.coeff(df(u1[t],t,0),1));
		}
		
		if(i>1)
		{
			dummy = eq_motion2_u1terms.coeff(df(u1[t],t,i-1),0);
			//cout << "dummy = " << dummy << endl;
			if(var_coeff0 == dummy)
			{
				vec_coeff.push_back(0);
			}
			
			if(var_coeff0 != dummy)
			{
				Symbolic new_term = var_coeff0 - dummy;
				//cout << "new term = " << new_term << endl;
				double cnt = new_term/df(u1[t],t,i-1);
				vec_coeff.push_back(cnt);
			}
		}
	}
	reverse(vec_coeff.begin(),vec_coeff.end());	// reverse the vector because we need the highest order at the first entry
	
	int n_vec = vec_coeff.size();
	vector<complex<double>> P;
	vector<complex<double>> vec_x0;
	vector<complex<double>> vec_ic_u1;
	vector<complex<double>> vec_ic_u2;
	int N = 100;
	// Populate complex vector P from vec_coeff
	for (int i = 0; i < n_vec; ++i)
	{
		if(vec_coeff[0] != 1.0)
		{
			complex<double> P_entry(divisiond(vec_coeff[i],vec_coeff[0]),0.0);
			P.push_back(P_entry);
		}
		else 
		{
			complex<double> P_entry(vec_coeff[i],0.0);
			P.push_back(P_entry);
		}
	}
	
	/*
		START OF IVP SOLUTION FOR HIGHER ORDER LINEAR DIFFERENTIAL EQUATION
	*/	

	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
	for(int i=0; i<n_vec-1; i++)
	{
		double real_part = distribution(generator);
		double imag_part = 0;
		complex<double> random_complex(real_part, imag_part);
		vec_x0.push_back(random_complex); 	
	}
	
	int n_Polynomial = P.size();
	int n_ic = vec_ic.size();
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_imagroots;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	if(n_ic != n_Polynomial-1)	
	{
		cerr << "Error: The number of initial conditions has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}
	cout << "\n*******************************************************************************" << endl;
	cout << "\n********************              For u_{1}(t)               ******************" << endl;
	cout << "\n*******************************************************************************" << endl;

	cout << "\nP: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}
	for (int i = 0; i < 2; ++i)
	{
		vec_ic_u1.push_back(vec_ic[i]);
	}
	for (int i = 2; i < 4; ++i)
	{
		double u1_ic_new = divisiond(k2*(real(vec_ic[i])-real(vec_ic[i-2])) - k1*real(vec_ic[i-2]),m1);
		complex<double> complex_u1_ic_new(u1_ic_new, 0.0);	
		vec_ic_u1.push_back(complex_u1_ic_new);
	}
	
	cout << "\nInitial conditions (u_{1} and u_{2}): " << endl;
	printComplexVector(vec_ic);
	cout << "\nInitial conditions (u_{1}): " << endl;
	printComplexVector(vec_ic_u1);
	cout << "\nInitial guess for the roots (generated randomly): " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		// To show the process of the Abert-Ehrlich
		//cout <<"\niteration: " << k << endl;
		//cout << "\nz_{i} new: " << endl;
		//printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
	// We use lround because there is an occurence if the root is obtained at very small decimal 
	// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

		if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),6), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),6));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),6),roundToDecimal(imag(vec_dummy[i]),6));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
	// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
	// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
	vector<std::complex<double>> vec_unique;
	vector<std::complex<double>> vec_duplicate;

	int m = vec_update.size();
	for (int i = 0 ; i < m ; ++i)
	{
		double a  = real(vec_update[i]);
		double b  = imag(vec_update[i]);

		std::complex<double> target(a, b);

		// Get number of occurrences
		long count = std::count(vec_update.begin(), vec_update.end(), target);

		//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
	
		if(count == 1)
		{
			vec_unique.push_back(vec_update[i]);
		}
		else if(count > 1)
		{
			vec_duplicate.push_back(vec_update[i]);
		}
	}
	//cout << "\nUnique vector:" << std::endl;
	//printComplexVector(vec_unique);
	//cout << "\nDuplicate vector:" << std::endl;
    	//printComplexVector(vec_duplicate);

	// End of splitting into duplicate and unique vectors

	int n_unique = vec_unique.size();
	int n_duplicate = vec_duplicate.size();
	
	// This is for the case when the roots are unique, no duplicate / repeated roots.
	if(n_duplicate == 0)
	{
		Symbolic t("t"), c("c");
		Symbolic general_solution, ivp_solution, df_solution;
		for(int i = 0; i < n;++i)
		{
			complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
			vec_imagroots.push_back(root) ;
		}
		//cout << "\nThe abs imag roots" << endl ;
		//printComplexVector(vec_imagroots);

		// 1. Sort the vector using a custom comparator
		
		std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
		{
		if (a.real() != b.real()) 
		{
			return a.real() < b.real();
		}
		else
		{	
			return a.imag() < b.imag();
		}
		});	
		setprecision(5);
		//    Use std::unique to move all non-duplicate elements to the front
		//    and return an iterator to the new logical end of the unique range.
		auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

		//    Erase the duplicate elements from the end of the vector.
		vec_imagroots.erase(last, vec_imagroots.end());

		// FInd a way to delete duplicate root / the complex conjugate.
		//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
		//printComplexVector(vec_imagroots);
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			//cout << "mu = " << mu << endl;
			double lambda = real(vec_imagroots[i]);
			if(mu != 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
				i=i+1;
			}
			else if(mu == 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
				
			}
		}
		cout << "\nThe general solution is: \nu_{1}(t) = "<< general_solution << endl;
		df_solution = general_solution;
		/*for(int i = 0; i < n;++i)
		{
			df_solution = df(df_solution,t);
			cout << "\ny^{(" << i << ")} = " << df_solution << endl;
			cout << "\ny^{(" << i << ")} (0)= " << df_solution[t==0] << endl;
		}*/

		vector<vector<double>> mat_A(n, vector<double>(n));
		vector<vector<double>> vec_b(n, vector<double>(1));
		for(int i = 0; i < n;++i)
		{
			Symbolic df_solution_ivp = df_solution[t==0] ;
			for(int j = 0; j < n;++j)
			{
				mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
			}	
			df_solution = df(df_solution,t);		
		}		
		for(int i = 0; i < n;++i)
		{
			vec_b[i][0] = real(vec_ic_u1[i]);
		}

		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
		//printVector(c_solution);
		
		for(int i = 0; i < n;++i)
		{
			if (abs(c_solution[i]) > 1e-5 )
			{
				c_solution[i] = c_solution[i];
			}
			if (abs(c_solution[i]) < 1e-5 )
			{
				c_solution[i] = 0;
			}

		}
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			double lambda = real(vec_imagroots[i]);
			if(mu != 0)
			{
				ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
				i=i+1;
			}
			else if(mu == 0)
			{
				ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
				
			}
		}
		cout << "\nThe initial value problem solution is: \nu_{1}(t) = "<< ivp_solution << endl;
		u1_solution = ivp_solution;
		//cout << "u1 solution = "<< u1_solution << endl;
		// The equation of motion for mass 2 in u1 terms, all are function of t
		Symbolic eq_motion2_withu1solution = m2*df(df(u2[t],t),t) + k2*(u2[t]-u1_solution);
	 	Symbolic eq_motion2_withoutu1solution = m2*df(df(u2[t],t),t) + k2*(u2[t]);
	 	
		cout << "\nThe equation of motion 2 with u_{1}(t):" << endl;
		cout << eq_motion2_withu1solution << " = 0 " << endl;

		// to get the coefficients from equation of motion 2, for the associated homogeneous equation with u_{1}(t)=0
		for(int i=1; i<=3 ; ++i)
		{
			Symbolic var_coeff = df(u2[t],t,i);
			Symbolic var_coeff0 = eq_motion2_withoutu1solution.coeff(var_coeff,0);
			//cout << "\ni = " << i << endl;		
			//cout << "var = " << var_coeff << endl;
			//cout << "var_coeff0 = " << var_coeff0 << endl;
			if(i==1)
			{
				vec_coeff_u2.push_back(var_coeff0.coeff(df(u2[t],t,0),1));
			}
			
			if(i>1)
			{
				dummy = eq_motion2_withoutu1solution.coeff(df(u2[t],t,i-1),0);
				//cout << "dummy = " << dummy << endl;
				if(var_coeff0 == dummy)
				{
					vec_coeff_u2.push_back(0);
				}
				
				if(var_coeff0 != dummy)
				{
					Symbolic new_term = var_coeff0 - dummy;
					//cout << "new term = " << new_term << endl;
					double cnt = new_term/df(u2[t],t,i-1);
					vec_coeff_u2.push_back(cnt);
				}
			}
		}
		
		reverse(vec_coeff_u2.begin(),vec_coeff_u2.end());	// reverse the vector because we need the highest order at the first entry
	
		int n_vec_u2 = vec_coeff_u2.size();
		P.clear();
		vec_x0.clear();
		// Populate complex vector P from vec_coeff
		for (int i = 0; i < n_vec_u2; ++i)
		{
			if(vec_coeff_u2[0] != 1.0)
			{
				complex<double> P_entry(divisiond(vec_coeff_u2[i],vec_coeff_u2[0]),0.0);
				P.push_back(P_entry);
			}
			else 
			{
				complex<double> P_entry(vec_coeff_u2[i],0.0);
				P.push_back(P_entry);
			}
		}
		cout << "\n*******************************************************************************" << endl;
		cout << "\n********************              For u_{2}(t)               ******************" << endl;
		cout << "\n*******************************************************************************" << endl;
		
		/*
		START OF IVP SOLUTION FOR HIGHER ORDER LINEAR DIFFERENTIAL EQUATION for u_{2}(t)
		*/	

		// 1. Obtain a seed:
		// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
		// provides a more robust seed than a fixed value.
		std::default_random_engine generator(
		std::chrono::system_clock::now().time_since_epoch().count());
		
		vec.clear();
	 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
		for(int i=0; i<n_vec_u2-1; i++)
		{
			double real_part = distribution(generator);
			double imag_part = 0;
			complex<double> random_complex(real_part, imag_part);
			vec_x0.push_back(random_complex); 	
		}

		for (int i = 0; i < 2; ++i)
		{
			vec_ic_u2.push_back(vec_ic[i+2]);
		}

		n_Polynomial = P.size();
		n_ic = vec_ic_u2.size();
		complex<double> nP2(n_Polynomial,0.0);
		n = vec_x0.size();
		P_derivative.clear();
		vec_update.clear();
		vec_imagroots.clear();
		vec_dummy.clear();
		vec_check.clear();
		complex<double> i_derivative2(1.0,0.0);
		for (int i = 0; i < n; ++i)
		{
			vec_dummy.push_back(vec_x0[i]);
		}
		
		if(n != n_Polynomial-1)	
		{
			cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
		}
		if(n_ic != n_Polynomial-1)	
		{
			cerr << "Error: The number of initial conditions has to be: the number of highest order of the derivative." << endl;
		}
		for (int i = 0; i < n_Polynomial - 1; ++i)
		{
			P_derivative.push_back((nP2 - i_derivative2)*P[i]);
			i_derivative2 = i_derivative2 + c1;
		}
		cout << "\nP: " << endl;
		printComplexVector(P);
		cout << "\nP': " << endl;
		printComplexVector(P_derivative);

		cout << "\nInitial conditions (u_{2}): " << endl;
		printComplexVector(vec_ic_u2);
		cout << "\nInitial guess for the roots (generated randomly): " << endl;
		printComplexVector(vec_dummy);

		for (int k = 0; k < N ; ++k)
		{
			vec_check.clear();
			for (int i = 0; i < n; ++i)
			{
				vec_check.push_back(vec_dummy[i]);
			}

			for (int i = 0; i < n; ++i)
			{
				//cout <<"\ni: " << i << endl;

				// Newton step
				complex<double> P_zi(0.0,0.0);
				complex<double> P_zi_derivative(0.0,0.0);
				complex<double> i_der(1.0,0.0);
				complex<double> i_der2(2.0,0.0);
				for (int j = 0; j < n_Polynomial ; ++j)
				{
					P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

					i_der = i_der + c1;
				}		
				//cout <<" P(z_{i}) : " << P_zi << endl;
				for (int j = 0; j < n_Polynomial - 1 ; ++j)
				{
					P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

					i_der2 = i_der2 + c1;
				}
				//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

				complex<double> N_zi =P_zi/P_zi_derivative;
				//cout <<" N(z_{i}) : " << N_zi << endl;

				// Shift
				complex<double> shift(0.0,0.0);
				for (int j = 0; j < n ; ++j)
				{
					if(i != j)
					{
						shift += c1/(vec_dummy[i]-vec_dummy[j]);
					}
					else if(i == j)
					{
						shift += 0;
					}
				}
				//cout <<" S(z_{i}) : " << shift << endl;

				vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
			}

			// To show the process of the Abert-Ehrlich
			//cout <<"\niteration: " << k << endl;
			//cout << "\nz_{i} new: " << endl;
			//printComplexVector(vec_dummy);

			complex<double> epsilon(0.0,0.0);		
			for (int i = 0; i < n; ++i)
			{
				epsilon += vec_dummy[i]-vec_check[i];
			}
			double diffnorm = moduluscomplex(epsilon);
			if( diffnorm < 1e-12)
			{
				cout << "\nRoots found at iteration: " << k << endl;
				k = N-1;		
			}
			
		}

		for(int i = 0; i < n;++i)
		{
		// We use lround because there is an occurence if the root is obtained at very small decimal 
		// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

			if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),6), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),6));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),6),roundToDecimal(imag(vec_dummy[i]),6));
			vec_update.push_back(root);
		}

		}

		cout << "\n************************************************************************" << endl;
		cout << "\nEnd of iteration" << endl;
		cout << "\nz_{i} final: " << endl;
		printComplexVector(vec_update);
			
		// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
		// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
		vector<std::complex<double>> vec_unique;
		vector<std::complex<double>> vec_duplicate;

		int m = vec_update.size();
		for (int i = 0 ; i < m ; ++i)
		{
			double a  = real(vec_update[i]);
			double b  = imag(vec_update[i]);

			std::complex<double> target(a, b);

			// Get number of occurrences
			long count = std::count(vec_update.begin(), vec_update.end(), target);

			//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
		
			if(count == 1)
			{
				vec_unique.push_back(vec_update[i]);
			}
			else if(count > 1)
			{
				vec_duplicate.push_back(vec_update[i]);
			}
		}
		//cout << "\nUnique vector:" << std::endl;
		//printComplexVector(vec_unique);
		//cout << "\nDuplicate vector:" << std::endl;
	    	//printComplexVector(vec_duplicate);

		// End of splitting into duplicate and unique vectors

		//int n_unique = vec_unique.size();
		int n_duplicate = vec_duplicate.size();
		
		// This is for the case when the roots are unique, no duplicate / repeated roots.
		if(n_duplicate == 0)
		{
			Symbolic t("t"), c("c");
			Symbolic general_solution, ivp_solution, df_solution;
			for(int i = 0; i < n;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}
			//cout << "\nThe abs imag roots" << endl ;
			//printComplexVector(vec_imagroots);

			// 1. Sort the vector using a custom comparator
			
			std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			setprecision(5);
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			// FInd a way to delete duplicate root / the complex conjugate.
			//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
			//printComplexVector(vec_imagroots);
			for(int i = 0; i < n;++i)
			{
				double mu = imag(vec_imagroots[i]);
				//cout << "mu = " << mu << endl;
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					
				}
			}
			Symbolic A("A");
			Symbolic A_solve = solve(m2*A*df(u1_solution,t,2)+(k2*A*u1_solution)-(k2*u1_solution),A).front().rhs;
			//cout << A_solve<< endl;
			A_solve = roundToDecimal(A_solve,2);
			Symbolic particular_solution_u2 = A_solve*u1_solution;
			
			general_solution = general_solution + particular_solution_u2;
			cout << "\nThe general solution is: \nu_{2}(t) = "<< general_solution << endl;
			df_solution = general_solution;

			vector<double> c_solution;

			for(int i = 0; i < n;++i)
			{
				//cout << "\ny^{(" << i << ")} = " << df_solution << endl;
				//cout << "\ny^{(" << i << ")} (0)= " << df_solution[t==0] << endl;
				Symbolic dfs_0 = df_solution[t==0];
				
				double cs = solve(dfs_0-real(vec_ic_u2[i]),c[i]).front().rhs;

				c_solution.push_back(cs);

				df_solution = df(df_solution,t);
			}

			for(int i = 0; i < n;++i)
			{
				if (abs(c_solution[i]) > 1e-5 )
				{
					c_solution[i] = c_solution[i];
				}
				if (abs(c_solution[i]) < 1e-5 )
				{
					c_solution[i] = 0;
				}

			}
			//printVector(c_solution);

			for(int i = 0; i < n;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
					
				}
			}
			cout << "\nThe initial value problem solution is: \nu_{2}(t) = "<< ivp_solution + particular_solution_u2 << endl;
		}
	}

	else if(n_duplicate != 0 )
	{
		
		Symbolic t("t"), c("c");
		Symbolic general_solution, ivp_solution, df_solution;
		if (n_unique != 0 && n_duplicate != n)
		{
			// We handle for the unique roots first
			for(int i = 0; i < n_unique;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			//cout << "\nThe abs imag unique roots after delete duplicate" << endl ;
			//printComplexVector(vec_imagroots);

			/*

				TEST 

			*/

			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			cout << "\nSorted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				//cout << "i = " << i << endl;
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "vec occurence = " << vec_occurence[i] << endl;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last2, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last3 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last3, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size() - n_unique;

			int index_i_continuing;
			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					
				}
				index_i_continuing = i+1;
			}
			int m_index = index_i_continuing;
			int i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
				}
				i_occurence = i_occurence+1;
			}
		
			/*

				END OF TEST 

			*/

			cout << "\nThe general solution is: \nu_{1}(t) = "<< general_solution << endl;
			df_solution = general_solution;
		
			vector<vector<double>> mat_A(n, vector<double>(n));
			vector<vector<double>> vec_b(n, vector<double>(1));
			for(int i = 0; i < n;++i)
			{
				Symbolic df_solution_ivp = df_solution[t==0] ;
				for(int j = 0; j < n;++j)
				{
					mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
				}	
				df_solution = df(df_solution,t);		
			}		
			for(int i = 0; i < n;++i)
			{
				vec_b[i][0] = real(vec_ic_u1[i]);
			}

			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
			//printMatrix(mat_A);
			//printVector(c_solution);

			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
				}
				index_i_continuing = i+1;
			}

			m_index = index_i_continuing;
			i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						ivp_solution += exp(lambda*t)*c_solution[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
				}
				i_occurence = i_occurence+1;
			}
			cout << "\nThe initial value problem solution is: \nu_{1}(t) = "<< ivp_solution << endl;
			u1_solution = ivp_solution;
			//cout << "u1 solution = "<< u1_solution << endl;
			// The equation of motion for mass 2 in u1 terms, all are function of t
			Symbolic eq_motion2_withu1solution = m2*df(df(u2[t],t),t) + k2*(u2[t]-u1_solution);
		 	Symbolic eq_motion2_withoutu1solution = m2*df(df(u2[t],t),t) + k2*(u2[t]);
	 	
			cout << "\nThe equation of motion 2 with u_{1}(t):" << endl;
			cout << eq_motion2_withu1solution << " = 0 " << endl;

			// to get the coefficients from equation of motion 2, for the associated homogeneous equation with u_{1}(t)=0
			for(int i=1; i<=3 ; ++i)
			{
				Symbolic var_coeff = df(u2[t],t,i);
				Symbolic var_coeff0 = eq_motion2_withoutu1solution.coeff(var_coeff,0);
				//cout << "\ni = " << i << endl;		
				//cout << "var = " << var_coeff << endl;
				//cout << "var_coeff0 = " << var_coeff0 << endl;
				if(i==1)
				{
					vec_coeff_u2.push_back(var_coeff0.coeff(df(u2[t],t,0),1));
				}
				
				if(i>1)
				{
					dummy = eq_motion2_withoutu1solution.coeff(df(u2[t],t,i-1),0);
					//cout << "dummy = " << dummy << endl;
					if(var_coeff0 == dummy)
					{
						vec_coeff_u2.push_back(0);
					}
					
					if(var_coeff0 != dummy)
					{
						Symbolic new_term = var_coeff0 - dummy;
						//cout << "new term = " << new_term << endl;
						double cnt = new_term/df(u2[t],t,i-1);
						vec_coeff_u2.push_back(cnt);
					}
				}
			}
			
			reverse(vec_coeff_u2.begin(),vec_coeff_u2.end());	// reverse the vector because we need the highest order at the first entry
		
			int n_vec_u2 = vec_coeff_u2.size();
			P.clear();
			vec_x0.clear();
			// Populate complex vector P from vec_coeff
			for (int i = 0; i < n_vec_u2; ++i)
			{
				if(vec_coeff_u2[0] != 1.0)
				{
					complex<double> P_entry(divisiond(vec_coeff_u2[i],vec_coeff_u2[0]),0.0);
					P.push_back(P_entry);
				}
				else 
				{
					complex<double> P_entry(vec_coeff_u2[i],0.0);
					P.push_back(P_entry);
				}
			}
			cout << "\n*******************************************************************************" << endl;
			cout << "\n********************              For u_{2}(t)               ******************" << endl;
			cout << "\n*******************************************************************************" << endl;
			
			/*
			START OF IVP SOLUTION FOR HIGHER ORDER LINEAR DIFFERENTIAL EQUATION for u_{2}(t)
			*/	

			// 1. Obtain a seed:
			// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
			// provides a more robust seed than a fixed value.
			std::default_random_engine generator(
			std::chrono::system_clock::now().time_since_epoch().count());
			
			vec.clear();
		 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
			for(int i=0; i<n_vec_u2-1; i++)
			{
				double real_part = distribution(generator);
				double imag_part = 0;
				complex<double> random_complex(real_part, imag_part);
				vec_x0.push_back(random_complex); 	
			}

			for (int i = 0; i < 2; ++i)
			{
				vec_ic_u2.push_back(vec_ic[i+2]);
			}

			n_Polynomial = P.size();
			n_ic = vec_ic_u2.size();
			complex<double> nP2(n_Polynomial,0.0);
			n = vec_x0.size();
			P_derivative.clear();
			vec_update.clear();
			vec_imagroots.clear();
			vec_dummy.clear();
			vec_check.clear();
			complex<double> i_derivative2(1.0,0.0);
			for (int i = 0; i < n; ++i)
			{
				vec_dummy.push_back(vec_x0[i]);
			}
			
			if(n != n_Polynomial-1)	
			{
				cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
			}
			if(n_ic != n_Polynomial-1)	
			{
				cerr << "Error: The number of initial conditions has to be: the number of highest order of the derivative." << endl;
			}
			for (int i = 0; i < n_Polynomial - 1; ++i)
			{
				P_derivative.push_back((nP2 - i_derivative2)*P[i]);
				i_derivative2 = i_derivative2 + c1;
			}
			cout << "\nP: " << endl;
			printComplexVector(P);
			cout << "\nP': " << endl;
			printComplexVector(P_derivative);

			cout << "\nInitial conditions (u_{2}): " << endl;
			printComplexVector(vec_ic_u2);
			cout << "\nInitial guess for the roots (generated randomly): " << endl;
			printComplexVector(vec_dummy);

			for (int k = 0; k < N ; ++k)
			{
				vec_check.clear();
				for (int i = 0; i < n; ++i)
				{
					vec_check.push_back(vec_dummy[i]);
				}

				for (int i = 0; i < n; ++i)
				{
					//cout <<"\ni: " << i << endl;

					// Newton step
					complex<double> P_zi(0.0,0.0);
					complex<double> P_zi_derivative(0.0,0.0);
					complex<double> i_der(1.0,0.0);
					complex<double> i_der2(2.0,0.0);
					for (int j = 0; j < n_Polynomial ; ++j)
					{
						P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

						i_der = i_der + c1;
					}		
					//cout <<" P(z_{i}) : " << P_zi << endl;
					for (int j = 0; j < n_Polynomial - 1 ; ++j)
					{
						P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

						i_der2 = i_der2 + c1;
					}
					//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

					complex<double> N_zi =P_zi/P_zi_derivative;
					//cout <<" N(z_{i}) : " << N_zi << endl;

					// Shift
					complex<double> shift(0.0,0.0);
					for (int j = 0; j < n ; ++j)
					{
						if(i != j)
						{
							shift += c1/(vec_dummy[i]-vec_dummy[j]);
						}
						else if(i == j)
						{
							shift += 0;
						}
					}
					//cout <<" S(z_{i}) : " << shift << endl;

					vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
				}

				// To show the process of the Abert-Ehrlich
				//cout <<"\niteration: " << k << endl;
				//cout << "\nz_{i} new: " << endl;
				//printComplexVector(vec_dummy);

				complex<double> epsilon(0.0,0.0);		
				for (int i = 0; i < n; ++i)
				{
					epsilon += vec_dummy[i]-vec_check[i];
				}
				double diffnorm = moduluscomplex(epsilon);
				if( diffnorm < 1e-12)
				{
					cout << "\nRoots found at iteration: " << k << endl;
					k = N-1;		
				}
				
			}

			for(int i = 0; i < n;++i)
			{
			// We use lround because there is an occurence if the root is obtained at very small decimal 
			// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

				if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
				{
					complex<double> root(roundToDecimal(real(vec_dummy[i]),6), 0.0);
					vec_update.push_back(root);
				}
				if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
				{
					complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),6));
					vec_update.push_back(root);
				}
				if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
				{
					complex<double> root(roundToDecimal(real(vec_dummy[i]),6),roundToDecimal(imag(vec_dummy[i]),6));
					vec_update.push_back(root);
				}

			}

			cout << "\n************************************************************************" << endl;
			cout << "\nEnd of iteration" << endl;
			cout << "\nz_{i} final: " << endl;
			printComplexVector(vec_update);
				
			// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
			// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
			vector<std::complex<double>> vec_unique;
			vector<std::complex<double>> vec_duplicate;

			int m = vec_update.size();
			for (int i = 0 ; i < m ; ++i)
			{
				double a  = real(vec_update[i]);
				double b  = imag(vec_update[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				long count = std::count(vec_update.begin(), vec_update.end(), target);

				//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
			
				if(count == 1)
				{
					vec_unique.push_back(vec_update[i]);
				}
				else if(count > 1)
				{
					vec_duplicate.push_back(vec_update[i]);
				}
			}
			//cout << "\nUnique vector:" << std::endl;
			//printComplexVector(vec_unique);
			//cout << "\nDuplicate vector:" << std::endl;
		    	//printComplexVector(vec_duplicate);

			// End of splitting into duplicate and unique vectors

			//int n_unique = vec_unique.size();
			int n_duplicate = vec_duplicate.size();
			
			// This is for the case when the roots are unique, no duplicate / repeated roots.
			if(n_duplicate == 0)
			{
				Symbolic t("t"), c("c");
				Symbolic general_solution, ivp_solution, df_solution;
				for(int i = 0; i < n;++i)
				{
					complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
					vec_imagroots.push_back(root) ;
				}
				//cout << "\nThe abs imag roots" << endl ;
				//printComplexVector(vec_imagroots);

				// 1. Sort the vector using a custom comparator
				
				std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
				{
				if (a.real() != b.real()) 
				{
					return a.real() < b.real();
				}
				else
				{	
					return a.imag() < b.imag();
				}
				});	
				setprecision(5);
				//    Use std::unique to move all non-duplicate elements to the front
				//    and return an iterator to the new logical end of the unique range.
				auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

				//    Erase the duplicate elements from the end of the vector.
				vec_imagroots.erase(last, vec_imagroots.end());

				// FInd a way to delete duplicate root / the complex conjugate.
				//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
				//printComplexVector(vec_imagroots);
				for(int i = 0; i < n;++i)
				{
					double mu = imag(vec_imagroots[i]);
					//cout << "mu = " << mu << endl;
					double lambda = real(vec_imagroots[i]);
					if(mu != 0)
					{
						general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
						i=i+1;
					}
					else if(mu == 0)
					{
						general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
						
					}
				}
				Symbolic A("A");
				Symbolic A_solve = solve(m2*A*df(u1_solution,t,2)+(k2*A*u1_solution)-(k2*u1_solution),A).front().rhs;
				//cout << A_solve<< endl;
				A_solve = roundToDecimal(A_solve,2);
				Symbolic particular_solution_u2 = A_solve*u1_solution;
				
				general_solution = general_solution + particular_solution_u2;
				cout << "\nThe general solution is: \nu_{2}(t) = "<< general_solution << endl;
				df_solution = general_solution;

				vector<double> c_solution;

				for(int i = 0; i < n;++i)
				{
					//cout << "\ny^{(" << i << ")} = " << df_solution << endl;
					//cout << "\ny^{(" << i << ")} (0)= " << df_solution[t==0] << endl;
					Symbolic dfs_0 = df_solution[t==0];
					
					double cs = solve(dfs_0-real(vec_ic_u2[i]),c[i]).front().rhs;

					c_solution.push_back(cs);

					df_solution = df(df_solution,t);
				}

				for(int i = 0; i < n;++i)
				{
					if (abs(c_solution[i]) > 1e-5 )
					{
						c_solution[i] = c_solution[i];
					}
					if (abs(c_solution[i]) < 1e-5 )
					{
						c_solution[i] = 0;
					}

				}
				//printVector(c_solution);

				for(int i = 0; i < n;++i)
				{
					double mu = imag(vec_imagroots[i]);
					double lambda = real(vec_imagroots[i]);
					if(mu != 0)
					{
						ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
						i=i+1;
					}
					else if(mu == 0)
					{
						ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
						
					}
				}
				cout << "\nThe initial value problem solution is: \nu_{2}(t) = "<< ivp_solution + particular_solution_u2 << endl;
			}
		}
		else if (n_duplicate == n)
		{
			//cout << "\nn unique == 0 "<< endl;
			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
		
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			//cout << "\nSorted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last2, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size();

			int m_index2 = 0;
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
				}
			}
			cout << "\nThe general solution is: \nu_{1}(t) = "<< general_solution << endl;
			df_solution = general_solution;
		
			vector<vector<double>> mat_A(n, vector<double>(n));
			vector<vector<double>> vec_b(n, vector<double>(1));
			for(int i = 0; i < n;++i)
			{
				Symbolic df_solution_ivp = df_solution[t==0] ;
				for(int j = 0; j < n;++j)
				{
					mat_A[i][j] = df_solution_ivp.coeff(c[j],1);
				}	
				df_solution = df(df_solution,t);		
			}		
			for(int i = 0; i < n;++i)
			{
				vec_b[i][0] = real(vec_ic_u1[i]);
			}

			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);
			//printMatrix(mat_A);
			//printVector(c_solution);

			m_index2 = 0;
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c_solution[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c_solution[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c_solution[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
				}
			}

		}
		cout << "\nThe initial value problem solution is: \nu_{1}(t) = "<< ivp_solution << endl;
		u1_solution = ivp_solution;
		//cout << "u1 solution = "<< u1_solution << endl;
		// The equation of motion for mass 2 in u1 terms, all are function of t
		Symbolic eq_motion2_withu1solution = m2*df(df(u2[t],t),t) + k2*(u2[t]-u1_solution);
	 	Symbolic eq_motion2_withoutu1solution = m2*df(df(u2[t],t),t) + k2*(u2[t]);
	 	
		cout << "\nThe equation of motion 2 with u_{1}(t):" << endl;
		cout << eq_motion2_withu1solution << " = 0 " << endl;

			// to get the coefficients from equation of motion 2, for the associated homogeneous equation with u_{1}(t)=0
			for(int i=1; i<=3 ; ++i)
			{
				Symbolic var_coeff = df(u2[t],t,i);
				Symbolic var_coeff0 = eq_motion2_withoutu1solution.coeff(var_coeff,0);
				//cout << "\ni = " << i << endl;		
				//cout << "var = " << var_coeff << endl;
				//cout << "var_coeff0 = " << var_coeff0 << endl;
				if(i==1)
				{
					vec_coeff_u2.push_back(var_coeff0.coeff(df(u2[t],t,0),1));
				}
				
				if(i>1)
				{
					dummy = eq_motion2_withoutu1solution.coeff(df(u2[t],t,i-1),0);
					//cout << "dummy = " << dummy << endl;
					if(var_coeff0 == dummy)
					{
						vec_coeff_u2.push_back(0);
					}
					
					if(var_coeff0 != dummy)
					{
						Symbolic new_term = var_coeff0 - dummy;
						//cout << "new term = " << new_term << endl;
						double cnt = new_term/df(u2[t],t,i-1);
						vec_coeff_u2.push_back(cnt);
					}
				}
			}
			
			reverse(vec_coeff_u2.begin(),vec_coeff_u2.end());	// reverse the vector because we need the highest order at the first entry
		
			int n_vec_u2 = vec_coeff_u2.size();
			P.clear();
			vec_x0.clear();
			// Populate complex vector P from vec_coeff
			for (int i = 0; i < n_vec_u2; ++i)
			{
				if(vec_coeff_u2[0] != 1.0)
				{
					complex<double> P_entry(divisiond(vec_coeff_u2[i],vec_coeff_u2[0]),0.0);
					P.push_back(P_entry);
				}
				else 
				{
					complex<double> P_entry(vec_coeff_u2[i],0.0);
					P.push_back(P_entry);
				}
			}
			cout << "\n*******************************************************************************" << endl;
			cout << "\n********************              For u_{2}(t)               ******************" << endl;
			cout << "\n*******************************************************************************" << endl;
			
			/*
			START OF IVP SOLUTION FOR HIGHER ORDER LINEAR DIFFERENTIAL EQUATION for u_{2}(t)
			*/	

			// 1. Obtain a seed:
			// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
			// provides a more robust seed than a fixed value.
			std::default_random_engine generator(
			std::chrono::system_clock::now().time_since_epoch().count());
			
			vec.clear();
		 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
			for(int i=0; i<n_vec_u2-1; i++)
			{
				double real_part = distribution(generator);
				double imag_part = 0;
				complex<double> random_complex(real_part, imag_part);
				vec_x0.push_back(random_complex); 	
			}

			for (int i = 0; i < 2; ++i)
			{
				vec_ic_u2.push_back(vec_ic[i+2]);
			}

			n_Polynomial = P.size();
			n_ic = vec_ic_u2.size();
			complex<double> nP2(n_Polynomial,0.0);
			n = vec_x0.size();
			P_derivative.clear();
			vec_update.clear();
			vec_imagroots.clear();
			vec_dummy.clear();
			vec_check.clear();
			complex<double> i_derivative2(1.0,0.0);
			for (int i = 0; i < n; ++i)
			{
				vec_dummy.push_back(vec_x0[i]);
			}
			
			if(n != n_Polynomial-1)	
			{
				cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
			}
			if(n_ic != n_Polynomial-1)	
			{
				cerr << "Error: The number of initial conditions has to be: the number of highest order of the derivative." << endl;
			}
			for (int i = 0; i < n_Polynomial - 1; ++i)
			{
				P_derivative.push_back((nP2 - i_derivative2)*P[i]);
				i_derivative2 = i_derivative2 + c1;
			}
			cout << "\nP: " << endl;
			printComplexVector(P);
			cout << "\nP': " << endl;
			printComplexVector(P_derivative);

			cout << "\nInitial conditions (u_{2}): " << endl;
			printComplexVector(vec_ic_u2);
			cout << "\nInitial guess for the roots (generated randomly): " << endl;
			printComplexVector(vec_dummy);

			for (int k = 0; k < N ; ++k)
			{
				vec_check.clear();
				for (int i = 0; i < n; ++i)
				{
					vec_check.push_back(vec_dummy[i]);
				}

				for (int i = 0; i < n; ++i)
				{
					//cout <<"\ni: " << i << endl;

					// Newton step
					complex<double> P_zi(0.0,0.0);
					complex<double> P_zi_derivative(0.0,0.0);
					complex<double> i_der(1.0,0.0);
					complex<double> i_der2(2.0,0.0);
					for (int j = 0; j < n_Polynomial ; ++j)
					{
						P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

						i_der = i_der + c1;
					}		
					//cout <<" P(z_{i}) : " << P_zi << endl;
					for (int j = 0; j < n_Polynomial - 1 ; ++j)
					{
						P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

						i_der2 = i_der2 + c1;
					}
					//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

					complex<double> N_zi =P_zi/P_zi_derivative;
					//cout <<" N(z_{i}) : " << N_zi << endl;

					// Shift
					complex<double> shift(0.0,0.0);
					for (int j = 0; j < n ; ++j)
					{
						if(i != j)
						{
							shift += c1/(vec_dummy[i]-vec_dummy[j]);
						}
						else if(i == j)
						{
							shift += 0;
						}
					}
					//cout <<" S(z_{i}) : " << shift << endl;

					vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
				}

				// To show the process of the Abert-Ehrlich
				//cout <<"\niteration: " << k << endl;
				//cout << "\nz_{i} new: " << endl;
				//printComplexVector(vec_dummy);

				complex<double> epsilon(0.0,0.0);		
				for (int i = 0; i < n; ++i)
				{
					epsilon += vec_dummy[i]-vec_check[i];
				}
				double diffnorm = moduluscomplex(epsilon);
				if( diffnorm < 1e-12)
				{
					cout << "\nRoots found at iteration: " << k << endl;
					k = N-1;		
				}
				
			}

			for(int i = 0; i < n;++i)
			{
			// We use lround because there is an occurence if the root is obtained at very small decimal 
			// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

				if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
				{
					complex<double> root(roundToDecimal(real(vec_dummy[i]),6), 0.0);
					vec_update.push_back(root);
				}
				if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
				{
					complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),6));
					vec_update.push_back(root);
				}
				if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
				{
					complex<double> root(roundToDecimal(real(vec_dummy[i]),6),roundToDecimal(imag(vec_dummy[i]),6));
					vec_update.push_back(root);
				}

			}

			cout << "\n************************************************************************" << endl;
			cout << "\nEnd of iteration" << endl;
			cout << "\nz_{i} final: " << endl;
			printComplexVector(vec_update);
				
			// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
			// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
			vector<std::complex<double>> vec_unique;
			vector<std::complex<double>> vec_duplicate;

			int m = vec_update.size();
			for (int i = 0 ; i < m ; ++i)
			{
				double a  = real(vec_update[i]);
				double b  = imag(vec_update[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				long count = std::count(vec_update.begin(), vec_update.end(), target);

				//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
			
				if(count == 1)
				{
					vec_unique.push_back(vec_update[i]);
				}
				else if(count > 1)
				{
					vec_duplicate.push_back(vec_update[i]);
				}
			}
			//cout << "\nUnique vector:" << std::endl;
			//printComplexVector(vec_unique);
			//cout << "\nDuplicate vector:" << std::endl;
		    	//printComplexVector(vec_duplicate);

			// End of splitting into duplicate and unique vectors

			//int n_unique = vec_unique.size();
			int n_duplicate = vec_duplicate.size();
			
			// This is for the case when the roots are unique, no duplicate / repeated roots.
			if(n_duplicate == 0)
			{
				Symbolic t("t"), c("c");
				Symbolic general_solution, ivp_solution, df_solution;
				for(int i = 0; i < n;++i)
				{
					complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
					vec_imagroots.push_back(root) ;
				}
				//cout << "\nThe abs imag roots" << endl ;
				//printComplexVector(vec_imagroots);

				// 1. Sort the vector using a custom comparator
				
				std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
				{
				if (a.real() != b.real()) 
				{
					return a.real() < b.real();
				}
				else
				{	
					return a.imag() < b.imag();
				}
				});	
				setprecision(5);
				//    Use std::unique to move all non-duplicate elements to the front
				//    and return an iterator to the new logical end of the unique range.
				auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

				//    Erase the duplicate elements from the end of the vector.
				vec_imagroots.erase(last, vec_imagroots.end());

				// FInd a way to delete duplicate root / the complex conjugate.
				//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
				//printComplexVector(vec_imagroots);
				for(int i = 0; i < n;++i)
				{
					double mu = imag(vec_imagroots[i]);
					//cout << "mu = " << mu << endl;
					double lambda = real(vec_imagroots[i]);
					if(mu != 0)
					{
						general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
						i=i+1;
					}
					else if(mu == 0)
					{
						general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
						
					}
				}
				Symbolic A("A");
				Symbolic A_solve = solve(m2*A*df(u1_solution,t,2)+(k2*A*u1_solution)-(k2*u1_solution),A).front().rhs;
				//cout << A_solve<< endl;
				A_solve = roundToDecimal(A_solve,2);
				Symbolic particular_solution_u2 = A_solve*u1_solution;
				
				general_solution = general_solution + particular_solution_u2;
				cout << "\nThe general solution is: \nu_{2}(t) = "<< general_solution << endl;
				df_solution = general_solution;

				vector<double> c_solution;

				for(int i = 0; i < n;++i)
				{
					//cout << "\ny^{(" << i << ")} = " << df_solution << endl;
					//cout << "\ny^{(" << i << ")} (0)= " << df_solution[t==0] << endl;
					Symbolic dfs_0 = df_solution[t==0];
					
					double cs = solve(dfs_0-real(vec_ic_u2[i]),c[i]).front().rhs;

					c_solution.push_back(cs);

					df_solution = df(df_solution,t);
				}

				for(int i = 0; i < n;++i)
				{
					if (abs(c_solution[i]) > 1e-5 )
					{
						c_solution[i] = c_solution[i];
					}
					if (abs(c_solution[i]) < 1e-5 )
					{
						c_solution[i] = 0;
					}

				}
				//printVector(c_solution);

				for(int i = 0; i < n;++i)
				{
					double mu = imag(vec_imagroots[i]);
					double lambda = real(vec_imagroots[i]);
					if(mu != 0)
					{
						ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i+1] *sin(mu*t) ;
						i=i+1;
					}
					else if(mu == 0)
					{
						ivp_solution += exp(lambda*t)*c_solution[i] *cos(mu*t) + exp(lambda*t)*c_solution[i] *sin(mu*t) ;
						
					}
				}
				cout << "\nThe initial value problem solution is: \nu_{2}(t) = "<< ivp_solution + particular_solution_u2 << endl;
			}

	}
	/*
		END OF IVP SOLUTION FOR HIGHER ORDER LINEAR DIFFERENTIAL EQUATION
	*/	
}

void higherorderlineardiffeq_nonhomogeneousequationsgeneralsolution(const vector<complex<double>> &P, const Symbolic &rhs_function, const Symbolic &t)
{
	Symbolic A("A"), Yt_particular, Yt_particular_total;
	Symbolic Yt_final, Yt_assumed;
	Symbolic general_solution;
		
	vector<complex<double>> vec_x0;
	int n_vec = P.size() ;
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
	for(int i=0; i<n_vec-1; i++)
	{
		double real_part = distribution(generator);
		double imag_part = 0;
		complex<double> random_complex(real_part, imag_part);
		vec_x0.push_back(random_complex); 	
	}

	int N = 100;
	int n_Polynomial = P.size();
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	int m_index3;
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_imagroots;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}

	cout << "\nP: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}

	cout << "\nInitial guess for the roots (generated randomly): " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		// To show the process of the Abert-Ehrlich
		//cout <<"\niteration: " << k << endl;
		//cout << "\nz_{i} new: " << endl;
		//printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
	// We use lround because there is an occurence if the root is obtained at very small decimal 
	// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

		if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2),roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
	// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
	// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
	vector<std::complex<double>> vec_unique;
	vector<std::complex<double>> vec_duplicate;

	int m = vec_update.size();
	for (int i = 0 ; i < m ; ++i)
	{
		double a  = real(vec_update[i]);
		double b  = imag(vec_update[i]);

		std::complex<double> target(a, b);

		// Get number of occurrences
		long count = std::count(vec_update.begin(), vec_update.end(), target);

		//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
	
		if(count == 1)
		{
			vec_unique.push_back(vec_update[i]);
		}
		else if(count > 1)
		{
			vec_duplicate.push_back(vec_update[i]);
		}
	}
	//cout << "\nUnique vector:" << std::endl;
	//printComplexVector(vec_unique);
	//cout << "\nDuplicate vector:" << std::endl;
    	//printComplexVector(vec_duplicate);

	// End of splitting into duplicate and unique vectors

	int n_unique = vec_unique.size();
	int n_duplicate = vec_duplicate.size();
	
	// This is for the case when the roots are unique, no duplicate / repeated roots.
	if(n_duplicate == 0)
	{
		Symbolic c("c");
		for(int i = 0; i < n;++i)
		{
			complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
			vec_imagroots.push_back(root) ;
		}
		//cout << "\nThe abs imag roots" << endl ;
		//printComplexVector(vec_imagroots);

		// 1. Sort the vector using a custom comparator
		
		std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
		{
		if (a.real() != b.real()) 
		{
			return a.real() < b.real();
		}
		else
		{	
			return a.imag() < b.imag();
		}
		});	
		setprecision(5);
		//    Use std::unique to move all non-duplicate elements to the front
		//    and return an iterator to the new logical end of the unique range.
		auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

		//    Erase the duplicate elements from the end of the vector.
		vec_imagroots.erase(last, vec_imagroots.end());

		// FInd a way to delete duplicate root / the complex conjugate.
		//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
		//printComplexVector(vec_imagroots);
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			//cout << "mu = " << mu << endl;
			double lambda = real(vec_imagroots[i]);
			int k = 1;
			if(mu != 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
				i=i+1;

			}
			else if(mu == 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
			}
				Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
				
		}
		cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;
		Yt_particular = A*Yt_assumed;
		m_index3 = n_Polynomial-1;
		for(int i = 0; i <= n_Polynomial-1; ++i)
		{
			//cout <<"P[i] = " << real(P[i]) << endl;
			//cout <<"mindex3 = " << m_index3 << endl;
			//cout << "Yt^{(n)} = " << df(Yt_particular,t,m_index3)*real(P[i]) << endl;
			Yt_particular_total += df(Yt_particular,t,m_index3)*real(P[i]) ;
			m_index3 = m_index3-1;
			//cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
		}
	}
	else if(n_duplicate != 0 )
	{
		Symbolic c("c");
		if (n_unique != 0 && n_duplicate != n)
		{
			// We handle for the unique roots first
			for(int i = 0; i < n_unique;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			//cout << "\nThe abs imag unique roots after delete duplicate" << endl ;
			//printComplexVector(vec_imagroots);

			/*

				TEST 

			*/

			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			cout << "\nSorted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				//cout << "i = " << i << endl;
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "vec occurence = " << vec_occurence[i] << endl;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last2, vec_duplicate.end());
			//cout << "\nDeleted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last3 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last3, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size() - n_unique;
			
			int index_i_continuing;
			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					
				}
				index_i_continuing = i+1;
			}
			int m_index = index_i_continuing;
			int i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				//cout <<" i = "<< i << endl;
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				//cout << "n occurence = " <<  n_occurence << endl;
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
					
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j <  m_index+n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
					
				}
				i_occurence = i_occurence + 1;
			}
			cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;

			Yt_particular = A*Yt_assumed;
			m_index3 = n_Polynomial-1;

			for(int i = 0; i <= n_Polynomial-1; ++i)
			{
				//cout <<"P[i] = " << real(P[i]) << endl;
				//cout <<"mindex3 = " << m_index3 << endl;
				//cout << "Yt^{(n)} = " << df(Yt_particular,t,m_index3)*real(P[i]) << endl;
				Yt_particular_total += df(Yt_particular,t,m_index3)*real(P[i]) ;
				m_index3 = m_index3-1;
				//cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
			}

			/*

				END OF TEST 

			*/			
		}
		else if (n_duplicate == n)
		{
			//cout << "\nn unique == 0 "<< endl;
			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
		
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			//cout << "\nSorted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last2, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size();

			int m_index2 = 0;
			
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				//cout << "n_occurence = " << n_occurence  << endl;
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
				}
				
			}
			cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;

			m_index3 = n_Polynomial-1;
			Yt_particular = A*Yt_assumed;
			for(int i = 0; i <= n_Polynomial-1; ++i)
			{
				//cout <<"P[i] = " << real(P[i]) << endl;
				//cout <<"mindex3 = " << m_index3 << endl;
				//cout << "Yt^{(n)} = " << df(Yt_particular,t,m_index3)*real(P[i]) << endl;
				Yt_particular_total += df(Yt_particular,t,m_index3)*real(P[i]) ;
				m_index3 = m_index3-1;
				//cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
			}
						
		}	
	}
		Symbolic A_solve;
		Symbolic F = rhs_function;
		cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
		cout << "\nY(t) assumed = "<< Yt_assumed << endl;

		// Case 1 : g(t) = a*exp(t)
		if (df(F,t,0) == df(F,t,1) &&  df(F,t,1) == df(F,t,2) && df(F,t,0) == df(F,t,2)) // meaning that the rhs_function is in the form of a*exp(t)
		{
			//double a = df(F,t,0)/exp(t);
			Symbolic lhs_nh = Yt_particular_total.coeff(exp(t),1);
			A_solve = solve(lhs_nh*exp(t)- rhs_function, A).front().rhs;
			Yt_final = A_solve*Yt_assumed ;		
		}
		// Case 2 : g(t) = a*sin(t) or g(t) = a*cos(t)
		if (df(F,t,0) == -df(F,t,2) &&  df(F,t,0) != df(F,t,1) && df(F,t,0) == df(F,t,4)) // meaning that the rhs_function is in the form of a*sin(t) or a* cos(t)
		{
			if (df(F/sin(t),t) == 0)  // meaning that the rhs_function is in the form of a*sin(t)
			{
				double a = F/sin(t);
				Symbolic lhs_nh = Yt_particular_total.coeff(A*sin(t),1);
				A_solve = a/lhs_nh;
				Yt_final = A_solve*Yt_assumed.coeff(sin(t),1)*sin(t) ;		
			}
			if (df(F/cos(t),t) == 0)  // meaning that the rhs_function is in the form of a*cos(t)
			{
				double a = F/cos(t);
				Symbolic lhs_nh = Yt_particular_total.coeff(A*cos(t),1);
				A_solve = a/lhs_nh;
				Yt_final = A_solve*Yt_assumed.coeff(cos(t),1)*cos(t) ;		
			}
			
		}
		
		//higherorderlineardiffeq_nonhomogeneousequations_undeterminedcoefficients(Yt_particular_total, rhs_function, t);
						
		//cout << "\nY(t) assumed = "<< Yt_assumed << endl;
		//cout << "\n"<< Yt_particular_total << " = " << rhs_function << endl;
		cout << "\nThe general solution of the nonhomogeneous equation is: \nY(t) = "<< Yt_final << endl;
			
		cout << "\nThe general solution of the differential equation is: \nY(t) = "<< general_solution + Yt_final<< endl;
}

void higherorderlineardiffeq_nonhomogeneousequations_undeterminedcoefficients(const Symbolic &Yt_particular_total, const Symbolic &rhs_function, const Symbolic &t)
{
	// a pain in the ass with equation list doesn't contain lhs
		Symbolic A("A");
		Symbolic A_solve;
		cout << "\n rhs = "<<  rhs_function << endl;

		if(rhs_function != 0 )
		{
			list<Equations> eq;
			list<Equations>::iterator i;
			UniqueSymbol a, b;
			// Case 1 : g(t) = a*exp(b*t)
			eq = (a*exp(b*t)).match(rhs_function, (a,b));
			
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i, b); // equation list doesn't contain lhs
			Symbolic lhs_nh = Yt_particular_total.coeff(exp(bp*t),1);
			A_solve = solve(lhs_nh*exp(bp*t)- rhs_function, A).front().rhs;
			//Yt_final = A_solve*Yt_assumed ;
			} catch(const SymbolicError &se) {}
			}

			// Case 2 : g(t) = a*exp(t)
			eq = (a*exp(t)).match(rhs_function, (a,b));
			
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i, b);  // equation list doesn't contain lhs
				
			Symbolic lhs_nh = Yt_particular_total.coeff(exp(t),1);
			A_solve = solve(lhs_nh*exp(t)- rhs_function, A).front().rhs;
			cout << "\n A solve = "<<  A_solve << endl;
			//Yt_final = A_solve*Yt_assumed ;
			} catch(const SymbolicError &se) {}
			}

		}
}

void higherorderlineardiffeq_nonhomogeneousequationsgeneralsolution(const vector<complex<double>> &P, const SymbolicMatrix &Matrix_A, const Symbolic &t)
{ //Code it on April 12th, 2026
	
	Symbolic A("A"), Yt_particular, Yt_particular_total;
	Symbolic Yt_solution, Yt_final, Yt_final_total, Yt_assumed;;
	Symbolic general_solution;
		
	vector<complex<double>> vec_x0;
	int n_vec = P.size() ;
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
	for(int i=0; i<n_vec-1; i++)
	{
		double real_part = distribution(generator);
		double imag_part = 0;
		complex<double> random_complex(real_part, imag_part);
		vec_x0.push_back(random_complex); 	
	}

	int N = 100;
	int n_Polynomial = P.size();
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	int m_index3;
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_imagroots;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}

	cout << "\nP: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}

	cout << "\nInitial guess for the roots (generated randomly): " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		// To show the process of the Abert-Ehrlich
		//cout <<"\niteration: " << k << endl;
		//cout << "\nz_{i} new: " << endl;
		//printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
	// We use lround because there is an occurence if the root is obtained at very small decimal 
	// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

		if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2),roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
	// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
	// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
	vector<std::complex<double>> vec_unique;
	vector<std::complex<double>> vec_duplicate;

	int m = vec_update.size();
	for (int i = 0 ; i < m ; ++i)
	{
		double a  = real(vec_update[i]);
		double b  = imag(vec_update[i]);

		std::complex<double> target(a, b);

		// Get number of occurrences
		long count = std::count(vec_update.begin(), vec_update.end(), target);

		//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
	
		if(count == 1)
		{
			vec_unique.push_back(vec_update[i]);
		}
		else if(count > 1)
		{
			vec_duplicate.push_back(vec_update[i]);
		}
	}
	//cout << "\nUnique vector:" << std::endl;
	//printComplexVector(vec_unique);
	//cout << "\nDuplicate vector:" << std::endl;
    	//printComplexVector(vec_duplicate);

	// End of splitting into duplicate and unique vectors

	int n_unique = vec_unique.size();
	int n_duplicate = vec_duplicate.size();
	
	// This is for the case when the roots are unique, no duplicate / repeated roots.
	if(n_duplicate == 0)
	{
		Symbolic c("c");
		for(int i = 0; i < n;++i)
		{
			complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
			vec_imagroots.push_back(root) ;
		}
		//cout << "\nThe abs imag roots" << endl ;
		//printComplexVector(vec_imagroots);

		// 1. Sort the vector using a custom comparator
		
		std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
		{
		if (a.real() != b.real()) 
		{
			return a.real() < b.real();
		}
		else
		{	
			return a.imag() < b.imag();
		}
		});	
		setprecision(5);
		//    Use std::unique to move all non-duplicate elements to the front
		//    and return an iterator to the new logical end of the unique range.
		auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

		//    Erase the duplicate elements from the end of the vector.
		vec_imagroots.erase(last, vec_imagroots.end());

		// FInd a way to delete duplicate root / the complex conjugate.
		//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
		//printComplexVector(vec_imagroots);
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			//cout << "mu = " << mu << endl;
			double lambda = real(vec_imagroots[i]);
			int k = 1;
			if(mu != 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
				i=i+1;

			}
			else if(mu == 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
			}
				Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
				
		}
		cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;
		Yt_particular = A*Yt_assumed;
		m_index3 = n_Polynomial-1;
		for(int i = 0; i <= n_Polynomial-1; ++i)
		{
			//cout <<"P[i] = " << real(P[i]) << endl;
			//cout <<"mindex3 = " << m_index3 << endl;
			//cout << "Yt^{(n)} = " << df(Yt_particular,t,m_index3)*real(P[i]) << endl;
			Yt_particular_total += df(Yt_particular,t,m_index3)*real(P[i]) ;
			m_index3 = m_index3-1;
			//cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
		}
	}
	else if(n_duplicate != 0 )
	{
		Symbolic c("c");
		if (n_unique != 0 && n_duplicate != n)
		{
			// We handle for the unique roots first
			for(int i = 0; i < n_unique;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			//cout << "\nThe abs imag unique roots after delete duplicate" << endl ;
			//printComplexVector(vec_imagroots);

			/*

				TEST 

			*/

			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			cout << "\nSorted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				//cout << "i = " << i << endl;
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "vec occurence = " << vec_occurence[i] << endl;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last2, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last3 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last3, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size() - n_unique;
			
			int index_i_continuing;
			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					
				}
				index_i_continuing = i+1;
			}
			int m_index = index_i_continuing;
			int i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
					
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j < n_occurence + m_index ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
					
				}
				i_occurence = i_occurence+1;
			}
			cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;

			Yt_particular = A*Yt_assumed;
			m_index3 = n_Polynomial-1;

			for(int i = 0; i <= n_Polynomial-1; ++i)
			{
				//cout <<"P[i] = " << real(P[i]) << endl;
				//cout <<"mindex3 = " << m_index3 << endl;
				//cout << "Yt^{(n)} = " << df(Yt_particular,t,m_index3)*real(P[i]) << endl;
				Yt_particular_total += df(Yt_particular,t,m_index3)*real(P[i]) ;
				m_index3 = m_index3-1;
				//cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
			}

			/*

				END OF TEST 

			*/			
		}
		else if (n_duplicate == n)
		{
			//cout << "\nn unique == 0 "<< endl;
			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
		
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			//cout << "\nSorted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last2, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size();

			int m_index2 = 0;
			
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				//cout << "n_occurence = " << n_occurence  << endl;
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
					Yt_assumed = exp(lambda*t)* cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t) * sin(mu*t) * pow(t,Symbolic(k));
				}
				
			}
			cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;

			m_index3 = n_Polynomial-1;
			Yt_particular = A*Yt_assumed;
			for(int i = 0; i <= n_Polynomial-1; ++i)
			{
				//cout <<"P[i] = " << real(P[i]) << endl;
				//cout <<"mindex3 = " << m_index3 << endl;
				//cout << "Yt^{(n)} = " << df(Yt_particular,t,m_index3)*real(P[i]) << endl;
				Yt_particular_total += df(Yt_particular,t,m_index3)*real(P[i]) ;
				m_index3 = m_index3-1;
				//cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
			}
						
		}	
	}
	
	int n_row = Matrix_A.rows();


	for (int i = 0; i < n_row ; ++i)
	{
		Symbolic rhs_function = Matrix_A[i][0];
		cout << "\n***********************************************************"<< endl;
		cout << "\nFor rhs = " << rhs_function << endl;
			
		Symbolic A_solve;
		Symbolic F = rhs_function;
		cout << "\nY(t) particular = "<< Yt_particular_total  << endl;
		cout << "\nY(t) assumed = "<< Yt_assumed << endl;

		// Case 1 : g(t) = a*exp(t)
		if (df(F,t,0) == df(F,t,1) &&  df(F,t,1) == df(F,t,2) && df(F,t,0) == df(F,t,2)) // meaning that the rhs_function is in the form of a*exp(t)
		{
			//double a = df(F,t,0)/exp(t);
			Symbolic lhs_nh = Yt_particular_total.coeff(exp(t),1);
			A_solve = solve(lhs_nh*exp(t)- rhs_function, A).front().rhs;
			Yt_final = A_solve*Yt_assumed ;		
		}
		// Case 2 : g(t) = a*sin(t) or g(t) = a*cos(t)
		if (df(F,t,0) == -df(F,t,2) &&  df(F,t,0) != df(F,t,1) && df(F,t,0) == df(F,t,4)) // meaning that the rhs_function is in the form of a*sin(t) or a* cos(t)
		{
			if (df(F/sin(t),t) == 0)  // meaning that the rhs_function is in the form of a*sin(t)
			{
				double a = F/sin(t);
				Symbolic lhs_nh = Yt_particular_total.coeff(A*sin(t),1);
				A_solve = a/lhs_nh;
				Yt_final = A_solve*Yt_assumed.coeff(sin(t),1)*sin(t) ;		
			}
			if (df(F/cos(t),t) == 0)  // meaning that the rhs_function is in the form of a*cos(t)
			{
				double a = F/cos(t);
				Symbolic lhs_nh = Yt_particular_total.coeff(A*cos(t),1);
				A_solve = a/lhs_nh;
				Yt_final = A_solve*Yt_assumed.coeff(cos(t),1)*cos(t) ;		
			}
			
		}
		Yt_final_total += Yt_final;
		cout << "\nThe general solution of this nonhomogeneous equation is: \nY_{i}(t) = "<< Yt_final << endl;
	}

	cout << "\nThe general solution of the differential equation is: \nY(t) = "<< general_solution + Yt_final_total<< endl;

}

void higherorderlineardiffeq_nonhomogeneousequations_variationofparameters(const vector<complex<double>> &P, const Symbolic &rhs_function, const Symbolic &t)
{
	Symbolic A("A");
	Symbolic general_solution, Yt_final;
	Symbolic c("c");
	vector<complex<double>> vec_x0;
	int n_vec = P.size() ;
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(5, 1.2); // mu = 5, sigma = 1.2
	for(int i=0; i<n_vec-1; i++)
	{
		double real_part = distribution(generator);
		double imag_part = 0;
		complex<double> random_complex(real_part, imag_part);
		vec_x0.push_back(random_complex); 	
	}

	int N = 100;
	int n_Polynomial = P.size();
	Matrix<Symbolic> W(n_Polynomial-1,n_Polynomial-1);
	Matrix<Symbolic> Mat_particular(n_Polynomial-1,1);
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	int m_index3;
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_imagroots;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}

	cout << "\nP: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}

	cout << "\nInitial guess for the roots (generated randomly): " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		// To show the process of the Abert-Ehrlich
		//cout <<"\niteration: " << k << endl;
		//cout << "\nz_{i} new: " << endl;
		//printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
	// We use lround because there is an occurence if the root is obtained at very small decimal 
	// if a root obtained is like this: 1.00000004575, and another root is : 0.9999999765,  it is hard to split them into duplicate and unique vector without lround

		if(abs(imag(vec_dummy[i])) < 1e-8 && abs(real(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(0.0,roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-8 && abs(imag(vec_dummy[i])) > 1e-8)
		{
			complex<double> root(roundToDecimal(real(vec_dummy[i]),2),roundToDecimal(imag(vec_dummy[i]),2));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
	// Splitting vec_update into unique vector(vector with unique element) and duplicate vector (vector with element that occurs more than 1)
	// Complex Equality: std::complex uses operator== which checks if both real and imaginary parts are equal.
	vector<std::complex<double>> vec_unique;
	vector<std::complex<double>> vec_duplicate;

	int m = vec_update.size();
	for (int i = 0 ; i < m ; ++i)
	{
		double a  = real(vec_update[i]);
		double b  = imag(vec_update[i]);

		std::complex<double> target(a, b);

		// Get number of occurrences
		long count = std::count(vec_update.begin(), vec_update.end(), target);

		//std::cout << "Element "  << i << "-th occurs " << count << " times." << std::endl;
	
		if(count == 1)
		{
			vec_unique.push_back(vec_update[i]);
		}
		else if(count > 1)
		{
			vec_duplicate.push_back(vec_update[i]);
		}
	}
	//cout << "\nUnique vector:" << std::endl;
	//printComplexVector(vec_unique);
	//cout << "\nDuplicate vector:" << std::endl;
    	//printComplexVector(vec_duplicate);

	// End of splitting into duplicate and unique vectors

	int n_unique = vec_unique.size();
	int n_duplicate = vec_duplicate.size();
	
	// This is for the case when the roots are unique, no duplicate / repeated roots.
	if(n_duplicate == 0)
	{
		Symbolic c("c");
		for(int i = 0; i < n;++i)
		{
			complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
			vec_imagroots.push_back(root) ;
		}
		//cout << "\nThe abs imag roots" << endl ;
		//printComplexVector(vec_imagroots);

		// 1. Sort the vector using a custom comparator
		
		std::sort(vec_imagroots.begin(), vec_imagroots.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
		{
		if (a.real() != b.real()) 
		{
			return a.real() < b.real();
		}
		else
		{	
			return a.imag() < b.imag();
		}
		});	
		setprecision(5);
		//    Use std::unique to move all non-duplicate elements to the front
		//    and return an iterator to the new logical end of the unique range.
		auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

		//    Erase the duplicate elements from the end of the vector.
		vec_imagroots.erase(last, vec_imagroots.end());

		// FInd a way to delete duplicate root / the complex conjugate.
		//cout << "\nThe abs imag roots after delete duplicate" << endl ; 
		//printComplexVector(vec_imagroots);
		for(int i = 0; i < n;++i)
		{
			double mu = imag(vec_imagroots[i]);
			//cout << "mu = " << mu << endl;
			double lambda = real(vec_imagroots[i]);
			
			if(mu != 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
				W[0][i] = (exp(lambda*t)*c[i] *cos(mu*t))/c[i]; 
				W[0][i+1] = (exp(lambda*t)*c[i+1] *sin(mu*t) )/c[i+1]; 
				
				i=i+1;

			}
			else if(mu == 0)
			{
				general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
				W[0][i] = (exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t))/c[i]; 

			}
				
		}
		cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;
		
		
	}
	else if(n_duplicate != 0 )
	{
		
		if (n_unique != 0 && n_duplicate != n)
		{
			// We handle for the unique roots first
			for(int i = 0; i < n_unique;++i)
			{
				complex<double> root(real(vec_unique[i]),abs(imag(vec_unique[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last, vec_imagroots.end());

			//cout << "\nThe abs imag unique roots after delete duplicate" << endl ;
			//printComplexVector(vec_imagroots);

			/*

				TEST 

			*/

			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			cout << "\nSorted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				//cout << "i = " << i << endl;
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "vec occurence = " << vec_occurence[i] << endl;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last2, vec_duplicate.end());
			//cout << "\nDeleted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last3 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last3, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size() - n_unique;
			
			int index_i_continuing;
			for(int i = 0; i < n_unique;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				if(mu != 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i+1] *sin(mu*t) ;
					W[0][i] = (exp(lambda*t)*c[i] *cos(mu*t))/c[i]; 
					W[0][i+1] = (exp(lambda*t)*c[i+1] *sin(mu*t) )/c[i+1]; 
					i=i+1;
				}
				else if(mu == 0)
				{
					general_solution += exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) ;
					W[0][i] = (exp(lambda*t)*c[i] *cos(mu*t) + exp(lambda*t)*c[i] *sin(mu*t) )/c[i]; 
				}
				index_i_continuing = i+1;
			}
			int m_index = index_i_continuing;
			int i_occurence = 0;
			for(int i = index_i_continuing; i < index_i_continuing + n_duplicatefinal;++i)
			{
				//cout <<" i = "<< i << endl;
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i_occurence];
				//cout << "n occurence = " <<  n_occurence << endl;
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index; j < (2*n_occurence) + m_index; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						W[0][j] = (exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)))/c[j]; 
						W[0][j+1] = (exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k)))/c[j+1]; 
						k = k+1;		
						j = j+1;				
					}
					m_index = (2*n_occurence) + m_index;
					
					
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index; j <  m_index+n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;
						W[0][j] = (exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k)))/c[j];		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index = m_index+n_occurence;
					
					
				}
				i_occurence = i_occurence + 1;
			}
			cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;

			m_index3 = n_Polynomial-1;

			/*

				END OF TEST 

			*/			
		}
		else if (n_duplicate == n)
		{
			//cout << "\nn unique == 0 "<< endl;
			vector<int> vec_occurence;
			
			// 1. Sort the vector using a custom comparator
		
			std::sort(vec_duplicate.begin(), vec_duplicate.end(), [](const std::complex<double>& a, const std::complex<double>& b) 
			{
			if (a.real() != b.real()) 
			{
				return a.real() < b.real();
			}
			else
			{	
				return a.imag() < b.imag();
			}
			});	
			//cout << "\nSorted Duplicate vector:" << endl;
    			//printComplexVector(vec_duplicate);
			int m_stop;

			for (int i = 0 ; i < n_duplicate ; ++i)
			{
				double a  = real(vec_duplicate[i]);
				double b  = imag(vec_duplicate[i]);

				std::complex<double> target(a, b);

				// Get number of occurrences
				int count = std::count(vec_duplicate.begin(), vec_duplicate.end(), target);

				//cout << "Element "  << i << "-th occurs " << count << " times." << endl;
				vec_occurence.push_back(count);
				m_stop = std::accumulate(vec_occurence.begin(), vec_occurence.end(), 0) ;
				//cout << "m stop = " << m_stop << endl;
				if (m_stop == n_duplicate )
				{
					i = n_duplicate-1;		
				}
			}
			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last = std::unique(vec_duplicate.begin(), vec_duplicate.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_duplicate.erase(last, vec_duplicate.end());
			cout << "\nDeleted Duplicate vector:" << endl;
    			printComplexVector(vec_duplicate);

			int n_duplicate2 = vec_duplicate.size();
			// Remove the complex conjugate and store the last final root/s in vec_imagroots
			for(int i = 0; i < n_duplicate2;++i)
			{
				complex<double> root(real(vec_duplicate[i]),abs(imag(vec_duplicate[i])));
				vec_imagroots.push_back(root) ;
			}

			//    Use std::unique to move all non-duplicate elements to the front
			//    and return an iterator to the new logical end of the unique range.
			auto last2 = std::unique(vec_imagroots.begin(), vec_imagroots.end());

			//    Erase the duplicate elements from the end of the vector.
			vec_imagroots.erase(last2, vec_imagroots.end());
			cout << "\nFinal root vector:" << endl;
    			printComplexVector(vec_imagroots);
			int n_duplicatefinal = vec_imagroots.size();

			int m_index2 = 0;
			
			for(int i = 0; i < n_duplicatefinal;++i)
			{
				double mu = imag(vec_imagroots[i]);
				double lambda = real(vec_imagroots[i]);
				int n_occurence = vec_occurence[i];
				//cout << "n_occurence = " << n_occurence  << endl;
				if(mu != 0)
				{					
					int k = 0;
					for(int j = m_index2; j < (2*n_occurence) + m_index2 ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k));
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						W[0][j]=(exp(lambda*t)*c[j] *cos(mu*t) * pow(t,Symbolic(k)))/c[j];
						W[0][j+1]=(exp(lambda*t)*c[j+1] *sin(mu*t) * pow(t,Symbolic(k)))/c[j+1];

						k = k+1;		
						j = j+1;				
					}
					m_index2 = (2*n_occurence) + m_index2;
				}
				else if(mu == 0)
				{
					int k = 0;
					for(int j = m_index2; j < n_occurence  ; j++)
					{
						//cout <<" j = "<< j << endl;
						general_solution += exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k))  ;
						W[0][j] = (exp(lambda*t)*c[j] * cos(mu*t) * pow(t,Symbolic(k)) + exp(lambda*t)*c[j] * sin(mu*t) * pow(t,Symbolic(k)))/c[j];		
						//cout << "\nThe general solution is: \ny(t) = "<< general_solution << endl;
						k = k+1;		
					}
					m_index2 = m_index2 + n_occurence;
				}
				
			}
			cout << "\nThe general solution of the homogeneous equation is: \ny(t) = "<< general_solution << endl;

			m_index3 = n_Polynomial-1;
						
		}	
		
	}
		
		for (int j = 0; j < m_index3; ++j)
		{
			for (int i = 1; i < m_index3; ++i)
			{
				W[i][j] = df(W[i-1][j],t);
			}
		}
		cout << "\nW(t):\n" << W <<endl;
					
		cout << "\ndet (W(t)) = " << W.determinant() <<endl;
		for(int j = 0; j < m_index3 ; ++j)
		{		
			Matrix<Symbolic> Wi(m_index3,m_index3);
			for (int j = 0; j < m_index3; ++j)
			{
				for (int i = 0; i < m_index3; ++i)
				{
					Wi[i][j] = W[i][j];
				}
			}
			for (int i = 0; i < m_index3; ++i)
			{
				Wi[i][j] = 0;
				if (i== m_index3-1)
				{
					Wi[i][j] = 1;
				}
			}
			cout << "\nW_{" << j+1 << "}(t):\n" << Wi <<endl;
		
			cout << "\ndet (W_{" << j+1 << "}(t)):\n" << Wi.determinant() <<endl;
			Mat_particular[j][0] = Wi.determinant();
		}
		//cout << "\nM(t):\n" << Mat_particular <<endl;
		for(int i = 0; i < m_index3 ; ++i)
		{	
			Mat_particular[i][0] = Mat_particular[i][0]*rhs_function/(W.determinant());
		}
		//cout << "\nM(t):\n" << Mat_particular <<endl;
		
		for(int i = 0; i < m_index3 ; ++i)
		{	
			//cout << "\nW_{i}(t):\n" << W[0][i] <<endl;
			Yt_final += W[0][i]*integrate(Mat_particular[i][0],t);
		}
		cout << "\nThe particular solution of the differential equation is: \nY(t) = "<< Yt_final << endl;
}

// Helper structure to compute combinations for polynomial shifting
long long binomialCoefficient(int n, int k) 
{
	if (k < 0 || k > n) return 0;
	if (k == 0 || k == n) return 1;
	long long res = 1;
	for (int i = 1; i <= k; ++i) 
	{
		res = res * (n - i + 1) / i;
	}
	return res;
}

/* 

	Initialize class to compute the series solution for homogeneous linear differential equation
	with constant coefficients 

*/

// Constructor initializes the terms (series degree+1) and base boundary conditions
HigherOrderODE_Homogeneous_PowerSeriesSolver::HigherOrderODE_Homogeneous_PowerSeriesSolver(const vector<double>& ode_coeffs, const vector<double>& init_conditions) 
{
	ode_coefficients = ode_coeffs;
	initial_values = init_conditions;
	order = ode_coeffs.size() - 1;
	
}

void HigherOrderODE_Homogeneous_PowerSeriesSolver::computeSeries(int terms) 
{
	if (terms < 2) 
	{
		return;
       	}

	coefficients.resize(terms,0.0);

	// Step 1: Assign initial conditions to the first m series terms
	// Recall that y^(k)(0) = k! * coefficients_k, so coefficients_k = y^(k)(0) / k!
	double factorial = 1.0;
		for (int k = 0; k <= order - 1 && k < terms; ++k) 
	{
		if (k > 0) 
		{
		factorial *= k;
		}
		coefficients[k] = initial_values[k] / factorial;
	}

	// Step 2: Use recurrence relation to compute subsequent coefficients
	double b_m = ode_coefficients[order];
	if (b_m == 0.0) 
	{
		std::cerr << "Error: Leading coefficient b_m cannot be zero." << endl;
	}

	for (int n = 0; n <= terms - order - 1; ++n) 
	{
		double sum_terms = 0.0;
		for (int k = 0; k < order; ++k) 
		{
			double mult = get_factorial_multiplier(n, k);
			sum_terms += ode_coefficients[k] * coefficients[n + k] * mult;
		}
		    
		double divisor = b_m * get_factorial_multiplier(n, order);
		coefficients[n + order] = -sum_terms / divisor;
	}
	// To show the differential equation nicely
	cout << "\n( " << ode_coefficients[0] << " ) y ";
	for(int i = 1; i <= order ; ++i)
	{
		cout <<  "+ ( " << ode_coefficients[i] << " ) y^(" << i << ")";
	}
	cout << " = 0 " << endl;

	cout << "\nInitial conditions:"<< endl;
	cout << "y(0) = " << initial_values[0] << endl;
	for(int i = 1; i < order ; ++i)
	{
		cout << "y^(" << i << ") (0) = " << initial_values[i] << endl;
	}
}

vector<double> HigherOrderODE_Homogeneous_PowerSeriesSolver::coefficientsvector(int terms) 
{
	if (terms < 2) 
	{
		return {};
       	}

	coefficients.resize(terms,0.0);

	// Step 1: Assign initial conditions to the first m series terms
	// Recall that y^(k)(0) = k! * coefficients_k, so coefficients_k = y^(k)(0) / k!
	double factorial = 1.0;
		for (int k = 0; k <= order - 1 && k < terms; ++k) 
	{
		if (k > 0) 
		{
		factorial *= k;
		}
		coefficients[k] = initial_values[k] / factorial;
	}

	// Step 2: Use recurrence relation to compute subsequent coefficients
	double b_m = ode_coefficients[order];
	if (b_m == 0.0) 
	{
		std::cerr << "Error: Leading coefficient b_m cannot be zero." << endl;
		return {};
	}

	for (int n = 0; n <= terms - order - 1; ++n) 
	{
		double sum_terms = 0.0;
		for (int k = 0; k < order; ++k) 
		{
			double mult = get_factorial_multiplier(n, k);
			sum_terms += ode_coefficients[k] * coefficients[n + k] * mult;
		}
		    
		double divisor = b_m * get_factorial_multiplier(n, order);
		coefficients[n + order] = -sum_terms / divisor;
	}

	return coefficients;
}

// Evaluates the power series at a specific value of x
double HigherOrderODE_Homogeneous_PowerSeriesSolver::evaluateAt(double x, int terms)  
{
	double current_x_power = 1.0;
	double sum = 0.0;
        
	for (int i = 0; i <= terms; ++i) 
	{
		sum += coefficients[i] * current_x_power;
		current_x_power *= x;
	}
	return sum;
}

void HigherOrderODE_Homogeneous_PowerSeriesSolver::printSeries() const 
{
	cout << "\nSeries solution: \n"<< endl;
	cout << "y(x) = ";
	bool first = true;
	for (size_t i = 0; i < coefficients.size(); ++i) 
	{
		if (coefficients[i] == 0.0) 
		{
			continue;
		}
		if (!first && coefficients[i] > 0) 
		{
			cout << " + ";
		}		
		if (coefficients[i] < 0) 
		{
			cout << " - ";
		}
		double abs_val = std::abs(coefficients[i]);
		if (i == 0) 
		{
			cout << abs_val;
		} 
		else if (i == 1) 
		{
			cout << abs_val << "x";
		} 
		else 
		{
			cout << abs_val << "x^" << i;
		}
 		first = false;
        }
	cout << " + ...\n";
}

/* 

	Initialize class to compute the series solution for homogeneous second order linear differential equation
	Solves p(x)y'' + q(x)*y' + r(x)*y = 0 with variable coefficients


*/


// Constructor initializes the terms (series degree+1) and base boundary conditions
SecondOrderODE_Homogeneous_PowerSeriesSolver::SecondOrderODE_Homogeneous_PowerSeriesSolver(const PolynomialDouble& P_input, 
	const PolynomialDouble& Q_input, const PolynomialDouble& R_input, double x0_input, double y0_input, double dy0_input) 
{
	P = P_input;
	Q = Q_input;
	R = R_input;
	x0 = x0_input;
	// Initial conditions: y(0) = c_0, y'(0) = c_1
	y0 = y0_input;
	dy0 = dy0_input;
}

void SecondOrderODE_Homogeneous_PowerSeriesSolver::computeSeries(int terms) 
{

	coefficients.resize(terms, 0.0);
	coefficients[0] = y0;
	coefficients[1] = dy0;

	 // 1. Shift variable coefficients around x0: P(t+x0), Q(t+x0), R(t+x0)
	PolynomialDouble P_shifted = P.shift_around(x0);
	PolynomialDouble Q_shifted = Q.shift_around(x0);
	PolynomialDouble R_shifted = R.shift_around(x0);

	// Ensure the point is ordinary (P(x0) cannot be zero)
	double p0 = P_shifted.get_coeff(0);
	if (std::abs(p0) < 1e-9) 
	{
		std::cerr << "Error: x0 is a singular point. P(x0) cannot be 0.\n";
	}

	// Iteratively determine c_n using the algebraic recurrence relation
	// The equation coefficient for X^m in P(X)y'' + Q(X)y' + R(X)y = 0 must equal 0
	for (int m = 0; m < terms - 2; ++m) 
	{
		double sum_terms = 0.0;

		// Contribution from P(X)y''
		for (int k = 0; k <= m; ++k) 
		{
			int n = k + 2;
			if (n < terms) 
			{
				sum_terms += P_shifted.get_coeff(m - k) * n * (n - 1) * coefficients[n];
			}
		}

		// Contribution from Q(X)y'
		for (int k = 0; k <= m; ++k) 
		{
			int n = k + 1;
			if (n < terms) 
			{
				sum_terms += Q_shifted.get_coeff(m - k) * n * coefficients[n];
			}
		}

		// Contribution from R(X)y
		for (int k = 0; k <= m; ++k) 
		{
			int n = k;
			if (n < terms) 
			{
				sum_terms += R_shifted.get_coeff(m - k) * coefficients[n];
			}
		}

		// The term containing the unknown coefficient c[m+2] is isolated:
		// p0 * (m+2) * (m+1) * c[m+2] + sum_terms_excluding_this_one = 0
		// Therefore, we can correct the sum by tracking how much c[m+2] contributed dynamically
		// and solving for it directly.
		
		// Alternatively, calculate next term explicitly by peeling off the c[m+2] multiplier:
		double known_sum = 0.0;
		
		// P(X)y'' parts up to c[m+1]
		for (int j = 1; j <= m; ++j) 
		{
			known_sum += P_shifted.get_coeff(j) * (m - j + 2) * (m - j + 1) * coefficients[m - j + 2];
		}
		// Q(X)y' parts up to c[m+1]
		for (int j = 0; j <= m; ++j) 
		{
			known_sum += Q_shifted.get_coeff(j) * (m - j + 1) * coefficients[m - j + 1];
		}
		// R(X)y parts up to c[m]
		for (int j = 0; j <= m; ++j) 
		{
			known_sum += R_shifted.get_coeff(j) * coefficients[m - j];
		}

		// Solve for c[m+2]
		coefficients[m + 2] = -known_sum / (p0 * (m + 2) * (m + 1));
	}

	// To show the differential equation nicely
	bool first = true;

	cout << "\n\nDifferential Equation: ( ";
	for (int i = 0; i < P.maxDegree()+1 ; ++i) 
	{
		if (P.sum_coeff() == 0.0) 
		{
			cout << "0";
			i = P.maxDegree();
			continue;
		}

		if (P.get_coeff(i) == 0.0) 
		{
			continue;
		}
		if (!first && P.get_coeff(i) > 0) 
		{
			cout << " + ";
		}		
		if (P.get_coeff(i) < 0) 
		{
			cout << " - ";
		}
		double abs_val = std::abs(P.get_coeff(i));

		if (i == 0) 
		{	
			cout << abs_val;
		} 
		else if (i == 1) 
		{
			cout << abs_val << "x";
		} 
		else 
		{
			cout << abs_val << "x^" << i;
		}
 		first = false;
		
        }
	cout << " ) y'' + ( ";
	first = true;
	for (int i = 0; i < Q.maxDegree()+1  ; ++i) 
	{
		if (Q.sum_coeff() == 0.0) 
		{
			cout << "0";
			i = Q.maxDegree();
		}
		if (Q.get_coeff(i) == 0.0) 
		{
			continue;
		}
		if (!first && Q.get_coeff(i) > 0) 
		{
			cout << " + ";
		}		
		if (Q.get_coeff(i) < 0) 
		{
			cout << " - ";
		}
		double abs_val = std::abs(Q.get_coeff(i));
		if (i == 0) 
		{
			cout << abs_val;
		} 
		else if (i == 1) 
		{
			cout << abs_val << "x";
		} 
		else 
		{
			cout << abs_val << "x^" << i;
		}
 		first = false;
		
        }
	cout << " ) y' + ( ";

	first = true;
	for (int i = 0; i < R.maxDegree()+1 ; ++i) 
	{
		if (R.sum_coeff() == 0.0) 
		{
			cout << "0";
			i = R.maxDegree();
		}
		if (R.get_coeff(i) == 0.0) 
		{
			continue;
		}
		if (!first && R.get_coeff(i) > 0) 
		{
			cout << " + ";
		}		
		if (R.get_coeff(i) < 0) 
		{
			cout << " - ";
		}
		double abs_val = std::abs(R.get_coeff(i));
		if (i == 0) 
		{
			cout << abs_val;
		} 
		else if (i == 1) 
		{
			cout << abs_val << "x";
		} 
		else 
		{
			cout << abs_val << "x^" << i;
		}
 		first = false;
		
        }
	cout << " ) y ";
	cout << "\nInitial conditions: y(0) = " << y0 << " , y'(0) = " << dy0 << endl;
	cout << "\nx_{0} = " << x0 << std::endl;

}

// Alternative 2
/*
// Computes the power series coefficients for y = c0 + c1*(x-x0) + c2*(x-x0)^2 + ...
// Up to a specified degree N using Taylor's method on: P(x)y'' + Q(x)y' + R(x)y = 0
void SecondOrderODE_Homogeneous_PowerSeriesSolver::computeSeries(int terms) 
{
	// 1. Shift variable coefficients so they are expanded around u = (x - x0)
	PolynomialDouble Pu = P.shiftToCenter(x0);
	PolynomialDouble Qu = Q.shiftToCenter(x0);
	PolynomialDouble Ru = R.shiftToCenter(x0);

	// Ensure x0 is an ordinary point (P(x0) != 0)
	if (std::abs(Pu.evaluateAt(0.0)) < 1e-12) 
	{
		std::cerr << "Error: x0 is a singular point. P(x0) cannot be 0.\n";
	}

	// y_derivatives[n] stores the n-th derivative of y evaluated at x0: y^(n)(x0)
	std::vector<double> y_derivatives(terms, 0.0);
	y_derivatives[0] = y0;
	y_derivatives[1] = dy0;

	// Helper lambdas to fetch coefficient of u^k from shifted polynomials safely
	auto getCoeff = [](const PolynomialDouble& poly, int k) 
	{
		return (k >= 0 && k < static_cast<int>(poly.coeffs.size())) ? poly.coeffs[k] : 0.0;
	};

	// 2. Compute higher order derivatives recursively up to N using Leibniz rule
	// Differentiating the ODE (n-2) times yields a system for y^(n)(x0)
	for (int n = 2; n <= terms; ++n) 
	{
		int m = n - 2; // number of differentiations applied to the entire ODE
        
		double sum_P = 0.0;
		for (int k = 1; k <= m; ++k) 
		{
			sum_P += binomialCoefficient(m, k) * getCoeff(Pu, k) * y_derivatives[m - k + 2];
		}

		double sum_Q = 0.0;
		for (int k = 0; k <= m; ++k) 
		{
			sum_Q += binomialCoefficient(m, k) * getCoeff(Qu, k) * y_derivatives[m - k + 1];
		}

		double sum_R = 0.0;
		for (int k = 0; k <= m; ++k) 
		{
			sum_R += binomialCoefficient(m, k) * getCoeff(Ru, k) * y_derivatives[m - k];
		}

		double P_val = getCoeff(Pu, 0); // P(x0)
		y_derivatives[n] = -(sum_P + sum_Q + sum_R) / P_val;
	}

    // 3. Convert derivatives to power series coefficients: a_n = y^(n)(x0) / n!
	coefficients.resize(terms, 0.0);
	coefficients[0] = y0;
	coefficients[1] = dy0;
//    std::vector<double> series_coeffs(N + 1, 0.0);
	double factorial = 1.0;
	for (int n = 0; n <= terms; ++n) 
	{
		if (n > 0) 
		{
			factorial *= n;
		}
		coefficients[n] = y_derivatives[n] / factorial;
	}

} */

// Alternative 3
/*
void SecondOrderODE_Homogeneous_PowerSeriesSolver::computeSeries(int terms) 
{
	// 1. Shift variable coefficients around x0: P(t+x0), Q(t+x0), R(t+x0)
	PolynomialDouble p_sh = P.shift(x0);
	PolynomialDouble q_sh = Q.shift(x0);
	PolynomialDouble r_sh = R.shift(x0);

	// Check if x0 is an ordinary point
	double p0 = p_sh.get_coeff(0);
	if (std::abs(p0) < 1e-9) 
	{
		std::cerr << "Error: x0 = " << x0 << " is a singular point. This solver requires an ordinary point.\n";
	}

        // Initialize coefficient array for y(t) = c_0 + c_1*t + c_2*t^2 + ...
	coefficients.resize(terms, 0.0);
	coefficients[0] = y0;
	if (terms > 1) 
	{
		coefficients[1] = dy0;
	}

	// 2. Iteratively solve the recurrence relation for higher-order terms
	for (int m = 0; m < terms - 2; ++m) 
	{
		double sum = 0.0;

		// Contribution from Q(t)*y'
		for (int j = 0; j <= m; ++j) 
		{
			sum += q_sh.get_coeff(m - j) * (j + 1) * coefficients[j + 1];
		}

		// Contribution from R(t)*y
		for (int j = 0; j <= m; ++j) 
		{
			sum += r_sh.get_coeff(m - j) * coefficients[j];
		}

		// Contribution from P(t)*y'' (excluding the leading p0 term)
		for (int j = 1; j <= m + 1; ++j) 
		{
		sum += p_sh.get_coeff(m - j + 2) * (j + 1) * j * coefficients[j + 1];
		}

		// Recurrence formula derived by matching t^m coefficients
		coefficients[m + 2] = -sum / (p0 * (m + 2) * (m + 1));
	}
} */

vector<double> SecondOrderODE_Homogeneous_PowerSeriesSolver::coefficientsvector(int terms) 
{
	if (terms < 2) 
	{
		return {};
       	}

	coefficients.resize(terms, 0.0);
	coefficients[0] = y0;
	coefficients[1] = dy0;

	PolynomialDouble P_shifted = P.shift_around(x0);
	PolynomialDouble Q_shifted = Q.shift_around(x0);
	PolynomialDouble R_shifted = R.shift_around(x0);

	// Ensure the point is ordinary (P(x0) cannot be zero)
	double p0 = P_shifted.get_coeff(0);
	if (std::abs(p0) < 1e-9) 
	{
		std::cerr << "Error: x is a singular point. P(x) cannot be 0.\n";
	}

	// Iteratively determine c_n using the algebraic recurrence relation
	// The equation coefficient for X^m in P(X)y'' + Q(X)y' + R(X)y = 0 must equal 0
	for (int m = 0; m < terms - 2; ++m) 
	{
		double sum_terms = 0.0;

		// Contribution from P(X)y''
		for (int k = 0; k <= m; ++k) 
		{
			int n = k + 2;
			if (n < terms) 
			{
				sum_terms += P_shifted.get_coeff(m - k) * n * (n - 1) * coefficients[n];
			}
		}

		// Contribution from Q(X)y'
		for (int k = 0; k <= m; ++k) 
		{
			int n = k + 1;
			if (n < terms) 
			{
				sum_terms += Q_shifted.get_coeff(m - k) * n * coefficients[n];
			}
		}

		// Contribution from R(X)y
		for (int k = 0; k <= m; ++k) 
		{
			int n = k;
			if (n < terms) 
			{
				sum_terms += R_shifted.get_coeff(m - k) * coefficients[n];
			}
		}

		// The term containing the unknown coefficient c[m+2] is isolated:
		// p0 * (m+2) * (m+1) * c[m+2] + sum_terms_excluding_this_one = 0
		// Therefore, we can correct the sum by tracking how much c[m+2] contributed dynamically
		// and solving for it directly.
		
		// Alternatively, calculate next term explicitly by peeling off the c[m+2] multiplier:
		double known_sum = 0.0;
		
		// P(X)y'' parts up to c[m+1]
		for (int j = 1; j <= m; ++j) 
		{
			known_sum += P_shifted.get_coeff(j) * (m - j + 2) * (m - j + 1) * coefficients[m - j + 2];
		}
		// Q(X)y' parts up to c[m+1]
		for (int j = 0; j <= m; ++j) 
		{
			known_sum += Q_shifted.get_coeff(j) * (m - j + 1) * coefficients[m - j + 1];
		}
		// R(X)y parts up to c[m]
		for (int j = 0; j <= m; ++j) 
		{
			known_sum += R_shifted.get_coeff(j) * coefficients[m - j];
		}

		// Solve for c[m+2]
		// Isolate and extract c[m+2] using the recurrence alignment rule
		coefficients[m + 2] = -known_sum / (p0 * (m + 2) * (m + 1));
	}
	
	return coefficients;
}

void SecondOrderODE_Homogeneous_PowerSeriesSolver::printCoefficients() const 
{
	// Output calculated power series coefficients to the console
	cout << "\nComputed Power Series Coefficients:\n";
	for (size_t i = 0; i < coefficients.size(); ++i) 
	{
		cout << "c_{" << i << "} = " << std::setw(10) << coefficients[i] << "\n";
	}
}

void SecondOrderODE_Homogeneous_PowerSeriesSolver::printSolution() const 
{
	cout << "\nSeries solution: \n"<< endl;
	cout << "y(x) = ";
	bool first = true;
	int n_terms = 0;
	for (size_t i = 0; i < coefficients.size(); ++i) 
	{
		if (std::abs(coefficients[i]) < 1e-9) 
		{
			continue;
		}
		if (!first && coefficients[i] > 0) 
		{
			cout << " + ";
		}
		if (coefficients[i] < 0) 
		{
			cout << " - ";
		}
		cout << std::abs(coefficients[i]);
		if (i > 0) 
		{
			if (x0==0)
			{
				cout << "*x";
			}
			else if (x0 != 0)
			{
				cout << "*(x - " << x0 << ")";
			}
			if (i > 1) 
			{
				cout << "^" << i;
			}
		}
        first = false;
	n_terms += 1;
        }
	
	if (n_terms >= int(coefficients.size()))
	{
		cout << " + ... \n";
	}
	else if (n_terms < int(coefficients.size()) )
	{
		cout << " " ;
	}
}


// Evaluates the power series at a specific value of x
double SecondOrderODE_Homogeneous_PowerSeriesSolver::evaluateAt(double x, int terms)  
{
	double result = 0.0;
        for (int i = terms - 1; i >= 0; --i) 
	{
		result = result * x + coefficients[i]; // Horner's method for numeric stability
	}
        return result;

}


// Function to find the minimum index shift for a regular singular point
int find_L(const PolynomialComplex& P, const PolynomialComplex& Q, const PolynomialComplex& R) 
{
	int L = 1e9;
	for (int i = 0; i < int(P.coeffs.size()); ++i) 
	{
		if (abs(P.coeffs[i]) > 1e-9) 
		{ 
			L = min(L, i - 2); 
			break; 
		}
	}
	for (int i = 0; i < int(Q.coeffs.size()); ++i) 
	{
		if (abs(Q.coeffs[i]) > 1e-9) 
		{ 
			L = min(L, i - 1); 
			break; 
		}
	}
	for (int i = 0; i < int(R.coeffs.size()); ++i) 
	{
		if (abs(R.coeffs[i]) > 1e-9) 
		{ 
			L = min(L, i); 
			break; 
		}
	}
	return L;
}

/* 

	Initialize class to compute the series solution for homogeneous second order linear differential equation
	Solves p(x)y'' + q(x)*y' + r(x)*y = 0 with variable coefficients
	Near a Regular Singular Points
	with Frobenius method

*/
// Small epsilon value to manage floating-point accuracy with complex numbers
const double EPSILON = 1e-7;

bool is_near_zero(complex<double> val) 
{
	return abs(val) < EPSILON;
}
void SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::classify_ode()
{
	PolynomialComplex Pc = P_unshifted; 
	PolynomialComplex Qc = Q_unshifted; 
	PolynomialComplex Rc = R_unshifted; 

	// 1. Check if it's an ordinary point
	if (!is_near_zero(Pc.evaluateAt(x0))) 
	{
		cout << "--> Result: Ordinary Point.\n";
		cout << "    Solve using a standard Power Series (Taylor Series).\n\n";
	}

	// 2. It is a Singular Point. Let's find limits using algebraic reduction.
	// We factor out (x - x0) from P(x), Q(x), and R(x) to compute the analytic limits.
    
	// Divide P once and twice
	auto [P_div1, P_rem1] = Pc.divide_by_linear(x0);
	auto [P_div2, P_rem2] = P_div1.divide_by_linear(x0);

	auto [Q_div1, Q_rem1] = Qc.divide_by_linear(x0);
	auto [R_div1, R_rem1] = Rc.divide_by_linear(x0);

	// Limit p0 = lim (x-x0)*Q(x)/P(x)
	// If Q has at least the same multiplicity root as P_minus_1_factor, the limit exists.
	p0_exists = is_near_zero(P_rem1); 
	complex<double> p0(0.0,0.0);
	if (p0_exists) 	
	{
		// If P has root of mult 1, then lim = Q(x0) / P_div1(x0)
		complex<double> p_denom = P_div1.evaluateAt(x0);
		if (!is_near_zero(p_denom)) {
			p0 = Qc.evaluateAt(x0) / p_denom;
		} 
		else 
		{
			// P has higher multiplicity root, check if Q cancels it out
			if (is_near_zero(Q_rem1)) 
			{
				p0 = Q_div1.evaluateAt(x0) / P_div2.evaluateAt(x0);
			} 
			else 
			{
				p0_exists = false;
			}
		}
	}

	// Limit q0 = lim (x-x0)^2*R(x)/P(x)
	q0_exists = false;
	complex<double> q0(0.0,0.0);
	complex<double> p_denom_2 = P_div1.evaluateAt(x0);
    
	if (!is_near_zero(p_denom_2)) 
	{
		// P has single root, (x-x0)^2 * R / P vanishes to 0 because of the extra (x-x0) up top
		q0_exists = true;
		q0 = 0.0;
	} 
	else 
	{
		// P has at least multiplicity 2 root
		complex<double> p_denom_3 = P_div2.evaluateAt(x0);
		if (!is_near_zero(p_denom_3)) 
		{
			q0_exists = true;
			q0 = Rc.evaluateAt(x0) / p_denom_3;
		}
	}
	// Classify the differential equation first
	cout << "Analyzing ODE at point x0 = " << x0 << "\n";
	if (p0_exists && q0_exists) 
	{
		// 3. Check for specific Euler-Cauchy structure: a*(x-x0)^2 y'' + b*(x-x0) y' + c y = 0
		// Centered around x0, this implies deg(P)<=2, deg(Q)<=1, deg(R)==0 relative to shift
		if (Pc.maxDegree() == 2 && Qc.maxDegree() == 1 && Rc.maxDegree() == 0) 
		{
			cout << "--> Result: Euler-Cauchy Equation.\n";
			cout << "    Solve exactly using the auxiliary algebraic characteristic equation.\n";
			EulerCauchy = true;
			Frobenius = false;
			// Indicial Equation for Euler-Cauchy: m^2 + (p0 - 1)m + q0 = 0
			/*complex<double> b = p0 - 1.0;
			complex<double> c = q0;
			complex<double> disc = b*b - 4.0*c;
			complex<double> m1 = (-b + sqrt(disc)) / 2.0;
			complex<double>m2 = (-b - sqrt(disc)) / 2.0;
			cout << "    Indicial complex roots: m1 = " << m1 << ", m2 = " << m2 << "\n\n";*/
		} 
		else 
		{
			cout << "--> Result: Regular Singular Point.\n";
			cout << "    Solve using the Frobenius method.\n";
			EulerCauchy = false;
			Frobenius = true;
			/*complex<double>b = p0 - 1.0;
			complex<double> c = q0;
			complex<double> disc = b*b - 4.0*c;
			complex<double> r1 = (-b + sqrt(disc)) / 2.0;
		    	complex<double> r2 = (-b - sqrt(disc)) / 2.0;
			cout << "    Frobenius Indicial roots: r1 = " << r1 << ", r2 = " << r2 << "\n\n";*/
		}
	} 
	else 
	{
		cout << "--> Result: Irregular Singular Point.\n";
		cout << "    Frobenius method fails. Solutions are divergent or non-Frobenius series.\n\n";
	}
}

// Evaluates one Frobenius series solution for a given root  up to max_terms
void SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::computeSeriesCoefficients() 
{
	if(Frobenius)
	{
		
		/*
			case r1 != r2 and abs(r1-r2) is an integer
			Computes coefficients for the standard Frobenius series: y = x^r * sum(a_n * x^n)
		*/
		
		// Convert to complex polynomial using the conversion constructor
		PolynomialComplex Pc = P; 
		PolynomialComplex Qc = Q; 
		PolynomialComplex Rc = R; 
		
		cout << "\nP(x) y'' + Q(x)y'+ R(x) y = 0" <<endl ;

		cout << "\nP(x) = ";
		Pc.print();
		cout << "\nQ(x) = " ; 
		Qc.print();
		cout << "\nR(x) = ";
		Rc.print() ;

		cout << "\nx0 = " << x0 ;
		//cout << "\nreal(r1 - r2) = " << N << "\n";
		cout << "\nr1 = " << r1 << "\n";
		cout << "r2 = " << r2 << "\n\n";
		
		coefficients_root1.resize(max_terms,0.0);
		coefficients_root2.resize(max_terms,0.0);
		
		if (integer_diff && !complex_roots) 
		{
			
			int N = static_cast<int>(std::round(real(r1) - real(r2)));

			coefficients_root1[0] = 1.0; // Arbitrary normalization factor
			coefficients_root2[0] = 1.0; // Arbitrary initialization

			cout << "Case 1: Different roots; roots differ by an integer." << endl;

			// Compute the first solution y1(x)
			for (int n = 1; n < max_terms; ++n) 
			{
				complex<double> Rn = 0.0;
				complex<double> current_r = r1 + static_cast<double>(n);

				// Collect contributions from previous terms
				for (int k = 0; k < n; ++k) 
				{
					complex<double> r_k = r1 + static_cast<double>(k);
					int shift = n - k;

					complex<double> P_part = Pc.get_coeff(shift + 2) * r_k * (r_k - 1.0);
					complex<double> Q_part = Qc.get_coeff(shift + 1) * r_k;
					complex<double> R_part = Rc.get_coeff(shift);

					Rn += coefficients_root1[k] * (P_part + Q_part + R_part);
				}

				// Indicial-like denominator for the current step
				complex<double> P2 = Pc.get_coeff(2);
				complex<double> Q1 = Qc.get_coeff(1);
				complex<double> R0 = Rc.get_coeff(0);
				complex<double> denom = P2 * current_r * (current_r - 1.0) + Q1 * current_r + R0;

				if (abs(denom) < 1e-9) 
				{
					// Indicial zero encountered (expected for r2 when roots differ by integer)
					coefficients_root1[n] = 0.0; 
				} 
				else 
				{
					coefficients_root1[n] = -Rn / denom;
				}
			}
				
			// Compute the second solution y2(x)
			// Extracting regular singular point values
			int L = find_L(P, Q, R);

			Complex p2 = Pc.get_coeff(L+2);
			Complex q1 = Qc.get_coeff(L + 1);
			Complex r0 = Rc.get_coeff(L );

			// Set up Frobenius parameter around smaller root: r = r2 + eps
			Dual2 r(r2, 1.0, 0.0);

			vector<Dual2> b(max_terms);

			// For distinct roots differing by integer, initialize a0 = r - r2 = eps
			b[0] = Dual2(0.0, 1.0, 0.0);

			// Lambda to evaluate Indicial polynomial of the form: P2*K*(K-1) + Q1*K + R0
			auto indicial_poly = [&](Dual2 K) 
			{
				return p2 * K * (K - Complex(1.0, 0.0)) + q1 * K + r0;
			};

			// Compute coefficients using Frobenius recurrence relations
			for (int n = 1; n < max_terms; ++n) 
			{
				Dual2 sum(0.0, 0.0, 0.0);
				for (int k = 0; k < n; ++k) 
				{
					Dual2 K = r + Complex(double(k), 0.0);
					Complex p_term = Pc.get_coeff(n - k + 2);
					Complex q_term = Qc.get_coeff(n - k + 1);
					Complex r_term = Rc.get_coeff(n - k);

					Dual2 bracket = p_term * K * (K - Complex(1.0, 0.0)) + q_term * K + r_term;
					sum = sum + b[k] * bracket;
				}
				Dual2 denom = indicial_poly(r + Complex(double(n), 0.0));
				b[n] = (Complex(-1.0, 0.0) * sum) / denom;
			}
			C_log = b[N].a; 
			for (int n = 0; n < max_terms; ++n) 
			{
				coefficients_root2[n] =b[n].b;
			}
			
		}
		if (noninteger_diff && !complex_roots) 
		{
			coefficients_root1[0] = 1.0; // Arbitrary scaling
			coefficients_root2[0] = 1.0; // Arbitrary initialization

			cout << "Case 2: Different roots; roots differ by non-integer." << endl;
			// Compute the first solution y1(x)
			for (int n = 1; n < max_terms; ++n) 
			{
				complex<double> Rn = 0.0;
				complex<double> current_r = r1 + complex<double>(n);

				// Collect contributions from previous terms / convolve with shifted indices
				for (int k = 0; k < n; ++k) 
				{
					complex<double> r_k = r1 + complex<double>(k);
					int shift = n - k;

					complex<double> P_part = Pc.get_coeff(shift + 2) * r_k * (r_k - 1.0);
					complex<double> Q_part = Qc.get_coeff(shift + 1) * r_k;
					complex<double> R_part = Rc.get_coeff(shift);

					Rn += coefficients_root1[k] * (P_part + Q_part + R_part);

				}

				// Indicial-like denominator for the current step
				complex<double> P2 = Pc.get_coeff(2);
				complex<double> Q1 = Qc.get_coeff(1);
				complex<double> R0 = Rc.get_coeff(0);
				complex<double> denom = P2 * current_r * (current_r - 1.0) + Q1 * current_r + R0;
				
				if (abs(denom) < 1e-9) 
				{
					// Indicial zero encountered (expected for r2 when roots differ by integer)
					coefficients_root1[n] = 0.0; 
				} 
				else 
				{
					coefficients_root1[n] = -Rn / denom;
				}
			}

			// Compute the second solution y2(x)
			for (int n = 1; n < max_terms; ++n) 
			{
				
				complex<double> Rn = 0.0;
				complex<double> current_r = r2 + complex<double>(n);

				// Collect contributions from previous terms / convolve with shifted indices
				for (int k = 0; k < n; ++k) 
				{
					complex<double> r_k = r2 + complex<double>(k);
					int shift = n - k;

					complex<double> P_part = Pc.get_coeff(shift + 2) * r_k * (r_k - 1.0);
					complex<double> Q_part = Qc.get_coeff(shift + 1) * r_k;
					complex<double> R_part = Rc.get_coeff(shift);
					
					Rn += coefficients_root2[k] * (P_part + Q_part + R_part);
				}
				
				// Indicial-like denominator for the current step
				complex<double> P2 = Pc.get_coeff(2);
				complex<double> Q1 = Qc.get_coeff(1);
				complex<double> R0 = Rc.get_coeff(0);
				complex<double> denom = P2 * current_r * (current_r - 1.0) + Q1 * current_r + R0;

				if (abs(denom) < 1e-9) 
				{
					// Indicial zero encountered (expected for r2 when roots differ by integer)
					coefficients_root2[n] = 0.0; 
				} 
				else 
				{
					coefficients_root2[n] = -Rn / denom;
				}
			}
			
		}
		if (repeated_roots && !complex_roots) 
		{
			complex<double> p0 = Pc.get_coeff(2);
			complex<double> q0 = Qc.get_coeff(1);
			complex<double> r0 = Rc.get_coeff(0);

			coefficients_root1[0] = 1.0; // Arbitrary scaling
			coefficients_root2[0] = 0.0; // Arbitrary initialization
			cout << "Case 3: Repeated roots." << endl;

			//  Setup recurrence using Dual numbers for automatic differentiation
			// We compute a_n as a function of r: a_n(r) and its derivative a_n'(r)
			vector<Dual> a(max_terms);
			a[0] = Dual(1.0, 0.0); // Set a_0(r) = 1, so da_0/dr = 0

			// Compute the first solution y1(x) and the second solution y2(x), only need 1 for loop since it is repeated roots
			for (int k = 1; k < max_terms; ++k) 
			{
				Dual sum(0.0, 0.0);
				
				// Collect contributions from previous terms / convolve with shifted indices
				for (int j = 0; j < k; ++j) 
				{
					int shift = k - j;
		    
					complex<double> p_k = Pc.get_coeff(shift + 2);
					complex<double> q_k = Qc.get_coeff(shift + 1);
					complex<double> r_k = Rc.get_coeff(shift);

					// Compute the linear operator inside recurrence for variable r using Dual numbers
					// F(r + j) = p_k*(r+j)*(r+j-1) + q_k*(r+j) + r_k
					Dual r_plus_j(r1 + complex<double>(j), 1.0); // Variable r has a derivative of 1.0
					Dual F_val = Dual(p_k) * r_plus_j * (r_plus_j - Dual(1.0)) + Dual(q_k) * r_plus_j + Dual(r_k);
		               
					sum = sum + F_val * a[j];
				}

				// Indicial function denominator evaluated at (r + k)
				Dual r_plus_k(r1 + Complex(k), 1.0);
				Dual Denom = Dual(p0) * r_plus_k * (r_plus_k - Dual(1.0)) + Dual(q0) * r_plus_k + Dual(r0);

				// a_k(r) = - sum / Denom
				a[k] = Dual(0.0) - (sum / Denom);

				coefficients_root1[k] = a[k].val;
				coefficients_root2[k] = a[k].der;
				
			}

		}
		if (complex_roots) 
		{
			cout << "Case 4: Complex roots." << endl;
			coefficients_root1[0] = 1.0; // Arbitrary scaling
			coefficients_root2[0] = 1.0; // Arbitrary initialization

			// Compute the first solution y1(x)
			vector<Dual1> a(max_terms);
			a[0] = Dual1(1.0,0.0); // Standard choice for a_0

			Dual1 r1_dual(r1, 1.0); // Seed the root into the Dual1 number
			Dual1 p2 = Pc.get_coeff(2);
			Dual1 q1 = Qc.get_coeff(1);
			Dual1 r0 = Rc.get_coeff(0);

			for (int n = 1; n < max_terms; ++n) 
			{
				Dual1 sum(0.0);
				for (int k = 0; k < n; ++k) 
				{
					complex<double> p_term = Pc.get_coeff(n - k + 2);
					complex<double> q_term = Qc.get_coeff(n - k + 1);
					complex<double> r_term = Rc.get_coeff(n - k);

					Dual1 k_r = r1_dual + double(k);
					Dual1 factor = p_term * k_r * (k_r - 1.0) + q_term * k_r + r_term;
					sum = sum + factor * a[k];
				}
				Dual1 n_r = r1_dual + double(n);
				Dual1 denom = p2 * n_r * (n_r - 1.0) + q1 * n_r + r0;
				
				//cout << "sum = " << sum.val << endl;
				a[n] = (Dual1(0.0, 0.0) - sum) / denom;
				coefficients_root1[n] = a[n].val;
			}
		
			// Compute the second solution y2(x)
			vector<Dual1> b(max_terms);
			b[0] = Dual1(1.0,0.0); // Standard choice for b_0

			Dual1 r2_dual(r2, 1.0); // Seed the root into the Dual1 number

			for (int n = 1; n < max_terms; ++n) 
			{
				Dual1 sum(0.0);
				for (int k = 0; k < n; ++k) 
				{
					complex<double> p_term = Pc.get_coeff(n - k + 2);
					complex<double> q_term = Qc.get_coeff(n - k + 1);
					complex<double> r_term = Rc.get_coeff(n - k);

					Dual1 k_r = r2_dual + double(k);
					Dual1 factor = p_term * k_r * (k_r - 1.0) + q_term * k_r + r_term;
					sum = sum + factor * b[k];
				}
				Dual1 n_r = r2_dual + double(n);
				Dual1 denom = p2 * n_r * (n_r - 1.0) + q1 * n_r + r0;
				
				//cout << "sum = " << sum.val << endl;
				b[n] = (Dual1(0.0, 0.0) - sum) / denom;
				coefficients_root2[n] = b[n].val;
			}

		}
	}
	else if(EulerCauchy)
	{
		// Convert to complex polynomial using the conversion constructor
		PolynomialComplex Pc = P; 
		PolynomialComplex Qc = Q; 
		PolynomialComplex Rc = R; 

		int N = static_cast<int>(real(r1 - r2)); 
		integer_diff = (N > 0 && abs((r1 - r2) - (double)N) < 1e-7);
		noninteger_diff = ( N > 0 && abs((r1 - r2) - (double)N) > 1e-7 );
		
		cout << "\nP(x) y'' + Q(x)y'+ R(x) y = 0" <<endl ;

		cout << "\nP(x) = ";
		Pc.print();
		cout << "\nQ(x) = " ; 
		Qc.print();
		cout << "\nR(x) = ";
		Rc.print() ;

		cout << "\nx0 = " << x0 ;
		//cout << "\nreal(r1 - r2) = " << N << "\n";
		cout << "\nr1 = " << r1 << "\n";
		cout << "r2 = " << r2 << "\n\n";
		
		
		if (integer_diff && !complex_roots) 
		{
			cout << "Case 1: Different roots; roots differ by an integer." << endl;
			
		}
		if (noninteger_diff && !complex_roots) 
		{
			cout << "Case 2: Different roots; roots differ by non-integer." << endl;
			
		}
		if (repeated_roots && !complex_roots) 
		{
			cout << "Case 3: Repeated roots." << endl;

		}
		if (complex_roots) 
		{
			cout << "Case 4: Complex roots." << endl;
			

		}

	}
}

complex<double> SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::get_indicial_value(complex<double> rho) // necessary or not?
{
	complex<double> one(1,0.0);
	return P.get_coeff(2) * rho * (rho - one) + P.get_coeff(1) * rho + R.get_coeff(0);
}

void SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::solve_indicial_equation() 
{
	// 1.All arbitrary polynomials already shifted around the point x0
	

	// 2. Identify valuations (lowest non-zero powers)
	int vP = P.valuation();
	int vQ = Q.valuation();
	int vR = R.valuation();

	// 3. Check Regular Singularity Condition
	// If P(x0) != 0, it's an ordinary point, not a singular point. 
	// For it to be a regular singular point: vQ >= vP - 1 and vR >= vP - 2
	if (vP == 0 || vQ < vP - 1 || vR < vP - 2) 
	{
	std::cerr << "Error: Point x0 = " << x0 << " is not a regular singular point!\n";
        
	}

	// 4. Extract dominant structural coefficients
	double leadP = P.coeffs[vP];
	double leadQ = (vQ == vP - 1) ? Q.coeffs[vQ] : 0.0;
	double leadR = (vR == vP - 2) ? R.coeffs[vR] : 0.0;

	// 5. Compute limits p0 and q0 exactly
	double p0 = leadQ / leadP;
	double q0 = leadR / leadP;

	// 6. Set up and solve the quadratic indicial equation: r^2 + (p0 - 1)r + q0 = 0
	complex<double> b(p0 - 1.0, 0.0);
	complex<double> c(q0, 0.0);
    
	complex<double> discriminant = b * b - 4.0 * c;
	complex<double> sqrt_disc = std::sqrt(discriminant);

	r1 = (-b + sqrt_disc) / 2.0;
	r2 = (-b - sqrt_disc) / 2.0;     

	// Ensure r1 has the larger real part
	if (r2.real() > r1.real()) 
	{
		swap(r1, r2);
	}
	// categorize the roots
	complex_roots = (abs(imag(sqrt_disc)) > 1e-10 );
	repeated_roots = (abs(real(sqrt_disc)) < 1e-10 && abs(imag(sqrt_disc)) < 1e-10);
	int N = static_cast<int>(real(r1 - r2)); 
	integer_diff = (N > 0 && abs((r1 - r2) - (double)N) < 1e-7);
	noninteger_diff = ( N > 0 && abs((r1 - r2) - (double)N) > 1e-7 );
	//repeated_roots = ( real(r1) - real(r2) < 1e-8);
}
// Solves the indicial equation: P_0*r*(r-1) + Q_0*r + R_0 = 0
std::pair<complex<double>, complex<double>> SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::solve_indicial_equation_inpair() {


	// 2. Identify valuations (lowest non-zero powers)
	int vP = P.valuation();
	int vQ = Q.valuation();
	int vR = R.valuation();

	// 3. Check Regular Singularity Condition
	// If P(x0) != 0, it's an ordinary point, not a singular point. 
	// For it to be a regular singular point: vQ >= vP - 1 and vR >= vP - 2
	if (vP == 0 || vQ < vP - 1 || vR < vP - 2) 
	{
		std::cerr << "Error: Point x0 = " << x0 << " is not a regular singular point!\n";
	}

	// 4. Extract dominant structural coefficients
	double leadP = P.coeffs[vP];
	double leadQ = (vQ == vP - 1) ? Q.coeffs[vQ] : 0.0;
	double leadR = (vR == vP - 2) ? R.coeffs[vR] : 0.0;

	// 5. Compute limits p0 and q0 exactly
	double p0 = leadQ / leadP;
	double q0 = leadR / leadP;

	// 6. Set up and solve the quadratic indicial equation: r^2 + (p0 - 1)r + q0 = 0
	complex<double> b(p0 - 1.0, 0.0);
	complex<double> c(q0, 0.0);
    
	complex<double> discriminant = b * b - 4.0 * c;
	complex<double> sqrt_disc = std::sqrt(discriminant);

	complex<double> r1 = (-b + sqrt_disc) / 2.0;
	complex<double> r2 = (-b - sqrt_disc) / 2.0;     

	// Ensure r1 has the larger real part
	if (r2.real() > r1.real()) 
	{
		swap(r1, r2);
	}
	return {r1, r2};
}

// Evaluates fundamental solution y_basis at target x
complex<double> SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::evaluate_series_y1(complex<double> x) 
{
		complex<double> t = x - x0;
		if (abs(t) <= 0.0 && std::floor(real(r1)) != real(r1)) 
		{
			return 0.0; // Guard against negative fractional bases
		}
		complex<double> sum(0.0,0.0);
		for (int n = max_terms - 1; n >= 0; --n) 
		{
			sum = sum * t + coefficients_root1[n];
		}
		return sum * std::pow(t, r1);
}
// This function is necessary since the formula to compute y2(x) depends on the type of roots.
complex<double> SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::evaluate_series_y2(complex<double> x) 
{
		complex<double> t = x - x0;
		if (abs(t) <= 0.0 && std::floor(real(r2)) != real(r2)) 
		{
			return 0.0; // Guard against negative fractional bases
		}
		complex<double> y1(0.0,0.0), y2(0.0, 0.0);
		for (int n = max_terms - 1; n >= 0; --n) 
		{
			y1 = y1 * t + coefficients_root1[n];
		}
		y1 = y1 * std::pow(t, r1);

		complex<double> sum(0.0,0.0);
		if(noninteger_diff  || complex_roots)
		{
			for (int n = max_terms - 1; n >= 0; --n) 
			{
				sum = sum * t + coefficients_root2[n];
			}
			y2 = sum * std::pow(t, r2);
		}
		if(integer_diff && !complex_roots)
		{
			for (int n = max_terms - 1; n >= 0; --n) 
			{
				sum = sum * t + coefficients_root2[n];
			}
			y2 = (C_log * log(t) * y1) + ( sum * std::pow(t, r2) );
		}
		if(repeated_roots)
		{
			for (int n = max_terms - 1; n >= 0; --n) 
			{
				sum = sum * t + coefficients_root2[n];
			}
			y2 = (log(t) * y1) + ( sum * std::pow(t, r2) );
		}
		
		return y2;
}

// Create evaluate_series_y2(const vector<complex<double>>& a, vector<complex<double>>& y1, complex<double> r2, complex<double> x) 
// Evaluates derivative of basis solution
complex<double> SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::evaluate_series_derivative(const vector<complex<double>>& a, complex<double> r_val, complex<double> x) 
{
		complex<double> t = x - x0;
		complex<double> sum(0.0,0.0);
		for (int n = max_terms - 1; n >= 0; --n) 
		{
			// Power rule tracking d/dt [a_n * t^(n+r)] = (n+r)*a_n * t^(n+r-1)
			sum = sum * t + a[n] * (complex<double>(n) + r_val);
		}
		return sum * std::pow(t, r_val - complex<double>(1));
}

// Solves Initial Value Problem given y(x_init) and y'(x_init)
void SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::solve_ivp(double x_init, double y_init, double dy_init, const vector<complex<double>>& test_points) 
{

	// Basis evaluations at configuration step
	complex<double> y1_0 = evaluate_series_y1(x_init);
	complex<double> y2_0 = evaluate_series_y2(x_init);
	complex<double> dy1_0 = evaluate_series_derivative(coefficients_root1, r1, x_init);
	complex<double> dy2_0 = evaluate_series_derivative(coefficients_root2, r2, x_init);

	// Solve Cramer's system: 
	// [ y1_0  y2_0 ] [ C1 ]  =  [ y_init  ]
	// [dy1_0 dy2_0 ] [ C2 ] = [ dy_init ]
	complex<double> det = y1_0 * dy2_0 - y2_0 * dy1_0;
	if (std::abs(det) < 1e-9) 
	{
		std::cerr << "Error: Wronskian determinant too small at initial point." << std::endl;
		return;
	}

	complex<double> C1 = (y_init * dy2_0 - y2_0 * dy_init) / det;
	complex<double> C2 = (y1_0 * dy_init - y_init * dy1_0) / det;

	cout << "\n--- IVP Solution Results ---" << endl;
 	cout << "Initial Conditions:\nt0 = " << x_init << ",\t y(t0) = " << y_init << ", \t y'(t0) = " << dy_init << endl;
	cout << "Linear Coefficients: C1 = " << C1 << ", C2 = " << C2 << "\n\n";
	cout << "x\t\t\ty(x)" << endl;
	cout << "-----------------------------------" << endl;
        
	for (complex<double> x : test_points) 
	{
		complex<double> y_val = C1 * evaluate_series_y1(x) + C2 * evaluate_series_y2(x);
		cout << x << "\t\t\t" << y_val << endl;
	}
}

void SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::printCoefficients() const 
{
	// Output calculated power series coefficients to the console
	cout << "\n\nComputed Power Series Coefficients\n\n";
	cout << "y1(x) coefficients (r = " << r1 << "):" << endl;
	for (int i = 0; i < int(coefficients_root1.size()); ++i) 
	{
		cout << "a_{" << i << "} = " << std::setw(1) << coefficients_root1[i] << "\n";
	}

	cout << "\n\ny2(x) coefficients (r = " << r2 << "):" << endl;
	if(integer_diff)
	{	
		cout << "Logarithmic constant C = " << C_log << endl;
	}
	for (int i = 0; i < int(coefficients_root1.size()); ++i) 
	{
		cout << "b_{" << i << "} = " << std::setw(1) << coefficients_root2[i] << "\n";
	}

	
}

void SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver::printSolution() const 
{
	if(Frobenius)
	{
		if(integer_diff )
		{
			cout << "\nSeries solution: \n"<< endl;
			cout << "y1(x) = x^" << r1 << " ( ";
			bool first = true;
			int n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(real(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << std::abs((coefficients_root1[i]));
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " " ;
			}

			cout << ")\ny2(x) = " << C_log << " * ln(x) y1(x) + x^" << r2 <<" ( ";

			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root2.size(); ++i) 
			{
				if (std::abs(real(coefficients_root2[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root2[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root2[i]) < 0) 
				{
					cout << " - ";
				}
				cout << std::abs((coefficients_root2[i]));
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root2.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root2.size()) )
			{
				cout << " " ;
			}
			cout << ")";
		}
		if(noninteger_diff)
		{
			cout << "\nSeries solution: \n"<< endl;
			cout << "y1(x) = x^" << r1 <<" ( ";
			bool first = true;
			int n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(real(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << std::abs((coefficients_root1[i]));
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " " ;
			}
			
			cout << ")\ny2(x) = x^" << r2 <<" ( ";

			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root2.size(); ++i) 
			{
				if (std::abs(real(coefficients_root2[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root2[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root2[i]) < 0) 
				{
					cout << " - ";
				}
				cout << std::abs((coefficients_root2[i]));
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root2.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root2.size()) )
			{
				cout << " " ;
			}
			cout <<" )";
		}
		if(repeated_roots)
		{
			cout << "\nSeries solution: \n"<< endl;
			cout << "y1(x) = x^"<< r1 << " ( ";
			bool first = true;
			int n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(real(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << std::abs((coefficients_root1[i]));
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " " ;
			}
			cout << ")\ny2(x) = y1(x) ln(x) + x^" << r1 <<" ( ";

			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root2.size(); ++i) 
			{
				if (std::abs(real(coefficients_root2[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root2[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root2[i]) < 0) 
				{
					cout << " - ";
				}
				cout << std::abs((coefficients_root2[i]));
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root2.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root2.size()) )
			{
				cout << " " ;
			}
			cout <<")";
		}
		if(complex_roots)
		{
			cout << "\nSeries solution: \n"<< endl;
			cout << "y1(x) = x^" << r1 <<" ( ";
			bool first = true;
			int n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(real(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root1[i]) < 0) 
				{
					cout << " + ";
				}
				cout << coefficients_root1[i];
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " " ;
			}
			
			cout << ")\ny2(x) = x^" << r2 <<" ( ";

			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root2.size(); ++i) 
			{
				if (std::abs(real(coefficients_root2[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root2[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root2[i]) < 0) 
				{
					cout << " + ";
				}
				cout << coefficients_root2[i];
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root2.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root2.size()) )
			{
				cout << " " ;
			}
			cout <<" )";

			// to show real-valued solutions

			cout << "\n\nReal-valued series solution: \n"<< endl;
			if (abs(x0)==0)
			{
				cout << "y1(x) = x^" << r1.real() <<"[ cos( " << r1.imag() << " ln (x) ) * ( ";
			}
			else if (abs(x0) != 0)
			{
				cout << "y1(x) = x^" << r1.real() <<"[ cos( " << r1.imag() << " ln (x-" << x0 << " ) ) * ( " ;
			}
			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(real(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << real(coefficients_root1[i]);
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " ) " ;
			}

			if (abs(x0)==0)
			{
				cout << " - sin( " << r1.imag() << " ln (x) ) * ( " ;
			}
			else if (abs(x0) != 0)
			{
				cout << " - sin( " << r1.imag() << " ln (x-" << x0 << " ) ) * ( " ;
			}
			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(imag(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && imag(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (imag(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << imag(coefficients_root1[i]);
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " ) " ;
			}

			cout << "]\ny2(x) = x^" << r1.real() <<" [ ";

			if (abs(x0)==0)
			{
				cout << "sin( " << r1.imag() << " ln (x) ) * ( ";
			}
			else if (abs(x0) != 0)
			{
				cout << "sin( " << r1.imag() << " ln (x-" << x0 << " ) ) * ( " ;
			}
			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(real(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && real(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (real(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << real(coefficients_root1[i]);
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " ) " ;
			}

			if (abs(x0)==0)
			{
				cout << " + cos( " << r1.imag() << " ln (x) ) * ( " ;
			}
			else if (abs(x0) != 0)
			{
				cout << " + cos( " << r1.imag() << " ln (x-" << x0 << " ) ) * ( " ;
			}
			first = true;
			n_terms = 0;
			for (size_t i = 0; i < coefficients_root1.size(); ++i) 
			{
				if (std::abs(imag(coefficients_root1[i])) < 1e-9) 
				{
					continue;
				}
				if (!first && imag(coefficients_root1[i]) > 0) 
				{
					cout << " + ";
				}
				if (imag(coefficients_root1[i]) < 0) 
				{
					cout << " - ";
				}
				cout << imag(coefficients_root1[i]);
				if (i > 0) 
				{
					if (abs(x0)==0)
					{
						cout << "*x";
					}
					else if (abs(x0) != 0)
					{
						cout << "*(x - " << x0 << ")";
					}
					if (i > 1) 
					{
						cout << "^" << i;
					}
				}
			first = false;
			n_terms += 1;
			}
			
			if (n_terms >= int(coefficients_root1.size()))
			{
				cout << " + ... \n";
			}
			else if (n_terms < int(coefficients_root1.size()) )
			{
				cout << " ) ]" ;
			}
		}
	}
	else if(EulerCauchy)
	{
		if(integer_diff || noninteger_diff)
		{
			if (abs(x0)==0)
			{
				cout << "\nGeneral solution: \n"<< endl;
				cout << "y(x) = C1 * x^" << r1 << " + C2 * x^" << r2 ;
			}
			else if (abs(x0)!=0)
			{
				cout << "\nGeneral solution: \n"<< endl;
				cout << "y(x) = C1 * (x - "<< x0 << ")" <<"^" << r1 << " + C2 * (x - " << x0 << ")^" << r2 ;
			}
		}
		
		if(repeated_roots )
		{
			if (abs(x0)==0)
			{
				cout << "\nGeneral solution: \n"<< endl;
				cout << "y(x) = ( C1 + C2 ln(x) ) x^" << r1  ;
			}
			else if (abs(x0)!=0)
			{
				cout << "\nGeneral solution: \n"<< endl;
				cout << "y(x) = ( C1 + C2 ln(x - " << x0 << ") ) (x - " << x0 << ")^" << r1  ;
			}

		}
		if(complex_roots )
		{
			if (abs(x0)==0)
			{
				cout << "\nGeneral solution: \n"<< endl;
				cout << "y(x) = [ C1 *  x^{" << r1.real() << "} cos ( "<< abs(r1.imag()) << " ln(x) ) ] + [ C2 *  x^{" << r1.real() << "} sin ( "<< abs(r1.imag()) << " ln(x) ) ] ";
			}
			else if (abs(x0)!=0)
			{
				cout << "\nGeneral solution: \n"<< endl;
				cout << "y(x) = [ C1 *  (x - " << x0 << ")^{" << r1.real() << "} cos ( "<< abs(r1.imag()) << " ln(x - " << x0 << ") ) ] + [ C2 *  (x - " << x0 << ")^{" << r1.real() << "} sin ( "<< abs(r1.imag()) << " ln(x - " << x0 << ") ) ] ";
			}
			
		}
	}
}

void secondorderlineardiffeq_derivativesvalueatx0(const Symbolic &diffeq, const Symbolic &y, const Symbolic &x, double x0, Symbolic y0, Symbolic dy0)
{
	cout << "\nThe differential equation :\n" <<  diffeq << endl;
	cout << "\nx0 = " << x0<< ", y(" << x0 << ") = " << y0 << ", y'(" << x0 << ") = " << dy0 << endl;
	
	Equations Eq_dy2_x0 = solve(diffeq,df(y[x],x,2));
	
	Equations rules = (  df(y[x],x,1) == dy0, y[x] == y0, x == x0);
	Symbolic dy2_x0_final = Eq_dy2_x0.front().rhs.subst_all(rules);

	cout << "\ny'' = " <<  Eq_dy2_x0.front().rhs << endl;
	cout << "\ny''(x0) = " << dy2_x0_final  << endl;	

	Symbolic dy3 = df(diffeq,x) ;
	Symbolic dy4 = df(diffeq,x,2);
	cout << "\nDifferentiate the differential equation with respect to " << x <<" :\n" <<  dy3 << endl;
	
	Equations Eq_dy3_x0 = solve(dy3,df(y[x],x,3));
	
	Equations rules2 = ( df(y[x],x,2) == dy2_x0_final, df(y[x],x,1) == dy0, y == y0, x == x0);
	Symbolic dy3_x0_final = Eq_dy3_x0.front().rhs.subst_all(rules2);

	cout << "\ny''' = " <<  Eq_dy3_x0.front().rhs << endl;
	cout << "\ny'''(x0) = " << dy3_x0_final  << endl;
	
	cout << "\nDifferentiate the differential equation again with respect to " << x << " :\n" <<  dy4 << endl;

	Equations Eq_dy4_x0 = solve(dy4,df(y[x],x,4));
	
	Equations rules3 = ( df(y[x],x,3) == dy3_x0_final , df(y[x],x,2) == dy2_x0_final, df(y[x],x,1) == dy0, y == y0, x == x0);
	Symbolic dy4_x0_final = Eq_dy4_x0.front().rhs.subst_all(rules3);

	cout << "\ny''' = " <<  Eq_dy4_x0.front().rhs << endl;
	cout << "\ny'''(x0) = " << dy4_x0_final  << endl;
}

#endif
#endif