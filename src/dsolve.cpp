/*
   
*/

#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, mu, C("C");
 
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d;
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
	}
	return dsol;
}

#endif
#endif