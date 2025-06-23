/*
   
*/

#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, mu, C("C"), d1;
 
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d;
		eq = (a*x*y + d).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = exp(integrate(-ap*x,x));
		dsol = (integrate(mu*dp,x))/(mu) + (C)/(mu);
		} catch(const SymbolicError &se) {}
		}
	}
	return dsol;
}

#endif
#endif