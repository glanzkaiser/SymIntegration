/*
   
*/

#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICSMLE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICSMLE_DEFINE

Symbolic mle(const Symbolic &fx, const Symbolic &x, const Symbolic &t, const Symbolic &s, int n)
{
	Symbolic mlesol, sum;
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d;
		//int n1;
		// Case 1 : a*θ/(b*x^(c*θ+d))
		eq = ((a*t)*(b*x^(-c*t-d))).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c),dp = rhs(*i, d);

		for (int i = 1; i < n+1 ; ++i)
		{
			sum += b*ln(x[i]);
		}
		mlesol = divisions(n,(b*sum));

		//if(df(rhs(*i, a), x) == 0) 
		//{
		//	return mlesol;
		//}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : θ/(x^(θ+1)) works amazingly
		eq = ((t)*(x^(-t-1))).match(fx, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		
		for (int i = 1; i < n+1 ; ++i)
		{
			sum += ln(x[i]);
		}

		mlesol = divisions(n,sum);

		} catch(const SymbolicError &se) {}
		}
		
		// Case 3 : (e^{-λ} * λ^{x} ) / x!
		eq = ((exp(-t))*(pow(t,x))/(factorialsym(x))).match(fx, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		
		for (int i = 1; i < n+1 ; ++i)
		{
			sum += x[i];
		}
		if(factorialsym(x) == 1) 
		{
			mlesol = divisions(sum,n);
		}
		} catch(const SymbolicError &se) {}
		}
	}
	return mlesol;
}

#endif
#endif