/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_INFINITESERIES_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_INFINITESERIES_DEFINE

Symbolic taylorseries(const Symbolic &f, const Symbolic &x, double x0, int N)
{
	Symbolic n("n"), sol, F;
	F=f;
 	for (int i = 0; i<N; ++i)
 	{
		sol += F[x==x0]/(factorial(i)) * pow(x-x0,Symbolic(i));
		F = df(F,x);
	}
	return sol;
}

#endif
#endif