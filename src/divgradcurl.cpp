/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DIVGRADCURL_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DIVGRADCURL_DEFINE

Symbolic div(const Symbolic &F1, const Symbolic &F2, const Symbolic &F3, const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol, i("i"), j("j"), k("k");
 	
	if(F1 != 0 || F2 != 0 || F3 != 0)
 	{
		sol = df(F1,x) + df(F2,y) + df(F3,z);
	}
	return sol;
}

Symbolic div(const Symbolic &F, const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol, i("i"), j("j"), k("k"), î("î"), ĵ("ĵ"), k̂("k̂");
 	
	Symbolic F1 = F.coeff(i,1);
	Symbolic F2 = F.coeff(j,1);
	Symbolic F3 = F.coeff(k,1);
	Symbolic F1_hat = F.coeff(î,1);
	Symbolic F2_hat = F.coeff(ĵ,1);
	Symbolic F3_hat = F.coeff(k̂,1);
	
	if(F1 != 0 || F2 != 0 || F3 != 0)
 	{
		sol = df(F1,x) + df(F2,y) + df(F3,z);
	}
	else if(F1_hat != 0 || F2_hat != 0 || F3_hat != 0)
 	{
		sol = df(F1_hat,x) + df(F2_hat,y) + df(F3_hat,z);
	}
	return sol;
}

Symbolic grad(const Symbolic &F, const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol;
 	
	if(F != 0)
 	{
		sol = (df(F,x), df(F,y),  df(F,z));
		sol = sol.transpose();
	}
	return sol;
}

Symbolic grad(const Symbolic &F, const Symbolic &x, const Symbolic &y)
{
	Symbolic sol;
 	
	if(F != 0)
 	{
		sol = (df(F,x), df(F,y));
		sol = sol.transpose();
	}
	return sol;
}

Symbolic curl(const Symbolic &F1, const Symbolic &F2, const Symbolic &F3, const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol, i("i"), j("j"), k("k");
 	
	if(F1 != 0 || F2 != 0 || F3 != 0)
 	{
		sol = (df(F3,y) - df(F2,z), -(df(F3,x) - df(F1,z)),  df(F2,x) - df(F1,y));
		sol = sol.transpose();
	}
	return sol;
}

Symbolic curl(const Symbolic &F, const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol, i("i"), j("j"), k("k"), î("î"), ĵ("ĵ"), k̂("k̂");
 	
	Symbolic F1 = F.coeff(i,1);
	Symbolic F2 = F.coeff(j,1);
	Symbolic F3 = F.coeff(k,1);
	Symbolic F1_hat = F.coeff(î,1);
	Symbolic F2_hat = F.coeff(ĵ,1);
	Symbolic F3_hat = F.coeff(k̂,1);
	
	if(F1 != 0 || F2 != 0 || F3 != 0)
 	{
		sol = (df(F3,y) - df(F2,z), -(df(F3,x) - df(F1,z)),  df(F2,x) - df(F1,y));
		sol = sol.transpose();
	}
	else if(F1_hat != 0 || F2_hat != 0 || F3_hat != 0)
 	{
		sol = (df(F3_hat,y) - df(F2_hat,z), -(df(F3_hat,x) - df(F1_hat,z)),  df(F2_hat,x) - df(F1_hat,y));
		sol = sol.transpose();
	}
	return sol;
}

#endif
#endif