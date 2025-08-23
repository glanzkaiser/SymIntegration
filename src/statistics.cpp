/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE

// Factorial and combinations
Symbolic factorial1(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	Symbolic result = 1;
	for (int i = 2; i <= n; ++i) 
	{
		result *= i;
	}
	return result;
}

Symbolic combinations1(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factorial1(n) / (factorial1(r) * factorial1(n - r));
}

Symbolic binomialpmf(int x, int n, double p)
{
 	Symbolic b = combinations1(n,x)*pow(p,x)*pow(1-p,n-x);
	
	return b;
}

Symbolic binomialcdf(int x, int n, double p)
{
 	Symbolic b;
	for (int i=0; i<=x;i++)
	{
		b += combinations1(n,i)*pow(p,i)*pow(1-p,n-i);
	}	
	return b;
}

Symbolic binomialmean(int x, int n, double p)
{
 	Symbolic b = n*p;

	return b;
}

Symbolic binomialvar(int x, int n, double p)
{
 	Symbolic b = n*p*(1-p);

	return b;
}

Symbolic binomialmgf(int x, int n, double p)
{
 	Symbolic t("t");
	Symbolic b = (1-p)+p*exp(t);

	return b;
}

Symbolic negativebinomialpmf(int x, int k, double p)
{
 	Symbolic b = combinations1(x-1,k-1)*pow(p,k)*pow(1-p,x-k);
	
	return b;
}

Symbolic negativebinomialmean(int x, int k, double p)
{
 	Symbolic b = k/p;

	return b;
}

Symbolic negativebinomialvar(int x, int k, double p)
{
 	Symbolic b = k*(1-p)/(p*p);

	return b;
}

Symbolic negativebinomialmgf(const Symbolic &x, const Symbolic &k, const Symbolic &p)
{
 	Symbolic t("t");

	Symbolic b =pow(p,k)*pow(1-((1-p)*exp(t)),-k);
	
	return b;
}

Symbolic geometricpmf(int x, double p)
{
 	Symbolic b = p*(pow(1-p,x-1));
	
	return b;
}

Symbolic geometricmean(int x, double p)
{
 	Symbolic b = 1/p;

	return b;
}

Symbolic geometricvar(int x, double p)
{
 	Symbolic b = (1-p)/(p*p);

	return b;
}

Symbolic geometricmgf(const Symbolic &x, const Symbolic &p)
{
 	Symbolic t("t");

	Symbolic b =p*pow(1-((1-p)*exp(t)),Symbolic(-1));
	
	return b;
}

Symbolic poissonpmf(int x, int k)
{
 	Symbolic b = exp(-k)*pow(k,x)/(factorial1(x));
	
	return b;
}

Symbolic poissoncdf(int x, int k)
{
 	Symbolic b;
	for (int i=0; i<=x;i++)
	{
		b += exp(-k)*pow(k,i)/(factorial1(i));
	}	
	return b;
}

Symbolic poissonmean(int x, int k)
{
 	Symbolic b = k;

	return b;
}

Symbolic poissonvar(int x, int k)
{
 	Symbolic b = k;

	return b;
}

Symbolic poissonmgf(const Symbolic &x, const Symbolic &k)
{
 	Symbolic t("t");

	Symbolic b =exp(k*(exp(t)-1));
	
	return b;
}

#endif
#endif