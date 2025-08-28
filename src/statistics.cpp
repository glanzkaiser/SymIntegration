/*
   
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define π 3.1415926535897f

double divisionint(double x, double y)
{
	return x/y;
}

double rising_pochhammer(double x, int n) 
{
	if (n < 0) 
	{
	// Handle negative n if needed, often defined using Gamma function
	// (x)_n = Gamma(x+n) / Gamma(x)
	// For now, let's assume n >= 0 for simplicity
		return 0.0; // Or throw an exception
	}
	if (n == 0) 
	{
		return 1.0;
	}
	double result = 1.0;
	for (int i = 0; i < n; ++i) 
	{
		result *= (x + i);
	}
	return result;
}

double falling_pochhammer(double x, int n) 
{
	if (n < 0) 
	{
	// Handle negative n if needed
		return 0.0; // Or throw an exception
	}
	if (n == 0) 
	{
		return 1.0;
	}
	double result = 1.0;
	for (int i = 0; i < n; ++i) 
	{
		result *= (x - i);
        }
	return result;
}

int fibonacciseries(int n)
{
	int t1 = 0, t2 = 1, nextTerm = 0;

    	cout << "Fibonacci Series: ";

    	for (int i = 1; i <= n; ++i) {
        // Prints the first two terms.
        if(i == 1) {
            cout << t1 << ", ";
            continue;
     	   }
        if(i == 2) {
            cout << t2 << ", ";
            continue;
        }
        nextTerm = t1 + t2;
        t1 = t2;
        t2 = nextTerm;
        
        cout << nextTerm << ", " ;
    	}
	cout << "\nFibonacci n-th number: "<< endl; 
	return nextTerm;
}

Symbolic divisions(Symbolic x, Symbolic y)
{
	return x/y;
}
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

// Discrete Distributions
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

Symbolic hypergeometricpmf(int x, int N, int n, int k)
{
 	Symbolic b = combinations1(k,x)*divisions(combinations1(N-k,n-x),combinations1(N,n));
	
	return b;
}

Symbolic hypergeometricmean(int x, int N, int n, int k)
{
 	Symbolic b = divisionint(n*k,N);

	return b;
}

Symbolic hypergeometricvar(int x, int N, int n, int k)
{
 	Symbolic b = n*divisionint(N-n,N-1)*divisionint(k,N)*(1-divisionint(k,N));

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

// Continuous Distributions

Symbolic uniformpdf(double x, double a, double b)
{
 	Symbolic u = divisions(1,b-a);

	return u;
}

Symbolic uniformcdf(double x1, double a, double b)
{
	return divisions(x1,b-a) - divisions(a,b-a) ;
}

Symbolic uniformmgf(double x, double a, double b)
{
 	Symbolic t("t");
	Symbolic u = (exp(b*t) - exp(a*t))/((b-a)*t);

	return u;
}

Symbolic uniformmean(double x, double a, double b)
{
 	Symbolic u = 0.5*(a+b);

	return u;
}

Symbolic uniformvar(double x, double a, double b)
{
 	Symbolic u = divisions(1,12)*(b-a)*(b-a);

	return u;
}

Symbolic normalpdf(double x, double μ, double σ)
{
 	Symbolic b = divisions(1,sqrt(2*π)*σ)*exp(-0.5*divisions(x-μ,σ)*divisions(x-μ,σ));

	return evalf(b,1,1);
}

// Function to compute the standard normal CDF
double normalCDF(double x) {
    // M_SQRT1_2 is 1/sqrt(2)
    // Some compilers/libraries might require defining M_SQRT1_2 explicitly
    // const double M_SQRT1_2 = 0.70710678118; 
    return 0.5 * erfc(-x * M_SQRT1_2);
}

Symbolic normalcdf(double x1, double μ, double σ)
{
	double x = ((x1 - μ) / σ);
	
 	Symbolic N = 0.5 * erfc(-x * M_SQRT1_2);

	return N;
}

Symbolic normalmgf(double x, double μ, double σ)
{
 	Symbolic t("t");
	Symbolic b = exp(μ*t + 0.5*σ*σ*t*t);

	return b;
}

Symbolic normalmean(double x, double μ, double σ)
{
 	Symbolic b = μ;

	return b;
}

Symbolic normalvar(double x, double μ, double σ)
{
 	Symbolic b = σ*σ;

	return b;
}

Symbolic gammapdf(double x, double α, double β)
{
 	Symbolic b = divisions(1,(Symbolic(β)^(Symbolic(α)))*tgamma(α))*(Symbolic(x)^(Symbolic(α-1)))*exp(-divisions(x,β));

	return evalf(b,1,1);
}

Symbolic gammacdf(double x1, double α, double β)
{	
	Symbolic x("x");
	Symbolic f = divisions(1,(Symbolic(β)^(Symbolic(α)))*tgamma(α))*(x^(Symbolic(α-1)))*exp(-divisions(x,β));

 	Symbolic b = integrate(f,x)[x==x1] - integrate(f,x)[x==0];

	return evalf(b,1,1);
}

Symbolic gammamgf(double x, double α, double β)
{
 	Symbolic t("t");
	Symbolic b = ((1-β*t)^(Symbolic(-α)));

	return b;
}

Symbolic gammamean(double x, double α, double β)
{
 	Symbolic b = α*β;

	return b;
}

Symbolic gammavar(double x, double α, double β)
{
 	Symbolic b = α*β*β;

	return b;
}

Symbolic exponentialpdf(double x, double λ)
{
 	Symbolic b = λ*exp(-λ*x);

	return evalf(b,1,1);
}

Symbolic exponentialcdf(double x1, double λ)
{	
	Symbolic x("x");
	Symbolic f =  λ*exp(-λ*x);

 	Symbolic b = 1 + integrate(f,x)[x==x1] ;

	return evalf(b,1,1);
}

Symbolic exponentialmgf(double x, double λ)
{
 	Symbolic t("t");
	Symbolic b = (1-(divisions(t,λ)))^(Symbolic(-1));

	return b;
}

Symbolic exponentialmean(double x, double λ)
{
 	Symbolic b = divisions(1,λ);

	return b;
}

Symbolic exponentialvar(double x, double λ)
{
 	Symbolic b = divisions(1,λ*λ);

	return b;
}

Symbolic betapdf(double x, double  α, double β)
{
 	Symbolic b = ( tgamma(α+β)/(tgamma(α)*tgamma(β)) ) * (x^(Symbolic(α-1))) * ((1-x)^(Symbolic(β-1)));

	return evalf(b,1,1);
}

Symbolic betacdf(double x1, double α, double β)
{	
	Symbolic x("x");
	int alpha  = α;
	int beta = β;
	Symbolic f = ( tgamma(α+β)/(tgamma(α)*tgamma(β)) ) * (pow(x,Symbolic(alpha-1))) * (pow(1-x,(Symbolic(beta-1))));

 	Symbolic b = integrate(f,x)[x==x1] ;

	return evalf(b,1,1);
}

Symbolic betamgf(double x, double α, double β)
{
 	Symbolic t("t");
	Symbolic b = hypergeometric_1F1(α, α+β,t,5);

	return b;
}

Symbolic betamean(double x, double  α, double β)
{
 	Symbolic b = divisions(α,α+β);

	return b;
}

Symbolic betavar(double x, double  α, double β)
{
 	Symbolic b = divisions(α*β,(α+β+1)*(α+β)*(α+β));

	return b;
}



#endif
#endif