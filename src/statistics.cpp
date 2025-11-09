/*
#include <boost/math/distributions/beta.hpp> 
#include <boost/math/distributions/gamma.hpp> 
#include <boost/math/distributions/students_t.hpp>    
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)
#include <boost/math/distributions/laplace.hpp> //Faster computing of laplace cdf and pdf
#include <boost/math/distributions/logistic.hpp> //Faster computing of logistic cdf and pdf

#include <vector>
#include <map>
#include <algorithm> // For std::max_element,  std::sort
#include <numeric> // For std::accumulate

#include <random> // For random number generation
#include <chrono>

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define π 3.141592653589793238462643383279502884f

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
double factoriald(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	double result = 1;
	for (int i = 2; i <= n; ++i) 
	{
		result *= i;
	}
	return result;
}

double combinationsd(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factoriald(n) / (factoriald(r) * factoriald(n - r));
}
// Generate random numbers
int randomnumberint(int a, int b, int n)
{
	// 1. Obtain a seed:
	// Using std::random_device for a non-deterministic seed if available,
	// otherwise falling back to a time-based seed.
	std::random_device rd;
	std::mt19937 generate(rd()); // Mersenne Twister engine seeded with rd

	// Alternatively, for a time-based seed (less truly random but often sufficient):
	// std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());

	// 2. Define a distribution:
	// For integers within a specific range (e.g., 1 to 100 inclusive)
	std::uniform_int_distribution<> distrib(a, b);
	
	for(int i=1; i<=n; i++)
	{
		cout << distrib(generate) << endl;
	}
	return distrib(generate);
}

double randomnumberreal(double a, double b, int n)
{
	// 1. Obtain a seed:
	// Using std::random_device for a non-deterministic seed if available,
	// otherwise falling back to a time-based seed.
	std::random_device rd;
	std::mt19937 generate(rd()); // Mersenne Twister engine seeded with rd

	// Alternatively, for a time-based seed (less truly random but often sufficient):
	// std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());

	// 2. Define a distribution:
	// For floating-point numbers within a specific range (e.g., 0.0 to 1.0)
	std::uniform_real_distribution<> real_distrib(a, b);

	for(int i=1; i<n; i++)
	{
		cout << real_distrib(generate) << endl;
	}
	return real_distrib(generate);
}

double randomnumbergamma(double alpha, double beta, int n)
{
	// 1. Obtain a seed:
	// otherwise falling back to a time-based seed.
	// Seed it with a time-based value for more randomness across runs.
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generate(seed); 

	// 2. Define the gamma distribution parameters
	// alpha (shape parameter) and beta (scale parameter)
	std::gamma_distribution<> distribution(alpha, beta);

	for(int i=1; i<n; i++)
	{
		cout << distribution(generate) << endl;
	}
	return distribution(generate);
}

// Generating random numbers in C++ can be achieved using either the older C-style rand() and srand() functions or the more modern C++11 <random> library. 
// The <random> library is generally preferred for its better statistical properties and more flexible control over distributions.
// This approach provides more robust and statistically sound random number generation.
// Using std::random_device for a non-deterministic seed (if available and entropy > 0)
// or std::chrono::system_clock::now().time_since_epoch().count() for a time-based seed.

std::vector<double> vrandn_bernoulli(double p, int n)
{
	// 1. Obtain a seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::vector<double> vec;
	// 2. Create a Bernoulli distribution
	// The constructor takes the probability 'p' of generating 'true'.
	std::bernoulli_distribution distribution(p);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(distribution(generator))); 
	}
	return vec;
}

std::vector<double> vrandn_binomial(double p, int n)
{
	// 1. Obtain a seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::vector<double> vec;

	// 2. Create a binomial distribution object
	std::binomial_distribution<> dist_binomial(n, p);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(dist_binomial(generator))); 
	}
	return vec;
}

std::vector<double> vrandn_normal(double mu, double sigma, int n)
{
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<double> vec;
 	std::normal_distribution<double> distribution(mu, sigma);
	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(distribution(generator))); 
	}
	return vec;
}

std::vector<double> vrandn_exponential(double lambda, int n)
{
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	// We seed it with the current time for better randomness.
	std::mt19937 generator(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	std::vector<double> vec;
 	// 2. Create an exponential distribution object
	// The constructor takes the 'lambda' (rate) parameter.
	// The mean of an exponential distribution is 1/lambda.
	std::exponential_distribution<double> distribution(lambda);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(distribution(generator))); 
	}
	return vec;
}

std::vector<double> vrandn_chisquared(double nu, int n)
{
	// 1. Obtain a seed:
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);

	std::vector<double> vec;
	// 2. Define the Chi-squared distribution
	// The constructor takes the degrees of freedom (n) as a parameter.
	std::chi_squared_distribution<float> chi_squared_dist(nu);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(chi_squared_dist(generator))); 
	}
	return vec;
}

// Instantiate an std::fisher_f_distribution: This distribution represents the F-distribution. 
// Its constructor takes two parameters: the numerator degrees of freedom (v1) and the denominator degrees of freedom (v2). Both must be positive.

std::vector<double> vrandn_fdist(double numerator_df, double denominator_df, int n)
{
	// 1. Obtain a seed:
	std::random_device rd; // Obtain a non-deterministic seed
	std::mt19937 gen(rd()); // Standard Mersenne Twister engine seeded with rd()

	std::vector<double> vec;
 	// 2. Instantiate the F-distribution
	std::fisher_f_distribution<> f_dist(numerator_df, denominator_df);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(f_dist(gen))); 
	}
	return vec;
}

std::vector<double> vrandn_tdist(double nu, int n)
{
	// 1. Obtain a seed:
	std::random_device rd; // Obtain a non-deterministic seed
	std::mt19937 gen(rd()); // Standard Mersenne Twister engine seeded with rd()

	std::vector<double> vec;
 	// 2. Define the Student's t-distribution.
	// The constructor takes the degrees of freedom (nu).
	std::student_t_distribution<> t_dist(nu);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(t_dist(gen))); 
	}
	return vec;
}

// Generating Erlang-distributed random numbers in C++ involves using the <random> library, 
// which provides facilities for generating pseudo-random numbers with various distributions. 
// Since the Erlang distribution is a special case of the Gamma distribution where the shape parameter (k or alpha) is an integer, 
// you can leverage the std::gamma_distribution to achieve this. 

std::vector<double> vrandn_erlang(double k, double lambda, int n)
{
	// 1. Obtain a seed:
	std::mt19937 engine(std::chrono::system_clock::now().time_since_epoch().count());

	std::vector<double> vec;
	const double beta = 1.0 / lambda; // Scale parameter for std::gamma_distribution
 	// 2. Create a gamma_distribution object
	// For Erlang, the shape parameter (alpha) is 'k' and the scale parameter (beta) is '1/lambda'.
	std::gamma_distribution<double> erlang_dist(k, beta);

	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(erlang_dist(engine))); 
	}
	return vec;
}

std::vector<double> vrandn_gamma(double alpha, double beta, int n)
{
	// 1. Obtain a seed:
	// otherwise falling back to a time-based seed.
	// Seed it with a time-based value for more randomness across runs.
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::mt19937 generate(seed); 

	// 2. Define the gamma distribution parameters
	// alpha (shape parameter) and beta (scale parameter)
	std::gamma_distribution<> distribution(alpha, beta);
	std::vector<double> vec;
 
	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(distribution(generate))); 
	}
	return vec;
}

// The C++ Standard Library does not directly offer a std::beta_distribution
// a Beta distribution can be simulated using two independently generated Gamma distributions.

double generateBeta(double alpha, double beta, std::mt19937& generator) {
	// Create two gamma distributions with scale parameter 1
	std::gamma_distribution<double> gamma_alpha(alpha, 1.0);
	std::gamma_distribution<double> gamma_beta(beta, 1.0);

	// Generate random numbers from each gamma distribution
	double x = gamma_alpha(generator);
	double y = gamma_beta(generator);

	// Calculate the beta-distributed random number
	return x / (x + y);
}

std::vector<double> vrandn_beta(double alpha, double beta, int n)
{
	// 1. Obtain a seed:
	std::random_device rd;
	std::mt19937 generate(rd()); // Mersenne Twister engine seeded with rd

	std::vector<double> vec;
 	
	for(int i=1; i<n; i++)
	{
		vec.push_back(static_cast<double>(generateBeta(alpha, beta, generate))); 
	}
	return vec;
}



// Discrete Distributions
double binomialpmf(int x, int n, double p)
{
 	Symbolic b = combinationsd(n,x)*pow(p,x)*pow(1-p,n-x);
	
	return b;
}

double binomialcdf(int x, int n, double p)
{
 	double b;
	for (int i=0; i<=x;i++)
	{
		b += combinationsd(n,i)*pow(p,i)*pow(1-p,n-i);
	}	
	return b;
}

double binomialmean(int x, int n, double p)
{
 	double mean = n*p;

	return mean;
}

double binomialvar(int x, int n, double p)
{
 	double var = n*p*(1-p);

	return var;
}

Symbolic binomialmgf(int x, int n, double p)
{
 	Symbolic t("t");
	Symbolic b = (1-p)+p*exp(t);

	return b;
}

double negativebinomialpmf(int x, int k, double p)
{
 	double b = combinationsd(x-1,k-1)*pow(p,k)*pow(1-p,x-k);
	
	return b;
}

double negativebinomialmean(int x, int k, double p)
{
 	double mean = k/p;

	return mean;
}

double negativebinomialvar(int x, int k, double p)
{
 	double var = k*(1-p)/(p*p);

	return var;
}

Symbolic negativebinomialmgf(const Symbolic &x, const Symbolic &k, const Symbolic &p)
{
 	Symbolic t("t");

	Symbolic mgf =pow(p,k)*pow(1-((1-p)*exp(t)),-k);
	
	return mgf;
}

double geometricpmf(int x, double p)
{
 	double pmf = p*(pow(1-p,x-1));
	
	return pmf;
}

double geometricmean(int x, double p)
{
 	double mean = 1/p;

	return mean;
}

double geometricvar(int x, double p)
{
 	double var = (1-p)/(p*p);

	return var;
}

Symbolic geometricmgf(const Symbolic &x, const Symbolic &p)
{
 	Symbolic t("t");

	Symbolic mgf =p*pow(1-((1-p)*exp(t)),Symbolic(-1));
	
	return mgf;
}

double hypergeometricpmf(int x, int N, int n, int k)
{
 	double pmf = combinationsd(k,x)*divisiond(combinationsd(N-k,n-x),combinationsd(N,n));
	
	return pmf;
}

double hypergeometricmean(int x, int N, int n, int k)
{
 	double mean = divisionint(n*k,N);

	return mean;
}

double hypergeometricvar(int x, int N, int n, int k)
{
 	double var = n*divisionint(N-n,N-1)*divisionint(k,N)*(1-divisionint(k,N));

	return var;
}

double poissonpmf(int x, int k)
{
 	double pmf = exp(-k)*pow(k,x)/(factorial1(x));
	
	return pmf;
}

double poissoncdf(int x, int k)
{
 	double cdf;
	for (int i=0; i<=x;i++)
	{
		cdf += exp(-k)*pow(k,i)/(factoriald(i));
	}	
	return cdf;
}

double poissonmean(int x, int k)
{
 	double mean = k;

	return mean;
}

double poissonvar(int x, int k)
{
 	double var = k;

	return var;
}

Symbolic poissonmgf(const Symbolic &x, const Symbolic &k)
{
 	Symbolic t("t");

	Symbolic mgf =exp(k*(exp(t)-1));
	
	return mgf;
}

// Continuous Distributions

double uniformpdf(double x, double a, double b)
{
 	double pdf = divisiond(1,b-a);

	return pdf;
}

double uniformcdf(double x1, double a, double b)
{
	return divisiond(x1,b-a) - divisiond(a,b-a) ;
}

Symbolic uniformmgf(double x, double a, double b)
{
 	Symbolic t("t");
	Symbolic mgf = (exp(b*t) - exp(a*t))/((b-a)*t);

	return mgf;
}

double uniformmean(double x, double a, double b)
{
 	double mean = 0.5*(a+b);

	return mean;
}

double uniformvar(double x, double a, double b)
{
 	double var = divisions(1,12)*(b-a)*(b-a);

	return var;
}

Symbolic normalpdf(double x, double mu, double sigma)
{
 	Symbolic pdf = divisiond(1,sqrt(2*π)*sigma)*exp(-0.5*divisions(x-mu,sigma)*divisiond(x-mu,sigma));

	return evalf(pdf,1,1);
}

// Function to compute the standard normal CDF
double normalCDF(double x) {
    // M_SQRT1_2 is 1/sqrt(2)
    // Some compilers/libraries might require defining M_SQRT1_2 explicitly
    // const double M_SQRT1_2 = 0.70710678118; 
    return 0.5 * erfc(-x * M_SQRT1_2);
}

double normalcdf(double x1, double mu, double sigma)
{
	double x = ((x1 - mu) / sigma);
	
 	double cdf = 0.5 * erfc(-x * M_SQRT1_2);

	return cdf;
}

Symbolic normalmgf(double x, double mu, double sigma)
{
 	Symbolic t("t");
	Symbolic mgf = exp(mu*t + 0.5*sigma*sigma*t*t);

	return mgf;
}

double normalmean(double x, double mu, double sigma)
{
 	double mean = mu;

	return mean;
}

double normalvar(double x, double mu, double sigma)
{
 	double var = sigma*sigma;

	return var;
}

double zquantile(double cdf)
{
	double z_quantile = r8_normal_01_cdf_inverse(cdf);

	return z_quantile;	
}

double gammapdf(double x, double alpha, double beta) // 16 digits precise with boost::math::gamma_distribution<> gamma_dist(alpha,beta)
{
 	//Symbolic b = divisions(1,(Symbolic(beta)^(Symbolic(alpha)))*tgamma(alpha))*(Symbolic(x)^(Symbolic(alpha-1)))*exp(-divisions(x,beta));
	//return evalf(b,1,1);

	//boost::math::gamma_distribution<> gamma_dist(alpha,beta);
	//double x_value = x;
	//double pdf_value = boost::math::pdf(gamma_dist, x_value);

	if (x < 0 || alpha <= 0 || beta <= 0) 
	{
	// Handle invalid input: Gamma distribution is defined for x > 0, alpha > 0, beta > 0
	return 0.0; 
	}

	// Calculate the Gamma function of alpha
	double gamma_alpha = std::tgamma(alpha);

	// Compute the PDF
	double pdf_value = (std::pow(beta, -alpha) * std::pow(x, alpha - 1) * std::exp(-x/beta)) / gamma_alpha;

	return pdf_value;
}

double gammacdf(double x, double alpha, double beta)
{	
	//boost::math::gamma_distribution<> gamma_dist(alpha,beta);
	//double x_value = x;
	//double cdf_value = boost::math::cdf(gamma_dist, x_value);

	// 13 decimal digit precision with boost::math::gamma_distribution
	
	double cdf_value = lowergamma(alpha,x/beta)/(tgamma(alpha));

	return cdf_value;
}

Symbolic gammamgf(double x, double alpha, double beta)
{
 	Symbolic t("t");
	Symbolic mgf = ((1-beta*t)^(Symbolic(-alpha)));

	return mgf;
}

double gammamean(double x, double alpha, double beta)
{
 	double mean = alpha*beta;

	return mean;
}

double gammavar(double x, double alpha, double beta)
{
 	double var = alpha*beta*beta;

	return var;
}

Symbolic exponentialpdf(double x, double lambda)
{
 	Symbolic pdf = lambda*exp(-lambda*x);

	return evalf(pdf,1,1);
}

Symbolic exponentialcdf(double x1, double lambda)
{	
	Symbolic x("x");
	Symbolic f =  lambda*exp(-lambda*x);

 	Symbolic cdf = 1 + integrate(f,x)[x==x1] ;

	return evalf(cdf,1,1);
}

Symbolic exponentialmgf(double x, double lambda)
{
 	Symbolic t("t");
	Symbolic mgf = (1-(divisions(t,lambda)))^(Symbolic(-1));

	return mgf;
}

double exponentialmean(double x, double lambda)
{
 	double mean = divisions(1,lambda);

	return mean;
}

double exponentialvar(double x, double lambda)
{
 	double var = divisiond(1,lambda*lambda);

	return var;
}

double betapdf(double x, double alpha, double beta)
{
 	//Symbolic b = ( tgamma(alpha+beta)/(tgamma(alpha)*tgamma(beta)) ) * (x^(Symbolic(alpha-1))) * ((1-x)^(Symbolic(beta-1)));
	//return evalf(b,1,1);

	//boost::math::beta_distribution<> my_beta(alpha,beta);
	//double pdf_value = boost::math::pdf(my_beta,x);

	if (x < 0.0 || x > 1.0 || alpha <= 0.0 || beta <= 0.0) 
	{
		// Handle invalid input: PDF is 0 outside [0,1] or parameters are invalid
		return 0.0;	
	}
	else 
	{
		// Calculate log of the Beta function using log-Gamma functions
		double log_beta_function = std::lgamma(alpha) + std::lgamma(beta) - std::lgamma(alpha + beta);

		// Calculate log of the numerator
		double log_numerator = (alpha - 1.0) * std::log(x) + (beta - 1.0) * std::log(1.0 - x);

		// Calculate the log of the PDF
		double log_pdf = log_numerator - log_beta_function;

		// Return the exponentiated value
		return std::exp(log_pdf);
	}
	
}

double betacdf(double x, double alpha, double beta)
{	
	/*
	boost::math::beta_distribution<> my_beta(alpha,beta);
	double cdf_value = boost::math::cdf(my_beta,x);*/

	// Precision with boost::math::beta_distribution is insanely exact
	double cdf_value = incbeta(x,alpha,beta);
	return cdf_value;
}

Symbolic betamgf(double x, double alpha, double beta)
{
 	Symbolic t("t");
	Symbolic mgf = hypergeometric_1F1(alpha, alpha+beta,t,5);

	return mgf;
}

double betamean(double x, double alpha, double beta)
{
 	double mean = divisiond(alpha,alpha+beta);

	return mean;
}

double betavar(double x, double alpha, double beta)
{
 	double var = divisiond(alpha*beta,(alpha+beta+1)*(alpha+beta)*(alpha+beta));

	return var;
}

double cauchypdf(double x)
{
 	double pdf = divisiond(1,π)*divisiond(1,x*x+1);

	return pdf;
}


double chisquaredpdf(double x, double r)
{
	int a = 0.5*r;
	int c = 0.5*r-1;
 	double pdf = divisionint(1,tgamma(0.5*r)*pow(2,a)) * pow(x,c) * exp(-0.5*x);

	return pdf;
}


double chisquaredcdf(double x, double r)
{
	//double degrees_of_freedom = r; 
	//boost::math::chi_squared_distribution<> chi_squared(degrees_of_freedom);

	/*The cumulative distribution function (CDF) for Chi-Squared distribution
	It has at least 16 decimal digit precision with boost::math::chi_squared_distribution
	*/
	double cdf_value =(1/(tgamma(r/2.0)))*lowergamma(r/2.0,x/2.0);

	return cdf_value;
}

Symbolic chisquaredmgf(double x, double r)
{
 	Symbolic t("t");
	Symbolic mgf= (1-2*t)^(Symbolic(-0.5*r));

	return mgf;
}

double chisquaredmean(double x, double r)
{
 	double mean = r;

	return mean;
}

double chisquaredvar(double x, double r)
{
 	double var = 2*r;

	return var;
}

double lowergamma(double s, double x)
{
	int sgn = 1;
	double l;
	for (int i = 0; i < 100 ; i++)
	{
		l += sgn*pow(x,i+s)/(factoriald(i)*(i+s));
		sgn = -sgn;
	}
	return l;
}

#include <math.h>

#define STOP 1.0e-8
#define TINY 1.0e-30

double incbeta(double x, double a, double b) 
{
	if (x < 0.0 ) 
	{
		return 1.0/0.0;
	}	
	if (x > 1.0) 
	{
		x = 1/(1+(a*x)/(b));
		cout << x << endl;
	}
	/*The continued fraction converges nicely for x < (a+1)/(a+b+2)*/
	if (x > (a+1.0)/(a+b+2.0)) 
	{
		return (1.0-incbeta(1.0-x,b,a)); /*Use the fact that beta is symmetrical.*/
	}

	/*Find the first part before the continued fraction.*/
	const double lbeta_ab = lgamma(a)+lgamma(b)-lgamma(a+b);
	const double front = exp(log(x)*a+log(1.0-x)*b-lbeta_ab) / a;

	/*Use Lentz's algorithm to evaluate the continued fraction.*/
	double f = 1.0, c = 1.0, d = 0.0;

	int i, m;
	for (i = 0; i <= 200; ++i) 
	{
		m = i/2;

		double numerator;
		if (i == 0) 
		{
			numerator = 1.0; /*First numerator is 1.0.*/
		} 
		else if (i % 2 == 0) 
		{
			numerator = (m*(b-m)*x)/((a+2.0*m-1.0)*(a+2.0*m)); /*Even term.*/
		} 
		else 
		{
			numerator = -((a+m)*(a+b+m)*x)/((a+2.0*m)*(a+2.0*m+1)); /*Odd term.*/
		}

		/*Do an iteration of Lentz's algorithm.*/
		d = 1.0 + numerator * d;
		if (fabs(d) < TINY) 
		{
			d = TINY;
		}
		d = 1.0 / d;

		c = 1.0 + numerator / c;
		if (fabs(c) < TINY) 
		{
			c = TINY;
		}
		const double cd = c*d;
		f *= cd;

		/*Check for stop.*/
		if (fabs(1.0-cd) < STOP) 
		{
			return front * (f-1.0);
		}
	}

	return 1.0/0.0; /*Needed more loops, did not converge.*/
}

double Fpdf(double x, double r1, double r2)
{
	//boost::math::fisher_f_distribution<> fisher_f(r1,r2);
	//double pdf_value = boost::math::pdf(fisher_f,x);

	// this function returns the same function as boost::math::fisher_f_distribution<> fisher_f(r1,r2); for F pdf
	double pdf_value = (tgamma(0.5*(r1+r2)))*(pow(r1/r2,0.5*r1))/(tgamma(0.5*r1)*tgamma(0.5*r2)) * (pow(x,(0.5*r1)-1)/(pow(1+(r1*x/r2),0.5*(r1+r2)))) ;

	return pdf_value;
}

double Fcdf(double x, double r1, double r2)
{
	//boost::math::fisher_f_distribution<> fisher_f(r1,r2);
	//double cdf_value = boost::math::cdf(fisher_f,x);

	/*The cumulative distribution function (CDF) for Fisher's F distribution
	It has at least 11 decimal digit precision with boost::math::fisher_f_distribution
	*/

	double z ;
	double cdf_value;
	
	z = (x*r1)/(r2 + r1*x);
	cdf_value = incbeta(z,r1/2.0, r2/2.0);


	return cdf_value;
}

double Fmean(double x, double r1, double r2)
{
	double mean = divisions(r2,r2-2) ;

	return mean;
}

double Fvar(double x, double r1, double r2)
{
	double var = 2*divisions(r2,r2-2)*divisions(r2,r2-2) * divisions(r1+r2-2,r1*(r2-4));

	return var;
}

double tpdf(double x, double r)
{
	//boost::math::students_t_distribution<> students_t(r);
	//double pdf_value = boost::math::pdf(students_t,x);

	/*The PDF for Student's t distribution
	It has 7 decimal digit precision with boost::math::students_t_distribution
	*/
	double z = (tgamma((r+1)/2.0))/(tgamma(r/2.0) * sqrt(π*r)) ;
	double pdf_value = z * pow((1 + (x*x)/(r)),-((r+1)/2.0));
	return pdf_value;
}

double tcdf(double x, double r)
{
	//boost::math::students_t_distribution<> students_t(r);
	//double cdf_value = boost::math::cdf(students_t,x);

	/*The cumulative distribution function (CDF) for Student's t distribution
	It has 15 decimal digit precision with boost::math::students_t_distribution
	*/
	double t = (x + sqrt(x * x + r)) / (2.0 * sqrt(x * x + r));
	double cdf_value = incbeta(t, r/2.0, r/2.0);

	return cdf_value;
}

double tmean(double x, double r)
{
	double mean = 0;

	return mean;
}

double tvar(double x, double r)
{
	double var = divisiond(r,r-2);

	return var;
}

Symbolic laplacepdf(double x, double θ)
{
	boost::math::laplace_distribution<> laplace(θ);
	double pdf_value = boost::math::pdf(laplace,x);

	return pdf_value;
}

Symbolic laplacecdf(double x, double θ)
{
	boost::math::laplace_distribution<> laplace(θ);
	double cdf_value = boost::math::cdf(laplace,x);

	return cdf_value;
}

Symbolic laplacemgf(double x, double θ)
{
 	Symbolic t("t");
	Symbolic mgf= exp(t*θ) * ( (1-t*t)^(Symbolic(-1)) );

	return mgf;
}

Symbolic laplacemean(double x, double θ)
{
	Symbolic mean = θ;

	return mean;
}

Symbolic laplacevar(double x, double θ)
{
	Symbolic var = 2;

	return var;
}

Symbolic logisticpdf(double x, double θ)
{
	boost::math::logistic_distribution<> logistic(θ);
	double pdf_value = boost::math::pdf(logistic,x);

	return pdf_value;
}

Symbolic logisticcdf(double x, double θ)
{
	boost::math::logistic_distribution<> logistic(θ);
	double cdf_value = boost::math::cdf(logistic,x);

	return cdf_value;
}

Symbolic logisticmean(double x, double θ)
{
	Symbolic mean = θ;

	return mean;
}

Symbolic logisticvar(double x, double θ)
{
	Symbolic var = divisions(π*π,3);

	return var;
}



// Function to find the mode of a vector
double findMode(const std::vector<double>& data) 
{
	if (data.empty()) 
	{
        // Handle empty input case, e.g., throw an exception or return a special value
	return 0; // Or some indicator of no mode
	}

	std::map<double, double> frequencyMap;
	for (double value : data) 
	{
		frequencyMap[value]++;
	}

	double mode = data[0]; // Initialize with the first element
	double maxFrequency = 0;

	for (const auto& pair : frequencyMap) 
	{
		if (pair.second > maxFrequency) 
		{
			maxFrequency = pair.second;
			mode = pair.first;
		}
	}
	return mode;
}

double descriptivestatistics(vector<double> vector_x)
{
	int n = vector_x.size();
	double mean, variance, stdev;
	cout << "\nSample size : " << n << endl; 
	int modeValue = findMode(vector_x);
 	cout << "\nMode : " << modeValue << endl; 

	std::sort(vector_x.begin(), vector_x.end());
	// Median: The middle value in a sorted dataset. If the dataset has an even number of elements, it's the average of the two middle values.
	if (n % 2 == 1) 
	{
		cout << "Median : " << vector_x[n / 2] << endl;
        } 
	else 
	{
		cout << "Median : " << (vector_x[n / 2 - 1] + vector_x[n / 2]) / 2.0 << endl;
	}	

	double sum = 0.0;
	for (double val : vector_x) 
	{
		sum += val;
	}
        mean = sum / n;
	
	cout << "Mean : " << mean << endl;

	// Assuming the vector containing the 'data' is already sorted
	double p = 0.25; // For the 1/4
	double index_double = p * (n - 1);
	size_t lower_index = static_cast<size_t>(std::floor(index_double));
	size_t upper_index = static_cast<size_t>(std::ceil(index_double));
	
	double quantile_value, quantile_value2;
	if (lower_index == upper_index) 
	{
		quantile_value = vector_x[lower_index];
	} 
	else 
	{
		double weight = index_double - lower_index;
		quantile_value = vector_x[lower_index] * (1.0 - weight) + vector_x[upper_index] * weight;
	}

	double p2 = 0.75; // For the 3/4
	double index_double2 = p2 * (n - 1);
	size_t lower_index2 = static_cast<size_t>(std::floor(index_double2));
	size_t upper_index2 = static_cast<size_t>(std::ceil(index_double2));
	if (lower_index2 == upper_index2) 
	{
		quantile_value2 = vector_x[lower_index2];
	} 
	else 
	{
		double weight2 = index_double2 - lower_index2;
		quantile_value2 = vector_x[lower_index2] * (1.0 - weight2) + vector_x[upper_index2] * weight2;
	}
	cout << "Quantile 1/4: " << quantile_value << endl;

	cout << "Quantile 3/4: " << quantile_value2 << endl;

	double sumSquaredDiff = 0.0;
	for (double val : vector_x) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	variance = sumSquaredDiff / (n - 1); // Sample variance

	cout << "Variance : " << variance << endl;
	stdev = sqrt(variance);
	cout << "Standard deviation : " << stdev;
	return stdev;
}
double rpearson(const SymbolicMatrix &A, int N)
{
	double sum_xy, sum_x, sum_y, x_bar, y_bar, sum_xsquared, sum_ysquared, Sxy, Sx, Sy, r_pearson;
	sum_xy = 0;
	sum_x = 0;
	sum_y = 0;
	sum_xsquared = 0;
	sum_ysquared = 0;
	for(int i=0; i < N; ++i)
	{
		sum_x = sum_x + A[i][0];
		sum_y = sum_y + A[i][1];
		sum_xy = sum_xy +(A[i][0]*A[i][1]);
		sum_xsquared = sum_xsquared + (A[i][0]*A[i][0]);
		sum_ysquared = sum_ysquared + (A[i][1]*A[i][1]);
		
	}	

	x_bar = sum_x/N;
	y_bar = sum_y/N;

	Sxy = (sum_xy/N) - (x_bar*y_bar);
	Sx = sqrt(sum_xsquared/N - (x_bar*x_bar));
	Sy = sqrt(sum_ysquared/N - (y_bar*y_bar));
	r_pearson = Sxy/(Sx*Sy);

	cout << "x bar: " << x_bar <<endl;
	cout << "y bar: " << y_bar <<endl;
	cout << "sum x: " << sum_x <<endl;
	cout << "sum y: " << sum_y <<endl;
	cout << "sum xy: " << sum_xy <<endl;
	cout << "sum x^2: " << sum_xsquared <<endl;
	cout << "sum y^2: " << sum_ysquared <<endl;

	cout << "\nPearson's correlation coefficient, r : " << endl;

	return r_pearson;
}

Symbolic regressionline(const SymbolicMatrix &A, int N)
{
	Symbolic x("x"), y("y");
	double sum_xy, sum_x, sum_y, x_bar, y_bar, sum_xsquared, Sxy, Sx;
	sum_xy = 0;
	sum_x = 0;
	sum_y = 0;
	sum_xsquared = 0;
	for(int i=0; i < N; ++i)
	{
		sum_x = sum_x + A[i][0];
		sum_y = sum_y + A[i][1];
		sum_xy = sum_xy +(A[i][0]*A[i][1]);
		sum_xsquared = sum_xsquared + (A[i][0]*A[i][0]);	
	}	

	x_bar = sum_x/N;
	y_bar = sum_y/N;

	Sxy = (sum_xy/N) - (x_bar*y_bar);
	Sx = sqrt(sum_xsquared/N - (x_bar*x_bar));
	
	Symbolic regression_line = (Sxy/(Sx*Sx))*(x-x_bar) + y_bar;

	return regression_line;
}

vector<vector<double>> multipleregression(vector<vector<double>> &X, vector<vector<double>> &y)
{
	dmat X_t = transpose(X);
	dmat XtX= multiply(X_t,X);
	dmat inv = inverse(XtX);
	dmat A = multiply(inv,X_t);
	dmat b = multiply(A,y);

	cout << endl;

	int n = b.size();
	for(int i = 0; i < n ; ++i)
	{
		cout <<"b["<< i <<"] = " << b[i][0] << endl;
	}
	return b;
}
// Function to calculate the mean of a vector
double calculateMean(vector<double> data) 
{
	return std::accumulate(data.begin(), data.end(), 0.0) / data.size();
}

// Function to calculate the covariance between two vectors
double calculateCovariance(vector<double> data1, vector<double> data2) 
{
	if (data1.size() != data2.size() || data1.empty()) 
	{
		return 0.0; // Handle error or return appropriate value
	}

	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sum_of_products = 0.0;

	for (size_t i = 0; i < data1.size(); ++i) 
	{
		sum_of_products += (data1[i] - mean1) * (data2[i] - mean2);
	}

	return sum_of_products / (data1.size() - 1); // Sample covariance
}

SymbolicMatrix covariancematrix(vector<vector<double>> matrix)
{
	vector<vector<double>> cov_matrix;
	
	//int R = matrix.size();
	int C = matrix[0].size();

	Matrix<Symbolic> B_mat(C,C);

	for(int i = 0; i < C; i++)
	{
		for(int j=0; j < C; j++)
		{
			B_mat[i][j] = calculateCovariance(getColumn(matrix, i),getColumn(matrix, j));	
		}
	}

	cout << "\nThe covariance matrix:" <<endl;
	
	return B_mat;
}

void hypothesistest(vector<double> data, double mu, double sigma, double alpha, double effect_size)
{
	double critical_value, critical_value2, z_quantile, z_value, z_beta, meanbeta, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n = data.size();
	double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	double sigma_sample = sqrt(variance);

	std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nHypothesis testing two-tailed" << endl;
	cout << "\nH0: μ = " << mu << endl;
	cout << "H1: μ != " << mu << endl;
	
	if (n < 30) // For small sample size we use t-statistic
	{
		t_value = (mean-mu)/(sigma_sample/sqrt(n));
		p_value = tcdf(t_value,n);

		critical_value = mu + t_quantile*(sigma_sample)/(sqrt(n)) ;
		meanbeta = mu + effect_size; // for right-tailed

		t_beta = (critical_value-meanbeta)/(sigma_sample/sqrt(n));
		beta = tcdf(t_beta,n-1);

		p_value = 1 - tcdf(t_value,n-1);
	}
	else if (n >= 30) // For large sample size we use Z-statistic
	{
		z_value = ((mu-effect_size)-mu)/(sigma/sqrt(n));
		
		z_quantile = zquantile(1-(0.5*alpha));

		critical_value = mu + z_quantile*(sigma)/(sqrt(n)) ; // upper bound
		critical_value2 = mu - z_quantile*(sigma)/(sqrt(n)) ; // lower bound

		meanbeta = mu + effect_size; // for two-tailed
		z_beta = (critical_value-meanbeta)/(sigma/sqrt(n));
		
		beta = (1 - 2*normalcdf(z_beta,0,1) );	
		power = 1-beta;
		p_value = 2*normalcdf(z_value,0,1);	
		
		cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
		cout << "The probability of a type II error (beta): " << beta << endl;
		cout << "Power: " << power << endl;
		cout << "\nCritical value (upper bound): " << critical_value << endl;
		cout << "Critical value (lower bound): " << critical_value2 << endl;
		cout << "Critical region : z < - " << z_quantile << " and z > " << z_quantile << endl;
		cout << "Computed z : " << z_value << endl;
		cout << "\n*Reject the null hypothesis when the sample average is greater than " << critical_value << " or less than " << critical_value2 << endl;
		cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << " or smaller than -" << z_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
	
	}
}

void hypothesistest_righttailed(vector<double> data, double mu, double sigma, double alpha, double effect_size)
{
	double critical_value, z_quantile, z_value, z_beta, meanbeta, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n = data.size();
	double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	double sigma_sample = sqrt(variance);

	std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nHypothesis testing right-tailed" << endl;
	cout << "\nH0: μ = " << mu << endl;
	cout << "H1: μ > " << mu << endl;
	
	if (n < 30) // For small sample size we use t-statistic
	{
		t_value = (mean-mu)/(sigma_sample/sqrt(n));
		p_value = tcdf(t_value,n);

		critical_value = mu + t_quantile*(sigma_sample)/(sqrt(n)) ;
		meanbeta = mu + effect_size; // for right-tailed

		t_beta = (critical_value-meanbeta)/(sigma_sample/sqrt(n));
		beta = tcdf(t_beta,n-1);

		p_value = 1 - tcdf(t_value,n-1);
	}
	else if (n >= 30) // For large sample size we use Z-statistic
	{
		z_value = ((mu+effect_size)-mu)/(sigma/sqrt(n));
		
		z_quantile = zquantile(1-alpha);

		critical_value = mu + z_quantile*(sigma)/(sqrt(n)) ;

		meanbeta = mu + effect_size; // for right-tailed
		z_beta = (critical_value-meanbeta)/(sigma/sqrt(n));
		beta = normalcdf(z_beta,0,1);	
		power = 1-beta;
		p_value = 1 - normalcdf(z_value,0,1);	
		
		cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
		cout << "The probability of a type II error (beta): " << beta << endl;
		cout << "Power: " << power << endl;
		cout << "\nCritical value: " << critical_value << endl;
		cout << "Critical region : z > " << z_quantile << endl;
		cout << "Computed z : " << z_value << endl;
		cout << "\n*Reject the null hypothesis when the sample average is greater than " << critical_value << endl;
		cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
	
	}
}


void ANOVA(vector<vector<double>> matrix)
{
	vector<vector<double>> anovamatrix;
	vector<double> total_column;
	vector<double> mean_column;
	vector<double> sst_column;
	vector<double> ssa_column;
	int C = matrix[0].size();
	int N = 0; // computing total data

	for(int i = 0 ; i < C ; ++i)
	{
		double sum = 0;
		int R = matrix.size();		
		for (int j = 0 ; j < R ; ++j)
		{
			sum += matrix[j][i]; // sum of the column
		
			if(matrix[j][i] !=0 )		
			{
				N = N+1;
			}
			else if(matrix[j][i] == 0 )		
			{
				N = N;
			}

		}
			total_column.push_back(sum);
			mean_column.push_back(sum/(R));
			//cout << "sum of column " << i << " = " << total_column[i] << endl;
			//cout << "mean sum of column " << i << " = " << mean_column[i] << endl;
			
	}
	
	double Y = 0.0;
	double y_bar = 0.0;
	for (double val : total_column) 
	{
		Y += val;
	}
        for (double val : mean_column) 
	{
		y_bar += val;
	}
	y_bar = y_bar/C;
	
	for(int i = 0 ; i < C ; ++i)
	{
		double sst_temp = 0.0;
		vector<double> v_col = getColumn(matrix,i);
		for (double val : v_col)
		{
			sst_temp += (val-y_bar)*(val-y_bar);
		}
		sst_column.push_back(sst_temp);
	}
	double SST = 0.0; // Total sum of squares
	for (double val : sst_column)
	{
		SST += val;
	}
	
	for(int i = 0 ; i < C ; ++i)
	{
		double ssa_temp = 0.0;
		vector<double> v_col2 = getColumn(matrix,i);
		int n = v_col2.size();

		ssa_temp = n*(mean_column[i] - y_bar)*(mean_column[i] - y_bar);
		
		ssa_column.push_back(ssa_temp);
	}
	double SSA = 0.0; // Treatment sum of squares
	for (double val : ssa_column)
	{
		SSA += val;
	}	
	double SSE = SST-SSA;
	
	int df_model = C-1;
	int df_total = N-1;
	int df_error = df_total - df_model;

	double s1_square = SSA/(df_model);
	double s1 = SSE/(df_error);
	double computed_f = s1_square/s1; 

	int r1 = df_model;
	int r2 = df_error;
	double p_value = 1 - Fcdf(computed_f,r1,r2);
	double r_square = 1 - SSE/SST;

	cout << "\nSource" << setw(10) << "DF" << setw(23) << "SS" << setw(23) << "MS" << setw(23) << "F" << setw(23) << "P" << endl;
	cout << "\nModel" << setw(10) <<  df_model << setw(23) << SSA << setw(23) << s1_square << setw(23) << computed_f << setw(23) << p_value << endl;
	cout << "\nError" << setw(10) << df_error << setw(23) << SSE << setw(23) << s1 << setw(23) << endl;
	cout << "\nTotal" << setw(10) << df_total << setw(23) << SST << endl;
	
	cout << "\nR-Square"<< setw(23) << "Grand Mean" << setw(23) << "Root MSE" << endl;
	cout << r_square << setw(23) << y_bar << setw(23) << sqrt(s1) << endl;
	
}
#endif
#endif