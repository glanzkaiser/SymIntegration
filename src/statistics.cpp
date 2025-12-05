/*
Thanks Freya, Berlin and Sentinel!
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
#include <iostream>
#include <string>
#include <fstream> // For file operations

#include <random> // For random number generation
#include <chrono>

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_DEFINE
#define π 3.141592653589793238462643383279502884f

#define EPSILON 0.00001

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
Symbolic factorialsym(Symbolic n) { // Dummy factorial symbolic for symbolic computational purpose like MLE
	
	Symbolic result = n*(n-1)*(n-2)*(n-3)*1;
	
	return result;
}

// Factorial and combinations
Symbolic factorials(int n) {
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

Symbolic combinationss(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factorials(n) / (factorials(r) * factorials(n - r));
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

vector<string> loadStringVector(const string& filename) 
{
	vector<string> vecx;
	ifstream inputFile(filename);

	if (!inputFile.is_open()) 
	{
		cerr << "Error: Could not open file " << filename << endl;
		return vecx; // Return empty matrix on error
	}

	if (inputFile.is_open()) 
	{
		string line;
		while (getline(inputFile,line)) 
		{ // Reads numbers separated by whitespace
			vecx.push_back(line);
		}
		inputFile.close();
		
	} 
	inputFile.close();
	return vecx;
}

void printStringVector(vector<string> vectorx)
{
	//int n = vectorx.size();
	for (const std::string& s : vectorx) 
	{
		cout << s << endl;
	}
}

void saveVectordouble(vector<double> vecx, const string& filename) 
{
	// Create an ofstream object and open the file for writing
	ofstream outputFile(filename);
	
	if (!outputFile.is_open()) 
	{
	cerr << "Error: Could not open the file for writing." << endl;
	}

	//  Iterate through the vector and write each element to the file
	for (double num : vecx) 	
	{
		outputFile << num << endl; // Write the number followed by a newline
	}

	// Close the file to release resources
	outputFile.close();
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
 	double pmf = exp(-k)*pow(k,x)/(factoriald(x));
	
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
	double z_quantile = r8_normal_01_cdf_inverse(cdf); // for double-precision floating point calculations

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
	double a = 0.5*r;
	double c = 0.5*r-1;
 	double pdf = divisiond(1,tgamma(0.5*r)*pow(2,a)) * pow(x,c) * exp(-0.5*x);

	return pdf;
}


double chisquaredcdf(double x, double r)
{
	//double degrees_of_freedom = r; 
	//boost::math::chi_squared_distribution<> chi_squared(degrees_of_freedom);

	/*The cumulative distribution function (CDF) for Chi-Squared distribution
	It has at least 16 decimal digit precision with boost::math::chi_squared_distribution
	*/
	
	// double cdf_value =(divisiond(1,(tgamma(0.5*r))))*lowergamma(0.5*r,0.5*x); 

	// We approximate the cdf using the Simpsons rule
	double b = x;
	double a = 0;
	int n = 108;
        double h = (b - a) / n;
        double sum = chisquaredpdf(a,r) + chisquaredpdf(b,r); // Endpoints
        for (int i = 1; i < n; i += 2) 
	{
            sum += 4 * chisquaredpdf(a + i * h,r); // Odd-indexed points
        }
        for (int i = 2; i < n; i += 2) 
	{
            sum += 2 * chisquaredpdf(a + i * h,r); // Even-indexed points
        }
        double cdf_value =  divisiond(sum*h,3.0) ;

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

double chisquaredquantile(double r, double cdf) 
{
	
	// Newton-Raphson
	double x0 = r; // nu is the best estimator for initial approximation
	double root;
	double h =( chisquaredcdf(x0,r) - cdf) / chisquaredpdf(x0,r);
	while (abs(h) >= EPSILON)
	{
		h = ( chisquaredcdf(x0,r) - cdf) / chisquaredpdf(x0,r);
		
		// x(i+1) = x(i) - f(x) / f'(x)  
		x0 = x0 - h;
	}
	root = x0;
	
	return root;

	/*
	// with Bisection method, slow but converges.
	double a = 0.00001;
	double b;
	if(r < 10 && r>0)
	{
		b = 10;
	}
	if(r < 10 && r>5 && cdf >0.5)
	{
		b = 21;
	}
	if(r < 30 && r > 20)
	{
		b = 100;
	}
	if(r >= 30)
	{
		b = 1000;
	}
	int N = 100;
	double c;
	for (int i = 0; i < N; ++i) 
	{
		c = (a + b) / 2;
		if (abs(divisiond(b - a, 2)) < EPSILON) 
		{
			return c; // Root found within tolerance
		}

		if ((chisquaredcdf(c,r) - cdf) * (chisquaredcdf(a,r) - cdf) < 0) 
		{
			b = c;
		} 
		else 
		{
			a = c;
		}
	}
	cerr << "Warning: Maximum iterations reached. Approximated root: " << c << endl;
	return c; // Return the best approximation found
	*/
}


double lowergamma(double s, double x)
{
	//int sgn = 1;
	double l;
	for (int i = 0; i < 101 ; i++)
	{
		l += divisiond(pow(-1,i)*pow(x,i+s),(factoriald(i)*(i+s)));
		//sgn = -sgn;
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

double Fquantile(double r1, double r2, double cdf)
{
	double x0 = 1;
	double root;
	double h =( Fcdf(x0,r1,r2) - cdf) / Fpdf(x0,r1,r2);
	while (abs(h) >= EPSILON)
	{
		h = ( Fcdf(x0,r1,r2) - cdf) / Fpdf(x0,r1,r2);
		
		// x(i+1) = x(i) - f(x) / f'(x)  
		x0 = x0 - h;
	}
	root = x0;
	return root;
}

double tpdf(double x, double r)
{
	//boost::math::students_t_distribution<> students_t(r);
	//double pdf_value = boost::math::pdf(students_t,x);

	/*The PDF for Student's t distribution
	It has 7 decimal digit precision with boost::math::students_t_distribution
	*/

	/*
	double pi = std::acos(-1.0); // Calculate pi
	double coefficient = std::tgamma(0.5 * (nu + 1.0)) /  (std::sqrt(nu * pi) * std::tgamma(0.5 * nu));
	double term = std::pow(1.0 + (x * x) / nu, -0.5 * (nu + 1.0));
	return coefficient * term;
	*/
	if (r <= 0) 
	{
		// Handle invalid degrees of freedom
		return 0.0; 
	}
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

double tquantile(double x0, double r, double cdf) // with Newton-Raphson method
{
	double root;
	double h =( tcdf(x0,r) - cdf) / tpdf(x0,r);
	while (abs(h) >= EPSILON)
	{
		h = ( tcdf(x0,r) - cdf) / tpdf(x0,r);
		
		// x(i+1) = x(i) - f(x) / f'(x)  
		x0 = x0 - h;
	}
	root = x0;
	return root;
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

vector<double> normalquantile(vector<double> vector_x, double scale)
{
	vector<double> nquantile;
	sort(vector_x.begin(), vector_x.end());
	double mean = calculateMean(vector_x);
	int n = vector_x.size();	
		
	for(int i = 0 ; i < n ; ++i)
	{
		nquantile.push_back(scale*(vector_x[i] - mean));
	} 
	return nquantile;
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

void confidenceinterval_onesampletwosides(vector<double> data, double sigma, double alpha)
{
	double mean, z_quantile, lower_bound, upper_bound;
	int n = data.size();

	mean = calculateMean(data);

	z_quantile = zquantile(1-(0.5*alpha));
	lower_bound = mean - z_quantile*(divisiond(sigma,sqrt(n)));
	upper_bound = mean + z_quantile*(divisiond(sigma,sqrt(n)));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval one sample two sides, σ known" << endl;
	cout << "Mean data = "<< mean << ", σ = " << sigma << ", n = " << n << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < μ1 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_onesampleonesided(vector<double> data, double sigma, double alpha)
{
	double mean, z_quantile, lower_bound, upper_bound;
	int n = data.size();

	mean = calculateMean(data);

	z_quantile = zquantile(1-alpha);
	lower_bound = mean - z_quantile*(divisiond(sigma,sqrt(n)));
	upper_bound = mean + z_quantile*(divisiond(sigma,sqrt(n)));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval one sample one-sided confidence bounds, σ known" << endl;
	cout << "Mean data = "<< mean << ", σ = " << sigma << ", n = " << n << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << "Upper one-sided bound: "<< upper_bound << endl;	
	cout << "Lower one-sided bound: "<< lower_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_onesampletwosides(vector<double> data, double alpha)
{
	double mean, sigma, t_quantile, lower_bound, upper_bound;
	int n = data.size();

	mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	sigma = sqrt(variance);

	t_quantile = tquantile(1,n-1,1-(0.5*alpha));
	lower_bound = mean - t_quantile*(divisiond(sigma,sqrt(n)));
	upper_bound = mean + t_quantile*(divisiond(sigma,sqrt(n)));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval one sample two sides, σ unknown" << endl;
	cout << "Mean data = "<< mean << ", σ = " << sigma << ", n = " << n << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < μ1 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_onesampleonesided(vector<double> data, double alpha)
{
	double mean, sigma, t_quantile, lower_bound, upper_bound;
	int n = data.size();

	mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	sigma = sqrt(variance);

	t_quantile = tquantile(1,n-1,1-(alpha));
	lower_bound = mean - t_quantile*(divisiond(sigma,sqrt(n)));
	upper_bound = mean + t_quantile*(divisiond(sigma,sqrt(n)));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval one sample one-sided confidence bounds, σ unknown" << endl;
	cout << "Mean data = "<< mean << ", σ = " << sigma << ", n = " << n << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << "Upper one-sided bound: "<< upper_bound << endl;	
	cout << "Lower one-sided bound: "<< lower_bound << endl;		
	cout << "\n********************************************************" << endl;
}

void predictioninterval(vector<double> data, double sigma, double alpha)
{
	double mean, z_quantile, lower_bound, upper_bound;
	int n = data.size();

	mean = calculateMean(data);

	z_quantile = zquantile(1-(0.5*alpha));
	lower_bound = mean - z_quantile*sigma*(sqrt(1+divisiond(1,n)));
	upper_bound = mean + z_quantile*sigma*(sqrt(1+divisiond(1,n)));

	cout << "\n********************************************************" << endl;
	cout << "\nPrediction interval for one sample data, σ known" << endl;
	cout << "Mean data = "<< mean << ", σ = " << sigma << ", n = " << n << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% prediction interval is :\n" << endl;	
	cout << lower_bound << " < X0 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void predictioninterval(vector<double> data, double alpha)
{
	double mean, sigma, t_quantile, lower_bound, upper_bound;
	int n = data.size();

	mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	sigma = sqrt(variance);

	t_quantile = tquantile(1,n-1,1-(0.5*alpha));
	lower_bound = mean - t_quantile*sigma*(sqrt(1+divisiond(1,n)));
	upper_bound = mean + t_quantile*sigma*(sqrt(1+divisiond(1,n)));

	cout << "\n********************************************************" << endl;
	cout << "\nPrediction interval for one sample data, σ unknown" << endl;
	cout << "Mean data = "<< mean << ", s = " << sigma << ", n = " << n << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% prediction interval is :\n" << endl;	
	cout << lower_bound << " < X0 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}


void confidenceinterval_knownsigma(vector<double> data1, vector<double> data2, double sigma1, double sigma2, double alpha)
{
	double mean1, mean2, z_quantile, lower_bound, upper_bound;
	int n1 = data1.size();
	int n2 = data2.size();

	mean1 = calculateMean(data1);
	/*double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}*/
	mean2 = calculateMean(data2);
	/*double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance1 = sumSquaredDiff1 / (n - 1); // Sample variance
	double sigma_sample1 = sqrt(variance1);
	double variance2 = sumSquaredDiff2 / (n - 1); // Sample variance
	double sigma_sample2 = sqrt(variance2);
	*/

	z_quantile = zquantile(1-(0.5*alpha));
	lower_bound = abs(mean2 - mean1) - z_quantile*(sqrt(pow(sigma1,2)/n1 + pow(sigma2,2)/n2));
	upper_bound = abs(mean2 - mean1) + z_quantile*(sqrt(pow(sigma1,2)/n1 + pow(sigma2,2)/n2));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval with population variances known" << endl;
	cout << "Mean data 1 = "<< mean1 << ", σ1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 <<  ", σ2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < μ1 - μ2 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_sameunknownsigma(vector<double> data1, vector<double> data2, double alpha)
{
	double mean1, mean2, Sp, t_quantile, lower_bound, upper_bound, sigma1, sigma2;
	int n1 = data1.size();
	int n2 = data2.size();

	mean1 = calculateMean(data1);
	double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	double variance1 = sumSquaredDiff1 / (n1 - 1); // Sample variance
	double s1 = sqrt(variance1);
	
	mean2 = calculateMean(data2);
	double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance2 = sumSquaredDiff2 / (n2 - 1); // Sample variance
	double s2 = sqrt(variance2);
	sigma1 = s1;
	sigma2 = s2;

	Sp = ((n1-1)*pow(sigma1,2) + (n2-1)*pow(sigma2,2))/(n1+n2-2);
	Sp = sqrt(Sp);
	t_quantile = tquantile(1,n1+n2-2,1-(0.5*alpha));
	lower_bound = abs(mean2 - mean1) - t_quantile*Sp*(sqrt(divisiond(1,n1) + divisiond(1,n2)));
	upper_bound = abs(mean2 - mean1) + t_quantile*Sp*(sqrt(divisiond(1,n1) + divisiond(1,n2)));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval with same population variances but both unknown" << endl;
	cout << "Mean data 1 = "<< mean1 << ", s1 = " << s1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << s2 << ", n2 = " << n2 << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < μ1 - μ2 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_unequalunknownsigma(vector<double> data1, vector<double> data2, double alpha)
{
	double mean1, mean2, t_quantile, lower_bound, upper_bound, sigma1, sigma2;
	int nu;
	int n1 = data1.size();
	int n2 = data2.size();

	mean1 = calculateMean(data1);
	double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	double variance1 = sumSquaredDiff1 / (n1 - 1); // Sample variance
	double s1 = sqrt(variance1);
	
	mean2 = calculateMean(data2);
	double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance2 = sumSquaredDiff2 / (n2 - 1); // Sample variance
	double s2 = sqrt(variance2);
	sigma1 = s1;
	sigma2 = s2;

	nu = pow((divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)),2)/(divisiond(pow(divisiond(sigma1*sigma1,n1),2),n1-1) + divisiond(pow(divisiond(sigma2*sigma2,n2),2),n2-1));
	t_quantile = tquantile(1,nu,1-(0.5*alpha));
	lower_bound = abs(mean2 - mean1) - t_quantile*(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));
	upper_bound = abs(mean2 - mean1) + t_quantile*(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval with unequal and unknown population variances" << endl;
	cout << "Mean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < μ1 - μ2 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_paired(vector<double> data1, vector<double> data2, double alpha)
{
	vector<double> d;
	double meand, t_quantile, lower_bound, upper_bound, sigmad;
	int n1 = data1.size();
	int n2 = data2.size();

	if (n1 != n2)
	{
		cerr << "Error: The samples size are different." << endl;
	}
	else if(n1==n2)
	{
		for (int i = 0; i < n1; ++i) 
		{
			d.push_back(data1[i]-data2[i]);
		}
		
		meand = calculateMean(d);
		double sumSquaredDiff = 0.0;
		for (double val : d) 
		{
			sumSquaredDiff += std::pow(val - meand, 2);
		}
		double variance = sumSquaredDiff / (n1 - 1); // Sample variance
		sigmad = sqrt(variance);

		t_quantile = tquantile(1,n1-1,1-(0.5*alpha));
		lower_bound = meand - t_quantile*(divisiond(sigmad,sqrt(n1)));
		upper_bound = meand + t_quantile*(divisiond(sigmad,sqrt(n1)));

		cout << "\n********************************************************" << endl;
		cout << "\nConfidence interval for paired observations" << endl;
		cout << "Mean d = "<< meand << ", sd = " << sigmad << ", n = " << n1 << endl;	
		cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
		cout << lower_bound << " < μD < " << upper_bound << endl;	
		cout << "\n********************************************************" << endl;
	}
}

void confidenceinterval_proportion(vector<string> data, const string &inputString, double alpha)
{
	int x_count = std::count(data.begin(), data.end(), inputString);
	
	double z_quantile, lower_bound, upper_bound, p, q;
	int n = data.size();

	p = divisiond(x_count,n);
	q = 1-p;

	z_quantile = zquantile(1-(0.5*alpha));
	// method 2
	double point_estimate = (p+divisiond(z_quantile*z_quantile,2*n))/(1+divisiond(z_quantile*z_quantile,n));
	double margin_of_error =  divisiond(z_quantile,1+divisiond(z_quantile*z_quantile,n))*(sqrt(divisiond(p*q,n) + divisiond(z_quantile*z_quantile,4*n*n)));
	 
	lower_bound = point_estimate - margin_of_error ;
	upper_bound = point_estimate + margin_of_error ;

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval for proportion" << endl;
	cout << "Data size =  "<< n << endl;	
	cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
	cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < p < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_differencebetweentwoproportions(vector<string> data1, vector<string> data2, const string &inputString, double alpha)
{
	int x_count = std::count(data1.begin(), data1.end(), inputString);
	int y_count = std::count(data2.begin(), data2.end(), inputString);
	
	double z_quantile, lower_bound, upper_bound, p1, q1, p2, q2;
	int n1 = data1.size();
	int n2 = data2.size();

	p1 = divisiond(x_count,n1);
	q1 = 1-p1;
	p2 = divisiond(y_count,n2);
	q2 = 1-p2;

	z_quantile = zquantile(1-(0.5*alpha));
	// method 2
	double point_estimate = (p1 - p2);
	double margin_of_error =  z_quantile*sqrt(divisiond(p1*q1,n1) +divisiond(p2*q2,n2));
	 
	lower_bound = point_estimate - margin_of_error ;
	upper_bound = point_estimate + margin_of_error ;

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval for difference between two proportions" << endl;
	cout << "Data 1 size =  "<< n1 << ", Data 2 size =  "<< n2 << endl;	
	cout << "Amount of "<< inputString << " in the data 1 = " << x_count << endl;		
	cout << "Amount of "<< inputString << " in the data 2 = " << y_count << endl;		
	cout << "Proportion for "<< inputString << ", p1 = " << p1 << ", q1 = " << q1 << endl;	
	cout << "Proportion for "<< inputString << ", p2 = " << p2 << ", q2 = " << q2 << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < p1 - p2 < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_variance(vector<double> data, double alpha)
{
	double mean, chi_quantile1, chi_quantile2, lower_bound, upper_bound;
	int n = data.size();
	mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance = sumSquaredDiff / (n - 1); // Sample variance

	double r = n-1;
	chi_quantile1 = chisquaredquantile(r,1-(0.5*alpha));
	chi_quantile2 = chisquaredquantile(r,0.5*alpha);
	cout << chi_quantile1 << endl;
	cout << chi_quantile2 << endl;
	// method 
	lower_bound = divisiond(r*variance,chi_quantile1);
	upper_bound = divisiond(r*variance,chi_quantile2);

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval for variance" << endl;
	cout << "Data size =  "<< n << endl;	
	cout << "Mean sample =  "<< mean << ", variance sample = " << variance << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " < σ^{2} < " << upper_bound << endl;	
	cout << "\n********************************************************" << endl;
}

void confidenceinterval_ratiotwovariances(vector<double> data1, vector<double> data2, double alpha)
{
	double mean1, mean2, F_quantile1, F_quantile2, lower_bound, upper_bound;
	int n1 = data1.size();
	int n2 = data2.size();

	mean1 = calculateMean(data1);
	double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	double variance1 = sumSquaredDiff1 / (n1 - 1); // Sample variance
	double s1 = sqrt(variance1);
	
	mean2 = calculateMean(data2);
	double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance2 = sumSquaredDiff2 / (n2 - 1); // Sample variance
	double s2 = sqrt(variance2);
	
	F_quantile1 = Fquantile(n1-1,n2-1,1-(0.5*alpha));
	F_quantile2 = Fquantile(n2-1,n1-1,1-(0.5*alpha));
	lower_bound = divisiond(s1*s1,s2*s2)*(divisiond(1,F_quantile1));
	upper_bound = divisiond(s1*s1,s2*s2)*(F_quantile2);

	cout << "\n********************************************************" << endl;
	cout << "\nConfidence interval for the ratio of two variances" << endl;
	cout << "Mean data 1 = "<< mean1 << ", s1 = " << s1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << s2 << ", n2 = " << n2 << endl;	
	cout << "\nThe " << 100*(1-alpha) << "\% confidence interval is :\n" << endl;	
	cout << lower_bound << " <  σ1^{2} /  σ2^{2} < " << upper_bound << endl;	
	cout << sqrt(lower_bound) << " <  σ1 /  σ2 < " << sqrt(upper_bound) << endl;	
	cout << "\n********************************************************" << endl;
}

void hypothesistest(vector<double> data, double mu, double sigma, double alpha, double effect_size)
{
	double critical_value, critical_value2, z_quantile, z_value, z_beta, meanbeta, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n = data.size();
	/*double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	double sigma_sample = sqrt(variance);*/

	std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing" << endl;
	cout << "\nH0: μ = " << mu << endl;
	cout << "H1: μ != " << mu << endl;
	
	if (n < 30) // For small sample size we use t-statistic
	{
		t_value = ((mu-effect_size)-mu)/(sigma/sqrt(n));

		t_quantile = tquantile(1,n-1,1-(0.5*alpha));

		critical_value = mu + abs(t_quantile)*(sigma)/(sqrt(n)) ;
		critical_value2 = mu - abs(t_quantile)*(sigma)/(sqrt(n)) ;
			
		meanbeta = mu + effect_size; // for two-tailed
		t_beta = (critical_value-meanbeta)/(sigma/sqrt(n));
		if (t_beta > 0)
		{
			t_beta = -t_beta;
		}
		else if (t_beta < 0)
		{
			t_beta = t_beta;
		}
		beta = (1 - 2*tcdf(t_beta,n-1) );
		power = 1-beta;
		p_value = 2*tcdf(t_value,n-1);

		cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
		cout << "The probability of a type II error (beta): " << beta << endl;
		cout << "Power: " << power << endl;
		cout << "\nCritical value (upper bound): " << critical_value << endl;
		cout << "Critical value (lower bound): " << critical_value2 << endl;
		cout << "Critical region : t < - " << t_quantile << " and t > " << t_quantile << endl;
		cout << "Computed t : " << t_value << endl;
		cout << "\n*Reject the null hypothesis when the sample average is greater than " << critical_value << " or less than " << critical_value2 << endl;
		cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << " or less than -" << t_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
	}
	else if (n >= 30) // For large sample size we use Z-statistic
	{
		z_value = ((mu-effect_size)-mu)/(sigma/sqrt(n));
		
		z_quantile = zquantile(1-(0.5*alpha));

		critical_value = mu + abs(z_quantile)*(sigma)/(sqrt(n)) ; // upper bound
		critical_value2 = mu - abs(z_quantile)*(sigma)/(sqrt(n)) ; // lower bound

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
		cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << " or less than -" << z_quantile << endl;
		
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
	/*double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	double sigma_sample = sqrt(variance);*/

	std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing" << endl;
	cout << "\nH0: μ = " << mu << endl;
	cout << "H1: μ > " << mu << endl;
	
	if (n < 30) // For small sample size we use t-statistic
	{
		t_value = ((mu+effect_size)-mu)/(sigma/sqrt(n));

		t_quantile = tquantile(1,n-1,1-alpha);

		critical_value = mu + abs(t_quantile)*(sigma)/(sqrt(n)) ;

		meanbeta = mu + effect_size; // for right-tailed
		t_beta = (critical_value-meanbeta)/(sigma/sqrt(n));
		if (t_beta > 0)
		{
			t_beta = t_beta;
		}
		else if (t_beta < 0)
		{
			t_beta = -t_beta;
		}

		beta = tcdf(t_beta,n-1);
		power = 1-beta;
		p_value = 1 - tcdf(t_value,n-1);

		cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
		cout << "The probability of a type II error (beta): " << beta << endl;
		cout << "Power: " << power << endl;
		cout << "\nCritical value: " << critical_value << endl;
		cout << "Critical region : t > " << t_quantile << endl;
		cout << "Computed t : " << t_value << endl;
		cout << "\n*Reject the null hypothesis when the sample average is greater than " << critical_value << endl;
		cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
	}
	else if (n >= 30) // For large sample size we use Z-statistic
	{
		z_value = ((mu+effect_size)-mu)/(sigma/sqrt(n));
		
		z_quantile = zquantile(1-alpha);

		critical_value = mu + abs(z_quantile)*(sigma)/(sqrt(n)) ;

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

void hypothesistest_lefttailed(vector<double> data, double mu, double sigma, double alpha, double effect_size)
{
	double critical_value, z_quantile, z_value, z_beta, meanbeta, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n = data.size();
	/*double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	
	double variance = sumSquaredDiff / (n - 1); // Sample variance
	double sigma_sample = sqrt(variance);*/

	std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed Hypothesis testing" << endl;
	cout << "\nH0: μ = " << mu << endl;
	cout << "H1: μ < " << mu << endl;
	
	if (n < 30) // For small sample size we use t-statistic
	{
		t_value = ((mu-effect_size)-mu)/(sigma/sqrt(n));

		t_quantile = tquantile(1,n-1,alpha);

		critical_value = mu - abs(t_quantile)*(sigma)/(sqrt(n)) ;

		meanbeta = mu - effect_size; // for left-tailed
		t_beta = (critical_value-meanbeta)/(sigma/sqrt(n));
		if (t_beta > 0)
		{
			t_beta = -t_beta;
		}
		else if (t_beta < 0)
		{
			t_beta = t_beta;
		}

		beta = 1 - tcdf(t_beta,n-1);
		power = 1-beta;
		p_value = tcdf(t_value,n-1);
		
		cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
		cout << "The probability of a type II error (beta): " << beta << endl;
		cout << "Power: " << power << endl;
		cout << "\nCritical value: " << critical_value << endl;
		cout << "Critical region : t < " << t_quantile << endl;
		cout << "Computed t : " << t_value << endl;
		cout << "\n*Reject the null hypothesis when the sample average is less than " << critical_value << endl;
		cout << "\n*Reject the null hypothesis when the computed t is less than " << t_quantile << endl;

		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
		
	}
	else if (n >= 30) // For large sample size we use Z-statistic
	{
		z_value = ((mu-effect_size)-mu)/(sigma/sqrt(n));
		
		z_quantile = zquantile(alpha);

		critical_value = mu - abs(z_quantile)*(sigma)/(sqrt(n)) ;

		meanbeta = mu - effect_size; // for right-tailed
		z_beta = (critical_value-meanbeta)/(sigma/sqrt(n));

		beta = 1 - normalcdf(z_beta,0,1);	
		power = 1-beta;
		p_value = normalcdf(z_value,0,1);	
		
		cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
		cout << "The probability of a type II error (beta): " << beta << endl;
		cout << "Power: " << power << endl;
		cout << "\nCritical value: " << critical_value << endl;
		cout << "Critical region : z < " << z_quantile << endl;
		cout << "Computed z : " << z_value << endl;
		cout << "\n*Reject the null hypothesis when the sample average is less than " << critical_value << endl;
		cout << "\n*Reject the null hypothesis when the computed z is less than " << z_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
	
	}
}

void hypothesistest_knownvariances_righttailed(vector<double> data1, vector<double> data2, double d0, double sigma1, double sigma2, double alpha, double delta)
{
	double critical_value, p_value, power, beta;
	double z_value, z_quantile, z_beta;

	int n1 = data1.size();
	int n2 = data2.size();

	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing, known variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 > " << d0 << endl;
	
	// Computation part
	z_value = ((mean1 - mean2)-d0)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	z_quantile = zquantile(1-alpha);

	critical_value = d0 + abs(z_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula

	z_beta = abs(z_quantile) - (delta)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));
	beta = normalcdf(-abs(z_beta),0,1);
	power = 1-beta;
	p_value = 1 - normalcdf(z_value,0,1);

	cout << "\nMean data 1 = "<< mean1 << ", σ1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", σ2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value: " << critical_value << endl;
	cout << "Critical region : z > " << z_quantile << endl;
	cout << "Computed z : " << z_value << endl;

	cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_knownvariances_lefttailed(vector<double> data1, vector<double> data2, double d0, double sigma1, double sigma2, double alpha, double delta)
{
	double critical_value, p_value, power, beta;
	double z_value, z_quantile, z_beta;

	int n1 = data1.size();
	int n2 = data2.size();

	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing, known variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 < " << d0 << endl;
	
	// Computation part
	z_value = ((mean1 - mean2)-d0)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	z_quantile = zquantile(alpha);

	critical_value = d0 - abs(z_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula

	z_beta = z_quantile + (delta)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));
	beta = normalcdf(-abs(z_beta),0,1);
	power = 1-beta;
	p_value = normalcdf(z_value,0,1);

	cout << "\nMean data 1 = "<< mean1 << ", σ1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", σ2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value: " << critical_value << endl;
	cout << "Critical region : z < " << z_quantile << endl;
	cout << "Computed z : " << z_value << endl;

	cout << "\n*Reject the null hypothesis when the computed z is less than " << z_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_knownvariances_twotailed(vector<double> data1, vector<double> data2, double d0, double sigma1, double sigma2, double alpha, double delta)
{
	double critical_value, critical_value2, p_value, power, beta;
	double z_value, z_quantile, z_beta;

	int n1 = data1.size();
	int n2 = data2.size();

	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing, known variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 != " << d0 << endl;
	
	// Computation part
	z_value = ((mean1 - mean2)-d0)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	z_quantile = zquantile(1-0.5*alpha);

	critical_value = d0 + abs(z_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula
	critical_value2 = d0 - abs(z_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula

	z_beta = abs(z_quantile) - (delta)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));
	//z_beta2 = -abs(z_quantile) + (delta)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2))); // the result is - z_beta
	beta = 2*normalcdf(-abs(z_beta),0,1) ;//   + (1-normalcdf(z_beta2,0,1));
	power = 1-beta;
	p_value = 2*normalcdf(-abs(z_value),0,1);

	cout << "\nMean data 1 = "<< mean1 << ", σ1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", σ2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value (upper bound): " << critical_value << endl;
	cout << "Critical value (lower bound): " << critical_value2 << endl;
	cout << "Critical region : z < - " << z_quantile << " and z > " << z_quantile << endl;
	cout << "Computed z : " << z_value << endl;
	cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << " or less than -" << z_quantile << endl;
				
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_equalunknownvariances_righttailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n1 = data1.size();
	int n2 = data2.size();
	int nu = n1 + n2 - 2;
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing, equal but unknown variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 > " << d0 << endl;
	
	// Computation part
	double sp = sqrt(divisiond((n1-1)*(variance1) + (n2-1)*(variance2),nu));
	t_value = ((mean1 - mean2)-d0)/(sp * sqrt(divisiond(1,n1) + divisiond(1,n2)));

	t_quantile = tquantile(1,n1+n2-2,1-alpha);

	critical_value = d0 + abs(t_quantile)*(sp)/(sqrt(n1+n2-2)) ; // my own formula
	
	t_beta = abs(t_quantile) - divisiond( delta , sp * sqrt(divisiond(1,n1) + divisiond(1,n2)) );
	
	beta = tcdf(-abs(t_beta),n1+n2-2);
	power = 1-beta;
	p_value = 1 - tcdf(t_value,n1+n2-2);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sqrt(variance1) << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sqrt(variance2) << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value: " << critical_value << endl;
	cout << "Critical region : t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_equalunknownvariances_lefttailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n1 = data1.size();
	int n2 = data2.size();
	int nu = n1 + n2 - 2;
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing, equal but unknown variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 < " << d0 << endl;
	
	// Computation part
	double sp = sqrt(divisiond((n1-1)*(variance1) + (n2-1)*(variance2),nu));
	t_value = ((mean1 - mean2)-d0)/(sp * sqrt(divisiond(1,n1) + divisiond(1,n2)));

	t_quantile = tquantile(1,n1+n2-2,alpha);

	critical_value = d0 - abs(t_quantile)*(sp)/(sqrt(n1+n2-2)) ; // my own formula

	t_beta = t_quantile + divisiond( delta , sp * sqrt(divisiond(1,n1) + divisiond(1,n2)) );
	
	beta = tcdf(-abs(t_beta),n1+n2-2);
	power = 1-beta;
	p_value = tcdf(t_value,n1+n2-2);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sqrt(variance1) << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sqrt(variance2) << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value: " << critical_value << endl;
	cout << "Critical region : t < " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is less than " << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_equalunknownvariances_twotailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, critical_value2, p_value, power, beta;
	double t_value, t_quantile, t_beta;

	int n1 = data1.size();
	int n2 = data2.size();
	int nu = n1 + n2 - 2;
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing, equal but unknown variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 != " << d0 << endl;
	
	// Computation part
	double sp = sqrt(divisiond((n1-1)*(variance1) + (n2-1)*(variance2),nu));
	t_value = ((mean1 - mean2)-d0)/(sp * sqrt(divisiond(1,n1) + divisiond(1,n2)));

	t_quantile = tquantile(1,n1+n2-2,1-(0.5*alpha));

	critical_value = d0 + abs(t_quantile)*(sp)/(sqrt(n1+n2-2)) ; // my own formula
	critical_value2 = d0 - abs(t_quantile)*(sp)/(sqrt(n1+n2-2)) ; // my own formula

	t_beta = abs(t_quantile) - divisiond( delta , sp * sqrt(divisiond(1,n1) + divisiond(1,n2)) );
	
	beta = 2*tcdf(-abs(t_beta),n1+n2-2);
	power = 1-beta;
	p_value = 2*tcdf(-abs(t_value),n1+n2-2);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sqrt(variance1) << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sqrt(variance2) << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value (upper bound): " << critical_value << endl;
	cout << "Critical value (lower bound): " << critical_value2 << endl;
	cout << "Critical region : t < - " << t_quantile << " and t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << " or less than -" << t_quantile << endl;
				
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_unequalunknownvariances_righttailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, p_value, power, beta;
	double t_value, t_quantile, t_beta, sigma1, sigma2;

	int n1 = data1.size();
	int n2 = data2.size();
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	sigma1 = sqrt(variance1);
	sigma2 = sqrt(variance2);
	
	int nu = divisiond(pow(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2),2), divisiond(pow(divisiond(sigma1*sigma1,n1),2),n1-1) + divisiond(pow(divisiond(sigma2*sigma2,n2),2),n2-1) );
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing, unequal and unknown variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 > " << d0 << endl;
	
	// Computation part
	t_value = ((mean1 - mean2)-d0)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	t_quantile = tquantile(1,nu,1-alpha);

	critical_value = d0 + abs(t_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula

	t_beta = abs(t_quantile) - divisiond(delta, sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));
	beta = tcdf(-abs(t_beta),nu);
	power = 1-beta;
	p_value = 1 - tcdf(t_value,nu);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value: " << critical_value << endl;
	cout << "Critical region : t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_unequalunknownvariances_lefttailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, p_value, power, beta;
	double t_value, t_quantile, t_beta, sigma1, sigma2;

	int n1 = data1.size();
	int n2 = data2.size();
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	sigma1 = sqrt(variance1);
	sigma2 = sqrt(variance2);
	
	int nu = divisiond(pow(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2),2), divisiond(pow(divisiond(sigma1*sigma1,n1),2),n1-1) + divisiond(pow(divisiond(sigma2*sigma2,n2),2),n2-1) );
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing, unequal and unknown variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 < " << d0 << endl;
	
	// Computation part
	t_value = ((mean1 - mean2)-d0)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	t_quantile = tquantile(1,nu,alpha);

	critical_value = d0 - abs(t_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula
	
	t_beta = t_quantile + divisiond(delta, sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)) );
	beta = tcdf(-abs(t_beta),nu);
	power = 1-beta;
	p_value = tcdf(t_value,nu);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value: " << critical_value << endl;
	cout << "Critical region : t < " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is less than " << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_unequalunknownvariances_twotailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, critical_value2, p_value, power, beta;
	double t_value, t_quantile, t_beta, sigma1, sigma2;

	int n1 = data1.size();
	int n2 = data2.size();
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	sigma1 = sqrt(variance1);
	sigma2 = sqrt(variance2);
	
	int nu = divisiond(pow(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2),2), divisiond(pow(divisiond(sigma1*sigma1,n1),2),n1-1) + divisiond(pow(divisiond(sigma2*sigma2,n2),2),n2-1) );
	
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing, unequal and unknown variances" << endl;
	cout << "\nH0: μ1 - μ2 = " << d0 << endl;
	cout << "H1: μ1 - μ2 != " << d0 << endl;
	
	// Computation part
	t_value = ((mean1 - mean2)-d0)/(sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)));

	t_quantile = tquantile(1,nu,1-(0.5*alpha));

	critical_value = d0 + abs(t_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula
	critical_value2 = d0 - abs(t_quantile)*(0.5*( (sigma1)/(sqrt(n1)) + (sigma2)/(sqrt(n2)) ) ); // my own formula
	
	t_beta = abs(t_quantile) - divisiond(delta, sqrt(divisiond(sigma1*sigma1,n1) + divisiond(sigma2*sigma2,n2)) );
	beta = 2*tcdf(-abs(t_beta),nu);
	power = 1-beta;
	p_value = 2*tcdf(-abs(t_value),nu);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value (upper bound): " << critical_value << endl;
	cout << "Critical value (lower bound): " << critical_value2 << endl;
	cout << "Critical region : t < - " << t_quantile << " and t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << " or less than -" << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_paired_righttailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, p_value, power, beta, sigmad ;
	double t_value, t_quantile, t_beta, sigma1, sigma2;
	vector<double> datad;
	int n1 = data1.size();
	int n2 = data2.size();
	if (n1 != n2)
	{
		cerr << "Error: The samples size are different." << endl;
	}
	for (int i = 0; i <n1 ; ++i)
	{
		datad.push_back(data1[i] - data2[i]);
	}
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double meand = calculateMean(datad);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	double sumSquaredDiffd = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	for (double val : datad) 
	{
		sumSquaredDiffd += std::pow(val - meand, 2);
	}
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	double varianced = sumSquaredDiffd / (n1 - 1);
	sigma1 = sqrt(variance1);
	sigma2 = sqrt(variance2);
	sigmad = sqrt(varianced);

	int nu = n1-1;
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nRight tailed hypothesis testing, paired observations" << endl;
	cout << "\nH0: μd = μ1 - μ2 = " << d0 << endl;
	cout << "H1: μd = μ1 - μ2 > " << d0 << endl;
	
	// Computation part
	t_value = (meand - d0)/(divisiond(sigmad,sqrt(n1)));

	t_quantile = tquantile(1,nu,1-alpha);

	critical_value = d0 + abs(t_quantile)*divisiond(sigmad,sqrt(n1)); // my own formula
	
	t_beta = abs(t_quantile) - divisiond(delta*sqrt(n1),sigmad);
	beta = tcdf(-abs(t_beta),nu);
	power = 1-beta;
	p_value = 1 - tcdf(t_value,nu);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "Mean d = "<< meand << ", sd = " << sigmad << ", nd = " << n1 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value (upper bound): " << critical_value << endl;
	cout << "Critical region : t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_paired_lefttailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, p_value, power, beta, sigmad ;
	double t_value, t_quantile, t_beta, sigma1, sigma2;
	vector<double> datad;
	int n1 = data1.size();
	int n2 = data2.size();
	if (n1 != n2)
	{
		cerr << "Error: The samples size are different." << endl;
	}
	for (int i = 0; i <n1 ; ++i)
	{
		datad.push_back(data1[i] - data2[i]);
	}
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double meand = calculateMean(datad);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	double sumSquaredDiffd = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	for (double val : datad) 
	{
		sumSquaredDiffd += std::pow(val - meand, 2);
	}
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	double varianced = sumSquaredDiffd / (n1 - 1);
	sigma1 = sqrt(variance1);
	sigma2 = sqrt(variance2);
	sigmad = sqrt(varianced);

	int nu = n1-1;
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing, paired observations" << endl;
	cout << "\nH0: μd = μ1 - μ2 = " << d0 << endl;
	cout << "H1: μd = μ1 - μ2 < " << d0 << endl;
	
	// Computation part
	t_value = (meand - d0)/(divisiond(sigmad,sqrt(n1)));

	t_quantile = tquantile(1,nu,alpha);

	critical_value = d0 - abs(t_quantile)*divisiond(sigmad,sqrt(n1)); // my own formula
	
	t_beta = abs(t_quantile) + divisiond(delta*sqrt(n1),sigmad);
	beta = tcdf(-abs(t_beta),nu);
	power = 1-beta;
	p_value = tcdf(t_value,nu);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "Mean d = "<< meand << ", sd = " << sigmad << ", nd = " << n1 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value (upper bound): " << critical_value << endl;
	cout << "Critical region : t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_paired_twotailed(vector<double> data1, vector<double> data2, double d0, double alpha, double delta)
{
	double critical_value, critical_value2, p_value, power, beta, sigmad ;
	double t_value, t_quantile, t_beta, sigma1, sigma2;
	vector<double> datad;
	int n1 = data1.size();
	int n2 = data2.size();
	if (n1 != n2)
	{
		cerr << "Error: The samples size are different." << endl;
	}
	for (int i = 0; i <n1 ; ++i)
	{
		datad.push_back(data1[i] - data2[i]);
	}
	double mean1 = calculateMean(data1);
	double mean2 = calculateMean(data2);
	double meand = calculateMean(datad);
	double sumSquaredDiff1 = 0.0;
	double sumSquaredDiff2 = 0.0;
	double sumSquaredDiffd = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	for (double val : datad) 
	{
		sumSquaredDiffd += std::pow(val - meand, 2);
	}
	double variance1 = sumSquaredDiff1 / (n1 - 1);
	double variance2 = sumSquaredDiff2 / (n2 - 1);
	double varianced = sumSquaredDiffd / (n1 - 1);
	sigma1 = sqrt(variance1);
	sigma2 = sqrt(variance2);
	sigmad = sqrt(varianced);

	int nu = n1-1;
	//std::sort(data.begin(), data.end());

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing, paired observations" << endl;
	cout << "\nH0: μd = μ1 - μ2 = " << d0 << endl;
	cout << "H1: μd = μ1 - μ2 != " << d0 << endl;
	
	// Computation part
	t_value = (meand - d0)/(divisiond(sigmad,sqrt(n1)));

	t_quantile = tquantile(1,nu,1-(0.5*alpha));

	critical_value = d0 + abs(t_quantile)*divisiond(sigmad,sqrt(n1)); // my own formula
	critical_value2 = d0 - abs(t_quantile)*divisiond(sigmad,sqrt(n1)); // my own formula
	
	t_beta = abs(t_quantile) - divisiond(delta*sqrt(n1),sigmad);
	beta = 2*tcdf(-abs(t_beta),nu);
	power = 1-beta;
	p_value = 2*tcdf(-abs(t_value),nu);

	cout << "\nMean data 1 = "<< mean1 << ", s1 = " << sigma1 << ", n1 = " << n1 << endl;	
	cout << "Mean data 2 = "<< mean2 << ", s2 = " << sigma2 << ", n2 = " << n2 << endl;	
	cout << "Mean d = "<< meand << ", sd = " << sigmad << ", nd = " << n1 << endl;	
	cout << "\nThe probability of a type I error (alpha): " << alpha << endl;
	cout << "The probability of a type II error (beta): " << beta << endl;
	cout << "Power: " << power << endl;
	cout << "\nCritical value (upper bound): " << critical_value << endl;
	cout << "Critical value (lower bound): " << critical_value2 << endl;
	cout << "Critical region : t < - " << t_quantile << " and t > " << t_quantile << endl;
	cout << "Computed t : " << t_value << endl;

	cout << "\n*Reject the null hypothesis when the computed t is greater than " << t_quantile << " or less than -" << t_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
}

void hypothesistest_proportion_righttailed(vector<string> data, const string &inputString, double p0, double alpha)
{
	double p_value, p, q, q0;
	double z_value, z_quantile;

	int x_count = std::count(data.begin(), data.end(), inputString);
	int n = data.size();
	p = divisiond(x_count,n);
	q = 1-p;
	q0 = 1 - p0;

	if(n < 30)
	{
		cout << "\n********************************************************" << endl;
		cout << "\nRight-tailed hypothesis testing on a single proportion" << endl;
		cout << "\nH0: p = " << p0 << endl;
		cout << "H1: p > " << p0 << endl;
		cout << "Data size =  "<< n << endl;	
		cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
		cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;	

		p_value = 1 - binomialcdf(x_count,n,p0);
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;	
	}
	else if(n >= 30)
	{
		cout << "\n********************************************************" << endl;
		cout << "\nRight-tailed hypothesis testing on a single proportion" << endl;
		cout << "\nH0: p = " << p0 << endl;
		cout << "H1: p > " << p0 << endl;
		cout << "Data size =  "<< n << endl;	
		cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
		cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;	

		z_quantile = zquantile(1-alpha);
		z_value = divisiond(p - p0,sqrt(divisiond(p0*q0,n))) ;
		p_value = 1 - normalcdf(z_value,0,1);

		cout << "Critical region : z > " << z_quantile << endl;
		cout << "Computed z : " << z_value << endl;
		cout << "\n********************************************************" << endl;

		cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;
	}
}

void hypothesistest_proportion_lefttailed(vector<string> data, const string &inputString, double p0, double alpha)
{
	double p_value, p, q, q0;
	double z_value, z_quantile;

	int x_count = std::count(data.begin(), data.end(), inputString);
	int n = data.size();
	p = divisiond(x_count,n);
	q = 1-p;
	q0 = 1 - p0;

	if(n < 30)
	{
		cout << "\n********************************************************" << endl;
		cout << "\nLeft-tailed hypothesis testing on a single proportion" << endl;
		cout << "\nH0: p = " << p0 << endl;
		cout << "H1: p < " << p0 << endl;
		cout << "Data size =  "<< n << endl;	
		cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
		cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;	
		p_value = binomialcdf(x_count,n,p0);
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;	
	}
	else if(n >=30)
	{
		cout << "\n********************************************************" << endl;
		cout << "\nLeft-tailed hypothesis testing on a single proportion" << endl;
		cout << "\nH0: p = " << p0 << endl;
		cout << "H1: p < " << p0 << endl;
		cout << "Data size =  "<< n << endl;	
		cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
		cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;	

		z_quantile = zquantile(alpha);
		z_value = divisiond(p - p0,sqrt(divisiond(p0*q0,n))) ;
		p_value = normalcdf(z_value,0,1);

		cout << "Critical region : z < " << z_quantile << endl;
		cout << "Computed z : " << z_value << endl;
		cout << "\n********************************************************" << endl;

		cout << "\n*Reject the null hypothesis when the computed z is less than " << z_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;	
	}
}

void hypothesistest_proportion_twotailed(vector<string> data, const string &inputString, double p0, double alpha)
{
	double p_value, p, q, q0;
	double z_value, z_quantile;

	int x_count = std::count(data.begin(), data.end(), inputString);
	int n = data.size();
	double x_compare = n*p0;
	p = divisiond(x_count,n);
	q = 1-p;
	q0 = 1 - p0;

	if(n < 30)
	{
		
		cout << "\n********************************************************" << endl;
		cout << "\nTwo-tailed hypothesis testing on a single proportion" << endl;
		cout << "\nH0: p = " << p0 << endl;
		cout << "H1: p != " << p0 << endl;
		cout << "Data size =  "<< n << endl;	
		cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
		cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;
		if(x_count < x_compare)
		{			
			p_value = 2*binomialcdf(x_count,n,p0);
		}
		else if(x_count >= x_compare)
		{
			p_value = 2*(1-binomialcdf(x_count,n,p0));
		}	
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;	
	}
	else if(n >= 30)
	{
		cout << "\n********************************************************" << endl;
		cout << "\nTwo-tailed hypothesis testing on a single proportion" << endl;
		cout << "\nH0: p = " << p0 << endl;
		cout << "H1: p != " << p0 << endl;
		cout << "Data size =  "<< n << endl;	
		cout << "Amount of "<< inputString << " in the data = " << x_count << endl;		
		cout << "Proportion for "<< inputString << ", p = " << p << ", q = " << q << endl;	

		z_quantile = zquantile(1-0.5*(alpha));
		z_value = divisiond(p - p0,sqrt(divisiond(p0*q0,n))) ;
		if(x_count < x_compare)
		{
			p_value = 2*normalcdf(z_value,0,1);
		}
		else if(x_count >= x_compare)
		{
			p_value = 2*(1-normalcdf(z_value,0,1));
		}

		cout << "Critical region : z < - " << z_quantile << " or > " << z_quantile << endl;
		cout << "Computed z : " << z_value << endl;
		cout << "\n********************************************************" << endl;

		cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << " or less than -" << z_quantile << endl;
		
		cout << "\nP-value: " << p_value << endl;
		cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
		cout << "\n********************************************************" << endl;	
	}
}

void hypothesistest_twoproportions_righttailed(vector<string> data1, vector<string> data2, const string &inputString, double alpha)
{
	double p_value, p1, p2, q1, q2, p_pooled, q_pooled;
	double z_value, z_quantile;

	int x1_count = std::count(data1.begin(), data1.end(), inputString);
	int n1 = data1.size();
	int x2_count = std::count(data2.begin(), data2.end(), inputString);
	int n2 = data2.size();
	p1 = divisiond(x1_count,n1);
	q1 = 1 - p1;
	p2 = divisiond(x2_count,n2);
	q2 = 1 - p2;
	p_pooled = divisiond(x1_count + x2_count, n1 + n2);
	q_pooled = 1 - p_pooled;

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing on two proportions" << endl;
	cout << "\nH0: p1 = p2" << endl;
	cout << "H1: p1 > p2" << endl;
	cout << "Data 1 size =  "<< n1 << ", Data 2 size =  "<< n2 << endl;	
	cout << "Amount of "<< inputString << " in the data 1 = " << x1_count << endl;		
	cout << "Amount of "<< inputString << " in the data 2 = " << x2_count << endl;		
	cout << "Proportion for "<< inputString << "\np1 = " << p1 << ", q1 = " << q1 << endl;	
	cout << "p2 = " << p2 << ", q2 = " << q2 << endl;	
	cout << "pooled estimate of p = " << p_pooled << ", pooled estimate of q = " << q_pooled << endl;	

	z_quantile = zquantile(1-alpha);
	z_value = divisiond(p1 - p2,sqrt(p_pooled*q_pooled*(divisiond(1,n1) + divisiond(1,n2)))) ;
	p_value = 1 - normalcdf(z_value,0,1);

	cout << "Critical region : z > " << z_quantile << endl;
	cout << "Computed z : " << z_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_twoproportions_lefttailed(vector<string> data1, vector<string> data2, const string &inputString, double alpha)
{
	double p_value, p1, p2, q1, q2, p_pooled, q_pooled;
	double z_value, z_quantile;

	int x1_count = std::count(data1.begin(), data1.end(), inputString);
	int n1 = data1.size();
	int x2_count = std::count(data2.begin(), data2.end(), inputString);
	int n2 = data2.size();
	p1 = divisiond(x1_count,n1);
	q1 = 1 - p1;
	p2 = divisiond(x2_count,n2);
	q2 = 1 - p2;
	p_pooled = divisiond(x1_count + x2_count, n1 + n2);
	q_pooled = 1 - p_pooled;

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing on two proportions" << endl;
	cout << "\nH0: p1 = p2" << endl;
	cout << "H1: p1 < p2" << endl;
	cout << "Data 1 size =  "<< n1 << ", Data 2 size =  "<< n2 << endl;	
	cout << "Amount of "<< inputString << " in the data 1 = " << x1_count << endl;		
	cout << "Amount of "<< inputString << " in the data 2 = " << x2_count << endl;		
	cout << "Proportion for "<< inputString << "\np1 = " << p1 << ", q1 = " << q1 << endl;	
	cout << "p2 = " << p2 << ", q2 = " << q2 << endl;	
	cout << "pooled estimate of p = " << p_pooled << ", pooled estimate of q = " << q_pooled << endl;	

	z_quantile = zquantile(alpha);
	z_value = divisiond(p1 - p2,sqrt(p_pooled*q_pooled*(divisiond(1,n1) + divisiond(1,n2)))) ;
	p_value = normalcdf(z_value,0,1);

	cout << "Critical region : z < " << z_quantile << endl;
	cout << "Computed z : " << z_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed z is less than " << z_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_twoproportions_twotailed(vector<string> data1, vector<string> data2, const string &inputString, double alpha)
{
	double p_value, p1, p2, q1, q2, p_pooled, q_pooled;
	double z_value, z_quantile;

	int x1_count = std::count(data1.begin(), data1.end(), inputString);
	int n1 = data1.size();
	int x2_count = std::count(data2.begin(), data2.end(), inputString);
	int n2 = data2.size();
	p1 = divisiond(x1_count,n1);
	q1 = 1 - p1;
	p2 = divisiond(x2_count,n2);
	q2 = 1 - p2;
	p_pooled = divisiond(x1_count + x2_count, n1 + n2);
	q_pooled = 1 - p_pooled;

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing on two proportions" << endl;
	cout << "\nH0: p1 = p2" << endl;
	cout << "H1: p1 != p2" << endl;
	cout << "Data 1 size =  "<< n1 << ", Data 2 size =  "<< n2 << endl;	
	cout << "Amount of "<< inputString << " in the data 1 = " << x1_count << endl;		
	cout << "Amount of "<< inputString << " in the data 2 = " << x2_count << endl;		
	cout << "Proportion for "<< inputString << "\np1 = " << p1 << ", q1 = " << q1 << endl;	
	cout << "p2 = " << p2 << ", q2 = " << q2 << endl;	
	cout << "pooled estimate of p = " << p_pooled << ", pooled estimate of q = " << q_pooled << endl;	

	z_quantile = zquantile(1-(0.5*alpha));
	z_value = divisiond(p1 - p2,sqrt(p_pooled*q_pooled*(divisiond(1,n1) + divisiond(1,n2)))) ;
	p_value = 2*normalcdf(-abs(z_value),0,1);

	cout << "Critical region : z < -" << z_quantile << " or > " << z_quantile << endl;
	cout << "Computed z : " << z_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed z is greater than " << z_quantile << " or less than -" << z_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_onesamplevariance_righttailed(vector<double> data, double sigma, double alpha)
{
	double p_value;
	double chi_value, chi_quantile;
	int n = data.size();
	int nu = n-1;
	double variance = sigma*sigma;

	double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance_sample = sumSquaredDiff / (n - 1);

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing on one-sample variance" << endl;
	cout << "\nH0: σ^{2} = " << sigma*sigma << endl;
	cout << "H1: σ^{2} > " << sigma*sigma << endl;
	cout <<  "Mean data = " << mean << ", s^{2} = " << variance_sample << ", n =  "<< n << endl;	

	chi_quantile = chisquaredquantile(nu,1-alpha);
	chi_value = divisiond(nu*variance_sample,variance) ;
	p_value = 1 - chisquaredcdf(chi_value,nu);

	cout << "Critical region : χ^{2} > " << chi_quantile << endl;
	cout << "Computed χ^{2} : " << chi_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed χ^{2} is greater than " << chi_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_onesamplevariance_lefttailed(vector<double> data, double sigma, double alpha)
{
	double p_value;
	double chi_value, chi_quantile;
	int n = data.size();
	int nu = n-1;
	double variance = sigma*sigma;

	double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance_sample = sumSquaredDiff / (n - 1);

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing on one-sample variance" << endl;
	cout << "\nH0: σ^{2} = " << sigma*sigma << endl;
	cout << "H1: σ^{2} < " << sigma*sigma << endl;
	cout <<  "Mean data = " << mean << ", s^{2} = " << variance_sample << ", n =  "<< n << endl;	

	chi_quantile = chisquaredquantile(nu,alpha);
	chi_value = divisiond(nu*variance_sample,variance) ;
	p_value = chisquaredcdf(chi_value,nu);

	cout << "Critical region : χ^{2} < " << chi_quantile << endl;
	cout << "Computed χ^{2} : " << chi_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed χ^{2} is less than " << chi_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_onesamplevariance_twotailed(vector<double> data, double sigma, double alpha)
{
	double p_value;
	double chi_value, chi_quantile, chi_quantile2;
	int n = data.size();
	int nu = n-1;
	double variance = sigma*sigma;

	double mean = calculateMean(data);
	double sumSquaredDiff = 0.0;
	for (double val : data) 
	{
		sumSquaredDiff += std::pow(val - mean, 2);
	}
	double variance_sample = sumSquaredDiff / (n - 1);

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing on one-sample variance" << endl;
	cout << "\nH0: σ^{2} = " << sigma*sigma << endl;
	cout << "H1: σ^{2} != " << sigma*sigma << endl;
	cout <<  "Mean data = " << mean << ", s^{2} = " << variance_sample << ", n =  "<< n << endl;	

	chi_quantile = chisquaredquantile(nu,0.5*alpha);
	chi_quantile2 = chisquaredquantile(nu,1 - (0.5*alpha));
	chi_value = divisiond(nu*variance_sample,variance) ;
	p_value = chisquaredcdf(chi_value,nu);

	cout << "Critical region : χ^{2} < " << chi_quantile << " or  χ^{2} > " << chi_quantile2 << endl;
	cout << "Computed χ^{2} : " << chi_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed χ^{2} is less than " << chi_quantile << " or  computed χ^{2} is greater than " << chi_quantile2 << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_twosamplesvariances_righttailed(vector<double> data1,  vector<double> data2, double alpha)
{
	double p_value;
	double F_value, F_quantile;
	int n1 = data1.size();
	int nu1 = n1-1;
	int n2 = data2.size();
	int nu2 = n2-1;

	double mean1 = calculateMean(data1);
	double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	double variance_sample1 = sumSquaredDiff1 / (n1 - 1);

	double mean2 = calculateMean(data2);
	double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance_sample2 = sumSquaredDiff2 / (n2 - 1);

	cout << "\n********************************************************" << endl;
	cout << "\nRight-tailed hypothesis testing on two-samples variances" << endl;
	cout << "\nH0: σ1^{2} = σ2^{2}" << endl;
	cout << "H1: σ1^{2} > σ2^{2}" << endl;
	cout <<  "Mean data 1 = " << mean1 << ", s1^{2} = " << variance_sample1 << ", n1 =  "<< n1 << endl;	
	cout <<  "Mean data 2 = " << mean2 << ", s2^{2} = " << variance_sample2 << ", n2 =  "<< n2 << endl;	

	F_quantile = Fquantile(nu1,nu2,1-alpha);
	F_value = divisiond(variance_sample1,variance_sample2) ;
	p_value = 1 - Fcdf(F_value,nu1, nu2);

	cout << "Critical region : f > " << F_quantile << endl;
	cout << "Computed f : " << F_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed f is greater than " << F_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_twosamplesvariances_lefttailed(vector<double> data1,  vector<double> data2, double alpha)
{
	double p_value;
	double F_value, F_quantile;
	int n1 = data1.size();
	int nu1 = n1-1;
	int n2 = data2.size();
	int nu2 = n2-1;

	double mean1 = calculateMean(data1);
	double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	double variance_sample1 = sumSquaredDiff1 / (n1 - 1);

	double mean2 = calculateMean(data2);
	double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance_sample2 = sumSquaredDiff2 / (n2 - 1);

	cout << "\n********************************************************" << endl;
	cout << "\nLeft-tailed hypothesis testing on two-samples variances" << endl;
	cout << "\nH0: σ1^{2} = σ2^{2}" << endl;
	cout << "H1: σ1^{2} < σ2^{2}" << endl;
	cout <<  "Mean data 1 = " << mean1 << ", s1^{2} = " << variance_sample1 << ", n1 =  "<< n1 << endl;	
	cout <<  "Mean data 2 = " << mean2 << ", s2^{2} = " << variance_sample2 << ", n2 =  "<< n2 << endl;	

	F_quantile = Fquantile(nu1,nu2,alpha);
	F_value = divisiond(variance_sample1,variance_sample2) ;
	p_value = Fcdf(F_value,nu1, nu2);

	cout << "Critical region : f < " << F_quantile << endl;
	cout << "Computed f : " << F_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed f is less than " << F_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_twosamplesvariances_twotailed(vector<double> data1,  vector<double> data2, double alpha)
{
	double p_value;
	double F_value, F_quantile, F_quantile2;
	int n1 = data1.size();
	int nu1 = n1-1;
	int n2 = data2.size();
	int nu2 = n2-1;

	double mean1 = calculateMean(data1);
	double sumSquaredDiff1 = 0.0;
	for (double val : data1) 
	{
		sumSquaredDiff1 += std::pow(val - mean1, 2);
	}
	double variance_sample1 = sumSquaredDiff1 / (n1 - 1);

	double mean2 = calculateMean(data2);
	double sumSquaredDiff2 = 0.0;
	for (double val : data2) 
	{
		sumSquaredDiff2 += std::pow(val - mean2, 2);
	}
	double variance_sample2 = sumSquaredDiff2 / (n2 - 1);

	cout << "\n********************************************************" << endl;
	cout << "\nTwo-tailed hypothesis testing on two-samples variances" << endl;
	cout << "\nH0: σ1^{2} = σ2^{2}" << endl;
	cout << "H1: σ1^{2} != σ2^{2}" << endl;
	cout <<  "Mean data 1 = " << mean1 << ", s1^{2} = " << variance_sample1 << ", n1 =  "<< n1 << endl;	
	cout <<  "Mean data 2 = " << mean2 << ", s2^{2} = " << variance_sample2 << ", n2 =  "<< n2 << endl;	

	F_quantile = Fquantile(nu1,nu2,1-(0.5*alpha));
	F_quantile2 = Fquantile(nu1,nu2,0.5*alpha);
	F_value = divisiond(variance_sample1,variance_sample2) ;
	p_value = 2*(1-Fcdf(F_value,nu1, nu2));

	cout << "Critical region : f > " << F_quantile << " and f < " << F_quantile2 << endl;
	cout << "Computed f : " << F_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed f is greater than " << F_quantile << " or less than " << F_quantile2 << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void hypothesistest_categoricaldata(vector<vector<double>> data, double alpha)
{
	double p_value, Total_sum;
	double chi_value, chi_quantile;
	int n = data.size();
	int C = data[0].size();
	int nu = (n-1)*(C-1);
	vector<double> total_column;
	vector<double> total_row;
	vector<vector<double>> expected_frequencies(n, vector<double>(C, 0.0));

	for(int i = 0 ; i < C ; ++i)
	{
		double sum = 0;
		for (int j = 0 ; j < n ; ++j)
		{
			sum += data[j][i]; // sum of the column
		}
			total_column.push_back(sum);
			cout << "sum of column " << i << " = " << total_column[i] << endl;	
	}
	for(int i = 0 ; i < n ; ++i)
	{
		double sum = 0;
		for (int j = 0 ; j < C ; ++j)
		{
			sum += data[i][j]; // sum of the column
		}
			total_row.push_back(sum);
			cout << "sum of row " << i << " = " << total_row[i] << endl;	
	}
	Total_sum = std::accumulate(total_row.begin(), total_row.end(), 0.0) ;
	
	for(int i = 0 ; i < n ; ++i)
	{
		for (int j = 0 ; j < C ; ++j)
		{
			expected_frequencies[i][j] = divisiond(total_row[i]*total_column[j],Total_sum); 
		}
	}
	for(int i = 0 ; i < n ; ++i)
	{
		for (int j = 0 ; j < C ; ++j)
		{
			chi_value += divisiond((data[i][j] - expected_frequencies[i][j])*(data[i][j] - expected_frequencies[i][j]), expected_frequencies[i][j]);
		}
	}
	chi_quantile = chisquaredquantile(nu,1 - alpha);
	p_value = 1 - chisquaredcdf(chi_value,nu);

	cout << "\n********************************************************" << endl;
	cout << "\nHypothesis testing on categorical data" << endl;
	cout << "\nH0: Independent" << endl;
	cout << "H1: Not independent" << endl;
	cout << "Number of rows =  "<< n << ", number of columns = " << C << endl;	
	cout << "Degrees of freedom =  "<< nu << endl;	
	
	cout << "\nObserved Frequencies:" << endl;
	printMatrix(data);
	cout << "\nExpected Frequencies:" << endl;
	printMatrix(expected_frequencies);
	
	cout << "Critical region : χ^{2} > " << chi_quantile << endl;
	cout << "Computed χ^{2} : " << chi_value << endl;
	cout << "\n********************************************************" << endl;

	cout << "\n*Reject the null hypothesis when the computed χ^{2} is greater than " << chi_quantile << endl;
		
	cout << "\nP-value: " << p_value << endl;
	cout << "\n*Reject null hypothesis if P-value <= "<< alpha << " , we fail to reject the null hypothesis if P-value > " << alpha << endl;
	cout << "\n********************************************************" << endl;
	
}

void countstring(vector<string> data, const string &inputString)
{
	int x_count = std::count(data.begin(), data.end(), inputString);
	
	cout << "The string \"" << inputString << "\" appears " << x_count << " times in the data." << endl;
}

void addstring(vector<string> &data, const string &inputString)
{
	data.push_back(inputString);
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



void regressionline(vector<vector<double>> matrix, double alpha)
{
	vector<double> x;
	vector<double> y;
	vector<double> y_fit;
	
	double x_mean, y_mean, b0,b1;
	double sum_xy, sum_x, sum_y, sum_xsquared, Sxy, Sxx, Syy;
	double SSE, SSR, SST, R_Squared, computed_f, p_value;
	double s, variance_estimate, b1_lowerconfidenceinterval, b1_upperconfidenceinterval, b0_lowerconfidenceinterval, b0_upperconfidenceinterval ;
	double fit, SE_fit, CI_lower, CI_upper, PI_lower, PI_upper;
	sum_xy = 0;
	sum_x = 0;
	sum_y = 0;
	sum_xsquared = 0;
	Sxx = 0;
	Syy = 0;
	Sxy = 0;
	int n = matrix.size(); // number of row
	double t_quantile = tquantile(1,n-2,1-(0.5*alpha));

	for(int i = 0; i < n; ++i)
	{
		x.push_back(matrix[i][0]);
		y.push_back(matrix[i][1]);
	}
	x_mean = calculateMean(x);
	y_mean = calculateMean(y);
	
	for(int i=0; i < n; ++i)
	{
		sum_x = sum_x + x[i];
		sum_y = sum_y + y[i];
		sum_xy = sum_xy +(x[i]*y[i]);
		sum_xsquared = sum_xsquared + (x[i]*x[i]);	
		Sxx = Sxx + (x[i] - x_mean)*(x[i] - x_mean);
		Syy = Syy + (y[i] - y_mean)*(y[i] - y_mean);
		Sxy = Sxy + (x[i] - x_mean)*(y[i] - y_mean);
	}	
	b1 = divisiond((n*sum_xy) - (sum_x*sum_y),(n*sum_xsquared)-(sum_x*sum_x));
	b0 = divisiond(sum_y-(b1*sum_x),n);
	SSE = Syy - (b1*Sxy);
	SSR = b1*Sxy;
	SST = SSR + SSE;
	R_Squared = divisiond(SSR,SST);
	variance_estimate = divisiond(SSE,n-2);
	s = sqrt(variance_estimate) ;
	computed_f = divisiond(SSR,variance_estimate);
	int df_model = 1;
	int df_total = n-1;
	int df_error = df_total - df_model;
	p_value = 1 - Fcdf(computed_f,df_model,df_error);
	
	for(int i = 0; i < n; ++i)
	{
		y_fit.push_back(b1*x[i] +b0);	
	}
	
	//for(int i = 0; i < n; ++i)
	//{
	//	SE_fit = SE_fit + sqrt(divisiond((y[i] - y_fit[i])*(y[i] - y_fit[i]),n-2));	
	//}
	
	b1_lowerconfidenceinterval = b1 - abs(t_quantile)*divisiond(sqrt(variance_estimate),sqrt(Sxx));
	b1_upperconfidenceinterval = b1 + abs(t_quantile)*divisiond(sqrt(variance_estimate),sqrt(Sxx));
	b0_lowerconfidenceinterval = b0 - abs(t_quantile)*divisiond(sqrt(variance_estimate),sqrt(n*Sxx))*sqrt(sum_xsquared);
	b0_upperconfidenceinterval = b0 + abs(t_quantile)*divisiond(sqrt(variance_estimate),sqrt(n*Sxx))*sqrt(sum_xsquared);

	cout << "\n********************************************************" << endl;
	cout << "\nSimple Linear Regression" << endl;
	cout << "\nSxx = " << Sxx << " , Syy = "<< Syy << " , Sxy = " << Sxy << endl;
	cout << "SSE = " << SSE << ", SSR = " << SSR << ", SST = " << SST << endl;
	cout << "R-Squared = " << 100*R_Squared << " %" << endl;
	cout << "Computed f = " << computed_f << ", P-value = " << p_value << endl;
	cout << "\nAn unbiased estimate of σ^{2} is : \ns^{2} = " << variance_estimate << endl; 	
	cout << "The estimated regression line is given by: \ny = " << b0 << " + " << b1 << "x" << endl;

	cout << "\n********************************************************" << endl;
	cout << "\nAnalysis of Variance" << endl;
	
	cout << "\nSource" << setw(10) << "\t\t" << "DF" << setw(10) << "\t" << "SS" << setw(10) << "\t\t" << "MS" << setw(10) << "\t\t" << "F" << setw(10) << "\t" << "P" << endl;
	cout << "Regression" << setw(10) << "\t" <<  df_model << setw(10) << "\t"  << SSR << setw(10) << "\t" << SSR << setw(10) << "\t" << computed_f << setw(10) << "\t" << p_value << endl;
	cout << "Residual Error" << setw(10) << "\t" << df_error << setw(10) << "\t" << SSE << setw(10) << "\t" << variance_estimate << endl;
	cout << "Total" << setw(10) << "\t\t" << df_total << setw(10) << "\t" << SST << endl;
	
	cout << "\nR-Squared"<< setw(10) << "\t" << "S" <<  endl;
	cout << 100*R_Squared << " %" << setw(10) << "\t" << s << endl;

	cout << "\nIf P value is less than the level of significance " << alpha << " then reject H0: β1 = 0"<< endl;
	cout << "\nRejection of H0 may suggest that the relationship is, indeed, linear."<< endl;
	cout << "The failure to reject H0: β1 = 0 suggests that there is no linear relationship between Y and x."<< endl;
	cout << "\n********************************************************" << endl;
	cout << "\nA "<< 100*(1-alpha) << "% confidence interval for the parameter β1 is : \n" << b1_lowerconfidenceinterval << " < β1 < " << b1_upperconfidenceinterval << endl; 	
	cout << "A "<< 100*(1-alpha) << "% confidence interval for the parameter β0 is : \n" << b0_lowerconfidenceinterval << " < β0 < " << b0_upperconfidenceinterval << endl; 	
	
	cout << setprecision(6) << "\nObs" << setw(10) << "y_data"  << "\t\t" << "Fit" << "\t\t"<< "SE Fit" << setw(21) << 100*(1-alpha) << "% CI" << setw(27) << 100*(1-alpha) << "% PI" << endl;
	for(int i=0; i < n; ++i)
	{
		fit = b1*x[i] +b0;
		SE_fit = s*sqrt(divisiond(1,n) + divisiond((x[i]- x_mean)*(x[i] - x_mean),Sxx));
		CI_lower = fit - abs(t_quantile)*sqrt(variance_estimate)*sqrt(divisiond(1,n) + divisiond((x[i]-x_mean)*(x[i]-x_mean),Sxx));
		CI_upper = fit + abs(t_quantile)*sqrt(variance_estimate)*sqrt(divisiond(1,n) + divisiond((x[i]-x_mean)*(x[i]-x_mean),Sxx));
		PI_lower = fit - abs(t_quantile)*sqrt(variance_estimate)*sqrt(1 + divisiond(1,n) + divisiond((x[i]-x_mean)*(x[i]-x_mean),Sxx));
		PI_upper = fit + abs(t_quantile)*sqrt(variance_estimate)*sqrt(1 + divisiond(1,n) + divisiond((x[i]-x_mean)*(x[i]-x_mean),Sxx)) ;
		cout << i+1 << "\t"<< y[i]  << "\t\t"<< fit << "\t\t" << SE_fit << setw(14) << "\t" << "("<< CI_lower <<" , " << CI_upper << ")" << "\t\t" << "("<< PI_lower <<" , " << PI_upper << ")" << endl;		
	}
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

vector<vector<double>> multipleregressionfull(vector<vector<double>> &X, vector<vector<double>> &y, double alpha)
{
	dmat X_t = transpose(X);
	dmat XtX= multiply(X_t,X);
	dmat inv = inverse(XtX);
	dmat A = multiply(inv,X_t);
	dmat b = multiply(A,y);
	vector<double> vec_b = getColumn(b,0);
	double SSE, SSR, SST, y_bar, MSR, MSE;
	double t_quantile, computed_f, p_value;
	int df_model, df_error, df_total;
	double RMSE, r_square, r_square_adj;
	double fit, SE_fit, CI_lower, CI_upper, PI_lower, PI_upper;

	cout << endl;
	int N = X.size();
	int n = b.size();
	for(int i = 0; i < n ; ++i)
	{
		cout <<"b["<< i <<"] = " << b[i][0] << endl;
	}
	
	cout << "\ny = " << b[0][0] ;
	for(int i = 1; i < n ; ++i)
	{
		cout << " +  " << b[i][0] << " x"<< i ;
	}
	
	for(int i = 0; i < N; ++i)
	{	
		y_bar += y[i][0];
	}	
	y_bar = divisiond(y_bar,N);

	for(int i = 0; i < N; ++i)
	{	
		double yi_hat = dot(getRow(X,i),vec_b);
		SSE += (yi_hat-y[i][0])*(yi_hat-y[i][0]);
		SST += (y[i][0] - y_bar)*(y[i][0] - y_bar);
	}
	SSR = SST - SSE;
	df_model = X[0].size() - 1 ;
	df_error = N - (df_model+1);
	df_total = N-1;
	MSR = divisiond(SSR,df_model);
	MSE = divisiond(SSE,df_error);
	computed_f = divisiond(MSR,MSE);
	p_value = 1 - Fcdf(computed_f,df_model,df_error);
	cout << "\n********************************************************" << endl;
	cout << "\nAnalysis of Variance" << endl;
	
	cout << "\nSource" << setw(10) << "\t\t" << "DF" << setw(10) << "\t" << "SS" << setw(10) << "\t\t" << "MS" << setw(10) << "\t\t" << "F" << setw(10) << "\t\t" << "P-value" << endl;
	cout << "Regression" << setw(10) << "\t" <<  df_model << setw(10) << "\t"  << SSR << setw(10) << "\t" << MSR << setw(10) << "\t" << computed_f << setw(10) << "\t" << p_value << endl;
	cout << "Residual Error" << setw(10) << "\t" << df_error << setw(10) << "\t" << SSE << setw(10) << "\t" << MSE << endl;
	cout << "Total" << setw(10) << "\t\t" << df_total << setw(10) << "\t" << SST << endl;
	
	cout << "\nRejection of H0 implies that the regression equation differs from a constant\n"<< endl;
	
	t_quantile = tquantile(1,df_error,1-(0.5*alpha));
	r_square = 1 - SSE/SST;
	r_square_adj = 1 - divisiond(SSE/df_error,SST/df_total);
	RMSE = sqrt(MSE);

	cout << "\nR-Square"<< setw(23) << "R-Square adjusted" << setw(23)  << "Grand Mean" << setw(23) << "Root MSE" << endl;
	cout << r_square << setw(23) << r_square_adj <<setw(23) << y_bar << setw(23) << RMSE << endl;

	// Compute the standard error from variable intercept to x1, x2, ..., xn
	// Compute the t value from variable intercept to x1, x2, ..., xn
	// Compute the P value from variable intercept to x1, x2, ..., xn
	dmat se = scalarmultiplication(inv,MSE);
	vector<double> v_stderror;
	vector<double> v_tvalue;
	vector<double> v_pvalue;
	for(int i = 0; i < n ; ++i)
	{
		v_stderror.push_back(sqrt(se[i][i]));
		v_tvalue.push_back(divisiond(vec_b[i],RMSE*sqrt(inv[i][i])));
		v_pvalue.push_back(1-tcdf(abs(v_tvalue[i]),df_error));
	}

	cout << "\nVariable" << setw(10) << "DF"  << "\t\t" << "Parameter Estimate" << "\t\t"<< "Standard Error" << setw(10) << "\t\t"<< " t value" << "\t\t" << "P-value" << endl;
	
	cout << "Intercept" << setprecision(6) << "\t"<< 1  << "\t\t"<< vec_b[0] << "\t\t\t" << v_stderror[0] << setw(14) << "\t\t" << v_tvalue[0] << "\t\t" << v_pvalue[0] << endl;			
	for(int i=0; i < n-1; ++i)
	{
		cout << "x" << i+1 << "\t\t"<< 1  << "\t\t"<< vec_b[i+1] << "\t\t\t" << v_stderror[i+1] << setw(14) << "\t\t" << v_tvalue[i+1]  << "\t\t" << v_pvalue[i+1] << endl;		
	}

	cout << "\nObs" << setw(11) << "y_data"  << "\t\t\t" << "Fit" << "\t\t\t"<< "SE Fit" << setw(11) << "\t"<< int(100*(1-alpha)) << "% CI" << setw(17) << "\t\t" << int(100*(1-alpha))<< "% PI" << setw(17)  << "\t\t" << "Residual" << endl;
	for(int i=0; i < N; ++i)
	{
		vector<vector<double>> x0 ;
		x0.push_back(getRow(X,i));
		dmat x0_t = transpose(x0);
		dmat r_mat = multiply(inv,x0_t);
		
		dvec v1 = getColumn(r_mat,0);
		dvec v2 = getColumn(x0_t,0);
		double std_error = dot(v1,v2);
		
		fit = dot(vec_b,v2);

		SE_fit = abs(t_quantile)*RMSE*sqrt( std_error );
		CI_lower = fit - abs(t_quantile)*RMSE*sqrt( std_error );
		CI_upper = fit + abs(t_quantile)*RMSE*sqrt( std_error );
		PI_lower = fit - abs(t_quantile)*RMSE*sqrt( 1 + std_error );
		PI_upper = fit + abs(t_quantile)*RMSE*sqrt( 1 + std_error );
		cout << i+1 << setprecision(6) << "\t"<< y[i][0]  << "\t\t"<< fit << "\t\t" << SE_fit << setw(14) << "\t" << "("<< CI_lower <<" , " << CI_upper << ")" << "\t\t" << "("<< PI_lower <<" , " << PI_upper << ")" << "\t\t" <<  y[i][0] - fit << endl;		
	}
	return b;
}

void ANOVA(vector<vector<double>> matrix)
{
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

void twoway_ANOVA(vector<vector<double>> matrix, int a, int b)
{
	vector<double> total_column;
	vector<double> mean_column;
	vector<double> total_row;
	vector<double> total_row_a;
	vector<double> mean_row;
	vector<double> sst_column;
	vector<double> ssa_column;

	double mean_total;
	double SST, SSA, SSB, SSAB, SSE, SSModel;

	int R = matrix.size();
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
	}
	for(int i = 0 ; i < R ; ++i)
	{
		double sum = 0;
		for (int j = 0 ; j < C ; ++j)
		{
			sum += matrix[i][j]; // sum of the row
		}
			total_row.push_back(sum);
	}
	int n = N/(a*b);
	int n_a = N/a; // number of samples per factor A
	for(int i = 0 ; i < R ;)
	{
		double sum = 0;
		for (int j = i ; j < i+n ; ++j)
		{
			sum += total_row[j]; // sum of the row for factor A level 1, 2, ..., a
		}
		total_row_a.push_back(sum);
		mean_row.push_back(sum/(n_a));		
		i = i+n;	
		
	}
	/*for(int i = 0 ; i < a ; ++i)
	{
		cout << "\nsum of row" << i << " = " << total_row_a[i] << endl;
		cout << "mean sum of row " << i << " = " << mean_row[i] << endl;
	}*/
	mean_total = std::accumulate(total_column.begin(), total_column.end(), 0.0) / N;
	
	// Compute SSA
	for(int i = 0 ; i < a ; ++i)
	{
		SSA += (mean_row[i] - mean_total)*(mean_row[i] - mean_total);
	}
	SSA = SSA*b*n;

	// Compute SSB	
	for(int i = 0 ; i < b ; ++i)
	{
		SSB += (mean_column[i] - mean_total)*(mean_column[i] - mean_total);
	}
	SSB = SSB*a*n;

	vector<vector<double>> mean_ij(a, vector<double>(C, 0.0));
	//cout << mean_ij.size() << endl;
	//cout << mean_ij[0].size() << endl;
	
	// To compute mean of the observations in the (ij)th cell
	int index = 0;
	for(int k = 0 ; k < a ; ++k)
	{
		for(int i = 0 ; i < C ; ++i)
		{
			double sum = 0;	
			for (int j = index ; j < index+n ; ++j)
			{
				sum += matrix[j][i]; 
			}
			mean_ij[k][i] = sum/n;
			//cout << "sum of  " << k << ","<< i << "cell = "<< mean_ij[k][i]<< endl;
		}
		index = index + n;
	}

	// Compute the marginal averages
	vector<double> column_marginalaverages;
	vector<double> row_marginalaverages;
	for(int i = 0 ; i < a ; ++i)
	{
		vector<double> row_marginal = getRow(mean_ij,i);
		row_marginalaverages.push_back(std::accumulate(row_marginal.begin(), row_marginal.end(), 0.0) / b);		
	}
	for(int i = 0 ; i < b ; ++i)
	{
		vector<double> column_marginal = getColumn(mean_ij,i);
		column_marginalaverages.push_back(std::accumulate(column_marginal.begin(), column_marginal.end(), 0.0) / a);		
	}

	// Show the two-way table of averages
	cout << "\n********************************************************" << endl;
	cout << "\nTwo-way table of averages" << endl;
	for (int i = 0 ; i < a ; ++i)
	{
		for (int j = 0 ; j < b ; ++j)
		{
			cout << std::fixed << setw(23) << setprecision(6) << mean_ij[i][j] ;
		}
		cout << setw(23) << " | " <<row_marginalaverages[i]<< endl;
	}
	for (int i = 0 ; i < b ; ++i)
	{
		cout << "--------------------------------" ;
	}
	cout << endl;
	for (int i = 0 ; i < b ; ++i)
	{
		cout << std::fixed << setw(23) << setprecision(6) <<  column_marginalaverages[i] ;
	}
	
	// Compute SSE
	index = 0;
	for(int k = 0 ; k < a ; ++k)
	{
		for(int i = 0 ; i < C ; ++i)
		{	
			for (int j = index ; j < index+n ; ++j)
			{
				SSE += (matrix[j][i] - mean_ij[k][i])*(matrix[j][i] - mean_ij[k][i]); 
			}
		}
		index = index + n;
	}
	
	// Compute SST
	for(int i = 0 ; i < R ; ++i)
	{	
		for (int j = 0 ; j < C ; ++j)
		{
			SST += (matrix[i][j] - mean_total)*(matrix[i][j] - mean_total); 
		}
	}

	SSAB = SST - SSA - SSB - SSE;
	SSModel = SST - SSE;
	int df_total = N-1;
	int df_error = a*b*(n-1);
	int df_model= df_total - df_error;
	double MS_Model = divisiond(SSModel,df_model);
	double MSE = divisiond(SSE,df_error);
	double computed_f = divisiond(MS_Model,MSE);
	double p_value = 1 - Fcdf(abs(computed_f),df_model,df_error);
	double r_square = 1 - SSE/SST;

	cout << "\n********************************************************" << endl;
	cout << "\nHypotheses for the Two-Factor Problem" << endl;
	cout << "\nH0': α_{1} = α_{2} = ... = α_{a} = 0" << endl;
	cout << "H1':  At least one of the α_{i} is not equal to zero" << endl;
	cout << "\nH0'': β_{1} = β_{2} = ... = β_{b} = 0" << endl;
	cout << "H1'':  At least one of the β_{j} is not equal to zero" << endl;
	cout << "\nH0''': (αβ)_{11} = (αβ)_{12} = ... = (αβ)_{ab} = 0" << endl;
	cout << "H1''':  At least one of the (αβ)_{ij} is not equal to zero" << endl;

	cout << "\n********************************************************" << endl;
	cout << "\nAnalysis of Variance" << endl;
	
	cout << "\nSource" << setw(10) << "\t" << "DF" << setw(10) << "\t" << "SS" << setw(10) << "\t" << "MS" << "\t\t" << "F" << setw(10) << "\t" << "P" << endl;
	cout << "Model" << setw(10) << "\t" <<  df_model << setw(10) << "\t" << SSModel << setw(5) << "\t" << MS_Model << setw(5) << "\t" << computed_f << setw(5) << "\t" << p_value << endl;
	cout << "Error" << setw(10) << "\t" << df_error << setw(10) << "\t" << SSE << setw(5) << "\t" << MSE << endl;
	cout << "Total" << setw(10) << "\t" << df_total << setw(10) << "\t" << SST << endl;
	
	cout << "\nR-Square"<< setw(23) << "Grand Mean" << setw(23) << "Root MSE" << endl;
	cout << r_square << setw(23) << mean_total << setw(23) << sqrt(MSE) << endl;

	int df_A = a-1;
	int df_B = b-1;
	int df_AB = (a-1)*(b-1);
	double MS_A = divisiond(SSA,df_A);
	double MS_B = divisiond(SSB,df_B);
	double MS_AB = divisiond(SSAB,df_AB);
	double computed_f_A = divisiond(MS_A,MSE);
	double computed_f_B = divisiond(MS_B,MSE);
	double computed_f_AB = divisiond(MS_AB,MSE);
	double p_value_A = 1 - Fcdf(computed_f_A,df_A,df_error);
	double p_value_B = 1 - Fcdf(computed_f_B,df_B,df_error);
	double p_value_AB = 1 - Fcdf(computed_f_AB,df_AB,df_error);
	cout << "\nAnalysis of Variance for the Two-Factor Experiment with n replications" << endl;
	
	cout << "\nSource" << setw(10) << "\t" << "DF" << setw(10) << "\t" << "Type III SS" << setw(10) << "\t" << "MS" << setw(10) << "\t\t" << "F" << setw(10) << "\t\t" << "P" << endl;
	cout << "A" << setprecision(6) << setw(10) << "\t" <<  df_A << setw(10) << "\t" << SSA << setw(10) << "\t" << MS_A << setw(10) << "\t" << computed_f_A << setw(10) << "\t" << p_value_A << endl;
	cout << "B" << setw(10) << "\t" << df_B << setw(10) << "\t" << SSB << setw(10) << "\t" <<  MS_B << setw(10) << "\t" << computed_f_B << setw(10) << "\t" << p_value_B << endl;
	cout << "A * B" << setw(10) << "\t" << df_AB << setw(10) << "\t" << SSAB << setw(10) << "\t" << MS_AB << setw(10) << "\t" << computed_f_AB << setw(5) << "\t\t" << p_value_AB << endl;
	
}

void threeway_ANOVA(vector<vector<double>> matrix, int a, int b, int c)
{
	vector<double> total_column;
	vector<double> mean_column;
	vector<double> total_row;
	vector<double> total_factor_A;
	vector<double> mean_factor_A;
	vector<double> mean_factor_B;
	vector<double> total_factor_C;
	vector<double> mean_factor_C;
	vector<double> sst_column;
	vector<double> ssa_column;

	double mean_total;
	double SST, SSA, SSB, SSC, SSAB, SSAC, SSBC, SSABC, SSE, SSModel;

	int R = matrix.size();
	int C = matrix[0].size();
	int N = 0; // computing total data

	for(int i = 0 ; i < C ; ++i)  // Compute the average for the j-th level of factor B 
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
	}
	for(int i = 0 ; i < R ; ++i)
	{
		double sum = 0;
		for (int j = 0 ; j < C ; ++j)
		{
			sum += matrix[i][j]; // sum of the row
		}
			total_row.push_back(sum);
	}

	int n = N/(a*b*c);
	int n_a = N/a; // number of samples per factor A
	for(int i = 0 ; i < R ;) // Compute the average for the i-th level of factor A 
	{
		double sum = 0;
		for (int j = i ; j < i+n ; ++j)
		{
			sum += total_row[j]; // sum of the row for factor A level 1, 2, ..., a
		}
		total_factor_A.push_back(sum);
		mean_factor_A.push_back(sum/(n_a));		
		i = i+n;	
		
	}
	int n_c = N/c;
	for(int i = 0 ; i < C ;)// Compute the average for the k-th level of factor C
	{
		double sum = 0;
		for (int j = i ; j < i+b ; ++j)
		{
			sum += total_column[j]; // sum of the columns for factor C level 1, 2, ..., c
		}
		total_factor_C.push_back(sum);
		mean_factor_C.push_back(sum/(n_c));		
		i = i+b;	
		
	}

	/*
	for(int i = 0 ; i < a ; ++i) // Compute the average for the i-th level of factor A 
	{
		cout << "\nsum of row" << i << " = " << total_factor_a[i] << endl;
		cout << "mean sum of row " << i << " = " << mean_factor_A[i] << endl;
	}
	for(int i = 0 ; i < C ; ++i) // Compute the average for the j-th level of factor B 
	{
		cout << "\nsum of column" << i << " = " << total_column[i] << endl;
		cout << "mean sum of column " << i << " = " << mean_column[i] << endl;
	}
	for(int i = 0 ; i < c ; ++i) // Compute the average for the k-th level of factor C
	{
		cout << "\nsum of factor C" << i << " = " << total_factor_C[i] << endl;
		cout << "mean sum of factor C " << i << " = " << mean_factor_C[i] << endl;
	}
	*/


	mean_total = std::accumulate(total_column.begin(), total_column.end(), 0.0) / N;
	
	// Compute SSA
	for(int i = 0 ; i < a ; ++i)
	{
		SSA += (mean_factor_A[i] - mean_total)*(mean_factor_A[i] - mean_total);
	}
	
	SSA = SSA*b*c*n;
	
	// Compute SSB	
	for(int i = 0 ; i < b ; ++i)
	{
		double sum = 0;
		for (int j = 0; j < c; ++j)
		{
			sum += mean_column[i+b*j] ; // sum of the columns for factor B level 1, 2, ..., b 
		}
		mean_factor_B.push_back(sum/(c));		
	}
	for(int i = 0 ; i < b ; ++i)
	{
		SSB += (mean_factor_B[i] - mean_total)*(mean_factor_B[i] - mean_total);
		
	}
	SSB = SSB*a*c*n;
	
	// Compute SSC
	for(int i = 0 ; i < c ; ++i)
	{
		SSC += (mean_factor_C[i] - mean_total)*(mean_factor_C[i] - mean_total);
	}
	SSC = SSC*a*b*n;
	
	// Compute SST
	for(int i = 0 ; i < R ; ++i)
	{	
		for (int j = 0 ; j < C ; ++j)
		{
			SST += (matrix[i][j] - mean_total)*(matrix[i][j] - mean_total); 
		}
	}
	vector<vector<double>> mean_ij(a, vector<double>(C, 0.0));
	
	// To compute mean of the observations in the (ij)th cell
	int index = 0;
	for(int k = 0 ; k < a ; ++k)
	{
		for(int i = 0 ; i < C ; ++i)
		{
			double sum = 0;	
			for (int j = index ; j < index+n ; ++j)
			{
				sum += matrix[j][i]; 
			}
			mean_ij[k][i] = sum/n;
		}
		index = index + n;
	}
	

	vector<vector<double>> mean_ij_final(a, vector<double>(b, 0.0));
	
	for(int k = 0 ; k < a ; ++k)
	{
		for(int i = 0 ; i < b ; ++i)
		{
			double sum = 0;	
			for (int j = 0 ; j < c ; ++j)
			{
				sum += mean_ij[k][i+b*j]; 
			}
			mean_ij_final[k][i] = sum/c;
		}
	}
	//printMatrix(mean_ij_final);

	vector<vector<double>> mean_jk(b, vector<double>(c, 0.0));
	for(int i = 0 ; i < b ; ++i)
	{
		for(int j = 0 ; j < c ; ++j)
		{
			vector<double> column_marginal = getColumn(matrix,i+j*b);		
			mean_jk[i][j] = std::accumulate(column_marginal.begin(), column_marginal.end(), 0.0) / R ; 
		}
	}
	//printMatrix(mean_jk);
	
	vector<vector<double>> mean_ik(a, vector<double>(c, 0.0));

	int n_ik = N/(a*c);
	int row_A = R/a;
	int col_C = C/c;
	for(int k_row = 0 ; k_row < a ; ++k_row)
	{
		for(int k_col = 0 ; k_col < c ; ++k_col)
		{
			double sum = 0;	
			for(int i = 0 ; i < row_A ; ++i)
			{
				for (int j = 0 ; j < col_C ; ++j)
				{
					sum += matrix[i+ k_row*row_A][j + k_col*col_C]; 
				}				
			}			
			mean_ik[k_row][k_col] = sum/n_ik;
		}		
	}
	//printMatrix(mean_ik);

	vector<vector<double>> mean_ijk(a, vector<double>(b*c, 0.0));

	int col_BC = b*c;
	for(int k_row = 0 ; k_row < a ; ++k_row)
	{
		for(int k_col = 0 ; k_col < col_BC ; ++k_col)
		{
			double sum = 0;	
			for(int i = 0 ; i < n ; ++i)
			{
				sum += matrix[i+ k_row*row_A][k_col]; 				
			}			
			mean_ijk[k_row][k_col] = sum/n;
		}		
	}
	//printMatrix(mean_ijk);

	// Compute SS(AB)
	for(int i = 0 ; i < a ; ++i)
	{	
		for (int j = 0 ; j < b ; ++j)
		{
			SSAB += (mean_ij_final[i][j] - mean_factor_A[i] - mean_factor_B[j] + mean_total)*(mean_ij_final[i][j] - mean_factor_A[i] - mean_factor_B[j] + mean_total); 
		}
	}
	SSAB = c*n*SSAB;
	
	// Compute SS(AC)
	for(int i = 0 ; i < a ; ++i)
	{	
		for (int j = 0 ; j < c ; ++j)
		{
			SSAC += (mean_ik[i][j] - mean_factor_A[i] - mean_factor_C[j] + mean_total)*(mean_ik[i][j] - mean_factor_A[i] - mean_factor_C[j] + mean_total); 
		}
	}
	SSAC = b*n*SSAC;

	// Compute SS(BC)
	for(int i = 0 ; i < b ; ++i)
	{	
		for (int j = 0 ; j < c ; ++j)
		{
			SSBC += (mean_jk[i][j] - mean_factor_B[i] - mean_factor_C[j] + mean_total)*(mean_jk[i][j] - mean_factor_B[i] - mean_factor_C[j] + mean_total); 
		}
	}
	SSBC = a*n*SSBC;

	// Compute SSE
	index = 0;
	for(int k = 0 ; k < a ; ++k)
	{
		for(int i = 0 ; i < col_BC ; ++i)
		{	
			for (int j = index ; j < index+n ; ++j)
			{
				SSE += (matrix[j][i] - mean_ijk[k][i])*(matrix[j][i] - mean_ijk[k][i]); 
			}
		}
		index = index + n;
	}

	// Compute SS(ABC)
	SSABC = SST - SSA - SSB - SSC - SSAB - SSAC- SSBC - SSE;
	
	
	// Compute the marginal averages
	vector<double> column_marginalaverages;
	vector<double> row_marginalaverages;
	for(int i = 0 ; i < a ; ++i)
	{
		vector<double> row_marginal = getRow(mean_ij_final,i);
		row_marginalaverages.push_back(std::accumulate(row_marginal.begin(), row_marginal.end(), 0.0) / b);		
	}
	for(int i = 0 ; i < b ; ++i)
	{
		vector<double> column_marginal = getColumn(mean_ij_final,i);
		column_marginalaverages.push_back(std::accumulate(column_marginal.begin(), column_marginal.end(), 0.0) / a);		
	}

	
	//	************************AB table of averages************************
	
	vector<double> column_marginalaverages_ij;
	vector<double> row_marginalaverages_ij;
	for(int i = 0 ; i < a ; ++i)
	{
		vector<double> row_marginal_ij = getRow(mean_ij_final,i);
		row_marginalaverages_ij.push_back(std::accumulate(row_marginal_ij.begin(), row_marginal_ij.end(), 0.0) / b);		
	}
	for(int i = 0 ; i < b ; ++i)
	{
		vector<double> column_marginal_ij = getColumn(mean_ij_final,i);
		column_marginalaverages_ij.push_back(std::accumulate(column_marginal_ij.begin(), column_marginal_ij.end(), 0.0) / a);		
	}
	// Show the two-way table of averages
	cout << "\n********************************************************" << endl;
	cout << "\nTwo-way table of averages, A (row) and B (column)" << endl;
	for (int i = 0 ; i < a ; ++i)
	{
		for (int j = 0 ; j < b ; ++j)
		{
			cout << std::fixed << setw(23) << setprecision(6) << mean_ij_final[i][j] ;
		}
		cout << setw(23) << " | " <<row_marginalaverages_ij[i]<< endl;
	}
	for (int i = 0 ; i < a ; ++i)
	{
		cout << "--------------------------------" ;
	}
	cout << endl;
	for (int i = 0 ; i < b ; ++i)
	{
		cout << std::fixed << setw(23) << setprecision(6) <<  column_marginalaverages_ij[i] ;
	}
	//	***********************************************************************

	//	************************AC table of averages************************
	vector<double> column_marginalaverages_ik;
	vector<double> row_marginalaverages_ik;
	for(int i = 0 ; i < a ; ++i)
	{
		vector<double> row_marginal_ik = getRow(mean_ik,i);
		row_marginalaverages_ik.push_back(std::accumulate(row_marginal_ik.begin(), row_marginal_ik.end(), 0.0) / c);		
	}
	for(int i = 0 ; i < c ; ++i)
	{
		vector<double> column_marginal_ik = getColumn(mean_ik,i);
		column_marginalaverages_ik.push_back(std::accumulate(column_marginal_ik.begin(), column_marginal_ik.end(), 0.0) / a);		
	}
	// Show the two-way table of averages
	cout << "\n\n********************************************************" << endl;
	cout << "\nTwo-way table of averages, A (row) and C (column)" << endl;
	for (int i = 0 ; i < a ; ++i)
	{
		for (int j = 0 ; j < c ; ++j)
		{
			cout << std::fixed << setw(23) << setprecision(6) << mean_ik[i][j] ;
		}
		cout << setw(23) << " | " <<row_marginalaverages_ik[i]<< endl;
	}
	for (int i = 0 ; i < a ; ++i)
	{
		cout << "--------------------------------" ;
	}
	cout << endl;
	for (int i = 0 ; i < c ; ++i)
	{
		cout << std::fixed << setw(23) << setprecision(6) <<  column_marginalaverages_ik[i] ;
	}
	//	***********************************************************************

	//	************************BC table of averages************************
	vector<double> column_marginalaverages_jk;
	vector<double> row_marginalaverages_jk;
	for(int i = 0 ; i < b ; ++i)
	{
		vector<double> row_marginal_jk = getRow(mean_jk,i);
		row_marginalaverages_jk.push_back(std::accumulate(row_marginal_jk.begin(), row_marginal_jk.end(), 0.0) / c);		
	}
	for(int i = 0 ; i < c ; ++i)
	{
		vector<double> column_marginal_jk = getColumn(mean_jk,i);
		column_marginalaverages_jk.push_back(std::accumulate(column_marginal_jk.begin(), column_marginal_jk.end(), 0.0) / b);		
	}
	// Show the two-way table of averages
	cout << "\n\n********************************************************" << endl;
	cout << "\nTwo-way table of averages, B (row) and C (column)" << endl;
	for (int i = 0 ; i < b ; ++i)
	{
		for (int j = 0 ; j < c ; ++j)
		{
			cout << std::fixed << setw(23) << setprecision(6) << mean_jk[i][j] ;
		}
		cout << setw(23) << " | " <<row_marginalaverages_jk[i]<< endl;
	}
	for (int i = 0 ; i < b ; ++i)
	{
		cout << "--------------------------------" ;
	}
	cout << endl;
	for (int i = 0 ; i < c ; ++i)
	{
		cout << std::fixed << setw(23) << setprecision(6) <<  column_marginalaverages_jk[i] ;
	}
	//	***********************************************************************
	
	SSModel = SST - SSE;
	int df_total = N-1;
	int df_error = a*b*c*(n-1);
	int df_model= df_total - df_error;
	double MS_Model = divisiond(SSModel,df_model);
	double MSE = divisiond(SSE,df_error);
	double computed_f = divisiond(MS_Model,MSE);
	double p_value = 1 - Fcdf(abs(computed_f),df_model,df_error);
	double r_square = 1 - SSE/SST;

	cout << "\n********************************************************" << endl;
	cout << "\nHypotheses for the Three-Factor Problem" << endl;
	cout << "\nH0^{(1)}: α_{1} = α_{2} = ... = α_{a} = 0" << endl;
	cout << "H1^{(1)}:  At least one of the α_{i} is not equal to zero" << endl;
	cout << "\nH0^{(2)}: β_{1} = β_{2} = ... = β_{b} = 0" << endl;
	cout << "H1^{(2)}:  At least one of the β_{j} is not equal to zero" << endl;
	cout << "\nH0^{(3)}: γ_{1} = γ_{2} = ... = γ_{c} = 0" << endl;
	cout << "H1^{(3)}:  At least one of the γ_{k} is not equal to zero" << endl;
	cout << "\nH0^{(4)}: (αβ)_{11} = (αβ)_{12} = ... = (αβ)_{ab} = 0" << endl;
	cout << "H1^{(4)}:  At least one of the (αβ)_{ij} is not equal to zero" << endl;
	cout << "\nH0^{(5)}: (αγ)_{11} = (αγ)_{12} = ... = (αγ)_{ac} = 0" << endl;
	cout << "H1^{(5)}:  At least one of the (αγ)_{ik} is not equal to zero" << endl;
	cout << "\nH0^{(6)}: (βγ)_{11} = (βγ)_{12} = ... = (βγ)_{bc} = 0" << endl;
	cout << "H1^{(6)}:  At least one of the (βγ)_{jk} is not equal to zero" << endl;
	cout << "\nH0^{(7)}: (αβγ)_{111} = (αβγ)_{112} = ... = (αβγ)_{abc} = 0" << endl;
	cout << "H1^{(7)}:  At least one of the (αβγ)_{ijk} is not equal to zero" << endl;

	cout << "\n********************************************************" << endl;
	cout << "\nAnalysis of Variance" << endl;
	
	cout << "\nSource" << setw(10) << "\t" << "DF" << setw(10) << "\t" << "SS" << setw(10) << "\t" << "MS" << "\t\t" << "F" << setw(10) << "\t" << "P" << endl;
	cout << "Model" << setw(10) << "\t" <<  df_model << setw(10) << "\t" << SSModel << setw(5) << "\t" << MS_Model << setw(5) << "\t" << computed_f << setw(5) << "\t" << p_value << endl;
	cout << "Error" << setw(10) << "\t" << df_error << setw(10) << "\t" << SSE << setw(5) << "\t" << MSE << endl;
	cout << "Total" << setw(10) << "\t" << df_total << setw(10) << "\t" << SST << endl;
	
	cout << "\nR-Square"<< setw(23) << "Grand Mean" << setw(23) << "Root MSE" << endl;
	cout << r_square << setw(23) << mean_total << setw(23) << sqrt(MSE) << endl;

	int df_A = a-1;
	int df_B = b-1;
	int df_C = c-1;
	int df_AB = (a-1)*(b-1);
	int df_AC = (a-1)*(c-1);
	int df_BC = (b-1)*(c-1);
	int df_ABC = (a-1)*(b-1)*(c-1);
	double MS_A = divisiond(SSA,df_A);
	double MS_B = divisiond(SSB,df_B);
	double MS_C = divisiond(SSC,df_C);
	double MS_AB = divisiond(SSAB,df_AB);
	double MS_AC = divisiond(SSAC,df_AC);
	double MS_BC = divisiond(SSBC,df_BC);
	double MS_ABC = divisiond(SSABC,df_ABC);
	double computed_f_A = divisiond(MS_A,MSE);
	double computed_f_B = divisiond(MS_B,MSE);
	double computed_f_C = divisiond(MS_C,MSE);
	double computed_f_AB = divisiond(MS_AB,MSE);
	double computed_f_AC = divisiond(MS_AC,MSE);
	double computed_f_BC = divisiond(MS_BC,MSE);
	double computed_f_ABC = divisiond(MS_ABC,MSE);
	double p_value_A = 1 - Fcdf(computed_f_A,df_A,df_error);
	double p_value_B = 1 - Fcdf(computed_f_B,df_B,df_error);
	double p_value_C = 1 - Fcdf(computed_f_C,df_C,df_error);
	double p_value_AB = 1 - Fcdf(computed_f_AB,df_AB,df_error);
	double p_value_AC = 1 - Fcdf(computed_f_AC,df_AC,df_error);
	double p_value_BC = 1 - Fcdf(computed_f_BC,df_BC,df_error);
	double p_value_ABC = 1 - Fcdf(computed_f_ABC,df_ABC,df_error);
	cout << "\nAnalysis of Variance for the Three-Factor Experiment with n replications" << endl;
	
	cout << "\nSource" << setw(10) << "\t" << "DF" << setw(10) << "\t" << "Type III SS" << setw(10) << "\t" << "MS" << setw(10) << "\t\t" << "F" << setw(10) << "\t\t" << "P" << endl;
	cout << "A" << setprecision(6) << setw(10) << "\t" <<  df_A << setw(10) << "\t" << SSA << setw(10) << "\t" << MS_A << setw(10) << "\t" << computed_f_A << setw(10) << "\t" << p_value_A << endl;
	cout << "B" << setw(10) << "\t" << df_B << setw(10) << "\t" << SSB << setw(10) << "\t" <<  MS_B << setw(10) << "\t" << computed_f_B << setw(10) << "\t" << p_value_B << endl;
	cout << "C" << setw(10) << "\t" << df_C << setw(10) << "\t" << SSC << setw(10) << "\t" <<  MS_C << setw(10) << "\t" << computed_f_C << setw(10) << "\t" << p_value_C << endl;
	cout << "A * B" << setw(10) << "\t" << df_AB << setw(10) << "\t" << SSAB << setw(10) << "\t" << MS_AB << setw(10) << "\t" << computed_f_AB << setw(5) << "\t\t" << p_value_AB << endl;
	cout << "A * C" << setw(10) << "\t" << df_AC << setw(10) << "\t" << SSAC << setw(10) << "\t" << MS_AC << setw(10) << "\t" << computed_f_AC << setw(5) << "\t\t" << p_value_AC << endl;
	cout << "B * C" << setw(10) << "\t" << df_BC << setw(10) << "\t" << SSBC << setw(10) << "\t" << MS_BC << setw(10) << "\t" << computed_f_BC << setw(5) << "\t\t" << p_value_BC << endl;
	cout << "A * B * C" << "\t" << df_ABC << setw(10) << "\t" << SSABC << setw(10) << "\t" << MS_ABC << setw(10) << "\t" << computed_f_ABC << setw(5) << "\t\t" << p_value_ABC << endl;
	cout << "Error" << setw(10) << "\t" << df_error << setw(10) << "\t" << SSE << setw(10) << "\t" << MSE << endl;
	cout << "-----------------------------------------------------------" << endl;
	cout << "Total" << setw(10) << "\t" << df_total << setw(10) << "\t" << SST << endl;
	
}
#endif
#endif