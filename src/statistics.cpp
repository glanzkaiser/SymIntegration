/*
   
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)
#include <boost/math/distributions/chi_squared.hpp> // Faster computing of chi squared cdf
#include <boost/math/distributions/beta.hpp> //Faster computing of beta cdf and pdf
#include <boost/math/distributions/fisher_f.hpp> //Faster computing of Fisher' F cdf and pdf
#include <boost/math/distributions/gamma.hpp> // Faster computing of gamma cdf
#include <boost/math/distributions/laplace.hpp> //Faster computing of laplace cdf and pdf
#include <boost/math/distributions/logistic.hpp> //Faster computing of logistic cdf and pdf
#include <boost/math/distributions/students_t.hpp> //Faster computing of students' t cdf and pdf

#include <random> // For random number generation
#include <chrono>

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
 	Symbolic mean = n*p;

	return mean;
}

Symbolic binomialvar(int x, int n, double p)
{
 	Symbolic var = n*p*(1-p);

	return var;
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
 	Symbolic mean = k/p;

	return mean;
}

Symbolic negativebinomialvar(int x, int k, double p)
{
 	Symbolic var = k*(1-p)/(p*p);

	return var;
}

Symbolic negativebinomialmgf(const Symbolic &x, const Symbolic &k, const Symbolic &p)
{
 	Symbolic t("t");

	Symbolic mgf =pow(p,k)*pow(1-((1-p)*exp(t)),-k);
	
	return mgf;
}

Symbolic geometricpmf(int x, double p)
{
 	Symbolic pmf = p*(pow(1-p,x-1));
	
	return pmf;
}

Symbolic geometricmean(int x, double p)
{
 	Symbolic mean = 1/p;

	return mean;
}

Symbolic geometricvar(int x, double p)
{
 	Symbolic var = (1-p)/(p*p);

	return var;
}

Symbolic geometricmgf(const Symbolic &x, const Symbolic &p)
{
 	Symbolic t("t");

	Symbolic mgf =p*pow(1-((1-p)*exp(t)),Symbolic(-1));
	
	return mgf;
}

Symbolic hypergeometricpmf(int x, int N, int n, int k)
{
 	Symbolic pmf = combinations1(k,x)*divisions(combinations1(N-k,n-x),combinations1(N,n));
	
	return pmf;
}

Symbolic hypergeometricmean(int x, int N, int n, int k)
{
 	Symbolic mean = divisionint(n*k,N);

	return mean;
}

Symbolic hypergeometricvar(int x, int N, int n, int k)
{
 	Symbolic var = n*divisionint(N-n,N-1)*divisionint(k,N)*(1-divisionint(k,N));

	return var;
}

Symbolic poissonpmf(int x, int k)
{
 	Symbolic pmf = exp(-k)*pow(k,x)/(factorial1(x));
	
	return pmf;
}

Symbolic poissoncdf(int x, int k)
{
 	Symbolic cdf;
	for (int i=0; i<=x;i++)
	{
		cdf += exp(-k)*pow(k,i)/(factorial1(i));
	}	
	return cdf;
}

Symbolic poissonmean(int x, int k)
{
 	Symbolic mean = k;

	return mean;
}

Symbolic poissonvar(int x, int k)
{
 	Symbolic var = k;

	return var;
}

Symbolic poissonmgf(const Symbolic &x, const Symbolic &k)
{
 	Symbolic t("t");

	Symbolic mgf =exp(k*(exp(t)-1));
	
	return mgf;
}

// Continuous Distributions

Symbolic uniformpdf(double x, double a, double b)
{
 	Symbolic pdf = divisions(1,b-a);

	return pdf;
}

Symbolic uniformcdf(double x1, double a, double b)
{
	return divisions(x1,b-a) - divisions(a,b-a) ;
}

Symbolic uniformmgf(double x, double a, double b)
{
 	Symbolic t("t");
	Symbolic mgf = (exp(b*t) - exp(a*t))/((b-a)*t);

	return mgf;
}

Symbolic uniformmean(double x, double a, double b)
{
 	Symbolic mean = 0.5*(a+b);

	return mean;
}

Symbolic uniformvar(double x, double a, double b)
{
 	Symbolic var = divisions(1,12)*(b-a)*(b-a);

	return var;
}

Symbolic normalpdf(double x, double μ, double σ)
{
 	Symbolic pdf = divisions(1,sqrt(2*π)*σ)*exp(-0.5*divisions(x-μ,σ)*divisions(x-μ,σ));

	return evalf(pdf,1,1);
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
	
 	Symbolic cdf = 0.5 * erfc(-x * M_SQRT1_2);

	return cdf;
}

Symbolic normalmgf(double x, double μ, double σ)
{
 	Symbolic t("t");
	Symbolic mgf = exp(μ*t + 0.5*σ*σ*t*t);

	return mgf;
}

Symbolic normalmean(double x, double μ, double σ)
{
 	Symbolic mean = μ;

	return mean;
}

Symbolic normalvar(double x, double μ, double σ)
{
 	Symbolic var = σ*σ;

	return var;
}

Symbolic gammapdf(double x, double α, double β)
{
 	//Symbolic b = divisions(1,(Symbolic(β)^(Symbolic(α)))*tgamma(α))*(Symbolic(x)^(Symbolic(α-1)))*exp(-divisions(x,β));
	//return evalf(b,1,1);

	boost::math::gamma_distribution<> gamma_dist(α,β);
	double x_value = x;
	double pdf_value = boost::math::pdf(gamma_dist, x_value);

	return pdf_value;
}

Symbolic gammacdf(double x, double α, double β)
{	
	boost::math::gamma_distribution<> gamma_dist(α,β);
	double x_value = x;
	double cdf_value = boost::math::cdf(gamma_dist, x_value);

	return cdf_value;
}

Symbolic gammamgf(double x, double α, double β)
{
 	Symbolic t("t");
	Symbolic mgf = ((1-β*t)^(Symbolic(-α)));

	return mgf;
}

Symbolic gammamean(double x, double α, double β)
{
 	Symbolic mean = α*β;

	return mean;
}

Symbolic gammavar(double x, double α, double β)
{
 	Symbolic var = α*β*β;

	return var;
}

Symbolic exponentialpdf(double x, double λ)
{
 	Symbolic pdf = λ*exp(-λ*x);

	return evalf(pdf,1,1);
}

Symbolic exponentialcdf(double x1, double λ)
{	
	Symbolic x("x");
	Symbolic f =  λ*exp(-λ*x);

 	Symbolic cdf = 1 + integrate(f,x)[x==x1] ;

	return evalf(cdf,1,1);
}

Symbolic exponentialmgf(double x, double λ)
{
 	Symbolic t("t");
	Symbolic mgf = (1-(divisions(t,λ)))^(Symbolic(-1));

	return mgf;
}

Symbolic exponentialmean(double x, double λ)
{
 	Symbolic mean = divisions(1,λ);

	return mean;
}

Symbolic exponentialvar(double x, double λ)
{
 	Symbolic var = divisions(1,λ*λ);

	return var;
}

Symbolic betapdf(double x, double α, double β)
{
 	//Symbolic b = ( tgamma(α+β)/(tgamma(α)*tgamma(β)) ) * (x^(Symbolic(α-1))) * ((1-x)^(Symbolic(β-1)));
	//return evalf(b,1,1);

	boost::math::beta_distribution<> my_beta(α,β);
	double pdf_value = boost::math::pdf(my_beta,x);

	return pdf_value;
}

Symbolic betacdf(double x, double α, double β)
{	
	/*Symbolic x("x");
	int alpha  = α;
	int beta = β;
	Symbolic f = ( tgamma(α+β)/(tgamma(α)*tgamma(β)) ) * (pow(x,Symbolic(alpha-1))) * (pow(1-x,(Symbolic(beta-1))));

 	Symbolic b = integrate(f,x)[x==x1] ;

	return evalf(b,1,1);*/

	boost::math::beta_distribution<> my_beta(α,β);
	double cdf_value = boost::math::cdf(my_beta,x);

	return cdf_value;
}

Symbolic betamgf(double x, double α, double β)
{
 	Symbolic t("t");
	Symbolic mgf = hypergeometric_1F1(α, α+β,t,5);

	return mgf;
}

Symbolic betamean(double x, double α, double β)
{
 	Symbolic mean = divisions(α,α+β);

	return mean;
}

Symbolic betavar(double x, double α, double β)
{
 	Symbolic var = divisions(α*β,(α+β+1)*(α+β)*(α+β));

	return var;
}

Symbolic cauchypdf(double x)
{
 	Symbolic pdf = divisionint(1,π)*divisionint(1,x*x+1);

	return pdf;
}


Symbolic chisquaredpdf(double x, double r)
{
	int a = 0.5*r;
	int c = 0.5*r-1;
 	Symbolic pdf = divisionint(1,tgamma(0.5*r)*pow(2,a)) * pow(x,c) * exp(-0.5*x);

	return evalf(pdf,1,1);
}


Symbolic chisquaredcdf(double x, double r)
{
	// We are using Boost to compute chi squared cdf
	// We try to use the lower and upper incomplete gamma function but it cannot produce the correct result
	double degrees_of_freedom = r; 
	boost::math::chi_squared_distribution<> chi_squared(degrees_of_freedom);

	double x_value = x; // The value at which to evaluate the CDF
	double cdf_value = boost::math::cdf(chi_squared, x_value);

	return cdf_value;
}

Symbolic chisquaredmgf(double x, double r)
{
 	Symbolic t("t");
	Symbolic mgf= (1-2*t)^(Symbolic(-0.5*r));

	return mgf;
}

Symbolic chisquaredmean(double x, double r)
{
 	Symbolic mean = r;

	return mean;
}

Symbolic chisquaredvar(double x, double r)
{
 	Symbolic var = 2*r;

	return var;
}

Symbolic Fpdf(double x, double r1, double r2)
{
	boost::math::fisher_f_distribution<> fisher_f(r1,r2);
	double pdf_value = boost::math::pdf(fisher_f,x);

	return pdf_value;
}

Symbolic Fcdf(double x, double r1, double r2)
{
	boost::math::fisher_f_distribution<> fisher_f(r1,r2);
	double cdf_value = boost::math::cdf(fisher_f,x);

	return cdf_value;
}

Symbolic Fmean(double x, double r1, double r2)
{
	Symbolic mean = divisions(r2,r2-2) ;

	return mean;
}

Symbolic Fvar(double x, double r1, double r2)
{
	Symbolic var = 2*divisions(r2,r2-2)*divisions(r2,r2-2) * divisions(r1+r2-2,r1*(r2-4));

	return var;
}

Symbolic tpdf(double x, double r)
{
	boost::math::students_t_distribution<> students_t(r);
	double pdf_value = boost::math::pdf(students_t,x);

	return pdf_value;
}

Symbolic tcdf(double x, double r)
{
	boost::math::students_t_distribution<> students_t(r);
	double cdf_value = boost::math::cdf(students_t,x);

	return cdf_value;
}

Symbolic tmean(double x, double r)
{
	Symbolic mean = 0;

	return mean;
}

Symbolic tvar(double x, double r)
{
	Symbolic var = divisions(r,r-2);

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

#include <vector>
#include <map>
#include <algorithm> // For std::max_element,  std::sort
#include <numeric> // For std::accumulate

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
	cout << "Standard deviation : " ;
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

#endif
#endif