// g++ -o result main.cpp -lcblas -lgsl
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector> // For std::vector (example container)
#include <algorithm> // For std::sort
#include <chrono>

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

using namespace std::chrono;
using namespace std;

#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>

double beta_mgf(double alpha, double beta, double t) 
{
	return boost::math::hypergeometric_1F1(alpha, alpha + beta, t);
}

long double H_1F1( double a, double b, double z)
{
	long double Fa = gsl_sf_hyperg_1F1(a, b, z);
	return Fa;
}


int main() 
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	double t = 0.5;
	double alpha = 2.0;
	double beta = 3.0;

	double result = H_1F1(alpha, beta, t);
	printf("GSL: Hypergeometric 1F1 (%.1f, %.1f, %.1f): %e\n", alpha, beta, t, result);

	//printf("Boost: Hypergeometric 1F1 (%.1f, %.1f, %.1f): %e\n", alpha, beta, t, boost::math::hypergeometric_1F1(alpha,  beta, t));
	
	double mgf_result = beta_mgf(alpha, beta, t);
	//printf("Boost: Hypergeometric 1F1 (%.1f, %.1f, %.1f): %e\n", alpha, beta, t, boost::math::hypergeometric_1F1(alpha,  beta, t));

	//cout << "Boost: Hypergeometric 1F1(" << alpha << ", " << beta<< ") at t=" << t << " is: " << boost::math::hypergeometric_1F1(alpha,  beta, t) << endl;

	cout << "Boost: MGF of Beta(" << alpha << ", " << beta<< ") at t=" << t << " is: " << beta_mgf(alpha,  beta, t) << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	
        return 0;
    }
