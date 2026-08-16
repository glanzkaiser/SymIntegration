// g++ -o result main.cpp -lsymintegration 
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include "symintegrationc++.h"

#include <chrono>

using namespace std::chrono;
using namespace std;

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/fisher_f.hpp>

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	double degrees_of_freedom = 11.0; // Replace with your desired degrees of freedom
	boost::math::chi_squared_distribution<> my_chi_squared(degrees_of_freedom);

	double x_value = 11.264; // The value at which to evaluate the CDF
	double cdf_value = boost::math::cdf(my_chi_squared, x_value);
	cout << "\nchisquaredcdf(11.264,11) = " << cdf_value << endl;

	boost::math::gamma_distribution<> my_gamma(degrees_of_freedom*0.5,2);
	double cdf_value2 = boost::math::cdf(my_gamma, x_value);
	cout << "\nchisquaredcdf(11.264,11) = " << cdf_value2 << endl;

	double alpha = 1, beta = 10;
	boost::math::beta_distribution<> my_beta(alpha,beta);
	double cdf_value3 = boost::math::cdf(my_beta,0.2);
	cout << "\nbetacdf(0.2;1,10) = " << cdf_value3 << endl;


	boost::math::fisher_f_distribution<> my_f(6,10);
	double cdf_value4 = boost::math::pdf(my_f,3.22);
	cout << "\nFcdf(0.2;1,10) = " << cdf_value4 << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}