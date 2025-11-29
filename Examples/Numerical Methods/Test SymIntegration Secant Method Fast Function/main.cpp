// Merci beaucoup Freya et Sentinel
// g++ main.cpp -o result -lsymintegration
#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include "symintegrationc++.h"

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f
#define pi  3.1415926535897

#include <chrono>
#include <string>
#include <bitset>

#include <iostream>
#include <cmath>
#include <iomanip>
#include <stdexcept>


using namespace std::chrono;
using namespace std;

// Define the function f(x) whose root we want to find
double f(double x) 
{
	// Example function: f(x) = x^3 + x - 1
	return pow(x, 3) + x - 1;
}
// Secant method function
double secant(double x0, double x1, double tolerance, int maxIterations) 
{
	double x_new, fx0, fx1;
	
	for (int i = 0; i < maxIterations; ++i) 
	{
		fx0 = f(x0);
		fx1 = f(x1);
		// Check for division by zero
		if (abs(fx1 - fx0) < 1e-10) 
		{
			throw runtime_error("Secant method: Division by zero (f(x1) == f(x0))");
		}

		// Apply the secant method formula: 
		// x_new = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
		x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0);

		// Check for convergence (stopping criterion)
		if (abs(x_new - x1) < tolerance) 
		{
			return x_new; // Root found within the desired tolerance
		}

		// Update values for the next iteration
		x0 = x1;
		x1 = x_new;
		}

		// If the loop finishes without converging, throw an exception
		throw runtime_error("Secant method: Did not converge within max iterations");
}

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	
	// Set output precision
	cout << fixed << setprecision(6);
	// Initial guesses and tolerance
	double initial_x0 = 0.0;
	double initial_x1 = 1.0;
	double tolerance = 0.0001;
	int maxIter = 20;

	try 
	{
		double root = secant(initial_x0, initial_x1, tolerance, maxIter);
		cout << "Root of the equation: " << root << endl;
	} 
	catch (const runtime_error& e) 
	{
		cerr << "Error: " << e.what() << endl;
	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}