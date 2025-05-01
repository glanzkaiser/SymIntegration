// For direct computation, the series representation can be used, but care must be taken regarding convergence and numerical stability:
// This direct computation method is suitable for |z| < 1. 
// For other cases, analytic continuation or other representations of the hypergeometric function are needed. 
// Libraries like Boost Math or specialized numerical software are recommended for robust and accurate computation
//  across a wider range of parameters.

#include <cmath>
#include <iostream>
#include <stdexcept>
 #include "symintegrationc++.h"

using namespace std;

double division(double x, double y)
{
	return x/y;
}

double pochhammer(double x, int n) {
	if (n < 0) 
	{
		throw invalid_argument("n must be non-negative");
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

Symbolic hypergeometric_2F1_direct(double a, double b, double c, const Symbolic &s, double tolerance = 1e-10, int max_iterations = 5) {
	Symbolic sum = 1;
	Symbolic term = 1;
	for (int n = 1; n <= max_iterations; ++n) 
	{
		term *= (pochhammer(a, n) * pochhammer(b, n) * (s^(n))) / (pochhammer(c, n) * tgamma(n + 1));
		sum += term;
	}
	return sum;
}


int main() {
	Symbolic x("x"), y_sin;

	int n = 3;
	
	y_sin = -(sin(x)^(n+1))*(cos(x))*((sin(x)^(2))^(-0.5*n-0.5));
	double a = 0.5;
	double b = (1-n)*0.5;
	double c = 1.5;
	Symbolic z = cos(x)*cos(x);

	try 
	{
		Symbolic result = hypergeometric_2F1_direct(a, b, c, z);
		cout << "Hypergeometric_2F1_direct(" << a << ", " << b << ", " << c << ", " << z << ") = " << result << endl;
	} 
	catch (const exception& e) {
		cerr << "Error: " << e.what() << endl;
	}

	cout << "\nfor n = 3" << endl;
	cout << "\nintegral (Sin^{n} x) = " << y_sin*hypergeometric_2F1_direct(a, b, c, z) << endl;
	
	n=5;
	y_sin = -(sin(x)^(n+1))*(cos(x))*((sin(x)^(2))^(-0.5*n-0.5));
	b = (1-n)*0.5;
	cout << "\nfor n = 5" << endl;
	cout << "\nintegral (Sin^{n} x) = " << y_sin*hypergeometric_2F1_direct(a, b, c, z) << endl;

	n=7;
	y_sin = -(sin(x)^(n+1))*(cos(x))*((sin(x)^(2))^(-0.5*n-0.5));
	b = (1-n)*0.5;
	cout << "\nfor n = 7" << endl;
	cout << "\nintegral (Sin^{n} x) = " << y_sin*hypergeometric_2F1_direct(a, b, c, z) << endl;
	
	n=8;
	y_sin = -(sin(x)^(n+1))*(cos(x))*((sin(x)^(2))^(-0.5*n-0.5));
	b = (1-n)*0.5;
	cout << "\nfor n = 8" << endl;
	cout << "\nintegral (Sin^{n} x) = " << y_sin*hypergeometric_2F1_direct(a, b, c, z) << endl;
	
	return 0;
}