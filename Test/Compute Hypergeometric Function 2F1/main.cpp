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

Symbolic pochhammer(Symbolic x, int n) {
	if (n < 0) 
	{
		throw invalid_argument("n must be non-negative");
	}
	if (n == 0) 
	{
		return 1.0;
	}
	Symbolic result = 1.0;
	for (int i = 0; i < n; ++i) 
	{
		result *= (x + i);
	}
	return result;
}

Symbolic hypergeometric_2F1_direct(Symbolic a, Symbolic b, Symbolic c, const Symbolic &s, double tolerance = 1e-10, int max_iterations = 3) {
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
	Symbolic x("x"), y_sin, y_cos, y_tan, y_cot, y_sec, y_csc;

	int n = 3;
	
	y_sin = -(sin(x)^(n+1))*(cos(x))*((sin(x))^(-n-1));
	Symbolic a = 0.5;
	Symbolic b = (1-n)*0.5;
	Symbolic c = 1.5;
	Symbolic z = cos(x)*cos(x);

	y_cos = -(sin(x))*(cos(x)^(n+1))*(csc(x)) / (n+1);
	Symbolic a_cos = 0.5;
	Symbolic b_cos = (1+n)*0.5;
	Symbolic c_cos = (n+3)*0.5;
	Symbolic z_cos = cos(x)*cos(x);

	y_tan = (tan(x)^(n+1)) / (n+1);
	Symbolic a_tan = 1;
	Symbolic b_tan = (1+n)*0.5;
	Symbolic c_tan = (n+3)*0.5;
	Symbolic z_tan = -tan(x)*tan(x);
	
	y_cot = -(cot(x)^(n+1)) / (n+1);
	Symbolic a_cot = 1;
	Symbolic b_cot = (1+n)*0.5;
	Symbolic c_cot = (n+3)*0.5;
	Symbolic z_cot = -cot(x)*cot(x);
	
	y_sec = (sin(x))*(sec(x)^(n+1))*((cos(x))^(n+1));
	Symbolic a_sec = 0.5;
	Symbolic b_sec = (1+n)*0.5;
	Symbolic c_sec = 1.5;
	Symbolic z_sec = sin(x)*sin(x);

	y_csc = -(cos(x))*(csc(x)^(n-1))*((sin(x))^(n-1)) ;
	Symbolic a_csc = 0.5;
	Symbolic b_csc = (1+n)*0.5;
	Symbolic c_csc = 1.5;
	Symbolic z_csc = cos(x)*cos(x);
	
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
	cout << "\nintegral (Cos^{n} x) = " << y_cos*hypergeometric_2F1_direct(a_cos, b_cos, c_cos, z_cos) << endl;
	cout << "\nintegral (Tan^{n} x) = " << y_tan*hypergeometric_2F1_direct(a_tan, b_tan, c_tan, z_tan) << endl;
	cout << "\nintegral (Cot^{n} x) = " << y_cot*hypergeometric_2F1_direct(a_cot, b_cot, c_cot, z_cot) << endl;
	cout << "\nintegral (Sec^{n} x) = " << y_sec*hypergeometric_2F1_direct(a_sec, b_sec, c_sec, z_sec) << endl;
	cout << "\nintegral (Csc^{n} x) = " << y_csc*hypergeometric_2F1_direct(a_csc, b_csc, c_csc, z_csc) << endl;

	/*
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
	*/	
	return 0;
}