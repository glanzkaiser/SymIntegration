// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector> // For std::vector (example container)
#include "symintegrationc++.h"
#include <algorithm> // For std::sort

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

using namespace std;

double division(double x, double y)
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

// Driver code
int main(int argc, char** argv)
{
	int n;

    	cout << "Enter the number of terms: ";
    	cin >> n;

	cout << "Rising pochhammer (6)n :  " << rising_pochhammer(6,n) << endl ;
	cout << "Falling pochhammer (6)n :  " << falling_pochhammer(6,n) << endl ;
	return 0;
}