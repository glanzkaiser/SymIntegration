/*
   
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_SPECIALFUNCTIONS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_SPECIALFUNCTIONS_DEFINE
#define π 3.1415926535897f

Symbolic hypergeometric_1F1(double a, double b, const Symbolic &s, int max_iterations) {
	Symbolic sum = 1;
	
	for (int i = 1; i <= max_iterations; ++i) 
	{
		sum += (rising_pochhammer(a, i) * (s^(i))) / (rising_pochhammer(b, i) * factorial(i));
	}
	return sum;
}

double hypergeometric_1F1(double a, double b, double s, int max_iterations) {
	Symbolic sum = 1;
	
	for (int i = 1; i <= max_iterations; ++i) 
	{
		sum += (rising_pochhammer(a, i) * (s^(Symbolic(i)))) / (rising_pochhammer(b, i) * factorial(i));
	}
	return sum;
}
#endif
#endif