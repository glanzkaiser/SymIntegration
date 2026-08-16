/*
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


// dsolve.h

#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DECLARE

#include "polynomial.h"
#include "rational.h"

double roundToDecimal(double, int); 

Symbolic dsolve(const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic dsolve(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic dsolvelogistic(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic dsolveseparable(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);

void secondorderlineardiffeq_dsolve(double, double, double, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_ivpsolution(double, double, double, const Symbolic &, const Symbolic &, double, double, double);
void secondorderlineardiffeq_springmasssystem(double, double, double, const Symbolic &, const Symbolic &, double, double);
void secondorderlineardiffeq_RLCserieselectriccircuit(double, double, double, double, double);

void wronskian(double, double, double, const Symbolic &, const Symbolic &, double);
Symbolic wronskian_resultonly(double, double, double, const Symbolic &, const Symbolic &);
void wronskian_fundamentalsetofsolutions(double, double, double, const Symbolic &, const Symbolic &);

void reductionoforder(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &); 

void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationsivpsolution(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double, double, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationsivpforcedvibrationssolution(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double, double, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &, const Symbolic &, const Symbolic &, Polynomial<double> &);
void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &, const Symbolic &, const Symbolic &, const SymbolicMatrix &, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationssolution_variationofparameters(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);

vector<complex<double>> higherorderlineardiffeq_vectorize(const Symbolic &,const Symbolic &,const Symbolic &, int);

void higherorderlineardiffeq_homogeneousequationsivpsolution(const vector<complex<double>> &,  const vector<complex<double>> &, const vector<complex<double>> &,  int);
void higherorderlineardiffeq_homogeneousequationsivpsolution(const vector<complex<double>> &,  const vector<complex<double>> &);
void higherorderlineardiffeq_twospringtwomasssystem(double, double, double, double, Symbolic &, const vector<complex<double>> &);
void higherorderlineardiffeq_nonhomogeneousequationsgeneralsolution(const vector<complex<double>> &, const Symbolic &, const Symbolic &);
void higherorderlineardiffeq_nonhomogeneousequationsgeneralsolution(const vector<complex<double>> &, const SymbolicMatrix &, const Symbolic &);
void higherorderlineardiffeq_nonhomogeneousequations_undeterminedcoefficients(const Symbolic &, const Symbolic &, const Symbolic &);
void higherorderlineardiffeq_nonhomogeneousequations_variationofparameters(const vector<complex<double>> &, const Symbolic &, const Symbolic &);

#ifndef HIGHERORDERODE_HOMOGENEOUS_POWERSERIESSOLVER_H
#define HIGHERORDERODE_HOMOGENEOUS_POWERSERIESSOLVER_H

class HigherOrderODE_Homogeneous_PowerSeriesSolver {
private:
	vector<double> ode_coefficients; // ODE coefficients [b_0, b_1, ..., b_m]
	vector<double> initial_values; // Initial conditions [y(0), y'(0), ..., y^(m-1)(0)]
	int order;                        // Order of the differential equation (m)
	vector<double> coefficients;

	// Helper function to calculate Pochhammer symbol / falling factorial factors
	double get_factorial_multiplier(int n, int k) const 
	{
		double multiplier = 1.0;
		for (int j = 1; j <= k; ++j) 
		{
			multiplier *= (n + j);
		}
		return multiplier;
	}
public:
	// Constructor initializes the terms (series degree+1) and base boundary conditions
	HigherOrderODE_Homogeneous_PowerSeriesSolver(const vector<double>& ode_coeffs, const vector<double>& init_conditions) ;

	// Computes coefficients up to order N using the recurrence relation:
	void computeSeries(int terms) ;
	vector<double> coefficientsvector(int terms) ;
	double evaluateAt(double x, int terms) ;
	void printSeries() const ;
};

#endif

// Helper structure to compute combinations for polynomial shifting
long long binomialCoefficient(int n, int k) 
{
	if (k < 0 || k > n) return 0;
	if (k == 0 || k == n) return 1;
	long long res = 1;
	for (int i = 1; i <= k; ++i) 
	{
		res = res * (n - i + 1) / i;
	}
	return res;
}

// Represents a polynomial coefficient (e.g., P(x) = p0 + p1*x + p2*x^2 + ...)
// Polynomialcoeff P({p0, p1, p2, p3, ...});
// A basic polynomial class to handle shifting and evaluation
class Polynomialcoeff {
public:
	vector<double> coeffs; // coeffs[i] is the coefficient of x^i

	Polynomialcoeff() : coeffs({0.0}) {}
	Polynomialcoeff(const std::vector<double>& c) : coeffs(c) {
		trim();
	}

	// Remove trailing zeros to maintain accurate degree
	void trim() 
	{
		while (coeffs.size() > 1 && std::abs(coeffs.back()) < 1e-9) 
		{
			coeffs.pop_back();
		}
	}

	// Shift polynomial to be centered around x0 using Taylor expansion / Horner's scheme
	// Returns a new polynomial P(X) where X = x - x0
	Polynomialcoeff shift_around(double x0) const 
	{
		int n = coeffs.size();
		vector<double> shifted(n, 0.0);
		vector<double> temp = coeffs;

		// Synthetic division (Horner's method variant) to find Taylor coefficients
		for (int i = 0; i < n; ++i) 
		{
			for (int j = n - 1; j > i; --j) 
			{
				temp[j - 1] += temp[j] * x0;
			}
			shifted[i] = temp[i];
		}
		return Polynomialcoeff(shifted);
	}

	// Alternative to shift_around, will produce same result anyway
	// Shift polynomial from variable x to variable (x - x0)
	// P(x) -> P(u + x0) where u = x - x0
	Polynomialcoeff shiftToCenter(double x0) const 
	{
		int deg = static_cast<int>(coeffs.size()) - 1;
		vector<double> shifted(deg + 1, 0.0);

		for (int i = 0; i <= deg; ++i) 
		{
			double c = coeffs[i];
			if (std::abs(c) < 1e-14) 
			{	
				continue;
			}
			// Expand c * (u + x0)^i using Binomial Theorem
			for (int k = 0; k <= i; ++k) 
			{
				shifted[k] += c * binomialCoefficient(i, k) * std::pow(x0, i - k);
			}
		}
		return Polynomialcoeff(shifted);
	}

	// Alternative to shift_around and shiftToCenter, will produce the same result
	// Shifts the polynomial P(x) into P(t + x0) where t = (x - x0)
	// This allows us to work directly with power series centered at x0
	Polynomialcoeff shift(double x0) const 
	{
		int deg = static_cast<int>(coeffs.size()) - 1;
		vector<double> shifted(deg + 1, 0.0);
        
		// Using binomial expansion to shift coordinates
		for (int i = 0; i <= deg; ++i) 
		{
			double c = coeffs[i];
			vector<double> term(i + 1, 0.0);
			term[0] = 1.0;
            
			// Generate (t + x0)^i
			for (int j = 0; j < i; ++j) 
			{
				vector<double> next_term(j + 2, 0.0);
				for (int k = 0; k <= j; ++k) 
				{
					next_term[k] += term[k] * x0;      // multiply by x0
					next_term[k + 1] += term[k];       // multiply by t
				}
				term = next_term;
			}

			for (int j = 0; j <= i; ++j) 
			{
				shifted[j] += c * term[j];
			}
		}
		return Polynomialcoeff(shifted);
	}

	double get_coeff(size_t degree) const 
	{
		if (degree < coeffs.size()) 
		{
			return coeffs[degree];
		}
		return 0.0;
	}
	
	// to obtain the sum of the polynomial coefficients
	double sum_coeff() const
	{
		int n = coeffs.size();
		double sum = 0;
		for (int i = 0; i < n; ++i) 
		{
			sum += coeffs[i];
		}
		return sum;

	}

	// Get the maximum stored degree
	int maxDegree() const 
	{
		return coeffs.empty() ? 0 : coeffs.size() - 1;
	}

	// Evaluates the polynomial at a specific value using Horner's method
	double evaluateAt(double x) const 
	{
		double result = 0.0;
		for (int i = coeffs.size() - 1; i >= 0; --i) 
		{
			result = result * x + coeffs[i];
		}
		return result;
	}
};

#ifndef SECONDORDERODE_HOMOGENEOUS_POWERSERIESSOLVER_H
#define SECONDORDERODE_HOMOGENEOUS_POWERSERIESSOLVER_H
class SecondOrderODE_Homogeneous_PowerSeriesSolver {
private:
	Polynomialcoeff P; // Coefficient for y''
	Polynomialcoeff Q; // Coefficient for y'
	Polynomialcoeff R; // Coefficient for y
	double x0;
 	double y0; // for y(0)
	double dy0; // for y'(0)
	vector<double> coefficients; // Power series coefficients [c_0, c_1, c_2, ...]

public:
	SecondOrderODE_Homogeneous_PowerSeriesSolver(
	const Polynomialcoeff& P_input, 
	const Polynomialcoeff& Q_input, 
	const Polynomialcoeff& R_input, 
	double x0_input, 
	double y0_input, 
	double dy0_input) ;

	// Initialize solver with ODE component polynomials and targeted precision order
//	SecondOrderODE_Homogeneous_PowerSeriesSolver(Polynomialcoeff p_coeff, Polynomialcoeff q_coeff, Polynomialcoeff r_coeff, int max_order)
//		: P(p_coeff), Q(q_coeff), R(r_coeff), terms(terms_input) {
//		coefficients.resize(terms, 0.0);
//		}
	void computeSeries(int terms) ;
	vector<double> coefficientsvector(int terms) ;
	double evaluateAt(double x, int terms) ;
	void printCoefficients() const ;
	void printSolution() const ;
	
};

#endif

Symbolic ivp(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);


#endif
#endif


#endif
