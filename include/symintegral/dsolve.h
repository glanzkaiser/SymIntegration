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

long long binomialCoefficient(int, int );

// Represents a polynomial coefficient (e.g., P(x) = p0 + p1*x + p2*x^2 + ...)
// PolynomialDouble P({p0, p1, p2, p3, ...});
// A basic polynomial class to handle shifting and evaluation
class PolynomialDouble {
public:
	vector<double> coeffs; // coeffs[i] is the coefficient of x^i
	
	PolynomialDouble() : coeffs({0.0}) {}
	PolynomialDouble(const std::vector<double>& c) : coeffs(c) {
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
	PolynomialDouble shift_around(double x0) const 
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
		return PolynomialDouble(shifted);
	}

	// Alternative to shift_around, will produce same result anyway
	// Shift polynomial from variable x to variable (x - x0)
	// P(x) -> P(u + x0) where u = x - x0
	PolynomialDouble shiftToCenter(double x0) const 
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
		return PolynomialDouble(shifted);
	}

	// Alternative to shift_around and shiftToCenter, will produce the same result
	// Shifts the polynomial P(x) into P(t + x0) where t = (x - x0)
	// This allows us to work directly with power series centered at x0
	PolynomialDouble shift(double x0) const 
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
		return PolynomialDouble(shifted);
	}

	// Returns the lowest power with a non-zero coefficient (Valuation)
	int valuation() const 
	{
		for (int i = 0; i < coeffs.size(); ++i) 
		{
			if (std::abs(coeffs[i]) > 1e-9)
			{ 
				return i;
			}
		}
		return coeffs.size(); // default fallback
	}

	double evaluate(double x) const 
	{
		double result = 0;
		for (int i = coeffs.size() - 1; i >= 0; --i) 
		{
			result = result * x + coeffs[i];
		}
		return result;
	}

	// Evaluates the i-th derivative at a given point (used for shifting / Taylor expansion)
	double evaluate_derivative(int order, double x0) const 
	{
		double val = 0.0;
		for (int i = order; i < coeffs.size(); ++i) 
		{
			double term = coeffs[i];
			for (int j = 0; j < order; ++j) 
			{
				term *= (i - j);
			}
			val += term * std::pow(x0, i - order);
		}
		return val;
	}

	// Shift polynomial to be centered around x0 -> returns a new Polynomial in terms of t = (x - x0)
	PolynomialDouble shiftFrobeniuswithderivative(double x0) const 
	{
		vector<double> shifted_coeffs;
		double factorial = 1.0;
		for (int i = 0; i < coeffs.size(); ++i) 
		{
			if (i > 0) 
			{
				factorial *= i;			
			}
			double c = evaluate_derivative(i, x0) / factorial;
			shifted_coeffs.push_back(c);
		}
		return PolynomialDouble(shifted_coeffs);
	}

	// Alternative to shift, shift_around and shiftToCenter, will produce the same result
	// Shifts the polynomial P(x) to p(t) where t = x - x0 -> P(t + x0)
	PolynomialDouble shiftFrobenius(double x0) const 
	{
		int n = coeffs.size();
		vector<double> shifted(n, 0.0);
        
		// Use Taylor expansion / Horner's scheme to shift center
		vector<double> temp = coeffs;
		for (int i = 0; i < n; ++i) 
		{
			shifted[i] = temp[0];
			for (size_t j = 0; j < temp.size() - 1; ++j) 
			{
				temp[j] = temp[j + 1] * (j + 1);
			}
			temp.pop_back();
			double factorial = 1.0;
			for (int k = 1; k <= i; ++k) 
			{	
				factorial *= k;
			}
			shifted[i] /= factorial;
		}

		// Exact shifting using binomial theorem for accuracy
		vector<double> result(n, 0.0);
		for (int i = 0; i < n; ++i) 
		{
			double c = coeffs[i];
			// Expand c * (t + x0)^i
			for (int j = 0; j <= i; ++j) 
			{
				double binomial = 1.0;
				for (int k = 0; k < j; ++k) 
				{
					binomial *= (i - k) / (k + 1.0);
				}
				// combination(i, j) * x0^(i-j) * t^j
				double term = c * std::pow(x0, i - j);
				// Compute combination manually or via standard loops
			}
		}
        
		// Alternative clean structural shift:
		vector<double> current = {1.0};
		vector<double> total(n, 0.0);
		for(int i = 0; i < coeffs.size(); ++i) 
		{
			for(int j = 0; j < current.size(); ++j) 
			{
				total[j] += coeffs[i] * current[j];
			}
			// Multiply current by (t + x0)
			vector<double> next(current.size() + 1, 0.0);
			for(int j = 0; j < current.size(); ++j) 
			{
				next[j] += current[j] * x0;
				next[j+1] += current[j];
			}
			current = next;
		}
		return PolynomialDouble(total);
	}

	double get_coeff(size_t degree) const 
	{
		if (degree < coeffs.size()) 
		{
			return coeffs[degree];
		}
		return 0.0;
	}
	// Getter to access the coefficients
	const vector<double>& getCoeffs() const 
	{ 
		return coeffs; 
	}
	
	// to obtain the absolute sum of the polynomial coefficients
	double sum_coeff() const
	{
		int n = coeffs.size();
		double sum = 0;
		for (int i = 0; i < n; ++i) 
		{
			sum += abs(coeffs[i]);
		}
		return sum;

	}

	// Get the maximum stored degree
	int maxDegree() const 
	{
		return coeffs.empty() ? 0 : coeffs.size() - 1;
	}

	// Evaluates the polynomial at a specific value using Horner's method
	 // Iterate backwards from the second-to-last coefficient down to index 0
	double evaluateAt(double x) const 
	{
		double result = 0.0;
		for (int i = coeffs.size() - 1; i >= 0; --i) 
		{
			result = result * x + coeffs[i];
		}
		return result;
	}

	//  Print function for visualization
	void print() const 	
	{
		bool first = true;
		for (int i = maxDegree(); i >= 0; --i) 
		{
			if (coeffs[i] == 0 && maxDegree() > 0) 
			{
				continue;
			}
			if (!first && coeffs[i] > 0) 
			{
				cout << " + ";
			}
			else if (coeffs[i] < 0) 
			{
				cout << " - ";
			}

			double absVal = std::abs(coeffs[i]);
			if (absVal != 1.0 || i == 0) 
			{
				cout << absVal;
			}

			if (i > 0) 
			{
				cout << "x";
				if (i > 1)
				{
					cout << "^" << i;
				}
			}
		
		first = false;
		}
		cout << endl;
	}

};

class PolynomialComplex{
public:
	// Coefficients stored from lowest degree to highest degree
	// coeffs[i] represents the coefficient for x^i
	vector<complex<double>> coeffs;

	PolynomialComplex() : coeffs({0.0}) {}

	// Constructor from std::vector<std::complex<double>>
	PolynomialComplex(const vector<complex<double>>& c) : coeffs(c) { trim(); }
	
	// Initializer list constructor for easy syntax: Polynomial p({{1,0}, {2,1}, {3,0}});
	PolynomialComplex(std::initializer_list<complex<double>> list) : coeffs(list) 
	{
		if (coeffs.empty()) 
		{
			coeffs.push_back(0.0);
		}
		trim();
	}

	int size()
	{
		return coeffs.size();	
	}

	bool empty()
	{
		coeffs.empty();
	}
	// Helper to remove trailing zeros and keep the vector compact
	void trim() 
	{
		while (!coeffs.empty() && std::norm(coeffs.back()) < 1e-15) 
		{
			coeffs.pop_back();
		}
	}

	// Conversion Constructor: Creates a complex polynomial from a double polynomial
	PolynomialComplex(const PolynomialDouble& other) 
	{
		const auto& double_coeffs = other.getCoeffs();
		
		// Reserve memory to prevent multiple allocations
		coeffs.reserve(double_coeffs.size()); 
		
		// Convert and copy each double to std::complex<double>
		for (double c : double_coeffs) 
		{
			coeffs.emplace_back(c, 0.0); // Real part = c, Imaginary part = 0.0
		}
	}

	// Shift function: Computes P(r + shift)
	// Uses binomial expansion to calculate the new coefficients
	PolynomialComplex shift(complex<double> s) const 
	{
		int deg = maxDegree();
		if (deg < 0) 
		{
			return PolynomialComplex();
		}

		vector<complex<double>> new_coeffs(deg + 1, 0.0);

		// Standard binomial expansion for shifting: 
		// For each term a_i * x^i, replace x with (x + s) -> a_i * \sum binom(i, k) * x^k * s^(i-k)
		for (int i = 0; i <= deg; ++i) 
		{
			complex<double> ai = coeffs[i];
			if (std::norm(ai) < 1e-15) 
			{
				continue;
			}

			// Compute binomial coefficients on the fly
			complex<double> s_pow = 1.0; 
			double binom = 1.0;

			for (int k = i; k >= 0; --k) 
			{
				new_coeffs[k] += ai * binom * s_pow;
		        
				// Update binomial for next iteration (moving down to k-1)
				if (k > 0) 
				{
					binom = binom * k / (i - k + 1);
					s_pow *= s;
				}
			}
		}

		return PolynomialComplex(new_coeffs);
	}

	// Shift polynomial P(x) to P(u + x0) using Taylor expansion / synthetic division
	PolynomialComplex shiftTaylor(complex<double> x0) const 
	{
		int n = maxDegree();
		vector<complex<double>> result = coeffs;
		for (int i = 0; i <= n; ++i) 
		{
			for (int j = n - 1; j >= i; --j) 
			{
				result[j] += result[j + 1] * x0;
			}
		}
		return PolynomialComplex(result);
	}

	// Basic Addition
	PolynomialComplex Add(const PolynomialComplex& other) const 
	{
		int max_size = std::max(coeffs.size(), other.coeffs.size());
		vector<complex<double>> result_coeffs(max_size, 0.0);

		for (size_t i = 0; i < max_size; ++i) 
		{
			result_coeffs[i] = get_coeff(i) + other.get_coeff(i);
		}
		return PolynomialComplex(result_coeffs);
	}

	// Safely get coefficient for x^pow
	complex<double> get_coeff(size_t pow) const 	
	{
		if (pow < coeffs.size()) 
		{
			return coeffs[pow];
		}
		return 0.0;
	}

	// Get the degree of the polynomial
	int maxDegree() const 
	{
		return coeffs.empty() ? 0 : coeffs.size() - 1;
	}

	// Horner's Method to evaluate the polynomial at a given complex point z
	complex<double> evaluateAt(const std::complex<double>& z) const 
	{
		complex<double> result = 0.0;
        
		// Loop backwards from the highest degree coefficient to the lowest
		for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it) 
		{
			result = result * z + *it;
		}
        
		return result;
	}
	
	// Durand-Kerner method to find all complex roots
	vector<complex<double>> find_roots_durand_kerner() const 
	{
		int max_iterations = 1000; 
		double tolerance = 1e-12;
		size_t n = maxDegree();
		if (n == 0) return {};

		// Normalize leading coefficient (make it monic)
		complex<double> leading_coeff = coeffs.back();
		
		// Initialize roots with distinct complex values (Aberth/Guggenheimer initial values)
		vector<complex<double>> roots(n);
		Complex center(0.4, 0.9); // Arbitrary complex starting point shifts
		for (size_t i = 0; i < n; ++i) 
		{
			double angle = 2.0 * M_PI * i / n;
			roots[i] = center + std::pow(complex<double>(0.9, 0.2), i); // Generates distinct initial guesses
		}

			bool converged = false;
		for (int iter = 0; iter < max_iterations && !converged; ++iter) 
		{
			converged = true;
			vector<Complex> old_roots = roots;

			for (size_t i = 0; i < n; ++i) 
			{
				// Calculate the denominator: Product(x_i - x_j) for j != i
				complex<double> denominator = 1.0;
				for (size_t j = 0; j < n; ++j) 
				{
					if (i != j) 
					{
						denominator *= (old_roots[i] - old_roots[j]);
					}
				}

				// Durand-Kerner formula step: x_i = x_i - P(x_i) / (a_n * denominator)
				complex<double> numerator = evaluateAt(old_roots[i]);
				complex<double> correction = numerator / (leading_coeff * denominator);
		        
				roots[i] = old_roots[i] - correction;

				// Check if any root is still changing above tolerance limit
				if (std::abs(correction) > tolerance) 
				{
					converged = false;
				}
			}
		}
		for (size_t i = 0; i < n; ++i) 
		{
			if(abs(real(roots[i])) < 1e-6 )
			{
				//roots[i] = (0,imag(roots[i]));
			}
			if(abs(imag(roots[i])) < 1e-6 )
			{
				//roots[i] = (real(roots[i]),0);
			}
		}
		return roots;
	}

	// Synthetic division: divides polynomial by (x - root)
	// Returns pair: {quotient, remainder}
	pair<PolynomialComplex, complex<double>> divide_by_linear(Complex root) const 
	{
	 	if (coeffs.empty()) 
		{
			return {PolynomialComplex({0}), 0};
		}
		
		int n = maxDegree();
		vector<complex<double>> q_coeffs(n > 0 ? n : 1, 0.0);
		Complex rem = coeffs.back();
		
		if (n > 0) 
		{
			q_coeffs[n - 1] = coeffs.back();
			for (int i = n - 2; i >= 0; --i) 
			{
				q_coeffs[i] = coeffs[i + 1] + q_coeffs[i + 1] * root;
			}
			rem = coeffs[0] + q_coeffs[0] * root;
		}
		
		return {PolynomialComplex(q_coeffs), rem};
	}

	// A pure polynomial is fundamentally continuous everywhere in C
	bool isAlwaysContinuous() const 
	{
		return true;
	}
	// Print utility
	void print() const 
	{
		for (size_t i = 0; i < coeffs.size(); ++i) 
		{
			cout << coeffs[i] << " * x^" << i;
			if (i < coeffs.size() - 1) 
			{
				cout << " + ";
			}
		}
		cout << endl;
	}
	/*
	// Print utility alternative
	void print() const 
	{
		if (maxDegree() < 0) 
		{
			cout << "0";
			return;
		}
		for (int i = maxDegree(); i >= 0; --i) 
		{
			complex<double> c = coeffs[i];
			if (std::norm(c) < 1e-15) 
			{
				continue;
			}

			if (i != maxDegree()) 
			{
				cout << " + ";
			}
			cout << "(" << c.real() << " + " << c.imag() << "i)";
			if (i > 0)
			{ 
				cout << "x^" << i;
			}
		}
		cout << endl;
	}*/
};


#ifndef SECONDORDERODE_HOMOGENEOUS_POWERSERIESSOLVER_H
#define SECONDORDERODE_HOMOGENEOUS_POWERSERIESSOLVER_H
class SecondOrderODE_Homogeneous_PowerSeriesSolver {
private:
	PolynomialDouble P; // Coefficient for y''
	PolynomialDouble Q; // Coefficient for y'
	PolynomialDouble R; // Coefficient for y
	double x0;
 	double y0; // for y(0)
	double dy0; // for y'(0)
	vector<double> coefficients; // Power series coefficients [c_0, c_1, c_2, ...]

public:
	SecondOrderODE_Homogeneous_PowerSeriesSolver(
	const PolynomialDouble& P_input, 
	const PolynomialDouble& Q_input, 
	const PolynomialDouble& R_input, 
	double x0_input, 
	double y0_input, 
	double dy0_input) ;

	// Initialize solver with ODE component polynomials and targeted precision order
//	SecondOrderODE_Homogeneous_PowerSeriesSolver(PolynomialDouble p_coeff, PolynomialDouble q_coeff, PolynomialDouble r_coeff, int max_order)
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

#ifndef SECONDORDERODE_HOMOGENEOUS_FROBENIUS_POWERSERIESSOLVER_H
#define SECONDORDERODE_HOMOGENEOUS_FROBENIUS_POWERSERIESSOLVER_H
class SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver {
private:
	 // Shifted polynomials
	PolynomialDouble P, P_unshifted; // Coefficient for y''
	PolynomialDouble Q, Q_unshifted; // Coefficient for y'
	PolynomialDouble R, R_unshifted; // Coefficient for y
	PolynomialComplex Pc, Qc, Rc;
	
	complex<double> x0;
	int max_terms;
	complex<double> r1, r2;      // Roots of indicial equation
	vector<complex<double>> coefficients_root1; // Power series coefficients [c_0, c_1, c_2, ...] for r1
	vector<complex<double>> coefficients_root2; // Power series coefficients [c_0, c_1, c_2, ...] for r2
	complex<double> C_log;
	bool integer_diff, complex_roots, repeated_roots, noninteger_diff; // to categorize the roots
	bool p0_exists, q0_exists, EulerCauchy, Frobenius; // to classify the ode
	complex<double> get_indicial_value(complex<double> rho);// necessary or not?

public:
	SecondOrderODE_Homogeneous_Frobenius_PowerSeriesSolver(PolynomialDouble P_in, PolynomialDouble Q_in, PolynomialDouble R_in, complex<double> x0_point, int terms) 
	: max_terms(terms), x0(x0_point)
	{
		P_unshifted = P_in;
		Q_unshifted = Q_in;
		R_unshifted = R_in;
		// Shift equations so that the regular singular point is at t = 0
		P = P_in.shift_around(x0.real());
		Q = Q_in.shift_around(x0.real());
		R = R_in.shift_around(x0.real());
		solve_indicial_equation();
		classify_ode() ;
	}
	void classify_ode() ;
	void computeSeriesCoefficients() ; 

	void solve_indicial_equation() ;
	std::pair<complex<double>, complex<double>> solve_indicial_equation_inpair() ;
	complex<double> evaluate_series_y1(complex<double> x) ;
	complex<double> evaluate_series_y2(complex<double> x) ;
	complex<double> evaluate_series_derivative(const vector<complex<double>>& a, complex<double> r_val, complex<double> x) ;
	void solve_ivp(double x_init, double y_init, double dy_init, const vector<complex<double>>& test_points);
	void printCoefficients() const ;
	void printSolution() const ;
	
};

#endif

#ifndef LAPLACE_TRANSFORM
#define LAPLACE_TRANSFORM
// Class to calculate both exact (analytical) and numerical Laplace transforms
class LaplaceTransform {
private:


public:
	LaplaceTransform();
	// Helper for numerical integration: Adaptive Simpson's Rule
	double simpsonStep(const std::function<double(double)>& f, double a, double b, double epsilon, double whole, int limit) const;
	double integrate(const std::function<double(double)>& f, double a, double b, double epsilon = 1e-6) const ;

	// 1. Constant: f(t) = 1  =>  F(s) = 1 / s
	double constant(double s) const ;
	// 2. Exponential: f(t) = e^(at)  =>  F(s) = 1 / (s - a)
	double exponential(double s, double a) const ;
	// 3. Polynomial / t-power: f(t) = t^n  =>  F(s) = n! / s^(n+1)
	double tPower(double s, int n) const ;
	// 4. Sine: f(t) = sin(wt)  =>  F(s) = w / (s^2 + w^2)
	double sine(double s, double w) const ;
	// 5. Cosine: f(t) = cos(wt)  =>  F(s) = s / (s^2 + w^2)
	double cosine(double s, double w) const ;
	// 6. Hyperbolic Sine: f(t) = sinh(at)  =>  F(s) = a / (s^2 - a^2)
	double sinh(double s, double a) const ;
	// 7. Hyperbolic Cosine: f(t) = cosh(at)  =>  F(s) = s / (s^2 - a^2)
	double cosh(double s, double a) const ;
	// 8. General Numerical Laplace Transform (for arbitrary custom functions)
	// Computes integral from 0 to upper_bound (approximating infinity) of f(t) * e^(-st)
	double numericalTransform(const std::function<double(double)>& f, double s, double upperBound = 50.0) const;
};

#endif

int find_L(const PolynomialComplex&, const PolynomialComplex& , const PolynomialComplex& );

bool isInteger_sym(Symbolic);
bool is_constant(const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic symbolicfactorial(const Symbolic &);
Symbolic SymbolicLaplaceTransform(const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic SymbolicInverseLaplaceTransform(const Symbolic &, const Symbolic &, const Symbolic &); 
Symbolic ivpsolve_laplacetransform(const Equation &, const Symbolic &, const Symbolic &, const double &, const double &);

Symbolic ivp(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_derivativesvalueatx0(const Symbolic &, const Symbolic &, const Symbolic &, double, Symbolic, Symbolic);

#endif
#endif


#endif
