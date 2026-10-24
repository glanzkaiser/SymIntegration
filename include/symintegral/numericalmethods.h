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

#include "symintegral/dual.h"
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_DECLARE

double divisiond(double, double);

Symbolic bisectionmethod(const Symbolic &, const Symbolic &,double, double, int);
Symbolic newtonmethod(const Symbolic &, const Symbolic &, const Symbolic &, int);
Symbolic eulermethod(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double);
Symbolic directionfield(const Symbolic &, const Symbolic &, const Symbolic &, double, double, double, double, double, double);

double NewtonRaphson(const Symbolic &, const Symbolic &, double);
double secantmethod(const Symbolic &, const Symbolic &, double,  double, double, int);

void AberthEhrlich(const vector<complex<double>> &,  const vector<complex<double>> &,  int);

double numericaldifferentiation(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation3pointoneside(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation3pointbothsides(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation5pointoneside(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation5pointbothsides(const Symbolic &, const Symbolic &, double, double);

double trapezoidalrule(const Symbolic &, const Symbolic &, double, double);
double simpsonsrule(const Symbolic &, const Symbolic &, double, double);
double simpsonsrule38(const Symbolic &, const Symbolic &, double, double);

double richardsonextrapolation(const Symbolic &, const Symbolic &, double, double, int, int);

void gradientdescent(const Symbolic &, const Symbolic &, double, double, int);

void choleskyDecomposition(vector<vector<double>>);
void LUDecomposition(vector<vector<double>>,vector<vector<double>>,vector<vector<double>>);

void JacobiMethod(const vector<vector<double>>& , const vector<double>& , int, double); 
bool GaussSeidel(const vector<vector<double>>&, const vector<double>& , vector<double>&, double, int); 
bool SORIterativeMethod(const vector<vector<double>>&, const vector<double>&, vector<double>&,  double, double, int ); 

/*

	Stochastic Optimizers

*/

#ifndef STOCHASTIC_OPTIMIZERS_H
#define STOCHASTIC_OPTIMIZERS_H

#include <vector>
#include <string>
#include <functional> // For std::function wrapper

// --- Standardized Structural Return Packages ---
struct DMCMetrics {
	std::vector<double> best_x;
	double best_energy;
	int total_samples;
};

struct STUNMetrics {
	std::vector<double> best_x;
	double best_energy;
	int total_steps;
	int accepted_moves;
};

struct PTMetrics {
	std::vector<double> best_x;
	double best_energy;
	int swap_attempts;
	int swap_successes;
};

// --- Benchmark Core Interface Structure ---
struct BenchmarkRecord {
	std::string solver_name;
	double best_energy;
	std::vector<double> best_x;
	double execution_time_ms;
};

// --- N-Dimensional Core Optimizer Interfaces (Accepting custom function hooks) ---
DMCMetrics run_direct_monte_carlo(
	int total_samples, 
	const std::vector<double>& lower_bounds, 
	const std::vector<double>& upper_bounds, 
	std::function<double(const std::vector<double>&)> obj_func);
STUNMetrics run_stochastic_tunneling(
	int max_steps, 
	double initial_step_size, 
	double gamma, 
	double temperature, 
	const std::vector<double>& initial_base, 
	std::function<double(const std::vector<double>&)> obj_func);
PTMetrics run_parallel_tempering(
	int total_cycles, 
	int local_steps, 
	const std::vector<double>& initial_base, 
	const std::vector<double>& temperatures, 
	std::function<double(const std::vector<double>&)> obj_func);

// --- Unified Cross-Platform Benchmark Suite Function ---
std::vector<BenchmarkRecord> run_unified_stochastic_benchmark(
	int dimensions,
	const std::vector<double>& lower_bounds,
	const std::vector<double>& upper_bounds,
	const std::vector<double>& initial_base,
	int budget_evaluations,
	std::function<double(const std::vector<double>&)> obj_func
);

#endif // STOCHASTIC_OPTIMIZERS_H

/*

	Iterative Techniques in Matrix Algebra

*/
#ifndef STRUCT_CRS_MATRIX
#define STRUCT_CRS_MATRIX

// Defining struct  in .h file (header file) will make it work when we call it from .cpp file from anywhere that use SymIntegration library
// Structure to store a sparse matrix in Compressed Row Storage (CRS) format
struct CRSMatrix {
	vector<double> values;     // Non-zero elements
	vector<int> col_indices;   // Column indices of non-zero elements
	vector<int> row_ptr;       // Row pointers
	int num_rows;
	int num_cols;
};

CRSMatrix denseToCRS(const vector<vector<double>>& dense);
vector<double> spmv(const CRSMatrix& A, const std::vector<double>& x);
double dotProduct(const vector<double>& u, const vector<double>& v);
vector<double> get_Jacobi_preconditioner(const CRSMatrix& A);
vector<double> ConjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iterations);
vector<double> PreconditionedConjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iterations);
#endif

/*

	Numerical Solutions of Nonlinear Systems of Equations

*/

using SystemFunction = std::function<vector<DualClass>(const vector<DualClass>&)>;

#ifndef NONLINEAR_NEWTONSOLVER_H
#define NONLINEAR_NEWTONSOLVER_H

class NonlinearNewtonSolver {
private:
	vector<complex<double>> roots;
	SystemFunction system;
	int max_iter ;
	double tol ;
	
public:
	NonlinearNewtonSolver(int max_iter_in, double tol_in)
	: max_iter(max_iter_in), tol(tol_in)
	{}

	// Gaussian Elimination with partial pivoting for complex systems
	vector<complex<double>> solveLinearSystem(vector<vector<complex<double>>> J, vector<complex<double>> F) ;
	vector<complex<double>> solve(SystemFunction func, const vector<complex<double>>& initial_guess)  ;
	void solve_printSolution(SystemFunction func, const vector<complex<double>>& initial_guess) ;
};


#endif

typedef vector<DualClass> (*NonlinearSystem)(const vector<DualClass>&);

#ifndef QUASINEWTONSOLVER_H
#define QUASINEWTONSOLVER_H

class QuasiNewtonSolver {
private:
	vector<complex<double>> roots;
	int max_iter ;
	double tol ;
	
public:
	QuasiNewtonSolver(int max_iter_in, double tol_in)
	: max_iter(max_iter_in), tol(tol_in)
	{}

	// Gaussian Elimination with partial pivoting for complex systems
	bool solveLinearSystem(vector<vector<complex<double>>> A, vector<complex<double>> b, vector<complex<double>>& x) ;
	void solve(NonlinearSystem system, const vector<complex<double>>& x0) ;
	vector<complex<double>> solve_getvector(NonlinearSystem system, const vector<complex<double>>& x0) ;
	void printSolution() const;
};


#endif

#ifndef STEEPESTDESCENTSOLVER_H
#define STEEPESTDESCENTSOLVER_H

class SteepestDescentSolver {
private:
	vector<complex<double>> roots;
	int max_iter ;
	double tol ;
	
public:
	using SystemFunction = std::function<vector<DualClass>(const vector<DualClass>&)>;

	SteepestDescentSolver(int max_iter_in, double tol_in)
	: max_iter(max_iter_in), tol(tol_in)
	{}

	vector<complex<double>> solve(SystemFunction system, const vector<complex<double>>& x0) ;
	void solve_printSolution(SystemFunction system, const vector<complex<double>>& x0) ;
};


#endif


#ifndef PSOMULTIROOTFINDER_H
#define PSOMULTIROOTFINDER_H

class PSOMultiRootFinder {
private:
	int num_vars;        // Number of complex variables
	int pop_size;        // Swarm population size
	int max_iter;        // Maximum iterations per run
	
	double lower_bound;  // Bound for initializing search space
	double upper_bound;  // Bound for initializing search space

	vector<ComplexVector> final_roots;

	struct Particle {
	vector<double> position; // Alternating [Real(x0), Imag(x0), Real(x1), Imag(x1)...]
	vector<double> velocity;
	vector<double> best_position;
	double best_fitness;
	};

	// Converts internal real representation back to std::complex<double>
	ComplexVector toComplexVector(const vector<double>& real_vec);

	// Evaluates system and returns the total residual magnitude squared
	double evaluateFitness(const vector<double>& real_vec, const ComplexVectorFunction& system);

public:
	PSOMultiRootFinder(int vars, int pop, int iter, double lb, double ub )
	: num_vars(vars), pop_size(pop), max_iter(iter), lower_bound(lb), upper_bound(ub) {}

	// Executes a single PSO run to locate a single global minimum root
	ComplexVector findOneRoot(const ComplexVectorFunction& system, mt19937& rng) ;

	// Multiroot manager using multi-start check routine to identify unique roots
	vector<ComplexVector> findMultipleRoots(const ComplexVectorFunction& system, int num_runs);
	void findMultipleRoots_printSolution(const ComplexVectorFunction& system, int num_runs);

};


#endif


#ifndef RK4MULTIROOTSOLVER_H
#define RK4MULTIROOTSOLVER_H
// Typedef for your system function wrapper using DualClassTemplate numbers
typedef vector<DualClassTemplate<complex<double>>> (*SystemFunc)(const vector<DualClassTemplate<complex<double>>>&);

class RK4MultiRootSolver {
private:
	vector<vector<complex<double>>> initial_guesses;
	SystemFunc targetSystem;
	int steps;

public:
	RK4MultiRootSolver(vector<vector<complex<double>>> initial_guesses_in, SystemFunc targetSystem_in, int steps_in)
	: initial_guesses(initial_guesses_in), targetSystem(targetSystem_in), steps(steps_in){}

	bool solveLinearSystem(vector<vector<complex<double>>>& A, vector<complex<double>>& b, vector<complex<double>>& x) ;
	vector<complex<double>> computeDerivatives(const vector<complex<double>>& x, SystemFunc system); 
	vector<complex<double>> solve(vector<complex<double>> x0, SystemFunc system);

	void printSolution();
};


#endif

#endif
#endif


#endif
