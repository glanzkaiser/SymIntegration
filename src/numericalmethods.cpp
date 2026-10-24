/*
   
*/
#include "symintegral/symintegrationc++.h"
#include "symintegral/dual.h"
#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_DEFINE

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()

using namespace std;

double divisiond(double x, double y)
{
	return x/y;
}

Symbolic bisectionmethod(const Symbolic &f, const Symbolic &x, double a, double b, int N)
{
 	Symbolic xa, xb, xp, fa, fp;
	
	float p = a + (b-a)/2 ;

	cout << setw(6) << "iteration" << "\t\t" << "a" << "\t\t\t" << "b" << "\t\t\t" << "p" << "\t\t\t" << "f(p)" << "\n";
	for (int i = 1; i <=N; i++)
	{
		p = a + (b-a)/2 ;
		xa = a;
		xp = p;
		fa = f[x==xa] ;
		fp = f[x==xp] ;
		double fafp = fa*fp;
	
		if (fafp> 0)
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t" << b << "\t\t\t" << p << "\t\t\t" <<  f[x==xp]  << "\n";
			a = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}	
		else
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t"<< b << "\t\t\t" << p << "\t\t\t" <<  f[x==xp]  << "\n";
			b = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}
		
		cout << endl;
	}	
	cout << "solution = " ;
	return p;
}

Symbolic newtonmethod(const Symbolic &f, const Symbolic &x, const Symbolic &x0, int N)
{
 	Symbolic fd, fp, fpd, pn, p0;

	p0 = x0;
	fd = df(f,x);	
	fp = f[x==p0] ;
	fpd = fd[x==p0] ;
	cout << "\nf(x) = " << f <<endl;
	cout << "f'(x) = " << fd <<endl;

	cout << endl;
	cout << setw(6) << "n" << "\t\t" << "p_{n}"  << "\n";
	cout << setprecision(14) << setw(6) << "0" << "\t\t" << p0 << "\n";	
	for (int i = 1; i <=N; i++)
	{
		fp = f[x==p0] ;
		fpd = fd[x==p0] ;
		pn = p0 - (fp/fpd);

		cout << setprecision(14) << setw(6) << i << "\t\t" << pn << "\n";
		double err = p0-pn;
		if (abs(err) < pow(10,-5))
		{
			cout << "The procedure was successful." << endl;			
			break;
		}
		p0 = pn;
	}
	cout << "solution = "<< endl;
	return pn;
}
#define EPSILON 0.00001
// Secant method function
double secantmethod(const Symbolic &fx, const Symbolic &x, double x0, double x1, double tolerance, int maxIterations) 
{
	double x_new, fx0, fx1;
	
	for (int i = 0; i < maxIterations; ++i) 
	{
		//fx0 = f(x0);
		//fx1 = f(x1);
		fx0 = fx[x==x0];
		fx1 = fx[x==x1];
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


double NewtonRaphson(const Symbolic &f, const Symbolic &x, double x0)
{
	Symbolic fd;
	fd = df(f,x);
	Equations rules = (  SymbolicConstant::e == exp(1), SymbolicConstant::i == sqrt(-1));
	fd = fd.subst_all(rules) ;

	double root;
	double h = f[x==x0].subst_all(rules)  / fd[x==x0].subst_all(rules) ;
	while (abs(h) >= EPSILON)
	{
		h = f[x==x0].subst_all(rules) / fd[x==x0].subst_all(rules) ;
		
		// x(i+1) = x(i) - f(x) / f'(x)  
		x0 = x0 - h;
	}
	root = x0;
	return root;
}

void AberthEhrlich(const vector<complex<double>> &P, const vector<complex<double>> &vec_x0,  int N)
{
// Done on March 23rd, 2026
	int n_Polynomial = P.size();
	complex<double> nP(n_Polynomial,0.0);
	int n = vec_x0.size();
	complex<double> root(0.0,0.0);
	complex<double> c1(1.0, 0.0); // means complex number with real part 1 and imag part 0
	complex<double> c0(0.0, 0.0);
	vector<complex<double>> P_derivative;
	vector<complex<double>> vec_update;
	vector<complex<double>> vec_dummy;
	vector<complex<double>> vec_check;
	complex<double> i_derivative(1.0,0.0);

	if(n != n_Polynomial-1)	
	{
		cerr << "Error: Initial guess has to be: the number of highest order of the derivative." << endl;
	}
	for (int i = 0; i < n_Polynomial - 1; ++i)
	{
		P_derivative.push_back((nP - i_derivative)*P[i]);
		i_derivative = i_derivative + c1;
	}
	cout << "P: " << endl;
	printComplexVector(P);
	cout << "\nP': " << endl;
	printComplexVector(P_derivative);
	//cout << "\n accumulate P: " << accumulate(P.begin(), P.end(), c0) << endl;
	//cout << "\n accumulate P': " << accumulate(P_derivative.begin(), P_derivative.end(), c0) << endl;
	
	for (int i = 0; i < n; ++i)
	{
		vec_dummy.push_back(vec_x0[i]);
	}

	cout << "\nInitial guess: " << endl;
	printComplexVector(vec_dummy);
	
	for (int k = 0; k < N ; ++k)
	{
		cout <<"\niteration: " << k << endl;

		vec_check.clear();
		for (int i = 0; i < n; ++i)
		{
			vec_check.push_back(vec_dummy[i]);
		}

		for (int i = 0; i < n; ++i)
		{
			//cout <<"\ni: " << i << endl;

			// Newton step
			complex<double> P_zi(0.0,0.0);
			complex<double> P_zi_derivative(0.0,0.0);
			complex<double> i_der(1.0,0.0);
			complex<double> i_der2(2.0,0.0);
			for (int j = 0; j < n_Polynomial ; ++j)
			{
				P_zi += P[j]*pow(vec_dummy[i], nP - i_der); 

				i_der = i_der + c1;
			}		
			//cout <<" P(z_{i}) : " << P_zi << endl;
			for (int j = 0; j < n_Polynomial - 1 ; ++j)
			{
				P_zi_derivative += P_derivative[j]*pow(vec_dummy[i], nP - i_der2); 

				i_der2 = i_der2 + c1;
			}
			//cout <<" P'(z_{i}) : " << P_zi_derivative << endl;

			complex<double> N_zi =P_zi/P_zi_derivative;
			//cout <<" N(z_{i}) : " << N_zi << endl;

			// Shift
			complex<double> shift(0.0,0.0);
			for (int j = 0; j < n ; ++j)
			{
				if(i != j)
				{
					shift += c1/(vec_dummy[i]-vec_dummy[j]);
				}
				else if(i == j)
				{
					shift += 0;
				}
			}
			//cout <<" S(z_{i}) : " << shift << endl;

			vec_dummy[i] = vec_dummy[i] - (N_zi)/(c1 - N_zi*shift);
		}

		cout << "\nz_{i} new: " << endl;
		printComplexVector(vec_dummy);

		complex<double> epsilon(0.0,0.0);		
		for (int i = 0; i < n; ++i)
		{
			epsilon += vec_dummy[i]-vec_check[i];
		}
		double diffnorm = moduluscomplex(epsilon);
		if( diffnorm < 1e-12)
		{
			cout << "\nRoots found at iteration: " << k << endl;
			k = N-1;		
		}
		
	}

	for(int i = 0; i < n;++i)
	{
		if(abs(imag(vec_dummy[i])) < 1e-12 && abs(real(vec_dummy[i])) > 1e-12)
		{
			complex<double> root(real(vec_dummy[i]), 0.0);
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) < 1e-12 && abs(imag(vec_dummy[i])) > 1e-12)
		{
			complex<double> root(0.0,imag(vec_dummy[i]));
			vec_update.push_back(root);
		}
		if(abs(real(vec_dummy[i])) > 1e-12 && abs(imag(vec_dummy[i])) > 1e-12)
		{
			complex<double> root(real(vec_dummy[i]),imag(vec_dummy[i]));
			vec_update.push_back(root);
		}

	}

	cout << "\n************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\nz_{i} final: " << endl;
	printComplexVector(vec_update);
		
}

Symbolic eulermethod(const Symbolic &f, const Symbolic &y, const Symbolic &x, const Symbolic &y0, const Symbolic &x0, const Symbolic &x1, double h)
{
 	Symbolic tangent, t("t");
	double t_now, y_now, f_now;
	double N1 = (x1-x0)/h;
	int N = N1;
	y_now = y0;
	t_now = x0;

	f_now = f[x==t_now, y==y_now] ;
	tangent = y_now + f_now*(t-t_now);

	cout << "\nf(x) = " << f <<endl;
	
	cout << endl;
	cout << setw(6) << "t" << "\t\t" << "Euler approximation y_{i}" << "\t\t\t" << "Tangent line"  << "\n";
	cout << setprecision(6) << setw(6) << x0 << "\t\t" << y_now << "\t\t\t\t\t" << tangent << "\n";	
	for (int i = x0; i < N;  i++)
	{
		t_now = t_now + h;
		y_now = y_now + f_now*h;
		f_now = f[x==t_now, y==y_now] ;
		tangent = y_now + f_now*(t-t_now);

		cout << setprecision(6) << setw(6) << t_now << "\t\t" << y_now << "\t\t\t\t\t" << tangent << "\n";
		
	}
	cout << "\nsolution = "<< endl;
	return y_now;
}

double numericaldifferentiation(const Symbolic &f, const Symbolic &x, double x0, double h)
{
	Symbolic df_symbolic = (f[x==x0+h]-f[x==x0])/(h);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double numericaldifferentiation3pointoneside(const Symbolic &f, const Symbolic &x, double x0, double h)
{
	Symbolic df_symbolic = (-3*f[x==x0] + 4*f[x==x0+h] - f[x==x0 + 2*h])/(2*h);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double numericaldifferentiation3pointbothsides(const Symbolic &f, const Symbolic &x, double x0, double h)
{
	Symbolic df_symbolic = (f[x==x0+h] - f[x==x0 - h])/(2*h);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double numericaldifferentiation5pointoneside(const Symbolic &f, const Symbolic &x, double x0, double h)
{
	Symbolic df_symbolic = (f[x==x0 - 2*h] - 8*f[x==x0-h] + 8*f[x==x0+h] - f[x==x0 + 2*h])/(12*h);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double numericaldifferentiation5pointbothsides(const Symbolic &f, const Symbolic &x, double x0, double h)
{
	Symbolic df_symbolic = (-25*f[x==x0] + 48*f[x==x0+h] -36*f[x==x0+2*h] + 16*f[x==x0 + 3*h] - 3*f[x==x0+4*h])/(12*h);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double trapezoidalrule(const Symbolic &f, const Symbolic &x, double a, double b)
{
	double h = (b-a)/1;
	Symbolic df_symbolic = (f[x==a] + f[x==b])*h/(2);

	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double simpsonsrule(const Symbolic &f, const Symbolic &x, double a, double b)
{
	double h = (b-a)/2;
	Symbolic df_symbolic = (f[x==a] + 4*f[x==a+h] + f[x==b])*h/(3);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double simpsonsrule38(const Symbolic &f, const Symbolic &x, double a, double b)
{
	double h = (b-a)/3;
	Symbolic df_symbolic = (f[x==a] + 3*f[x==a+h] + 3*f[x==a+2*h] + f[x==b])*3*h/(8);
	
	double df_numeric = df_symbolic[SymbolicConstant::e == exp(1)] ;
	return df_numeric;
}

double richardsonextrapolation(const Symbolic &f, const Symbolic &x, double x0, double h, int level, int order_of_error)
{
	std::vector<double> approx(level);
	std::vector<double> current_h_value(level);
	std::vector<vector<double>> squareMatrix(level, vector<double>(level));

	for (int i = 0; i< level ; i++)
	{
		current_h_value[i] = h / std::pow(2, i);
		approx[i] = (1/(2*current_h_value[i]))*(evalf(f[x==x0+current_h_value[i]],x,1) - evalf(f[x==x0-current_h_value[i]],x,1));
		squareMatrix[i][0] = approx[i] ;
	}
	int j1 = 0;
	for (int k = 1; k < level; ++k) 
	{
		for (int i = 0; i < level - k; ++i) 
		{
			// Burden Faires Numerical Analysis book subsection 4.2 centered difference formula with O(h^(order-of_error))
			approx[i] = approx[i + 1] + (approx[i + 1] - approx[i]) / (std::pow(order_of_error, k) - 1);
			
			squareMatrix[i+1+j1][k] = approx[i] ;			
		}
		j1 = j1+1;
	}
	cout << "Richardson's Extrapolation table:\n" << endl;
	// Display the Richardson's extrapolation in matrix display
	for(int i = 0; i<level; i++)
	{
		cout << setw(23);
		for(int j=0; j<level; j++)
		{
			cout << setprecision(10) << squareMatrix[i][j] << setw(23);
		}
		cout << endl;
	}
	cout << endl;
	cout << "\nBest approximation :" << endl;
	return squareMatrix[level-1][level-1] ;
}

void gradientdescent(const Symbolic &f, const Symbolic &x, double x0, double alpha, int epochs)
{
	double x_iter = x0;           // Starting point (initial guess)
	Symbolic der_f = df(f,x);
	for (int i = 0; i < epochs; ++i) 
	{
	double gradient = der_f[x==x_iter];
        
	// Gradient Descent Formula
	x_iter = x_iter - (alpha * gradient);

	cout << "Iteration " << i + 1 << ": \tx = " << x_iter 
		<< ", \tf(x) = " << f[x==x_iter] << endl;
	}
}

Symbolic directionfield(const Symbolic &f, const Symbolic &tf, const Symbolic &yf, double x_min, double x_max,  double y_min, double y_max,  double step_size, double k)
{
	// Open a file to write data for Gnuplot
	ofstream dataFile("direction_field.dat");
	
	double dx, dy, magnitude;
	for (double y = y_min; y <= y_max; y += step_size) 
	{
		for (double x = x_min; x <= x_max; x += step_size) 
			{
				if (std::abs(y) > 0.001) 
				{ // Avoid division by zero for this example
				double slope = f[tf==x, yf==y];
				// Calculate components of a unit vector in the direction of the slope
				dx = 1.0;
				dy = slope;
				magnitude = sqrt(dx * dx + dy * dy);
				dx = k * dx / magnitude;
				dy = k * dy / magnitude;
				
				// Write starting point (x,y) and vector components (dx, dy) to file
				dataFile << x << " \t \t  " << y << " \t  \t " << dx << " \t  \t " << dy << endl;
				}
			}	
	}

	dataFile.close();
	return 0;
}

void choleskyDecomposition(vector<vector<double>> matrix)
{
	int n = matrix.size();

	// to store the lower triangular matrix
	vector<vector<double>> lower(n, vector<double>(n, 0));

	// Decomposing a matrix into Lower Triangular
	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j <= i; j++) 
		{
			double sum = 0;
			// summation for diagonals
			if (j == i) 
			{
				for (int k = 0; k < j; k++)
				{
					sum += pow(lower[j][k], 2);
					
				}
				lower[j][j] = sqrt(matrix[j][j] - sum);
			} 
			else 
			{
		        // Evaluating L(i, j) using L(j, j)
				for (int k = 0; k < j; k++)
				{
					sum += (lower[i][k] * lower[j][k]);
					
				}
				lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
			}
		}
	}
	
	cout << "\nA = " << endl;
	// Displaying Lower Triangular Matrix
	for (int i = 0; i < n; i++) 
	{
	// Lower Triangular
		for (int j = 0; j < n; j++)
		{
			cout << setw(10) << lower[i][j] << setw(10);
		}
		cout << endl;
	}

	cout<<endl;

	cout << "\nA^{T} = " << endl;

	// Displaying Transpose of Lower Triangular Matrix
	for (int i = 0; i < n; i++) 
	{        
		// Lower Triangular
		for (int j = 0; j < n; j++)
		{
			cout <<  setw(10) << lower[j][i] <<  setw(10);
		}
		cout << endl;
	}
}

// Function to perform LU decomposition
void LUDecomposition(vector<vector<double>> &A, vector<vector<double>> &L, vector<vector<double>> &U) 
{
	int n = A.size();	
	
	// Initialize L with ones on the diagonal and zeros above
	// Initialize U with zeros below the diagonal
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			L[i][j] = (i == j) ? 1.0 : 0.0; // Ones on diagonal for L
			U[i][j] = 0.0;
		}
	}
	
	// Doolittle's algorithm
	for (int i = 0; i < n; ++i) 
	{
		// Calculate U elements
		for (int j = i; j < n; ++j) 
		{
			double sum = 0.0;
			for (int k = 0; k < i; ++k) 
			{
				sum += L[i][k] * U[k][j];
			}
			U[i][j] = A[i][j] - sum;
		}
        

		// Calculate L elements (below the diagonal)
		for (int j = i + 1; j < n; ++j) 
		{
			double sum = 0.0;
			for (int k = 0; k < i; ++k) 
			{
				sum += L[j][k] * U[k][i];
			}
			L[j][i] = (A[j][i] - sum) / U[i][i]; // U[i][i] cannot be zero
		}
	}
}

/*

	Stochastic Optimizers

*/


// --- Helper function to keep candidates strictly within bounds ---
namespace 
{
	vector<double> clamp_to_bounds(const vector<double>& x, const vector<double>& low, const vector<double>& high) 
	{
		vector<double> clamped = x;
		for (size_t i = 0; i < x.size(); ++i) 
		{
			if (clamped[i] < low[i]) clamped[i] = low[i];
			if (clamped[i] > high[i]) clamped[i] = high[i];
		}
		return clamped;
	}
}

// ============================================================================
// 1. N-DIMENSIONAL DIRECT MONTE CARLO SOLVER
// ============================================================================
DMCMetrics run_direct_monte_carlo(int total_samples, const vector<double>& lower_bounds, const vector<double>& upper_bounds, std::function<double(const std::vector<double>&)> obj_func) 
{
	size_t dim = lower_bounds.size();
	random_device rd; mt19937 gen(rd());
    
	vector<uniform_real_distribution<double>> distributions;
	distributions.reserve(dim);
	for (size_t d = 0; d < dim; ++d) 
	{
		distributions.emplace_back(lower_bounds[d], upper_bounds[d]);
	}

	DMCMetrics metrics;
	metrics.total_samples = total_samples;
	metrics.best_energy = numeric_limits<double>::infinity();
	metrics.best_x.assign(dim, 0.0);

	for (int i = 0; i < total_samples; i++) 
	{
		vector<double> current_x(dim);
		for (size_t d = 0; d < dim; ++d) 
		{
			current_x[d] = distributions[d](gen);
 		}
        
	double current_e = obj_func(current_x);
		if (current_e < metrics.best_energy) 
		{
			metrics.best_energy = current_e;
			metrics.best_x = current_x;
		}
	}
	return metrics;
}

// ============================================================================
// 2. N-DIMENSIONAL STOCHASTIC TUNNELING (STUN) SOLVER
// ============================================================================
STUNMetrics run_stochastic_tunneling(int max_steps, double initial_step_size, double gamma, double temperature, const vector<double>& initial_base, std::function<double(const std::vector<double>&)> obj_func) 
{
	size_t dim = initial_base.size();
	random_device rd; mt19937 gen(rd());
	uniform_real_distribution<double> rand_u(0.0, 1.0);
	normal_distribution<double> rand_n(0.0, initial_step_size); 

	vector<double> current_x = initial_base;
	double current_e = obj_func(current_x);
	double f_0 = current_e;

	STUNMetrics metrics;
	metrics.best_x = current_x;
	metrics.best_energy = current_e;
	metrics.total_steps = max_steps;
	metrics.accepted_moves = 0;

	for (int step = 0; step < max_steps; step++) 
	{
		vector<double> proposal(dim);
		for (size_t d = 0; d < dim; ++d) 
		{
			proposal[d] = current_x[d] + rand_n(gen);
		}
        
		double proposal_e = obj_func(proposal);

		double stun_current = (current_e >= f_0) ? (1.0 - exp(-gamma * (current_e - f_0))) : (current_e - f_0);
		double stun_proposal = (proposal_e >= f_0) ? (1.0 - exp(-gamma * (proposal_e - f_0))) : (proposal_e - f_0);
		double delta_stun = stun_proposal - stun_current;

		if (delta_stun <= 0.0 || rand_u(gen) < exp(-delta_stun / temperature)) 
		{
			current_x = proposal;
			current_e = proposal_e;
			metrics.accepted_moves++;

			if (current_e < f_0) 
			{
				f_0 = current_e;
				metrics.best_energy = current_e;
				metrics.best_x = current_x;
			}
		}
	}
	return metrics;
}

// ============================================================================
// 3. N-DIMENSIONAL PARALLEL TEMPERING SOLVER
// ============================================================================
struct PrivateReplica 
{
	vector<double> x;  
	double energy;     
	double temperature;
};

PTMetrics run_parallel_tempering(int total_cycles, int local_steps, const vector<double>& initial_base, const vector<double>& temperatures, std::function<double(const std::vector<double>&)> obj_func) 
{
	int num_replicas = temperatures.size();
	size_t dim = initial_base.size();
	random_device rd; mt19937 gen(rd());
	uniform_real_distribution<double> rand_u(0.0, 1.0);
	normal_distribution<double> rand_n(0.0, 0.1); 

	vector<PrivateReplica> replicas(num_replicas);
	for (int i = 0; i < num_replicas; i++) 
	{
		replicas[i].x = initial_base;
		for (size_t d = 0; d < dim; ++d) 
		{
			replicas[i].x[d] += i * 0.02;
		}
		replicas[i].energy = obj_func(replicas[i].x);
		replicas[i].temperature = temperatures[i];
	}

	PTMetrics metrics;
	metrics.best_x = replicas[0].x;
	metrics.best_energy = replicas[0].energy;
	metrics.swap_attempts = 0;
	metrics.swap_successes = 0;

	for (int cycle = 0; cycle < total_cycles; cycle++) 
	{
		for (int i = 0; i < num_replicas; i++) 
		{
			for (int step = 0; step < local_steps; step++) 
			{
				vector<double> proposal(dim);
				for (size_t d = 0; d < dim; ++d) 
				{
					proposal[d] = replicas[i].x[d] + rand_n(gen);
				}
					
				double proposal_energy = obj_func(proposal);
				double delta_e = proposal_energy - replicas[i].energy;
					
				if (delta_e <= 0 || rand_u(gen) < exp(-delta_e / replicas[i].temperature)) 
				{
					replicas[i].x = proposal;
					replicas[i].energy = proposal_energy;
					    
					if (replicas[i].energy < metrics.best_energy) 
					{
						metrics.best_energy = replicas[i].energy;
						metrics.best_x = proposal;
					}
				}
			}
		}
        
		for (int i = 0; i < num_replicas - 1; i++) 
		{
		 	int r1 = i; int r2 = i + 1;
			double beta1 = 1.0 / replicas[r1].temperature;
			double beta2 = 1.0 / replicas[r2].temperature;
			double swap_delta = (beta1 - beta2) * (replicas[r1].energy - replicas[r2].energy);
			    
			metrics.swap_attempts++;
			if (swap_delta >= 0 || rand_u(gen) < exp(swap_delta)) 
			{
				metrics.swap_successes++;
				swap(replicas[r1].x, replicas[r2].x);
				swap(replicas[r1].energy, replicas[r2].energy);
			}
		}
	}
	return metrics;
}

// ============================================================================
// UNIFIED CROSS-PLATFORM BENCHMARK SUITE IMPLEMENTATION
// ============================================================================
vector<BenchmarkRecord> run_unified_stochastic_benchmark(
	int dimensions, const vector<double>& lower_bounds, const vector<double>& upper_bounds,
	const vector<double>& initial_base, int budget_evaluations, std::function<double(const std::vector<double>&)> obj_func) 
{
	vector<BenchmarkRecord> records;
	chrono::high_resolution_clock::time_point start, end;

	// Apply strict axis boundary safety checking for optimization neighborhood generations
	auto bounded_obj_func = [&](const vector<double>& x) 
	{
		return obj_func(clamp_to_bounds(x, lower_bounds, upper_bounds));
	};

	// 1. Benchmark Direct Monte Carlo
	start = chrono::high_resolution_clock::now();
	DMCMetrics dmc_res = run_direct_monte_carlo(budget_evaluations, lower_bounds, upper_bounds, bounded_obj_func);
	end = chrono::high_resolution_clock::now();
	records.push_back({
		"Direct Monte Carlo", dmc_res.best_energy, dmc_res.best_x,
		chrono::duration<double, milli>(end - start).count()
	});

	// 2. Benchmark Stochastic Tunneling (STUN)
	start = chrono::high_resolution_clock::now();
	STUNMetrics stun_res = run_stochastic_tunneling(budget_evaluations, 0.1, 2.0, 0.05, initial_base, bounded_obj_func);
	end = chrono::high_resolution_clock::now();
	records.push_back({
		"Stochastic Tunneling", stun_res.best_energy, stun_res.best_x,
		chrono::duration<double, milli>(end - start).count()
	});

	// 3. Benchmark Parallel Tempering
	vector<double> pt_temps = {0.05, 0.22, 0.93, 4.0};
	int num_replicas = pt_temps.size();
	int local_steps = 10;
	int pt_cycles = budget_evaluations / (num_replicas * local_steps);
	if (pt_cycles < 1) pt_cycles = 1;

	start = chrono::high_resolution_clock::now();
	PTMetrics pt_res = run_parallel_tempering(pt_cycles, local_steps, initial_base, pt_temps, bounded_obj_func);
	end = chrono::high_resolution_clock::now();
	records.push_back({
		"Parallel Tempering", pt_res.best_energy, pt_res.best_x,
		chrono::duration<double, milli>(end - start).count()
	});

	return records;
}


/*

	Iterative Techniques in Matrix Algebra

*/

// Function to solve Ax = b using Jacobi Iteration
void JacobiMethod(const vector<vector<double>>& A, const vector<double>& b, int maxIterations, double tolerance) 
{
	int n = b.size();
	vector<double> x(n, 0.0);      // Current iteration values (initialized to 0)
	vector<double> x_old(n, 0.0);  // Previous iteration values

	cout << std::fixed << std::setprecision(6);
	cout << "Starting Jacobi Iteration...\n\n";

	for (int k = 1; k <= maxIterations; ++k) 
	{
		// Save current results to x_old before updating, simultaneous updates
		x_old = x;

		for (int i = 0; i < n; ++i) 
		{
			double sum = 0.0;
			for (int j = 0; j < n; ++j) 
			{
				if (i != j) 
				{
					sum += A[i][j] * x_old[j];
				}
			}
			// Apply Jacobi formula
			x[i] = (b[i] - sum) / A[i][i];
		}

		// Check for convergence (L2 norm of the difference / standard Euclidean difference)
		double diffNorm = 0.0;
		for (int i = 0; i < n; ++i) 
		{
			diffNorm += std::pow(x[i] - x_old[i], 2);
		}
		diffNorm = std::sqrt(diffNorm);

		// Print progress
  		cout << "Iteration " << k << ": ";
		for (int i = 0; i < n; ++i) 
		{
			cout << "x[" << i << "]=" << x[i] << "  \t ";
		}
		cout << "(Error: " << diffNorm << ")\n";

		// Stop if the solution has converged
		if (diffNorm < tolerance) 
		{
			cout << "\nConverged in " << k << " iterations.\n";
			return;
		}
	}
	cout << "\nReached maximum iterations without full convergence.\n";
}

// Function to perform the Gauss-Seidel Method
bool GaussSeidel(
	const vector<vector<double>>& A, 
	const vector<double>& b, 
	vector<double>& x, 
	double tolerance , 
	int maxIterations 
) 
{
	int n = b.size();
    
	cout << std::fixed << std::setprecision(6);
	cout << "Starting Gauss-Seidel Method \n";

	// Check if diagonal elements are zero
	for (int i = 0; i < n; ++i) 
	{
		if (std::abs(A[i][i]) < 1e-12) 
		{
			std::cerr << "Error: Diagonal element A[" << i << "][" << i << "] is close to zero." << std::endl;
			return false;
		}
	}

	cout << "Iterative Steps:\n";
    
	for (int iter = 1; iter <= maxIterations; ++iter) 
	{
		bool converged = true;
        	double diffNorm = 0.0;
		for (int i = 0; i < n; ++i) 
		{
			double sum = b[i];
            
			for (int j = 0; j < n; ++j) 
			{
				if (i != j) 
				{
					sum -= A[i][j] * x[j]; // Uses updated x[j] if j < i, and old x[j] if j > i
				}
			}
            
			double newValue = sum / A[i][i];

			diffNorm += std::pow(newValue - x[i], 2);

 			// Check convergence criteria based on absolute change
			// alternative: std::abs(newValue - x[i]) > tolerance
			if (std::abs(newValue - x[i]) > tolerance) 
			{
				converged = false;
			}
            
			x[i] = newValue; // Instant update
		}
		
		diffNorm = std::sqrt(diffNorm);

		// Print current iteration values
		cout << "Iteration " << iter << ": ";
		for (int i = 0; i < n; ++i) 
		{
			cout << "x[" << i << "]=" << x[i] << "  \t ";
		}
		cout << "(Error: " << diffNorm << ")\n"; // alternative use max_error

		if (converged) 
		{
			cout << "\nConverged successfully in " << iter << " iterations.\n";
			return true;
		}
	}

	cout << "\nWarning: Reached maximum iterations without complete convergence.\n";
	return false;
}

// Function to solve Ax = b using SOR method
bool SORIterativeMethod(
	const vector<vector<double>>& A, // Coefficient matrix
	const vector<double>& b,              // Right-hand side vector
	vector<double>& x,                    // Initial guess / output solution
	double omega,                              // Relaxation factor (1 < omega < 2)
	double tolerance,                          // Convergence threshold
	int maxIterations                          // Iteration safety cap
) 
{
	int n = A.size();
	vector<double> x_old = x;

	cout << std::fixed << std::setprecision(6);
	cout << "Starting SOR Method with omega = " << omega << "\n\n";

	for (int iter = 1; iter <= maxIterations; ++iter) 
	{
		double max_error = 0.0;

		for (int i = 0; i < n; ++i) 
		{
			double sum = 0.0;

			// Calculate the summation part of the equation
			for (int j = 0; j < n; ++j) 
			{
				if (j != i) 
				{
					sum += A[i][j] * x[j]; // Uses newly updated values for j < i automatically
				}
			}

		// Gauss-Seidel intermediate step
		double gs_value = (b[i] - sum) / A[i][i];

		// Apply the SOR relaxation formula
		x[i] = (1.0 - omega) * x_old[i] + omega * gs_value;

		// Track the maximum absolute change for convergence criteria
		max_error = std::max(max_error, std::abs(x[i] - x_old[i]));

		}
		// Check for convergence (L2 norm of the difference / standard Euclidean difference)
		double diffNorm = 0.0;
		for (int i = 0; i < n; ++i) 
		{
			diffNorm += std::pow(x[i] - x_old[i], 2);
		}
		diffNorm = std::sqrt(diffNorm);

		
		// Print progress (optional)
  		cout << "Iteration " << iter << ": ";
		for (int i = 0; i < n; ++i) 
		{
			cout << "x[" << i << "]=" << x[i] << "  \t ";
		}
		cout << "(Error: " << diffNorm << ")\n"; // alternative use max_error

		// Check for convergence
		if (diffNorm < tolerance) 
		{
			cout << "\nConverged in " << iter << " iterations.\n";
			return true;
		}

		// Update old values for the next iteration
		x_old = x;
	}

	cout << "\nReached maximum iterations without strict convergence.\n";
	return false;
}

// Functions to solve Ax = b using Conjugate Gradient with CRS format
// Function to convert a dense 2D vector to CRS format
CRSMatrix denseToCRS(const vector<vector<double>>& dense) 
{
	CRSMatrix crs;
	crs.num_rows = dense.size();
	crs.num_cols = crs.num_rows > 0 ? dense[0].size() : 0;
    
	crs.row_ptr.push_back(0); // First element is always 0

	for (int i = 0; i < crs.num_rows; ++i) 
	{
		for (int j = 0; j < crs.num_cols; ++j) 
		{
			if (dense[i][j] != 0.0) 
			{
				crs.values.push_back(dense[i][j]);
				crs.col_indices.push_back(j);
			}
		}
		crs.row_ptr.push_back(crs.values.size());
	}
	return crs;
}

// Helper function: Sparse Matrix-Vector Multiplication (y = A * x)
vector<double> spmv(const CRSMatrix& A, const std::vector<double>& x) 
{
	vector<double> y(A.num_rows, 0.0);
	for (int i = 0; i < A.num_rows; ++i) 
	{
		double sum = 0.0;
		int row_start = A.row_ptr[i];
		int row_end = A.row_ptr[i + 1];
		for (int k = row_start; k < row_end; ++k) 
		{
			sum += A.values[k] * x[A.col_indices[k]];
		}
		y[i] = sum;
	}
	return y;
}

// Helper function: Vector dot product (u . v)
double dotProduct(const vector<double>& u, const vector<double>& v) 
{
	double dot = 0.0;
	for (size_t i = 0; i < u.size(); ++i) 
	{
		dot += u[i] * v[i];
	}
	return dot;
}

// Extracts the diagonal elements for the Jacobi preconditioner
vector<double> get_Jacobi_preconditioner(const CRSMatrix& A) 
{
	vector<double> inv_M(A.num_rows, 0.0);
	for (int i = 0; i < A.num_rows; ++i) 
	{
		for (int j = A.row_ptr[i]; j < A.row_ptr[i + 1]; ++j) 
		{
			if (A.col_indices[j] == i) 
			{
				if (std::abs(A.values[j]) < 1e-12) 
				{
					throw std::runtime_error("Zero or near-zero diagonal element encountered.");
				}
				inv_M[i] = 1.0 / A.values[j]; // Store inverse of diagonal
				break;
			}
		}
	}
	return inv_M;
}

// Conjugate Gradient Solver for CRS Matrix
vector<double> ConjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iterations) 
{
	cout << "Starting Conjugate Gradient Method.\n\n";

	int n = A.num_rows;
	vector<double> x(n, 0.0); // Initial guess x_0 = 0
    
	// r = b - A * x
	vector<double> Ax = spmv(A, x);
	vector<double> r(n);
	for (int i = 0; i < n; ++i) 
	{
		r[i] = b[i] - Ax[i];
	}
    
	// p = r
	vector<double> p = r;
	double rsold = dotProduct(r, r);
    
	if (std::sqrt(rsold) < tolerance) 
	{
		return x; // Initial guess is already close enough
	}

	for (int iter = 0; iter < max_iterations; ++iter) 
	{
		vector<double> Ap = spmv(A, p);
		double alpha = rsold / dotProduct(p, Ap);
        
		// Update x and r
		for (int i = 0; i < n; ++i) 
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}
        
		double rsnew = dotProduct(r, r);
        	
		// Print progress (optional)
  		cout << "Iteration " << iter + 1 << ": ";
		for (int i = 0; i < n; ++i) 
		{
			cout << "x[" << i << "]=" << x[i] << "  \t ";
		}
		cout << "( || r || : " << rsnew << ")\n"; 

		// Check convergence
		if (std::sqrt(rsnew) < tolerance) 
		{
			cout << "Converged in " << iter + 1 << " iterations.\n";
			return x;
		}
        
		// Update direction p
		double beta = rsnew / rsold;
		for (int i = 0; i < n; ++i) 
		{
			p[i] = r[i] + beta * p[i];
		}
		rsold = rsnew;
	}
    
	cout << "Warning: Maximum iterations reached without full convergence.\n";
	return x;
}

// Preconditioned Conjugate Gradient (PCG) Solver
vector<double> PreconditionedConjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iter) 
{
	int n = A.num_rows;
	vector<double> x(n, 0.0); // Initial guess x_0 = 0
	vector<double> r(n);
	vector<double> z(n);
	vector<double> p(n);
	vector<double> Ap;

	//cout << std::fixed << std::setprecision(8);
	cout << "Starting Preconditioned Conjugate Gradient Method.\n\n";

	// Get the inverse diagonal for Jacobi preconditioning
	vector<double> inv_M = get_Jacobi_preconditioner(A);

	// Compute initial residual: r = b - A * x
	Ap = spmv(A, x);
	for (int i = 0; i < n; ++i) 
	{
		r[i] = b[i] - Ap[i];
	}

	// Check if initial guess is already accurate enough
	double r_norm = std::sqrt(dotProduct(r, r));
	if (r_norm < tolerance) 
	{
		cout << "Converged in " << 0 << " iterations.\n";
	}

	// Apply Preconditioner: z = M^-1 * r
	for (int i = 0; i < n; ++i) 
	{
		z[i] = inv_M[i] * r[i];
	}

	// Initial search direction: p = z
	p = z;

	double rho = dotProduct(r, z);

	for (int iter = 1; iter <= max_iter; ++iter) 
	{
		// Ap = A * p
		Ap = spmv(A, p);

		// Compute step length alpha = (r . z) / (p . Ap)
		double p_Ap = dotProduct(p, Ap);
		double alpha = rho / p_Ap;

		// Update solution and residual
		for (int i = 0; i < n; ++i) 
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}
		
		// Check convergence
		r_norm = std::sqrt(dotProduct(r, r));

		// Print progress (optional)
  		cout << "Iteration " << iter << ": ";
		for (int i = 0; i < n; ++i) 
		{
			cout << "x[" << i << "]=" << x[i] << "  \t ";
		}
		cout << "( || r || : " << r_norm << ")\n"; 

		
        	if (r_norm < tolerance) 
		{
			cout << "Converged in " << iter  << " iterations.\n";
			return x;
		}

		// Apply Preconditioner to new residual: z = M^-1 * r
		for (int i = 0; i < n; ++i) 
		{
			z[i] = inv_M[i] * r[i];
		}

		// Compute beta = (r_new . z_new) / (r_old . z_old)
		double rho_new = dotProduct(r, z);
		double beta = rho_new / rho;
		rho = rho_new;

		// Update search direction: p = z + beta * p
		for (int i = 0; i < n; ++i) 
		{
			p[i] = z[i] + beta * p[i];
		}
	}
	cout << "Warning: Maximum iterations reached without full convergence.\n";  // Did not converge within max_iter
	return x;
}


// Newton's Method to solve nonlinear system
using SystemFunction = function<vector<DualClass>(const vector<DualClass>&)>;
vector<complex<double>> NonlinearNewtonSolver::solve(SystemFunction func, const vector<complex<double>>& initial_guess) 
{	
	size_t n = initial_guess.size();
	vector<complex<double>> x = initial_guess;
	cout << std::fixed << std::setprecision(10);

	cout << "\n" ;
	cout << "========================================================\n";
	cout << "  Newton' Method Multivariable Complex Root Solver \n";
	cout << "========================================================\n";


	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// Seed independent variables into DualClass vector
		vector<DualClass> x_dual(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_dual[i] = DualClass(x[i], i, n);
		}

		// Evaluate system functions & automatically compute gradients (Jacobian matrix)
		vector<DualClass> f_dual = func(x_dual);
            
		vector<complex<double>> F(n);
		vector<vector<complex<double>>> J(n, vector<complex<double>>(n));

		double residual_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			F[i] = f_dual[i].val;
			residual_norm += norm(F[i]); // magnitude squared
			for (size_t j = 0; j < n; ++j) 
			{
				J[i][j] = f_dual[i].der[j];
			}
		}
		residual_norm = sqrt(residual_norm);

		if (residual_norm < tol) 
		{
			cout << "\n[Success] Converged in " << iter << " iterations." << endl;
			return x;
		}

		// Solve linear system: J * delta_x = F
		vector<complex<double>> delta_x = solveLinearSystem(J, F);

		// Update formula: x_{k+1} = x_k - delta_x
		// Update the state variables
		double delta_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			x[i] -= delta_x[i];
			delta_norm += std::norm(delta_x[i]);
		}
		delta_norm = std::sqrt(delta_norm);

		cout << "Iter " << iter << " | Residual Norm: " << residual_norm << " | Update Norm: " << delta_norm << "\n";
	}
	cout << "\n[Warning] Reached max iterations without reaching target tolerance." << endl;
	return x;
}

void NonlinearNewtonSolver::solve_printSolution(SystemFunction func, const vector<complex<double>>& initial_guess) 
{	
	size_t n = initial_guess.size();
	roots = initial_guess;
	cout << std::fixed << std::setprecision(10);

	cout << "\n" ;
	cout << "========================================================\n";
	cout << "  Newton' Method Multivariable Complex Root Solver \n";
	cout << "========================================================\n";

	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// Seed independent variables into DualClass vector
		vector<DualClass> x_dual(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_dual[i] = DualClass(roots[i], i, n);
		}

		// Evaluate system functions & automatically compute gradients (Jacobian matrix)
		vector<DualClass> f_dual = func(x_dual);
            
		vector<complex<double>> F(n);
		vector<vector<complex<double>>> J(n, vector<complex<double>>(n));

		double residual_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			F[i] = f_dual[i].val;
			residual_norm += norm(F[i]); // magnitude squared
			for (size_t j = 0; j < n; ++j) 
			{
				J[i][j] = f_dual[i].der[j];
			}
		}
		residual_norm = sqrt(residual_norm);

		if (residual_norm < tol) 
		{
			cout << "\n[Success] Converged in " << iter << " iterations." << endl;
			break;
		}

		// Solve linear system: J * delta_x = F
		vector<complex<double>> delta_x = solveLinearSystem(J, F);

		// Update formula: x_{k+1} = x_k - delta_x
		// Update the state variables
		double delta_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			roots[i] -= delta_x[i];
			delta_norm += std::norm(delta_x[i]);
		}
		delta_norm = std::sqrt(delta_norm);

		cout << "Iter " << iter << " | Residual Norm: " << residual_norm << " | Update Norm: " << delta_norm << "\n";
		if (iter == max_iter-1) 
		{
			cout << "\n[Warning] Reached max iterations without reaching target tolerance." << endl;
		}
	}
	
	cout << "\nComputed Roots:" << endl;
	cout << "---------------------------------------------------------\n";
	for (size_t i = 0; i < roots.size(); ++i) 
	{
		cout << "x" << i + 1 << " = " << roots[i].real() ;
		if (roots[i].imag() != 0.0) 
		{
			cout << (roots[i].imag() >= 0 ? " + " : " - ") << abs(roots[i].imag()) << "i" << endl;
		}
		cout << "\n";
	}
	cout << "---------------------------------------------------------\n";
}

vector<complex<double>> NonlinearNewtonSolver::solveLinearSystem(vector<vector<complex<double>>> J, vector<complex<double>> F) 
{ // Gaussian Elimination with Partial Pivoting optimized for Complex Matrices
	size_t n = J.size();
	for (size_t i = 0; i < n; ++i) 
	{
		size_t max_row = i;
		for (size_t k = i + 1; k < n; ++k) 
		{
			if (abs(J[k][i]) > abs(J[max_row][i])) 
			{
				max_row = k;
			}
 		}
		swap(J[i], J[max_row]);
		swap(F[i], F[max_row]);

		if (abs(J[i][i]) < 1e-15) 
		{
			cerr << "Error: Matrix is singular or near-singular!" << endl;
			return vector<complex<double>>(n, 0.0);
		}

		for (size_t k = i + 1; k < n; ++k) 
		{
			complex<double> factor = J[k][i] / J[i][i];
			for (size_t j = i; j < n; ++j) 
			{
				J[k][j] -= factor * J[i][j];
			}
			F[k] -= factor * F[i];
		}
	}

	// Back substitution
	vector<complex<double>> delta_x(n);
	for (int i = n - 1; i >= 0; --i) 
	{
		complex<double> sum = 0.0;
		for (size_t j = i + 1; j < n; ++j) 
		{
			sum += J[i][j] * delta_x[j];
		}
		delta_x[i] = (F[i] - sum) / J[i][i];
	}
	return delta_x;
}


// Quasi-Newton Method to solve nonlinear system
// Linear system solver (complex Gaussian elimination with partial pivoting)
bool QuasiNewtonSolver::solveLinearSystem(vector<vector<complex<double>>> A, vector<complex<double>> b, vector<complex<double>>& x) 
{
	size_t n = A.size();
	x.assign(n, 0.0);

	for (size_t i = 0; i < n; ++i) 
	{
		// Find pivot element
		size_t maxRow = i;
		for (size_t k = i + 1; k < n; ++k) 
		{
			if (abs(A[k][i]) > abs(A[maxRow][i])) 
			{
				maxRow = k;
			}
		}
		swap(A[i], A[maxRow]);
		swap(b[i], b[maxRow]);

		if (abs(A[i][i]) < 1e-12) 
		{
			return false; // Matrix is singular
		}

		// Elimination phase
		for (size_t k = i + 1; k < n; ++k) 
		{
			complex<double> factor = A[k][i] / A[i][i];
			for (size_t j = i; j < n; ++j) 
			{
				A[k][j] -= factor * A[i][j];
			}
			b[k] -= factor * b[i];
		}
	}

	// Back substitution phase
	for (int i = n - 1; i >= 0; --i) 
	{
		complex<double> sum = 0.0;
		for (size_t j = i + 1; j < n; ++j) 
		{
			sum += A[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / A[i][i];
	}
	return true;
}

vector<complex<double>> QuasiNewtonSolver::solve_getvector(NonlinearSystem system, const vector<complex<double>>& x0) 
{
	size_t n = x0.size();
	vector<complex<double>> x = x0;

	// Step 1: Compute Initial Function Evaluation and Exact Jacobian via Automatic Differentiation
	vector<DualClass> x_dual(n);
	for (size_t i = 0; i < n; ++i) 
	{
		vector<complex<double>> d(n, 0.0);
		d[i] = 1.0; // Seed the derivative vector
		x_dual[i] = DualClass(x[i], d);
	}

	vector<DualClass> F_dual = system(x_dual);
	vector<complex<double>> F(n);
	vector<vector<complex<double>>> J(n, vector<complex<double>>(n));

	for (size_t i = 0; i < n; ++i) 
	{
		F[i] = F_dual[i].val;
		for (size_t j = 0; j < n; ++j) 
		{
			J[i][j] = F_dual[i].der[j]; // Extract exact analytical derivatives
		}
	}

	// Iterative Quasi-Newton Updates
	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// Check for early convergence
		double f_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			f_norm += norm(F[i]);
		}
		f_norm = sqrt(f_norm);
		if (f_norm < tol) 
		{
			cout << "\nConverged successfully in " << iter << " iterations.\n";
			return x;
		}

		// Solve linear system: J * delta_x = -F
		vector<complex<double>> neg_F(n);
		for (size_t i = 0; i < n; ++i) 
		{
			neg_F[i] = -F[i];
		}

		vector<complex<double>> delta_x;
		if (!solveLinearSystem(J, neg_F, delta_x)) 
		{
			cerr << "Error: Singular Jacobian matrix met during iterations.\n";
			return x;
		}

		// Update Variable Strategy: x_new = x + delta_x
		vector<complex<double>> x_new(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_new[i] = x[i] + delta_x[i];
		}

		// Evaluate system at x_new
		for (size_t i = 0; i < n; ++i) 
		{
			x_dual[i] = DualClass(x_new[i]);
		}
		vector<DualClass> F_new_dual = system(x_dual);
		vector<complex<double>> F_new(n);
		for (size_t i = 0; i < n; ++i) 
		{
			F_new[i] = F_new_dual[i].val;
		}

		// Broyden Rank-One Update Steps
		vector<complex<double>> y(n);
		double delta_x_norm_sq = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			y[i] = F_new[i] - F[i];
			delta_x_norm_sq += norm(delta_x[i]);
		}

		// Compute rank-one vector update component: (y - J * delta_x) / ||delta_x||^2
		vector<complex<double>> J_delta_x(n, 0.0);
		for (size_t i = 0; i < n; ++i) 
		{
			for (size_t j = 0; j < n; ++j) 
			{
				J_delta_x[i] += J[i][j] * delta_x[j];
			}
		}

		for (size_t i = 0; i < n; ++i) 
		{
			complex<double> numerator = y[i] - J_delta_x[i];
			for (size_t j = 0; j < n; ++j) 
			{
				J[i][j] += (numerator * conj(delta_x[j])) / delta_x_norm_sq;
			}
		}

		// Propagate forward values
		x = x_new;
		F = F_new;
	}

	cout << "\nWarning: Reached maximum allowable iterations.\n";
	return x;
}


void QuasiNewtonSolver::solve(NonlinearSystem system, const vector<complex<double>>& x0) 
{
	size_t n = x0.size();
	roots.resize(n, 0.0);	
	roots = x0;

	cout << "\n" ;
	cout << fixed << setprecision(10);
	cout << "========================================================\n";
	cout << "  Quasi-Newton Multivariable Complex Root Solver \n";
	cout << "========================================================\n";

	// Step 1: Compute Initial Function Evaluation and Exact Jacobian via Automatic Differentiation
	vector<DualClass> x_dual(n);
	for (size_t i = 0; i < n; ++i) 
	{
		vector<complex<double>> d(n, 0.0);
		d[i] = 1.0; // Seed the derivative vector
		x_dual[i] = DualClass(roots[i], d);
	}

	vector<DualClass> F_dual = system(x_dual);
	vector<complex<double>> F(n);
	vector<vector<complex<double>>> J(n, vector<complex<double>>(n));

	for (size_t i = 0; i < n; ++i) 
	{
		F[i] = F_dual[i].val;
		for (size_t j = 0; j < n; ++j) 
		{
			J[i][j] = F_dual[i].der[j]; // Extract exact analytical derivatives
		}
	}

	// Iterative Quasi-Newton Updates
	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// Check for early convergence
		double f_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			f_norm += norm(F[i]);
		}
		f_norm = sqrt(f_norm);
		if (f_norm < tol) 
		{
			cout << "\nConverged successfully in " << iter << " iterations.\n";
			break;
		}

		// Solve linear system: J * delta_x = -F
		vector<complex<double>> neg_F(n);
		for (size_t i = 0; i < n; ++i) 
		{
			neg_F[i] = -F[i];
		}

		vector<complex<double>> delta_x;
		if (!solveLinearSystem(J, neg_F, delta_x)) 
		{
			cerr << "Error: Singular Jacobian matrix met during iterations.\n";
			
		}

		// Update Variable Strategy: x_new = x + delta_x
		vector<complex<double>> x_new(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_new[i] = roots[i] + delta_x[i];
		}

		// Evaluate system at x_new
		for (size_t i = 0; i < n; ++i) 
		{
			x_dual[i] = DualClass(x_new[i]);
		}
		vector<DualClass> F_new_dual = system(x_dual);
		vector<complex<double>> F_new(n);
		for (size_t i = 0; i < n; ++i) 
		{
			F_new[i] = F_new_dual[i].val;
		}

		// Broyden Rank-One Update Steps
		vector<complex<double>> y(n);
		double delta_x_norm_sq = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			y[i] = F_new[i] - F[i];
			delta_x_norm_sq += norm(delta_x[i]);
		}

		// Compute rank-one vector update component: (y - J * delta_x) / ||delta_x||^2
		vector<complex<double>> J_delta_x(n, 0.0);
		for (size_t i = 0; i < n; ++i) 
		{
			for (size_t j = 0; j < n; ++j) 
			{
				J_delta_x[i] += J[i][j] * delta_x[j];
			}
		}

		for (size_t i = 0; i < n; ++i) 
		{
			complex<double> numerator = y[i] - J_delta_x[i];
			for (size_t j = 0; j < n; ++j) 
			{
				J[i][j] += (numerator * conj(delta_x[j])) / delta_x_norm_sq;
			}
		}
		cout << "Iter " << iter << " | Residual Norm: " << f_norm << " | Update Norm: " << delta_x_norm_sq << "\n";
		// Propagate forward values
		roots = x_new;
		F = F_new;

		if (iter == max_iter-1) 
		{
			cout << "\nWarning: Reached maximum allowable iterations.\n";
		}
	}

}

void QuasiNewtonSolver::printSolution() const
{
	// Print calculated complex roots
	cout << setprecision(10) << fixed;
	cout << "\nComputed Roots:" << endl;
	cout << "---------------------------------------------------------\n";
	for (size_t i = 0; i < roots.size(); ++i) 
	{
		cout << "x" << i + 1 << " = " << roots[i].real() ;
		if (roots[i].imag() != 0.0) 
		{
			cout << (roots[i].imag() >= 0 ? " + " : " - ") << abs(roots[i].imag()) << "i" << endl;
		}
		cout << "\n";
	}
	cout << "---------------------------------------------------------\n";
	
}

// Steepest Descent to solve nonlinear system

vector<complex<double>> SteepestDescentSolver::solve(SystemFunction system, const vector<complex<double>>& x0) 
{
	size_t n = x0.size();
	vector<complex<double>> x = x0;

	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// 1. Evaluate system with DualClass variables to get values and Jacobian
		vector<DualClass> x_dual(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_dual[i] = DualClass(x[i], i, n);
		}

		vector<DualClass> f_dual = system(x_dual);
		
		// Extract F(x) and Jacobian J(x)
		vector<complex<double>> F(n);
		vector<vector<complex<double>>> J(n, vector<complex<double>>(n));
		double g_val = 0.0;

		for (size_t i = 0; i < n; ++i) 
		{
			F[i] = f_dual[i].val;
			g_val += std::norm(F[i]); // sum of squares of magnitudes
			for (size_t j = 0; j < n; ++j) 
			{
				J[i][j] = f_dual[i].der[j];
			}
		}
		g_val *= 0.5;

		// Check convergence
		if (std::sqrt(2.0 * g_val) < tol) 
		{
			cout << "Converged in " << iter << " iterations.\n";
			return x;
		}

		// 2. Compute the Complex Gradient: grad(g) = J^H * F
		vector<complex<double>> grad(n, 0.0);
		double grad_norm_sq = 0.0;
		for (size_t j = 0; j < n; ++j) 
		{
			for (size_t i = 0; i < n; ++i) 
			{
				grad[j] += std::conj(J[i][j]) * F[i];
			}
			grad_norm_sq += std::norm(grad[j]);
		}

		if (std::sqrt(grad_norm_sq) < 1e-12) 
		{
			cout << "Gradient vanished near a local minimum at iteration " << iter << ".\n";
			return x;
		}

		// 3. Backtracking Line Search (Armijo Rule)
		double alpha = 1.0;
		double c = 1e-4;
		double beta = 0.5;
		vector<complex<double>> x_new(n);

		while (alpha > 1e-5) 
		{
			for (size_t i = 0; i < n; ++i) 
			{
				x_new[i] = x[i] - alpha * grad[i];
			}

			// Evaluate g(x_new)
			vector<DualClass> x_new_dual(n);
			for (size_t i = 0; i < n; ++i) x_new_dual[i] = DualClass(x_new[i], i, n);
			vector<DualClass> f_new_dual = system(x_new_dual);

			double g_new = 0.0;
			for (size_t i = 0; i < n; ++i) 
			{
				g_new += std::norm(f_new_dual[i].val);
			}
			g_new *= 0.5;

			if (g_new <= g_val - c * alpha * grad_norm_sq) 
			{
				break;
			}
			alpha *= beta;
		}

		// Update variables
		x = x_new;
	} // the for loop

	cout << "Reached maximum iterations without perfect convergence.\n";
	return x;
}

void SteepestDescentSolver::solve_printSolution(SystemFunction system, const vector<complex<double>>& x0) 
{
	size_t n = x0.size();
	roots = x0;
	cout << fixed << setprecision(10);
	cout << "========================================================\n";
	cout << "  Steepest Descent Nonlinear Solver \n";
	cout << "========================================================\n";
	cout << "\nTolerance required: " << tol << endl;
	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// 1. Evaluate system with DualClass variables to get values and Jacobian
		vector<DualClass> x_dual(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_dual[i] = DualClass(roots[i], i, n);
		}

		vector<DualClass> f_dual = system(x_dual);

		// Check for early convergence
		double f_norm = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			f_norm += norm(f_dual[i].val);
		}
		f_norm = sqrt(f_norm);

		// Extract F(x) and Jacobian J(x)
		vector<complex<double>> F(n);
		vector<vector<complex<double>>> J(n, vector<complex<double>>(n));
		double g_val = 0.0;
		for (size_t i = 0; i < n; ++i) 
		{
			F[i] = f_dual[i].val;
			g_val += std::norm(F[i]); // sum of squares of magnitudes
			for (size_t j = 0; j < n; ++j) 
			{
				J[i][j] = f_dual[i].der[j];
			}
		}
		g_val *= 0.5;

		// Check convergence
		//if (std::sqrt(2.0 * g_val) < tol) 
		if (g_val < tol) 
		{
			cout << "Converged in " << iter << " iterations.\n";
			break;
		}

		// 2. Compute the Complex Gradient: grad(g) = J^H * F
		vector<complex<double>> grad(n, 0.0);
		double grad_norm_sq = 0.0;
		for (size_t j = 0; j < n; ++j) 
		{
			for (size_t i = 0; i < n; ++i) 
			{
				grad[j] += std::conj(J[i][j]) * F[i];
			}
			grad_norm_sq += std::norm(grad[j]);
		}

		if (std::sqrt(grad_norm_sq) < 1e-12) 
		{
			cout << "Gradient vanished near a local minimum at iteration " << iter << ".\n";
		}

		// 3. Backtracking Line Search (Armijo Rule)
		double alpha = 1.0;
		double c = 1e-4;
		double beta = 0.5;
		vector<complex<double>> x_new(n);

		while (alpha > 1e-5) 
		{
			for (size_t i = 0; i < n; ++i) 
			{
				x_new[i] = roots[i] - alpha * grad[i];
			}

			// Evaluate g(x_new)
			vector<DualClass> x_new_dual(n);
			for (size_t i = 0; i < n; ++i) x_new_dual[i] = DualClass(x_new[i], i, n);
			vector<DualClass> f_new_dual = system(x_new_dual);

			double g_new = 0.0;
			for (size_t i = 0; i < n; ++i) 
			{
				g_new += std::norm(f_new_dual[i].val);
			}
			g_new *= 0.5;

			if (g_new <= g_val - c * alpha * grad_norm_sq) 
			{
				break;
			}
			alpha *= beta;
		}
		
		// Update variables
		roots = x_new;
		
		cout << "Iter " << iter << " | 0.5||F(x)||^{2}: " << g_val <<  "\t | \t " << "||F(x)||: " << f_norm << endl; 

		if (iter == max_iter-1)
		{
			cout << "Reached maximum iterations without perfect convergence.\n";
		}
	} // the for loop

	// Output results
	cout << "\nComputed Roots:\n";
	cout << std::fixed << std::setprecision(10);
	for (size_t i = 0; i < roots.size(); ++i) 
	{
		cout << "x" << i + 1 << " = " << roots[i].real() 
		<< (roots[i].imag() >= 0 ? " + " : " - ") 
		<< std::abs(roots[i].imag()) << "i\n";
	}
	
	
}

// PSO to solve nonlinear system
// Converts internal real representation back to std::complex<double>
ComplexVector PSOMultiRootFinder::toComplexVector(const vector<double>& real_vec) 
{
	ComplexVector comp_vec(num_vars);
	for (int i = 0; i < num_vars; ++i) 
	{
		comp_vec[i] = Complex(real_vec[2 * i], real_vec[2 * i + 1]);
	}
	return comp_vec;
}
// Evaluates system and returns the total residual magnitude squared
double PSOMultiRootFinder::evaluateFitness(const vector<double>& real_vec, const ComplexVectorFunction& system) 
{
	ComplexVector x = toComplexVector(real_vec);
	ComplexVector residuals = system(x);
	double total_error = 0.0;
	for (const auto& res : residuals) 
	{
		total_error += norm(res); // norm(z) yields |z|^2
	}
	return total_error;
}
// Executes a single PSO run to locate a single global minimum root
ComplexVector PSOMultiRootFinder::findOneRoot(const ComplexVectorFunction& system, mt19937& rng) 
{
	int dim = 2 * num_vars; // 2 dimensions per complex number
        uniform_real_distribution<double> dist_pos(lower_bound, upper_bound);
        uniform_real_distribution<double> dist_vel((lower_bound - upper_bound) * 0.1, (upper_bound - lower_bound) * 0.1);
        uniform_real_distribution<double> dist_r(0.0, 1.0);

	vector<Particle> swarm(pop_size);
	vector<double> gbest_position(dim);
	double gbest_fitness = numeric_limits<double>::max();

	// Initialization
	for (int i = 0; i < pop_size; ++i) 
	{
		swarm[i].position.resize(dim);
		swarm[i].velocity.resize(dim);
		for (int d = 0; d < dim; ++d) 
		{
			swarm[i].position[d] = dist_pos(rng);
			swarm[i].velocity[d] = dist_vel(rng);
		}
		swarm[i].best_position = swarm[i].position;
		swarm[i].best_fitness = evaluateFitness(swarm[i].position, system);

		if (swarm[i].best_fitness < gbest_fitness) 
		{
			gbest_fitness = swarm[i].best_fitness;
			gbest_position = swarm[i].best_position;
		}
	}
	// PSO Optimization Parameters
	const double w = 0.729;  // Inertia weight
	const double c1 = 1.494; // Cognitive weight
	const double c2 = 1.494; // Social weight
	// Optimization Loop
	for (int iter = 0; iter < max_iter; ++iter) 
	{
		for (int i = 0; i < pop_size; ++i) 
		{
			for (int d = 0; d < dim; ++d) 
			{
				double r1 = dist_r(rng);
				double r2 = dist_r(rng);

				swarm[i].velocity[d] = w * swarm[i].velocity[d] +
					c1 * r1 * (swarm[i].best_position[d] - swarm[i].position[d]) +
					c2 * r2 * (gbest_position[d] - swarm[i].position[d]);
                    
				swarm[i].position[d] += swarm[i].velocity[d];
			}

			double fitness = evaluateFitness(swarm[i].position, system);
			if (fitness < swarm[i].best_fitness) 
			{
				swarm[i].best_fitness = fitness;
				swarm[i].best_position = swarm[i].position;
			}

			if (fitness < gbest_fitness) 
			{
				gbest_fitness = fitness;
				gbest_position = swarm[i].position;
			}
		}
	}
	return toComplexVector(gbest_position);
}

// Multiroot manager using multi-start check routine to identify unique roots
vector<ComplexVector> PSOMultiRootFinder::findMultipleRoots(const ComplexVectorFunction& system, int num_runs) 
{
	double fitness_threshold = 1e-6;
	double uniqueness_threshold = 1e-3;
	random_device rd;
	mt19937 rng(rd());
	vector<ComplexVector> unique_roots;

	for (int run = 0; run < num_runs; ++run) 
	{
		ComplexVector root = findOneRoot(system, rng);
		double fitness = 0.0;
		for (const auto& res : system(root)) 
		{
			fitness += norm(res);
		}

		// Keep only valid roots converging below our fitness target
		if (fitness < fitness_threshold) 
		{
			bool is_unique = true;
			for (const auto& existing : unique_roots) 
			{
				double distance = 0.0;
				for (size_t i = 0; i < root.size(); ++i) 
				{
					distance += norm(root[i] - existing[i]);
				}
				if (sqrt(distance) < uniqueness_threshold) 
				{
					is_unique = false;
					break;
				}
			}
			if (is_unique) 
			{
				unique_roots.push_back(root);
			}
		}
	}
	return unique_roots;
}

void PSOMultiRootFinder::findMultipleRoots_printSolution(const ComplexVectorFunction& system, int num_runs) 
{
	double fitness_threshold = 1e-6;
	double uniqueness_threshold = 1e-3;
	random_device rd;
	mt19937 rng(rd());

	cout << "========================================================\n";
	cout << "  PSO Multivariable Complex Root Finder  \n";
	cout << "========================================================\n";
	
	cout << "\nSwarm population: " << pop_size << endl;
	cout << "Maximum iteration: " << max_iter << endl;
	cout << "Lower bound: " << lower_bound << endl;
	cout << "Upper bound: " << upper_bound << endl;
	for (int run = 0; run < num_runs; ++run) 
	{
		ComplexVector root = findOneRoot(system, rng);
		double fitness = 0.0;
		for (const auto& res : system(root)) 
		{
			fitness += norm(res);
		}

		// Keep only valid roots converging below our fitness target
		if (fitness < fitness_threshold) 
		{
			bool is_unique = true;
			for (const auto& existing : final_roots) 
			{
				double distance = 0.0;
				for (size_t i = 0; i < root.size(); ++i) 
				{
					distance += norm(root[i] - existing[i]);
				}
				if (sqrt(distance) < uniqueness_threshold) 
				{
					is_unique = false;
					break;
				}
			}
			if (is_unique) 
			{
				final_roots.push_back(root);
			}
		}
	}

	// Output formatting
	cout << fixed << setprecision(10);
	cout << "\nFound " << final_roots.size() << " unique root(s):\n";
	for (size_t i = 0; i < final_roots.size(); ++i) 
	{
		cout << "----------------------------------------\n";
		cout << "Root #" << i + 1 << ":\n";
		for (size_t j = 0; j < final_roots[0].size(); ++j) 
		{
			cout << "x" << j + 1 << " = " << final_roots[i][j].real() ;

			cout << (final_roots[i][j].imag() >= 0 ? " + " : " - ") << abs(final_roots[i][j].imag()) << "i" << endl;
			
		}

        
		// Output verification residuals
		ComplexVector res = system(final_roots[i]);
		double f_norm=0;
		for (size_t j = 0; j < final_roots[0].size(); ++j) 
		{
			f_norm += abs(res[j]) ;
		}		
		cout << "Residual Norm = " << f_norm << "\n";
	}
	cout << "----------------------------------------\n";
	
}

// RK4 Multiroot solver

bool RK4MultiRootSolver::solveLinearSystem(vector<vector<complex<double>>>& A, vector<complex<double>>& b, vector<complex<double>>& x) 
{
	size_t n = A.size();
	x.assign(n, 0);

	for (size_t i = 0; i < n; ++i) 
	{
	  	// Pivot selection
		size_t maxRow = i;
		for (size_t k = i + 1; k < n; ++k) 
		{
			if (std::abs(A[k][i]) > std::abs(A[maxRow][i]))
			{ 
				maxRow = k;
			}
		}
		std::swap(A[i], A[maxRow]);
		std::swap(b[i], b[maxRow]);

		if (std::abs(A[i][i]) < 1e-12) 
		{
			return false; // Singular matrix
		}

		for (size_t k = i + 1; k < n; ++k) 
		{
			complex<double> c = -A[k][i] / A[i][i];
			for (size_t j = i; j < n; ++j) 
			{
				A[k][j] += c * A[i][j];
			}
			b[k] += c * b[i];
		}
	}

	// Back substitution
	for (int i = n - 1; i >= 0; --i) 
	{
		x[i] = b[i];
		for (size_t j = i + 1; j < n; ++j) 
		{
			x[i] -= A[i][j] * x[j];
		}
		x[i] /= A[i][i];
	}
	return true;
}


// Continuous Newton function derivative: dx/dt = -J^-1 * F(x)
vector<complex<double>> RK4MultiRootSolver::computeDerivatives(const vector<complex<double>>& x, SystemFunc system) 
{
        size_t n = x.size();
	vector<DualClassTemplate<complex<double>>> dual_x(n);
        
	// Initialize dual variables for automatic differentiation
	for (size_t i = 0; i < n; ++i) 
	{
		vector<complex<double>> der(n, 0.0);
		der[i] = 1.0;
		dual_x[i] = DualClassTemplate<complex<double>>(x[i], der);
	}

	// Evaluate system equations via DualClassTemplate auto-differentiation
	auto dual_res = system(dual_x);

	vector<vector<complex<double>>> J(n, vector<complex<double>>(n));
	vector<complex<double>> neg_F(n);

	for (size_t i = 0; i < n; ++i) 
	{
		neg_F[i] = -dual_res[i].val;
		for (size_t j = 0; j < n; ++j) 
		{
			J[i][j] = (j < dual_res[i].der.size()) ? dual_res[i].der[j] : 0.0;
		}
	}

	vector<complex<double>> dxdt(n);
	if (!solveLinearSystem(J, neg_F, dxdt)) 
	{
		return vector<complex<double>>(n, 0.0); // Fail-safe
	}
	return dxdt;
}

vector<complex<double>> RK4MultiRootSolver::solve(vector<complex<double>> x0, SystemFunc system) 
{
	double dt = 0.1;
	
	vector<complex<double>> x = x0;
	size_t n = x.size();

	for (int step = 0; step < steps; ++step) 
	{
		auto k1 = computeDerivatives(x, system);
            
		vector<complex<double>> x_k2(n);
		for (size_t i = 0; i < n; ++i)
		{ 
			x_k2[i] = x[i] + 0.5 * dt * k1[i];
		}
		auto k2 = computeDerivatives(x_k2, system);

		vector<complex<double>> x_k3(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_k3[i] = x[i] + 0.5 * dt * k2[i];
		}
		auto k3 = computeDerivatives(x_k3, system);

		vector<complex<double>> x_k4(n);
		for (size_t i = 0; i < n; ++i) 
		{
			x_k4[i] = x[i] + dt * k3[i];
		}
		auto k4 = computeDerivatives(x_k4, system);

		// Update state
		for (size_t i = 0; i < n; ++i) 
		{
			x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
		}
	}
	return x;

}

void RK4MultiRootSolver::printSolution() 
{
	cout << "========================================================\n";
	cout << "  RK4 Multiroot Solver  \n";
	cout << "========================================================\n";
	cout << "\nStep size: " << steps << endl;
	cout << endl;// Output formatting
	cout << fixed << setprecision(10);
	int root_count = 1;
	for (const auto& guess : initial_guesses) 
	{
		cout << "Initial Guess " << root_count << ":\n";
		for (size_t i = 0; i < guess.size(); ++i) 
		{
			cout << "  x" << i+1 << " = " << guess[i] << "\n";
		}
	

		// Find root using continuous Newton + RK4 Integration
		auto root = solve(guess, targetSystem);

		cout << "Computed Convergence Root:\n";
		for (size_t i = 0; i < root.size(); ++i) 
		{
			// Clean up extremely small floating values around 0 for clean output

			double real_part = std::abs(root[i].real()) < 1e-9 ? 0.0 : root[i].real();
			double imag_part = std::abs(root[i].imag()) < 1e-9 ? 0.0 : root[i].imag();
			cout << "  x" << i+1 << " = (" << real_part << ", " << imag_part << ")\n";
		}

		cout << "-------------------------------------------\n";

		root_count++;
	}
}
#endif

#endif