/*
   
*/
#include "symintegral/symintegrationc++.h"

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

// Conjugate Gradient Solver for CRS Matrix
vector<double> conjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iterations) 
{
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
#endif
#endif