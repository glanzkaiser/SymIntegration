/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_DEFINE

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()

using namespace std;

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
	cout << "solution = " << p << endl;
	return 0;
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


#endif
#endif