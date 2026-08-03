/*
   
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_SPECIALFUNCTIONS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_SPECIALFUNCTIONS_DEFINE
#define π 3.1415926535897f

Symbolic hypergeometric_1F1(double a, double b, const Symbolic &s, int max_iterations) 
{
	Symbolic sum = 1;
	
	for (int i = 1; i <= max_iterations; ++i) 
	{
		sum += (rising_pochhammer(a, i) * (s^(i))) / (rising_pochhammer(b, i) * factorial(i));
	}
	return sum;
}

double hypergeometric_1F1(double a, double b, double s, int max_iterations) 
{
	Symbolic sum = 1;
	
	for (int i = 1; i <= max_iterations; ++i) 
	{
		sum += (rising_pochhammer(a, i) * (s^(Symbolic(i)))) / (rising_pochhammer(b, i) * factorial(i));
	}
	return sum;
}

void probabilists_Hermitepolynomials(int n, const Symbolic &x)
{
	Symbolic He_prev = 1;
	cout << "He_{" << 0 << "} (x) = " << He_prev << endl;
	Symbolic He_current = x;
	Symbolic He_next;
	cout << "He_{" << 1 << "} (x) = " << He_current << endl;
	for(int i = 2; i <= n; ++i)
	{
		He_next = (x*He_current)- ((i-1)*He_prev);
		cout << "H_{" << i << "} (x) = " << He_next << endl;
		He_prev = He_current;
		He_current = He_next;
	}
}

void physicists_Hermitepolynomials(int n, const Symbolic &x)
{
	Symbolic H_prev = 1;
	cout << "H_{" << 0 << "} (x) = " << H_prev << endl;
	Symbolic H_current = 2*x;
	Symbolic H_next;
	cout << "H_{" << 1 << "} (x) = " << H_current << endl;
	for(int i = 2; i <= n; ++i)
	{
		H_next = (2*x*H_current)- (2*(i-1)*H_prev);
		cout << "H_{" << i << "} (x) = " << H_next << endl;
		H_prev = H_current;
		H_current = H_next;
	}
}

void Jacobianpolynomials(double alpha, double beta, int n, const Symbolic &x)
{
	Symbolic J_next;
	cout << "α = " << alpha << ", β = " << beta << endl;
	for(int i = 0; i <= n; ++i)
	{
		for(int j = 0; j <= i; ++j)
		{
			J_next += combinationsd(i+alpha,i-j)*combinationsd(i+j+alpha+beta,j)*pow(0.5*(x-1),Symbolic(j));
		}	
		cout << "P_{" << i << "} (x) = " << J_next << endl;
		J_next =  0;
	}
}

void firstkind_Chebyshevpolynomials(int n, const Symbolic &x)
{
	Symbolic T_prev = 1;
	cout << "T_{" << 0 << "} (x) = " << T_prev << endl;
	Symbolic T_current = x;
	Symbolic T_next;
	cout << "T_{" << 1 << "} (x) = " << T_current << endl;
	for(int i = 2; i <= n; ++i)
	{
		T_next = (2*x*T_current)- (T_prev);
		cout << "T_{" << i << "} (x) = " << T_next << endl;
		T_prev = T_current;
		T_current = T_next;
	}
}

void secondkind_Chebyshevpolynomials(int n, const Symbolic &x)
{
	Symbolic U_prev = 1;
	cout << "U_{" << 0 << "} (x) = " << U_prev << endl;
	Symbolic U_current = 2*x;
	Symbolic U_next;
	cout << "U_{" << 1 << "} (x) = " << U_current << endl;
	for(int i = 2; i <= n; ++i)
	{
		U_next = (2*x*U_current)- (U_prev);
		cout << "U_{" << i << "} (x) = " << U_next << endl;
		U_prev = U_current;
		U_current = U_next;
	}
}

void Gegenbauerpolynomials(double alpha, int n, const Symbolic &x)
{
	Symbolic C_prev = 1;
	cout << "C_{" << 0 << "} (x) = " << C_prev << endl;
	Symbolic C_current = 2*alpha*x;
	Symbolic C_next;
	cout << "C_{" << 1 << "} (x) = " << C_current << endl;
	for(int i = 2; i <= n; ++i)
	{
		C_next = ( (2*(i+alpha-1)*x*C_current) - ((i+2*alpha-2)*C_prev) )/(i);
		cout << "C_{" << i << "} (x) = " << C_next << endl;
		C_prev = C_current;
		C_current = C_next;
	}
}

void Laguerrepolynomials(double alpha, int n, const Symbolic &x)
{
	Symbolic L_prev = 1;
	cout << "L_{" << 0 << "} (x) = " << L_prev << endl;
	Symbolic L_current = -x+alpha+1;
	Symbolic L_next;
	cout << "L_{" << 1 << "} (x) = " << L_current << endl;
	for(int i = 2; i <= n; ++i)
	{
		L_next = ( ((2*(i-1)+alpha+1-x)*L_current) - ((i-1+alpha)*L_prev) )/(i);
		cout << "L_{" << i << "} (x) = " << L_next << endl;
		L_prev = L_current;
		L_current = L_next;
	}
}
#endif
#endif