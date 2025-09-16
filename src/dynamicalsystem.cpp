/*
   
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DYNAMICALSYSTEM_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DYNAMICALSYSTEM_DEFINE
#define π 3.1415926535897f

Symbolic lagrangian(const Symbolic &x1, const Symbolic &y1, const Symbolic &x2, const Symbolic &y2, const Symbolic &θ1, const Symbolic &dθ1, const Symbolic &ddθ1,  const Symbolic &θ2, const Symbolic &dθ2, const Symbolic &ddθ2 ) {
	Symbolic dx1, dx2, dy1, dy2, V, T, m1("m1"), m2("m2"), g("g");
	Symbolic L10, L11, L12, L20, L21, L22, Final_Eq1, Final_Eq2;
	dx1 = df(x1,θ1)*dθ1;
	dy1 = df(y1,θ1)*dθ1;
	dx2 = df(x2,θ1)*dθ1 + df(x2,θ2)*dθ2;
	dy2 = df(y2,θ1)*dθ1 + df(y2,θ2)*dθ2;

	V = m1*g*y1 + m2*g*y2;
	T = 0.5*m1*(dx1*dx1 + dy1*dy1) + 0.5*m2*(dx2*dx2 + dy2*dy2);
	T = simplify(T,θ1);
	T = simplify(T,θ2);
	Symbolic L = T-V;

	cout <<"\nPotential energy = " << V << endl;
	cout <<"Kinetic energy = " << T << endl;
	cout << "Lagrangian = " << L << endl;

	L10 = df(L,θ1) ;
	L11 = df(L,dθ1);
	L12 = df(df(L,dθ1),θ1)*dθ1 + df(df(L,dθ1),θ2)*dθ2 + df(L,dθ1).coeff(dθ1,1)*ddθ1 + df(L,dθ1).coeff(dθ2,1)*ddθ2;
	L20 = df(L,θ2) ;
	L21 = df(L,dθ2);
	L22 = df(df(L,dθ2),θ1)*dθ1 + df(df(L,dθ2),θ2)*dθ2 + df(L,dθ2).coeff(dθ1,1)*ddθ1 + df(L,dθ2).coeff(dθ2,1)*ddθ2;
	Final_Eq1 = L12 - L10;
	Final_Eq2 = L22 - L20;

	cout << "\ndL / dθ1 = " << L10 << endl;
	cout << "\ndL / dθ1' = " << L11 << endl;
	cout << "\nd/dt (dL / dθ1') = " << L12 << endl;

	cout << "\ndL / dθ2 = " << L20 << endl;
	cout << "\ndL / dθ2' = " << L21 << endl;
	cout << "\nd/dt (dL / dθ2') = " << L22 << endl;

	cout << "\n\nFirst Equation of motion = \n" << Final_Eq1 << endl;
	cout << "\nSecond Equation of motion = " << endl;
	return Final_Eq2;
}

#endif
#endif