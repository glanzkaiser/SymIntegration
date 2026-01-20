/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_GEOMETRYANDVECTORS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_GEOMETRYANDVECTORS_DEFINE

Symbolic cylindricaltocartesian(const Symbolic &r, const Symbolic &t, const Symbolic &z)
{
	Symbolic sol;
 	
	if(r != 0 || t != 0 || z != 0)
 	{
		sol = (r*cos(t), r*sin(t),  z);
		sol = sol.transpose();
	}
	return sol;
}

Symbolic cartesiantocylindrical(const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol, r;
 	double a = y, b = x;
	r = sqrt(x*x+y*y);
	if(x != 0 || y != 0 || z != 0)
 	{
		sol = (r, atanf(a/b),  z);
		sol = sol.transpose();
	}
	return sol;
}

Symbolic cartesiantospherical(const Symbolic &x, const Symbolic &y, const Symbolic &z)
{
	Symbolic sol, rho;
 	rho = sqrt(x*x+y*y+z*z);
	double a = y, b = x;
	if(x != 0 || y != 0 || z != 0)
 	{
		sol = (rho, atanf(a/b),  acos(z/rho));
		sol = sol.transpose();
	}
	return sol;
}

Symbolic sphericaltocartesian(const Symbolic &r, const Symbolic &t, const Symbolic &p)
{
	Symbolic sol;

	if(r != 0 || t != 0 || p != 0)
 	{
		sol = (r*sin(p)*cos(t), r*sin(p)*sin(t),  r*cos(p));
		sol = sol.transpose();
	}
	return sol;
}

double areaoftriangle(vector<double> &x1, vector<double> &x2, vector<double> &x3)
{
	double A;

	if (x1.size() != x2.size() || x1.size() !=  x3.size())
	{
		throw std::invalid_argument("Vectors size have to be the same.");
	}
	vector<double> a = subtract(x2,x1);
	vector<double> b = subtract(x3,x1);

	double norm_a = norm(a);
	double norm_b = norm(b);
 	double ab = dot(a,b);

	double theta = acosf(divisiond(ab,norm_a*norm_b));
	
	A = 0.5*norm_a*norm_b*sinf(theta);

	return A;
}

#endif
#endif