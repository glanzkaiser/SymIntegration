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

#endif
#endif