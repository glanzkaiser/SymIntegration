/*
   
*/

#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, mu, C("C");
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, r;
		// Case 1 : ay' + ty = b
		eq = (a*x*y + d).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = exp(integrate(-ap*x,x));
		dsol = (integrate(mu*dp,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : aty' + by = ct^2	/ y' = -(a/t)y + b*t	
		eq = ((a/x)*y + d*x).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = x^(-ap);
		dsol = (integrate(mu*dp*x,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : ay' + by = c*exp(d*t)	/ y' = -ay + b*exp(c*t) 	
		eq = (a*y + b*exp(c*x)).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(-ap*x);
		dsol = (integrate(mu*bp*exp(cp*x),x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : ay' + by = c*t+d	/ y' = -ay + b*t + c 	
		eq = (a*y + b*x+c).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(-ap*x);
		dsol = (integrate(mu*(bp*x+cp),x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

	}
	return dsol;
}

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x, const Symbolic &z)
{
	Symbolic dsol, mu, gs,  C("C");
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d;
		// Will work for this model: y' + ry/b = r/c	/ y' = r/b - ry/c 	
		mu = exp(-fx.coeff(y*z,1)*x*z);
		gs = z*(fx.coeff(z,1)).coeff(y,0);
		dsol = (integrate(mu*(gs),x))/(mu) + (C)/(mu);
		
		/* // Case 1 : y' + ry/b = r/c	/ y' = r/b - ry/c 	Useless this cannot be matched. Still keep it for now maybe will be of use one day.
		eq = (b*z - c*z*y).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i,a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(cp*z*x);
		dsol = (integrate(mu*(bp*z),x))/(mu) + (C)/(mu);
		
		return dsol;
		
		} catch(const SymbolicError &se) {}
		}*/
	}
	return dsol;
}

Symbolic dsolveseparable(const Symbolic &fdy, const Symbolic &fdx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, dsol_y, dsol_x, C("C");
 	
	if(fdx == 0 && fdy !=0)
 	{
		dsol = C;
	}
	else if(fdx != 0 && fdy != 0 && fdy.coeff(x,1) == 0 && fdx.coeff(y,1) == 0)
 	{
		dsol = integrate(fdy,y) - integrate(fdx,x)  - C; 		
		
	}
	else if(fdy.coeff(x,1) !=0 && fdx.coeff(y,1) !=0 && fdy.coeff(x,2) ==0 && fdx.coeff(y,2) ==0 )
	{
		dsol_y = fdy/x;
		dsol_x = fdx/x ;		
		dsol = dsol_x / dsol_y;
		//dsol = dsol[y*(x^-1)==x];
		//dsol = dsol -(x*dsol_y[y*(x^-1)==x] / dsol_y[y*(x^-1)==x] );
		//dsol = 1/dsol;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x] - x*dsol_y[y*(x^-1)==x];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}
	else if(fdy.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 )
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}
	else if(fdx.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 && fdy.coeff(x,1) !=0 ) // is this necessary?
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}
	return dsol;
}

#endif
#endif