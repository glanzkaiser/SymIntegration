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

#include "symintegral/symintegrationc++.h"


#ifdef  SYMBOLIC_DEFINE
#ifndef SYMBOLIC_CPLUSPLUS_SIMPLIFY_DEFINE
#define SYMBOLIC_CPLUSPLUS_SIMPLIFY_DEFINE

#include <cmath>

bool isPerfectSquare(int n) {
	if (n < 0) 
	{
	return false;
	}
	int int_sqrt = static_cast<int>(sqrt(n));
	return (int_sqrt * int_sqrt == n);
}

Symbolic simplify(const Symbolic &f, const Symbolic &x)
{
	Symbolic simplified_function;
	if(f != 0 ) 	
	{	
	// Tons of trigonometry identities= very CUPULU
	Equations rules = (  (1 - (sin(x)^(2))) == (cos(x)^(2)), 
					(1 - (cos(x)^(2))) == (sin(x)^(2)),
					((csc(x)^(2)) - 1) == (cot(x)^(2)),
					(1/cos(x)) == (sec(x)),
					(1/sin(x)) == (csc(x)),
					(tan(x)) == (sin(x)/cos(x)),
					(cot(x)) == (cos(x)/sin(x)),
					(sin(x)*cos(x)) == (0.5*sin(2*x)),
					(sin(x)*sin(x)) == (0.5 - 0.5*cos(2*x)),
					(cos(x)*cos(x)) == (0.5 + 0.5*cos(2*x)),
					(1/csc(x)) == (sin(x)),
					(1/sec(x)) == (cos(x)),
					(1/(cos(x)^(2))) == (sec(x)^(2)),
					(1/(sin(x)^(2))) == (csc(x)^(2)),
					 (0.5*(exp(x) - exp(-x))) == (sinh(x)),
					(0.5*(exp(x) + exp(-x))) == (cosh(x)),
					(sinh(x)/cosh(x)) == (tanh(x)),
					(cosh(x)/sinh(x)) == (coth(x)),  
					(1/cosh(x)) == (sech(x)),
					(1/sinh(x)) == (csch(x)),
					(  (cosh(x)^2) - (sinh(x)^2) ) == 1,
					1 + (sinh(x)^2) == (cosh(x)^2), 
					(sin(acos(x))) == sqrt(1-x*x),
					(cos(asin(x))) == sqrt(1-x*x),
					sinh(x) + cosh(x) == exp(x),
					cosh(x) - sinh(x) == exp(-x),
					2*sinh(x)*cosh(x) == sinh(2*x),
					(cosh(x)^(2)) + (sinh(x)^(2)) == cosh(2*x) );
	simplified_function = f.subst_all(rules) ;
	
	}
	else if (f == 0) 
	{
	return 0;
	}

	return simplified_function; 
}

Symbolic simplify(const Symbolic &f, const Symbolic &x, const Symbolic &y)
{
	Symbolic simplified_function;
	if(f != 0 ) 	
	{	
	// Tons of trigonometry identities= very PLUTUTUT
	Equations rules = (  cos(x)*cos(y) + sin(x)*sin(y) == cos(x-y),
					 cos(x)*cos(y) - sin(x)*sin(y) == cos(x+y),
					 sin(x)*cos(y) + cos(x)*sin(y) == sin(x+y), 
					 sin(x)*cos(y) - cos(x)*sin(y) == sin(x-y) );
	simplified_function = f.subst_all(rules) ;
	
	}
	else if (f == 0) 
	{
	return 0;
	}

	return simplified_function; 
}

Symbolic simplify1(const Symbolic &f, const Symbolic &x) // TEST ONLY
{
	Symbolic simplified_function;
	
	list<Equations> eq;
	list<Equations>::iterator i;
	UniqueSymbol k, a, b, c;
	eq = ( ((sin(x))^(b)) + ((sin(x))^(c)) ).match(f, (b,c));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
		 	Symbolic bp = rhs(*i, b), cp = rhs(*i, c);
			int b_int = bp, c_int = cp;
		
			if(b_int %2 == 0 && c_int % 2 == 0 )
			{
				if(b_int < c_int)
				{
					simplified_function =  (sin(x)^(bp)) * ( (sin(x)^(bp - bp))  + (sin(x)^(cp - bp))  );	
				}
				else if(c_int < b_int)
				{
					simplified_function =  (sin(x)^(cp)) * ( (sin(x)^(bp - cp))  + (sin(x)^(cp - cp))  );	
				}
			}
			else 
			{
				return f;
			}
		
	} catch(const SymbolicError &se) {}

	return simplified_function; 
}

#endif
#endif

