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
// THANKS SENTINEL!!! and Freya too

#include "symintegral/symintegrationc++.h"


#ifdef  SYMBOLIC_DEFINE
#ifndef SYMBOLIC_CPLUSPLUS_INTEGRATE_DEFINE
#define SYMBOLIC_CPLUSPLUS_INTEGRATE_DEFINE

#define PI 3.1415926535897f

Symbolic integrate(const Symbolic &f,const Symbolic &x)
{
 list<Equations> eq;
 list<Equations>::iterator i;
 UniqueSymbol a, b, c, d, e,a1,a2,a3,a4,a5;
 Symbolic π("π");
 eq = (exp(a*x)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0) return f/rhs(*i, a);
  } catch(const SymbolicError &se) {}

 eq = (exp(a*ln(x))).match(f, (a,b)); // to handle integral of e^{a * ln(x)} but it doesn't work still symbolic error
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
   if(df(rhs(*i, a), x) == 0) return ( x^(ap+1) ) / (ap+1) ;
  } catch(const SymbolicError &se) {}

 eq = (a*exp(b*ln(x))*c*x).match(f, (a,b,c)); // to handle integral of a*e^{b * ln(x)}*c*x
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i, c);
   if(df(rhs(*i, a), x) == 0) return (ap*cp*(x^(bp+2)) ) / (bp+2) ;
  } catch(const SymbolicError &se) {}

 eq = (exp(a*x*x)).match(f, (a,b)); // to handle integral of e^{ax^2}
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic ap = rhs(*i, a);
   if(df(rhs(*i, a), x) == 0) return (sqrt(π)*erf(x*sqrt(-ap)) ) / (2*sqrt(-ap)) ;
  } catch(const SymbolicError &se) {}

 eq = (b*exp(a*x)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
    if(df(bp, x) == 0) return f/ap;
   }
  } catch(const SymbolicError &se) {}

 eq = (b*x*exp(a*x)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
    if(df(bp, x) == 0) return f/ap - bp*exp(ap*x)/(ap^2);
   }
  } catch(const SymbolicError &se) {}

 eq = ((x^b)*exp(x)).match(f, (a,b)); // case for x^b * exp(x)
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic bp = rhs(*i, b);
   if(bp.type() == typeid(Numeric)
      && Number<void>(bp)->numerictype() == typeid(int)
      && Number<void>(bp)>Number<int>(0))
   {
    int n = CastPtr<const Number<int> >(bp)->n, sgn = 1;
    Symbolic integral, nf = 1;
    for(; n>=0; nf*=n, --n, sgn=-sgn) integral += sgn*nf*(x^n)*exp(x);
    return integral;
   }
  } catch(const SymbolicError &se) {}

 eq = ( cos(a*x)*sin(b*x) ).match(f, (a,b)); // case for cos(a*x) * sin(b*x)
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
   if(df(rhs(*i, a), x) == 0)
   {
	if(ap == bp || bp == ap || bp == -ap || ap == -bp)
	   {
	    return (sin(bp*x)*sin(bp*x))/(2*bp) ;
	    //return -(cos(bp*x)*cos(bp*x))/(2*bp) ;
	   }
	if(ap == bp && ap ==1)
	  {
	  return 0.5*sin(x)*sin(x);
	  }
	return (ap*sin(ap*x)*sin(bp*x))/(ap*ap - bp*bp) + (bp*cos(ap*x)*cos(bp*x))/(ap*ap - bp*bp) ;
   }
  } catch(const SymbolicError &se) {}

 eq = ( cos(a*x)*cos(b*x) ).match(f, (a,b)); // case for cos(a*x) * cos(b*x)
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
   if(df(rhs(*i, a), x) == 0)
   {
	if(ap == bp || bp == ap || bp == -ap || ap == -bp)
	   {
	    return (2*bp*x + sin(2*bp*x)) / (4*bp) ; //weirdly, it cannot be returned from here, but from 'Symbolic Power::integrate(const Symbolic &s) const' in functions.cpp instead
	   }
	return ap*sin(ap*x)*cos(bp*x)/(ap*ap - bp*bp) - bp*sin(bp*x)*cos(ap*x)/(ap*ap - bp*bp); 	
	//return 0.5 * ( (1/(ap+bp)) * (sin((ap+bp)*x))  + (1/(ap-bp)) * (sin((ap-bp)*x)) );
   }
  } catch(const SymbolicError &se) {}

 eq = ( sin(a*x)*sin(b*x) ).match(f, (a,b)); // case for sin(a*x) * sin(b*x)
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
   if(df(rhs(*i, a), x) == 0)
   {
	if(ap == bp || bp == ap || bp == -ap || ap == -bp)
	   {
	    return 0.5*x - (sin(2*bp*x) / (4*bp)) ;
	   }
	return (bp*sin(ap*x)*cos(bp*x))/(ap*ap - bp*bp) - (ap*cos(ap*x)*sin(bp*x))/(ap*ap - bp*bp) ;
   }
  } catch(const SymbolicError &se) {}

 eq = ((x^b)*exp(a*x)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
    if(bp.type() == typeid(Numeric)
       && Number<void>(bp)->numerictype() == typeid(int)
       && Number<void>(bp)>Number<int>(0))
    {
     int n = CastPtr<const Number<int> >(bp)->n, sgn = 1;
     Symbolic integral, nf = 1/ap;
     for(; n>=0; nf*=(n)/ap, --n, sgn=-sgn)
      integral += sgn*nf*(x^n)*exp(ap*x);
     return integral;
    }
   }
  } catch(const SymbolicError &se) {}

 eq = (c*(x^b)*exp(a*x)).match(f, (a,b,c));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0 && df(rhs(*i, c), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
    if(bp.type() == typeid(Numeric)
       && Number<void>(bp)->numerictype() == typeid(int)
       && Number<void>(bp)>Number<int>(0))
    {
     int n = CastPtr<const Number<int> >(bp)->n, sgn = 1;
     Symbolic integral, nf = 1/ap;
     for(; n>=0; nf*=(n)/ap, --n, sgn=-sgn)
      integral += sgn*nf*(x^n)*exp(ap*x);
     return rhs(*i, c)*integral;
    }
   }
  } catch(const SymbolicError &se) {}

 eq = (b*exp((a+c)*x)).match(f, (a,b,c));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0 && df(rhs(*i, c), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i, c);
    if(df(bp, x) == 0) return f/(ap+cp);
    if(df(bp, x, 5) == 0)
     return f/(ap+cp) - integrate(df(bp,x)*exp((ap+cp)*x)/(ap+cp), x);
   }
  } catch(const SymbolicError &se) {}

 eq = (b*x*exp(a*x*x)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
    if(df(bp, x) == 0) return f/(2*ap);
   }
  } catch(const SymbolicError &se) {}

 eq = (1/sqrt(a*x+b)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0 && df(rhs(*i, b), x) == 0)
   {
    Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
    return 2*sqrt(ap*x+bp)/ap;
   }
  } catch(const SymbolicError &se) {}

 return f.integrate(x);
}

Symbolic fractionintegrate(const Symbolic &fnum, const Symbolic &fdenom, const Symbolic &x)
{
  list<Equations> eq;
  list<Equations>::iterator i;
  UniqueSymbol anum, bnum, cnum;
  Symbolic integral_sol;
  double anump, bnump, cnump, ap, bp, cp;
  ap = fdenom.coeff(x,2);
  bp = fdenom.coeff(x,1);
  cp = fdenom.coeff(x,0); 

  if (fnum.coeff(x,2)!=0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)!=0 ) // for (a1x^2 + b1x + c1)/(ax^2+bx+c)
  {
	eq = ( (anum*x*x+bnum*x+cnum)).match(fnum, (anum,bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, anum), b1 = rhs(*i, bnum), c1 = rhs(*i, cnum); 
 	anump= a1;
 	bnump= b1;
 	cnump = c1;
 	} catch(const SymbolicError &se) {}

	double D = sqrt(bp*bp-(4*ap*cp))*(2*ap*ap*cnump-2*ap*anump*cp-ap*bp*bnump+anump*bp*bp)/(2*ap*ap*(4*ap*cp-(bp*bp))) ;
 
	integral_sol = ((ap*bnump-anump*bp)/(2*ap*ap) - D)*ln(x+(4*ap*ap*cp*((ap*bnump-anump*bp)/(2*ap*ap) - D) - ap*bp*bp*((ap*bnump-anump*bp)/(2*ap*ap) - D) + ap*bp*cnump-2*ap*bnump*cp+anump*bp*cp) / (2*ap*ap*cnump-2*ap*anump*cp-ap*bp*bnump+anump*bp*bp)) + ((ap*bnump-anump*bp)/(2*ap*ap) + D)*ln(x+(4*ap*ap*cp*((ap*bnump-anump*bp)/(2*ap*ap) + D) - ap*bp*bp*((ap*bnump-anump*bp)/(2*ap*ap) + D) + ap*bp*cnump-2*ap*bnump*cp+anump*bp*cp) / (2*ap*ap*cnump-2*ap*anump*cp-ap*bp*bnump+anump*bp*bp)) + (anump*x)/ap;   
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)!=0 ) // for (a1x+b1)/(ax^2+bx+c)
  {
	eq = ( (bnum*x+cnum)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
 	cnump = b1;
 	} catch(const SymbolicError &se) {}

	double D = (2*ap*cnump - bnump*bp)*sqrt(-4*ap*cp+(bp*bp))/(2*ap*(4*ap*cp-(bp*bp)));
	double denom1 = 2*ap*cnump - bnump*bp;	
 
	integral_sol = ((bnump/(2*ap)) - D)*ln(x + (4*ap*cp*((bnump/(2*ap))-D) - 2*bnump*cp -(bp*bp)*(bnump/(2*ap) - D) + bp*cnump)/(denom1)) + ((bnump/(2*ap)) + D)*ln(x + (4*ap*cp*((bnump/(2*ap))+D) - 2*bnump*cp -(bp*bp)*(bnump/(2*ap) + D) + bp*cnump)/(denom1)) ;   
  }
  
  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)!=0 ) // for (a1x)/(ax^2+bx+c)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	double D = (bp*sqrt(bp*bp-4*ap*cp))/(2*ap*(4*ap*cp-(bp*bp))) ;
	integral_sol = bnump*( (-D + (1/(2*ap)))*ln(x+(-4*ap*cp*(-D + 1/(2*ap)) +bp*bp*(-D+ 1/(2*ap)) + 2*cp)/(bp)) + (D + (1/(2*ap)))*ln(x+(-4*ap*cp*(D + 1/(2*ap)) +bp*bp*(D+ 1/(2*ap)) + 2*cp)/(bp)) );
  }
  
  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)==0 && fdenom.coeff(x,0)!=0 ) // for (a1x + b1)/(ax^2+c)
  {
	eq = ( (bnum*x+cnum)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
 	cnump = b1; 
 	} catch(const SymbolicError &se) {}

	double D = cnump*sqrt(-ap*ap*ap*cp)/(2*ap*ap*cp);
	integral_sol = (bnump/(2*ap) - D)*ln(x+(2*ap*cp*(bnump/(2*ap) - D) - bnump*cp)/(ap*cnump)) + (bnump/(2*ap) + D)*ln(x+(2*ap*cp*(bnump/(2*ap) + D) - bnump*cp)/(ap*cnump));
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)==0 && fdenom.coeff(x,0)!=0 ) // for (a1x)/(ax^2+c)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*( ln(ap*x*x+cp)/(2*ap) );
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)==0 ) // for (a1x + b1)/(ax^2+bx)
  {
	eq = ( (bnum*x+cnum)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
 	cnump = b1; 
 	} catch(const SymbolicError &se) {}

	double D = (bp*cnump+(bp*(ap*cnump-bnump*bp))/(ap))/(2*ap*cnump-bnump*bp);
	integral_sol = cnump*ln(x)/bp - ((ap*cnump-bnump*bp)*ln(x+D))/(ap*bp);
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)==0 ) // for (a1x)/(ax^2+bx)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*( ln(ap*x+bp)/(ap) );
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)==0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)!=0 ) // for (a1x + b1)/(bx+c)
  {
	eq = ( (bnum*x+cnum)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
 	cnump = b1; 
 	} catch(const SymbolicError &se) {}

	integral_sol = (bnump*x/bp) - ((bnump*cp-bp*cnump)*ln(bp*x+cp))/(bp*bp);
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)==0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)!=0 ) // for (a1x)/(bx+c)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*(x/bp - cp*ln(bp*x+cp)/(bp*bp) );
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)==0 && fdenom.coeff(x,0)==0 ) // for (a1x+b1)/(ax^2)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
	cnump = b1; 
 	} catch(const SymbolicError &se) {}

	integral_sol = (bnump*ln(x) - (cnump)/(x))/(ap);
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)!=0 && fdenom.coeff(x,1)==0 && fdenom.coeff(x,0)==0 ) // for (a1x)/(ax^2)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*ln(x)/ap;
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)==0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)==0 ) // for (a1x+b1)/(bx)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
	cnump = b1; 
 	} catch(const SymbolicError &se) {}

	integral_sol = (bnump*x +cnump*ln(x))/(bp);
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)==0 && fdenom.coeff(x,1)!=0 && fdenom.coeff(x,0)==0 ) // for (a1x)/(bx)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*x/bp;
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)!=0 && fdenom.coeff(x,2)==0 && fdenom.coeff(x,1)==0 && fdenom.coeff(x,0)!=0 ) // for (a1x+b1)/(c)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum), b1 = rhs(*i, cnum); 
 	bnump= a1;
	cnump = b1; 
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*x*x/(2*cp) + cnump*x/cp;
  }

  if (fnum.coeff(x,2)==0 && fnum.coeff(x,1)!=0 && fnum.coeff(x,0)==0 && fdenom.coeff(x,2)==0 && fdenom.coeff(x,1)==0 && fdenom.coeff(x,0)!=0 ) // for (a1x)/(c)
  {
	eq = ( (bnum*x)).match(fnum, (bnum,cnum)); 
 	for(i=eq.begin(); i!=eq.end(); ++i) 
 	try {
 	Symbolic a1 = rhs(*i, bnum); 
 	bnump= a1;
 	} catch(const SymbolicError &se) {}

	integral_sol = bnump*x*x/(2*cp);
  }

 return integral_sol;
}

Symbolic integrate(const Symbolic &f,const Symbolic &x,
                   const Symbolic &a, const Symbolic &b)
{
 Symbolic I = integrate(f,x);
 return I[x == b] - I[x == a];
}

Symbolic integrate(const Symbolic &s,const Symbolic &x,unsigned int i)
{
 Symbolic r = s;
 while(i-- > 0) r = integrate(r,x);
 return r;
}


#endif
#endif

