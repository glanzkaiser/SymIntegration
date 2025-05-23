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

Symbolic integrate(const Symbolic &f,const Symbolic &x)
{
 list<Equations> eq;
 list<Equations>::iterator i;
 UniqueSymbol a, b, c;

 eq = (exp(a*x)).match(f, (a,b));
 for(i=eq.begin(); i!=eq.end(); ++i)
  try {
   if(df(rhs(*i, a), x) == 0) return f/rhs(*i, a);
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

