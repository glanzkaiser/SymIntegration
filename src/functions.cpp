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
#include <cmath>

using namespace std;
#ifdef  SYMBOLIC_DEFINE
#ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS_DEFINE
#define SYMBOLIC_CPLUSPLUS_FUNCTIONS_DEFINE
#define SYMBOLIC_CPLUSPLUS_FUNCTIONS

//////////////////////////////////////
// Implementation of Sin            //
//////////////////////////////////////

Sin::Sin(const Sin &s) : Symbol(s) {}

Sin::Sin(const Symbolic &s) : Symbol(Symbol("sin")[s]) {}

Asin::Asin(const Asin &s) : Symbol(s) {}

Asin::Asin(const Symbolic &s) : Symbol(Symbol("asin")[s]) {}

Acos::Acos(const Acos &s) : Symbol(s) {}

Acos::Acos(const Symbolic &s) : Symbol(Symbol("acos")[s]) {}

Simplified Sin::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 const Symbolic &ac = Symbol(Symbol("acos")[s]) ;
 //const Symbolic &x = parameters.front();
 if(s == 0) return Number<int>(0);
 //if(s == ac) return Number<int>(8);
 //if(parameters.front() == Symbol("acos")[s])) return Number<int>(8);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Sin(-s);
 // if(p ->factors.front() == acos(s)) return sqrt(1-(parameters.front()*parameters.front()));
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(sin(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Sin::df(const Symbolic &s) const
{ return cos(parameters.front()) * parameters.front().df(s); }

Symbolic Sin::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return -cos(x) ;
 if(df(s) == 0) return *this * s;
 return -cos(parameters.front()) * (1/ parameters.front().df(s) );
}

//////////////////////////////////////
// Implementation of Asin            //
//////////////////////////////////////

//Asin::Asin(const Asin &s) : Symbol(s) {}

//Asin::Asin(const Symbolic &s) : Symbol(Symbol("asin")[s]) {}

Simplified Asin::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Asin(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(asin(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Asin::df(const Symbolic &s) const
{ return parameters.front().df(s) / (sqrt(1 - (parameters.front())*(parameters.front())) ); }

Symbolic Asin::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x*asin(x) + sqrt(1-x*x) ;
 if(df(s) == 0) return *this * s;
 return ( sqrt(1-(parameters.front())*(parameters.front()) ) /  parameters.front().df(s) ) + ( parameters.front()*asin(parameters.front()) /  parameters.front().df(s) );
}

//////////////////////////////////////
// Implementation of Cos            //
//////////////////////////////////////

Cos::Cos(const Cos &s) : Symbol(s) {}

Cos::Cos(const Symbolic &s) : Symbol(Symbol("cos")[s]) {}

Simplified Cos::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(1);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return Cos(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(cos(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Cos::df(const Symbolic &s) const
{ return -sin(parameters.front()) * parameters.front().df(s); }

Symbolic Cos::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return sin(x);
 if(df(s) == 0) return *this * s;
 return sin(parameters.front()) * ( 1 / (parameters.front().df(s)) ) ;
}

//////////////////////////////////////
// Implementation of Acos            //
//////////////////////////////////////

//Acos::Acos(const Acos &s) : Symbol(s) {}

//Acos::Acos(const Symbolic &s) : Symbol(Symbol("acos")[s]) {}

Simplified Acos::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<double>(1.5707963267948966);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return SymbolicConstant::pi - Acos(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(acos(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Acos::df(const Symbolic &s) const
{ return - parameters.front().df(s) / (sqrt(1 - (parameters.front())*(parameters.front())) ); }

Symbolic Acos::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x*acos(x) - sqrt(1-x*x) ;
 if(df(s) == 0) return *this * s;
 return ( - sqrt(1-(parameters.front())*(parameters.front()) ) /  parameters.front().df(s) ) + ( parameters.front()*acos(parameters.front()) /  parameters.front().df(s) );
}


//////////////////////////////////////
// Implementation of Tan           //
//////////////////////////////////////

Tan::Tan(const Tan &s) : Symbol(s) {}

Tan::Tan(const Symbolic &s) : Symbol(Symbol("tan")[s]) {}

Simplified Tan::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Tan(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(tan(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Tan::df(const Symbolic &s) const
{ return tan(parameters.front()) * tan(parameters.front()) * parameters.front().df(s) + parameters.front().df(s) ; }

Symbolic Tan::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return -ln(cos(x)) * (1 / parameters.front().df(s));
 if(parameters.front().coeff(s,1) != 0 && parameters.front().coeff(s,0) ==0) return -( ln(cos((parameters.front().coeff(s,1))*s)) ) * (1 / (parameters.front().coeff(s,1)) ) ;

 if(df(s) == 0) return *this * s;
 return ln( tan(parameters.front()) * tan(parameters.front()) + 1) * ( 1 / (2*parameters.front().df(s)) ) ;
}

//////////////////////////////////////
// Implementation of Atan           //
//////////////////////////////////////

Atan::Atan(const Atan &s) : Symbol(s) {}

Atan::Atan(const Symbolic &s) : Symbol(Symbol("atan")[s]) {}

Simplified Atan::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Atan(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(atan(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Atan::df(const Symbolic &s) const
{ return parameters.front().df(s) / ( 1 + (parameters.front())*(parameters.front()) ); }

Symbolic Atan::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x*atan(x) - 0.5*ln(1+x*x) ;
 if(df(s) == 0) return *this * s;
 return ( (x/parameters.front().df(s))*atan(parameters.front()) - ln(parameters.front()*parameters.front()+1) / (2*parameters.front().df(s)) + parameters.front().coeff(s,0)*(atan(parameters.front())) / (parameters.front().df(s))  );
}


//////////////////////////////////////
// Implementation of Cot           //
//////////////////////////////////////

Cot::Cot(const Cot &s) : Symbol(s) {}

Cot::Cot(const Symbolic &s) : Symbol(Symbol("cot")[s]) {}

Simplified Cot::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Symbol(Symbol("Inf"));
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Cot(-s);
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(double))
 {
	return Number<double>(1) / Number<double>(tan(CastPtr<const Number<double> >(s)->n));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(int))
 {
	return Number<double>(1) / Number<double>(tan(CastPtr<const Number<int> >(s)->n));
 }
 return *this;
}

Symbolic Cot::df(const Symbolic &s) const
{ return -cot(parameters.front()) * cot(parameters.front()) * parameters.front().df(s) - parameters.front().df(s) ; }

Symbolic Cot::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return ln(sin(x)) * (1 / parameters.front().df(s));
 if(parameters.front().coeff(s,1) != 0 && parameters.front().coeff(s,0) ==0) return ( ln(sin((parameters.front().coeff(s,1))*s)) ) * (1 / (parameters.front().coeff(s,1)) ) ;

 if(df(s) == 0) return *this * s;
 return -ln(tan(parameters.front()) * tan(parameters.front()) + 1) * ( 1 / (2*parameters.front().df(s)) ) + ln(tan(parameters.front())) * (1 / (parameters.front().df(s)) ) ;
}

//////////////////////////////////////
// Implementation of Acot           //
//////////////////////////////////////

Acot::Acot(const Acot &s) : Symbol(s) {}

Acot::Acot(const Symbolic &s) : Symbol(Symbol("acot")[s]) {}

Simplified Acot::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<double>(1.5707963267948966);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Acot(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(acot(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Acot::df(const Symbolic &s) const
{ return -parameters.front().df(s) / ( 1 + (parameters.front())*(parameters.front()) ); }

Symbolic Acot::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x*acot(x) + 0.5*ln(1+x*x) ;
 if(df(s) == 0) return *this * s;
 return ( (x/x.df(s)) * acot(parameters.front()) + ln(parameters.front()*parameters.front()+1) / (2*parameters.front().df(s)) + parameters.front().coeff(s,0)*(acot(parameters.front())) / (parameters.front().df(s))  );
}

//////////////////////////////////////
// Implementation of Sec            //
//////////////////////////////////////

Sec::Sec(const Sec &s) : Symbol(s) {}

Sec::Sec(const Symbolic &s) : Symbol(Symbol("sec")[s]) {}

Simplified Sec::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(1);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return Sec(-s);
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(double))
 {
	return Number<double>(1) / Number<double>(cos(CastPtr<const Number<double> >(s)->n));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(int))
 {
	return Number<double>(1) / Number<double>(cos(CastPtr<const Number<int> >(s)->n));
 }
 return *this;
}

Symbolic Sec::df(const Symbolic &s) const
{ return parameters.front().df(s) * tan(parameters.front()) * sec(parameters.front()); }

Symbolic Sec::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return -0.5*(ln(sin(x)-1)) + 0.5*(ln(sin(x)+1)) ;
 if(parameters.front().coeff(s,1) != 0 && parameters.front().coeff(s,0) ==0) return (-0.5*(ln(sin((parameters.front().coeff(s,1))*s)-1)) + 0.5*(ln(sin((parameters.front().coeff(s,1))*s)+1)) ) * (1 / (parameters.front().df(s)) ) ;
 if(df(s) == 0) return *this * s;
 return ln(tan(parameters.front()) + sec(parameters.front()) ) * (1 / (parameters.front().df(s)) ) ;
}

//////////////////////////////////////
// Implementation of Asec           //
//////////////////////////////////////

Asec::Asec(const Asec &s) : Symbol(s) {}

Asec::Asec(const Symbolic &s) : Symbol(Symbol("asec")[s]) {}

Simplified Asec::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Symbol(Symbol("Inf"));
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return SymbolicConstant::pi - Asec(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(asec(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Asec::df(const Symbolic &s) const
{ return (parameters.front().df(s)) / (parameters.front()*parameters.front()*sqrt(1-(1/(parameters.front()*parameters.front())))) ; }

Symbolic Asec::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) 
 {
	cout << "( for |x^2| > 1 )"<< endl ;
	return x*asec(x) - acosh(x) ;
 }
 if(df(s) == 0) return *this * s;
 if(parameters.front().coeff(s,0)==0 && parameters.front().coeff(s,1) !=0)
 {
	cout << "( for |x^2| > 1/"<< parameters.front().coeff(s,1) * parameters.front().coeff(s,1) << " )"<< endl ;
	return (x/x.df(s))*asec(parameters.front()) - acosh(parameters.front()) ;
 }
 return Integral(*this,s);
}


//////////////////////////////////////
// Implementation of Csc            //
//////////////////////////////////////

Csc::Csc(const Csc &s) : Symbol(s) {}

Csc::Csc(const Symbolic &s) : Symbol(Symbol("csc")[s]) {}

Simplified Csc::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(1);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Csc(-s);
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(double))
 {
	return Number<double>(1) / Number<double>(sin(CastPtr<const Number<double> >(s)->n));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(int))
 {
	return Number<double>(1) / Number<double>(sin(CastPtr<const Number<int> >(s)->n));
 }
 return *this;
}

Symbolic Csc::df(const Symbolic &s) const
{ return -parameters.front().df(s) * cot(parameters.front()) * csc(parameters.front()); }

Symbolic Csc::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return 0.5*(ln(cos(x)-1)) - 0.5*(ln(cos(x)+1)) ;
 if(parameters.front().coeff(s,1) != 0 && parameters.front().coeff(s,0) ==0) return (0.5*(ln(cos((parameters.front().coeff(s,1))*s)-1)) - 0.5*(ln(cos((parameters.front().coeff(s,1))*s)+1)) ) * (1 / (parameters.front().df(s)) ) ;
 if(df(s) == 0) return *this * s;
 return -ln(cot(parameters.front()) + csc(parameters.front()) ) * (1 / (parameters.front().df(s)) ) ;
}

//////////////////////////////////////
// Implementation of Acsc          //
//////////////////////////////////////

Acsc::Acsc(const Acsc &s) : Symbol(s) {}

Acsc::Acsc(const Symbolic &s) : Symbol(Symbol("acsc")[s]) {}

Simplified Acsc::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Symbol(Symbol("Inf"));
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return - Acsc(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(acsc(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Acsc::df(const Symbolic &s) const
{ return -(parameters.front().df(s)) / (parameters.front()*parameters.front()*sqrt(1-(1/(parameters.front()*parameters.front())))) ; }

Symbolic Acsc::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) 
 {
	cout << "( for |x^2| > 1 )"<< endl ;
	return x*acsc(x) + acosh(x) ;
 }
 if(df(s) == 0) return *this * s;
 if(parameters.front().coeff(s,0)==0 && parameters.front().coeff(s,1) !=0)
 {
	cout << "( for |x^2| > 1/"<< parameters.front().coeff(s,1) * parameters.front().coeff(s,1) << " )"<< endl ;
	return (x/x.df(s))*acsc(parameters.front()) + acosh(parameters.front()) ;
 }
 return Integral(*this,s);
}

//////////////////////////////////////
// Implementation of Sinh           //
//////////////////////////////////////

Sinh::Sinh(const Sinh &s) : Symbol(s) {}

Sinh::Sinh(const Symbolic &s) : Symbol(Symbol("sinh")[s]) {}

Simplified Sinh::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Sinh(-s);
  if(p->factors.front() == 1) return 0.5*(exp(s) - exp(-s));
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(sinh(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Sinh::df(const Symbolic &s) const
{ return cosh(parameters.front()) * parameters.front().df(s); }

Symbolic Sinh::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return cosh(x);
 if(df(s) == 0) return *this * s;
 return cosh(parameters.front()) * (1/ parameters.front().df(s)) ;
}

//////////////////////////////////////
// Implementation of Asinh           //
//////////////////////////////////////

Asinh::Asinh(const Asinh &s) : Symbol(s) {}

Asinh::Asinh(const Symbolic &s) : Symbol(Symbol("asinh")[s]) {}

Simplified Asinh::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Asinh(-s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(asinh(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Asinh::df(const Symbolic &s) const
{ return parameters.front().df(s) / sqrt(parameters.front()*parameters.front() + 1); }

Symbolic Asinh::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x*asinh(x) - sqrt(x*x+1);
 if(df(s) == 0) return *this * s;
 return ( x/x.df(s)) * asinh(parameters.front()) - (sqrt(parameters.front()*parameters.front() + 1)) / (parameters.front().df(s)) + (parameters.front().coeff(s,0) * asinh(parameters.front())) / (parameters.front().df(s))  ;
}

//////////////////////////////////////
// Implementation of Cosh           //
//////////////////////////////////////

Cosh::Cosh(const Cosh &s) : Symbol(s) {}

Cosh::Cosh(const Symbolic &s) : Symbol(Symbol("cosh")[s]) {}

Simplified Cosh::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(1);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return Cosh(-s);
  if(p->factors.front() == 1) return 0.5*(exp(s) + exp(-s));
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(cosh(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Cosh::df(const Symbolic &s) const
{ return sinh(parameters.front()) * parameters.front().df(s); }

Symbolic Cosh::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return sinh(x);
 if(df(s) == 0) return *this * s;
 return sinh(parameters.front()) * (1/ parameters.front().df(s)) ;
}


//////////////////////////////////////
// Implementation of Acosh           //
//////////////////////////////////////

Acosh::Acosh(const Acosh &s) : Symbol(s) {}

Acosh::Acosh(const Symbolic &s) : Symbol(Symbol("acosh")[s]) {}

Simplified Acosh::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Symbol(Symbol("Undefined")[s]);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return Acosh(s);
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(acosh(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Acosh::df(const Symbolic &s) const
{ return parameters.front().df(s) / sqrt(parameters.front()-1) * sqrt(parameters.front() + 1); }

Symbolic Acosh::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return Integral(*this,s);
 if(df(s) == 0) return *this * s;
 return Integral(*this,s);
}


//////////////////////////////////////
// Implementation of Tanh           //
//////////////////////////////////////

Tanh::Tanh(const Tanh &s) : Symbol(s) {}

Tanh::Tanh(const Symbolic &s) : Symbol(Symbol("tanh")[s]) {}

Simplified Tanh::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(1);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Tanh(-s);
  if(p->factors.front() == 1) return (exp(s) - exp(-s))/(exp(s) + exp(-s));
 }
 if(s.type() == typeid(Numeric) &&
    Number<void>(s).numerictype() == typeid(double))
  return Number<double>(tanh(CastPtr<const Number<double> >(s)->n));
 return *this;
}

Symbolic Tanh::df(const Symbolic &s) const
{ return parameters.front().df(s) - parameters.front().df(s) * ( tanh(parameters.front()) * tanh(parameters.front()) ); }

Symbolic Tanh::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x - ln(tanh(x) + 1);
 if(df(s) == 0) return *this * s;
 return ( x/x.df(s)) - ( ln(tanh(parameters.front()) + 1) ) * (1/ parameters.front().df(s)) ;
}

//////////////////////////////////////
// Implementation of Coth           //
//////////////////////////////////////

Coth::Coth(const Coth &s) : Symbol(s) {}

Coth::Coth(const Symbolic &s) : Symbol(Symbol("coth")[s]) {}

Simplified Coth::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(1);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Coth(-s);
  if(p->factors.front() == 1) return (exp(s) + exp(-s))/(exp(s) - exp(-s));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(double))
 {
	return Number<double>(1) / Number<double>(tanh(CastPtr<const Number<double> >(s)->n));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(int))
 {
	return Number<double>(1) / Number<double>(tanh(CastPtr<const Number<int> >(s)->n));
 }
 return *this;
}

Symbolic Coth::df(const Symbolic &s) const
{ return - parameters.front().df(s) * (1 / ( sinh(parameters.front()) * sinh(parameters.front()) ) ); }

Symbolic Coth::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return x - ln(tanh(x) + 1) + ln(tanh(x));
 if(df(s) == 0) return *this * s;
 return ( x/x.df(s)) - ( ( ln(tanh(parameters.front()) + 1) ) * (1/ parameters.front().df(s)) ) + ( ( ln(tanh(parameters.front())) ) * (1/ parameters.front().df(s)) );
}

//////////////////////////////////////
// Implementation of Sech           //
//////////////////////////////////////

Sech::Sech(const Sech &s) : Symbol(s) {}

Sech::Sech(const Symbolic &s) : Symbol(Symbol("sech")[s]) {}

Simplified Sech::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return Sech(-s);
  if(p->factors.front() == 1) return 2 / (exp(s) + exp(-s));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(double))
 {
	return Number<double>(1) / Number<double>(cosh(CastPtr<const Number<double> >(s)->n));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(int))
 {
	return Number<double>(1) / Number<double>(cosh(CastPtr<const Number<int> >(s)->n));
 }
 return *this;
}

Symbolic Sech::df(const Symbolic &s) const
{ return ( - parameters.front().df(s) ) * tanh(parameters.front()) * sech(parameters.front()); }

Symbolic Sech::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return 2 * atan(tanh(0.5*x));
 if(df(s) == 0) return *this * s;
 return 2 * (atan(tanh(0.5*(parameters.front()))))  * (1/ parameters.front().df(s)) ;
}


//////////////////////////////////////
// Implementation of Csch           //
//////////////////////////////////////

Csch::Csch(const Csch &s) : Symbol(s) {}

Csch::Csch(const Symbolic &s) : Symbol(Symbol("csch")[s]) {}

Simplified Csch::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 if(s.type() == typeid(Product))
 {
  CastPtr<const Product> p(s);
  if(p->factors.front() == -1) return -Csch(-s);
  if(p->factors.front() == 1) return 2 / (exp(s) - exp(-s));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(double))
 {
	return Number<double>(1) / Number<double>(sinh(CastPtr<const Number<double> >(s)->n));
 }
 if(s.type() == typeid(Numeric) && Number<void>(s).numerictype() == typeid(int))
 {
	return Number<double>(1) / Number<double>(sinh(CastPtr<const Number<int> >(s)->n));
 }
 return *this;
}

Symbolic Csch::df(const Symbolic &s) const
{ return ( - parameters.front().df(s) ) * coth(parameters.front()) * csch(parameters.front()); }

Symbolic Csch::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return 2 * atan(tanh(0.5*x));
 if(df(s) == 0) return *this * s;
 return (ln(tanh(0.5*(parameters.front()))))  * (1/ parameters.front().df(s)) ;
}

//////////////////////////////////////
// Implementation of Erf           //
//////////////////////////////////////

Erf::Erf(const Erf &s) : Symbol(s) {}

Erf::Erf(const Symbolic &s) : Symbol(Symbol("erf")[s]) {}

Simplified Erf::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s == 0) return Number<int>(0);
 return *this;
}

Symbolic Erf::df(const Symbolic &s) const
{ return Derivative(*this,s); }

Symbolic Erf::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.front();
 if(x == s) return Integral(*this,s);
 if(df(s) == 0) return Integral(*this,s);
 return Integral(*this,s);
}

//////////////////////////////////////
// Implementation of Log            //
//////////////////////////////////////

Log::Log(const Log &s) : Symbol(s) {}

Log::Log(const Symbolic &s1,const Symbolic &s2)
: Symbol(Symbol("log")[s1,s2]) {}

void Log::print(ostream &o) const
{
 if(parameters.size() == 2 && parameters.front() == SymbolicConstant::e)
 {
  Log l = *this;
  l.name = "ln";
  l.parameters.pop_front();
  l.print(o);
 }
 else
  Symbol::print(o);
}

Simplified Log::simplify() const
{
 // log_a(b)
 const Symbolic &a = parameters.front().simplify();
 const Symbolic &b = parameters.back().simplify();
 if(b == 1) return Number<int>(0);
 if(b == a) return Number<int>(1);
 if(b.type() == typeid(Power))
 {
  CastPtr<const Power> p = b;
  if(p->parameters.front() == a)
   return p->parameters.back();
 }
 if(b.type() == typeid(Numeric)                     &&
    Number<void>(b).numerictype() == typeid(double) &&
    CastPtr<const Number<double> >(b)->n > 0.0)
  return Product(Number<double>(log(CastPtr<const Number<double> >(b)->n)),
                 Power(ln(a),-1)).simplify();
 return *this;
}

// d/ds log_a(b) = d/ds (ln(b) / ln(a))
//               = (1/b db/ds - log_a(b) / a da/ds) / ln(a)
//               = (a db/ds - b log_a(b) da/ds) / (a b ln(a))
Symbolic Log::df(const Symbolic &s) const
{
 const Symbolic &a = parameters.front();
 const Symbolic &b = parameters.back();
 return (a * b.df(s) - b * *this * a.df(s)) / (a * b * ln(a));
}

Symbolic Log::integrate(const Symbolic &s) const
{
 const Symbolic &x = parameters.back();
 const Symbolic &a = parameters.front();
 // int(log_a(x)) = (x ln(x) - x) / ln(a)
 if(x == s && a.df(s) == 0)
  return (x * *this - x) / ln(a);
 if(df(s) == 0) return *this * s;
 return Integral(*this,s);
}

//////////////////////////////////////
// Implementation of Power          //
//////////////////////////////////////

Power::Power(const Power &s) : Symbol(s) {}

Power::Power(const Symbolic &s,const Symbolic &p) : Symbol("pow")
{ parameters.push_back(s); parameters.push_back(p); }

void Power::print(ostream &o) const
{
  if(*this == SymbolicConstant::i)
  { o << SymbolicConstant::i_symbol; return; }

  int parens1 = parameters.front().type() == typeid(Symbol)
             || parameters.front().type() == typeid(Sin)
	     || parameters.front().type() == typeid(Cos)
             || parameters.front().type() == typeid(Sinh)
             || parameters.front().type() == typeid(Cosh)
             || parameters.front().type() == typeid(Log)
             || parameters.front().type() == typeid(Derivative);
  int parens2 = parameters.back().type() == typeid(Symbol)
             || parameters.back().type() == typeid(Sin)
             || parameters.back().type() == typeid(Cos)
             || parameters.back().type() == typeid(Sinh)
             || parameters.back().type() == typeid(Cosh)
             || parameters.back().type() == typeid(Log)
             || parameters.back().type() == typeid(Derivative);
  parens1 = !parens1;
  parens2 = !parens2;
  if(parens1) o << "(";
  parameters.front().print(o);
  if(parens1) o << ")";
  o << "^";
  if(parens2) o << "(";
  parameters.back().print(o);
  if(parens2) o << ")";
}


Simplified Power::simplify() const
{
 list<Symbolic>::iterator i, j;
 const Symbolic &b = parameters.front().simplify();
 const Symbolic &n = parameters.back().simplify();
 if(n == 0) return Number<int>(1);
 if(n == 1) return b;
 if(b == 0) return Number<int>(0);
 if(b == 1) return Number<int>(1);

 if(b.type() == typeid(Power))
 {
  CastPtr<const Power> p = b;
  return (p->parameters.front() ^ (p->parameters.back() * n)).simplify();
 }
 if(b.type() == typeid(Numeric) &&
    n.type() == typeid(Numeric) &&
    Number<void>(n).numerictype() == typeid(int))
 {
  int i = CastPtr<const Number<int> >(n)->n;
  int inv = (i<0) ? 1 : 0;
  Number<void> r = Number<int>(1), x = b;
  i = (inv) ? -i : i;
  while(i != 0)
  {
   if(i & 1)  r = r * x;
   x = x * x;
   i >>= 1;
  }
  if(inv) return (Number<int>(1) / r);
  return r;
 }
 if(n.type() == typeid(Log))
 {
  CastPtr<const Log> l(n);
  if(l->parameters.front() == b)
   return l->parameters.back();
 }
 if(b.type() == typeid(Numeric) &&
    n.type() == typeid(Numeric) &&
    Number<void>(b).numerictype() == typeid(double))
 {
  double nd, bd = CastPtr<const Number<double> >(b)->n;
  if(Number<void>(n).numerictype() == typeid(int))
   nd = CastPtr<const Number<int> >(n)->n;
  else if(Number<void>(n).numerictype() == typeid(double))
   nd = CastPtr<const Number<double> >(n)->n;
  else if(Number<void>(n).numerictype() == typeid(Rational<Number<void> >))
   nd = double(CastPtr<const Number<Rational<Number<void> > > >(n)->n);
  else return Power(b,n);
  if(bd >= 0.0 || int(nd) == nd)
   return Number<double>(pow(bd,nd));
 }
 if(b.type() == typeid(Numeric) &&
    n.type() == typeid(Numeric) &&
    Number<void>(n).numerictype() == typeid(double))
 {
  double bd, nd = CastPtr<const Number<double> >(n)->n;
  if(Number<void>(b).numerictype() == typeid(int))
   bd = CastPtr<const Number<int> >(b)->n;
  else if(Number<void>(b).numerictype() == typeid(double))
   bd = CastPtr<const Number<double> >(b)->n;
  else if(Number<void>(b).numerictype() == typeid(Rational<Number<void> >))
   bd = double(CastPtr<const Number<Rational<Number<void> > > >(b)->n);

  else return Power(b,n);
  if(bd >= 0.0 || int(nd) == nd)
   return Number<double>(pow(bd,nd));
  if(bd == -1 && floor(nd) == nd && nd != 0.5 && nd != -0.5) return Number<int>(1);
  else if (bd == -1 && floor(nd) != nd && nd != 0.5 && nd != -0.5) return Number<int>(-1);
 }
 return Power(b,n);
}

Expanded Power::expand() const
{
 Symbolic b = parameters.front().expand();
 Symbolic n = parameters.back().expand();
 if(b.type() != typeid(Sum)     &&
    b.type() != typeid(Product) &&
    b.type() != typeid(Numeric))
  return Power(b,n);

 // a^(b+c) == a^b a^c  when b and c commute
 if(n.type() == typeid(Sum))
 {
  CastPtr<const Sum> s(*n);
  Product p;
  list<Symbolic>::const_iterator k, k1;
  for(k=s->summands.begin();k!=s->summands.end();++k)
   for(++(k1=k);k1!=s->summands.end();++k1)
    if(!k->commute(*k1)) return Power(b,n);
  for(k=s->summands.begin();k!=s->summands.end();++k)
   p.factors.push_back(Power(b,*k));
  return p;
 }

 // (a*b)^c == a^c b^c  when a and b commute
 if(b.type() == typeid(Product))
 {
  CastPtr<const Product> p(*b);
  Product r;
  list<Symbolic>::const_iterator k, k1;
  for(k=p->factors.begin();k!=p->factors.end();++k)
   for(++(k1=k);k1!=p->factors.end();++k1)
    if(!k->commute(*k1)) return Power(b,n);
  for(k=p->factors.begin();k!=p->factors.end();++k)
   r.factors.push_back(Power(*k,n));
  return r;
 }

 if(n.type() != typeid(Numeric) ||
    Number<void>(n).numerictype() != typeid(int))
  return Power(b,n);

 int ae = Symbolic::auto_expand;
 int i = CastPtr<const Number<int> >(n)->n;
 if(i<0) return *this;
 int sgn = (i>=0) ? 1 : -1;
 Symbolic r = 1, x = b;
 i *= sgn;

 Symbolic::auto_expand = 0;
 while(i != 0)
 {
  if(i & 1) r = Product(r,x).expand();
  i >>= 1;
  if(i != 0) x = Product(x,x).expand();
 }
 Symbolic::auto_expand = ae;

 return Power(r,Number<int>(sgn));
}

Symbolic Power::subst(const Symbolic &x,const Symbolic &y,int &n) const
{
 if(x.type() == typeid(Power))
 {
  CastPtr<const Power> p(x);
  if(parameters.front() == p->parameters.front()                  &&
     parameters.back().type() == typeid(Numeric)                  &&
     p->parameters.back().type() == typeid(Numeric)               &&
     Number<void>(parameters.back()).numerictype() == typeid(int) &&
     Number<void>(p->parameters.back()).numerictype() == typeid(int))
  {
   int s = CastPtr<const Number<int> >(parameters.back())->n;
   int t = CastPtr<const Number<int> >(p->parameters.back())->n;
   if((s > 0) && (t > 0) && (s >= t))
   {
    ++n;
    return Power(parameters.front(),s - t).subst(x,y,n) * y;
   }
  }
 }
 return Symbol::subst(x,y,n);
}

// d/ds (a^b) = d/ds exp(b ln(a))
//            = (ln(a) db/ds + b/a da/ds) a^b
Symbolic Power::df(const Symbolic &s) const
{
 const Symbolic &a = parameters.front();
 const Symbolic &b = parameters.back();
 return (ln(a)*b.df(s) + (b/a)*a.df(s)) * *this;
}

// Pochammer function
Symbolic pochhammer(Symbolic x, int n) {
	if (n < 0) 
	{
		throw invalid_argument("n must be non-negative");
	}
	if (n == 0) 
	{
		return 1.0;
	}
	Symbolic result = 1.0;
	for (int i = 0; i < n; ++i) 
	{
		result *= (x + i);
	}
	return result;
}

// Hypergeometric 2F1 function
Symbolic hypergeometric_2F1_direct(Symbolic a, Symbolic b, Symbolic c, const Symbolic &s, int max_iterations = 5) {
	Symbolic sum = 1;
	Symbolic term = 1;
	for (int n = 1; n <= max_iterations; ++n) 
	{
		term *= (pochhammer(a, n) * pochhammer(b, n) * (s^(n))) / (pochhammer(c, n) * tgamma(n + 1));
		sum += term;
	}
	return sum;
}

// Factorial and combinations
Symbolic factorial(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	Symbolic result = 1;
	for (int i = 2; i <= n; ++i) 
	{
		result *= i;
	}
	return result;
}

Symbolic combinations(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factorial(n) / (factorial(r) * factorial(n - r));
}

Symbolic Power::integrate(const Symbolic &s) const
{
 const Symbolic &a = parameters.front();
 const Symbolic &b = parameters.back();
 int bpower = parameters.back().coeff(s,0);
 if(a == 1 - cos(s))
 {
	list<Equations> eq, eq2;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (1-cos(s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	if(b.df(s) == 0 && b==-1) 
	{
		return (-1)*(tan(0.5*s));
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return (-1)*(2*tan(0.5*s)) + (-1)*(tan(0.5*s)^(3));
	}
	} catch(const SymbolicError &se) {}	
 }
 if(a == (sin(s)^(-bpower+1))*(1-cos(s)) ) // for sin(x) * (1-cos(x))^(-b)
 {
	if(b.df(s) == 0 && b ==-1) 
	{
		return ln(cos(s)-1);
	}
	if(b.df(s) == 0 && b !=-1 && bpower < 0) 
	{
		return cos(s)/(bpower*((1-cos(s))^(bpower)) - ((1-cos(s))^(bpower))) - 1/(bpower*((1-cos(s))^(bpower)) - ((1-cos(s))^(bpower)));
	}
		
 }
 if(a.type() == typeid(Cos))
 {
	list<Equations> eq, eq2;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (cos(c*s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	Symbolic ap = rhs(*i, c);
	if(b.df(s) == 0 && b==2) 
	{
		return (2*ap*s + sin(2*ap*s)) / (4*ap);
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return sin(ap*s) / (ap*cos(ap*s)) ;
	}
	} catch(const SymbolicError &se) {}	
	
	eq2 = (cos(s)).match(a, (c,d));
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	try {
	if(b.df(s) == 0 && b==2) 
	{
		return 0.5*s + 0.5*(sin(s) * cos(s));
	}	
	if(b.df(s) == 0 && b==-1) 
	{
		return -0.5*ln(sin(s)-1) + 0.5*( ln(sin(s) + 1) );
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return sin(s) / cos(s);
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 != 0) // odd power case
	{
		Symbolic sgn = 1;
		Symbolic integral; 
		int j =1;
		int m = bpower-(0.5*(bpower+1));
		for(int i = 1 ; i <= bpower ; i = i+2)
		{
			integral += sgn*combinations(m,j-1)*((sin(s))^(i))/(i);
			sgn = -sgn;
			j = j+1;
		} 
		return integral;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 == 0) //  even power case
	{
		Symbolic c = 1;
		Symbolic d = parameters.back().coeff(s,0);
		Symbolic integral; 
		for(int i = bpower ; i >= 2 ; i = i-2)
		{
			integral += c*sin(s)*((cos(s))^(i-1))/(d) ;
			d = d*(i-2);
			c = c*(i-1);
			if (i==4)
			{
				integral += ( (c*s)/(d) );
			}
		} 
		return integral;
	}
	} catch(const SymbolicError &se) {}	
 }
 if(a.type() == typeid(Sin))
 {
	list<Equations> eq, eq2, eq3;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (sin(c*s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	Symbolic ap = rhs(*i, c);
	if(b.df(s) == 0 && b==2) 
	{
		return 0.5*s - ( sin(2*ap*s) / (4*ap) );
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return - cos(ap*s) / (ap*sin(ap*s)) ;
	}
	} catch(const SymbolicError &se) {}	
	eq2 = (sin(s)).match(a, (c,d));
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	try {
	if(b.df(s) == 0 && b==2) 
	{
		return 0.5*s - 0.5*(sin(s) * cos(s));
	}
	if(b.df(s) == 0 && b==-1) 
	{
		return 0.5*ln(cos(s)-1) - 0.5*( ln(cos(s) + 1) );
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return -cos(s) / sin(s);
	}
	// Too slow, we are not using hypergeometric function
	/*
	if(b.df(s) == 0 && b!=0) 
	{
		Symbolic y_sin = -(sin(s)^(b+1))*(cos(s))*((sin(s)^(2))^(-0.5*b-0.5));
		Symbolic a_sin = 0.5;
		Symbolic b_sin = 0.5*(1-b);
		Symbolic c_sin = 1.5;
		Symbolic z_sin = cos(s)*cos(s);
		return y_sin*hypergeometric_2F1_direct(a_sin, b_sin, c_sin, z_sin) ;
	}
	*/
	if(b.df(s) == 0 && b!=0 && bpower % 2 != 0) // odd power case
	{
		Symbolic sgn = -1;
		Symbolic integral; 
		int j =1;
		int m = bpower-(0.5*(bpower+1));
		for(int i = 1 ; i <= bpower ; i = i+2)
		{
			integral += sgn*combinations(m,j-1)*((cos(s))^(i))/(i);
			sgn = -sgn;
			j = j+1;
		} 
		return integral;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 == 0) //  even power case
	{
		Symbolic c = 1;
		Symbolic d = parameters.back().coeff(s,0);
		Symbolic integral; 
		for(int i = bpower ; i >= 2 ; i = i-2)
		{
			integral -= c*cos(s)*((sin(s))^(i-1))/(d) ;
			d = d*(i-2);
			c = c*(i-1);
			if (i==4)
			{
				integral += ( (c*s)/(d) );
			}
		} 
		return integral;
	}
	// NOT YET
	/*eq3 = ( (sin(s)^c) * (cos(s)^d) ).match(f, (c,d)); // case for cos(x)^n * sin(x)^m
	for(i=eq3.begin(); i!=eq3.end(); ++i)
	try {
	Symbolic m = rhs(*i, c), n = rhs(*i, d);
	if(df(rhs(*i, a), s) == 0)
	{
	if(ap == bp || bp == ap || bp == -ap || ap == -bp)
	{
	    return ;
	}
	return ;
	}*/

	} catch(const SymbolicError &se) {}	
 }
 if(a.type() == typeid(Tan))
 {
	list<Equations> eq, eq2;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (tan(c*s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	Symbolic ap = rhs(*i, c);
	if(b.df(s) == 0 && b==2) 
	{
		return (-ap*s + (sin(ap*s) / cos(ap*s))) / ap ;
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return (-ap*s - (cos(ap*s) / sin(ap*s))) / ap ;
	}
	} catch(const SymbolicError &se) {}	
	eq2 = (tan(s)).match(a, (c,d));
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	try {
	if(b.df(s) == 0 && b==2) 
	{
		return -s + (sin(s) / cos(s)) ;
	}
	if(b.df(s) == 0 && b==-1) 
	{
		return ln(sin(s));
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return -s - (cos(s) / sin(s)) ;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 != 0) //  odd power case
	{
		Symbolic integral;
		Symbolic integral_front;
		Symbolic sgn = 1;
		vector<int> v={1};
		vector<int> v_temp={1};
		vector<int> mc={1}; // to store the middle coefficient
		vector<int> mc_temp={1}; // to store the temporary / new middle coefficient
		int d0 = 2;
		int k = 1, l=2;
		// For the coefficient at the numerator 
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{
			if (i == 1)
			{
				v[0] = 1; 
				v.assign({v[0]});
			}	

			if (i ==2)
			{		
				mc[0]=1;
				mc.assign({mc[0]});
				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
							
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}
			}
			if (i == 3)
			{ 
				mc[0] = mc[0] + 1; // 
				mc.assign({mc[0]});	

				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}					
			}
			if (i == 4)
			{ 
				mc[0] = mc[0] + 1;	
				mc.assign({mc[0]});
				mc.push_back(mc[0]);
				
				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}					
			}
			if (i >= 5)
			{ 
				mc_temp[0] = mc[0]+1;
				mc_temp.assign({mc_temp[0]});
				for(int ic = 1 ; ic < l  ; ic = ic+1)
				{
					mc_temp[ic] = mc[ic-1] + mc[ic];
				}
				mc[l] = (mc_temp[0]);
				mc[l+1] = (mc_temp[0]);

				mc[0]=mc_temp[0];	
				mc.assign({mc[0]});

				for(int ic = 1 ; ic < l ; ic = ic+1)
				{
					mc.push_back(mc_temp[ic]);
				}
				mc.push_back(mc_temp[0]);
				mc.push_back(0);

				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}

				l = l+1;
			}

			d0 = d0*k;
			k = k+1;	
		} 	
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{			
			integral += sgn*v[i-1]*(cos(s)^(2*(i-1)));
			sgn=-sgn;
			
		} 	
		
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{	
			integral_front = sgn;			
			sgn = -sgn;
		} 	
		return integral_front*ln(cos(s)) + (integral)/(d0*(cos(s)^(bpower-1))) ;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 == 0) //  even power case
	{
		Symbolic integral;
		Symbolic integral_front;
		Symbolic sgn = -1;
		for(int i = 2 ; i <= (bpower)  ; i = i+2)
		{	
			integral = (-1)*integral;			
			integral += (sin(s)^(i-1)) / ((i-1)*(cos(s)^(i-1)));
		} 	
		
		for(int i = 1 ; i <= (bpower)/2  ; i = i+1)
		{	
			integral_front = sgn;			
			sgn = -sgn;
		} 	
		return integral_front*s + integral;
	}
	} catch(const SymbolicError &se) {}	
 }
 if(a.type() == typeid(Cot))
 {
	list<Equations> eq, eq2;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (cot(c*s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	Symbolic ap = rhs(*i, c);
	if(b.df(s) == 0 && b==2) 
	{
		return (-ap*s - (cos(ap*s) / sin(ap*s))) / ap ;
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return (-ap*s + (sin(ap*s) / cos(ap*s))) / ap ;
	}
	} catch(const SymbolicError &se) {}	
	eq2 = (cot(s)).match(a, (c,d));
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	try {
	if(b.df(s) == 0 && b==2) 
	{
		return -s - (cos(s) / sin(s)) ;
	}
	if(b.df(s) == 0 && b==-1) 
	{
		return -ln(cos(s));
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return -s + (sin(s) / cos(s)) ;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 != 0) //  odd power case
	{
		Symbolic integral;
		Symbolic integral_front;
		Symbolic sgn = -1;
		vector<int> v={1};
		vector<int> v_temp={1};
		vector<int> mc={1}; // to store the middle coefficient
		vector<int> mc_temp={1}; // to store the temporary / new middle coefficient
		int d0 = 2;
		int k = 1, l=2;
		// For the coefficient at the numerator 
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{
			if (i == 1)
			{
				v[0] = 1; 
				v.assign({v[0]});
			}	

			if (i ==2)
			{		
				mc[0]=1;
				mc.assign({mc[0]});
				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
							
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}
			}
			if (i == 3)
			{ 
				mc[0] = mc[0] + 1; // 
				mc.assign({mc[0]});	

				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}					
			}
			if (i == 4)
			{ 
				mc[0] = mc[0] + 1;	
				mc.assign({mc[0]});
				mc.push_back(mc[0]);
				
				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}					
			}
			if (i >= 5)
			{ 
				mc_temp[0] = mc[0]+1;
				mc_temp.assign({mc_temp[0]});
				for(int ic = 1 ; ic < l  ; ic = ic+1)
				{
					mc_temp[ic] = mc[ic-1] + mc[ic];
				}
				mc[l] = (mc_temp[0]);
				mc[l+1] = (mc_temp[0]);

				mc[0]=mc_temp[0];	
				mc.assign({mc[0]});

				for(int ic = 1 ; ic < l ; ic = ic+1)
				{
					mc.push_back(mc_temp[ic]);
				}
				mc.push_back(mc_temp[0]);
				mc.push_back(0);

				v_temp[0] = v[0]*(i-1); 
				v_temp.assign({v_temp[0]});
				for(int j = 1 ; j < i-1  ; j = j+1)
				{
					v_temp.push_back(v[j-1]*((2*i)-j) + v_temp[0]*mc[j-1] );	
				}
				v_temp.push_back((v[i-2]*(i+1) ) + (v_temp[0]) );
				// Assign the temporary vector to vector v for future use
				v[0] = v_temp[0];
				v.assign({v[0]});
				for(int j = 1 ; j < i  ; j = j+1)
				{
					v.push_back(v_temp[j]);	
				}

				l = l+1;
			}

			d0 = d0*k;
			k = k+1;	
		} 	
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{			
			integral += sgn*v[i-1]*(sin(s)^(2*(i-1)));
			sgn=-sgn;
			
		} 	
		
		for(int i = 1 ; i <= (bpower-1)/2  ; i = i+1)
		{	
			integral_front = sgn;			
			sgn = -sgn;
		} 	
		return integral_front*ln(sin(s)) + (integral)/(d0*(sin(s)^(bpower-1))) ;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 == 0) //  even power case
	{
		Symbolic integral;
		Symbolic integral_front;
		Symbolic sgn = -1;
		for(int i = 2 ; i <= (bpower)  ; i = i+2)
		{	
			integral = (-1)*integral;			
			integral += - (cos(s)^(i-1)) / ((i-1)*(sin(s)^(i-1)));
		} 	
		
		for(int i = 1 ; i <= (bpower)/2  ; i = i+1)
		{	
			integral_front = sgn;			
			sgn = -sgn;
		} 	
		return integral_front*s + integral;
	}
	} catch(const SymbolicError &se) {}	
 }
 if(a.type() == typeid(Sec))
 {
	list<Equations> eq, eq2;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (sec(c*s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	Symbolic ap = rhs(*i, c);
	if(b.df(s) == 0 && b==2) 
	{
		return sin(ap*s) / (ap*cos(ap*s)) ;
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return (0.5*ap*s + (0.5*( sin(ap*s) * cos(ap*s) ) ) ) / ap ;
	}
	} catch(const SymbolicError &se) {}	
	eq2 = (sec(s)).match(a, (c,d));
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	try {
	if(b.df(s) == 0 && b==2) 
	{
		return sin(s) / cos(s) ;
	}
	if(b.df(s) == 0 && b==-1) 
	{
		return sin(s);
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return 0.5*s + 0.5*(sin(s) * cos(s)) ;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 != 0) // odd power case
	{
		int k = 1, l=2;
		vector<int> v={1};
		vector<int> mc={1}; // to store the middle coefficient
		Symbolic sgn = -1;
		Symbolic integral_numerator, integral_denominator; 
		Symbolic d0 = 2, d1 = 2;
		int j =1;
		int m = bpower-(0.5*(bpower+1));
		int last_coeff;
		int first_coeff = 3;
		
		// For the coefficient at the numerator of sine with odd power
		for(int i = 1 ; i < (bpower-1)/2  ; i = i+1)
		{
			if (i >= 2)
			{ 
				for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					mc[ic-1] = v[ic-1] + v[ic];
				}
			}
			k = k*l;
			last_coeff = v[j-1]; 
			v[0] = v[0]*(first_coeff+2*(i-1)); 
			v.assign({v[0]});			
					
			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v.push_back(mc[ic-1]*(first_coeff+2*(i-1)));	
				}
			d1 = d1*(2*i);
			d0 = d0*(2+2*i);
			v.push_back(last_coeff*(first_coeff+2*(i-1))+k);	
			l = l+2;	
			j = j+1;	
		} 	
		int j_num = 0 ;
		for(int i = bpower-2 ; i >= 1 ; i = i-2)
		{
			integral_numerator += sgn*v[j_num]*((sin(s))^(i));
			sgn = -sgn;
			j_num = j_num+1;
		} 
		sgn = 1;
		j = 1;
		for(int i = bpower-1 ; i >= 0 ; i = i-2)
		{
			integral_denominator += sgn*d0*combinations(m,j-1)*((sin(s))^(i));
			sgn = -sgn;
			j = j+1;
		} 
		d1 = d1*(bpower-1);
		return (-v[0]*ln(sin(s)-1))/(d1) + (v[0]*ln(sin(s)+1))/(d1) + (integral_numerator)/(integral_denominator);
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 == 0) //  even power case
	{
		int v[999];
		v[0] = 1;		 
		Symbolic integral; 
		Symbolic d0 = 1;	
		int k = 1, l = 2, c=1;
		for(int i = 1 ; i < bpower  ; i = i+2) // to compute the denominator
		{
			d0 *= i;		
		}
		for(int i = 1 ; i < bpower  - 1; i = i+2)
		{
			c = v[0];
			int arrsec[999]; // make the size of the array as big as possible
			
			for(int j = 0 ; j < i-1  ; j = j+1)
			{	
				arrsec[j] = v[j]; 
			}			
			int d;
			for(int j = 1 ; j < i  ; j = j+1)
			{		
				d = arrsec[j-1];	
				v[j] = d*l;		
			}
			
			v[0] = c*k;
			v[1] = c*l;
			
			k= k + 2;
			l = l + 2;
		} 	
			
		int j_d = 1 ;
		for(int i = 1 ; i < (0.5*bpower)+1; i = i+1)
		{
			integral += ( v[i-1] * sin(s) ) / ( d0*(cos(s)^(bpower-j_d)) ) ;
			j_d = j_d+2;
		} 
		return integral;
	}
	} catch(const SymbolicError &se) {}	
 }
 if(a.type() == typeid(Csc))
 {
	list<Equations> eq, eq2;
	list<Equations>::iterator i;
	UniqueSymbol c, d;

	eq = (csc(c*s)).match(a, (c,d));
	for(i=eq.begin(); i!=eq.end(); ++i)
	try {
	Symbolic ap = rhs(*i, c);
	if(b.df(s) == 0 && b==2) 
	{
		return -cos(ap*s) / (ap*sin(ap*s)) ;
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return (0.5*ap*s - (0.5*( sin(ap*s) * cos(ap*s) ) ) ) / ap ;
	}
	} catch(const SymbolicError &se) {}	
	eq2 = (csc(s)).match(a, (c,d));
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	try {
	if(b.df(s) == 0 && b==2) 
	{
		return -cos(s) / sin(s) ;
	}
	if(b.df(s) == 0 && b==-1) 
	{
		return -cos(s);
	}
	if(b.df(s) == 0 && b==-2) 
	{
		return 0.5*s - 0.5*(sin(s) * cos(s)) ;
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 != 0) // odd power case
	{
		int k = 1, l=2;
		vector<int> v={1};
		vector<int> mc={1}; // to store the middle coefficient
		Symbolic sgn = -1;
		Symbolic integral_numerator, integral_denominator; 
		Symbolic d0 = 2, d1 = 2;
		int j =1;
		int m = bpower-(0.5*(bpower+1));
		int last_coeff;
		int first_coeff = 3;
		
		// For the coefficient at the numerator of cosine with odd power
		for(int i = 1 ; i < (bpower-1)/2  ; i = i+1)
		{
			if (i >= 2)
			{ 
				for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					mc[ic-1] = v[ic-1] + v[ic];
				}
			}
			k = k*l;
			last_coeff = v[j-1]; 
			v[0] = v[0]*(first_coeff+2*(i-1)); 
			v.assign({v[0]});			
					
			for(int ic = 1 ; ic < i  ; ic = ic+1)
				{
					v.push_back(mc[ic-1]*(first_coeff+2*(i-1)));	
				}
			d1 = d1*(2*i);
			d0 = d0*(2+2*i);
			v.push_back(last_coeff*(first_coeff+2*(i-1))+k);	
			l = l+2;	
			j = j+1;	
		} 	
		int j_num = 0 ;
		for(int i = bpower-2 ; i >= 1 ; i = i-2)
		{
			integral_numerator += sgn*v[j_num]*((cos(s))^(i));
			sgn = -sgn;
			j_num = j_num+1;
		} 
		sgn = 1;
		j = 1;
		for(int i = bpower-1 ; i >= 0 ; i = i-2)
		{
			integral_denominator += sgn*d0*combinations(m,j-1)*((cos(s))^(i));
			sgn = -sgn;
			j = j+1;
		} 
		d1 = d1*(bpower-1);
		return (v[0]*ln(cos(s)-1))/(d1) - (v[0]*ln(cos(s)+1))/(d1) - (integral_numerator)/(integral_denominator);
	}
	if(b.df(s) == 0 && b!=0 && bpower % 2 == 0) //  even power case
	{
		int v[999];
		v[0] = 1;		 
		Symbolic integral; 
		Symbolic d0 = 1;	
		int k = 1, l = 2, c=1;
		for(int i = 1 ; i < bpower  ; i = i+2) // to compute the denominator
		{
			d0 *= i;		
		}
		for(int i = 1 ; i < bpower  - 1; i = i+2)
		{
			c = v[0];
			int arrsec[999]; // make the size of the array as big as possible
			
			for(int j = 0 ; j < i-1  ; j = j+1)
			{	
				arrsec[j] = v[j]; 
			}			
			int d;
			for(int j = 1 ; j < i  ; j = j+1)
			{		
				d = arrsec[j-1];	
				v[j] = d*l;		
			}
			
			v[0] = c*k;
			v[1] = c*l;
			
			k= k + 2;
			l = l + 2;
		} 	
			
		int j_d = 1 ;
		for(int i = 1 ; i < (0.5*bpower)+1; i = i+1)
		{
			integral -= ( v[i-1] * cos(s) ) / ( d0*(sin(s)^(bpower-j_d)) ) ;
			j_d = j_d+2;
		} 
		return integral;
	}
	} catch(const SymbolicError &se) {}	
 }

 if(b == -1 && parameters.front().coeff(s,2) == 0) return ln(parameters.front()) * (1/ parameters.front().df(s));
 if(b == -1 && parameters.front().coeff(s,2) != 0 && parameters.front().coeff(s,1) == 0 && parameters.front().coeff(s,0) == 0 ) return -(1 / (s*parameters.front().coeff(s,2)) ) ;
 if(b == -1 && parameters.front().coeff(s,2) != 0)
 {
  double a1 = parameters.front().coeff(s,2);
  double b1 = parameters.front().coeff(s,1);
  double c1 = parameters.front().coeff(s,0);
  double D_inv = sqrt(-1/(4*a1*c1-b1*b1));
  
  return - D_inv * ln(s + (-4*a1*c1*D_inv + b1*b1*D_inv + b1)/(2*a1)) + D_inv * ln(s + (4*a1*c1*D_inv - b1*b1*D_inv + b1)/(2*a1)) ;
 } 
 if(a == s && b.df(s) == 0 )
 {
  if(b == -1) return ln(a);
  return (a^(b+1)) / (b+1);
 }
 if(a == SymbolicConstant::e && b == s)
  return *this;
 if(a == SymbolicConstant::e && parameters.back().coeff(s,1) != 0 && parameters.back().coeff(s,0) != 0)
  return exp(parameters.back()) * (1/ parameters.back().df(s)); 
 if(df(s) == 0) return *this * s;
 return (a^(b+1)) / (b+1) ; 
}

PatternMatches
Power::match_parts(const Symbolic &s, const list<Symbolic> &p) const
{
 PatternMatches l;
 if(s.type() == typeid(Power))
 {
  CastPtr<const Power> pw(s);
  l = pw->parameters.front().match(parameters.front(), p);
  if(!l.empty()                                                   &&
     parameters.back().type() == typeid(Numeric)                  &&
     pw->parameters.back().type() == typeid(Numeric)              &&
     Number<void>(parameters.back()).numerictype() == typeid(int) &&
     Number<void>(pw->parameters.back()).numerictype() == typeid(int))
  {
   int s = CastPtr<const Number<int> >(parameters.back())->n;
   int t = CastPtr<const Number<int> >(pw->parameters.back())->n;
   if((s < 0) || (t < 0) || (s < t))
    pattern_match_FALSE(l);
  }
 }

 if(parameters.back().type() == typeid(Numeric) &&
    Number<void>(parameters.back()).numerictype() == typeid(int))
 {
  int n = abs(CastPtr<const Number<int> >(parameters.back())->n);
  Product pr;
  while(n-->0) pr.factors.insert(pr.factors.end(), parameters.back());
  pattern_match_OR(l, pr.match_parts(s, p));
 }

 pattern_match_OR(l, Symbol::match_parts(s, p));
 return l;
}

//////////////////////////////////////
// Implementation of Derivative     //
//////////////////////////////////////

Derivative::Derivative(const Derivative &d) : Symbol(d) { }

Derivative::Derivative(const Symbolic &s1,const Symbolic &s2)
: Symbol("df")
{
 if(s1.type() == typeid(Derivative))
 {
  parameters = CastPtr<const Derivative>(s1)->parameters;
  parameters.push_back(s2);
 }
 else
 {
  parameters.push_back(s1);
  parameters.push_back(s2);
 }
}

Symbolic Derivative::subst(const Symbolic &x,const Symbolic &y,int &n) const
{
 if(*this == x) return y;

 list<Symbolic>::const_iterator i;
 list<Symbolic>::iterator j;

 if(x.type() == type() &&
    parameters.front() == CastPtr<const Derivative>(x)->parameters.front())
 {
  // make a copy of this
  CastPtr<Derivative> l(*this);
  CastPtr<const Derivative> d(x);
  i = d->parameters.begin();
  for(++i;i!=d->parameters.end();++i)
  {
   j = l->parameters.begin();
   for(++j;j!=l->parameters.end();++j)
    if(*j == *i) break;
   if(j == l->parameters.end()) break;
   l->parameters.erase(j);
  }
  if(i == d->parameters.end())
  {
   ++n;
   Symbolic newdf = y;
   j = l->parameters.begin();
   for(++j;j!=l->parameters.end();++j)
    newdf = ::df(newdf,*j);
   return newdf;
  }
 }

 i = parameters.begin();
 Symbolic dy = i->subst(x,y,n);

 if(dy == *i)
 {
  Derivative d(*this);
  for(j=d.parameters.begin(),++j;j!=d.parameters.end();++j)
   *j = j->subst(x,y,n);
  return d;
 }

 for(++i;i!=parameters.end();++i)
  dy = dy.df(i->subst(x,y,n));

 return dy;
}

Symbolic Derivative::df(const Symbolic &s) const
{
 list<Symbolic>::const_iterator i;

 if(parameters.front().type() == typeid(Symbol))
 {
  Symbolic result;
  CastPtr<const Symbol> sym(parameters.front());
  for(i=sym->parameters.begin();i!=sym->parameters.end();++i)
   result = result + Derivative(*this,*i) * i->df(s);

  return result;
 }

 if(parameters.front().df(s) != 0)
 {
  Derivative d(*this);
  d.parameters.push_back(s);
  return d;
 }

 return 0;
}

int Derivative::compare(const Symbolic &s) const
{
 list<Symbolic>::const_iterator i;
 list<Symbolic>::iterator j;

 if(s.type() != type()) return 0;
 // make a copy of s
 CastPtr<Derivative> d(*s);
 if(d->parameters.size() != parameters.size()) return 0;
 if(d->parameters.front() != parameters.front()) return 0;
 for(i=parameters.begin(),++i;i!=parameters.end();++i)
 {
  for(j=d->parameters.begin(),++j;j!=d->parameters.end();++j)
   if(*i == *j) break;
  if(j == d->parameters.end()) return 0;
  d->parameters.erase(j);
 }
 return 1;
}

Symbolic Derivative::integrate(const Symbolic &s) const
{
 int n = 0, n1;
 list<Symbolic>::const_iterator i, i1 = parameters.end();
 list<Symbolic>::iterator j;

 for(i=parameters.begin();i!=parameters.end();++i,++n)
  if(*i == s) { i1 = i; n1 = n; }

 if(i1 != parameters.end())
 {
  // make a copy of *this
  CastPtr<Derivative> d(*this);
  for(j=d->parameters.begin();n1!=0;++j,--n1);
  d->parameters.erase(j);
  if(d->parameters.size() == 1) return d->parameters.front();
  return *d;
 }

 if(parameters.front().df(s) == 0) return *this * s;
 return Integral(*this,s);
}

//////////////////////////////////////
// Implementation of Integral       //
//////////////////////////////////////

Integral::Integral(const Integral &d) : Symbol(d) { }

Integral::Integral(const Symbolic &s1,const Symbolic &s2) : Symbol("int")
{
 if(s1.type() == typeid(Integral))
 {
  parameters = CastPtr<const Integral>(s1)->parameters;
  parameters.push_back(s2);
 }
 else
 {
  parameters.push_back(s1);
  parameters.push_back(s2);
 }
}

Symbolic Integral::subst(const Symbolic &x,const Symbolic &y,int &n) const
{
 if(*this == x) { ++n; return y; }

 list<Symbolic>::const_iterator i = parameters.begin();
 Symbolic dy = i->subst(x,y,n);
 for(++i;i!=parameters.end();++i)
  dy = ::integrate(dy, i->subst(x,y,n));
 return dy;
}

Symbolic Integral::df(const Symbolic &s) const
{
 int n = 0, n1;
 list<Symbolic>::const_iterator i, i1 = parameters.end();
 list<Symbolic>::iterator j;

 for(i=parameters.begin();i!=parameters.end();++i,++n)
  if(*i == s) { i1 = i; n1 = n; }

 if(i1 != parameters.end())
 {
  // make a copy of *this
  CastPtr<Integral> in(*this);
  for(j=in->parameters.begin();n1!=0;++j,--n1);
  in->parameters.erase(j);
  if(in->parameters.size() == 1) return in->parameters.front();
  return *in;
 }

 if(parameters.front().df(s) == 0) return 0;
 return Derivative(*this,s);
}

Symbolic Integral::integrate(const Symbolic &s) const
{
 if(parameters.front().df(s) != 0)
 {
  Integral i(*this);
  i.parameters.push_back(s);
  return i;
 }

 return *this * s;
}

//////////////////////////////////////
// Implementation of Rows           //
//////////////////////////////////////

Rows::Rows(const Rows &s) : Symbol(s) {}

Rows::Rows(const Symbolic &s) : Symbol(Symbol("rows")[s]) {}

Simplified Rows::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Rows(s);
 return Symbolic(CastPtr<const SymbolicMatrix>(s)->rows());
}

//////////////////////////////////////
// Implementation of Columns        //
//////////////////////////////////////

Columns::Columns(const Columns &s) : Symbol(s) {}

Columns::Columns(const Symbolic &s) : Symbol(Symbol("columns")[s]) {}

Simplified Columns::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Columns(s);
 return Symbolic(CastPtr<const SymbolicMatrix>(s)->cols());
}

//////////////////////////////////////
// Implementation of Row            //
//////////////////////////////////////

Row::Row(const Row &s) : Symbol(s), row(s.row) {}

Row::Row(const Symbolic &s,int r) : Symbol(Symbol("row")[s,r]), row(r) {}

Simplified Row::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Row(s, row);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s)->operator[](row));
}

//////////////////////////////////////
// Implementation of Column         //
//////////////////////////////////////

Column::Column(const Column &s) : Symbol(s), column(s.column) {}

Column::Column(const Symbolic &s,int c)
 : Symbol(Symbol("column")[s,c]), column(c) {}

Simplified Column::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Column(s, column);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s)->operator()(column));
}

//////////////////////////////////////
// Implementation of Transpose      //
//////////////////////////////////////

Transpose::Transpose(const Transpose &s) : Symbol(s) {}

Transpose::Transpose(const Symbolic &s) : Symbol(Symbol("transpose")[s]) {}

Simplified Transpose::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Transpose(s);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s)->transpose());
}

//////////////////////////////////////
// Implementation of Trace          //
//////////////////////////////////////

Trace::Trace(const Trace &s) : Symbol(s) {}

Trace::Trace(const Symbolic &s) : Symbol(Symbol("tr")[s]) {}

Simplified Trace::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Trace(s);
 return CastPtr<const SymbolicMatrix>(s)->trace();
}

//////////////////////////////////////
// Implementation of Determinant    //
//////////////////////////////////////

Determinant::Determinant(const Determinant &s) : Symbol(s) {}

Determinant::Determinant(const Symbolic &s) : Symbol(Symbol("det")[s]) {}

Simplified Determinant::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Determinant(s);
 return CastPtr<const SymbolicMatrix>(s)->determinant();
}

//////////////////////////////////////
// Implementation of Vec            //
//////////////////////////////////////

Vec::Vec(const Vec &s) : Symbol(s) {}

Vec::Vec(const Symbolic &s) : Symbol(Symbol("vec")[s]) {}

Simplified Vec::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(SymbolicMatrix)) return Vec(s);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s)->vec());
}

//////////////////////////////////////
// Implementation of Kronecker      //
//////////////////////////////////////

Kronecker::Kronecker(const Kronecker &s) : Symbol(s) {}

Kronecker::Kronecker(const Symbolic &s1, const Symbolic &s2)
 : Symbol(Symbol("kron")[s1,s2]) {}

Simplified Kronecker::simplify() const
{
 const Symbolic &s1 = parameters.front().simplify();
 const Symbolic &s2 = parameters.back().simplify();
 if(s1.type() != typeid(SymbolicMatrix)) return Kronecker(s1,s2);
 if(s2.type() != typeid(SymbolicMatrix)) return Kronecker(s1,s2);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s1)
                       ->kron(*CastPtr<const SymbolicMatrix>(s2)));
}

//////////////////////////////////////
// Implementation of DirectSum      //
//////////////////////////////////////

DirectSum::DirectSum(const DirectSum &s) : Symbol(s) {}

DirectSum::DirectSum(const Symbolic &s1, const Symbolic &s2)
 : Symbol(Symbol("dsum")[s1,s2]) {}

Simplified DirectSum::simplify() const
{
 const Symbolic &s1 = parameters.front().simplify();
 const Symbolic &s2 = parameters.back().simplify();
 if(s1.type() != typeid(SymbolicMatrix)) return DirectSum(s1,s2);
 if(s2.type() != typeid(SymbolicMatrix)) return DirectSum(s1,s2);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s1)
                       ->dsum(*CastPtr<const SymbolicMatrix>(s2)));
}

//////////////////////////////////////
// Implementation of Hadamard       //
//////////////////////////////////////

Hadamard::Hadamard(const Hadamard &s) : Symbol(s) {}

Hadamard::Hadamard(const Symbolic &s1, const Symbolic &s2)
 : Symbol(Symbol("hadamard")[s1,s2]) {}

Simplified Hadamard::simplify() const
{
 const Symbolic &s1 = parameters.front().simplify();
 const Symbolic &s2 = parameters.back().simplify();
 if(s1.type() != typeid(SymbolicMatrix)) return Hadamard(s1,s2);
 if(s2.type() != typeid(SymbolicMatrix)) return Hadamard(s1,s2);
 return SymbolicMatrix(CastPtr<const SymbolicMatrix>(s1)
                        ->hadamard(*CastPtr<const SymbolicMatrix>(s2)));
}

//////////////////////////////////////
// Implementation of Gamma          //
//////////////////////////////////////

Gamma::Gamma(const Gamma &s) : Symbol(s) {}

Gamma::Gamma(const Symbolic &s) : Symbol(Symbol("gamma")[s]) {}

Simplified Gamma::simplify() const
{
 const Symbolic &s = parameters.front().simplify();
 if(s.type() != typeid(Numeric)) return Gamma(s);
 CastPtr<Numeric> n(s);
 if(n->numerictype() == typeid(int))
 {
  int i = int(s) - 1;
  if(i >= 0)
  {
   Symbolic f = 1;
   while(i > 0) f *= i--;
   return f;
  }
 }
 return Gamma(s);
}

#endif
#endif

