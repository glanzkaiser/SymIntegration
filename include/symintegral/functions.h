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


// functions.h

#ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS

#include <cmath>
using namespace std;

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS_FORWARD
#define SYMBOLIC_CPLUSPLUS_FUNCTIONS_FORWARD

class Sin;
class Asin;
class Cos;
class Acos;
class Tan;
class Atan;
class Cot;
class Acot;
class Sec;
class Asec;
class Csc;
class Acsc;
class Sinh;
class Asinh;
class Cosh;
class Acosh;
class Tanh;
class Coth;
class Sech;
class Csch;
class Log;
class Erf;
class Power;
class Derivative;
class Integral;
class Rows;
class Columns;
class Row;
class Column;
class Transpose;
class Trace;
class Determinant;
class Vec;
class Kronecker;
class DirectSum;
class Hadamard;
class Gamma;

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMBOLIC_CPLUSPLUS_FUNCTIONS
#ifndef SYMBOLIC_CPLUSPLUS_FUNCTIONS_DECLARE
#define SYMBOLIC_CPLUSPLUS_FUNCTIONS_DECLARE

class Sin: public Symbol
{
 public: Sin(const Sin&);
         Sin(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Asin: public Symbol
{
 public: Asin(const Asin&);
         Asin(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Cos: public Symbol
{
 public: Cos(const Cos&);
         Cos(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Acos: public Symbol
{
 public: Acos(const Acos&);
         Acos(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Tan: public Symbol
{
 public: Tan(const Tan&);
         Tan(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Atan: public Symbol
{
 public: Atan(const Atan&);
         Atan(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Cot: public Symbol
{
 public: Cot(const Cot&);
         Cot(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Acot: public Symbol
{
 public: Acot(const Acot&);
         Acot(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Sec: public Symbol
{
 public: Sec(const Sec&);
         Sec(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Asec: public Symbol
{
 public: Asec(const Asec&);
         Asec(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Csc: public Symbol
{
 public: Csc(const Csc&);
         Csc(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};


class Acsc: public Symbol
{
 public: Acsc(const Acsc&);
         Acsc(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Sinh: public Symbol
{
 public: Sinh(const Sinh&);
         Sinh(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Asinh: public Symbol
{
 public: Asinh(const Asinh&);
         Asinh(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Cosh: public Symbol
{
 public: Cosh(const Cosh&);
         Cosh(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Acosh: public Symbol
{
 public: Acosh(const Acosh&);
         Acosh(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Tanh: public Symbol
{
 public: Tanh(const Tanh&);
         Tanh(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Coth: public Symbol
{
 public: Coth(const Coth&);
         Coth(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Sech: public Symbol
{
 public: Sech(const Sech&);
         Sech(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Csch: public Symbol
{
 public: Csch(const Csch&);
         Csch(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Erf: public Symbol
{
 public: Erf(const Erf&);
         Erf(const Symbolic&);

         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Log: public Symbol
{
 public: Log(const Log&);
         Log(const Symbolic&, const Symbolic&);

         void print(ostream&) const;
         Simplified simplify() const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Power: public Symbol
{
 public: Power(const Power&);
         Power(const Symbolic&,const Symbolic&);

         void print(ostream&) const;
         Simplified simplify() const;
         Expanded expand() const;
         Symbolic subst(const Symbolic &x,const Symbolic &y,int &n) const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;
         PatternMatches match_parts(const Symbolic &s, const list<Symbolic>&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Derivative: public Symbol
{
 public: Derivative(const Derivative&);
         Derivative(const Symbolic&,const Symbolic&);

         Symbolic subst(const Symbolic &x,const Symbolic &y,int &n) const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;
         int compare(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Integral: public Symbol
{
 public: Integral(const Integral&);
         Integral(const Symbolic&,const Symbolic&);

         Symbolic subst(const Symbolic &x,const Symbolic &y,int &n) const;
         Symbolic df(const Symbolic&) const;
         Symbolic integrate(const Symbolic&) const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Rows: public Symbol
{
 public: Rows(const Rows&);
         Rows(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Columns: public Symbol
{
 public: Columns(const Columns&);
         Columns(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Row: public Symbol
{
 private: int row;
 public:  Row(const Row&);
          Row(const Symbolic&,int);

          Simplified simplify() const;

          Cloning *clone() const { return Cloning::clone(*this); }
};

class Column: public Symbol
{
 private: int column;
 public:  Column(const Column&);
          Column(const Symbolic&,int);

          Simplified simplify() const;

          Cloning *clone() const { return Cloning::clone(*this); }
};

class Transpose: public Symbol
{
 public: Transpose(const Transpose&);
         Transpose(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Trace: public Symbol
{
 public: Trace(const Trace&);
         Trace(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Determinant: public Symbol
{
 public: Determinant(const Determinant&);
         Determinant(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Vec: public Symbol
{
 public: Vec(const Vec&);
         Vec(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Kronecker: public Symbol
{
 public: Kronecker(const Kronecker&);
         Kronecker(const Symbolic&,const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class DirectSum: public Symbol
{
 public: DirectSum(const DirectSum&);
         DirectSum(const Symbolic&,const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Hadamard: public Symbol
{
 public: Hadamard(const Hadamard&);
         Hadamard(const Symbolic&,const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

class Gamma: public Symbol
{
 public: Gamma(const Gamma&);
         Gamma(const Symbolic&);

         Simplified simplify() const;

         Cloning *clone() const { return Cloning::clone(*this); }
};

#endif
#endif


#endif
