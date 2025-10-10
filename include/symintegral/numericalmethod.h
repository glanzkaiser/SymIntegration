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


// dsolve.h

#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_DECLARE

Symbolic bisectionmethod(const Symbolic &, const Symbolic &,double, double, int);
Symbolic newtonmethod(const Symbolic &, const Symbolic &, const Symbolic &, int);
Symbolic eulermethod(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double);
Symbolic directionfield(const Symbolic &, const Symbolic &, const Symbolic &, double, double, double, double, double, double);

double numericaldifferentiation(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation3pointoneside(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation3pointbothsides(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation5pointoneside(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation5pointbothsides(const Symbolic &, const Symbolic &, double, double);

double trapezoidalrule(const Symbolic &, const Symbolic &, double, double);
double simpsonsrule(const Symbolic &, const Symbolic &, double, double);
double simpsonsrule38(const Symbolic &, const Symbolic &, double, double);

double richardsonextrapolation(const Symbolic &, const Symbolic &, double, double, int, int);
#endif
#endif


#endif
