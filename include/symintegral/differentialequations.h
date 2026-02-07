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

#ifndef SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS
#ifndef SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS_DECLARE

Symbolic picarditeration(const Symbolic &, double, Symbolic &, Symbolic &, int);
Symbolic picarditeration_mathematicalinduction(const Symbolic &, double, Symbolic &, Symbolic &);

void adjointequation(const Symbolic &, const Symbolic &, const Symbolic &, Symbolic &);

#endif
#endif


#endif
