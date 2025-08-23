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

#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS
#ifndef SYMINTEGRATION_CPLUSPLUS_STATISTICS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_STATISTICS_DECLARE

Symbolic binomialpmf(int, int, double);
Symbolic binomialcdf(int, int, double);
Symbolic binomialmean(int, int, double);
Symbolic binomialvar(int, int, double);
Symbolic binomialmgf(int, int, double);

Symbolic negativebinomialpmf(int, int, double);
Symbolic negativebinomialcdf(int, int, double);
Symbolic negativebinomialmean(int, int, double);
Symbolic negativebinomialvar(int, int, double);
Symbolic negativebinomialmgf(const Symbolic &, const Symbolic &, const Symbolic &);

Symbolic geometricpmf(int, double);
Symbolic geometricmean(int, double);
Symbolic geometricvar(int, double);
Symbolic geometricmgf(const Symbolic &, const Symbolic &);

Symbolic poissonpmf(int, int);
Symbolic poissoncdf(int, int);
Symbolic poissonmean(int, int);
Symbolic poissonvar(int, int);
Symbolic poissonmgf(const Symbolic &, const Symbolic &);

#endif
#endif


#endif
