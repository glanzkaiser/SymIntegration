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

int fibonacciseries(int);
double rising_pochhammer(double, int);
double falling_pochhammer(double, int);

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

Symbolic hypergeometricpmf(int, int, int, int);
Symbolic hypergeometricmean(int, int, int, int);
Symbolic hypergeometricvar(int, int, int, int);

Symbolic poissonpmf(int, int);
Symbolic poissoncdf(int, int);
Symbolic poissonmean(int, int);
Symbolic poissonvar(int, int);
Symbolic poissonmgf(const Symbolic &, const Symbolic &);

Symbolic uniformpdf(double, double, double);
Symbolic uniformcdf(double, double, double);
Symbolic uniformmgf(double, double, double);
Symbolic uniformmean(double, double, double);
Symbolic uniformvar(double, double, double);

Symbolic normalpdf(double, double, double);
Symbolic normalcdf(double, double, double);
Symbolic normalmgf(double, double, double);
Symbolic normalmean(double, double, double);
Symbolic normalvar(double, double, double);

Symbolic gammapdf(double, double, double);
Symbolic gammacdf(double, double, double);
Symbolic gammamgf(double, double, double);
Symbolic gammamean(double, double, double);
Symbolic gammavar(double, double, double);

Symbolic exponentialpdf(double, double);
Symbolic exponentialcdf(double, double);
Symbolic exponentialmgf(double, double);
Symbolic exponentialmean(double, double);
Symbolic exponentialvar(double, double);

Symbolic betapdf(double, double, double);
Symbolic betacdf(double, double, double);
Symbolic betamgf(double, double, double);
Symbolic betamean(double, double, double);
Symbolic betavar(double, double, double);

Symbolic chisquaredpdf(double, double);
Symbolic chisquaredcdf(double, double);
Symbolic chisquaredmgf(double, double);
Symbolic chisquaredmean(double, double);
Symbolic chisquaredvar(double, double);

Symbolic cauchypdf(double);

Symbolic Fpdf(double, double, double);
Symbolic Fcdf(double, double, double);
Symbolic Fmean(double, double, double);
Symbolic Fvar(double, double, double);

Symbolic tpdf(double, double);
Symbolic tcdf(double, double);
Symbolic tmean(double, double);
Symbolic tvar(double, double);

Symbolic laplacepdf(double, double);
Symbolic laplacecdf(double, double);
Symbolic laplacemgf(double, double);
Symbolic laplacemean(double, double);
Symbolic laplacevar(double, double);

Symbolic logisticpdf(double, double);
Symbolic logisticcdf(double, double);
Symbolic logisticmgf(double, double);
Symbolic logisticmean(double, double);
Symbolic logisticvar(double, double);

double descriptivestatistics(vector<double>);

double rpearson(const SymbolicMatrix&, int);
Symbolic regressionline(const SymbolicMatrix&, int);

int randomnumberint(int, int, int);
double randomnumberreal(double, double, int);

#endif
#endif


#endif
