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

#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DECLARE

#include "polynomial.h"
#include "rational.h"

double roundToDecimal(double, int); 

Symbolic dsolve(const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic dsolve(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic dsolvelogistic(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
Symbolic dsolveseparable(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);

void secondorderlineardiffeq_dsolve(double, double, double, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_ivpsolution(double, double, double, const Symbolic &, const Symbolic &, double, double, double);
void secondorderlineardiffeq_springmasssystem(double, double, double, const Symbolic &, const Symbolic &, double, double);
void secondorderlineardiffeq_RLCserieselectriccircuit(double, double, double, double, double);

void wronskian(double, double, double, const Symbolic &, const Symbolic &, double);
Symbolic wronskian_resultonly(double, double, double, const Symbolic &, const Symbolic &);
void wronskian_fundamentalsetofsolutions(double, double, double, const Symbolic &, const Symbolic &);

void reductionoforder(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &); 

void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationsivpsolution(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double, double, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationsivpforcedvibrationssolution(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double, double, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &, const Symbolic &, const Symbolic &, Polynomial<double> &);
void secondorderlineardiffeq_nonhomogeneousequationssolution(const Symbolic &, const Symbolic &, const Symbolic &, const SymbolicMatrix &, const Symbolic &, const Symbolic &);
void secondorderlineardiffeq_nonhomogeneousequationssolution_variationofparameters(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);

vector<complex<double>> higherorderlineardiffeq_vectorize(const Symbolic &,const Symbolic &,const Symbolic &, int);

void higherorderlineardiffeq_homogeneousequationsivpsolution(const vector<complex<double>> &,  const vector<complex<double>> &, const vector<complex<double>> &,  int);
void higherorderlineardiffeq_homogeneousequationsivpsolution(const vector<complex<double>> &,  const vector<complex<double>> &);
void higherorderlineardiffeq_twospringtwomasssystem(double, double, double, double, Symbolic &, const vector<complex<double>> &);
void higherorderlineardiffeq_nonhomogeneousequationsgeneralsolution(const vector<complex<double>> &, const Symbolic &, const Symbolic &);
void higherorderlineardiffeq_nonhomogeneousequationsgeneralsolution(const vector<complex<double>> &, const SymbolicMatrix &, const Symbolic &);
void higherorderlineardiffeq_nonhomogeneousequations_undeterminedcoefficients(const Symbolic &, const Symbolic &, const Symbolic &);
void higherorderlineardiffeq_nonhomogeneousequations_variationofparameters(const vector<complex<double>> &, const Symbolic &, const Symbolic &);

Symbolic ivp(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);

#endif
#endif


#endif
