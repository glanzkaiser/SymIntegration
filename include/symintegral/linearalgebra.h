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

#include <random>
#include <iostream>
#include <vector>

#ifndef SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA
#ifndef SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DECLARE

vector<vector<double>> loadMatrixFromFile(const string&);
vector<double> loadVectorFromFile(const string&);
void printMatrix(vector<vector<double>>);
void printVector(vector<double>);
vector<double> getColumn(vector<vector<double>>, int); 
vector<double> getRow(vector<vector<double>>, int); 

vector<vector<double>> multiply(vector<vector<double>>  &, vector<vector<double>> &);
vector<vector<double>> add(vector<vector<double>>  &, vector<vector<double>> &);
vector<vector<double>> transpose(const vector<vector<double>> &);
void scalarmultiplication_alt(vector<vector<double>> &, double);
vector<vector<double>> scalarmultiplication(vector<vector<double>>, double);

void getCofactor(vector<vector<double>> &, vector<vector<double>> &, int, int, int);
double determinant1(vector<vector<double>> &, int);
vector<vector<double>> adjugate(vector<vector<double>> &);
vector<vector<double>> inverse(vector<vector<double>> &);

double determinant_alt(vector<vector<double>>);
long double determinant(const vector<vector<double>> &);
double norm(const vector<double> &);
double dot(const vector<double> &, const vector<double> &);
double distance(const vector<double> &, const vector<double> &);
double distance(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const vector<double> &);

double angle(const vector<double> &, const vector<double> &);
double degtorad(double);
double radtodeg(double);
vector<double> orthogonalprojection(vector<double>, vector<double>);
vector<double> add(vector<double> &, vector<double> &);
vector<double> subtract(vector<double> &, vector<double> &);
vector<double> scalarmultiplication(vector<double> &, double);
vector<double> crossproduct(vector<double> &, vector<double> &);
double scalartripleproduct(vector<double> &, vector<double> &, vector<double> &);

SymbolicMatrix vectorequation(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &);
void vectorequationdecomp(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, vector<double> &, vector<double> &, vector<double> &);

Symbolic gaussianelimination(const SymbolicMatrix&, int);
Symbolic gaussianeliminationtest(const SymbolicMatrix&, int, int);

void solve_nhsystem(vector<vector<double>> &,vector<vector<double>> &, vector<double> &);
void LUsolve_nhsystem(vector<vector<double>> &,vector<vector<double>> &, vector<double> &);

#endif
#endif


#endif
