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


#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHODS_DECLARE

double divisiond(double, double);

Symbolic bisectionmethod(const Symbolic &, const Symbolic &,double, double, int);
Symbolic newtonmethod(const Symbolic &, const Symbolic &, const Symbolic &, int);
Symbolic eulermethod(const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, const Symbolic &, double);
Symbolic directionfield(const Symbolic &, const Symbolic &, const Symbolic &, double, double, double, double, double, double);

double NewtonRaphson(const Symbolic &, const Symbolic &, double);
double secantmethod(const Symbolic &, const Symbolic &, double,  double, double, int);

void AberthEhrlich(const vector<complex<double>> &,  const vector<complex<double>> &,  int);

double numericaldifferentiation(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation3pointoneside(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation3pointbothsides(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation5pointoneside(const Symbolic &, const Symbolic &, double, double);
double numericaldifferentiation5pointbothsides(const Symbolic &, const Symbolic &, double, double);

double trapezoidalrule(const Symbolic &, const Symbolic &, double, double);
double simpsonsrule(const Symbolic &, const Symbolic &, double, double);
double simpsonsrule38(const Symbolic &, const Symbolic &, double, double);

double richardsonextrapolation(const Symbolic &, const Symbolic &, double, double, int, int);

void gradientdescent(const Symbolic &, const Symbolic &, double, double, int);

void choleskyDecomposition(vector<vector<double>>);
void LUDecomposition(vector<vector<double>>,vector<vector<double>>,vector<vector<double>>);

void JacobiMethod(const vector<vector<double>>& , const vector<double>& , int, double); 
bool GaussSeidel(const vector<vector<double>>&, const vector<double>& , vector<double>&, double, int); 
bool SORIterativeMethod(const vector<vector<double>>&, const vector<double>&, vector<double>&,  double, double, int ); 

/*

	Iterative Techniques in Matrix Algebra

*/
#ifndef STRUCT_CRS_MATRIX
#define STRUCT_CRS_MATRIX

// Defining struct  in .h file (header file) will make it work when we call it from .cpp file from anywhere that use SymIntegration library
// Structure to store a sparse matrix in Compressed Row Storage (CRS) format
struct CRSMatrix {
	vector<double> values;     // Non-zero elements
	vector<int> col_indices;   // Column indices of non-zero elements
	vector<int> row_ptr;       // Row pointers
	int num_rows;
	int num_cols;
};

CRSMatrix denseToCRS(const vector<vector<double>>& dense);
vector<double> spmv(const CRSMatrix& A, const std::vector<double>& x);
double dotProduct(const vector<double>& u, const vector<double>& v);
vector<double> get_Jacobi_preconditioner(const CRSMatrix& A);
vector<double> ConjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iterations);
vector<double> PreconditionedConjugateGradient(const CRSMatrix& A, const vector<double>& b, double tolerance, int max_iterations);
#endif

#endif
#endif


#endif
