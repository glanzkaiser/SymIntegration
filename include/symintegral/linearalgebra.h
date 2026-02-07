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

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

vector<vector<double>> loadMatrixFromFile(const string&);
vector<vector<complex<double>>> loadComplexMatrixFromFile(const string&);
vector<complex<double>> loadComplexVectorFromFile(const string&); 
vector<double> loadVectorFromFile(const string&);
void printMatrix(vector<vector<double>>);
void printVector(vector<double>);
vector<double> createVector(int, double); 
vector<vector<double>> createMatrix(int, int, double);
vector<vector<double>> createIdentityMatrix(int);
vector<vector<double>> createMatrixFromColumnVectors(const vector<vector<double>> &);
vector<vector<double>> VandermondeMatrix(const vector<double> &, int);
vector<double> getColumn(vector<vector<double>>, int); 
vector<double> getRow(vector<vector<double>>, int); 
vector<vector<double>> addColumn(vector<vector<double>>, vector<double>, int);
vector<vector<double>> addRow(vector<vector<double>>, vector<double>, int);
vector<vector<double>> deleteColumn(vector<vector<double>> &, int);
vector<vector<double>> deleteRow(vector<vector<double>> &, int);
int MaxElementIndex(vector<double>);

void printComplexMatrix(const ComplexMatrix &);
void printComplexVector(const ComplexVector &);
vector<vector<double>> ComplextoRealMatrix(vector<vector<complex<double>>> &);
complex<double> complexdivision(complex<double>, double);
vector<complex<double>> complexvecrand_normal(double, double, int);
vector<complex<double>> complexvecrand_normal_zeroimaginary(double, double, int);
vector<complex<double>> complexvecrand_uniform(double, double, int);
double complexnorm(const ComplexVector &);
double moduluscomplex(complex<double>);
complex<double> complexdotproduct (const ComplexVector &, const ComplexVector &);
complex<double> conjugate(complex<double>);
vector<complex<double>> addComplexVectors(const ComplexVector &, const ComplexVector &);
vector<complex<double>> subtractComplexVectors(const ComplexVector &, const ComplexVector &);
vector<complex<double>> scalarmultiplicationComplexVector(const ComplexVector &, double);
vector<complex<double>> complexnumbermultiplicationComplexVector(const vector<complex<double>>&, complex<double> );
vector<complex<double>> getcolumnComplexMatrix(const vector<vector<complex<double>>>& , int);

ComplexMatrix addComplexMatrices(const ComplexMatrix &, const ComplexMatrix &); 
ComplexMatrix subtractComplexMatrices(const ComplexMatrix &, const ComplexMatrix &); 
ComplexMatrix multiplyComplexMatrices(const ComplexMatrix &, const ComplexMatrix &);
ComplexMatrix scalarmultiplicationComplexMatrix(const ComplexMatrix &, double);
ComplexMatrix complexnumbermultiplicationComplexMatrix(const ComplexMatrix &, complex<double> );
ComplexVector multiplycomplexmatrixvector(const ComplexMatrix &, const ComplexVector &); 
ComplexMatrix createIdentityComplexMatrix(int);
ComplexMatrix createSubMatrix(vector<vector<complex<double>>> &, int, int, int, int);
ComplexMatrix TransposeComplexMatrix(ComplexMatrix &);
ComplexMatrix ComplexMatrixInverse(ComplexMatrix &);
complex<double> complexquadraticmultiplication(vector<vector<complex<double>>> &, vector<complex<double>> &);

vector<vector<double>> PermutationMatrixMax(vector<double> &);
vector<double> multiplymatrixvector(vector<vector<double>>  &, vector<double> &);
vector<vector<double>> multiply(vector<vector<double>>  &, vector<vector<double>> &);
vector<vector<double>> matrixsubtraction(vector<vector<double>>  &, vector<vector<double>> &);
vector<vector<double>> subtract(vector<vector<double>>  &, vector<vector<double>> &);
vector<vector<double>> add(vector<vector<double>>  &, vector<vector<double>> &);
vector<vector<double>> transpose(const vector<vector<double>> &);
void scalarmultiplication_alt(vector<vector<double>> &, double);
vector<vector<double>> scalarmultiplication(vector<vector<double>>, double);
double quadraticmultiplication(vector<vector<double>>  &, vector<double> &);

void spanningtest(vector<vector<double>> &);
void basistest(vector<vector<double>> &);
void linearindependencetest(vector<vector<double>> &);
void coordinatevector(vector<vector<double>> &, vector<vector<double>> &);
void basistransition(vector<vector<double>> &, vector<vector<double>> &);
void basistransition_withcoordinatevector(vector<vector<double>> &, vector<vector<double>> &, vector<double> &);
void rowspacebasis(vector<vector<double>> &);
void columnspacebasis(vector<vector<double>> &);
void homogeneouslinearsystembasis(vector<vector<double>> &);
vector<vector<double>> normalizehomogeneouslinearsystembasis_matrix(vector<vector<double>> &); 

vector<double> reflection_yaxis(vector<double> &); 
vector<double> reflection_xaxis(vector<double> &); 
vector<double> reflection_linex(vector<double> &); 
vector<double> reflection_xyplane(vector<double> &); 
vector<double> reflection_xzplane(vector<double> &); 
vector<double> reflection_yzplane(vector<double> &); 
vector<double> rotation_ccw(vector<double> &, double); 
vector<double> contractiondilation(vector<double> &, double); 
vector<double> compressionexpansion_xdirection(vector<double> &, double); 
vector<double> compressionexpansion_ydirection(vector<double> &, double); 
vector<double> shear_xdirection(vector<double> &, double); 
vector<double> shear_ydirection(vector<double> &, double); 
vector<double> orthogonalprojection_xaxis(vector<double> &); 
vector<double> orthogonalprojection_yaxis(vector<double> &); 
vector<double> orthogonalprojection_xyplane(vector<double> &); 
vector<double> orthogonalprojection_xzplane(vector<double> &); 
vector<double> orthogonalprojection_yzplane(vector<double> &); 
vector<double> rotation3d_xaxis_ccw(vector<double> &, double); 
vector<double> rotation3d_yaxis_ccw(vector<double> &, double); 
vector<double> rotation3d_zaxis_ccw(vector<double> &, double); 
vector<double> compressionexpansion3d_xdirection(vector<double> &, double); 
vector<double> compressionexpansion3d_ydirection(vector<double> &, double); 
vector<double> compressionexpansion3d_zdirection(vector<double> &, double); 

vector<double> MarkovChain(vector<vector<double>> &, vector<double> &, int);

vector<vector<double>> penrose(vector<double> &, vector<double> &);
vector<vector<complex<double>>> penroseComplex(vector<complex<double>> &, vector<complex<double>> &);

vector<vector<double>> gramschmidt(vector<vector<double>> &);
void QRDecomposition(vector<vector<double>> &, vector<vector<double>> &, vector<vector<double>> &); 
void QRDecompositionComplex(vector<vector<complex<double>>> &, vector<vector<complex<double>>> &, vector<vector<complex<double>>> &); 
void QRDecompositionComplexHouseholder(vector<vector<complex<double>>> &, vector<vector<complex<double>>> &, vector<vector<complex<double>>> &); 

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

void gaussianelimination(const vector<vector<double>> &);
Symbolic gaussianelimination(const SymbolicMatrix&, int);
Symbolic gaussianeliminationtest(const SymbolicMatrix&, int, int);
void GaussJordanComplexMatrix(ComplexMatrix &);
void GaussJordanComplexMatrixTEST(ComplexMatrix &);// this can be used to compute complex eigenvectors

void solve_nhsystem(vector<vector<double>> &,vector<vector<double>> &, vector<double> &);
void solve_nhsystem_resultsonly(vector<vector<double>> &,vector<vector<double>> &, vector<double> &);
void LUsolve_nhsystem(vector<vector<double>> &,vector<vector<double>> &, vector<double> &);

SymbolicMatrix IVPSolution_firstorderdiffeq(vector<vector<complex<double>>> &, vector<complex<double>> &, double, double);

vector<double> PowerMethod(vector<vector<double>> &, int);
vector<double> EigenvaluesEigenvectorsApproximation(vector<vector<double>> &, int);
vector<double> RayleighQuotientIteration(vector<vector<double>>&, int, double, double);
vector<double> diagonalization(vector<vector<double>> &);
vector<double> SVD(vector<vector<double>> &);

vector<vector<complex<double>>> QRDecompositionHouseholderReflections(vector<vector<complex<double>>> &);
vector<vector<complex<double>>> HessenbergDecomposition(vector<vector<complex<double>>> &);
vector<vector<complex<double>>> HessenbergDecompositionResultOnly(vector<vector<complex<double>>> &);
vector<complex<double>> QRAlgorithmwithShifts(vector<vector<complex<double>>>& , int);
vector<complex<double>> FrancisDoubleStepQRReal(vector<vector<complex<double>>>& , int);
vector<complex<double>> FrancisDoubleStepQR(vector<vector<complex<double>>>& , int);
vector<complex<double>> ComplexRayleighQuotientIteration(vector<vector<complex<double>>> & , int, int, double , double);
vector<vector<complex<double>>> SpectralDecomposition(vector<vector<complex<double>>> &, double, double);

#endif
#endif


#endif
