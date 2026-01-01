/*
    SymIntegration is branching from SymbolicC++ 3.35
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
// THANKS SENTINEL!!! and Freya too

// g++ -o result main.cpp -lsymintegration

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

#include <complex>
#include <iomanip> // for std::setw

#include <iostream>
#include <vector>
#include <cmath> // For std::abs, std::atan2, M_PI

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	//try 
	//{
	// Define two 2x2 complex matrices
	ComplexMatrix A = 
	{
		{Complex(1, 2), Complex(3, -1)},
		{Complex(0, -1), Complex(4, 5)}
	};

	ComplexMatrix B = 
	{
		{Complex(2, 0), Complex(1, 3)},
		{Complex(-1, 4), Complex(0, -2)}
	};

	ComplexVector x = 
	{
		{Complex(1, 1)},
		{Complex(0, 1)},
		{Complex(3, -1)}
	};

	ComplexVector y = 
	{
		{Complex(1, 1)},
		{Complex(2, 0)},
		{Complex(0, 4)}
	};

	cout << "Matrix A:\n";
	printComplexMatrix(A);

	cout << "\nMatrix B:\n";
	printComplexMatrix(B);

	cout << "\nVector x:\n";
	printComplexVector(x);

	cout << "\nVector y:\n";
	printComplexVector(y);

	complex<double> complexdot1 = complexdotproduct(x,y);
	cout <<"\nx . y = " << complexdot1 << endl;

	complex<double> complexdot2 = complexdotproduct(y,x);
	cout <<"\ny . x = " << complexdot2 << endl;

	double normx = complexnorm(x);
	cout <<"\n|| x || = " << normx<< endl;

	double normy = complexnorm(y);
	cout <<"\n|| y || = " << normy<< endl;

	ComplexVector addvec = addComplexVectors(x,y);
	cout <<"\nx + y : " << endl;
	printComplexVector(addvec);

	ComplexVector subtractvec = subtractComplexVectors(x,y);
	cout <<"\nx - y : " << endl;
	printComplexVector(subtractvec);

	ComplexVector scalarproductvec = scalarmultiplicationComplexVector(x,3);
	cout <<"\nx * 3 : " << endl;
	printComplexVector(scalarproductvec);

	// Addition
	ComplexMatrix sum = addComplexMatrices(A, B);
	cout << "\nA + B:\n";
	printComplexMatrix(sum);

	// Subtraction
	ComplexMatrix subtraction = subtractComplexMatrices(A, B);
	cout << "\nA - B:\n";
	printComplexMatrix(subtraction);

	// Multiplication
	ComplexMatrix product = multiplyComplexMatrices(A, B);
	cout << "\nA * B:\n";
	printComplexMatrix(product);

	// Scalar Multiplication
	ComplexMatrix scalarproduct = scalarmultiplicationComplexMatrix(A, 3);
	cout << "\nA * 3:\n";
	printComplexMatrix(scalarproduct);

	//std::complex<double> c1(4.0, 3.0);
	//cout << conjugate(c1) << endl;

	//} catch (const std::exception &e) 
	//{
	//	cerr << "Error: " << e.what() << "\n";
	//}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
