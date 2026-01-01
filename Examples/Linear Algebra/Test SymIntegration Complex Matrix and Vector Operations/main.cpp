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

	ComplexMatrix A = loadComplexMatrixFromFile("matrix.txt");
	ComplexVector x = loadComplexVectorFromFile("vector.txt");
	
	cout << "Matrix A:\n";
	printComplexMatrix(A);

	ComplexMatrix scalarproduct = scalarmultiplicationComplexMatrix(A, 3);
	cout << "\nA * 3:\n";
	printComplexMatrix(scalarproduct);
	
	cout << "\nVector x:\n";
	printComplexVector(x);

	ComplexVector Ax = multiplycomplexmatrixvector(A,x);
	cout << "\nA * x:" << endl;
	printComplexVector(Ax);

	std::complex<double> input(-1.0, 1.0);

	// Use the overloaded std::sqrt function for complex numbers
	std::complex<double> result = std::sqrt(input);

	std::complex<double> result2 = complexdivision(input,3);

	// Output the result
	std::cout << "The square root of " << input << " is " << result << std::endl;

	std::cout << "complex number /3 is " << result2 << std::endl;

	vector<complex<double>> y = complexvecrand_normal(2, 0.5, 4);
	printComplexVector(y);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
