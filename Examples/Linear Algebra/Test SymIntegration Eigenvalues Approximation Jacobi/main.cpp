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

#include <iostream>
#include <vector>
#include <cmath> // For std::abs, std::atan2, M_PI


// Function to perform Jacobi iteration
void jacobiEigen(vector<vector<double>>& A, vector<double>& eigenvalues, vector<vector<double>>& eigenvectors, double tolerance, int maxIterations) 
{
	int n = A.size();
	eigenvectors.assign(n, vector<double>(n, 0.0));
	for (int i = 0; i < n; ++i) 
	{
		eigenvectors[i][i] = 1.0; // Initialize eigenvectors as identity matrix
	}

	for (int iter = 0; iter < maxIterations; ++iter) 
	{
		// Find the largest off-diagonal element
		int p = 0, q = 1;
		double maxOffDiagonal = 0.0;
		for (int i = 0; i < n; ++i) 
		{
			for (int j = i + 1; j < n; ++j) 
			{
				if (std::abs(A[i][j]) > std::abs(maxOffDiagonal)) 
				{
					maxOffDiagonal = A[i][j];
					p = i;
					q = j;
				}
			}
 		}

		// Check for convergence
		if (abs(maxOffDiagonal) < tolerance) 
		{
			break;
		}

		// Calculate rotation angle
		double theta;
		if (A[p][p] == A[q][q]) 
		{
			theta = (maxOffDiagonal > 0) ? M_PI / 4.0 : -M_PI / 4.0;
		} 
		else 
		{
			theta = 0.5 * std::atan2(2.0 * A[p][q], A[p][p] - A[q][q]);
		}

		double cos_theta = std::cos(theta);
		double sin_theta = std::sin(theta);

		// Apply rotation to A
		double App = A[p][p];
		double Aqq = A[q][q];
		double Apq = A[p][q];

		A[p][p] = App * cos_theta * cos_theta - 2 * Apq * sin_theta * cos_theta + Aqq * sin_theta * sin_theta;
		A[q][q] = App * sin_theta * sin_theta + 2 * Apq * sin_theta * cos_theta + Aqq * cos_theta * cos_theta;
		A[p][q] = 0.0; // Zero out off-diagonal elements
		A[q][p] = 0.0;

		for (int i = 0; i < n; ++i) 
		{
			if (i != p && i != q) 
			{
				double Aip = A[i][p];
				double Aiq = A[i][q];
				A[i][p] = Aip * cos_theta - Aiq * sin_theta;
				A[p][i] = A[i][p]; // Maintain symmetry
				A[i][q] = Aip * sin_theta + Aiq * cos_theta;
				A[q][i] = A[i][q]; // Maintain symmetry
			}
		}

		// Apply rotation to eigenvectors
		for (int i = 0; i < n; ++i) 
		{
			double Vip = eigenvectors[i][p];
			double Viq = eigenvectors[i][q];
			eigenvectors[i][p] = Vip * cos_theta - Viq * sin_theta;
			eigenvectors[i][q] = Vip * sin_theta + Viq * cos_theta;
		}
	}

	// Extract eigenvalues from the diagonal of A
	eigenvalues.resize(n);
	for (int i = 0; i < n; ++i) 
	{
		eigenvalues[i] = A[i][i];
	}
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	dmat A = loadMatrixFromFile("matrix.txt");

	dvec eigenvalues;
	dmat eigenvectors;
	double tolerance = 1e-9;
	int maxIterations = 100;

	jacobiEigen(A, eigenvalues, eigenvectors, tolerance, maxIterations);

	std::cout << "Eigenvalues:" << std::endl;
	for (double val : eigenvalues) 
	{
		std::cout << val << std::endl;
	}

	std::cout << "\nEigenvectors:" << std::endl;
	printMatrix(eigenvectors);	

	//cout << pow(2.23849,3) - 5*pow(2.23849,2) + (8*2.23849) - 4 << endl;
	//cout << pow(2.12095,3) - 5*pow(2.12095,2) + (8*2.12095) - 4 << endl;
	//cout << pow(0.64056,3) - 5*pow(0.64056,2) + (8*0.64056) - 4 << endl;
	//cout << pow(1,3) - 5*pow(1,2) + (8*1) - 4 << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
