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

// Helper function for dot product of two vectors
double dot_product(const vector<double>& v1, const vector<double>& v2) 
{
	double result = 0.0;
	for (size_t i = 0; i < v1.size(); ++i) 
	{
		result += v1[i] * v2[i];
	}
	return result;
}

// Helper function for vector subtraction (v1 - v2)
vector<double> vector_subtract(const vector<double>& v1, const vector<double>& v2) 
{
	vector<double> result(v1.size());
	for (size_t i = 0; i < v1.size(); ++i) 
	{
		result[i] = v1[i] - v2[i];
	}
	return result;
}

// Helper function for scalar-vector multiplication (scalar * v)
vector<double> scalar_multiply(double scalar, const vector<double>& v) 
{
	vector<double> result(v.size());
	for (size_t i = 0; i < v.size(); ++i) 
	{
		result[i] = scalar * v[i];
	}
	return result;
}

// Helper function to get a column from a matrix
vector<double> get_column(const vector<vector<double>>& A, int col_idx) 
{
	vector<double> column(A.size());
	for (size_t i = 0; i < A.size(); ++i) 
	{
		column[i] = A[i][col_idx];
	}
	return column;
}

// Function to perform QR decomposition using Modified Gram-Schmidt
void qr_decomposition_gram_schmidt(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) 
{
	int m = A.size(); // Rows
	int n = A[0].size(); // Columns

	Q.assign(m, std::vector<double>(n, 0.0));
	R.assign(n, std::vector<double>(n, 0.0));

	for (int i = 0; i < n; ++i) 
	{
		vector<double> a_i = get_column(A, i);
		vector<double> u_i = a_i;

		for (int j = 0; j < i; ++j) 
		{
			vector<double> q_j = get_column(Q, j);
			// r[j][i] = dot(q_j, a_i)
			R[j][i] = dot_product(q_j, a_i);
			// u_i = u_i - r[j][i] * q_j
			u_i = vector_subtract(u_i, scalar_multiply(R[j][i], q_j));
		}

		// Compute norm of u_i (which becomes the diagonal element R[i][i])
		double norm_u_i = std::sqrt(dot_product(u_i, u_i));
		R[i][i] = norm_u_i;
	
		// Normalize u_i to get q_i: q_i = u_i / norm_u_i
		for (int k = 0; k < m; ++k) 
		{
			Q[k][i] = u_i[k] / norm_u_i;
		}
	}
}
int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	dmat A = loadMatrixFromFile("matrix.txt");
	int n = A.size();
	int C = A[0].size();

	dmat Q(n, vector<double>(C));
	dmat R(n, vector<double>(C));

	qr_decomposition_gram_schmidt(A, Q, R);
	
	cout << "\nMatrix Q : " << endl;
	printMatrix(Q);
	cout << "\nMatrix R : " << endl; 
	printMatrix(R);
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
