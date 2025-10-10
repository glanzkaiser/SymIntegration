// g++ -o result main.cpp -lsymbolicc++
// Merci beaucoup Freya et Sentinel

#include<bits/stdc++.h>
#include<iostream>
#include "symintegrationc++.h"
#include<vector>
#include <chrono>
#include <algorithm> // For std::next_permutation
#include <string>
#include <iomanip> // For std::setprecision
using namespace std::chrono;
using namespace std;

#define R 2 // number of rows
#define C 7 // number of columns

// Function to perform Gaussian Elimination
void gaussianElimination(vector<vector<double>>& augmentedMatrix, int n) {
	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (std::abs(augmentedMatrix[k][i]) > std::abs(augmentedMatrix[pivotRow][i])) 
			{
			pivotRow = k;
			}
		}
		swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);

		// Check for singular matrix (no unique solution)
		if (std::abs(augmentedMatrix[i][i]) < 1e-9) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
	 		return;
		}

		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
			for (int j = i; j <= n; ++j) 
			{ // Iterate up to n for the constant term
			augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
			}
		}
	}

	// Back Substitution
	vector<double> solution(n);
	for (int i = n - 1; i >= 0; --i) 
	{
		double sum = 0.0;
		for (int j = i + 1; j < n; ++j) 
		{
			sum += augmentedMatrix[i][j] * solution[j];
		}
		solution[i] = (augmentedMatrix[i][n] - sum) / augmentedMatrix[i][i];
	}

	// Print Solution
	cout << "Solution:" << endl;
	for (int i = 0; i < n; ++i) 
	{
	cout << "x" << i + 1 << " = " << fixed << setprecision(4) << solution[i] << endl;
	}
}

void comb(int N, int K)
{
	string bitmask(K, 1); // K leading 1's
	bitmask.resize(N, 0); // N-K trailing 0's
 
	// print integers and permute bitmask
	do 
	{
	for (int i = 0; i < N; ++i) // [0..N-1] integers
	{
		if (bitmask[i]) 
		{
			cout << " " << i;
		}
	}
	cout << endl;
	} 
	while (prev_permutation(bitmask.begin(), bitmask.end()));
}

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	int n = 2;
 	// Construct a symbolic matrix of size 2 X 3
	Matrix<Symbolic> B_mat(2,3);
	B_mat[0][0] = 0.0333;		B_mat[0][1] = 0.2;	B_mat[0][2] = 1;
	B_mat[1][0] = 0.1;		B_mat[1][1] = 0.4;	B_mat[1][2] = 1;
	cout << "The matrix A:\n" << B_mat <<endl;
	cout << endl;
	vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));

	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = B_mat[i][j];
		}
	}
	gaussianElimination(augmentedMatrix, n);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}