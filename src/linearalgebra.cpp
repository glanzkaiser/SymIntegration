/*
   Thanks Freya the Goddess, RK, Sentinel, Berlin, Mother Mary.
*/
#include "symintegral/symintegrationc++.h"

using namespace std;

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DEFINE

#include <string>
#include<iostream>
#include<vector>
#include <fstream>
#include <cmath>
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()
#include <map>
#include <algorithm> // For std::max_element,  std::sort
#include <numeric> // For std::accumulate
#include <random> // For random number generation

// Function to load a matrix from a text file
//template <typename T>
vector<vector<double>> loadMatrixFromFile(const string& filename) 
{
	vector<vector<double>> matrix;
	ifstream inputFile(filename);

	if (!inputFile.is_open()) 
	{
	cerr << "Error: Could not open file " << filename << endl;
	return matrix; // Return empty matrix on error
	}

	string line;
	while (getline(inputFile, line)) 
	{
	if (line.empty()) 
	{ // Skip empty lines
            continue;
        }
	istringstream iss(line);
	vector<double> row;
	double value;
        while (iss >> value) 
	{
		row.push_back(value);
	}
	matrix.push_back(row);
	}

	inputFile.close();
	return matrix;
}

// Function to print a matrix
void printMatrix(vector<vector<double>> matrix) 
{
	 for (const auto& row : matrix) 
	{
		for (double val : row) 
		{
			cout << std::fixed << setw(17) << setprecision(6) << val << setw(17);
		}
		cout << endl;
	}
}

void printVector(vector<double> vectorx) 
{
	 int n = vectorx.size();
	for (int i=0 ; i<n ; i++) 
	{
		cout << setprecision(6) << vectorx[i] << endl;		
	}
}

// Function to extract a column vector
vector<double> getColumn(vector<vector<double>> matrix, int columnIndex) 
{
	vector<double> columnVector;

	int R = matrix.size();
	
	// Check for an empty matrix or invalid column index
	if (matrix.empty() || columnIndex < 0) 
	{
		return columnVector; // Return an empty vector
	}
	for(int i=0; i<R; i++)
	{
		columnVector.push_back(static_cast<double>(matrix[i][columnIndex])); 
	}
	// Iterate through each row and add the element at columnIndex to the columnVector
	/*for (const auto& row : matrix) 
	{
		if (columnIndex < row.size()) 
		{ // Ensure the row has the specified column
			columnVector.push_back(row[columnIndex]);
		} 
		else 
		{
			// Handle cases where rows might have different lengths or the column is out of bounds for a specific row
			// You might choose to throw an exception, add a default value, or simply skip
			// For this example, we'll skip if the column doesn't exist in a row.
 			// Consider your specific error handling requirements here.
		}
	}*/
	return columnVector;
}

// Function to extract a row vector
vector<double> getRow(vector<vector<double>> matrix, int rowIndex) 
{
	vector<double> rowVector;

	int C = matrix[0].size();
	
	// Check for an empty matrix or invalid row index
	if (matrix.empty() || rowIndex < 0) 
	{
		return rowVector; // Return an empty vector
	}
	for(int i=0; i<C; i++)
	{
		rowVector.push_back(static_cast<double>(matrix[rowIndex][i])); 
	}
	
	return rowVector;
}

vector<vector<double>> multiply(vector<vector<double>> &matrixA, vector<vector<double>> &matrixB) 
{

	if (matrixA.empty() || matrixB.empty() || matrixB[0].empty()) 
	{
        	// Handle empty matrices or invalid dimensions
		return {};
	}

	size_t rowsA = matrixA.size();
	size_t colsA = matrixA[0].size();
	size_t rowsB = matrixB.size();
	size_t colsB = matrixB[0].size();

	if (colsA != rowsB) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions are incompatible for multiplication." << std::endl;
		return {};
	}

	// Initialize result matrix with appropriate dimensions
	vector<vector<double>> result(rowsA, vector<double>(colsB, 0.0));

	// Perform matrix multiplication
	for (size_t i = 0; i < rowsA; ++i) 
	{
		for (size_t j = 0; j < colsB; ++j) 
		{
			for (size_t k = 0; k < colsA; ++k) 
			{ // or rowsB, they are equal
				result[i][j] += matrixA[i][k] * matrixB[k][j];
			}
		}
	}
	return result;
}

vector<vector<double>> add(vector<vector<double>> &matrixA, vector<vector<double>> &matrixB) 
{

	if (matrixA.empty() || matrixB.empty() || matrixB[0].empty()) 
	{
        	// Handle empty matrices or invalid dimensions
		return {};
	}

	size_t rowsA = matrixA.size();
	size_t colsA = matrixA[0].size();
	size_t rowsB = matrixB.size();
	size_t colsB = matrixB[0].size();

	if (colsA != colsB && rowsA != rowsB) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions are incompatible for multiplication." << std::endl;
		return {};
	}

	// Initialize result matrix with appropriate dimensions
	vector<vector<double>> result(rowsA, vector<double>(colsB, 0.0));

	// Perform matrix addition
	for (size_t i = 0; i < rowsA; ++i) 
	{
		for (size_t j = 0; j < colsB; ++j) 
		{
			result[i][j] += matrixA[i][j] + matrixB[i][j];	
		}
	}
	return result;
}

Symbolic gaussianelimination(const SymbolicMatrix &A, int n)
{
	Matrix<Symbolic> B_mat(n,n+1);
	vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));
	
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = A[i][j];
			B_mat[i][j] = A[i][j];
		}
	}
	cout << "A :\n" << B_mat <<endl;
	
	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (abs(augmentedMatrix[k][i]) > abs(augmentedMatrix[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);

		// Check for singular matrix (no unique solution)
		if (abs(augmentedMatrix[i][i]) < 1e-9) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
	 		return 0;
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

	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			B_mat[i][j] = augmentedMatrix[i][j];
		}
	}
	
	cout << "A (in row reduced echelon form) :\n" << B_mat << endl;

	// Print Solution
	cout << "Solution:" << endl;
	for (int i = 0; i < n; ++i) 
	{
	cout << "x" << i + 1 << " = " << fixed << setprecision(5) << solution[i] << endl;
	}
	cout << endl;
	return 0;
}

void solve_nhsystem(vector<vector<double>> &A, vector<vector<double>> &b, vector<double> &x) 
{
	int n = A.size();
	vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = A[i][j];
			if(j==n)
			{
				augmentedMatrix[i][j] = b[i][0];
			}
		}
	}
	cout << "Augmented Matrix:" << endl;	
	printMatrix(augmentedMatrix);

	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (abs(augmentedMatrix[k][i]) > abs(augmentedMatrix[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);

		// Check for singular matrix (no unique solution)
		if (abs(augmentedMatrix[i][i]) < 1e-9) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
	 		
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
	
	cout << "\nAugmented Matrix in reduced row form:" << endl;
	printMatrix(augmentedMatrix);

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
		x.push_back(solution[i]);
	}	
	
	// Print Solution
	cout << "\nSolution:" << endl;
	for (int i = 0; i < n; ++i) 
	{
	cout << "x" << i + 1 << " = " << fixed << setprecision(5) << x[i] << endl;
	}
	cout << endl;

}

Symbolic gaussianeliminationtest(const SymbolicMatrix &A, int C, int R)
{
	Matrix<Symbolic> B_mat(R,C);
	
	for(int i = 0; i<R; i++)
	{
		for(int j=0; j<C; j++)
		{
			B_mat[i][j] = A[i][j];
		}
	}
	cout << "A :\n" << B_mat <<endl;
	
	for (int k=0; k<R; k++)
	{
		int i_max = k;
	
		if (i_max != k)
		{
			for (int m=0; m<C; k++)
			{
				Symbolic temp = B_mat[k][m];
				B_mat[k][m] = B_mat[i_max][m];
				B_mat[i_max][m] = temp;
			}
		}
		for (int i=k+1; i<R; i++)
		{
			 // factor f to set current row kth element to 0 and subsequently remaining kth column to 0 
			Symbolic f = B_mat[i][k]/B_mat[k][k];
			// subtract fth multiple of corresponding kth row element
			for (int j=k+1; j<C; j++)
			{
				B_mat[i][j] -= B_mat[k][j]*f;
			}
			B_mat[i][k] = 0;
		}
		
	}
	
	cout << "A (in reduced row form) :\n" << B_mat <<endl;
	for (int i = 0; i < R; i++)
		{
			double pivot = B_mat[i][i];	
			if (B_mat[i][i] != 0)
			{
				pivot = B_mat[i][i];					
			}
			else if (B_mat[i][i] == 0)
			{
				int j = 0;
				while (j < C)
				{
					if (B_mat[i][j] !=0)
					{
						pivot = B_mat[i][j];
					}
					j++;
				}					
			}
			for (int j = 0; j < C; j++) 
			{	
				B_mat[i][j] = B_mat[i][j]/pivot;
			}
		}

	cout << "A (in row reduced echelon form) :\n" << B_mat << endl;

	return 0;
}

#endif
#endif