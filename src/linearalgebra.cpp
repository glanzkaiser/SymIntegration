/*
   Thanks Freya the Goddess, RK, Sentinel, Berlin, Mother Mary.
*/
#include "symintegral/symintegrationc++.h"

using namespace std;

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DEFINE

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f
#define pi  3.1415926535897

#include <string>
#include<iostream>
#include<vector>
#include <fstream>
#include <cmath>
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()
#include <map>
#include <algorithm> // For std::max_element,  std::sort, std::reverse
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

vector<double> loadVectorFromFile(const string& filename) 
{
	vector<double> vecx;
	ifstream inputFile(filename);

	if (!inputFile.is_open()) 
	{
		cerr << "Error: Could not open file " << filename << endl;
		return vecx; // Return empty matrix on error
	}

	if (inputFile.is_open()) 
	{
		double num;
		while (inputFile >> num) 
		{ // Reads numbers separated by whitespace
			vecx.push_back(num);
		}
		inputFile.close();
		
	} 
	inputFile.close();
	return vecx;
}

// Function to print a matrix
void printMatrix(vector<vector<double>> matrix) 
{
	 for (const auto& row : matrix) 
	{
		for (double val : row) 
		{
			cout << std::fixed << setw(23) << setprecision(6) << val << setw(23);
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

vector<double> createVector(int n, double k)
{
	vector<double> resultVector(n, k);
	
	return resultVector;
}

vector<vector<double>> createMatrix(int R, int C, double k)
{
	vector<vector<double>> resultMatrix(R, vector<double>(C, k));
	
	return resultMatrix;
}

vector<vector<double>> addRow(vector<vector<double>> matrix, vector<double> vector_x, int k)
{
	int R = matrix.size();
	int C = matrix[0].size();
	vector<vector<double>> newMatrix(R, vector<double>(C));

	// Check if the number of new column elements matches the number of rows
	if (vector_x.size() != matrix.size()) 
	{
		cerr << "Error: Number of new row elements does not match the number of rows in the matrix." << endl;
	}
	for (int i = 0; i < R; ++i) 
	{
		for (int j = 0; j < C; ++j) 
		{
			newMatrix[i][j] = matrix[i][j];
		}
	}

	// Insert another row at a specific index (e.g., after the first row)
	// data.begin() + k points to the second element
	newMatrix.insert(newMatrix.begin() + k, vector_x);
	
	return newMatrix;
}

vector<vector<double>> addColumn(vector<vector<double>> matrix, vector<double> vector_x, int k)
{
	int R = matrix.size();
	int C = matrix[0].size();
	vector<vector<double>> newMatrix(R, vector<double>(C+1));

	// Check if the number of new column elements matches the number of rows
	if (vector_x.size() != matrix.size()) 
	{
		cerr << "Error: Number of new column elements does not match the number of columns in the matrix." << endl;
		 // Indicate an error
	}
	// Add the new column elements at index k
	for (int i = 0; i < R; ++i) 
	{
		for (int j = 0; j < C+1; ++j) 
		{
			if (j < k)
			{
				newMatrix[i][j] = matrix[i][j];
			}
			if (j == k)
			{
				newMatrix[i][j] = vector_x[i];
			}
			if (j > k)
			{
				newMatrix[i][j] = matrix[i][j-1];
			}
		}
	}

	return newMatrix;
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

vector<vector<double>> transpose(const vector<vector<double>> &matrix) // Done beautifully
{
	if (matrix.empty() || matrix[0].empty() ) 
	{
        	// Handle empty matrices or invalid dimensions
		return {};
	}	
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();

	vector<vector<double>> transposed_matrix(cols, vector<double>(rows, 0.0));
	for (size_t i = 0; i < rows; ++i) 
	{
		for (size_t j = 0; j < cols; ++j) 
		{
			transposed_matrix[j][i] = matrix[i][j];
			if(rows==1)
			{
				transposed_matrix[j][0] = matrix[0][j];
			}
			if(cols==1)
			{
				transposed_matrix[0][i] = matrix[i][0];
			}
		}
	}
	return transposed_matrix;
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
		cerr << "Error: Matrix dimensions are incompatible for multiplication." << endl;
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
		cerr << "Error: Matrix dimensions are incompatible." << std::endl;
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

void scalarmultiplication_alt(vector<vector<double>> &matrix, double scalar) 
{
	for (size_t i = 0; i < matrix.size(); ++i) 
	{
		for (size_t j = 0; j < matrix[i].size(); ++j) 
		{
			matrix[i][j] *= scalar; // Multiply each element by the scalar
		}
	}
}

vector<vector<double>> scalarmultiplication(vector<vector<double>> matrix, double scalar) 
{
	int R = matrix.size();
	int C = matrix[0].size();
	vector<vector<double>> vectorresult(R, vector<double>(C));
	for (size_t i = 0; i < matrix.size(); ++i) 
	{
		for (size_t j = 0; j < matrix[i].size(); ++j) 
		{
			matrix[i][j] *= scalar; // Multiply each element by the scalar
		}
	}
	for (size_t i = 0; i < matrix.size(); ++i) 
	{
		for (size_t j = 0; j < matrix[i].size(); ++j) 
		{
			vectorresult[i][j] = matrix[i][j];
		}
	}
	return vectorresult;
}

vector<double> add(vector<double> &vectorA, vector<double> &vectorB) 
{

	if (vectorA.empty() || vectorB.empty() )
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorA.size();
	int m = vectorB.size();
	
	if (n != m) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vectors dimensions are incompatible." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);

	// Perform vector addition
	for (int i = 0; i < n; ++i) 
	{
		result[i] += vectorA[i] + vectorB[i];	
	}
	return result;
}

vector<double> subtract(vector<double> &vectorA, vector<double> &vectorB) 
{

	if (vectorA.empty() || vectorB.empty() )
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorA.size();
	int m = vectorB.size();
	
	if (n != m) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vectors dimensions are incompatible." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);

	// Perform vector subtraction
	for (int i = 0; i < n; ++i) 
	{
		result[i] += vectorA[i] - vectorB[i];	
	}
	return result;
}

vector<double> scalarmultiplication(vector<double> &vectorA, double k) 
{

	if (vectorA.empty()  )
	{
        	// Handle empty vector
		return {};
	}

	int n = vectorA.size();
	
	// Initialize result vector with appropriate dimensions
	vector<double> result(n);

	// Perform scalar multiplicaiton for vector
	for (int i = 0; i < n; ++i) 
	{
		result[i] = vectorA[i]*k;	
	}
	return result;
}

vector<double> crossproduct(vector<double> &vectorA, vector<double> &vectorB) 
{

	if (vectorA.empty() || vectorB.empty() )
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorA.size();
	int m = vectorB.size();
	
	if (n != m && n!=3 && m!=3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vectors dimensions are incompatible or not in 3-space." << endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	
	result[0] = vectorA[1]*vectorB[2] - vectorA[2]*vectorB[1];	
	result[1] = -(vectorA[0]*vectorB[2] - vectorA[2]*vectorB[0]);
	result[2] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0];

	return result;
}

double scalartripleproduct(vector<double> &vectorA, vector<double> &vectorB, vector<double> &vectorC) 
{

	if (vectorA.empty() || vectorB.empty() || vectorC.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorA.size();
	int m = vectorB.size();
	int l = vectorC.size();

	if (n!=3 && m!=3 && l !=3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vectors dimensions are incompatible or not in 3-space." << endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> cresult(n);
	
	cresult[0] = vectorB[1]*vectorC[2] - vectorB[2]*vectorC[1];	
	cresult[1] = -(vectorB[0]*vectorC[2] - vectorB[2]*vectorC[0]);
	cresult[2] = vectorB[0]*vectorC[1] - vectorB[1]*vectorC[0];

	double finalresult;

	for (int i=0 ; i<n ; i++) 
	{
		finalresult += vectorA[i] * cresult[i] ;		
	}
	return finalresult;
}

double determinant_alt(vector<vector<double>> matrix) // alternative code to compute determinant
{
	int n = matrix.size();

	// Base cases
	if (n == 1) 
	{
		return matrix[0][0];
	}
	if (n == 2) 
	{
		return (matrix[0][0] * matrix[1][1]) - (matrix[0][1] * matrix[1][0]);
	}

	double result = 0.0;
	for (int p = 0; p < n; ++p) 
	{ // Iterate through the first row
		vector<vector<double>> submatrix;
		// Construct the submatrix by excluding the first row and the p-th column
		for (int i = 1; i < n; ++i) 
		{
			vector<double> row;
			for (int j = 0; j < n; ++j) 
			{
				if (j != p) 
				{
					row.push_back(matrix[i][j]);
				}
			}
		submatrix.push_back(row);
		}
	// Add the term to the determinant using cofactor expansion
	result += matrix[0][p] * std::pow(-1, p) * determinant_alt(submatrix);
	}
	return result;
}

long double determinant(const vector<vector<double>> &matrix) 
{
	int n = matrix.size();
	if (n == 0) 
	{
		return 0.0;
	}
	if (n == 1) 
	{
		return matrix[0][0];
	}

	long double det = 0.0;
	if (n == 2) 
	{
		return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
	}

	for (int p = 0; p < n; ++p) 
	{
		vector<vector<double>> submatrix(n - 1, vector<double>(n - 1));
		for (int i = 1; i < n; ++i) 
		{
			int q = 0;
			for (int j = 0; j < n; ++j) 
			{
				if (j == p) 
				{
					continue;
				}
			submatrix[i - 1][q++] = matrix[i][j];
			}
		}
		det += (p % 2 == 0 ? 1 : -1) * matrix[0][p] * determinant(submatrix);
	}
	return det;
}

// Function to get the submatrix for calculating minors
void getCofactor(vector<vector<double>> &mat, vector<vector<double>> &temp, int p, int q, int n) 
{
	int i = 0, j = 0;
	for (int row = 0; row < n; row++) 
	{
		for (int col = 0; col < n; col++) 
		{
			// Copy elements except for the given row and column
			if (row != p && col != q) 
			{
				temp[i][j++] = mat[row][col];
				if (j == n - 1) 
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

// Function to calculate the determinant of a matrix
double determinant1(vector<vector<double>> &mat, int n) 
{
	if (n == 1) 
	{
		return mat[0][0];
	}

	double det = 0;
	vector<vector<double>> temp(n - 1, vector<double>(n - 1));
	int sign = 1;

	for (int f = 0; f < n; f++) 
	{
		getCofactor(mat, temp, 0, f, n);
		det += sign * mat[0][f] * determinant1(temp, n - 1);
		sign = -sign;
	}
	return det;
}

// Function to calculate the adjugate of a matrix
vector<vector<double>> adjugate(vector<vector<double>> &mat) 
{
	int n = mat.size();
	int c = mat[0].size();
	if (n == 0 || c != n) 
	{
		// Handle non-square or empty matrices
		return {}; 
	}

	vector<vector<double>> adj(n, vector<double>(n));

	if (n == 1) 
	{
		adj[0][0] = 1; // Adjugate of a 1x1 matrix [a] is [1]
		return adj;
	}

	vector<vector<double>> temp(n - 1, vector<double>(n - 1));

	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < n; j++) 
		{
		getCofactor(mat, temp, i, j, n);
		// Calculate cofactor and then transpose it directly into adj
		adj[j][i] = std::pow(-1, i + j) * determinant1(temp, n - 1); 
		}
	}
	return adj;
}

// Function to calculate the inverse of a matrix
vector<vector<double>> inverse(vector<vector<double>> &matrix) 
{
	int n = matrix.size();
	int c = matrix[0].size();
	if (n == 0 || c != n) 
	{
		cerr << "Error: Matrix must be square and non-empty." << endl;
		return {};
	}

	double det = determinant(matrix);
	if (std::abs(det) < 1e-9) 
	{ // Check for near-zero determinant
		cerr << "Error: Matrix is singular, inverse does not exist." << endl;
		return {};
	}

	vector<vector<double>> adjugate_matrix = adjugate(matrix);
	vector<vector<double>> inverse_matrix(n, vector<double>(n));

	for (int i = 0; i < n; ++i) 
	{
	for (int j = 0; j < n; ++j) 
	{
		inverse_matrix[i][j] = adjugate_matrix[i][j] / det;
        }
	}
	return inverse_matrix;
}

double norm(const vector<double> &vectorx)
{
	double result;
	int n = vectorx.size();
		
	for (int i=0 ; i<n ; i++) 
	{
		result += vectorx[i]*vectorx[i] ;		
	}
	return sqrt(result);
}

double dot(const vector<double> &vectorx, const vector<double> &vectory)
{
	double result;
	int n = vectorx.size();
	int m = vectory.size();
	if (n !=m)
	{
		cerr << "Error: Vectors are of different dimension" << endl;
		return {};
	}
	for (int i=0 ; i<n ; i++) 
	{
		result += vectorx[i] * vectory[i] ;		
	}
	return result;
}

double distance(const vector<double> &vectorx, const vector<double> &vectory)
{
	double result;
	int n = vectorx.size();
	int m = vectory.size();
	if (n !=m)
	{
		cerr << "Error: Vectors are of different dimension" << endl;
		return {};
	}
	for (int i=0 ; i<n ; i++) 
	{
		result += (vectorx[i] - vectory[i])*(vectorx[i] - vectory[i]) ;		
	}
	return sqrt(result);
}

double distance(const Symbolic &fx, const Symbolic &x, const Symbolic &y,  const Symbolic &z, const vector<double> &vectorx)
{
	double result;
	int n = vectorx.size();
	if (n == 2)
	{
		double a = fx.coeff(x,1);
		double b = fx.coeff(y,1);

		double num = fx[x==vectorx[0],y==vectorx[1]];
		double denom = sqrt(a*a + b*b);
		result = abs(num)/denom;

		return result;
	}
	if (n == 3)
	{
		double a = fx.coeff(x,1);
		double b = fx.coeff(y,1);
		double c = fx.coeff(z,1);

		Equations rules = (  x == vectorx[0], y == vectorx[1], z == vectorx[2] );
		double num = fx.subst_all(rules);
		double denom = sqrt(a*a + b*b + c*c);
		result = abs(num)/denom;

		return result;
	}
	if (n!=2 && n !=3)
	{
		result = 0;	
	}
	return result;
}

SymbolicMatrix vectorequation(const Symbolic &fx, const Symbolic &x, const Symbolic &y,  const Symbolic &z)
{
	Matrix<Symbolic> B_mat(3,1);
	Symbolic t("t"), s("s");

	double a = fx.coeff(x,1);
	double b = fx.coeff(y,1);
	double c = fx.coeff(z,1);
	double d = fx.coeff(x,0).coeff(y,0).coeff(z,0);
	if (fx.coeff(x,1) != 0)
	{
		B_mat[0][0] = (-d - (b*t) - (c*s))/a;	
		B_mat[1][0] = t;
		B_mat[2][0] = s;
		if (fx.coeff(y,1) == 0 && fx.coeff(z,1) !=0)
		{		
			B_mat[0][0] = (-d - (c*s))/a;
			B_mat[1][0] = 0;
			B_mat[2][0] = s;
		}
		else if (fx.coeff(z,1) == 0 && fx.coeff(y,1) !=0)
		{	
			B_mat[0][0] = (-d - (b*t))/a;
			B_mat[1][0] = t;
			B_mat[2][0] = 0;
		}
	}
	if (fx.coeff(x,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = (-d - (c*s))/b;
		B_mat[2][0] = s;
	}
	if (fx.coeff(x,1) == 0 && fx.coeff(z,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = -d/b;
		B_mat[2][0] = 0;
	}
	if (fx.coeff(x,1) == 0 && fx.coeff(y,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = 0;
		B_mat[2][0] = -d/c;
	}
	if (fx.coeff(y,1) == 0 && fx.coeff(z,1) == 0)
	{
		B_mat[0][0] = -d/a;	
		B_mat[1][0] = 0;
		B_mat[2][0] = 0;
	}
	if (fx.coeff(x,1) == 0 && fx.coeff(y,1) == 0 && fx.coeff(z,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = 0;
		B_mat[2][0] = 0;
		cout << "The function is zero" << endl;
	}
	return B_mat;
}

void vectorequationdecomp(const Symbolic &fx, const Symbolic &x, const Symbolic &y,  const Symbolic &z, vector<double> &v1, vector<double> &v2, vector<double> &v3)
{
	Matrix<Symbolic> B_mat(3,1);
	Symbolic t("t"), s("s");

	double a = fx.coeff(x,1);
	double b = fx.coeff(y,1);
	double c = fx.coeff(z,1);
	double d = fx.coeff(x,0).coeff(y,0).coeff(z,0);
	if (fx.coeff(x,1) != 0)
	{
		B_mat[0][0] = (-d - (b*t) - (c*s))/a;	
		B_mat[1][0] = t;
		B_mat[2][0] = s;
		if (fx.coeff(y,1) == 0 && fx.coeff(z,1) !=0)
		{		
			B_mat[0][0] = (-d - (c*s))/a;
			B_mat[1][0] = 0;
			B_mat[2][0] = s;
		}
		else if (fx.coeff(z,1) == 0 && fx.coeff(y,1) !=0)
		{	
			B_mat[0][0] = (-d - (b*t))/a;
			B_mat[1][0] = t;
			B_mat[2][0] = 0;
		}
	}
	if (fx.coeff(x,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = (-d - (c*s))/b;
		B_mat[2][0] = s;
	}
	if (fx.coeff(x,1) == 0 && fx.coeff(z,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = -d/b;
		B_mat[2][0] = 0;
	}
	if (fx.coeff(x,1) == 0 && fx.coeff(y,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = 0;
		B_mat[2][0] = -d/c;
	}
	if (fx.coeff(y,1) == 0 && fx.coeff(z,1) == 0)
	{
		B_mat[0][0] = -d/a;	
		B_mat[1][0] = 0;
		B_mat[2][0] = 0;
	}
	if (fx.coeff(x,1) == 0 && fx.coeff(y,1) == 0 && fx.coeff(z,1) == 0)
	{
		B_mat[0][0] = 0;	
		B_mat[1][0] = 0;
		B_mat[2][0] = 0;
		cout << "The function is zero" << endl;
	}
	double c1 = B_mat[0][0].coeff(s,0).coeff(t,0);
	double c2 = B_mat[0][0].coeff(t,1);
	double c3 = B_mat[0][0].coeff(s,1);
	
	v1[0] = c1;
	v1[1] = 0;
	v1[2] = 0;
	v2[0] = c2;
	v2[1] = 1;
	v2[2] = 0;
	v3[0] = c3;
	v3[1] = 0;
	v3[2] = 1;

}

double angle(const vector<double> &vectorx, const vector<double> &vectory)
{
	int n = vectorx.size();
	int m = vectory.size();
	if (n !=m)
	{
		cerr << "Error: Vectors are of different dimension" << endl;
		return {};
	}
	double dotresult, normx,normy ;
	for (int i=0 ; i<n ; i++) 
	{
		dotresult += vectorx[i] * vectory[i] ;		
		normx += vectorx[i]*vectorx[i] ;		
		normy += vectory[i]*vectory[i] ;
	}
	normx = sqrt(normx);
	normy = sqrt(normy);

	return acos(dotresult/(normx*normy));
}

double degtorad(double x)
{
	return x*DEGTORAD;
}

double radtodeg(double x)
{
	return x*RADTODEG;
}

vector<double> orthogonalprojection(vector<double> vectorx, vector<double> vectory)
{
	int n = vectorx.size();
	int m = vectory.size();
	vector<double> w1 = vectory;

	if (n !=m)
	{
		cerr << "Error: Vectors are of different dimension" << endl;
		return {};
	}
	double dotresult, norma ;
	for (int i=0 ; i<n ; i++) 
	{
		dotresult += vectorx[i] * vectory[i] ;			
		norma += vectory[i]*vectory[i] ;
	}
	norma = sqrt(norma);
	double scalar = dotresult/(norma*norma);

	for (int i = 0; i < n; ++i) 
	{
		w1[i] *= scalar; // Multiply each element by the scalar	
	}
	return w1;
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
	reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	for (int i = 0; i < n; ++i) 
	{
		cout << "x" << i + 1 << " = " << fixed << setprecision(5) << x[i] << endl;
	}
	cout << endl;

}

void LUsolve_nhsystem(vector<vector<double>> &A, vector<vector<double>> &b, vector<double> &x) 
{
	int n = A.size();
	// Declare the matrix L and U with the size.
	dmat L(n, vector<double>(n));
	dmat U(n, vector<double>(n));

	// Initialize L with ones on the diagonal and zeros above
	// Initialize U with zeros below the diagonal	
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			L[i][j] = (i == j) ? 1.0 : 0.0; // Ones on diagonal for L
			U[i][j] = 0.0;
		}
	}
	
	// Doolittle's algorithm
	for (int i = 0; i < n; ++i) 
	{
		// Calculate U elements
		for (int j = i; j < n; ++j) 
		{
			double sum = 0.0;
			for (int k = 0; k < i; ++k) 
			{
				sum += L[i][k] * U[k][j];
			}
			U[i][j] = A[i][j] - sum;
		}
        

		// Calculate L elements (below the diagonal)
		for (int j = i + 1; j < n; ++j) 
		{
			double sum = 0.0;
			for (int k = 0; k < i; ++k) 
			{
				sum += L[j][k] * U[k][i];
			}
			L[j][i] = (A[j][i] - sum) / U[i][i]; // U[i][i] cannot be zero
		}
	}

	
	// Lz = b
	vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = L[i][j];
			if(j==n)
			{
				augmentedMatrix[i][j] = b[i][0];
			}
		}
	}
	//cout << "Augmented Matrix for Lz = b :" << endl;	
	//printMatrix(augmentedMatrix);
	
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
	
	//cout << "\nAugmented Matrix for Lz = b in reduced row form:" << endl;
	//printMatrix(augmentedMatrix);

	// Back Substitution
	vector<double> z(n);
	for (int i = n - 1; i >= 0; --i) 
	{
		double sum = 0.0;
		for (int j = i + 1; j < n; ++j) 
		{
			sum += augmentedMatrix[i][j] * z[j];
		}
		z[i] = (augmentedMatrix[i][n] - sum) / augmentedMatrix[i][i];
		
	}	
	
	// Print Solution
	/*cout << "\nSolution for Lz = b:" << endl;
	for (int i = 0; i < n; ++i) 
	{
		cout << "z" << i + 1 << " = " << fixed << setprecision(5) << z[i] << endl;
	}
	cout << endl;*/
	
	// Ux = z
	vector<vector<double>> augmentedMatrix2(n, vector<double>(n + 1));
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix2[i][j] = U[i][j];
			if(j==n)
			{
				augmentedMatrix2[i][j] = z[i];
			}
		}
	}
	//cout << "Augmented Matrix for Ux = z :" << endl;	
	//printMatrix(augmentedMatrix);

	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (abs(augmentedMatrix2[k][i]) > abs(augmentedMatrix2[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(augmentedMatrix2[i], augmentedMatrix2[pivotRow]);

		// Check for singular matrix (no unique solution)
		if (abs(augmentedMatrix2[i][i]) < 1e-9) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
	 		
		}

		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			double factor = augmentedMatrix2[k][i] / augmentedMatrix2[i][i];
			for (int j = i; j <= n; ++j) 
			{ // Iterate up to n for the constant term
				augmentedMatrix2[k][j] -= factor * augmentedMatrix2[i][j];
			}
		}
	}
	
	//cout << "\nAugmented Matrix for Ux = z in reduced row form:" << endl;
	//printMatrix(augmentedMatrix);

	// Back Substitution
	vector<double> solution(n);
	for (int i = n - 1; i >= 0; --i) 
	{
		double sum = 0.0;
		for (int j = i + 1; j < n; ++j) 
		{
			sum += augmentedMatrix2[i][j] * solution[j];
		}
		solution[i] = (augmentedMatrix2[i][n] - sum) / augmentedMatrix2[i][i];
		x.push_back(solution[i]);
	}	
	
	// Print Solution
	//cout << "\nSolution for Ux = z:" << endl;
	reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	/*for (int i = 0; i < n; ++i) 
	{
		cout << "x" << i + 1 << " = " << fixed << setprecision(5) << x[i] << endl;
	}
	cout << endl;
	*/
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