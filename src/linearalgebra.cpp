/*
   Thanks Freya the Goddess, RK, Sentinel, Berlin, Mother Mary.
*/
// Search for  // NOT YET to see which function we haven't finish yet

#include "symintegral/symintegrationc++.h"

using namespace std;

// Type alias for convenience
using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

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
#include <algorithm> // For std::max_element,  std::sort, std::reverse ,  std::for_each
#include <numeric> // For std::accumulate
#include <random> // For random number generation
#include <complex>

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

ComplexMatrix loadComplexMatrixFromFile(const string& filename) 
{
	vector<vector<complex<double>>> matrix;
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
	vector<complex<double>> row;
	complex<double> value;
        while (iss >> value) 
	{
		row.push_back(value);
	}
	matrix.push_back(row);
	}

	inputFile.close();
	return matrix;
}

ComplexVector loadComplexVectorFromFile(const string& filename) 
{
	vector<complex<double>> vecx;
	ifstream inputFile(filename);

	if (!inputFile.is_open()) 
	{
		cerr << "Error: Could not open file " << filename << endl;
		return vecx; // Return empty matrix on error
	}

	if (inputFile.is_open()) 
	{
		complex<double> num;
		while (inputFile >> num) 
		{ // Reads numbers separated by whitespace
			vecx.push_back(num);
		}
		inputFile.close();
		
	} 
	inputFile.close();
	return vecx;
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

// Function to print a complex matrix
void printComplexMatrix(const ComplexMatrix &mat) 
{
	for (const auto &row : mat) 
	{
		for (const auto &val : row) 
		{
			cout << setw(12) << setprecision(6) << val << " ";
		}
		cout << "\n";
	}
}

// Function to print a complex vector
void printComplexVector(const ComplexVector &vec) 
{
	int n = vec.size();
	for (int i = 0; i< n; ++i) 
	{
		cout << setprecision(6) << vec[i] ;
		cout << "\n";
	}
}

vector<vector<double>> ComplextoRealMatrix(vector<vector<complex<double>>> &A) 
{
	int R = A.size();
	int C = A[0].size();

	vector<vector<double>> matrix(R, vector<double>(C, 0.0));
	for (int i = 0; i < R; ++i) 
	{
		for (int j = 0; j < C; ++j) 
		{
			matrix[i][j] = real(A[i][j]);
		}
	}
	return matrix;
}

complex<double> complexdivision(complex<double> a, double b)
{
	double acomplex = divisiond(real(a),b);
	double bcomplex = divisiond(imag(a),b);
	
	complex<double> result(acomplex, bcomplex);
	
	return result;
}

vector<complex<double>> complexvecrand_normal(double mu, double sigma, int n)
{
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(mu, sigma);
	for(int i=0; i<n; i++)
	{
		double real_part = distribution(generator);
		double imag_part = distribution(generator);
		complex<double> random_complex(real_part, imag_part);
		vec.push_back(random_complex); 	
	}
	return vec;
}

vector<complex<double>> complexvecrand_normal_zeroimaginary(double mu, double sigma, int n)
{
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
 	std::normal_distribution<double> distribution(mu, sigma);
	for(int i=0; i<n; i++)
	{
		double real_part = distribution(generator);
		double imag_part = 0;
		complex<double> random_complex(real_part, imag_part);
		vec.push_back(random_complex); 	
	}
	return vec;
}

vector<complex<double>> complexvecrand_uniform(double a, double b, int n)
{
	// 1. Obtain a seed:
	// Seeding with std::chrono::system_clock::now().time_since_epoch().count()
	// provides a more robust seed than a fixed value.
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	std::vector<complex<double>> vec;
	// 2. Define the distribution for floating-point numbers (e.g., uniform distribution)
	// Range [a, b)
	std::uniform_real_distribution<double> distribution(a, b);	
	for(int i=0; i<n; i++)
	{
		double real_part = distribution(generator);
		double imag_part = distribution(generator);
		complex<double> random_complex(real_part, imag_part);
		vec.push_back(random_complex); 	
	}
	return vec;
}


complex<double> conjugate(complex<double> a)
{
	complex<double> result(real(a), -imag(a));
	// The global function std::imag() takes the std::complex object as an argument

	return result;
}
double moduluscomplex(complex<double> c)
{
	double result;
	double a = real(c);
	double b = imag(c);
	result = sqrt(a*a + b*b);

	return result;
}
// Function to compute complex Euclidean norm of a vector
double complexnorm (const ComplexVector &x)
{
	int n = x.size();
	double result;
	for (int i = 0; i< n; ++i) 
	{
		double a = real(x[i]);
		double b = imag(x[i]);

		result += (a*a) + (b*b) ;
	}
	result = sqrt(result);
	return result;
}

// Function to compute complex Euclidean inner product / complex dot product
complex<double> complexdotproduct (const ComplexVector &x, const ComplexVector &y)
{
	if (x.size() != y.size()) 
	{
		throw std::invalid_argument("Vector dimensions must match for addition.");
	}
	int n = x.size();
	complex<double> result;
	for (int i = 0; i< n; ++i) 
	{
		result += x[i]*conjugate(y[i]);
	}
	return result;
}

// Function to add two complex vectors
ComplexVector addComplexVectors(const ComplexVector &x, const ComplexVector &y) 
{
	if (x.size() != y.size()) 
	{
		throw std::invalid_argument("Vector dimensions must match for addition.");
	}
	int n = x.size();
	ComplexVector result(n);
	for (int i = 0; i< n; ++i) 
	{
		result[i] = x[i] + y[i];
	}
	return result;
}

// Function to subtract two complex vectors
ComplexVector subtractComplexVectors(const ComplexVector &x, const ComplexVector &y) 
{
	if (x.size() != y.size()) 
	{
		throw std::invalid_argument("Vector dimensions must match for addition.");
	}
	int n = x.size();
	ComplexVector result(n);
	for (int i = 0; i< n; ++i) 
	{
		result[i] = x[i] - y[i];
	}
	return result;
}

// Function to make scalar multiplicaiton on a complex vector
ComplexVector scalarmultiplicationComplexVector(const ComplexVector &x, double k) 
{
	int n = x.size();
	ComplexVector result(n);
	for (int i = 0; i< n; ++i) 
	{
		result[i] = x[i]*k;
	}
	return result;
}

// function for complex number-vector multiplication ((a+bi) * v)
ComplexVector complexnumbermultiplicationComplexVector(const vector<complex<double>>& v, complex<double> scalar) 
{
	vector<complex<double>> result(v.size());
	for (size_t i = 0; i < v.size(); ++i) 
	{
		result[i] = scalar * v[i];
	}
	return result;
}

// function to get a column from a complex matrix
vector<complex<double>> getcolumnComplexMatrix(const vector<vector<complex<double>>>& A, int col_idx) 
{
	vector<complex<double>> column(A.size());
	for (size_t i = 0; i < A.size(); ++i) 
	{
		column[i] = A[i][col_idx];
	}
	return column;
}

// Function to add two complex matrices
ComplexMatrix addComplexMatrices(const ComplexMatrix &A, const ComplexMatrix &B) 
{
	if (A.size() != B.size() || A[0].size() != B[0].size()) 
	{
		throw std::invalid_argument("Matrix dimensions must match for addition.");
	}

	ComplexMatrix result(A.size(), std::vector<Complex>(A[0].size()));
	for (size_t i = 0; i < A.size(); ++i) 
	{
		for (size_t j = 0; j < A[i].size(); ++j) 
		{
			result[i][j] = A[i][j] + B[i][j];
		}
	}
	return result;
}

// Function to subtact two complex matrices
ComplexMatrix subtractComplexMatrices(const ComplexMatrix &A, const ComplexMatrix &B) 
{
	if (A.size() != B.size() || A[0].size() != B[0].size()) 
	{
		throw std::invalid_argument("Matrix dimensions must match for addition.");
	}

	ComplexMatrix result(A.size(), std::vector<Complex>(A[0].size()));
	for (size_t i = 0; i < A.size(); ++i) 
	{
		for (size_t j = 0; j < A[i].size(); ++j) 
		{
			result[i][j] = A[i][j] - B[i][j];
		}
	}
	return result;
}

// Function to multiply two complex matrices
ComplexMatrix multiplyComplexMatrices(const ComplexMatrix &A, const ComplexMatrix &B) 
{
	if (A[0].size() != B.size()) 
	{
		throw std::invalid_argument("Invalid dimensions for matrix multiplication.");
	}

	ComplexMatrix result(A.size(), vector<Complex>(B[0].size(), Complex(0, 0)));
	for (size_t i = 0; i < A.size(); ++i) 
	{
		for (size_t j = 0; j < B[0].size(); ++j) 
		{
			for (size_t k = 0; k < A[0].size(); ++k) 
			{
				result[i][j] += A[i][k] * B[k][j];
			}
		}
	}
	return result;
}

// Function to multiply a complex matrix with a constant
ComplexMatrix scalarmultiplicationComplexMatrix(const ComplexMatrix &A, double k) 
{
	ComplexMatrix result(A.size(), std::vector<Complex>(A[0].size()));
	for (size_t i = 0; i < A.size(); ++i) 
	{
		for (size_t j = 0; j < A[i].size(); ++j) 
		{
			result[i][j] = A[i][j]*k;
		}
	}
	return result;
}


// Function to multiply a complex matrix with a complex number
ComplexMatrix complexnumbermultiplicationComplexMatrix(const ComplexMatrix &A, complex<double> k) 
{
	ComplexMatrix result(A.size(), std::vector<Complex>(A[0].size()));
	for (size_t i = 0; i < A.size(); ++i) 
	{
		for (size_t j = 0; j < A[i].size(); ++j) 
		{
			result[i][j] = A[i][j]*k;
		}
	}
	return result;
}

ComplexVector multiplycomplexmatrixvector(const ComplexMatrix &A, const ComplexVector &x) 
{
	int R = A.size();
	int C = A[0].size();
	
	int n = x.size();
	if (C != n) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions and the vector size are incompatible for multiplication." << endl;
		return {};
	}
	ComplexVector result(x.size());
	for (int i = 0; i < R; ++i)
	{
		for (int j = 0; j < C; ++j)
		{
			result[i] += A[i][j]*x[j];
		}
	}

	return result;
}

vector<vector<complex<double>>> createIdentityComplexMatrix(int n)
{
	vector<vector<complex<double>>> resultMatrix;
	resultMatrix.assign(n, std::vector<complex<double>>(n, 0.0));
	
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			if(i == j)
			{
				resultMatrix[i][j] = 1;
			}
		}
	}
	return resultMatrix;
}

vector<vector<complex<double>>> createSubMatrix(vector<vector<complex<double>>> &A, int i1, int i2, int j1, int j2)
{
	int n = A.size();
	int m = A[0].size();

	if (i1 > i2 || j1 > j2 || i1 >= n || j1 >= m) 
	{
		throw std::invalid_argument("Wrong input.");
	}
	int R = i2 - i1 + 1;
	int C= j2 - j1 + 1;
	
	vector<vector<complex<double>>> resultMatrix;
	resultMatrix.assign(R, std::vector<complex<double>>(C, 0.0));
	
	int i1_index = i1;
	for (int i = 0; i < R; ++i) 
	{
		int j1_index = j1;
		for (int j = 0; j < C; ++j) 
		{
			resultMatrix[i][j] = A[i1_index][j1_index];
			j1_index = j1_index + 1;
		}
		i1_index = i1_index + 1;
	}
	return resultMatrix;
}

vector<vector<complex<double>>> TransposeComplexMatrix(vector<vector<complex<double>>> &A)
{
	int rows = A.size();
	int cols = A[0].size();

	vector<vector<complex<double>>> transposed_matrix;
	transposed_matrix.assign(rows, std::vector<complex<double>>(cols, 0.0));
	for (int i = 0; i < rows; ++i) 
	{
		for (int j = 0; j < cols; ++j) 
		{
			transposed_matrix[j][i] = A[i][j];
			if(rows==1)
			{
				transposed_matrix[j][0] = A[0][j];
			}
			if(cols==1)
			{
				transposed_matrix[0][i] = A[i][0];
			}
		}
	}

	return transposed_matrix;
}

vector<vector<complex<double>>> ComplexMatrixInverse(vector<vector<complex<double>>> &A)
{
	int n = A.size();
	int m = A[0].size();
	
	if (n != m) 
	{
		throw std::invalid_argument("Matrix is not square.");
	}

	vector<vector<complex<double>>> A_rref;
	vector<vector<complex<double>>> I_complex = createIdentityComplexMatrix(n);
	A_rref.assign(n, std::vector<complex<double>>(m, 0.0));
	
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			A_rref[i][j] = A[i][j];
		}
	}
	
	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (moduluscomplex(A_rref[k][i]) > moduluscomplex(A_rref[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(A_rref[i], A_rref[pivotRow]);
		swap(I_complex[i], I_complex[pivotRow]);
	
		// Check for singular matrix (no unique solution)
		if (moduluscomplex(A_rref[i][i]) < 1e-12) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
		}
			
		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			complex<double> factor = A_rref[k][i] / A_rref[i][i];
			for (int j = 0; j < n; ++j) 
			{ // Iterate up to n for the constant term
				A_rref[k][j] -= factor * A_rref[i][j];
				I_complex[k][j] -= factor * I_complex[i][j];
			}
		}
	}
	
	// make the leading 1
	for (int i = 0; i < n; ++i) 
	{
		double realpart = lround(real(A_rref[i][i]));
		if (realpart != 1.0 )
		{
			complex<double> pivot = A_rref[i][i];
			for (int j = 0; j < m; ++j) 
			{ 
				A_rref[i][j] = A_rref[i][j]/pivot;
				I_complex[i][j] = I_complex[i][j]/pivot;
			}	
		}
	}

	int f = 2;
	// Backward elimination
	// make zeros above all leading 1 / make matrix A into reduced row echelon form
	for (int i = n-2; i >= 0; --i) 
	{
		for (int k = 1; k < f; ++k)
		{
			double realpart = real(A_rref[i][i+k]);
			//cout << realpart << endl;		
			if (realpart != 0.0  )
			{
				complex<double> pivot = A_rref[i][i+k];
				for (int j = 0; j < m; ++j) 
				{ 
					A_rref[i][j] = A_rref[i][j] - (pivot * A_rref[i+k][j]);
					I_complex[i][j] = I_complex[i][j] - (pivot * I_complex[i+k][j]);
				}
			}
		}
	f = f+1;
	}

	return I_complex;
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

vector<vector<double>> createIdentityMatrix(int n)
{
	vector<vector<double>> resultMatrix(n, vector<double>(n, 0.0));
	for (int i = 0; i < n ;++i)
	{
		for (int j = 0; j < n ;++j)
		{
			if (i != j)
			{
				resultMatrix[i][j] = 0.0 ;
			}
			else if (i == j)
			{
				resultMatrix[i][j] = 1.0 ;
			}
		}
	}
	return resultMatrix;
}

// Function to create a matrix from a vector of column vectors
vector<vector<double>> createMatrixFromColumnVectors(const vector<vector<double>>& columnVectors) 
{

	if (columnVectors.empty()) 
	{
		return {}; // Return an empty matrix if no column vectors are provided
	}

	// Determine the number of rows (height) from the first column vector's size
	size_t numRows = columnVectors[0].size();
	// Determine the number of columns (width) from the number of column vectors
	size_t numCols = columnVectors.size();

	// Initialize the matrix with the correct dimensions
	vector<vector<double>> matrix(numRows, vector<double>(numCols));

	// Populate the matrix by iterating through column vectors and their elements
	for (size_t j = 0; j < numCols; ++j) 
	{ // Iterate through columns
		if (columnVectors[j].size() != numRows) 
		{
			// Handle error: column vectors must have consistent sizes
			// For simplicity, this example assumes consistent sizes.
			// In a real application, you might throw an exception or return an error.
			cerr << "Error: Column vector " << j << " has an inconsistent size." << endl;
			return {};
		}
		for (size_t i = 0; i < numRows; ++i) 	
		{ // Iterate through rows
			matrix[i][j] = columnVectors[j][i];
		}
	}
	return matrix;
}

vector<vector<double>> VandermondeMatrix(const vector<double>& vector_x, int n)
{
	int R = vector_x.size();
	int C = n;
	vector<vector<double>> VMatrix(R, vector<double>(C));

	// Add the new column elements at index k
	for (int i = 0; i < R; ++i) 
	{
		for (int j = 0; j < C; ++j) 
		{
			double vx = vector_x[i];
			VMatrix[i][j] = pow(vx,j);
		}
	}

	return VMatrix;
}


int MaxElementIndex(vector<double> vector_x)
{
	// Find iterator to the maximum element
	auto maxIt = std::max_element(vector_x.begin(), vector_x.end());

	// Calculate index from iterator
	size_t index = std::distance(vector_x.begin(), maxIt);

	// Output results, if using auto, then write *maxIt
	//cout << "Maximum value: " << *maxIt << "\n";
	//cout << "Location (index): " << index << "\n";
	//double h = *maxIt;

	int result = index;

	return result;
}

vector<vector<double>> PermutationMatrixMax(vector<double>& vector_x)
{
	int n = vector_x.size();

	vector<vector<double>> PermutationMatrix(n, vector<double>(n,0.0));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			if (i == j)
			{
				PermutationMatrix[i][j] = 1.0;
			}
			else if (i != j)
			{
				PermutationMatrix[i][j] = 0.0;
			}
		}
	}
	// Find iterator to the maximum element
	auto maxIt = std::max_element(vector_x.begin(), vector_x.end());

	// Calculate index from iterator
	size_t index = std::distance(vector_x.begin(), maxIt);

	int max_vector_index = index;
	swap(PermutationMatrix[0],PermutationMatrix[max_vector_index]);

	return PermutationMatrix;

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

vector<vector<double>> deleteRow(vector<vector<double>>& matrix, int rowIndex) 
{
	int R = matrix.size();
	int C = matrix[0].size(); 
	vector<vector<double>> matrix_final(R-1, vector<double>(C, 0.0));

	if (matrix.empty()) 
	{
		cerr << "Error: Matrix is empty" << endl;
		return {};
	}
	// Ensure the rowIndex is valid
	if (rowIndex < 0 || rowIndex >= R) 
	{
		cerr << "Error: Invalid row index." << endl;
		return {};
	}

	else
	{
		matrix.erase(matrix.begin() + rowIndex);
	}
	for (int i = 0; i < R-1; ++i) 
	{
		for (int j = 0; j < C; ++j) 
		{
				matrix_final[i][j] = matrix[i][j];
		}
	}
	return matrix_final;
}

vector<vector<double>> deleteColumn(vector<vector<double>>& matrix, int columnIndex) 
{
	int R = matrix.size();
	int C = matrix[0].size(); 
	vector<vector<double>> matrix_final(R, vector<double>(C-1, 0.0));
	
	if (matrix.empty()) 
	{
		cerr << "Error: Matrix is empty" << endl;
		return {};
	}
	// Ensure the columnIndex is valid for at least the first row
	if (columnIndex < 0 || columnIndex >= C) 
	{
		cerr << "Error: Invalid column index." << endl;
		return {};
	}

	// Iterate through each row and erase the element at the specified column index
	for (auto& row : matrix) 
	{
		row.erase(row.begin() + columnIndex);
	}
	for (int i = 0; i < R; ++i) 
	{
		for (int j = 0; j < C-1; ++j) 
		{
				matrix_final[i][j] = matrix[i][j];
		}
	}
	return matrix_final;
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
vector<double> multiplymatrixvector(vector<vector<double>> &matrixA, vector<double> &vectorX) 
{
	int R = matrixA.size();
	int C = matrixA[0].size();
	
	int n = vectorX.size();
	if (C != n) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions and the vector size are incompatible for multiplication." << endl;
		return {};
	}
	vector<double> result(R, 0.0);
	for (int i = 0; i < R; ++i)
	{
		for (int j = 0; j < C; ++j)
		{
			result[i] += matrixA[i][j]*vectorX[j];
		}
	}
	return result;
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

vector<vector<double>> matrixsubtraction(vector<vector<double>> &matrixA, vector<vector<double>> &matrixB) 
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

	if (rowsA != rowsB && colsA != colsB) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions are incompatible for subtraction." << endl;
		return {};
	}

	// Initialize result matrix with appropriate dimensions
	vector<vector<double>> result(rowsA, vector<double>(colsA, 0.0));

	// Perform matrix multiplication
	for (size_t i = 0; i < rowsA; ++i) 
	{
		for (size_t j = 0; j < colsA; ++j) 
		{
			result[i][j] = matrixA[i][j] - matrixB[i][j];
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

	if (rowsA != rowsB && colsA != colsB) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions are incompatible for addition." << endl;
		return {};
	}

	// Initialize result matrix with appropriate dimensions
	vector<vector<double>> result(rowsA, vector<double>(colsA, 0.0));

	// Perform matrix multiplication
	for (size_t i = 0; i < rowsA; ++i) 
	{
		for (size_t j = 0; j < colsA; ++j) 
		{
			result[i][j] = matrixA[i][j] + matrixB[i][j];
		}
	}
	return result;
}

vector<vector<double>> subtract(vector<vector<double>> &matrixA, vector<vector<double>> &matrixB) 
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

	if (rowsA != rowsB && colsA != colsB) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions are incompatible for subtraction." << endl;
		return {};
	}

	// Initialize result matrix with appropriate dimensions
	vector<vector<double>> result(rowsA, vector<double>(colsA, 0.0));

	// Perform matrix multiplication
	for (size_t i = 0; i < rowsA; ++i) 
	{
		for (size_t j = 0; j < colsA; ++j) 
		{
			result[i][j] = matrixA[i][j] - matrixB[i][j];
		}
	}
	return result;
}

double quadraticmultiplication(vector<vector<double>> &matrixA, vector<double> &vectorX)
{
	int R = matrixA.size();
	int C = matrixA[0].size();
	int n = vectorX.size();
	double quadraticresult = 0;

	if (C != n) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions and the vector size are incompatible for multiplication." << endl;
		return {};
	}
	vector<double> result(R, 0.0);
	for (int i = 0; i < R; ++i)
	{
		for (int j = 0; j < C; ++j)
		{
			result[i] += matrixA[i][j]*vectorX[j];
		}
	}
	
	for (int i = 0; i < n; ++i)
	{
		quadraticresult += vectorX[i]*result[i];
	}

	return quadraticresult;
}

complex<double> complexquadraticmultiplication(vector<vector<complex<double>>> &matrixA, vector<complex<double>> &vectorX)
{
	int R = matrixA.size();
	int C = matrixA[0].size();
	int n = vectorX.size();
	complex<double> quadraticresult;

	if (C != n) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Matrix dimensions and the vector size are incompatible for multiplication." << endl;
		return {};
	}
	vector<complex<double>> result(R);
	vector<complex<double>> X_conjugate(R);
	for (int i = 0; i < R; ++i)
	{
		for (int j = 0; j < C; ++j)
		{
			result[i] += matrixA[i][j]*vectorX[j];
		}
	}
	
	for (int i = 0; i < n; ++i)
	{
		X_conjugate[i] = conjugate(vectorX[i]);
	}
	for (int i = 0; i < n; ++i)
	{
		quadraticresult += X_conjugate[i]*result[i];
	}

	return quadraticresult;
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

void spanningtest(vector<vector<double>> &matrix)
{
	int n = matrix.size();

	cout << "\nMatrix:" << endl;
	printMatrix(matrix);
	if(determinant(matrix) == 0 )
	{
		cout << "\n The column vectors do not span R^{" << n<< "}" << endl;
	}
	else if(determinant(matrix) != 0 )
	{
		cout << "\n The column vectors span R^{" << n<< "}" << endl;
	}
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
	if (std::abs(det) < 1e-12) 
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
	double result = 0;
	int n = vectorx.size();
		
	for (int i=0 ; i<n ; i++) 
	{
		result += vectorx[i]*vectorx[i] ;		
	}
	return sqrt(result);
}

double dot(const vector<double> &vectorx, const vector<double> &vectory)
{
	double result = 0;
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

void linearindependencetest(vector<vector<double>> &A) 
{
	int n = A.size();
	vector<vector<double>> b(n, vector<double>(1,0.0));
	cout << "\nMatrix:" << endl;
	printMatrix(A);
	
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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
	
	// Back Substitution
	vector<double> solution(n);
	vector<double> x;
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
	
	reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	
	double x_accumulate = std::accumulate(x.begin(), x.end(), 0.0) ;
	if(x_accumulate == 0)
	{
		cout << "\nThe column vectors are linearly independent." << endl;
	}
	else if (x_accumulate != 0)
	{
		cout << "\nThe column vectors are not linearly dependent." << endl;
	}

}

void basistest(vector<vector<double>> &A) 
{
	int n = A.size();
	vector<vector<double>> b(n, vector<double>(1,0.0));
	cout << "\nMatrix:" << endl;
	printMatrix(A);
	if(determinant(A) == 0 )
	{
		cout << "\nThe column vectors do not span R^{" << n<< "}" << endl;
	}
	else if(determinant(A) != 0 )
	{
		cout << "\nThe column vectors span R^{" << n<< "}" << endl;
	}

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
	cout << "\nAugmented Matrix:" << endl;	
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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
	
	// Back Substitution
	vector<double> solution(n);
	vector<double> x;
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
	reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	/*for (int i = 0; i < n; ++i) 
	{
		cout << "x" << i + 1 << " = " << fixed << setprecision(5) << x[i] << endl;
	}*/

	double x_accumulate = std::accumulate(x.begin(), x.end(), 0.0) ;
	if(determinant(A) != 0 && x_accumulate ==0)
	{
		cout << "\nThe column vectors are the basis in R^{" << n<< "}" << endl;
	}
	else
	{
		cout << "\nThe column vectors are not the basis in R^{" << n<< "}" << endl;
	}

}

void coordinatevector(vector<vector<double>> &A, vector<vector<double>> &b) 
{
	int n = A.size();
	cout << "\nMatrix:" << endl;
	printMatrix(A);
	
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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
	
	// Back Substitution
	vector<double> solution(n);
	vector<double> x;
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
	
	reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	
	cout <<"\n(v)_{S}:" << endl;
	printVector(x);
	
}

void basistransition(vector<vector<double>> &oldbasis, vector<vector<double>> &newbasis) 
{
	int n = newbasis.size();
	int m = oldbasis.size();
	vector<vector<double>> transitionMatrix(n, vector<double>(n));
	vector<vector<double>> transitionMatrix_old(n, vector<double>(n));
	vector<vector<double>> transitionMatrix_full(n, vector<double>(2*n));
	
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			transitionMatrix[i][j] = newbasis[i][j];
			transitionMatrix_old[i][j] = oldbasis[i][j];
			transitionMatrix_full[i][j] = newbasis[i][j];
		}
		int k = 0;
		for(int j=n; j<2*n; j++)
		{
			transitionMatrix_full[i][j] = oldbasis[i][k];
			k = k+1;
		}
		
	}
	cout << "\n****************************************************************************************************************************" << endl;
	cout << "\nNew basis:" << endl;
	printMatrix(transitionMatrix);

	cout << "\nOld basis:" << endl;
	printMatrix(transitionMatrix_old);

	cout << "\nTransition matrix from old basis to new basis:" << endl;
	printMatrix(transitionMatrix_full);

	if(m != n)
	{
		cerr << "Error: The old basis and new basis have different size" << endl;	
	}		
	else if(m == n)
	{
		// Forward Elimination
		for (int i = 0; i < n; ++i) 
		{
			// Partial Pivoting (optional but recommended for stability)
			int pivotRow = i;
			for (int k = i + 1; k < n; ++k) 
			{
				if (abs(transitionMatrix_full[k][i]) > abs(transitionMatrix_full[pivotRow][i])) 
				{
					pivotRow = k;
				}
			}
			swap(transitionMatrix_full[i], transitionMatrix_full[pivotRow]);

			// Check for singular matrix (no unique solution)
			if (abs(transitionMatrix_full[i][i]) < 1e-12) // Using a small epsilon
			{ 
				cout << "No unique solution or infinite solutions exist." << endl;
		 		
			}

			// Eliminate elements below the pivot
			for (int k = i + 1; k < n; ++k) 
			{
				double factor = transitionMatrix_full[k][i] / transitionMatrix_full[i][i];
				for (int j = i; j < 2*n; ++j) 
				{ // Iterate up to n for the constant term
					transitionMatrix_full[k][j] -= factor * transitionMatrix_full[i][j];
				}
			}
		}
		
		// eliminate -1 in the leading 1 at every row
		for (int i = 0; i < n; ++i) 
		{
			if (transitionMatrix_full[i][i] < 0 || transitionMatrix_full[i][i] == -1)
			{
				for (int j = i; j < 2*n; ++j) 
				{ // Iterate up to n for the constant term
					transitionMatrix_full[i][j] = -1* transitionMatrix_full[i][j];
				}
			}
		}
		// make the leading 1
		for (int i = 0; i < n; ++i) 
		{
			double pivot = transitionMatrix_full[i][i];
			for (int j = i; j < 2*n; ++j) 
			{ // Iterate up to n for the constant term
				transitionMatrix_full[i][j] = transitionMatrix_full[i][j]/pivot;
			}
		}
		// make zeros above all leading 1
		for (int i = 0; i < n; ++i) 
		{
			if (transitionMatrix_full[i][i+1] != 0 )
			{
				double pivot = transitionMatrix_full[i][i+1];
				for (int j = i; j < 2*n; ++j) 
				{ // Iterate up to n for the constant term
					transitionMatrix_full[i][j] = transitionMatrix_full[i][j] - (pivot* transitionMatrix_full[i+1][j]);
				}
			}
		}
		cout << "\nTransition matrix in reduced row form:" << endl;
		printMatrix(transitionMatrix_full);
	}
	cout << "\n****************************************************************************************************************************" << endl;
	
}

void basistransition_withcoordinatevector(vector<vector<double>> &oldbasis, vector<vector<double>> &newbasis, vector<double> &w) 
{
	int n = newbasis.size();
	int m = oldbasis.size();
	int n_w = w.size();
	vector<vector<double>> transitionMatrix(n, vector<double>(n));
	vector<vector<double>> transitionMatrix_old(n, vector<double>(n));
	vector<vector<double>> transitionMatrix_full(n, vector<double>(2*n));
	vector<vector<double>> identityMatrix(n, vector<double>(n, 0.0));
	vector<vector<double>> transitionMatrix_fullwithidentity(n, vector<double>(2*n));
	vector<vector<double>> transitionMatrix_fullwithidentity2(n, vector<double>(2*n));
	
	for(int i = 0; i<n; i++)
	{
		identityMatrix[i][i] = 1.0;
		for(int j=0; j<n; j++)
		{
			transitionMatrix[i][j] = newbasis[i][j];
			transitionMatrix_old[i][j] = oldbasis[i][j];
			transitionMatrix_full[i][j] = newbasis[i][j];	
			transitionMatrix_fullwithidentity[i][j] = newbasis[i][j];		
			transitionMatrix_fullwithidentity2[i][j] = oldbasis[i][j];		
		}
		int k = 0;
		for(int j=n; j<2*n; j++)
		{
			transitionMatrix_full[i][j] = oldbasis[i][k];
			transitionMatrix_fullwithidentity[i][j] = identityMatrix[i][k];		
			transitionMatrix_fullwithidentity2[i][j] = identityMatrix[i][k];	
			k = k+1;
		}
		
	}
	cout << "\n****************************************************************************************************************************" << endl;
	cout << "\nNew basis (B):" << endl;
	printMatrix(transitionMatrix);

	cout << "\nOld basis (B'):" << endl;
	printMatrix(transitionMatrix_old);

	cout << "\nTransition matrix from old basis to new basis:" << endl;
	printMatrix(transitionMatrix_full);

	cout << "\nTransition matrix from elementary basis to new basis:" << endl;
	printMatrix(transitionMatrix_fullwithidentity);

	cout << "\nTransition matrix from elementary basis to old basis:" << endl;
	printMatrix(transitionMatrix_fullwithidentity2);

	if(m != n && n_w != n)
	{
		cerr << "Error: The old basis and new basis have different size, or/and the vector w has different size with the basis" << endl;	
	}		
	else if(m == n)
	{
		// Forward Elimination
		for (int i = 0; i < n; ++i) 
		{
			// Partial Pivoting (optional but recommended for stability)
			int pivotRow = i;	
			
			for (int k = i + 1; k < n; ++k) 
			{
				if (abs(transitionMatrix_full[k][i]) > abs(transitionMatrix_full[pivotRow][i])) 
				{
					pivotRow = k;
				}
			}
			swap(transitionMatrix_full[i], transitionMatrix_full[pivotRow]);

			// Check for singular matrix (no unique solution)
			if (abs(transitionMatrix_full[i][i]) < 1e-12) // Using a small epsilon
			{ 
				cout << "No unique solution or infinite solutions exist." << endl;
		 		
			}

			// Eliminate elements below the pivot
			for (int k = i + 1; k < n; ++k) 
			{
				double factor = transitionMatrix_full[k][i] / transitionMatrix_full[i][i];
				for (int j = i; j < 2*n; ++j) 
				{ 
					transitionMatrix_full[k][j] -= factor * transitionMatrix_full[i][j];
				}
			}

		}

		// Forward Elimination
		for (int i = 0; i < n; ++i) 
		{
			// Partial Pivoting (optional but recommended for stability)
			int pivotRow = i;	
			
			for (int k = i + 1; k < n; ++k) 
			{
				if (abs(transitionMatrix_fullwithidentity[k][i]) > abs(transitionMatrix_fullwithidentity[pivotRow][i])) 
				{
					pivotRow = k;
				}
			}
			swap(transitionMatrix_fullwithidentity[i], transitionMatrix_fullwithidentity[pivotRow]);

			// Check for singular matrix (no unique solution)
			if (abs(transitionMatrix_fullwithidentity[i][i]) < 1e-12) // Using a small epsilon
			{ 
				cout << "No unique solution or infinite solutions exist." << endl;
		 		
			}

			// Eliminate elements below the pivot
			for (int k = i + 1; k < n; ++k) 
			{
				double factor = transitionMatrix_fullwithidentity[k][i] / transitionMatrix_fullwithidentity[i][i];
				for (int j = i; j < 2*n; ++j) 
				{ 
					transitionMatrix_fullwithidentity[k][j] -= factor * transitionMatrix_fullwithidentity[i][j];
				}
			}

		}

		// Forward Elimination
		for (int i = 0; i < n; ++i) 
		{
			// Partial Pivoting (optional but recommended for stability)
			int pivotRow = i;	
			
			for (int k = i + 1; k < n; ++k) 
			{
				if (abs(transitionMatrix_fullwithidentity2[k][i]) > abs(transitionMatrix_fullwithidentity2[pivotRow][i])) 
				{
					pivotRow = k;
				}
			}
			swap(transitionMatrix_fullwithidentity2[i], transitionMatrix_fullwithidentity2[pivotRow]);

			// Check for singular matrix (no unique solution)
			if (abs(transitionMatrix_fullwithidentity2[i][i]) < 1e-12) // Using a small epsilon
			{ 
				cout << "No unique solution or infinite solutions exist." << endl;
		 		
			}

			// Eliminate elements below the pivot
			for (int k = i + 1; k < n; ++k) 
			{
				double factor = transitionMatrix_fullwithidentity2[k][i] / transitionMatrix_fullwithidentity2[i][i];
				for (int j = i; j < 2*n; ++j) 
				{ 
					transitionMatrix_fullwithidentity2[k][j] -= factor * transitionMatrix_fullwithidentity2[i][j];
				}
			}

		}
		
		// eliminate -1 in the leading 1 at every row
		for (int i = 0; i < n; ++i) 
		{
			if (transitionMatrix_full[i][i] < 0 || transitionMatrix_full[i][i] == -1)
			{
				for (int j = i; j < 2*n; ++j) 
				{ // Iterate up to n for the constant term
					transitionMatrix_full[i][j] = -1* transitionMatrix_full[i][j];
				}
			}
		}

		for (int i = 0; i < n; ++i) 
		{
			if (transitionMatrix_fullwithidentity[i][i] < 0 || transitionMatrix_fullwithidentity[i][i] == -1)
			{
				for (int j = i; j < 2*n; ++j) 
				{ // Iterate up to n for the constant term
					transitionMatrix_fullwithidentity[i][j] = -1* transitionMatrix_fullwithidentity[i][j];
				}
			}
		}

		for (int i = 0; i < n; ++i) 
		{
			if (transitionMatrix_fullwithidentity2[i][i] < 0 || transitionMatrix_fullwithidentity2[i][i] == -1)
			{
				for (int j = i; j < 2*n; ++j) 
				{ // Iterate up to n for the constant term
					transitionMatrix_fullwithidentity2[i][j] = -1* transitionMatrix_fullwithidentity2[i][j];
				}
			}
		}

		// make the leading 1
		for (int i = 0; i < n; ++i) 
		{
			double pivot = transitionMatrix_full[i][i];
			for (int j = i; j < 2*n; ++j) 
			{
				transitionMatrix_full[i][j] = transitionMatrix_full[i][j]/pivot;
			}
		}
		
		for (int i = 0; i < n; ++i) 
		{
			double pivot = transitionMatrix_fullwithidentity[i][i];
			for (int j = i; j < 2*n; ++j) 
			{
				transitionMatrix_fullwithidentity[i][j] = transitionMatrix_fullwithidentity[i][j]/pivot;
			}
		}

		for (int i = 0; i < n; ++i) 
		{
			double pivot = transitionMatrix_fullwithidentity2[i][i];
			for (int j = i; j < 2*n; ++j) 
			{
				transitionMatrix_fullwithidentity2[i][j] = transitionMatrix_fullwithidentity2[i][j]/pivot;
			}
		}

		// make zeros above all leading 1
		for (int i = 0; i < n-1; ++i) 
		{
			if (transitionMatrix_full[i][i+1] != 0 )
			{
				double pivot = transitionMatrix_full[i][i+1];
				for (int j = i; j < 2*n; ++j) 
				{ 
					transitionMatrix_full[i][j] = transitionMatrix_full[i][j] - (pivot * transitionMatrix_full[i+1][j]);
				}
			}
		}

		for (int i = 0; i < n-1; ++i) 
		{
			if (transitionMatrix_fullwithidentity[i][i+1] != 0 )
			{
				double pivot = transitionMatrix_fullwithidentity[i][i+1];
				
				for (int j = i; j < 2*n; ++j) 
				{ 
					transitionMatrix_fullwithidentity[i][j] = transitionMatrix_fullwithidentity[i][j] - (pivot* transitionMatrix_fullwithidentity[i+1][j]);
				}
			}
		}

		for (int i = 0; i < n-1; ++i) 
		{
			if (transitionMatrix_fullwithidentity2[i][i+1] != 0 )
			{
				double pivot = transitionMatrix_fullwithidentity2[i][i+1];
				
				for (int j = i; j < 2*n; ++j) 
				{ 
					transitionMatrix_fullwithidentity2[i][j] = transitionMatrix_fullwithidentity2[i][j] - (pivot* transitionMatrix_fullwithidentity2[i+1][j]);
				}
			}
		}

		// Create P_{E -> B}
		vector<vector<double>> P_EB(n, vector<double>(2*n));
		for (int i = 0; i < n; ++i) 
		{
			for (int j = 0; j < 2*n; ++j) 
			{ 
				P_EB[i][j] = transitionMatrix_fullwithidentity[i][j];
			}
		}
		// Create P_{E -> B'}
		vector<vector<double>> P_EBold(n, vector<double>(2*n));
		for (int i = 0; i < n; ++i) 
		{
			for (int j = 0; j < 2*n; ++j) 
			{ 
				P_EBold[i][j] = transitionMatrix_fullwithidentity2[i][j];
			}
		}
		// Delete the first two columns for P_{E -> B}
		for (int i = 0; i < n; ++i) 
		{
			for (auto& row : P_EB) 
			{
				row.erase(row.begin() + 0);
			}
		}
		// Delete the first two columns for P_{E -> B'}
		for (int i = 0; i < n; ++i) 
		{
			for (auto& row : P_EBold) 
			{
				row.erase(row.begin() + 0);
			}
		}
		// Create P_{B -> B'}
		vector<vector<double>> P_BoldBnew(n, vector<double>(2*n));
		for (int i = 0; i < n; ++i) 
		{
			for (int j = 0; j < 2*n; ++j) 
			{ 
				P_BoldBnew[i][j] = transitionMatrix_full[i][j];
			}
		}
		// Delete the first two columns for P_{B -> B'}
		for (int i = 0; i < n; ++i) 
		{
			for (auto& row : P_BoldBnew) 
			{
				row.erase(row.begin() + 0);
			}
		}

		vector<vector<double>> P_BnewBold = inverse(P_BoldBnew);
		// Matrix vector multiplication
		vector<double> wB;
		for (int i = 0; i < n; ++i) 
		{
			double sum = 0;
			for (int j = 0; j < n; ++j) 
			{ 
				sum += P_EB[i][j]*w[j];
			}
			wB.push_back(sum);
		}
		vector<double> wBnew;
		for (int i = 0; i < n; ++i) 
		{
			double sum = 0;
			for (int j = 0; j < n; ++j) 
			{ 
				sum += P_BnewBold[i][j]*wB[j];
			}
			wBnew.push_back(sum);
		}
		
		cout << "\n[ I | transition from B' to B]:" << endl;
		printMatrix(transitionMatrix_full);

		cout << "\n[ I | transition from E to B]:" << endl;
		printMatrix(transitionMatrix_fullwithidentity);

		cout << "\nFrom matrix [ B | B' ],\nP_{B' -> B}:" << endl;
		printMatrix(P_BoldBnew);
		
		cout << "\nFrom matrix [ B' | B ],\nP_{B -> B'}:" << endl;
		printMatrix(P_BnewBold);
		
		cout << "\nP_{E -> B}:" << endl;
		printMatrix(P_EB);
		
		cout << "\nP_{E -> B'}:" << endl;
		printMatrix(P_EBold);

		cout << "\nw:" << endl;
		printVector(w);
		
		cout << "\nw_{B}:" << endl;
		printVector(wB);

		cout << "\nw_{B'}:" << endl;
		printVector(wBnew);
		
	}
	cout << "\n****************************************************************************************************************************" << endl;
	
}

void rowspacebasis(vector<vector<double>> &A) 
{
	int r = A.size();
	int c = A[0].size();
	vector<vector<double>> R(r, vector<double>(c));
	
	for(int i = 0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			R[i][j] = A[i][j];
		}		
	}
	cout << "\n****************************************************************************************************************************" << endl;
	cout << "\nA:" << endl;
	printMatrix(A);
	
		
	// Forward Elimination
	for (int i = 0; i < r; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < r; ++k) 
		{
			if (abs(R[k][i]) > abs(R[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(R[i], R[pivotRow]);

		

		// Eliminate elements below the pivot
		for (int k = i + 1; k < r; ++k) 
		{
			if (R[i][i] != 0)
			{
				double factor = R[k][i] / R[i][i];
				
				for (int j = i; j < c; ++j) 
				{ 
					R[k][j] -= factor * R[i][j];
				
				}
			}
			else if (R[i][i] ==0)
			{
				i=i+1;
			}
		}
	}

	// Eliminate row that is a linear combination of other row
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if(R[j][i] != 0)
				{
					double pivot = divisiond(R[i][i],R[j][i]);
					for(int k=i; k<c; ++k)
					{
						R[j][k] = pivot*R[j][k] - R[i][k];
					}
				}
				else if(R[j][i] == 0)
				{
					i=i+1;
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i+1; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int k = i+1; k < r; ++k) 
					{ 
						if(R[k][j] != 0)
						{
							double pivot = R[i][j]/R[k][j];
							
							for(int m=i; m<c; ++m)
							{
								R[k][m] = pivot*R[k][m] - R[i][m];
							}
						}
						else if(R[k][j] == 0)
						{
							i=i+1;
						}
					}
				}
				else if (R[i][j] == 0)
				{
					i = i+1;
				}
			}
		}
	}
	
		
	// eliminate -1 in the leading 1 at every row
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] < 0 || R[i][i] == -1)
		{
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = -1* R[i][j];
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i+1; j < c; ++j) 
			{ 
				if (R[i][j] < 0 || R[i][j] == -1)
				{
					for(int k=j; k<c; k++)
					{
						R[i][k] = -1* R[i][k];
					}
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					i = i+1;
				}
			}
		}
	}
	
	// make the leading 1
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			double pivot = R[i][i];
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = R[i][j]/pivot;
			}	
		}
		else if(R[i][i] ==0)
		{
			for (int j = i+1; j < c; ++j) 
			{
				double pivot = R[i][j];
				 if (R[i][j] != 0)
				{
					
					for (int k = j; k < c; ++k) 
					{
						R[i][k] = R[i][k]/pivot;
					}				
					j = c-1;
				}	
				else if (R[i][j] ==0)
				{
					i = i+1;
				}
			}
		}
	}
	
	// Another forward elimination
	for (int i = 0; i < r-1; ++i) 
	{	
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if (R[j][i] == 0)
				{
					
				}
				else if(R[j][i] !=0)
				{
					double pivot = R[i][i]/R[j][i];
					for(int k=i; k<c; ++k)
					{
						R[j][k] = R[j][k] - pivot*R[i][k];
					}
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int m = i+1; m < r; ++m) 
					{ 
						if (R[m][j] == 0)
						{
							
						}
						else if(R[m][j] !=0)
						{
							double pivot = R[i][j]/R[m][j];
							for(int k=j; k<c; ++k)
							{
								R[m][k] = R[m][k] - pivot*R[i][k];
							}
							
						}
					}
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					
				}
			}
		}
	}
	int rank = 0;
	vector<vector<double>> rowspace;
		
	// Find the basis for the row space of A
	for (int i = 0; i < r; ++i) 
	{	
		if (R[i][i] != 0)
		{
			rank = rank + 1;
			rowspace.push_back(getRow(R,i));
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					rank = rank + 1;
					rowspace.push_back(getRow(R,i));

					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					rank = rank;
				}
			}
		}
	}
	cout << "\nR:" << endl;
	printMatrix(R);
	//cout << "\nrank(A):" << rank << endl;
	cout << "\nBasis for the row space of the matrix A:" << endl;
	printMatrix(rowspace);
	cout << "\n****************************************************************************************************************************" << endl;
	
}

void columnspacebasis(vector<vector<double>> &A) 
{
	int r = A.size();
	int c = A[0].size();
	vector<vector<double>> R(r, vector<double>(c));
	
	for(int i = 0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			R[i][j] = A[i][j];
		}		
	}
	cout << "\n****************************************************************************************************************************" << endl;
	cout << "\nA:" << endl;
	printMatrix(A);
	
		
	// Forward Elimination
	for (int i = 0; i < r; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < r; ++k) 
		{
			if (abs(R[k][i]) > abs(R[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(R[i], R[pivotRow]);

		

		// Eliminate elements below the pivot
		for (int k = i + 1; k < r; ++k) 
		{
			if (R[i][i] != 0)
			{
				double factor = R[k][i] / R[i][i];
				
				for (int j = i; j < c; ++j) 
				{ 
					R[k][j] -= factor * R[i][j];
				
				}
			}
			else if (R[i][i] ==0)
			{
				i=i+1;
			}
		}
	}

	// Eliminate row that is a linear combination of other row
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if(R[j][i] != 0)
				{
					double pivot = divisiond(R[i][i],R[j][i]);
					for(int k=i; k<c; ++k)
					{
						R[j][k] = pivot*R[j][k] - R[i][k];
					}
				}
				else if(R[j][i] == 0)
				{
					i=i+1;
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i+1; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int k = i+1; k < r; ++k) 
					{ 
						if(R[k][j] != 0)
						{
							double pivot = R[i][j]/R[k][j];
							
							for(int m=i; m<c; ++m)
							{
								R[k][m] = pivot*R[k][m] - R[i][m];
							}
						}
						else if(R[k][j] == 0)
						{
							i=i+1;
						}
					}
				}
				else if (R[i][j] == 0)
				{
					i = i+1;
				}
			}
		}
	}
	
		
	// eliminate -1 in the leading 1 at every row
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] < 0 || R[i][i] == -1)
		{
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = -1* R[i][j];
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i+1; j < c; ++j) 
			{ 
				if (R[i][j] < 0 || R[i][j] == -1)
				{
					for(int k=j; k<c; k++)
					{
						R[i][k] = -1* R[i][k];
					}
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					i = i+1;
				}
			}
		}
	}
	
	// make the leading 1
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			double pivot = R[i][i];
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = R[i][j]/pivot;
			}	
		}
		else if(R[i][i] ==0)
		{
			for (int j = i+1; j < c; ++j) 
			{
				double pivot = R[i][j];
				 if (R[i][j] != 0)
				{
					
					for (int k = j; k < c; ++k) 
					{
						R[i][k] = R[i][k]/pivot;
					}				
					j = c-1;
				}	
				else if (R[i][j] ==0)
				{
					i = i+1;
				}
			}
		}
	}
	
	// Another forward elimination
	for (int i = 0; i < r-1; ++i) 
	{	
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if (R[j][i] == 0)
				{
					
				}
				else if(R[j][i] !=0)
				{
					double pivot = R[i][i]/R[j][i];
					for(int k=i; k<c; ++k)
					{
						R[j][k] = R[j][k] - pivot*R[i][k];
					}
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int m = i+1; m < r; ++m) 
					{ 
						if (R[m][j] == 0)
						{
							
						}
						else if(R[m][j] !=0)
						{
							double pivot = R[i][j]/R[m][j];
							for(int k=j; k<c; ++k)
							{
								R[m][k] = R[m][k] - pivot*R[i][k];
							}
							
						}
					}
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					
				}
			}
		}
	}
	int rank = 0;
	vector<vector<double>> columnspace;
		
	// Find the basis for the column space of A
	for (int i = 0; i < r; ++i) 
	{	
		if (R[i][i] != 0)
		{
			rank = rank + 1;
			columnspace.push_back(getColumn(A,i));
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					rank = rank + 1;
					columnspace.push_back(getColumn(A,j));

					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					rank = rank;
				}
			}
		}
	}
	cout << "\nR:" << endl;
	printMatrix(R);
	
	vector<vector<double>> columnspace_transpose = transpose(columnspace);
	
	//cout << "\nrank(A):" << rank << endl;
	cout << "\nBasis for the column space of the matrix A:" << endl;
	printMatrix(columnspace_transpose);
	cout << "\n****************************************************************************************************************************" << endl;
	
}

void homogeneouslinearsystembasis(vector<vector<double>> &A) 
{
	int r = A.size();
	int c = A[0].size();
	vector<vector<double>> R(r, vector<double>(c));
	
	for(int i = 0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			R[i][j] = A[i][j];
		}		
	}
	cout << "\n****************************************************************************************************************************" << endl;
	cout << "\nA:" << endl;
	printMatrix(A);
	
		
	// Forward Elimination
	for (int i = 0; i < r; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < r; ++k) 
		{
			if (abs(R[k][i]) > abs(R[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(R[i], R[pivotRow]);

		

		// Eliminate elements below the pivot
		for (int k = i + 1; k < r; ++k) 
		{
			if (R[i][i] != 0)
			{
				double factor = R[k][i] / R[i][i];
				
				for (int j = i; j < c; ++j) 
				{ 
					R[k][j] -= factor * R[i][j];
				
				}
			}
			else if (R[i][i] ==0)
			{
				i=i+1;
			}
		}
	}

	// Eliminate row that is a linear combination of other row
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if(R[j][i] != 0)
				{
					double pivot = divisiond(R[i][i],R[j][i]);
					for(int k=i; k<c; ++k)
					{
						R[j][k] = pivot*R[j][k] - R[i][k];
					}
				}
				else if(R[j][i] == 0)
				{
					i=i+1;
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i+1; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int k = i+1; k < r; ++k) 
					{ 
						if(R[k][j] != 0)
						{
							double pivot = R[i][j]/R[k][j];
							
							for(int m=i; m<c; ++m)
							{
								R[k][m] = pivot*R[k][m] - R[i][m];
							}
						}
						else if(R[k][j] == 0)
						{
							i=i+1;
						}
					}
				}
				else if (R[i][j] == 0)
				{
					i = i+1;
				}
			}
		}
	}
	
	// make the leading 1
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			double pivot = R[i][i];
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = R[i][j]/pivot;
			}	
		}
		else if(R[i][i] ==0)
		{
			for (int j = i+1; j < c; ++j) 
			{
				double pivot = R[i][j];
				 if (R[i][j] != 0)
				{	
					for (int k = j; k < c; ++k) 
					{
						R[i][k] = R[i][k]/pivot;
					}				
					j = c-1;
				}	
				else if (R[i][j] ==0)
				{
					
				}
			}
		}
	}

	// make zeros above all leading 1
	for (int i = 0; i < r-1; ++i) 
	{
		if (R[i][i+1] != 0 )
		{
			double pivot = R[i][i+1];
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = R[i][j] - (pivot * R[i+1][j]);
			}
		}
	}
	
	
	// Another forward elimination
	for (int i = 0; i < r-1; ++i) 
	{	
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if (R[j][i] == 0)
				{
					
				}
				else if(R[j][i] !=0)
				{
					double pivot = R[i][i]/R[j][i];
					for(int k=i; k<c; ++k)
					{
						R[j][k] = R[j][k] - pivot*R[i][k];
					}
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int m = i+1; m < r; ++m) 
					{ 
						if (R[m][j] == 0)
						{
							
						}
						else if(R[m][j] !=0)
						{
							double pivot = R[i][j]/R[m][j];
							for(int k=j; k<c; ++k)
							{
								R[m][k] = R[m][k] - pivot*R[i][k];
							}
							
						}
					}
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					
				}
			}
		}
	}

	// Computing rank and nullity
	int rank = 0;
	int dim = c;
		
	// Find the basis for the row space of A
	for (int i = 0; i < r; ++i) 
	{	
		if (R[i][i] != 0)
		{
			rank = rank + 1;
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					rank = rank + 1;
					
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					rank = rank;
				}
			}
		}
	}
	int null = dim-rank;
	// Compute the basis for the solution space of Ax = 0
	vector<vector<double>> solutionbasis;
	vector<vector<double>> identitybasis(null,vector<double>(null,0.0));

	for(int i = 0; i < null; ++i)
	{
		identitybasis[i][i] = 1.0;	
	}

	for(int i = 0; i < rank; ++i)
	{
		solutionbasis.push_back(getRow(R,i));	
	}
	// Multiply by -1
	solutionbasis = scalarmultiplication(solutionbasis,-1);

	// Delete column 1,2, ..., rank(A)
	for(int i = 0; i < rank; ++i)
	{
		// Iterate through each row and erase the first column with the amount of rank(A)
		for (auto& row : solutionbasis) 
		{
			row.erase(row.begin() + 0);
		}	
	}
	//cout << "\nIdentity basis:" << endl;
	//printMatrix(identitybasis);
	for(int i = 0; i < null; ++i)
	{
		solutionbasis.push_back(getRow(identitybasis,i));	
	}

	cout << "\nRank(A) = " << rank << endl;
	cout << "nullity(A) = " << null << endl;
	cout << "dim(A) = " << dim << endl;
	cout << "\nR:" << endl;
	printMatrix(R);
	cout << "\nBasis for the solution space of Ax = 0 :" << endl;
	printMatrix(solutionbasis);
	//vector<vector<double>> columnspace_transpose = transpose(columnspace);
	
	//cout << "\nrank(A):" << rank << endl;
	cout << "\n****************************************************************************************************************************" << endl;
	
}


vector<vector<double>> normalizehomogeneouslinearsystembasis_matrix(vector<vector<double>> &A) 
{
	int r = A.size();
	int c = A[0].size();
	vector<vector<double>> R(r, vector<double>(c));
	vector<vector<double>> result;
	
	for(int i = 0; i<r; i++)
	{
		for(int j=0; j<c; j++)
		{
			R[i][j] = A[i][j];
		}		
	}
		
	// Forward Elimination
	for (int i = 0; i < r; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < r; ++k) 
		{
			if (abs(R[k][i]) > abs(R[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(R[i], R[pivotRow]);

		

		// Eliminate elements below the pivot
		for (int k = i + 1; k < r; ++k) 
		{
			if (R[i][i] != 0)
			{
				double factor = R[k][i] / R[i][i];
				
				for (int j = i; j < c; ++j) 
				{ 
					R[k][j] -= factor * R[i][j];
				
				}
			}
			else if (R[i][i] ==0)
			{
				i=i+1;
			}
		}
	}

	// Eliminate row that is a linear combination of other row
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if(R[j][i] != 0)
				{
					double pivot = divisiond(R[i][i],R[j][i]);
					for(int k=i; k<c; ++k)
					{
						R[j][k] = pivot*R[j][k] - R[i][k];
					}
				}
				else if(R[j][i] == 0)
				{
					i=i+1;
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i+1; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int k = i+1; k < r; ++k) 
					{ 
						if(R[k][j] != 0)
						{
							double pivot = R[i][j]/R[k][j];
							
							for(int m=i; m<c; ++m)
							{
								R[k][m] = pivot*R[k][m] - R[i][m];
							}
						}
						else if(R[k][j] == 0)
						{
							i=i+1;
						}
					}
				}
				else if (R[i][j] == 0)
				{
					i = i+1;
				}
			}
		}
	}
	
	// make the leading 1
	for (int i = 0; i < r; ++i) 
	{
		if (R[i][i] != 0)
		{
			double pivot = R[i][i];
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = R[i][j]/pivot;
			}	
		}
		else if(R[i][i] ==0)
		{
			for (int j = i+1; j < c; ++j) 
			{
				double pivot = R[i][j];
				 if (R[i][j] != 0)
				{	
					for (int k = j; k < c; ++k) 
					{
						R[i][k] = R[i][k]/pivot;
					}				
					j = c-1;
				}	
				else if (R[i][j] ==0)
				{
					
				}
			}
		}
	}

	// make zeros above all leading 1
	for (int i = 0; i < r-1; ++i) 
	{
		if (R[i][i+1] != 0 )
		{
			double pivot = R[i][i+1];
			for (int j = i; j < c; ++j) 
			{ 
				R[i][j] = R[i][j] - (pivot * R[i+1][j]);
			}
		}
	}
	
	
	// Another forward elimination
	for (int i = 0; i < r-1; ++i) 
	{	
		if (R[i][i] != 0)
		{
			for (int j = i+1; j < r; ++j) 
			{ 
				if (R[j][i] == 0)
				{
					
				}
				else if(R[j][i] !=0)
				{
					double pivot = R[i][i]/R[j][i];
					for(int k=i; k<c; ++k)
					{
						R[j][k] = R[j][k] - pivot*R[i][k];
					}
				}
			}
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					for (int m = i+1; m < r; ++m) 
					{ 
						if (R[m][j] == 0)
						{
							
						}
						else if(R[m][j] !=0)
						{
							double pivot = R[i][j]/R[m][j];
							for(int k=j; k<c; ++k)
							{
								R[m][k] = R[m][k] - pivot*R[i][k];
							}
							
						}
					}
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					
				}
			}
		}
	}

	// Computing rank and nullity
	int rank = 0;
	int dim = c;
		
	// Find the basis for the row space of A
	for (int i = 0; i < r; ++i) 
	{	
		if (R[i][i] != 0)
		{
			rank = rank + 1;
		}
		else if (R[i][i] == 0)
		{
			for (int j = i; j < c; ++j) 
			{ 
				if (R[i][j] != 0)
				{
					rank = rank + 1;
					
					j = c-1;
				}
				else if (R[i][j] == 0)
				{
					rank = rank;
				}
			}
		}
	}
	int null = dim-rank;
	// Compute the basis for the solution space of Ax = 0
	vector<vector<double>> solutionbasis;
	vector<vector<double>> identitybasis(null,vector<double>(null,0.0));

	for(int i = 0; i < null; ++i)
	{
		identitybasis[i][i] = 1.0;	
	}

	for(int i = 0; i < rank; ++i)
	{
		solutionbasis.push_back(getRow(R,i));	
	}
	// Multiply by -1
	solutionbasis = scalarmultiplication(solutionbasis,-1);

	// Delete column 1,2, ..., rank(A)
	for(int i = 0; i < rank; ++i)
	{
		// Iterate through each row and erase the first column with the amount of rank(A)
		for (auto& row : solutionbasis) 
		{
			row.erase(row.begin() + 0);
		}	
	}
	for(int i = 0; i < null; ++i)
	{
		solutionbasis.push_back(getRow(identitybasis,i));
	}
	// normalize the solution basis
	for(int i = 0; i < null; ++i)
	{	
		vector<double> vector_v = getColumn(solutionbasis,i);
		double ssv = 0; // sum of squares of the new orthonormal basis.
		for (int j = 0; j < c ; ++j)
		{
			ssv += vector_v[j]*vector_v[j];
		}
		double norm = std::sqrt(ssv);
		// Normalize v (v = v / norm)
		vector_v = scalarmultiplication(vector_v,divisiond(1,norm));
		
		result.push_back(vector_v);
	}
	return result;	 // as row vectors
}


vector<double> reflection_yaxis(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][0] = -1;

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> reflection_xaxis(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][1] = -1;

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> reflection_linex(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 0;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 1;
			}
		}
	}
	
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> orthogonalprojection_xaxis(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][1] = 0;

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> orthogonalprojection_yaxis(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][0] = 0;

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> rotation_ccw(vector<double> &vectorx, double angle) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = cosf(angle);
			}
			else if (i != j )
			{
				standardMatrix[i][j] = sinf(angle);
			}
		}
	}
	standardMatrix[0][1] = -sinf(angle);

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> contractiondilation(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	// works for all R^n

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = k;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> compressionexpansion_xdirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][0] = k;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> compressionexpansion_ydirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][1] = k;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> shear_xdirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}
	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][1] = k;

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> shear_ydirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	if (n != 2) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 2." << std::endl;
		return {};
	}
	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][0] = k;

	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> reflection_xyplane(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[2][2] = -1;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> reflection_xzplane(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][1] = -1;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> reflection_yzplane(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][0] = -1;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> orthogonalprojection_xyplane(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[2][2] = 0;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> orthogonalprojection_xzplane(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][1] = 0;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> orthogonalprojection_yzplane(vector<double> &vectorx) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][0] = 0;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> rotation3d_xaxis_ccw(vector<double> &vectorx, double angle) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	if (vectorx[0] >= 0) // rotation about positive x-axis
	{
		for(int i = 0; i < n; i++)
		{
			for(int j=0; j < n; j++)
			{	
				if (i == j )
				{
					standardMatrix[i][j] = 1;
				}
				else if (i != j )
				{
					standardMatrix[i][j] = 0;
				}
			}
		}
		standardMatrix[1][1] = cosf(angle);
		standardMatrix[2][2] = cosf(angle);
		standardMatrix[1][2] = -sinf(angle);
		standardMatrix[2][1] = sinf(angle);
	}
	else if (vectorx[0] < 0) // rotation about negative x-axis
	{
		for(int i = 0; i < n; i++)
		{
			for(int j=0; j < n; j++)
			{	
				if (i == j )
				{
					standardMatrix[i][j] = 1;
				}
				else if (i != j )
				{
					standardMatrix[i][j] = 0;
				}
			}
		}
		standardMatrix[1][1] = cosf(angle);
		standardMatrix[2][2] = cosf(angle);
		standardMatrix[1][2] = sinf(angle);
		standardMatrix[2][1] = -sinf(angle);
	}
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> rotation3d_yaxis_ccw(vector<double> &vectorx, double angle) 
{
	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	if (vectorx[1] >= 0) // rotation about positive y-axis
	{
		for(int i = 0; i < n; i++)
		{
			for(int j=0; j < n; j++)
			{	
				if (i == j )
				{
					standardMatrix[i][j] = 1;
				}
				else if (i != j )
				{
					standardMatrix[i][j] = 0;
				}
			}
		}
		standardMatrix[0][0] = cosf(angle);
		standardMatrix[0][2] = sinf(angle);
		standardMatrix[2][0] = -sinf(angle);
		standardMatrix[2][2] = cosf(angle);
	}
	else if (vectorx[1] < 0) // rotation about negative y-axis
	{
		for(int i = 0; i < n; i++)
		{
			for(int j=0; j < n; j++)
			{	
				if (i == j )
				{
					standardMatrix[i][j] = 1;
				}
				else if (i != j )
				{
					standardMatrix[i][j] = 0;
				}
			}
		}
		standardMatrix[0][0] = cosf(angle);
		standardMatrix[0][2] = -sinf(angle);
		standardMatrix[2][0] = sinf(angle);
		standardMatrix[2][2] = cosf(angle);
	} 
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> rotation3d_zaxis_ccw(vector<double> &vectorx, double angle) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	if (vectorx[2] >= 0) // rotation about positive z-axis
	{
		for(int i = 0; i < n; i++)
		{
			for(int j=0; j < n; j++)
			{	
				if (i == j )
				{
					standardMatrix[i][j] = 1;
				}
				else if (i != j )
				{
					standardMatrix[i][j] = 0;
				}
			}
		}
		standardMatrix[0][0] = cosf(angle);
		standardMatrix[0][1] = -sinf(angle);
		standardMatrix[1][0] = sinf(angle);
		standardMatrix[1][1] = cosf(angle);
	}
	if (vectorx[2] < 0) // rotation about negative z-axis
	{
		for(int i = 0; i < n; i++)
		{
			for(int j=0; j < n; j++)
			{	
				if (i == j )
				{
					standardMatrix[i][j] = 1;
				}
				else if (i != j )
				{
					standardMatrix[i][j] = 0;
				}
			}
		}
		standardMatrix[0][0] = cosf(angle);
		standardMatrix[0][1] = sinf(angle);
		standardMatrix[1][0] = -sinf(angle);
		standardMatrix[1][1] = cosf(angle);
	}
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> compressionexpansion3d_xdirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[0][0] = k;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> compressionexpansion3d_ydirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[1][1] = k;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<double> compressionexpansion3d_zdirection(vector<double> &vectorx, double k) 
{

	if (vectorx.empty())
	{
        	// Handle empty vectors or invalid dimensions
		return {};
	}

	int n = vectorx.size();
	
	if (n != 3) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimension has to be 3." << std::endl;
		return {};
	}

	// Initialize result vector with appropriate dimensions
	vector<double> result(n);
	// Initialize the standard matrix for reflection about the y-axis
	vector<vector<double>> standardMatrix(n, vector<double>(n));

	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{	
			if (i == j )
			{
				standardMatrix[i][j] = 1;
			}
			else if (i != j )
			{
				standardMatrix[i][j] = 0;
			}
		}
	}
	standardMatrix[2][2] = k;
	// Perform matrix-vector multiplication
	for (int i = 0; i < n; i++) 
	{
		for(int j=0; j < n; j++)
		{
			result[i] += standardMatrix[i][j]*vectorx[j] ;	
		}
	}
	return result;
}

vector<vector<double>> penrose(vector<double>& vectorx, vector<double>& vectory) 
{
	int n = vectorx.size();
	int m = vectory.size();
	if (n != m) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimensions has to be the same." << std::endl;
		return {};
	}
	vector<vector<double>> penrosematrix(n, vector<double>(n, 0.0)); // initialize penrose matrix with size n x n and entries all zero
	double denom = dot(vectorx,vectorx);
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			penrosematrix[i][j] = vectory[i]*vectorx[j];
		}
	}
	penrosematrix = scalarmultiplication(penrosematrix,divisiond(1,denom));

	return penrosematrix;
}

vector<vector<complex<double>>> penroseComplex(vector<complex<double>>& vectorx, vector<complex<double>>& vectory) 
{
	int n = vectorx.size();
	int m = vectory.size();
	if (n != m) 
	{
		// Dimensions are incompatible for multiplication
		cerr << "Error: Vector dimensions has to be the same." << std::endl;
		return {};
	}
	vector<vector<complex<double>>> penrosematrix;
	penrosematrix.assign(n, std::vector<complex<double>>(n, 0.0)); // initialize penrose matrix with size n x n and entries all zero
	complex<double> denom = complexdotproduct(vectorx,vectorx);
	complex<double> z(1,0);
	denom = z/denom;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			penrosematrix[i][j] = vectory[i]*vectorx[j];
		}
	}
	penrosematrix = complexnumbermultiplicationComplexMatrix(penrosematrix,denom);

	return penrosematrix;
}

vector<double> MarkovChain(vector<vector<double>>& P, vector<double>& x0, int k) 
{
	int n = P.size(); // Rows	
	int m = n;
	int nv = x0.size(); // Rows	
	if (n != m) 
	{
		cerr << "Error: Matrix is not square." << std::endl;
		return {};
	}
	if (n != nv) 
	{
		cerr << "Error: Vector size is not compatible with matrix size." << std::endl;
		return {};
	}
	cout << "\nP = " << endl;
	printMatrix(P);
	cout << "\nx_{0} = " << endl;
	printVector(x0);
	
	vector<vector<double>> Pk(n, vector<double>(n));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			Pk[i][j] = P[i][j];
		}
	}
	for (int i = 0; i < k; ++i) 
	{		
		Pk = multiply(Pk,P);	
	}

	vector<double> steadystate_vector = multiplymatrixvector(Pk,x0);	

	cout << "\nState vector, k = "<< k<< " : "<< endl;
	printVector(steadystate_vector);

	return steadystate_vector; 
}

vector<vector<double>> gramschmidt(vector<vector<double>> &A)
{
	int n = A.size();
	int c = A[0].size();
	vector<vector<double>> orthogonal_basis;

	cout << "\nA:" << endl;
	printMatrix(A);
	cout << endl;
	vector<double> v1 = getColumn(A,0);
	orthogonal_basis.push_back(v1);
		
	for (int i =1; i < c ; ++i)
	{
		if(i <= 1)
		{
			vector<double> u = getColumn(A,i);
			vector<double> v1 = getRow(orthogonal_basis,i-1);
			double num = 0;
			double denom = 0;
			for (int l = 0 ; l< n; ++l)
			{
				num += u[l]*v1[l];
				denom += v1[l]*v1[l];
			}
			double k = divisiond(num,denom);
			vector<double> v_process = scalarmultiplication(v1,k);
			
			vector<double> v_new = subtract(u,v_process);
			
			orthogonal_basis.push_back(v_new);
		}
		if(i > 1)
		{
			vector<double> u = getColumn(A,i);
			vector<vector<double>> proj_matrix_final;
			for (int j = 0; j< i; ++j)
			{
				vector<double> v_now = getRow(orthogonal_basis,j);	
				double num = 0;
				double denom = 0;
				for (int l = 0 ; l< n; ++l)
				{
					num += u[l]*v_now[l];
					denom += v_now[l]*v_now[l];
				}
				double k_now = divisiond(num,denom);
				
				vector<double> v_processnow = scalarmultiplication(v_now,k_now);
				
				proj_matrix_final.push_back(v_processnow);
			}
			vector<double> v_total;
			for (int k = 0; k < n; ++k)
			{
				vector<double> v_col =getColumn(proj_matrix_final,k); 
				double sum_col = std::accumulate(v_col.begin(), v_col.end(), 0.0);
				v_total.push_back(sum_col);
			}
			vector<double> v_new2 = subtract(u,v_total);
			orthogonal_basis.push_back(v_new2);
		}		
	}
	vector<vector<double>> orthogonal_basis_final = transpose(orthogonal_basis); 
	cout << "\nOrthogonal basis:" << endl;
	printMatrix(orthogonal_basis_final);

	vector<vector<double>> orthonormal_basis;
	for (int i = 0; i < c; ++i)
	{
		vector<double> basis = getRow(orthogonal_basis,i);
		double norm_denom = norm(basis);
		vector<double> vector_orthonormal = scalarmultiplication(basis,divisiond(1,norm_denom));
		orthonormal_basis.push_back(vector_orthonormal);
	}
	vector<vector<double>> orthonormal_basis_final = transpose(orthonormal_basis); 
	
	cout << "\nOrthonormal basis:" << endl;
	printMatrix(orthonormal_basis_final);

	return orthonormal_basis_final;
}

void QRDecomposition(vector<vector<double>> &A, vector<vector<double>> &Q, vector<vector<double>> &R)
{
	int n = A.size();
	int c = A[0].size();
	vector<vector<double>> orthogonal_basis;

	cout << "\nA:" << endl;
	printMatrix(A);
	cout << endl;
	vector<double> v1 = getColumn(A,0);
	orthogonal_basis.push_back(v1);
		
	for (int i =1; i < c ; ++i)
	{
		if(i <= 1)
		{
			vector<double> u = getColumn(A,i);
			vector<double> v1 = getRow(orthogonal_basis,i-1);
			double num = 0;
			double denom = 0;
			for (int l = 0 ; l< n; ++l)
			{
				num += u[l]*v1[l];
				denom += v1[l]*v1[l];
			}
			double k = divisiond(num,denom);
			vector<double> v_process = scalarmultiplication(v1,k);
			
			vector<double> v_new = subtract(u,v_process);
			
			orthogonal_basis.push_back(v_new);
		}
		if(i > 1)
		{
			vector<double> u = getColumn(A,i);
			vector<vector<double>> proj_matrix_final;
			for (int j = 0; j< i; ++j)
			{
				vector<double> v_now = getRow(orthogonal_basis,j);	
				double num = 0;
				double denom = 0;
				for (int l = 0 ; l< n; ++l)
				{
					num += u[l]*v_now[l];
					denom += v_now[l]*v_now[l];
				}
				double k_now = divisiond(num,denom);
				
				vector<double> v_processnow = scalarmultiplication(v_now,k_now);
				
				proj_matrix_final.push_back(v_processnow);
			}
			vector<double> v_total;
			for (int k = 0; k < n; ++k)
			{
				vector<double> v_col =getColumn(proj_matrix_final,k); 
				double sum_col = std::accumulate(v_col.begin(), v_col.end(), 0.0);
				v_total.push_back(sum_col);
			}
			vector<double> v_new2 = subtract(u,v_total);
			orthogonal_basis.push_back(v_new2);
		}		
	}
	vector<vector<double>> orthogonal_basis_final = transpose(orthogonal_basis); 
	
	vector<vector<double>> orthonormal_basis;
	for (int i = 0; i < c; ++i)
	{
		vector<double> basis = getRow(orthogonal_basis,i);
		double norm_denom = norm(basis);
		vector<double> vector_orthonormal = scalarmultiplication(basis,divisiond(1,norm_denom));
		orthonormal_basis.push_back(vector_orthonormal);
	}
	vector<vector<double>> orthonormal_basis_final = transpose(orthonormal_basis); 
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			Q[i][j] = orthonormal_basis_final[i][j];
		}
	}
	//vector<vector<double>> R_final(n, vector<double>(c,0.0));
	for (int i = 0; i < n; ++i)
	{
		vector<double> q = getColumn(Q,i);
		for (int j = 0; j < c; ++j)
		{
			vector<double> u_R = getColumn(A,j);
			double inner_product = 0;
				
			for (int l = 0 ; l< n; ++l)
			{
				inner_product += u_R[l]*q[l];
			}
			R[i][j] = inner_product;
		}
	}

	
	cout << "\nQ:" << endl;
	printMatrix(Q);

	cout << "\nR:" << endl;
	printMatrix(R);

}

// Function to perform QR decomposition using Modified Gram-Schmidt for complex matrix
void QRDecompositionComplex(vector<vector<complex<double>>>& A, vector<vector<complex<double>>>& Q, vector<vector<complex<double>>>& R) 
{
	int m = A.size(); // Rows
	int n = A[0].size(); // Columns

	Q.assign(m, std::vector<complex<double>>(n, 0.0));
	R.assign(n, std::vector<complex<double>>(n, 0.0));

	for (int i = 0; i < n; ++i) 
	{
		vector<complex<double>> a_i = getcolumnComplexMatrix(A, i);
		vector<complex<double>> u_i = a_i;

		for (int j = 0; j < i; ++j) 
		{
			vector<complex<double>> q_j = getcolumnComplexMatrix(Q, j);
			// r[j][i] = dot(q_j, a_i)
			R[j][i] = complexdotproduct(q_j, a_i);
			// u_i = u_i - r[j][i] * q_j
			u_i = subtractComplexVectors(u_i, complexnumbermultiplicationComplexVector(q_j, R[j][i]));
		}

		// Compute norm of u_i (which becomes the diagonal element R[i][i])
		double norm_u_i = complexnorm(u_i);
		R[i][i] = norm_u_i;
	
		// Normalize u_i to get q_i: q_i = u_i / norm_u_i
		for (int k = 0; k < m; ++k) 
		{
			Q[k][i] = u_i[k] / norm_u_i;
		}
	}
}

// Function to perform QR decomposition using Householder Reflections
void QRDecompositionComplexHouseholder(vector<vector<complex<double>>>& A, vector<vector<complex<double>>>& Q, vector<vector<complex<double>>>& R) 
{
	int m = A.size(); // Rows
	int n = A[0].size(); // Columns

	Q =createIdentityComplexMatrix(n);
	R.assign(n, std::vector<complex<double>>(n, 0.0));

	int iter;
	if (m == n) // iterations is the number of column-1 of matrix A if m=n
	{
		iter = m-1;		
	}
	if (m < n) // iterations depend on the number of column of matrix A if m<n
	{
		iter = m;	
	}	
	for (int i = 0; i < iter; ++i) 
	{		
	
		// compute the next vector
		vector<complex<double>> x;  
		vector<complex<double>> w(n, 0.0);  
		vector<complex<double>> v;  
			
		for (int j = 0; j < n ; ++j) 
		{		
			x.push_back(A[j][i]);
		}
		if (i > 0)
		{
			for (int j = 0; j < i ; ++j) 
			{		
				x[j] = 0;
			}
			
		}
		double real_x= real(x[i]); // check again should we use x[i] or x[0], because x[0]  will always be 0 after iteration 1.
		double norm_x = complexnorm(x); 
		// w = (sgn) ||x_{i}|| e_{i}
		if(real_x < 0 )
		{
			w[i] = norm_x;
		}
		else if(real_x >= 0 )
		{
			w[i] = -norm_x;
		}
		// v = x_{i} - ||x_{i}|| e_{i}
		v = subtractComplexVectors(x,w);

		vector<vector<complex<double>>> mat_v;
		vector<vector<complex<double>>> mat_v_transpose;
		mat_v.assign(n, std::vector<complex<double>>(1, 0.0));
		mat_v_transpose.assign(1, std::vector<complex<double>>(n, 0.0));

		for (int i = 0; i < n ; ++i) 
		{		
			mat_v[i][0] = v[i];
			mat_v_transpose[0][i] = v[i];
		}

		vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
		vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		vector<vector<complex<double>>> H_hat; 
		vector<vector<complex<double>>> R_temp; 
		
		double norm_v = divisiond(1,real(vTv[0][0]));

		vector<vector<complex<double>>> P;
		P = scalarmultiplicationComplexMatrix(vvT,norm_v);
		
		P = scalarmultiplicationComplexMatrix(P,2);
		H_hat = subtractComplexMatrices(I,P);
		
		R_temp = multiplyComplexMatrices(H_hat,A);

		// Rounding down those with small decimal number to 0
		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < m ; ++j)
			{
				double real_value = real(R_temp[i][j]) ;
				if(abs(real_value) < 1e-12)
				{
					R_temp[i][j] = 0;
				}
			}
		}

		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < m ; ++j)
			{
				R[i][j] = R_temp[i][j];
			}
		}
		
		Q = multiplyComplexMatrices(Q,H_hat);
		// Rounding down those with small decimal number to 0
		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < m ; ++j)
			{
				double real_value = real(Q[i][j]) ;
				if(abs(real_value) < 1e-12)
				{
					//Q[i][j] = 0;
				}
			}
		}
		
	}
}

void GaussJordanComplexMatrix(vector<vector<complex<double>>> &A)
{
	int n = A.size();
	int m = A[0].size();
	
	if (n != m) 
	{
		throw std::invalid_argument("Matrix is not square.");
	}

	vector<vector<complex<double>>> A_rref;
	vector<vector<complex<double>>> I_complex = createIdentityComplexMatrix(n);
	A_rref.assign(n, std::vector<complex<double>>(m, 0.0));
	
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			A_rref[i][j] = A[i][j];
		}
	}
	cout << "A:" << endl;
	printComplexMatrix(A_rref);
	
	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (moduluscomplex(A_rref[k][i]) > moduluscomplex(A_rref[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(A_rref[i], A_rref[pivotRow]);
		swap(I_complex[i], I_complex[pivotRow]);
		cout << "Swap row\nA_{rref} :" << endl;
		printComplexMatrix(A_rref);
		cout << "A^{-1} :" << endl;
		printComplexMatrix(I_complex);
	
		// Check for singular matrix (no unique solution)
		if (moduluscomplex(A_rref[i][i]) < 1e-12) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
		}
			
		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			complex<double> factor = A_rref[k][i] / A_rref[i][i];
			for (int j = 0; j < n; ++j) 
			{ // Iterate up to n for the constant term
				A_rref[k][j] -= factor * A_rref[i][j];
				I_complex[k][j] -= factor * I_complex[i][j];
			}
			cout << "Eliminate elements below the pivot, i = " << i << ", k = " << k << "\nA_{rref} :" << endl;
			printComplexMatrix(A_rref);
			cout << "A^{-1} :" << endl;
			printComplexMatrix(I_complex);
	
		}
	}
	cout << "After forward elimination: \nA_{rref} :" << endl;
	printComplexMatrix(A_rref);
	cout << "A^{-1} :" << endl;
	printComplexMatrix(I_complex);
	
	// make the leading 1
	for (int i = 0; i < n; ++i) 
	{
		double realpart = lround(real(A_rref[i][i]));
		if (realpart != 1.0 )
		{
			complex<double> pivot = A_rref[i][i];
			for (int j = 0; j < m; ++j) 
			{ 
				A_rref[i][j] = A_rref[i][j]/pivot;
				I_complex[i][j] = I_complex[i][j]/pivot;
			}	
		}
	}
	cout << "After creating leading 1:\nA_{rref} :" << endl;
	printComplexMatrix(A_rref);
	cout << "A^{-1} :" << endl;
	printComplexMatrix(I_complex);

	int f = 2;
	// Backward elimination
	// make zeros above all leading 1 / make matrix A into reduced row echelon form
	for (int i = n-2; i >= 0; --i) 
	{
		for (int k = 1; k < f; ++k)
		{
			double realpart = real(A_rref[i][i+k]);
			//cout << realpart << endl;		
			if (realpart != 0.0  )
			{
				complex<double> pivot = A_rref[i][i+k];
				for (int j = 0; j < m; ++j) 
				{ 
					A_rref[i][j] = A_rref[i][j] - (pivot * A_rref[i+k][j]);
					I_complex[i][j] = I_complex[i][j] - (pivot * I_complex[i+k][j]);
				}
			}
			cout << "Eliminate elements above the pivot, i = " << i << ", k = " << k << "\nA_{rref} :" << endl;
			printComplexMatrix(A_rref);
			cout << "A^{-1} :" << endl;
			printComplexMatrix(I_complex);
		}
	f = f+1;
	}
	
	cout << "After backward elimination: \nA_{rref} :" << endl;
	printComplexMatrix(A_rref);
	cout << "A^{-1} :" << endl;
	printComplexMatrix(I_complex);
	
	cout << "AA^{-1} :" << endl;
	cmat AI = multiplyComplexMatrices(A,I_complex);
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			double realpart = real(A_rref[i][j]);
			if (abs(realpart) < 1e-12) // Using a small epsilon
			{ 
				AI[i][j] = 0;
			}
		}
	}
	printComplexMatrix(AI);
}

// this can be used to compute complex eigenvectors
void GaussJordanComplexMatrixTEST(vector<vector<complex<double>>> &A) //NOT YET
{
	int n = A.size();
	int m = A[0].size();
	
	if (n != m) 
	{
		throw std::invalid_argument("Matrix is not square.");
	}

	vector<vector<complex<double>>> A_ref;
	A_ref.assign(n, std::vector<complex<double>>(m, 0.0));
	
	for (int i = 0; i < n; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			A_ref[i][j] = A[i][j];
		}
	}
	cout << "A:" << endl;
	printComplexMatrix(A_ref);
	
	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (moduluscomplex(A_ref[k][i]) > moduluscomplex(A_ref[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(A_ref[i], A_ref[pivotRow]);
		cout << "Swap row\nA_{ref} :" << endl;
		printComplexMatrix(A_ref);
	
		// Check for singular matrix (no unique solution)
		if (moduluscomplex(A_ref[i][i]) < 1e-12) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
		}
			
		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			complex<double> factor = A_ref[k][i] / A_ref[i][i];
			//cout << "Pivot, a = " << real(factor) << ", b = " << imag(factor) << endl;
			
			for (int j = 0; j < n; ++j) 
			{ // Iterate up to n for the constant term
				A_ref[k][j] -= factor * A_ref[i][j];
			}
			cout << "Eliminate elements below the pivot, i = " << i << ", k = " << k << "\nA_{ref} :" << endl;
			printComplexMatrix(A_ref);
	
		}
	}
	
	// make the leading 1
	for (int i = 0; i < n; ++i) 
	{
		double realpart = lround(real(A_ref[i][i]));
		if (realpart != 1.0 )
		{
			complex<double> pivot = A_ref[i][i];
			for (int j = 0; j < m; ++j) 
			{ 
				A_ref[i][j] = A_ref[i][j]/pivot;
			}	
		}
	}
	cout << "After creating leading 1:\nA_{ref} :" << endl;
	printComplexMatrix(A_ref);

	vector<complex<double>> eigenvector;
	for (int i = 0; i < n; ++i) 
	{
		double realpart = real(A_ref[i][i]);
		if (realpart == 1.0 )
		{
			for (int j = 1; j < m; ++j) 
			{ 
				eigenvector.push_back(-A_ref[i][j]);
			}	
		}
	}
	eigenvector.push_back(complex<double>(1.0,0.0));
	
	cout << "\nEigenvector :" << endl;
	printComplexVector(eigenvector);
}


void gaussianelimination(const vector<vector<double>> &A)
{
	int n = A.size();
	vector<vector<double>> B_mat(n,vector<double>(n));
	vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));

	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = A[i][j];
			B_mat[i][j] = A[i][j];
		}
	}
	cout << "A :" << endl;
	printMatrix(B_mat);
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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

	// Back Substitution
	vector<double> solution(n,1.0);
	for (int i = n - 1; i >= 0; --i) 
	{
		double sum = 0.0;
			
		if(augmentedMatrix[i][i] !=0)
		{
			if (i >=1)
			{
				for (int k = i - 1; k >= 0; --k) 
				{
					if(augmentedMatrix[i][k] == 0)
					{
						solution[i] = (augmentedMatrix[i][n-1] ) / augmentedMatrix[i][i];
					}
					else if(augmentedMatrix[i][k] != 0)
					{
						for(int k = i-1 ; k >=0 ; --k)
						{
							sum += augmentedMatrix[i][k] * solution[k];
						}
						solution[i] = (augmentedMatrix[i][n-1] - sum) / augmentedMatrix[i][i];
					}
					
				}
			}
			else if(i < 1)
			{
				for (int j = i + 1; j < n-1; ++j) 
				{
					sum += augmentedMatrix[i][j] * solution[j];				
				}
				solution[i] = (augmentedMatrix[i][n-1] - sum) / augmentedMatrix[i][i];
			}
		}
		
	}	

	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			B_mat[i][j] = augmentedMatrix[i][j];
		}
	}
	
	cout << "A (in row reduced echelon form) :" << endl;
	printMatrix(B_mat);

	// Print Solution
	cout << "\nSolution:" << endl;
	for (int i = 0; i < n-1; ++i) 
	{
	cout << "x" << i + 1 << " = " << fixed << setprecision(5) << solution[i] << endl;
	}
	cout << endl;

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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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

void solve_nhsystem_resultsonly(vector<vector<double>> &A, vector<vector<double>> &b, vector<double> &x) 
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
		{ 
			
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
		x.push_back(solution[i]);
	}	
	
	// Print Solution
	reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	for (int i = 0; i < n; ++i) 
	{
	//	cout << "x" << i + 1 << " = " << fixed << setprecision(5) << x[i] << endl;
	}

}

void solve_homogeneoussystem(vector<vector<double>> &A, vector<double> &x) // not yet finished
{
	int n = A.size();
	vector<vector<double>> Matrix(n, vector<double>(n));
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			Matrix[i][j] = A[i][j];
			
		}
	}
	cout << "\nMatrix:" << endl;	
	printMatrix(Matrix);

	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (abs(Matrix[k][i]) > abs(Matrix[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(Matrix[i], Matrix[pivotRow]);

		// Check for singular matrix (no unique solution)
		if (abs(Matrix[i][i]) < 1e-12) // Using a small epsilon
		{ 
			//cout << "No unique solution or infinite solutions exist." << endl;
		}

		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			double factor = Matrix[k][i] / Matrix[i][i];
			for (int j = i; j < n; ++j) 
			{ // Iterate up to n for the constant term
				Matrix[k][j] -= factor * Matrix[i][j];
			}
		}
	}

	// make the leading 1
	for (int i = 0; i < n; ++i) 
	{
		if (Matrix[i][i] != 0)
		{
			double pivot = Matrix[i][i];
			for (int j = i; j < n; ++j) 
			{ 
				Matrix[i][j] = Matrix[i][j]/pivot;
			}	
		}
		else if(Matrix[i][i] ==0)
		{
			for (int j = i+1; j < n; ++j) 
			{
				double pivot = Matrix[i][j];
				 if (Matrix[i][j] != 0)
				{	
					for (int k = j; k < n; ++k) 
					{
						Matrix[i][k] = Matrix[i][k]/pivot;
					}				
					j = n-1;
				}	
				else if (Matrix[i][j] ==0)
				{
					
				}
			}
		}
	}

	// make zeros above all leading 1
	for (int i = 0; i < n-1; ++i) 
	{
		if (Matrix[i][i+1] != 0 )
		{
			double pivot = Matrix[i][i+1];
			for (int j = i; j < n; ++j) 
			{ 
				Matrix[i][j] = Matrix[i][j] - (pivot * Matrix[i+1][j]);
			}
		}
	}
	
	cout << "\nMatrix in reduced row form:" << endl;
	printMatrix(Matrix);

	/*// Back Substitution
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
	}	*/
	
	// Print Solution
	cout << "\nSolution:" << endl;
	/*reverse(x.begin(),x.end());	// reverse the vector because we are using back substitution
	for (int i = 0; i < n; ++i) 
	{
		cout << "x" << i + 1 << " = " << fixed << setprecision(5) << x[i] << endl;
	}*/
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
		if (abs(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
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
		if (abs(augmentedMatrix2[i][i]) < 1e-12) // Using a small epsilon
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

SymbolicMatrix IVPSolution_firstorderdiffeq(vector<vector<complex<double>>> &A, vector<complex<double>> &vec_b, double a, double b) 
// Done on Full Moon January 3rd, 2026
{	
	int n = A.size(); // Rows	
	int m = n;

	cout << "\nMatrix A : " << endl;
	printComplexMatrix(A);
	int N = 100;
	int iterations = 100;

	vector<complex<double>> vector_eigenvalues(n);		
	vector<vector<complex<double>>> matrix_eigenvectors;	
	vector<complex<double>> v;
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(n, 0.0));

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	int index = 0;
	//  N =  number of the outermost iteration to obtain the correct eigenvalues.
	for (int iter = 0; iter < N; ++iter) 
	{
		
		// 1. Obtain a seed:
		std::default_random_engine generator(
		std::chrono::system_clock::now().time_since_epoch().count());
		
		// 2. Define the distribution for floating-point numbers (e.g., uniform distribution)
		// Range [a, b)
		std::uniform_real_distribution<double> distribution(a, b);	
		std::uniform_real_distribution<double> distribution2(0, b);
		double mu = distribution(generator); 
		double sigma = distribution2(generator);
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		// mu and sigma are also random number that are generated with uniform distribution
		v = complexvecrand_normal_zeroimaginary(mu, sigma, m); // the eigenvector won't have imaginary parts
		double norm_v = divisiond(1,complexnorm(v));
		v = scalarmultiplicationComplexVector(v,norm_v);
		
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		complex<double> cdotproduct = complexdotproduct(v,v);
		complex<double> lambda = complexquadraticmultiplication(A,v) / cdotproduct;
		
		// Solve (A - λ^{k-1}I)w = v^{k-1}
		// Then multiply the reverse to obtain RQ = A^{k}
		for (int i = 0; i < iterations; ++i) 
		{		
			vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
			vector<vector<complex<double>>> lambda_I = complexnumbermultiplicationComplexMatrix(I,lambda); 
			vector<vector<complex<double>>> A_lambda_I = subtractComplexMatrices(A_original,lambda_I);
			vector<vector<complex<double>>> A_lambda_I_inverse = ComplexMatrixInverse(A_lambda_I);
		
			vector<complex<double>> w = multiplycomplexmatrixvector(A_lambda_I_inverse,v);

			double norm_w = divisiond(1,complexnorm(w));
			w = scalarmultiplicationComplexVector(w,norm_w);
			//cout <<"\nw_{normalize} = (A - λ^{k-1}I) * v^{k-1}: "<<endl;
			//printComplexVector(w);
		
			for (int i = 0; i < n ; ++i) 
			{		
				v[i] = w[i];
			}
		}
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		complex<double> dotproduct_final = complexdotproduct(v,v);
		complex<double> lambda_final = complexquadraticmultiplication(A,v) / dotproduct_final;
		complex<double> eigenvalue(lround(real(lambda_final)), lround(imag(lambda_final)));
		
		vector<complex<double>> A_x =  multiplycomplexmatrixvector(A,v);  
		vector<complex<double>> lambda_x =  complexnumbermultiplicationComplexVector(v,eigenvalue); 
		
		// to obtain the whole eigenvalues and eigenvectors
		vector<complex<double>> difference = subtractComplexVectors(A_x,lambda_x);
		double diffnorm = complexnorm(difference);
		if (abs(diffnorm) < 1e-12 )
		{
			for (int i = 0; i < n ; ++i)		
			{		
				if(eigenvalue != vector_eigenvalues[i])
				{
					if(i == n-1)	
					{
						vector_eigenvalues[index] = eigenvalue;
						matrix_eigenvectors.push_back(v);
						index= index+1;
					}
				}			
				else if(eigenvalue == vector_eigenvalues[i])
				{
					break;
				}			
			}		
		}	
		if (index == n)
		{
			iter = N-1;
		}
		
	}
	cout <<"\nEigenvalues : "<<endl;	
	printComplexVector(vector_eigenvalues);
	cout <<"\nEigenvectors (represented by column): "<<endl;	
	vector<vector<complex<double>>> matrix_eigenvectors_transpose = TransposeComplexMatrix(matrix_eigenvectors);
	printComplexMatrix(matrix_eigenvectors_transpose);
	
	// To obtain the constants c_{1}, c_{2}, ..., c_{n}
	vector<vector<complex<double>>> augmentedMatrix;
	augmentedMatrix.assign(n, std::vector<complex<double>>(n+1, 0.0));
	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = matrix_eigenvectors_transpose[i][j];
			if(j==n)
			{
				augmentedMatrix[i][j] = vec_b[i];
			}
		}
	}
	
	cout << "\nAugmented Matrix:" << endl;	
	printComplexMatrix(augmentedMatrix);


	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (moduluscomplex(augmentedMatrix[k][i]) > moduluscomplex(augmentedMatrix[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);
	
		// Check for singular matrix (no unique solution)
		if (moduluscomplex(augmentedMatrix[i][i]) < 1e-12) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
		}
			
		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			complex<double> factor = augmentedMatrix[k][i] /augmentedMatrix[i][i];
			for (int j = 0; j <= n; ++j) 
			{ // Iterate up to n for the constant term
				augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
			}
		}
	}
	//cout << "After forward elimination: \nA_{rref} :" << endl;
	//printComplexMatrix(augmentedMatrix);
	
	
	// make the leading 1
	for (int i = 0; i < n; ++i) 
	{
		double realpart = real(augmentedMatrix[i][i]);
		if (realpart != 1.0 )
		{
			complex<double> pivot = augmentedMatrix[i][i];
			for (int j = 0; j <= n; ++j) 
			{ 
				augmentedMatrix[i][j] = augmentedMatrix[i][j]/pivot;
			}	
		}
	}
	//cout << "After creating leading 1:\nA_{rref} :" << endl;
	//printComplexMatrix(augmentedMatrix);

	int f = 2;
	// Backward elimination
	// make zeros above all leading 1 / make matrix A into reduced row echelon form
	for (int i = n-2; i >= 0; --i) 
	{
		for (int k = 1; k < f; ++k)
		{
			double realpart = real(augmentedMatrix[i][i+k]);	
			if (realpart != 0.0  )
			{
				complex<double> pivot = augmentedMatrix[i][i+k];
				for (int j = 0; j <= n; ++j) 
				{ 
					augmentedMatrix[i][j] = augmentedMatrix[i][j] - (pivot * augmentedMatrix[i+k][j]);
				}
			}
		}
		f = f+1;
	}
	
	cout << "\nAfter backward elimination:" << endl;
	printComplexMatrix(augmentedMatrix);

	vector<complex<double>> solution;
	for (int i = 0; i <n; ++i) 
	{
		solution.push_back(augmentedMatrix[i][n]);
	}	
	//printComplexVector(solution);

	Matrix<Symbolic> B_mat(n,1);

	Symbolic x("x");
	for(int i = 0; i < n; i++)
	{
		Symbolic sum = 0;
		for(int j = 0; j < n; j++)
		{
			// we take the real part only
			sum += real(solution[j])*real(matrix_eigenvectors_transpose[i][j])*exp(real(vector_eigenvalues[j])*x);	
		}
		B_mat[i][0] = sum;
	}
	cout << "\n The solution satisfying the initial condition is : " << endl;
	cout <<B_mat << endl;
	return B_mat;
}

vector<vector<complex<double>>> Arnoldi(vector<vector<complex<double>>> &A, double a, double b) // NOT YET
// Done on Full Moon + 3: January 6th, 2026
{	
	int n = A.size(); // Rows	
	int m = n;

	cout << "\nMatrix A : " << endl;
	printComplexMatrix(A);
	int iterations = n-1;

	vector<vector<complex<double>>> matrix_basisvectors;
	matrix_basisvectors.assign(n, std::vector<complex<double>>(n, 0.0));	
	vector<complex<double>> q;
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(n, 0.0));

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
		
		// 1. Obtain a seed:
		std::default_random_engine generator(
		std::chrono::system_clock::now().time_since_epoch().count());
		
		// 2. Define the distribution for floating-point numbers (e.g., uniform distribution)
		// Range [a, b)
		std::uniform_real_distribution<double> distribution(a, b);	
		std::uniform_real_distribution<double> distribution2(0, b);
		double mu = distribution(generator); 
		double sigma = distribution2(generator);
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		// mu and sigma are also random number that are generated with uniform distribution
		q= complexvecrand_normal_zeroimaginary(mu, sigma, m); // the vector won't have imaginary parts
		double norm_q = divisiond(1,complexnorm(q));
		q = scalarmultiplicationComplexVector(q,norm_q); // normalize to obtain the first basis vector

		// Populate the first column of the matrix with the first basis vector
		for (int i = 0; i < n ; ++i) 
		{		
			matrix_basisvectors[i][0] = q[i];
		}
		int col = 1;
		for (int i = 0; i < iterations; ++i) 
		{		
			// compute the next vector
			vector<complex<double>> v_k =  multiplycomplexmatrixvector(A,q);  
			// compute projection coefficient
			complex<double> h_coefficient = complexdotproduct(q,v_k);

			// subtract component in q_{j}'s direction
			vector<complex<double>> h_q = complexnumbermultiplicationComplexVector(q,h_coefficient);
			v_k = subtractComplexVectors(v_k,h_q);

			double norm_v = divisiond(1,complexnorm(v_k));
			// normalize to obtain the next orthonormal vector
			vector<complex<double>> new_q = scalarmultiplicationComplexVector(v_k,norm_v); 

			for (int i = 0; i < n ; ++i) 
			{		
				q[i] = new_q[i];
				matrix_basisvectors[i][col] = q[i];
			}

			cout <<"\nNext orthonormal vector : "<<endl;	
			printComplexVector(new_q);
			col = col+1;	
		}
	cout <<"\nMatrix orthonormal vector : "<<endl;	
	printComplexMatrix(matrix_basisvectors);

	vector<vector<complex<double>>> matrix_basisvectors_transpose = TransposeComplexMatrix(matrix_basisvectors);
	cout <<"\nMatrix orthonormal vector transpose: "<<endl;	
	printComplexMatrix(matrix_basisvectors_transpose);

	vector<vector<complex<double>>> AP = multiplyComplexMatrices(A,matrix_basisvectors);
	vector<vector<complex<double>>> PTAP = multiplyComplexMatrices(matrix_basisvectors_transpose,AP);
	cout <<"\nP^{T}AP: "<<endl;	
	printComplexMatrix(PTAP);

	return matrix_basisvectors;
}

vector<vector<complex<double>>> HessenbergDecomposition(vector<vector<complex<double>>> &A) 
// Done on Full Moon + 3: January 6th, 2026
{	
	int n = A.size(); // Rows	
	int m = n;

	if (n != m) 
	{
		cerr << "Error: Matrix is not square." << std::endl;
		return {};
	}
	cout << "\nMatrix A : " << endl;
	printComplexMatrix(A);
	
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(n, 0.0));
	vector<vector<complex<double>>> P_Hessenberg = createIdentityComplexMatrix(n);
	vector<vector<complex<double>>> AP_Hessenberg; 
	vector<vector<complex<double>>> PAP_Hessenberg; 
	vector<vector<complex<double>>> P_Hessenberg_transpose; 

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
		
	int row = 1;	
	int col = 0;
	int n_x = n-1;
	for (int i = 0; i < n-2; ++i) 
	{		
		// compute the next vector
		vector<complex<double>> x;  
		vector<complex<double>> w(n-row, 0.0);  
		vector<complex<double>> v;  
			
		for (int j = row; j <= n_x ; ++j) 
		{		
			x.push_back(A[j][col]);
		}
		double real_x= real(x[0]);
		double norm_x = complexnorm(x); 
		if(real_x < 0 )
		{
			w[0] = norm_x;
		}
		else if(real_x >= 0 )
		{
			w[0] = -norm_x;
		}
		v = subtractComplexVectors(w,x);
		
		cout << "\nx: " << endl;
		printComplexVector(x);
		cout << "\nw: " << endl;
		printComplexVector(w);
		cout << "\nv: " << endl;
		printComplexVector(v);
		vector<vector<complex<double>>> mat_v;
		vector<vector<complex<double>>> mat_v_transpose;
		mat_v.assign(n-row, std::vector<complex<double>>(1, 0.0));
		mat_v_transpose.assign(1, std::vector<complex<double>>(n-row, 0.0));

		for (int i = 0; i < n-row ; ++i) 
		{		
			mat_v[i][0] = v[i];
			mat_v_transpose[0][i] = v[i];
		}

		vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
		vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n-row);
		vector<vector<complex<double>>> H_hat; 
		vector<vector<complex<double>>> AH; 
		vector<vector<complex<double>>> HAH; 
		vector<vector<complex<double>>> H_new = createIdentityComplexMatrix(n);
		
		double norm_v = divisiond(1,real(vTv[0][0]));

		vector<vector<complex<double>>> P;
		P = scalarmultiplicationComplexMatrix(vvT,norm_v);
		
		P = scalarmultiplicationComplexMatrix(P,2);
		H_hat = subtractComplexMatrices(I,P);
		
		int i1 = 0;
		for (int i = row; i < n ; ++i) 
		{	
			int j1 = 0;	
			for (int j = row; j < n ; ++j)
			{
				H_new[i][j] = H_hat[i1][j1];
				
				j1 = j1 + 1;
			}
			i1 = i1+1;
		}

		AH = multiplyComplexMatrices(A,H_new);
		HAH = multiplyComplexMatrices(H_new,AH);

		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < n ; ++j)
			{
				double real_value = real(HAH[i][j]) ;
				if(abs(real_value) < 1e-12)
				{
					HAH[i][j] = 0;
				}
			}
		}

		cout << "\nvTv: " << endl;
		printComplexMatrix(vTv);
		cout << "\nvvT: " << endl;
		printComplexMatrix(vvT);
		cout << "\nP: " << endl;
		printComplexMatrix(P);
		cout << "\nI - 2P: " << endl;
		printComplexMatrix(H_hat);
		cout << "\nH_{i}: " << endl;
		printComplexMatrix(H_new);
		cout << "\nH_{i} A H_{i}: " << endl;
		printComplexMatrix(HAH);

		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < n ; ++j)
			{
				A[i][j] = HAH[i][j];
			}
		}

		
		P_Hessenberg = multiplyComplexMatrices(P_Hessenberg,H_new);
		cout << "\nHouseholder Matrix: " << endl;
		printComplexMatrix(P_Hessenberg);


		row = row + 1;
		col = col+1;	
	}
	
	P_Hessenberg_transpose = TransposeComplexMatrix(P_Hessenberg);

	AP_Hessenberg = multiplyComplexMatrices(A_original,P_Hessenberg);
	PAP_Hessenberg = multiplyComplexMatrices(P_Hessenberg_transpose,AP_Hessenberg);
		
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			double real_value = real(PAP_Hessenberg[i][j]) ;
			if(abs(real_value) < 1e-12)
			{
				PAP_Hessenberg[i][j] = 0;
			}
		}
	}
	cout << "\n****************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\n****************************************************************************" << endl;
	cout << "\nP^{T} A P: " << endl;
	printComplexMatrix(PAP_Hessenberg);

	return PAP_Hessenberg;
}

vector<vector<complex<double>>> HessenbergDecompositionResultOnly(vector<vector<complex<double>>> &A) 
// Done on Full Moon + 3: January 6th, 2026
{	
	int n = A.size(); // Rows	
	int m = n;

	if (n != m) 
	{
		cerr << "Error: Matrix is not square." << std::endl;
		return {};
	}
	
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(n, 0.0));
	vector<vector<complex<double>>> P_Hessenberg = createIdentityComplexMatrix(n);
	vector<vector<complex<double>>> AP_Hessenberg; 
	vector<vector<complex<double>>> PAP_Hessenberg; 
	vector<vector<complex<double>>> P_Hessenberg_transpose; 

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
		
	int row = 1;	
	int col = 0;
	int n_x = n-1;
	for (int i = 0; i < n-2; ++i) 
	{		
		// compute the next vector
		vector<complex<double>> x;  
		vector<complex<double>> w(n-row, 0.0);  
		vector<complex<double>> v;  
			
		for (int j = row; j <= n_x ; ++j) 
		{		
			x.push_back(A[j][col]);
		}
		double real_x= real(x[0]);
		double norm_x = complexnorm(x); 
		if(real_x < 0 )
		{
			w[0] = norm_x;
		}
		else if(real_x >= 0 )
		{
			w[0] = -norm_x;
		}
		v = subtractComplexVectors(w,x);
		
		
		vector<vector<complex<double>>> mat_v;
		vector<vector<complex<double>>> mat_v_transpose;
		mat_v.assign(n-row, std::vector<complex<double>>(1, 0.0));
		mat_v_transpose.assign(1, std::vector<complex<double>>(n-row, 0.0));

		for (int i = 0; i < n-row ; ++i) 
		{		
			mat_v[i][0] = v[i];
			mat_v_transpose[0][i] = v[i];
		}

		vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
		vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n-row);
		vector<vector<complex<double>>> H_hat; 
		vector<vector<complex<double>>> AH; 
		vector<vector<complex<double>>> HAH; 
		vector<vector<complex<double>>> H_new = createIdentityComplexMatrix(n);
		
		double norm_v = divisiond(1,real(vTv[0][0]));

		vector<vector<complex<double>>> P;
		P = scalarmultiplicationComplexMatrix(vvT,norm_v);
		
		P = scalarmultiplicationComplexMatrix(P,2);
		H_hat = subtractComplexMatrices(I,P);
		
		int i1 = 0;
		for (int i = row; i < n ; ++i) 
		{	
			int j1 = 0;	
			for (int j = row; j < n ; ++j)
			{
				H_new[i][j] = H_hat[i1][j1];
				
				j1 = j1 + 1;
			}
			i1 = i1+1;
		}

		AH = multiplyComplexMatrices(A,H_new);
		HAH = multiplyComplexMatrices(H_new,AH);

		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < n ; ++j)
			{
				double real_value = real(HAH[i][j]) ;
				if(abs(real_value) < 1e-12)
				{
					HAH[i][j] = 0;
				}
			}
		}

		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < n ; ++j)
			{
				A[i][j] = HAH[i][j];
			}
		}

		
		P_Hessenberg = multiplyComplexMatrices(P_Hessenberg,H_new);
		
		row = row + 1;
		col = col+1;	
	}
	
	P_Hessenberg_transpose = TransposeComplexMatrix(P_Hessenberg);

	AP_Hessenberg = multiplyComplexMatrices(A_original,P_Hessenberg);
	PAP_Hessenberg = multiplyComplexMatrices(P_Hessenberg_transpose,AP_Hessenberg);
		
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			double real_value = real(PAP_Hessenberg[i][j]) ;
			if(abs(real_value) < 1e-12)
			{
				PAP_Hessenberg[i][j] = 0;
			}
		}
	}
	
	return PAP_Hessenberg;
}


vector<vector<complex<double>>> SpectralDecomposition(vector<vector<complex<double>>> &A, double a, double b) 
// Done on Full Moon + 2: January 5th, 2026
{	
	int n = A.size(); // Rows	
	int m = n;

	cout << "\nMatrix A : " << endl;
	printComplexMatrix(A);
	int N = 100;
	int iterations = 100;

	vector<complex<double>> vector_eigenvalues(n);		
	vector<vector<complex<double>>> matrix_eigenvectors;	
	vector<complex<double>> v;
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(n, 0.0));

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	int index = 0;
	//  N =  number of the outermost iteration to obtain the correct eigenvalues.
	for (int iter = 0; iter < N; ++iter) 
	{
		
		// 1. Obtain a seed:
		std::default_random_engine generator(
		std::chrono::system_clock::now().time_since_epoch().count());
		
		// 2. Define the distribution for floating-point numbers (e.g., uniform distribution)
		// Range [a, b)
		std::uniform_real_distribution<double> distribution(a, b);	
		std::uniform_real_distribution<double> distribution2(0, b);
		double mu = distribution(generator); 
		double sigma = distribution2(generator);
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		// mu and sigma are also random number that are generated with uniform distribution
		v = complexvecrand_normal_zeroimaginary(mu, sigma, m); // the eigenvector won't have imaginary parts
		double norm_v = divisiond(1,complexnorm(v));
		v = scalarmultiplicationComplexVector(v,norm_v);
		
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		complex<double> cdotproduct = complexdotproduct(v,v);
		complex<double> lambda = complexquadraticmultiplication(A,v) / cdotproduct;
		
		// Solve (A - λ^{k-1}I)w = v^{k-1}
		// Rayleigh quotient iteration
		for (int i = 0; i < iterations; ++i) 
		{		
			vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
			vector<vector<complex<double>>> lambda_I = complexnumbermultiplicationComplexMatrix(I,lambda); 
			vector<vector<complex<double>>> A_lambda_I = subtractComplexMatrices(A_original,lambda_I);
			vector<vector<complex<double>>> A_lambda_I_inverse = ComplexMatrixInverse(A_lambda_I);
		
			vector<complex<double>> w = multiplycomplexmatrixvector(A_lambda_I_inverse,v);

			double norm_w = divisiond(1,complexnorm(w));
			w = scalarmultiplicationComplexVector(w,norm_w);
			//cout <<"\nw_{normalize} = (A - λ^{k-1}I) * v^{k-1}: "<<endl;
			//printComplexVector(w);
		
			for (int i = 0; i < n ; ++i) 
			{		
				v[i] = w[i];
			}
		}
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		complex<double> dotproduct_final = complexdotproduct(v,v);
		complex<double> lambda_final = complexquadraticmultiplication(A,v) / dotproduct_final;
		complex<double> eigenvalue(lround(real(lambda_final)), lround(imag(lambda_final)));
		
		vector<complex<double>> A_x =  multiplycomplexmatrixvector(A,v);  
		vector<complex<double>> lambda_x =  complexnumbermultiplicationComplexVector(v,eigenvalue); 
		
		// to obtain the whole eigenvalues and eigenvectors
		vector<complex<double>> difference = subtractComplexVectors(A_x,lambda_x);
		double diffnorm = complexnorm(difference);
		if (abs(diffnorm) < 1e-12 )
		{
			for (int i = 0; i < n ; ++i)		
			{		
				if(eigenvalue != vector_eigenvalues[i])
				{
					if(i == n-1)	
					{
						vector_eigenvalues[index] = eigenvalue;
						matrix_eigenvectors.push_back(v);
						index= index+1;
					}
				}			
				else if(eigenvalue == vector_eigenvalues[i])
				{
					break;
				}			
			}		
		}	
		if (index == n)
		{
			iter = N-1;
		}
		
	}
	cout <<"\nEigenvalues : "<<endl;	
	printComplexVector(vector_eigenvalues);
	cout <<"\nEigenvectors (represented by column): "<<endl;	
	vector<vector<complex<double>>> matrix_eigenvectors_transpose = TransposeComplexMatrix(matrix_eigenvectors);
	printComplexMatrix(matrix_eigenvectors_transpose);
	
	// Spectral decomposition computation starts here
	for (int i = 0; i < n; ++i)
	{
		vector<complex<double>> vi = getcolumnComplexMatrix(matrix_eigenvectors_transpose,i);
	
		vector<vector<complex<double>>> mat_u;	
		vector<vector<complex<double>>> mat_u_transpose;
		mat_u.assign(n, std::vector<complex<double>>(1, 0.0));
		mat_u_transpose.assign(1, std::vector<complex<double>>(n, 0.0));

		for (int i = 0; i < n ; ++i) 
		{		
			mat_u[i][0] = vi[i];
			mat_u_transpose[0][i] = vi[i];
		}
		vector<vector<complex<double>>> SDi = multiplyComplexMatrices(mat_u,mat_u_transpose);

		cout <<"\nSpectral decomposition for λ_{" << i+1 <<"} = " << vector_eigenvalues[i] << " * " <<endl;	
		printComplexMatrix(SDi);
		
	}

	return matrix_eigenvectors_transpose;
}


vector<double> PowerMethod(vector<vector<double>>& A, int iterations) 
{
	int n = A.size();
	double mu = 2; 
	double sigma = 0.5;
	// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
	vector<double> v = vrandn_normal(mu, sigma, n);

	cout << "v initial vector = " << endl;
	printVector(v);
	// Repeatedly multiply the matrix A by the current vector v to get v_{k+1} = A v_{k}
	// Then normalize v_{k+1} at each step to prevent the values from growing too large.
	// The normalized vector converges to the dominant eigenvector, and the corresponding eigenvalue can be calculated from Av = λv
	for (int i = 0; i < iterations; ++i) 
	{
		vector<double> w = multiplymatrixvector(A,v);
		double ssw = 0; // sum of squares of w.
		for (int j = 0; j < n ; ++j)
		{
			ssw += w[j]*w[j];
		}
		double norm = std::sqrt(ssw);
		// Normalize w (w = w / norm)
		vector<double> w_normalize = scalarmultiplication(w,divisiond(1,norm));
		
		for (int j = 0; j < n ; ++j)
		{
			v[j] = w_normalize[j];
		}
	}
	// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
	double lambda = divisiond(quadraticmultiplication(A,v),dot(v,v));
	cout << "\nEigenvalue = " << lambda << endl;
	cout << "\nEigenvector corresponds to " << lambda << " :" << endl;
	printVector(v);
	
	return v; // returns the eigenvector
}

vector<double> EigenvaluesEigenvectorsApproximation(vector<vector<double>>& A, int iterations) 
// The power method is used first to compute the largest eigenvalue
// Method of deflation to approximate all the eigenvalues
// Inverse iteration method to compute the corresponding eigenvector for each eigenvalue

{
	int n = A.size();
	int m = n;
	vector<double> vector_eigenvalues;	
	vector<vector<double>> matrix_eigenvectors;	
	vector<double> v;
	vector<vector<double>> A_original(n, vector<double>(n,0.0));	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	cout <<"\nA: "<<endl;
	printMatrix(A);
	for (int k = 0; k < n; ++k)
	{
		double mu = 2; 
		double sigma = 0.5;
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		v = vrandn_normal(mu, sigma, m);

		//cout <<"\nA: "<<endl;
		//printMatrix(A);
		
		//cout << "\nv initial vector = " << endl;
		//printVector(v);
		// Repeatedly multiply the matrix A by the current vector v to get v_{k+1} = A v_{k}
		// Then normalize v_{k+1} at each step to prevent the values from growing too large.
		// The normalized vector converges to the dominant eigenvector, and the corresponding eigenvalue can be calculated from Av = λv
		for (int i = 0; i < iterations; ++i) 
		{
			vector<double> w = multiplymatrixvector(A,v);
			double ssw = 0; // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w[j]*w[j];
			}
			double norm = std::sqrt(ssw);
			// Normalize w (w = w / norm)
			vector<double> w_normalize = scalarmultiplication(w,divisiond(1,norm));
			
			for (int j = 0; j < n ; ++j)
			{
				v[j] = w_normalize[j];
			}
		}
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		double lambda = divisiond(quadraticmultiplication(A,v),dot(v,v));
		//cout << "\nEigenvalue = " << lambda << endl;

		vector_eigenvalues.push_back(lambda);

		//cout << "\nEigenvector corresponds to " << lambda << " :" << endl;
		//printVector(v);
		vector<vector<double>> P = PermutationMatrixMax(v);
		//cout << "\nP = " << endl;
		//printMatrix(P);

		vector<double> wi = multiplymatrixvector(P,v);
		//cout << "\nw_{1} = " << endl;
		//printVector(wi);

		vector<vector<double>> R_inv(m, vector<double>(m));
		for (int i = 0; i < m ; ++i)
		{
			R_inv[i][0] = wi[i];
		}
		for (int i = 1; i < m ; ++i)
		{
			for (int j = 1; j < m ; ++j)
			{
				if (i != j)
				{
					R_inv[i][j] = 0;
				}
				else if (i == j)
				{
					R_inv[i][j] = 1;
				}
			}
		}
		vector<vector<double>> R= inverse(R_inv);
		vector<vector<double>> RP = multiply(R,P);
		vector<vector<double>> PiRi = multiply(P,R_inv);
		vector<vector<double>> APiRi = multiply(A,PiRi);
		vector<vector<double>> RPAPiRi = multiply(RP,APiRi);
		
		/*
		cout << "\nR^{-1} = " << endl;
		printMatrix(R_inv);

		cout << "\nR = " << endl;
		printMatrix(R);

		cout << "\nRP = " << endl;
		printMatrix(RP);

		cout << "\nP^{-1}R^{-1} = " << endl;
		printMatrix(PiRi);

		cout << "\nAP^{-1}R^{-1} = " << endl;
		printMatrix(APiRi);

		cout << "\nB = RPAP^{-1}R^{-1}  = " << endl;
		printMatrix(RPAPiRi);
		*/

		A = deleteColumn(RPAPiRi,0); 
		A = deleteRow(RPAPiRi,0); 
		m = m-1;
	}
	
	// to determine the whole eigenvectors from known eigenvalue
	// we use inverse iteration with shift from power method
	for (int i = 0; i < n ;++i)
	{
		double eigen = vector_eigenvalues[i];
		vector<double> vector_eigen;
		cout << "\nEigenvalue(" << i+1 << ") = " << eigen << endl;
		
		vector<vector<double>> I = createIdentityMatrix(n);
		// why we use subtraction to 0.1? because we are using inverse iteration method and the inverse of (A - λI) does not exist 
		vector<vector<double>> lambda_I = scalarmultiplication(I,eigen-0.1); 
		//cout << "\nlambda I = " << endl;
		//printMatrix(lambda_I);
		vector<vector<double>> A_lambda_I = matrixsubtraction(A_original,lambda_I);
		vector<vector<double>> A_lambda_I_inverse = inverse(A_lambda_I);
		//cout << "\nA - lambda I" << endl;
		//printMatrix(A_lambda_I);
		double mu = 2; 
		double sigma = 0.5;
		vector_eigen = vrandn_normal(mu, sigma, n);
		double ssv = 0; // sum of squares of v.
		for (int j = 0; j < n ; ++j)
		{
			ssv += vector_eigen[j]*vector_eigen[j];
		}
		double norm = std::sqrt(ssv);
		// Normalize v (v = v / norm)
		vector_eigen = scalarmultiplication(vector_eigen,divisiond(1,norm));
			
		//cout << "\nv initial vector = " << endl;
		//printVector(vector_eigen);
		
		for (int k = 0; k < iterations; ++k) 
		{
			vector<double> w_eigen = multiplymatrixvector(A_lambda_I_inverse,vector_eigen);
			double ssw = 0; // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w_eigen[j]*w_eigen[j];
			}
			double norm_w = std::sqrt(ssw);
			// Normalize w (w = w / norm)
			vector<double> w_eigen_normalize = scalarmultiplication(w_eigen,divisiond(1,norm_w));
			
			for (int j = 0; j < n ; ++j)
			{
				vector_eigen[j] = w_eigen_normalize[j];
			}
			
		}
		cout << "\nCorresponding Eigenvector: " << endl;
		printVector(vector_eigen);
		matrix_eigenvectors.push_back(vector_eigen);
	}

	matrix_eigenvectors = transpose(matrix_eigenvectors);
	cout << "\nEigenvectors (in column vectors): " << endl;
	printMatrix(matrix_eigenvectors);

	return vector_eigenvalues; // returns the eigenvector
}

vector<vector<complex<double>>> QRDecompositionHouseholderReflections(vector<vector<complex<double>>>& A) 
{ // Done on January 9th, 2026

	int n = A.size(); // Rows	
	int m = A[0].size();
	
	if(m > n)
	{
		cerr << "Error: The column dimension is bigger than the row dimension" << endl;	
	}

	cout << "\nMatrix A : " << endl;
	printComplexMatrix(A);
	
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(m, 0.0));
	vector<vector<complex<double>>> R_final;	
	R_final.assign(n, std::vector<complex<double>>(m, 0.0));
	vector<vector<complex<double>>> Q = createIdentityComplexMatrix(n);
	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < m ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	int iter;
	if (m == n) // iterations is the number of column-1 of matrix A if m=n
	{
		iter = m-1;		
	}
	if (m < n) // iterations depend on the number of column of matrix A if m<n
	{
		iter = m;	
	}	
	for (int i = 0; i < iter; ++i) 
	{		

		cout << "\nIteration : " << i+1 << endl;
	
		// compute the next vector
		vector<complex<double>> x;  
		vector<complex<double>> w(n, 0.0);  
		vector<complex<double>> v;  
			
		for (int j = 0; j < n ; ++j) 
		{		
			x.push_back(A[j][i]);
		}
		if (i > 0)
		{
			for (int j = 0; j < i ; ++j) 
			{		
				x[j] = 0;
			}
			
		}
		double real_x= real(x[i]); // check again should we use x[i] or x[0], because x[0]  will always be 0 after iteration 1.
		double norm_x = complexnorm(x); 
		// w = (sgn) ||x_{i}|| e_{i}
		if(real_x < 0 )
		{
			w[i] = norm_x;
		}
		else if(real_x >= 0 )
		{
			w[i] = -norm_x;
		}
		// v = x_{i} - ||x_{i}|| e_{i}
		v = subtractComplexVectors(x,w);
		
		cout << "\nx: " << endl;
		printComplexVector(x);
		cout << "\nw: " << endl;
		printComplexVector(w);
		cout << "\nv: " << endl;
		printComplexVector(v);
		vector<vector<complex<double>>> mat_v;
		vector<vector<complex<double>>> mat_v_transpose;
		mat_v.assign(n, std::vector<complex<double>>(1, 0.0));
		mat_v_transpose.assign(1, std::vector<complex<double>>(n, 0.0));

		for (int i = 0; i < n ; ++i) 
		{		
			mat_v[i][0] = v[i];
			mat_v_transpose[0][i] = v[i];
		}

		vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
		vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		vector<vector<complex<double>>> H_hat; 
		vector<vector<complex<double>>> R; 
		
		double norm_v = divisiond(1,real(vTv[0][0]));

		vector<vector<complex<double>>> P;
		P = scalarmultiplicationComplexMatrix(vvT,norm_v);
		
		P = scalarmultiplicationComplexMatrix(P,2);
		H_hat = subtractComplexMatrices(I,P);
		
		R = multiplyComplexMatrices(H_hat,A);

		// Rounding down those with small decimal number to 0
		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < m ; ++j)
			{
				double real_value = real(R[i][j]) ;
				if(abs(real_value) < 1e-12)
				{
					R[i][j] = 0;
				}
			}
		}

		/*
		cout << "\nvTv: " << endl;
		printComplexMatrix(vTv);
		cout << "\nvvT: " << endl;
		printComplexMatrix(vvT);
		cout << "\nP: " << endl;
		printComplexMatrix(P);
		cout << "\nI - 2P: " << endl;
		printComplexMatrix(H_hat);*/
		cout << "\nR: " << endl;
		printComplexMatrix(R);

		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < m ; ++j)
			{
				A[i][j] = R[i][j];
				R_final[i][j] = R[i][j];
			}
		}
		
		Q = multiplyComplexMatrices(Q,H_hat);
		// Rounding down those with small decimal number to 0
		for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < m ; ++j)
			{
				double real_value = real(Q[i][j]) ;
				if(abs(real_value) < 1e-12)
				{
					Q[i][j] = 0;
				}
			}
		}
		cout << "\nHouseholder Matrix / Q: " << endl;
		printComplexMatrix(Q);

	}
	
	
	cout << "\n****************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\n****************************************************************************" << endl;
	
	cout << "\nQ: " << endl;
	printComplexMatrix(Q);
	cout << "\nR: " << endl;
	printComplexMatrix(R_final);

	cout << "\nQ * R: " << endl;
	vector<vector<complex<double>>> AQR = multiplyComplexMatrices(Q,R_final);
	printComplexMatrix(AQR);
	return Q;
	
}

vector<complex<double>> QRAlgorithmwithShifts(vector<vector<complex<double>>>& A, int iterations) 
{

	int n = A.size(); // Rows	
	int C = A[0].size();
	
	vector<complex<double>> vector_eigenvalues;	
	vector<vector<complex<double>>> matrix_eigenvectors;	
	vector<complex<double>> v;
	vector<vector<complex<double>>> A_original(n, vector<complex<double>>(n));
	vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
	vector<vector<complex<double>>> Q_final(n, vector<complex<double>>(n));
	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	
	// Repeatedly decompose complex matrix A to QR, A^{k-1} - μI = QR
	// Then multiply the reverse to obtain A^{k} = RQ + μI
	
	vector<vector<complex<double>>> A_loop(n, vector<complex<double>>(n));
	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_loop[i][j] = A_original[i][j];
		}
	}
	for (int i = 0; i < iterations; ++i) 
	{
		double a_complex = real(A_loop[n-1][n-1]);
		double b_complex = imag(A_loop[n-1][n-1]);
		complex<double> mu(a_complex,b_complex);
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		vector<vector<complex<double>>> mu_I = complexnumbermultiplicationComplexMatrix(I,mu); 
		vector<vector<complex<double>>> A_mu_I = subtractComplexMatrices(A_loop,mu_I);

		cmat Q(n, vector<complex<double>>(C));
		cmat R(n, vector<complex<double>>(C));

		QRDecompositionComplexHouseholder(A_mu_I, Q, R);

		vector<vector<complex<double>>> RQ = multiplyComplexMatrices(R,Q);

		cout <<"\nIteration: "<< i + 1 << endl;
		cout << "\nμ = " << mu << endl;
		cout <<"\nA: "<<endl;
		printComplexMatrix(A_loop);
		cout << "\nμI = " << endl;
		printComplexMatrix(mu_I);
		cout << "\nA - μI = " << endl;
		printComplexMatrix(A_mu_I);
		cout <<"\nQ: "<<endl;
		printComplexMatrix(Q);
		cout <<"\nR: "<<endl;
		printComplexMatrix(R);
		cout << "\nRQ = " << endl;
		printComplexMatrix(RQ);
			
			for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < n ; ++j)
			{
				A_loop[i][j] = RQ[i][j] + mu_I[i][j];
				Q_final[i][j] = Q[i][j];
			}
		}
		cout << "\nA_{new} = " << endl;
		printComplexMatrix(A_loop);
	}
	double sum_lowerdiagonal;
	for (int i = 1; i < n ; ++i) 
	{		
		for (int j = 0; j < i ; ++j) 
		{
			sum_lowerdiagonal += real(A_loop[i][j]);	
		}	
	}
			
	if(abs(sum_lowerdiagonal) < 1e-12)
	{
		for (int i = 0; i < n ; ++i) 
		{
			vector_eigenvalues.push_back(A_loop[i][i]);
		}
		
	}
	
	cout << "\n****************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\n****************************************************************************" << endl;
	
	// to determine the whole eigenvectors from known eigenvalue
	// we use inverse iteration with shift from power method
	for (int i = 0; i < n ;++i)
	{
		complex<double> eigen = vector_eigenvalues[i];
		vector<complex<double>> vector_eigen;
		cout << "\nEigenvalue(" << i+1 << ") = " << eigen << endl;
		
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		// why we use subtraction to 0.1? because we are using inverse iteration method and the inverse of (A - λI) does not exist 
		vector<vector<complex<double>>> lambda_I = scalarmultiplicationComplexMatrix(I,real(eigen)-0.1); 
		//cout << "\nlambda I = " << endl;
		//printMatrix(lambda_I);
		vector<vector<complex<double>>> A_lambda_I = subtractComplexMatrices(A_original,lambda_I);
		vector<vector<complex<double>>> A_lambda_I_inverse = ComplexMatrixInverse(A_lambda_I);
		//cout << "\nA - lambda I" << endl;
		//printMatrix(A_lambda_I);
		double mu = 2; 
		double sigma = 0.5;
		vector_eigen = complexvecrand_normal_zeroimaginary(mu, sigma, n);
		complex<double> ssv(0,0); // sum of squares of v.
		for (int j = 0; j < n ; ++j)
		{
			ssv += vector_eigen[j]*vector_eigen[j];
		}
		double norm = std::sqrt(real(ssv));
		// Normalize v (v = v / norm)
		vector_eigen = scalarmultiplicationComplexVector(vector_eigen,divisiond(1,norm));
			
		//cout << "\nv initial vector = " << endl;
		//printVector(vector_eigen);
		
		for (int k = 0; k < iterations; ++k) 
		{
			vector<complex<double>> w_eigen = multiplycomplexmatrixvector(A_lambda_I_inverse,vector_eigen);
			complex<double> ssw(0,0); // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w_eigen[j]*w_eigen[j];
			}
			double norm_w = std::sqrt(real(ssw));
			// Normalize w (w = w / norm)
			vector<complex<double>> w_eigen_normalize = scalarmultiplicationComplexVector(w_eigen,divisiond(1,norm_w));
			
			for (int j = 0; j < n ; ++j)
			{
				vector_eigen[j] = w_eigen_normalize[j];
			}
			
		}
		cout << "\nCorresponding Eigenvector: " << endl;
		printComplexVector(vector_eigen);
		matrix_eigenvectors.push_back(vector_eigen);
	}

	matrix_eigenvectors = TransposeComplexMatrix(matrix_eigenvectors);
	cout << "\nEigenvectors (in column vectors): " << endl;
	printComplexMatrix(matrix_eigenvectors);

	cout << "\nEigenvalues: " << endl;
	printComplexVector(vector_eigenvalues);

	return vector_eigenvalues; // returns the eigenvalues
}

// Homework
vector<complex<double>> QRAlgorithmDoubleShifts(vector<vector<complex<double>>>& A, int iterations) 
{

	int n = A.size(); // Rows	
	int C = A[0].size();
	
	vector<complex<double>> vector_eigenvalues;	
	vector<vector<complex<double>>> matrix_eigenvectors;	
	vector<complex<double>> v;
	vector<vector<complex<double>>> A_original(n, vector<complex<double>>(n));
	vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
	vector<vector<complex<double>>> Q_final(n, vector<complex<double>>(n));
	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	
	// Repeatedly decompose complex matrix A to QR, A^{k-1} - μI = QR
	// Then multiply the reverse to obtain A^{k} = RQ + μI
	
	vector<vector<complex<double>>> A_loop(n, vector<complex<double>>(n));
	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_loop[i][j] = A_original[i][j];
		}
	}
	for (int i = 0; i < iterations; ++i) 
	{
		double a_complex = real(A_loop[n-1][n-1]);
		double b_complex = imag(A_loop[n-1][n-1]);
		complex<double> mu(a_complex,b_complex);
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		vector<vector<complex<double>>> mu_I = complexnumbermultiplicationComplexMatrix(I,mu); 
		vector<vector<complex<double>>> A_mu_I = subtractComplexMatrices(A_loop,mu_I);

		cmat Q(n, vector<complex<double>>(C));
		cmat R(n, vector<complex<double>>(C));

		QRDecompositionComplexHouseholder(A_mu_I, Q, R);

		vector<vector<complex<double>>> RQ = multiplyComplexMatrices(R,Q);

		cout <<"\nIteration: "<< i + 1 << endl;
		cout << "\nμ = " << mu << endl;
		cout <<"\nA: "<<endl;
		printComplexMatrix(A_loop);
		cout << "\nμI = " << endl;
		printComplexMatrix(mu_I);
		cout << "\nA - μI = " << endl;
		printComplexMatrix(A_mu_I);
		cout <<"\nQ: "<<endl;
		printComplexMatrix(Q);
		cout <<"\nR: "<<endl;
		printComplexMatrix(R);
		cout << "\nRQ = " << endl;
		printComplexMatrix(RQ);
			
			for (int i = 0; i < n ; ++i) 
		{		
			for (int j = 0; j < n ; ++j)
			{
				A_loop[i][j] = RQ[i][j] + mu_I[i][j];
				Q_final[i][j] = Q[i][j];
			}
		}
		cout << "\nA_{new} = " << endl;
		printComplexMatrix(A_loop);
	}
	double sum_lowerdiagonal;
	for (int i = 1; i < n ; ++i) 
	{		
		for (int j = 0; j < i ; ++j) 
		{
			sum_lowerdiagonal += real(A_loop[i][j]);	
		}	
	}
			
	if(abs(sum_lowerdiagonal) < 1e-12)
	{
		for (int i = 0; i < n ; ++i) 
		{
			vector_eigenvalues.push_back(A_loop[i][i]);
		}
		
	}
	
	cout << "\n****************************************************************************" << endl;
	cout << "\nEnd of iteration" << endl;
	cout << "\n****************************************************************************" << endl;
	
	// to determine the whole eigenvectors from known eigenvalue
	// we use inverse iteration with shift from power method
	for (int i = 0; i < n ;++i)
	{
		complex<double> eigen = vector_eigenvalues[i];
		vector<complex<double>> vector_eigen;
		cout << "\nEigenvalue(" << i+1 << ") = " << eigen << endl;
		
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		// why we use subtraction to 0.1? because we are using inverse iteration method and the inverse of (A - λI) does not exist 
		vector<vector<complex<double>>> lambda_I = scalarmultiplicationComplexMatrix(I,real(eigen)-0.1); 
		//cout << "\nlambda I = " << endl;
		//printMatrix(lambda_I);
		vector<vector<complex<double>>> A_lambda_I = subtractComplexMatrices(A_original,lambda_I);
		vector<vector<complex<double>>> A_lambda_I_inverse = ComplexMatrixInverse(A_lambda_I);
		//cout << "\nA - lambda I" << endl;
		//printMatrix(A_lambda_I);
		double mu = 2; 
		double sigma = 0.5;
		vector_eigen = complexvecrand_normal_zeroimaginary(mu, sigma, n);
		complex<double> ssv(0,0); // sum of squares of v.
		for (int j = 0; j < n ; ++j)
		{
			ssv += vector_eigen[j]*vector_eigen[j];
		}
		double norm = std::sqrt(real(ssv));
		// Normalize v (v = v / norm)
		vector_eigen = scalarmultiplicationComplexVector(vector_eigen,divisiond(1,norm));
			
		//cout << "\nv initial vector = " << endl;
		//printVector(vector_eigen);
		
		for (int k = 0; k < iterations; ++k) 
		{
			vector<complex<double>> w_eigen = multiplycomplexmatrixvector(A_lambda_I_inverse,vector_eigen);
			complex<double> ssw(0,0); // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w_eigen[j]*w_eigen[j];
			}
			double norm_w = std::sqrt(real(ssw));
			// Normalize w (w = w / norm)
			vector<complex<double>> w_eigen_normalize = scalarmultiplicationComplexVector(w_eigen,divisiond(1,norm_w));
			
			for (int j = 0; j < n ; ++j)
			{
				vector_eigen[j] = w_eigen_normalize[j];
			}
			
		}
		cout << "\nCorresponding Eigenvector: " << endl;
		printComplexVector(vector_eigen);
		matrix_eigenvectors.push_back(vector_eigen);
	}

	matrix_eigenvectors = TransposeComplexMatrix(matrix_eigenvectors);
	cout << "\nEigenvectors (in column vectors): " << endl;
	printComplexMatrix(matrix_eigenvectors);

	cout << "\nEigenvalues: " << endl;
	printComplexVector(vector_eigenvalues);

	return vector_eigenvalues; // returns the eigenvalues
}

vector<complex<double>> FrancisDoubleStepQRReal(vector<vector<complex<double>>>& A, int iterations) 
// Done on January 18th, 2026 took almost 1 week to comprehend and finish this
{

	int n = A.size(); // Rows	
	int C = A[0].size();
	if (C != n) 
	{
		cerr << "Error: Matrix is not square." << endl;
		return {};
	}
	vector<complex<double>> vector_eigenvalues;	
	vector<vector<complex<double>>> matrix_eigenvectors;	
	
	vector<vector<double>> I = createIdentityMatrix(n);

	cout << "\nA = " << endl;
	printComplexMatrix(A);
	
	vector<vector<complex<double>>> H = HessenbergDecompositionResultOnly(A);
	vector<vector<double>> H_loop = ComplextoRealMatrix(H);
	
	cout << "\nH_{0} = " << endl;
	printComplexMatrix(H);
	cout << "\nH_{0} = " << endl;
	printMatrix(H_loop);
	int iter = 1;
	int q_index, p_index;
	int p = n;
	double s,t,x,y,z;		
	for (int pindex = 0; pindex < iterations; pindex++ )
	{
		cout << "\nIteration: " << iter << endl;
		cout << "\nH = " << endl;
		printMatrix(H_loop);
		p_index = p-1;
		q_index = p-2;
		
		s = H_loop[q_index][q_index] + H_loop[p_index][p_index];
		t = (H_loop[q_index][q_index] * H_loop[p_index][p_index]) - (H_loop[q_index][p_index] * H_loop[p_index][q_index]);
		cout << "\ns = " << s << endl;	
		cout << "t = " << t << endl;	
		x = (H_loop[0][0]*H_loop[0][0]) + (H_loop[0][1]*H_loop[1][0]) - s*H_loop[0][0] + t;
		y = H_loop[1][0]*(H_loop[0][0] + H_loop[1][1] - s);
		z = H_loop[1][0]*H_loop[2][1];
		cout << "x = " << x << endl;	
		cout << "y = " << y << endl;	
		cout << "z = " << z << endl;	

		vector<double> e1(3,0.0);
		e1[0] = 1.0;
		vector<double> vec_x(3,0.0) ;
		vec_x[0] = x;
		vec_x[1] = y;
		vec_x[2] = z;

		double norm_x = norm(vec_x);

		if(vec_x[0] > 0 )
		{
			norm_x = -1*norm_x;		
		}		
		else if(vec_x[0] <= 0 )
		{
			norm_x = norm_x;		
		}		

		vector<double> vec_w = scalarmultiplication(e1,norm_x);
		vector<double> vec_u = add(vec_x,vec_w);

		cout << "\nx = " << endl;
		printVector(vec_x);
		cout << "\nw = " << endl;
		printVector(vec_w);
		cout << "\nu = " << endl;
		printVector(vec_u);
		

		vector<vector<double>> mat_u; 
		vector<vector<double>> mat_u_transpose; 
		mat_u.assign(3, std::vector<double>(1, 0.0));
		mat_u_transpose.assign(1, std::vector<double>(3, 0.0));
		
		for (int i = 0; i < 3 ; ++i) 
		{		
			mat_u[i][0] = vec_u[i];
			mat_u_transpose[0][i] =vec_u[i];
		}
		cout << "\nu = " << endl;
		printMatrix(mat_u);
		cout << "\nuT = " << endl;
		printMatrix(mat_u_transpose);

		vector<vector<double>> uTu = multiply(mat_u_transpose,mat_u);
		vector<vector<double>> uuT = multiply(mat_u,mat_u_transpose);
		vector<vector<double>> P0; 
		P0.assign(3,std::vector<double>(3, 0.0));
		for (int i = 0; i < 3 ; ++i) 
		{		
			for (int j = 0; j < 3 ; ++j) 
			{		
				P0[i][j] = uuT[i][j];
			}
		}
		
		cout << "\nuTu = " << endl;
		printMatrix(uTu);
		cout << "\nuuT = " << endl;
		printMatrix(uuT);

		
		double norm_u = divisiond(1,uTu[0][0]); // ||u||^{2}
		cout << "\n||u|| = " << norm_u << endl;
		// compute 2 u u^{t} / ||u||
		P0 = scalarmultiplication(P0,norm_u);
		P0 = scalarmultiplication(P0,2); 
		//cout << "\nP_{" << iter <<"}:" << endl;
		//printMatrix(P0);

		vector<vector<double>> P0_final; 
		P0_final.assign(n,std::vector<double>(n, 0.0));
	
		for (int i = 0; i < 3 ; ++i) 
		{		
			for (int j = 0; j < 3 ; ++j) 
			{		
				P0_final[i][j] = P0[i][j];
			}
		}
		// compute I - (2 u u^{t} / ||u|| )
		P0_final = subtract(I,P0_final);
		
		cout << "\nP_{" << iter  - 1<<"}:" << endl;
		printMatrix(P0_final);

		vector<vector<double>> B = multiply(H_loop,P0_final); 
		B = multiply(P0_final,B); 

		cout << "\nB_{" << iter  - 1<<"}:" << endl;
		printMatrix(B);
		vector<vector<double>> H_new = multiply(I,B); 
		// "Bulge chasing" sequence coding
	
		int index = 1;
		int index2 = 0;
		for (int k = 0; k <= p-2; ++k)
		{
			vector<vector<double>> P1; 
			vector<vector<double>> P_last; 
			vector<double> vec_x1(3);
			if(k < p-3)
			{			
				vec_x1[0] = H_new[index][k];
				vec_x1[1] = H_new[index+1][k];
				vec_x1[2] = H_new[index+2][k];
				
				double norm_x1 = norm(vec_x1);
			
				if(real(vec_x1[0]) > 0 )
				{
					norm_x1 = -1*norm_x1;		
				}		
				else if(real(vec_x1[0]) <= 0 )
				{
					norm_x1 = norm_x1;		
				}		

				cout << "\nx_{" << k+1 <<"}:" << endl;
				printVector(vec_x1);

				vector<double> e1_3d(3,0.0);
				e1_3d[0] = 1.0;

				vector<double> vec_w1 = scalarmultiplication(e1_3d,norm_x1);
				vector<double> vec_v1 = add(vec_x1,vec_w1);

				vector<vector<double>> mat_v; 
				vector<vector<double>> mat_v_transpose; 
				mat_v.assign(3, std::vector<double>(1, 0.0));
				mat_v_transpose.assign(1, std::vector<double>(3, 0.0));

				for (int i = 0; i < 3 ; ++i) 
				{		
					mat_v[i][0] = vec_v1[i];
					mat_v_transpose[0][i] =vec_v1[i];
				}

				vector<vector<double>> vTv = multiply(mat_v_transpose,mat_v);
				vector<vector<double>> vvT = multiply(mat_v,mat_v_transpose);
				
				double norm_v1 = divisiond(1,vTv[0][0]); // ||u||^{2}

				// compute 2 u u^{t} / ||u||
				P1 = scalarmultiplication(vvT,norm_v1);				
				P1 = scalarmultiplication(P1,2);

				vector<vector<double>> P1_final;
				P1_final.assign(n, std::vector<double>(n, 0.0));
				int i1 = 0;
				for (int i = 1 + index2; i <= 3 + index2; ++i)
				{
					int j1=0;
					for (int j = 1 + index2; j <= 3 + index2; ++j)
					{
						P1_final[i][j] = P1[i1][j1];
						j1 = j1+1;
					}
					i1= i1 + 1;
				}
				vector<vector<double>> P1_f = subtract(I,P1_final);
				cout << "\nP_{" << k+1 <<"} reflector:" << endl;
				printMatrix(P1_f);

				H_new = multiply(H_new,P1_f);
				H_new = multiply(P1_f,H_new);

				cout << "\nH^{" << k+1 <<"}:" << endl;
				printMatrix(H_new);
				
			}
			else if (k == p-3)
			{
				vector<double> vec_x2(2);
				vec_x2[0] = H_new[index][k];
				vec_x2[1] = H_new[index+1][k];
				
				double norm_x2 = norm(vec_x2);
			
				if(vec_x2[0] > 0 )
				{
					norm_x2 = -1*norm_x2;		
				}		
				else if(vec_x2[0] <= 0 )
				{
					norm_x2 = norm_x2;		
				}		

				cout << "\nx_{" << k+1 <<"}:" << endl;
				printVector(vec_x2);

				vector<double> e1_2d(2,0.0);
				e1_2d[0] = 1.0;

				vector<double> vec_w2 = scalarmultiplication(e1_2d,norm_x2);
				vector<double> vec_v2 = add(vec_x2,vec_w2);

				vector<vector<double>> mat_v; 
				vector<vector<double>> mat_v_transpose; 
				mat_v.assign(2, std::vector<double>(1, 0.0));
				mat_v_transpose.assign(1, std::vector<double>(2, 0.0));

				for (int i = 0; i < 2 ; ++i) 
				{		
					mat_v[i][0] = vec_v2[i];
					mat_v_transpose[0][i] =vec_v2[i];
				}

				vector<vector<double>> vTv = multiply(mat_v_transpose,mat_v);
				vector<vector<double>> vvT = multiply(mat_v,mat_v_transpose);
				
				double norm_v2 = divisiond(1,vTv[0][0]); // ||u||^{2}

				P_last = scalarmultiplication(vvT,norm_v2);
				
				P_last = scalarmultiplication(P_last,2);
				//cout << "\nP_{" << k+1 <<"} last:" << endl;
				//printMatrix(P_last);
				vector<vector<double>> P_last_final;
				P_last_final.assign(n, std::vector<double>(n, 0.0));
				
				vector<vector<double>> P_last_f = subtract(I,P_last_final);
			
				P_last_f[n-2][n-2] = 1 - P_last[0][0] ;
				P_last_f[n-2][n-1] =0 - P_last[0][1] ;
				P_last_f[n-1][n-2] = 0 - P_last[1][0] ;
				P_last_f[n-1][n-1] = 1 - P_last[1][1] ;

				cout << "\nP_{" << k+1 <<"} reflector:" << endl;
				printMatrix(P_last_f);
				H_new = multiply(H_new,P_last_f);
				H_new = multiply(P_last_f,H_new);

				cout << "\nH^{" << k+1 <<"}:" << endl;
				printMatrix(H_new);
				
			}
			else if (k == p-2)
			{
				vector<double> vec_x2(2);
				vec_x2[0] = H_new[index2][k];
				vec_x2[1] = H_new[index2+1][k];
				
				double norm_x2 = norm(vec_x2);
			
				if(vec_x2[0] > 0 )
				{
					norm_x2 = -1*norm_x2;		
				}		
				else if(vec_x2[0] <= 0 )
				{
					norm_x2 = norm_x2;		
				}		

				cout << "\nx_{" << k+1 <<"}:" << endl;
				printVector(vec_x2);

				vector<double> e1_2d(2,0.0);
				e1_2d[0] = 1.0;

				vector<double> vec_w2 = scalarmultiplication(e1_2d,norm_x2);
				vector<double> vec_v2 = add(vec_x2,vec_w2);

				vector<vector<double>> mat_v; 
				vector<vector<double>> mat_v_transpose; 
				mat_v.assign(2, std::vector<double>(1, 0.0));
				mat_v_transpose.assign(1, std::vector<double>(2, 0.0));

				for (int i = 0; i < 2 ; ++i) 
				{		
					mat_v[i][0] = vec_v2[i];
					mat_v_transpose[0][i] =vec_v2[i];
				}

				vector<vector<double>> vTv = multiply(mat_v_transpose,mat_v);
				vector<vector<double>> vvT = multiply(mat_v,mat_v_transpose);
				
				double norm_v2 = divisiond(1,vTv[0][0]); // ||u||^{2}

				P_last = scalarmultiplication(vvT,norm_v2);
				
				P_last = scalarmultiplication(P_last,2);
				cout << "\nP_{" << k+1 <<"} last:" << endl;
				printMatrix(P_last);
				vector<vector<double>> P_last_final;
				P_last_final.assign(n, std::vector<double>(n, 0.0));
				
				vector<vector<double>> P_last_f = subtract(I,P_last_final);
			
				P_last_f[n-1][n-1] = 1 - P_last[0][0] ;

				cout << "\nP_{" << k+1 <<"} reflector:" << endl;
				printMatrix(P_last_f);
				H_new = multiply(H_new,P_last_f);
				H_new = multiply(P_last_f,H_new);

				cout << "\nH^{" << k+1 <<"}:" << endl;
				printMatrix(H_new);
				
			}
			index = index+1;
			index2=index2+1;
		}
		H_loop = multiply(I,H_new);
		iter = iter+1;
		
	}
	cout << "\n****************************************************************************" << endl;
	cout << "\nEnd of " << iter - 1 << " iterations" << endl;
	cout << "\n****************************************************************************" << endl;
	
	cout << "\nH updated :" << endl;
	printMatrix(H_loop);
	
	for (int i = 0; i < n ; ++i) 
	{		
		complex<double> a(1,0);
		complex<double> b = - H_loop[i+1][i+1] - H_loop[i][i];
		complex<double> c = (H_loop[i+1][i+1]*H_loop[i][i]) - (H_loop[i+1][i]*H_loop[i][i+1]);
			
		complex<double> σ1 = (-b + sqrt( (b*b) - (complex<double>(4,0)*a*c)) )/(complex<double>(2,0)*a);
		complex<double> σ2 = (-b - sqrt( (b*b) - (complex<double>(4,0)*a*c)) )/(complex<double>(2,0)*a);
		//cout << "\nσ1 :" <<   σ1 << endl;
		//cout << "\nσ2 :" <<   σ2 << endl;
		vector_eigenvalues.push_back(σ1);
		vector_eigenvalues.push_back(σ2);
		i = i+1;		
	}
	//matrix_eigenvectors = TransposeComplexMatrix(matrix_eigenvectors);
	//cout << "\nEigenvectors (in column vectors): " << endl;
	//printComplexMatrix(matrix_eigenvectors);

	cout << "\nEigenvalues: " << endl;
	printComplexVector(vector_eigenvalues);

	return vector_eigenvalues; // returns the eigenvalues
}


vector<complex<double>> FrancisDoubleStepQR(vector<vector<complex<double>>>& A, int iterations) 
{

	int n = A.size(); // Rows	
	int C = A[0].size();
	if (C != n) 
	{
		cerr << "Error: Matrix is not square." << endl;
		return {};
	}
	vector<complex<double>> vector_eigenvalues;	
	vector<vector<complex<double>>> matrix_eigenvectors;	
	vector<complex<double>> v;
	vector<vector<complex<double>>> A_original(n, vector<complex<double>>(n));
	vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
	vector<vector<complex<double>>> Q_final(n, vector<complex<double>>(n));
	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	cout << "\nA = " << endl;
	printComplexMatrix(A);
	// Repeatedly decompose complex matrix A to QR, A^{k-1} - μI = QR
	// Then multiply the reverse to obtain A^{k} = RQ + μI
	
	vector<vector<complex<double>>> H = HessenbergDecompositionResultOnly(A);
	vector<vector<complex<double>>> H_loop;
	H_loop.assign(n, std::vector<complex<double>>(n, 0.0));

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			H_loop[i][j] = H[i][j];
		}
	}

	cout << "\nH_{0} = " << endl;
	printComplexMatrix(H);
	int iter = 1;
	int p = n;
	int q_index, p_index;
	for (int i = 0; i < iterations; i++)
	{
		cout << "\nIteration: " << iter << endl;
		cout << "\nH = " << endl;
		printComplexMatrix(H_loop);
		
		/*
		// this method use the eigenvalues from 2 x 2 trailing matrix at the bottom right
		// the weakness is if the eigenvalues are complex, then it will be hard since Francis method basic is using real arithmetic
		complex<double> a(1,0);
		complex<double> b = - H_loop[n-2][n-2] - H_loop[n-1][n-1];
		complex<double> c = (H_loop[n-1][n-1]*H_loop[n-2][n-2]) - (H_loop[n-1][n-2]*H_loop[n-2][n-1]);
	
		complex<double> σ1 = (-b + sqrt( (b*b) - (complex<double>(4,0)*a*c)) )/(complex<double>(2,0)*a);
		complex<double> σ2 = (-b - sqrt( (b*b) - (complex<double>(4,0)*a*c)) )/(complex<double>(2,0)*a);
		cout << "σ1 = " << σ1 << endl;
		cout << "σ2 = " << σ2 << endl;
	
		vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
		vector<vector<complex<double>>> I_σ1 = complexnumbermultiplicationComplexMatrix(I,σ1);
		vector<vector<complex<double>>> I_σ2 = complexnumbermultiplicationComplexMatrix(I,σ2);
		vector<vector<complex<double>>> H_σ1 = subtractComplexMatrices(H_loop,I_σ1);
		vector<vector<complex<double>>> H_σ2 = multiplyComplexMatrices(H_loop,I_σ2);
		vector<complex<double>> e1(n,0.0);
		e1[0] = 1.0;
		vector<complex<double>> vec_x = multiplycomplexmatrixvector(G,e1);
		vector<vector<complex<double>>> G = multiplyComplexMatrices(H_σ1,H_σ2);*/

		p_index = p-1;
		q_index = p-2;
		
		complex<double> s(real(H_loop[q_index][q_index] + H_loop[p_index][p_index]), imag(H_loop[q_index][q_index] + H_loop[p_index][p_index]) );
		complex<double> t(real( (H_loop[q_index][q_index] * H_loop[p_index][p_index]) - (H_loop[q_index][p_index] * H_loop[p_index][q_index])) , imag((H_loop[q_index][q_index] * H_loop[p_index][p_index]) - (H_loop[q_index][p_index] * H_loop[p_index][q_index])) );
		cout << "\ns = " << s << endl;	
		cout << "t = " << t << endl;	
		complex<double> x(real((H_loop[0][0]*H_loop[0][0]) + (H_loop[0][1]*H_loop[1][0]) - s*H_loop[0][0] + t), imag((H_loop[0][0]*H_loop[0][0]) + (H_loop[0][1]*H_loop[1][0]) - s*H_loop[0][0] + t) );
		complex<double> y(real(H_loop[1][0]*(H_loop[0][0] + H_loop[1][1] - s)), imag(H_loop[1][0]*(H_loop[0][0] + H_loop[1][1] - s)) );
		complex<double> z(real(H_loop[1][0]*H_loop[2][1]) , imag(H_loop[1][0]*H_loop[2][1]) );
		cout << "x = " << x << endl;	
		cout << "y = " << y << endl;	
		cout << "z = " << z << endl;	

		vector<complex<double>> e1(n,0.0);
		e1[0] = 1.0;
		vector<complex<double>> vec_x(n,0.0);
		vec_x[0] = x;
		vec_x[1] = y;
		vec_x[2] = z;
		double norm_x = complexnorm(vec_x);
		
		if(real(vec_x[0]) > 0 )
		{
			norm_x = -1*norm_x;		
		}		
		else if(real(vec_x[0]) <= 0 )
		{
			norm_x = norm_x;		
		}		

		vector<complex<double>> vec_w = scalarmultiplicationComplexVector(e1,norm_x);
		vector<complex<double>> vec_u = addComplexVectors(vec_x,vec_w);

		cout << "\nx = " << endl;
		printComplexVector(vec_x);
		cout << "\nw = " << endl;
		printComplexVector(vec_w);
		cout << "\nu = " << endl;
		printComplexVector(vec_u);
		

		vector<vector<complex<double>>> mat_u; 
		vector<vector<complex<double>>> mat_u_transpose; 
		mat_u.assign(n, std::vector<complex<double>>(1, 0.0));
		mat_u_transpose.assign(1, std::vector<complex<double>>(n, 0.0));

		for (int i = 0; i < n ; ++i) 
		{		
			mat_u[i][0] = vec_u[i];
			mat_u_transpose[0][i] =vec_u[i];
		}

		vector<vector<complex<double>>> uTu = multiplyComplexMatrices(mat_u_transpose,mat_u);
		vector<vector<complex<double>>> uuT = multiplyComplexMatrices(mat_u,mat_u_transpose);
		vector<vector<complex<double>>> P0; 
		
		double norm_u = divisiond(1,real(uTu[0][0])); // ||u||^{2}

		P0 = scalarmultiplicationComplexMatrix(uuT,norm_u);
		
		P0 = scalarmultiplicationComplexMatrix(P0,2);
		P0 = subtractComplexMatrices(I,P0);
		cout << "\nP_{" << i <<"}:" << endl;
		printComplexMatrix(P0);
		
		vector<vector<complex<double>>> PHP = multiplyComplexMatrices(H_loop,P0);
		vector<vector<complex<double>>> P0_transpose = TransposeComplexMatrix(P0);
		PHP = multiplyComplexMatrices(P0_transpose,PHP);
		cout << "\nP_{" << i <<"}^{T} H P_{" << i <<"}:" << endl;
		printComplexMatrix(PHP);	
		vector<vector<complex<double>>> H_new = multiplyComplexMatrices(I,PHP);
		
		// "Bulge chasing" sequence coding
		int p = n;
		int index = 1;
		int index2 = 0;
		for (int k = 0; k <= p-2; ++k)
		{
			vector<vector<complex<double>>> P1; 
			vector<vector<complex<double>>> P_last; 
			vector<complex<double>> vec_x1(3);
			if(k < p-3)
			{			
				vec_x1[0] = H_new[index][k];
				vec_x1[1] = H_new[index+1][k];
				vec_x1[2] = H_new[index+2][k];
				
				double norm_x1 = complexnorm(vec_x1);
			
				if(real(vec_x1[0]) > 0 )
				{
					norm_x1 = -1*norm_x1;		
				}		
				else if(real(vec_x1[0]) <= 0 )
				{
					norm_x1 = norm_x1;		
				}		

				cout << "\nx_{" << k+1 <<"}:" << endl;
				printComplexVector(vec_x1);

				vector<complex<double>> e1_3d(3,0.0);
				e1_3d[0] = 1.0;

				vector<complex<double>> vec_w1 = scalarmultiplicationComplexVector(e1_3d,norm_x1);
				vector<complex<double>> vec_v1 = addComplexVectors(vec_x1,vec_w1);

				vector<vector<complex<double>>> mat_v; 
				vector<vector<complex<double>>> mat_v_transpose; 
				mat_v.assign(3, std::vector<complex<double>>(1, 0.0));
				mat_v_transpose.assign(1, std::vector<complex<double>>(3, 0.0));

				for (int i = 0; i < 3 ; ++i) 
				{		
					mat_v[i][0] = vec_v1[i];
					mat_v_transpose[0][i] =vec_v1[i];
				}

				vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
				vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
				
				double norm_v1 = divisiond(1,real(vTv[0][0])); // ||u||^{2}

				P1 = scalarmultiplicationComplexMatrix(vvT,norm_v1);
				
				P1 = scalarmultiplicationComplexMatrix(P1,2);
				vector<vector<complex<double>>> P1_final;
				P1_final.assign(n, std::vector<complex<double>>(n, 0.0));
				int i1 = 0;
				for (int i = 1 + index2; i <= 3 + index2; ++i)
				{
					int j1=0;
					for (int j = 1 + index2; j <= 3 + index2; ++j)
					{
						P1_final[i][j] = P1[i1][j1];
						j1 = j1+1;
					}
					i1= i1 + 1;
				}
				vector<vector<complex<double>>> P1_f = subtractComplexMatrices(I,P1_final);
				cout << "\nP_{" << k+1 <<"} reflector:" << endl;
				printComplexMatrix(P1_f);

				H_new = multiplyComplexMatrices(H_new,P1_f);
				H_new = multiplyComplexMatrices(P1_f,H_new);

				cout << "\nH^{" << k+1 <<"}:" << endl;
				printComplexMatrix(H_new);
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < n; ++j)
					{
						if (abs(real(H_new[i][j])) < 1e-10)
						{
							H_new[i][j] = 0;
						}
					}
				}
				H_new[n-1][k]= 0;
				cout << "\nH^{" << k+1 <<"}:" << endl;
				printComplexMatrix(H_new);
			}
			else if (k == p-3)
			{
				vector<complex<double>> vec_x2(2);
				vec_x2[0] = H_new[index][k];
				vec_x2[1] = H_new[index+1][k];
				
				double norm_x2 = complexnorm(vec_x2);
			
				if(real(vec_x2[0]) > 0 )
				{
					norm_x2 = -1*norm_x2;		
				}		
				else if(real(vec_x2[0]) <= 0 )
				{
					norm_x2 = norm_x2;		
				}		

				cout << "\nx_{" << k+1 <<"}:" << endl;
				printComplexVector(vec_x2);

				vector<complex<double>> e1_2d(2,0.0);
				e1_2d[0] = 1.0;

				vector<complex<double>> vec_w2 = scalarmultiplicationComplexVector(e1_2d,norm_x2);
				vector<complex<double>> vec_v2 = addComplexVectors(vec_x2,vec_w2);

				vector<vector<complex<double>>> mat_v; 
				vector<vector<complex<double>>> mat_v_transpose; 
				mat_v.assign(2, std::vector<complex<double>>(1, 0.0));
				mat_v_transpose.assign(1, std::vector<complex<double>>(2, 0.0));

				for (int i = 0; i < 2 ; ++i) 
				{		
					mat_v[i][0] = vec_v2[i];
					mat_v_transpose[0][i] =vec_v2[i];
				}

				vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
				vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
				
				double norm_v2 = divisiond(1,real(vTv[0][0])); // ||u||^{2}

				P_last = scalarmultiplicationComplexMatrix(vvT,norm_v2);
				
				P_last = scalarmultiplicationComplexMatrix(P_last,2);
				//cout << "\nP_{" << k+1 <<"} last:" << endl;
				//printComplexMatrix(P_last);
				vector<vector<complex<double>>> P_last_final;
				P_last_final.assign(n, std::vector<complex<double>>(n, 0.0));
				
				vector<vector<complex<double>>> P_last_f = subtractComplexMatrices(I,P_last_final);
			
				P_last_f[n-2][n-2] = complex<double>(1,0) - P_last[0][0] ;
				P_last_f[n-2][n-1] = complex<double>(0,0) - P_last[0][1] ;
				P_last_f[n-1][n-2] = complex<double>(0,0) - P_last[1][0] ;
				P_last_f[n-1][n-1] = complex<double>(1,0) - P_last[1][1] ;

				cout << "\nP_{" << k+1 <<"} reflector:" << endl;
				printComplexMatrix(P_last_f);
				H_new = multiplyComplexMatrices(H_new,P_last_f);
				H_new = multiplyComplexMatrices(P_last_f,H_new);

				cout << "\nH^{" << k+1 <<"}:" << endl;
				printComplexMatrix(H_new);
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < n; ++j)
					{
						if (abs(real(H_new[i][j])) < 1e-10)
						{
							H_new[i][j] = 0;
						}
					}
				}
				cout << "\nH^{" << k+1 <<"}:" << endl;
				printComplexMatrix(H_new);
			}
			else if (k == p-2)
			{
				vector<complex<double>> vec_x2(2);
				vec_x2[0] = H_new[index2][k];
				vec_x2[1] = H_new[index2+1][k];
				
				double norm_x2 = complexnorm(vec_x2);
			
				if(real(vec_x2[0]) > 0 )
				{
					norm_x2 = -1*norm_x2;		
				}		
				else if(real(vec_x2[0]) <= 0 )
				{
					norm_x2 = norm_x2;		
				}		

				cout << "\nx_{" << k+1 <<"}:" << endl;
				printComplexVector(vec_x2);

				vector<complex<double>> e1_2d(2,0.0);
				e1_2d[0] = 1.0;

				vector<complex<double>> vec_w2 = scalarmultiplicationComplexVector(e1_2d,norm_x2);
				vector<complex<double>> vec_v2 = addComplexVectors(vec_x2,vec_w2);

				vector<vector<complex<double>>> mat_v; 
				vector<vector<complex<double>>> mat_v_transpose; 
				mat_v.assign(2, std::vector<complex<double>>(1, 0.0));
				mat_v_transpose.assign(1, std::vector<complex<double>>(2, 0.0));

				for (int i = 0; i < 2 ; ++i) 
				{		
					mat_v[i][0] = vec_v2[i];
					mat_v_transpose[0][i] =vec_v2[i];
				}

				vector<vector<complex<double>>> vTv = multiplyComplexMatrices(mat_v_transpose,mat_v);
				vector<vector<complex<double>>> vvT = multiplyComplexMatrices(mat_v,mat_v_transpose);
				
				double norm_v2 = divisiond(1,real(vTv[0][0])); // ||u||^{2}

				P_last = scalarmultiplicationComplexMatrix(vvT,norm_v2);
				
				P_last = scalarmultiplicationComplexMatrix(P_last,2);
				cout << "\nP_{" << k+1 <<"} last:" << endl;
				printComplexMatrix(P_last);
				vector<vector<complex<double>>> P_last_final;
				P_last_final.assign(n, std::vector<complex<double>>(n, 0.0));
				
				vector<vector<complex<double>>> P_last_f = subtractComplexMatrices(I,P_last_final);
			
				P_last_f[n-1][n-1] = complex<double>(1,0) - P_last[0][0] ;

				cout << "\nP_{" << k+1 <<"} reflector:" << endl;
				printComplexMatrix(P_last_f);
				H_new = multiplyComplexMatrices(H_new,P_last_f);
				H_new = multiplyComplexMatrices(P_last_f,H_new);

				cout << "\nH^{" << k+1 <<"}:" << endl;
				printComplexMatrix(H_new);
				for (int i = 0; i < n; ++i)
				{
					for (int j = 0; j < n; ++j)
					{
						if (abs(real(H_new[i][j])) < 1e-10)
						{
							H_new[i][j] = 0;
						}
					}
				}
				cout << "\nH^{" << k+1 <<"}:" << endl;
				printComplexMatrix(H_new);
			}
			index = index+1;
			index2 = index2+1;
		}
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				H_loop[i][j] = H_new[i][j];
			}
		}
		
		iter = iter+1;
	}
	cout << "\n****************************************************************************" << endl;
	cout << "\nEnd of " << iter - 1 << " iterations." << endl;
	cout << "\n****************************************************************************" << endl;
	
	cout << "\nH updated :" << endl;
	printComplexMatrix(H_loop);
	
	for (int i = 0; i < n ; ++i) 
	{		
		complex<double> a1(1,0);
		complex<double> b1 = - H_loop[i+1][i+1] - H_loop[i][i];
		complex<double> c1 = (H_loop[i+1][i+1]*H_loop[i][i]) - (H_loop[i+1][i]*H_loop[i][i+1]);
			
		complex<double> λ1 = (-b1 + sqrt( (b1*b1) - (complex<double>(4,0)*a1*c1)) )/(complex<double>(2,0)*a1);
		complex<double> λ2 = (-b1 - sqrt( (b1*b1) - (complex<double>(4,0)*a1*c1)) )/(complex<double>(2,0)*a1);
		vector_eigenvalues.push_back(λ1);
		vector_eigenvalues.push_back(λ2);
		i = i+1;		
	}
	//matrix_eigenvectors = TransposeComplexMatrix(matrix_eigenvectors);
	//cout << "\nEigenvectors (in column vectors): " << endl;
	//printComplexMatrix(matrix_eigenvectors);

	cout << "\nEigenvalues: " << endl;
	printComplexVector(vector_eigenvalues);

	return vector_eigenvalues; // returns the eigenvalues
}

vector<double> RayleighQuotientIteration(vector<vector<double>>& A, int iterations, double a, double b) 
{

	int n = A.size(); // Rows	
	int m = n;
	
	vector<double> vector_eigenvalues;	
	vector<vector<double>> matrix_eigenvectors;	
	vector<double> v;
	vector<vector<double>> A_original(n, vector<double>(n));	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	
	// 1. Obtain a seed:
	std::default_random_engine generator(
        std::chrono::system_clock::now().time_since_epoch().count());
	
	// 2. Define the distribution for floating-point numbers (e.g., uniform distribution)
	// Range [a, b)
	std::uniform_real_distribution<double> distribution(a, b);	
	std::uniform_real_distribution<double> distribution2(0, b);
	double mu = distribution(generator); 
	double sigma = distribution2(generator);
	// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
	// mu and sigma are also random number that are generated with uniform distribution
	v = vrandn_normal(mu, sigma, m);
	double norm_v = divisiond(1,norm(v));
	//cout << "\nv initial vector = " << endl;
	//printComplexVector(v);
	v = scalarmultiplication(v,norm_v);
	cout << "\nv initial vector_{normalize} = " << endl;
	printVector(v);
	
	// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
	double cdotproduct = dot(v,v);
	double lambda = quadraticmultiplication(A,v) / cdotproduct;
	cout << "\nλ^{0} = " << lambda << endl;

	// Solve (A - λ^{k-1}I)w = v^{k-1}
	// Then multiply the reverse to obtain RQ = A^{k}
	for (int i = 0; i < iterations; ++i) 
	{		
		vector<vector<double>> I = createIdentityMatrix(n);
		vector<vector<double>> lambda_I = scalarmultiplication(I,lambda); 
		vector<vector<double>> A_lambda_I = matrixsubtraction(A_original,lambda_I);
		vector<vector<double>> A_lambda_I_inverse = inverse(A_lambda_I);
	
		vector<double> w = multiplymatrixvector(A_lambda_I_inverse,v);

		double norm_w = divisiond(1,norm(w));
		w = scalarmultiplication(w,norm_w);
		//cout <<"\nw_{normalize} = (A - λ^{k-1}I) * v^{k-1}: "<<endl;
		//printVector(w);
	
		for (int i = 0; i < n ; ++i) 
		{		
			v[i] = w[i];
		}
	}
	// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
	double dotproduct = dot(v,v);
	double lambda_final = divisiond(quadraticmultiplication(A,v),dotproduct);
	cout << "\nEigenvalue = " << lambda_final << endl;
	cout <<"\nCorresponding Eigenvector: "<<endl;
	printVector(v);
	vector_eigenvalues.push_back(lambda);

	//cout << "\nEigenvector corresponds to " << lambda << " :" << endl;
	//printComplexVector(v);
	/* */
	return vector_eigenvalues; // returns the eigenvector
}

vector<complex<double>> ComplexRayleighQuotientIteration(vector<vector<complex<double>>> &A, int iterations, int N, double a, double b) 
{

	int n = A.size(); // Rows	
	int m = n;

	cout << "\nMatrix A : " << endl;
	printComplexMatrix(A);

	vector<complex<double>> vector_eigenvalues(n);		
	vector<vector<complex<double>>> matrix_eigenvectors;	
	vector<complex<double>> v;
	vector<vector<complex<double>>> A_original;	
	A_original.assign(n, std::vector<complex<double>>(n, 0.0));

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	int index = 0;
	//  N =  number of the outermost iteration to obtain the correct eigenvalues.
	for (int iter = 0; iter < N; ++iter) 
	{
		cout <<"\n************************************************************************" << endl;	
		cout <<"\nIter: "<< iter << endl;	
	
		// 1. Obtain a seed:
		std::default_random_engine generator(
		std::chrono::system_clock::now().time_since_epoch().count());
		
		// 2. Define the distribution for floating-point numbers (e.g., uniform distribution)
		// Range [a, b)
		std::uniform_real_distribution<double> distribution(a, b);	
		std::uniform_real_distribution<double> distribution2(0, b);
		double mu = distribution(generator); 
		double sigma = distribution2(generator);
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		// mu and sigma are also random number that are generated with uniform distribution
		v = complexvecrand_normal(mu, sigma, m);
		double norm_v = divisiond(1,complexnorm(v));

		v = scalarmultiplicationComplexVector(v,norm_v);
		cout << "\nv initial vector_{normalize} = " << endl;
		printComplexVector(v);
		
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		complex<double> cdotproduct = complexdotproduct(v,v);
		complex<double> lambda = complexquadraticmultiplication(A,v) / cdotproduct;
		cout << "\nλ^{0} = " << lambda << endl;

		// Solve (A - λ^{k-1}I)w = v^{k-1}
		// Then multiply the reverse to obtain RQ = A^{k}
		for (int i = 0; i < iterations; ++i) 
		{		
			vector<vector<complex<double>>> I = createIdentityComplexMatrix(n);
			vector<vector<complex<double>>> lambda_I = complexnumbermultiplicationComplexMatrix(I,lambda); 
			vector<vector<complex<double>>> A_lambda_I = subtractComplexMatrices(A_original,lambda_I);
			vector<vector<complex<double>>> A_lambda_I_inverse = ComplexMatrixInverse(A_lambda_I);
		
			vector<complex<double>> w = multiplycomplexmatrixvector(A_lambda_I_inverse,v);

			double norm_w = divisiond(1,complexnorm(w));
			w = scalarmultiplicationComplexVector(w,norm_w);
			//cout <<"\nw_{normalize} = (A - λ^{k-1}I) * v^{k-1}: "<<endl;
			//printComplexVector(w);
		
			for (int i = 0; i < n ; ++i) 
			{		
				v[i] = w[i];
			}
		}
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		complex<double> dotproduct_final = complexdotproduct(v,v);
		complex<double> lambda_final = complexquadraticmultiplication(A,v) / dotproduct_final;
		complex<double> eigenvalue(lround(real(lambda_final)), lround(imag(lambda_final)));
		cout << "\nEigenvalue = " << eigenvalue << endl;
		cout <<"\nCorresponding Eigenvector: "<<endl;
		printComplexVector(v);

		vector<complex<double>> A_x =  multiplycomplexmatrixvector(A,v);  
		vector<complex<double>> lambda_x =  complexnumbermultiplicationComplexVector(v,eigenvalue); 
		cout <<"\nAx: "<<endl;
		printComplexVector(A_x);
		cout <<"\nλx: "<<endl;
		printComplexVector(lambda_x);	
		
		// to obtain the whole eigenvalues and eigenvectors
		vector<complex<double>> difference = subtractComplexVectors(A_x,lambda_x);
		double diffnorm = complexnorm(difference);
		if (abs(diffnorm) < 1e-12 )
		{
			for (int i = 0; i < n ; ++i)		
			{		
				if(eigenvalue != vector_eigenvalues[i])
				{
					if(i == n-1)	
					{
						cout << "\nEigenvalue accepted= " << eigenvalue << endl;
						vector_eigenvalues[index] = eigenvalue;
						matrix_eigenvectors.push_back(v);
						index= index+1;
					}
				}			
				else if(eigenvalue == vector_eigenvalues[i])
				{
					break;
				}			
			}		
		}	
		if (index == n)
		{
			iter = N-1;
		}
		
	}
	cout <<"\nEigenvalues : "<<endl;	
	printComplexVector(vector_eigenvalues);
	cout <<"\nEigenvectors (represented by row): "<<endl;	
	printComplexMatrix(matrix_eigenvectors);

	//vector<vector<complex<double>>> matrix_eigenvectors_transpose = TransposeComplexMatrix(matrix_eigenvectors);
	//printComplexMatrix(matrix_eigenvectors_transpose);
	
	return vector_eigenvalues; // returns the eigenvector
}

vector<double> diagonalization(vector<vector<double>>& A) 
{
	int iterations = 100;
	int n = A.size();
	int m = n;
	vector<double> vector_eigenvalues;	
	vector<vector<double>> matrix_eigenvectors;	
	vector<double> v;
	vector<vector<double>> A_original(n, vector<double>(n,0.0));	
	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	cout <<"\nA: "<<endl;
	printMatrix(A);
	for (int k = 0; k < n; ++k)
	{
		double mu = 2; 
		double sigma = 0.5;
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		v = vrandn_normal(mu, sigma, m);

		// Repeatedly multiply the matrix A by the current vector v to get v_{k+1} = A v_{k}
		// Then normalize v_{k+1} at each step to prevent the values from growing too large.
		// The normalized vector converges to the dominant eigenvector, and the corresponding eigenvalue can be calculated from Av = λv
		for (int i = 0; i < iterations; ++i) 
		{
			vector<double> w = multiplymatrixvector(A,v);
			double ssw = 0; // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w[j]*w[j];
			}
			double norm = std::sqrt(ssw);
			// Normalize w (w = w / norm)
			vector<double> w_normalize = scalarmultiplication(w,divisiond(1,norm));
			
			for (int j = 0; j < n ; ++j)
			{
				v[j] = w_normalize[j];
			}
		}
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		double lambda = divisiond(quadraticmultiplication(A,v),dot(v,v));
	
		vector_eigenvalues.push_back(lambda);

		vector<vector<double>> P = PermutationMatrixMax(v);

		vector<double> wi = multiplymatrixvector(P,v);
		//cout << "\nw_{1} = " << endl;
		//printVector(wi);

		vector<vector<double>> R_inv(m, vector<double>(m));
		for (int i = 0; i < m ; ++i)
		{
			R_inv[i][0] = wi[i];
		}
		for (int i = 1; i < m ; ++i)
		{
			for (int j = 1; j < m ; ++j)
			{
				if (i != j)
				{
					R_inv[i][j] = 0;
				}
				else if (i == j)
				{
					R_inv[i][j] = 1;
				}
			}
		}
		vector<vector<double>> R= inverse(R_inv);
		vector<vector<double>> RP = multiply(R,P);
		vector<vector<double>> PiRi = multiply(P,R_inv);
		vector<vector<double>> APiRi = multiply(A,PiRi);
		vector<vector<double>> RPAPiRi = multiply(RP,APiRi);
		
		A = deleteColumn(RPAPiRi,0); 
		A = deleteRow(RPAPiRi,0); 
		m = m-1;
	}
	
	// to determine the whole eigenvectors from known eigenvalue
	// we use inverse iteration with shift from power method
	for (int i = 0; i < n ;++i)
	{
		double eigen = vector_eigenvalues[i];
		vector<double> vector_eigen;
		cout << "\nEigenvalue(" << i+1 << ") = " << eigen << endl;
		
		vector<vector<double>> I = createIdentityMatrix(n);
		// why we use subtraction to 0.1? because we are using inverse iteration method and the inverse of (A - λI) does not exist 
		vector<vector<double>> lambda_I = scalarmultiplication(I,eigen-0.1); 
		vector<vector<double>> A_lambda_I = matrixsubtraction(A_original,lambda_I);
		vector<vector<double>> A_lambda_I_inverse = inverse(A_lambda_I);
		
		double mu = 2; 
		double sigma = 0.5;
		vector_eigen = vrandn_normal(mu, sigma, n);
		double ssv = 0; // sum of squares of w.
		for (int j = 0; j < n ; ++j)
		{
			ssv += vector_eigen[j]*vector_eigen[j];
		}
		double norm = std::sqrt(ssv);
		// Normalize v (v = v / norm)
		vector_eigen = scalarmultiplication(vector_eigen,divisiond(1,norm));
		
		for (int k = 0; k < iterations; ++k) 
		{
			vector<double> w_eigen = multiplymatrixvector(A_lambda_I_inverse,vector_eigen);
			double ssw = 0; // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w_eigen[j]*w_eigen[j];
			}
			double norm_w = std::sqrt(ssw);
			// Normalize w (w = w / norm)
			vector<double> w_eigen_normalize = scalarmultiplication(w_eigen,divisiond(1,norm_w));
			
			for (int j = 0; j < n ; ++j)
			{
				vector_eigen[j] = w_eigen_normalize[j];
			}
			
		}
		//cout << "\nCorresponding Eigenvector: " << endl;
		//printVector(vector_eigen);
		matrix_eigenvectors.push_back(vector_eigen);
	}

	matrix_eigenvectors = transpose(matrix_eigenvectors);
	cout << "\nP: " << endl;
	printMatrix(matrix_eigenvectors);
	vector<vector<double>> P_mat (n, vector<double>(n,0.0));
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n ; ++j)
		{
			P_mat[i][j] = matrix_eigenvectors[i][j];
		}		
	}
	vector<vector<double>> P_mat_inv = inverse(P_mat);
	vector<vector<double>> AP = multiply(A_original,P_mat);
	vector<vector<double>> PiAP = multiply(P_mat_inv,AP);

	cout << "\nP^{-1}: " << endl;
	printMatrix(P_mat_inv);
	cout << "\nP^{-1}AP: " << endl;	
	printMatrix(PiAP);
	
	return vector_eigenvalues; // returns the eigenvector
}

vector<double>SVD(vector<vector<double>>& A0) 
{
	int iterations = 100;
	int row_A = A0.size();
	int col_A = A0[0].size();
	vector<double> vector_eigenvalues;	
	vector<vector<double>> matrix_eigenvectors;	
	vector<double> v;
	vector<vector<double>> A;
	vector<vector<double>> A_transpose = transpose(A0);
	A = multiply(A_transpose,A0);
	int n = A.size();
	int m = n;
	vector<vector<double>> A_original(n, vector<double>(n,0.0));
	
	// SVD decomposition matrices	
	vector<vector<double>> sigma(row_A, vector<double>(col_A,0.0));	
	vector<vector<double>> U;	
	vector<vector<double>> V_transpose(col_A, vector<double>(col_A,0.0));	

	for (int i = 0; i < n ; ++i) 
	{		
		for (int j = 0; j < n ; ++j)
		{
			A_original[i][j] = A[i][j];
		}
	}
	
	cout <<"\nA: "<<endl;
	printMatrix(A0);
	cout <<"\nA^{T}A: "<<endl;
	printMatrix(A);
	for (int k = 0; k < n; ++k)
	{
		double mu = 2; 
		double sigma = 0.5;
		// start with a random initial vector v that is normally distributed with mean = mu and standard deviation = sigma
		v = vrandn_normal(mu, sigma, m);

		// Repeatedly multiply the matrix A by the current vector v to get v_{k+1} = A v_{k}
		// Then normalize v_{k+1} at each step to prevent the values from growing too large.
		// The normalized vector converges to the dominant eigenvector, and the corresponding eigenvalue can be calculated from Av = λv
		for (int i = 0; i < iterations; ++i) 
		{
			vector<double> w = multiplymatrixvector(A,v);
			double ssw = 0; // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w[j]*w[j];
			}
			double norm = std::sqrt(ssw);
			// Normalize w (w = w / norm)
			vector<double> w_normalize = scalarmultiplication(w,divisiond(1,norm));
			
			for (int j = 0; j < n ; ++j)
			{
				v[j] = w_normalize[j];
			}
		}
		// Calculate eigenvalue (lambda = (v^T * A * v) / (v^T * v))
		double lambda = divisiond(quadraticmultiplication(A,v),dot(v,v));
		vector_eigenvalues.push_back(lambda);

		vector<vector<double>> P = PermutationMatrixMax(v);

		vector<double> wi = multiplymatrixvector(P,v);
		//cout << "\nw_{1} = " << endl;
		//printVector(wi);

		vector<vector<double>> R_inv(m, vector<double>(m));
		for (int i = 0; i < m ; ++i)
		{
			R_inv[i][0] = wi[i];
		}
		for (int i = 1; i < m ; ++i)
		{
			for (int j = 1; j < m ; ++j)
			{
				if (i != j)
				{
					R_inv[i][j] = 0;
				}
				else if (i == j)
				{
					R_inv[i][j] = 1;
				}
			}
		}
		vector<vector<double>> R= inverse(R_inv);
		vector<vector<double>> RP = multiply(R,P);
		vector<vector<double>> PiRi = multiply(P,R_inv);
		vector<vector<double>> APiRi = multiply(A,PiRi);
		vector<vector<double>> RPAPiRi = multiply(RP,APiRi);
		
		A = deleteColumn(RPAPiRi,0); 
		A = deleteRow(RPAPiRi,0); 
		m = m-1;
	}
	
	// to determine the whole eigenvectors from known eigenvalue
	// we use inverse iteration with shift from power method
	for (int i = 0; i < n ;++i)
	{
		double eigen = vector_eigenvalues[i];
		vector<double> vector_eigen;
		//cout << "\nEigenvalue(" << i+1 << ") = " << eigen << endl;
		
		vector<vector<double>> I = createIdentityMatrix(n);
		// why we use subtraction to 0.1? because we are using inverse iteration method and the inverse of (A - λI) does not exist 
		vector<vector<double>> lambda_I = scalarmultiplication(I,eigen-0.1); 
		vector<vector<double>> A_lambda_I = matrixsubtraction(A_original,lambda_I);
		vector<vector<double>> A_lambda_I_inverse = inverse(A_lambda_I);
		
		double mu = 2; 
		double sigma = 0.5;
		vector_eigen = vrandn_normal(mu, sigma, n);
		double ssv = 0; // sum of squares of w.
		for (int j = 0; j < n ; ++j)
		{
			ssv += vector_eigen[j]*vector_eigen[j];
		}
		double norm = std::sqrt(ssv);
		// Normalize v (v = v / norm)
		vector_eigen = scalarmultiplication(vector_eigen,divisiond(1,norm));
		
		for (int k = 0; k < iterations; ++k) 
		{
			vector<double> w_eigen = multiplymatrixvector(A_lambda_I_inverse,vector_eigen);
			double ssw = 0; // sum of squares of w.
			for (int j = 0; j < n ; ++j)
			{
				ssw += w_eigen[j]*w_eigen[j];
			}
			double norm_w = std::sqrt(ssw);
			// Normalize w (w = w / norm)
			vector<double> w_eigen_normalize = scalarmultiplication(w_eigen,divisiond(1,norm_w));
			
			for (int j = 0; j < n ; ++j)
			{
				vector_eigen[j] = w_eigen_normalize[j];
			}
			
		}
		matrix_eigenvectors.push_back(vector_eigen);
	}

	matrix_eigenvectors = transpose(matrix_eigenvectors);
	
	V_transpose = transpose(matrix_eigenvectors);
	
	for (int i = 0; i < row_A ; ++i) 
	{		
		for (int j = 0; j < col_A ; ++j)
		{
			if( i == j)
			{
				sigma[i][j] = sqrt(vector_eigenvalues[i]);
			}
			else if( i != j)
			{
				sigma[i][j] = 0;
			}
		}
	}
	for (int i = 0; i < col_A ; ++i) 
	{		
		double k = divisiond(1,sqrt(vector_eigenvalues[i]));
		vector<double> vi = getColumn(matrix_eigenvectors,i);
		vector<double> ui = multiplymatrixvector(A0,vi);
		ui = scalarmultiplication(ui,k);
		U.push_back(ui);
	}
	

	vector<vector<double>> new_U = normalizehomogeneouslinearsystembasis_matrix(U);
	int nu = new_U.size();
	for (int i = 0; i < nu; ++i) 
	{	
		U.push_back(getRow(new_U,i));
	}
	U = transpose(U);

	cout << "\nU: " << endl;	
	printMatrix(U);
	cout << "\nΣ: " << endl;	
	printMatrix(sigma);
	cout << "\nV^{T}: " << endl;	
	printMatrix(V_transpose);
	
	vector<vector<double>> sigmaV = multiply(sigma,V_transpose);
	vector<vector<double>> UsV = multiply(U,sigmaV);
	cout << "\nUΣV^{T}: " << endl;	
	printMatrix(UsV);
	return vector_eigenvalues; // returns the eigenvector
}

#endif
#endif