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
#ifndef SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS_DEFINE

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f
#define pi  3.1415926535897

#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#include <string>
#include <fstream> // For file operations
#include <chrono>
#include <string>
#include<iostream>
#include<vector>
#include <cmath>
#include <bits/stdc++.h> //for setw(6) 
#include <iomanip> // to declare the manipulator of setprecision()
#include <map>
#include <algorithm> // For std::max_element,  std::sort, std::reverse ,  std::for_each
#include <numeric> // For std::accumulate
#include <random> // For random number generation
#include <complex>

vector<string> loadStringVector(const string& filename) 
{
	vector<string> vecx;
	ifstream inputFile(filename);

	if (!inputFile.is_open()) 
	{
		cerr << "Error: Could not open file " << filename << endl;
		return vecx; // Return empty matrix on error
	}

	if (inputFile.is_open()) 
	{
		string line;
		while (getline(inputFile,line)) 
		{ // Reads numbers separated by whitespace
			vecx.push_back(line);
		}
		inputFile.close();
		
	} 
	inputFile.close();
	return vecx;
}

vector<vector<string>> loadStringMatrix(const string& filename) 
{
	vector<vector<string>> matrix;
	ifstream inputFile(filename);

	if (!inputFile.is_open()) 
	{
		cerr << "Error: Could not open file " << filename << endl;
		return matrix; // Return empty matrix on error
	}

	if (inputFile.is_open()) 
	{
		string line;
		while (getline(inputFile, line)) 
		{
		if (line.empty()) 
		{ // Skip empty lines
		    continue;
		}
		istringstream iss(line);
		vector<string> row;
		string value;
		while (iss >> value) 
		{
			row.push_back(value);
		}
		matrix.push_back(row);
		}
		inputFile.close();
	} 
	return matrix;
}

void printStringVector(vector<string> vectorx)
{
	//int n = vectorx.size();
	for (const std::string& s : vectorx) 
	{
		cout << s << "";
	}
	cout << endl;
}

void printStringMatrix(vector<vector<string>> matrix) 
{
	 for (const auto& row : matrix) 
	{
		for (const auto& s : row) 
		{
			cout << s << " ";
		}
		cout << endl;
	}
}

vector<vector<string>> transposeStringMatrix(vector<vector<string>> &matrix) // Done beautifully
{
	if (matrix.empty() || matrix[0].empty() ) 
	{
        	// Handle empty matrices or invalid dimensions
		return {};
	}	
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();

	vector<vector<string>> transposed_matrix(cols, vector<string>(rows, ""));
	for (size_t i = 0; i < rows; ++i) 
	{
		for (size_t j = 0; j < cols; ++j) 
		{
			transposed_matrix[j][i] = matrix[i][j];	
		}
	}
	return transposed_matrix;
}

void stringmatrix_addColumn(vector<vector<string>>& data_table, const vector<vector<string>>& new_column_data) 
{
	// Check if the number of rows is compatible
	if (data_table.size() != new_column_data.size()) 
	{
		cerr << "Error: Number of rows in both vectors must be the same to add as a column." << endl;
		return;
	}

	// Iterate through each row and append the corresponding element(s)
	for (size_t i = 0; i < data_table.size(); ++i) 
	{
		// Append all elements from the inner vector of new_column_data[i] 
		// to the inner vector of data_table[i]
		data_table[i].insert(data_table[i].end(), new_column_data[i].begin(), new_column_data[i].end()); 
	}
}

// Function to extract a column vector
// This approach is efficient as it uses standard C++ algorithms like std::transform and std::back_inserter. 
// It also includes basic bounds checking to ensure safe access to the inner vectors. 
vector<string> stringmatrix_getColumn(vector<vector<string>>& matrix, int columnIndex) 
{
/**
 * Extracts a specific column from a 2D vector of strings.
 *
 * @param matrix The 2D vector (vector of vectors) to extract from.
 * @param columnIndex The index of the column to retrieve (0-based).
 * @return A new vector<string> containing the elements of the specified column.
 */
	vector<string> columnVector;
	int R = matrix.size();
	// Check for an empty matrix or invalid column index
	if (matrix.empty() || columnIndex < 0) 
	{
		return columnVector; // Return an empty vector
	}

	// Reserve space to avoid frequent reallocations if the matrix is large
	columnVector.reserve(R);

	// Use std::transform to efficiently iterate and extract the element at columnIndex
	std::transform(matrix.begin(), matrix.end(), std::back_inserter(columnVector),
        [columnIndex](const vector<string>& row) 
	{
	int R_row = row.size();
	// Check if the row has enough elements to prevent out-of-bounds access
	if (columnIndex < R_row) 
	{
		return row[columnIndex];
	} 
	else 
	{
		// Return an empty string or handle the error as appropriate for your use case
		return std::string("");
	}
	});

	return columnVector;	
}

void stringmatrix_deleteColumn(vector<vector<string>>& matrix, int columnIndex) 
{
	int C = matrix[0].size();
	if (matrix.empty()) 
	{
		cerr << "Error: Matrix is empty" << endl;
	}
	// Ensure the rowIndex is valid
	if (columnIndex < 0 || columnIndex >= C) 
	{
		cerr << "Error: Invalid column index." << endl;		
	}
	// Iterate through each row in the 2D vector
	for (auto& row : matrix) 
	{
		int R = row.size();
		// Check if the column index is valid for the current row
		if (columnIndex >= 0 && columnIndex < R) 
		{
			// Erase the element at the specified column index using an iterator
			row.erase(row.begin() + columnIndex);
		}
	}
}

void stringmatrix_addRow(vector<vector<string>>& matrix, vector<string>& newRow, int rowIndex) 
{
	int R = matrix.size();
	
	if (matrix.empty()) 
	{
		cerr << "Error: Matrix is empty" << endl;
	}
	// Ensure the rowIndex is valid
	if (rowIndex < 0 || rowIndex >= R) 
	{
		cerr << "Error: Invalid row index." << endl;		
	}
	else if (rowIndex >= 0 && rowIndex < R) 
	{
		// Inserting at matrix.begin() + matrix.size() is equivalent to push_back().
		// Alternative: If adding at the end, use matrix.push_back(newRow)
		matrix.insert(matrix.begin() + rowIndex, newRow); 
	}
	
}

// Function to extract a row vector
vector<string> stringmatrix_getRow(vector<vector<string>>& matrix, int rowIndex) 
{
	vector<string> rowVector;

	int C = matrix[0].size();
	
	// Check for an empty matrix or invalid row index
	if (matrix.empty() || rowIndex < 0) 
	{
		return rowVector; // Return an empty vector
	}
	for(int i=0; i<C; i++)
	{
		rowVector.push_back(static_cast<string>(matrix[rowIndex][i])); 
	}
	
	return rowVector;
}

void stringmatrix_deleteRow(vector<vector<string>>& matrix, int rowIndex) 
{
	int R = matrix.size();
	
	if (matrix.empty()) 
	{
		cerr << "Error: Matrix is empty" << endl;
	}
	// Ensure the rowIndex is valid
	if (rowIndex < 0 || rowIndex >= R) 
	{
		cerr << "Error: Invalid row index." << endl;		
	}
	else if (rowIndex >= 0 && rowIndex < R) 
	{
		matrix.erase(matrix.begin() + rowIndex);
	}
	
}

void stringmatrix_removeduplicaterows_sort(vector<vector<string>>& matrix) 
{
	// This method changes the original order of the rows.
	// 1. Sort the vector. std::vector has an overloaded operator< 
	//    that works on inner vectors and strings by default.
	std::sort(matrix.begin(), matrix.end());

	// 2. Use std::unique to move all non-duplicate elements to the front
	//    and return an iterator to the new logical end of the unique range.
	auto last = std::unique(matrix.begin(), matrix.end());

	// 3. Erase the duplicate elements from the end of the vector.
	matrix.erase(last, matrix.end());
}

vector<vector<string>> transpose_stringmatrix(const vector<vector<string>>& data) 
{
	if (data.empty()) 
	{
		return {};
	}

	// Assuming all inner vectors have the same size (like a table)
	size_t num_rows = data.size();
	size_t num_cols = data[0].size();
    
	// Create a new 2D vector for the "split columns" (transposed data)
	vector<vector<string>> split_columns(num_cols, vector<string>(num_rows));

	for (size_t i = 0; i < num_rows; ++i) 
	{
		for (size_t j = 0; j < num_cols; ++j) 
		{
			split_columns[j][i] = data[i][j];
		}
	}
	return split_columns;
}

vector<string> stringvector_groupeverynchars(const vector<string>& input_vector, int n) 
{
/*
To group a vector<string> every 'n' characters, it is necessary to clarify what "grouping" entails in this context. It could mean:

    Creating new strings by concatenating parts of existing strings:
    This involves taking 'n' characters from the first string, then 'n' from the next, and so on, to form new strings in a new vector.
    Splitting existing strings into segments of 'n' characters:
    This would result in each original string being broken down into multiple smaller strings, which are then collected into a new vector. 

This code is splitting existing strings into segments of 'n' characters
*/
	vector<string> result_vector;
	for (const string& s : input_vector) 
	{
		int l = s.length();
		for (int i = 0; i < l ; i += n) 
		{
			result_vector.push_back(s.substr(i, n));
		}
	}
    return result_vector;
}


#endif
#endif