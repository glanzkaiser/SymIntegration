// g++ -o result main.cpp -lsymbolicc++
// Merci beaucoup Freya et Sentinel

#include<bits/stdc++.h>
#include<iostream>
#include "symintegrationc++.h"
#include<vector>
#include <chrono>
#include <algorithm> // For std::next_permutation
#include <string>
using namespace std::chrono;
using namespace std;

#define R 2 // number of rows
#define C 7 // number of columns

// Function to load a matrix from a text file
template <typename T>
vector<vector<T>> loadMatrixFromFile(const string& filename) 
{
	vector<vector<T>> matrix;
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
	vector<T> row;
	T value;
        while (iss >> value) 
	{
		row.push_back(value);
	}
	matrix.push_back(row);
	}

	inputFile.close();
	return matrix;
}

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	string filename = "matrix.txt";
	
	vector<vector<double>> doubleMatrix = loadMatrixFromFile<double>(filename);
	if (!doubleMatrix.empty()) 
	{
		cout << "\nLoaded Double Matrix:" << endl;
		for (const auto& row : doubleMatrix) 
		{
			for (double val : row) 
			{
				cout << val << " ";
				
				
			}
			cout << endl;
		}
	}
	// Construct a symbolic matrix of size R X C
	Matrix<Symbolic> B_mat(R,C);
	for(int i = 0; i < R; i++)
	{
		for(int j=0; j < C; j++)
		{
			B_mat[i][j] = doubleMatrix[i][j] ;
		}
	}
	cout << "\nThe matrix A:\n" << B_mat <<endl;
	
	vector<double> obj_function = {3,10,1,2,3}; // 3x1 + 10x2+ x3 + 2x4 + 3x5

	cout << "Simplex Method :\n" << simplexmethod(B_mat, obj_function, C, R, 5) <<endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}