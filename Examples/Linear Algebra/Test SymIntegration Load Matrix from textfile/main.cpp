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

#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>

#define π 3.1415926535897f

using namespace std;
using namespace SymbolicConstant;
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

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

int main() {
	// Example usage:
	// Create a sample matrix file (e.g., "matrix.txt")
	// 1 2 3
	// 4 5 6
	// 7 8 9

	string filename = "matrix.txt";
	vector<std::vector<int>> intMatrix = loadMatrixFromFile<int>(filename);

	if (!intMatrix.empty()) 
	{
	cout << "Loaded Integer Matrix:" << endl;
	for (const auto& row : intMatrix) 
	{
		for (int val : row) 
		{
			cout << val << " ";
		}
		cout << endl;
	}
	}
	int n = 3;
	Matrix<Symbolic> B_mat(n,n);
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
	for(int i = 0; i < n; i++)
	{
		for(int j=0; j < n; j++)
		{
			B_mat[i][j] = doubleMatrix[i][j] ;
		}
	}
	cout << "A :\n" << B_mat <<endl;

    return 0;
}
