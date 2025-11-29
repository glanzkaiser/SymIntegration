// g++ -o result main.cpp -lsymintegration
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

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	//string filename = "matrix.txt";
	
	dmat doubleMatrix = loadMatrixFromFile("matrix.txt");
	dmat vecb = loadMatrixFromFile("vectorb.txt");
	int n = doubleMatrix.size();
	int R = doubleMatrix.size();
	int C = doubleMatrix[0].size();

	// Declare the matrix L and U with the size.
	dmat L(n, vector<double>(n));
	dmat U(n, vector<double>(n));
	
	cout << "Data:" << endl;
	for(int i = 0; i < R; i++)
	{
		for(int j=0; j < C; j++)
		{
			cout << setw(10) << doubleMatrix[i][j] << setw(10);
		}
		cout << endl;
	}
	LUDecomposition(doubleMatrix, L, U);
	
	cout << "\nMatrix A : " << endl;
	printMatrix(doubleMatrix);
	cout << "\nMatrix L : " << endl;
	printMatrix(L);
	cout << "\nMatrix U : " << endl; 
	printMatrix(U);
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}