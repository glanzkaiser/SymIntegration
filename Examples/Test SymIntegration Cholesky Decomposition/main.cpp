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

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	//string filename = "matrix.txt";
	
	vector<vector<double>> doubleMatrix = loadMatrixFromFile("matrix.txt");
	int R = doubleMatrix.size();
	int C = doubleMatrix[0].size();
	cout << "Data:" << endl;
	for(int i = 0; i < R; i++)
	{
		for(int j=0; j < C; j++)
		{
			cout << setw(10) << doubleMatrix[i][j] << setw(10);
		}
		cout << endl;
	}
	choleskyDecomposition(doubleMatrix);
	
	//printMatrix(doubleMatrix);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}