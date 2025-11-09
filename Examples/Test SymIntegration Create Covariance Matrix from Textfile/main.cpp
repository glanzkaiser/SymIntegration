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
	
	vector<vector<double>> doubleMatrix2 = loadMatrixFromFile("matrix.txt");
	int R = doubleMatrix2.size();
	int C = doubleMatrix2[0].size();
	cout << "Data:" << endl;
	for(int i = 0; i < R; i++)
	{
		for(int j=0; j < C; j++)
		{
			cout << setw(10) << doubleMatrix2[i][j] << setw(10);
		}
		cout << endl;
	}
	SymbolicMatrix y = covariancematrix(doubleMatrix2);
	
	cout << y << endl;
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}