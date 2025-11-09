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

	for(int i = 0; i < R; i++)
	{
		for(int j=0; j < C; j++)
		{
			cout << doubleMatrix2[i][j] << endl;
		}
	}


	Matrix<Symbolic> B_mat(R,C);
	for(int i = 0; i < R; i++)
	{
		for(int j=0; j < C; j++)
		{
			B_mat[i][j] = doubleMatrix2[i][j] ;
		}
	}
	cout << "\nThe matrix A:\n" << B_mat <<endl;
	cout << "\nRows of matrix A:\n" << B_mat.rows() <<endl;
	cout << "\nColumns of matrix A:\n" << B_mat.cols() <<endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}