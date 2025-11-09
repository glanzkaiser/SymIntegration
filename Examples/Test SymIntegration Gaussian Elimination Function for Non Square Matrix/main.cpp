// g++ -o result main.cpp -lsymintegration
// Merci beaucoup Freya et Sentinel

#include<bits/stdc++.h>
#include<iostream>
#include "symintegrationc++.h"
#include<vector>
#include <chrono>

using namespace std::chrono;
using namespace std;

#define R 3 // number of rows
#define C 6 // number of columns

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

 	// Construct a symbolic matrix of size 3 X 6
	Matrix<Symbolic> B_mat(3,6);
	B_mat[0][0] = 1;		B_mat[0][1] = 1;	B_mat[0][2] = 1;	B_mat[0][3] = 0;	B_mat[0][4] = 0;	B_mat[0][5] = 12;
	B_mat[1][0] = 2;		B_mat[1][1] = 1;	B_mat[1][2] = 0;	B_mat[1][3] = 1;	B_mat[1][4] = 0;	B_mat[1][5] = 16;
	B_mat[2][0] = -40;		B_mat[2][1] = -30;	B_mat[2][2] = 0;	B_mat[2][3] = 0;	B_mat[2][4] = 1;	B_mat[2][5] = 0;
	cout << "The matrix A:\n" << B_mat <<endl;
	cout << endl;
	
	
	cout << "Gaussian Elimination :\n" << gaussianeliminationtest(B_mat,C,R) <<endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}