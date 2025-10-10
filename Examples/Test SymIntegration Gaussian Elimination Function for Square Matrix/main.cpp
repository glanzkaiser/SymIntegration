// g++ -o result main.cpp -lsymbolicc++
// Merci beaucoup Freya et Sentinel

#include<bits/stdc++.h>
#include<iostream>
#include "symintegrationc++.h"
#include<vector>
#include <chrono>
#include <string>
#include <iomanip> // For std::setprecision
using namespace std::chrono;
using namespace std;

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

 	// Construct a symbolic matrix of size 3 X 6
	Matrix<Symbolic> B_mat(2,3);
	B_mat[0][0] = 0.0333;		B_mat[0][1] = 0.2;	B_mat[0][2] = 1;
	B_mat[1][0] = 0.1;		B_mat[1][1] = 0.4;	B_mat[1][2] = 1;
	cout << "The matrix A:\n" << B_mat <<endl;
	cout << endl;
	
	
	cout << "Gaussian Elimination :\n" << gaussianelimination(B_mat,2) <<endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}