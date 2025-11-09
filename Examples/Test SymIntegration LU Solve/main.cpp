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
	
	cout << "Data:" << endl;
	printMatrix(doubleMatrix);
	
	dvec x;
	LUsolve_nhsystem(doubleMatrix,vecb,x);
	cout << "\nx : " << endl;
	printVector(x);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}