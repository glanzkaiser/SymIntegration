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
	
	vector<double> inputVector = loadVectorFromFile("matrix.txt");

	
	vector<vector<vector<double>>> Matrix3D = Create3DMatrixfromVector(inputVector, 2, 3, 3);
	cout << "\nMatrix A : " << endl;
	print3DMatrix(Matrix3D);
	
	vector<vector<vector<double>>> upsampleMatrix = CNN_2DUpsample3DMatrix(Matrix3D,2,2);
	cout << "\nUpsample Matrix A : " << endl;
	print3DMatrix(upsampleMatrix);
	
	vector<vector<vector<double>>> avgupsampleMatrix = CNN_2DaverageUpsample3DMatrix(Matrix3D,2,2);
	cout << "\nAverage Upsample Matrix A : " << endl;
	print3DMatrix(avgupsampleMatrix);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}