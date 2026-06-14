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
	
	vector<vector<double>> inputMatrix = loadMatrixFromFile("matrix.txt");
	vector<vector<double>> kernelMatrix = loadMatrixFromFile("kernel.txt");
	
	cout << "\nMatrix A : " << endl;
	printMatrix(inputMatrix);
	cout << "\nKernel : " << endl;
	printMatrix(kernelMatrix);
	CNN_2DConvolutionOperation(inputMatrix, kernelMatrix) ;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}