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

	vector<vector<double>> kernelMatrix180 = rotate180(kernelMatrix);
	vector<vector<double>> inputMatrix180 = rotate180(inputMatrix);
	cout << "\nrot180(Kernel) : " << endl;
	printMatrix(kernelMatrix180);


	vector<vector<double>> OutputMap = CNN_2DConvolutionOperation(inputMatrix, kernelMatrix, 1);
	cout << "\nOutput Map from function : " << endl;
	printMatrix(OutputMap);

	vector<vector<double>> OutputMap2 = CNN_2DConvolutionOperation(inputMatrix, kernelMatrix180, 1);
	cout << "\nOutput Map from function rot180: " << endl;
	printMatrix(OutputMap2);

	vector<vector<double>> OutputMap3 = CNN_2DConvolutionOperation(inputMatrix180, kernelMatrix, 1);
	cout << "\nOutput Map from function rot180: " << endl;
	printMatrix(OutputMap3);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}