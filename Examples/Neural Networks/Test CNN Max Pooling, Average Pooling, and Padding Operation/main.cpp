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
	
	vector<vector<int>> inputMatrix = loadIntMatrixFromFile("matrix.txt");

	cout << "\nMatrix A : " << endl;
	printIntMatrix(inputMatrix);
	
	int kernelSize = 2;
	int stride = 2;
	vector<vector<int>> maxMatrix = CNN_2DmaxPooling(inputMatrix, kernelSize, stride) ;
	cout << "\nMax Pooling : " << endl;
	printIntMatrix(maxMatrix);

	vector<vector<double>> avgMatrix = CNN_2DaveragePooling(inputMatrix, kernelSize, stride) ;
	cout << "\nAverage Pooling : " << endl;
	printMatrix(avgMatrix);

	vector<vector<int>> padMatrix = CNN_2DpadBorder(inputMatrix, 1);
	cout << "\nPadding Matrix A : " << endl;
	printIntMatrix(padMatrix);


	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}