// g++ -o result main.cpp -lsymintegration
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

	vector<vector<double>> MatrixA = loadMatrixFromFile("matrixA.txt");
	vector<vector<double>> MatrixB = loadMatrixFromFile("matrixB.txt");

	cout << "A:"<<endl;
	printMatrix(MatrixA);
	cout << "\nB:"<<endl;
	printMatrix(MatrixB);

	vector<vector<double>> AB = multiply(MatrixA, MatrixB);
	cout << "\nA*B:"<<endl;
	printMatrix(AB);
	
		
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}