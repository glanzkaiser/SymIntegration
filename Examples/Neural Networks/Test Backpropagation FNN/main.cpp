// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	
	vector<vector<double>> weights = loadMatrixFromFile("weights.txt"); // [num_neurons][num_inputs]
	vector<vector<double>> weightshidden = loadMatrixFromFile("weightshidden.txt"); // [num_hidden][num_outputs]
	vector<double> bias1 = loadVectorFromFile("bias1.txt");
	vector<double> bias2 = loadVectorFromFile("bias2.txt");
	vector<double> input = loadVectorFromFile("input.txt");
	vector<double> output = loadVectorFromFile("output.txt");

	cout << "Input:" << endl;
	printVector(input);
	cout << "\nInitial Weights :" << endl;
	printMatrix(weights);
	cout << "\nInitial Hidden Layer Weights :" << endl;
	printMatrix(weightshidden);
	cout << "\nBiases:" << endl;
	printVector(bias1);
	cout << "\nActual Output:" << endl;
	printVector(output);

	cout << "\nFeed-forward Neural Network with 1 hidden layer:" << endl;
	FNN_1hiddenlayer(input,weights,weightshidden,bias1, bias2, output);

	cout << "\nFeed-forward Neural Network with 1 hidden layer type 2:" << endl;
	FNN_1hiddenlayer(input,output, 2);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}