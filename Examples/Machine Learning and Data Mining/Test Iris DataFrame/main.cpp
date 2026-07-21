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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	vector<vector<double>> IrisMatrix = loadMatrixFromFile("iris.csv");
	
	cout << "6 Points from Iris dataset:" << endl;
	printMatrix(IrisMatrix);

	vector<vector<double>> LKM = LinearKernelMatrix(IrisMatrix);
	cout << "\nLinear Kernel Matrix:" << endl;
	printMatrix(LKM);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}