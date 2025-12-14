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

	
	dvec veca = loadVectorFromFile("veca.txt");
	dvec vecb = loadVectorFromFile("vecb.txt");
	dvec b = loadVectorFromFile("b.txt");
	
	dmat columnVectors = {veca, vecb, b};
	dmat newmatrix = createMatrixFromColumnVectors(columnVectors);

	printMatrix(newmatrix);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}