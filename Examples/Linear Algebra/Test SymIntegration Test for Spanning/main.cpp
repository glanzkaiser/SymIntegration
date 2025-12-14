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
	
	dvec veca = loadVectorFromFile("vectora.txt");
	dvec vecb = loadVectorFromFile("vectorb.txt");
	dvec vecc = loadVectorFromFile("vectorc.txt");
	
	/*dmat A (3, vector<double>(3));
	A = addColumn(A,veca,0);
	A = addColumn(A,vecb,1);
	A = addColumn(A,vecc,2);
	deleteColumn(A,3);
	deleteColumn(A,3);
	deleteColumn(A,3);
	deleteRow(A,1);
	cout << "Matrix:" << endl;
	printMatrix(A);*/
	
	dmat columnVectors = {veca, vecb, vecc};
	dmat newmatrix = createMatrixFromColumnVectors(columnVectors);

	spanningtest(newmatrix);

	dmat mat = loadMatrixFromFile("Mat.txt");
	spanningtest(mat);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}