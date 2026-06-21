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

	dmat MatrixXorig = loadMatrixFromFile("MatrixX.txt");
	dvec y= loadVectorFromFile("vectory.txt");
	
	dvec ones = vector<double>(6,1.0) ;
	dmat X = addColumn(MatrixXorig,ones,0);
	dmat X_transpose = transpose(X);
	cout << "\nMatrix X:" << endl;
	printMatrix(X);
	cout << "\nMatrix X^{T}:" << endl;
	printMatrix(X_transpose);
	
	cout <<"\nVector y:" << endl;
	printVector(y);

	// (X^{T} *X)^{-1} * X^{T} * y
	dmat op1 = multiply(X_transpose,X);
	dmat op2 = inverse(op1);
	dmat op3 = multiply(op2,X_transpose);
	dvec B = multiplymatrixvector(op3,y);

	cout << "\nB : " << endl;
	printVector(B);

	dvec xnew = loadVectorFromFile("Xnew.txt");
	cout << "\nXnew : " << endl;
	printVector(xnew);
	
	cout <<"\nXnew * B = " << dot(xnew,B) << endl;
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}