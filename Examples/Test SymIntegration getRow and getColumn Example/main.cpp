// g++ -o result main.cpp -lsymbolicc++
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

	vector<vector<double>> MatrixA = loadMatrixFromFile("matrix.txt");

	int R = MatrixA.size();
	int C = MatrixA[0].size();
	 			
	cout << "Column-1" << endl;
	for (int i = 0; i<R; i++)
	{	
		// take the first column
		cout << getColumn(MatrixA,0)[i] << endl; 
	}

	cout << "Row-2" << endl;
	for (int i = 0; i<C; i++)
	{
		// take the second row
		cout << getRow(MatrixA,1)[i] << endl;
	}

	cout << "vector row 1:"<<endl;
	vector<double> vectora = getRow(MatrixA,0) ;
	printVector(vectora);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}