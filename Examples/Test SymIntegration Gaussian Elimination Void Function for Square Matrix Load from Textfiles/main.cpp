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

	vector<vector<double>> doubleMatrix2 = loadMatrixFromFile("matrix.txt");
	vector<vector<double>> vector1 = loadMatrixFromFile("vectorb.txt");
	vector<double> sol;

	solve_nhsystem(doubleMatrix2,vector1,sol);

	//cout << sol[0] << endl; // To get the first solution (x1).
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}