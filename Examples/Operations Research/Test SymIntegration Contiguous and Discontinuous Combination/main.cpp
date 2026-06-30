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

	int n = 5; // last index from 0
	int k = 4;
    
	// Generate combinations (assuming 0-based indexing for 0 to n)
	
	vector<vector<int>> contiguous_combinations = contiguousCombinations(n,k,0);
	cout << "All contiguous combinations of " << k << " from " << n+1 << " items:\n";
	cout <<"contiguous matrix:" << endl;
	printIntMatrix(contiguous_combinations);

	vector<vector<int>> discontiguous_combinations = discontiguousCombinations(n,k,0);
	cout << "All discontiguous combinations of " << k << " from " << n+1 << " items:\n";
	cout <<"discontiguous matrix:" << endl;
	printIntMatrix(discontiguous_combinations);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}