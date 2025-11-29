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

	dvec vecu = loadVectorFromFile("vectorx.txt");
	dvec vecv = loadVectorFromFile("vectory.txt");
	dvec vecw = loadVectorFromFile("vectorz.txt");
		
	cout << "\nVector u : " << endl;
	printVector(vecu);
	
	cout << "\nVector v : " << endl;
	printVector(vecv);
	
	cout << "\nVector w : " << endl;
	printVector(vecw);
	
	cout << "\n Scalar triple product  u . (v x w) : " << endl;
	double stp = scalartripleproduct(vecu,vecv,vecw);
	cout << stp << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}