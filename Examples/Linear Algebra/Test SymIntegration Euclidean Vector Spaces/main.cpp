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

	dvec vecx = loadVectorFromFile("vectorx.txt");
	dvec vecy = loadVectorFromFile("vectory.txt");
		
	cout << "\nVector x : " << endl;
	printVector(vecx);
	
	cout << "\nVector y : " << endl;
	printVector(vecy);
	
	cout << "\nnorm(x) : " << norm(vecx) << endl;
	cout << "\nnorm(y) : " << norm(vecy) << endl;
	
	cout << "\ndot(x,y) : " << dot(vecx,vecy) << endl;
	
	cout << "\nangle(x,y) in degree: " << radtodeg(angle(vecx,vecy)) << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}