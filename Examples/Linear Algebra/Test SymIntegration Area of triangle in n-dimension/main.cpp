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
	
	dvec a = loadVectorFromFile("a.txt");
	dvec b = loadVectorFromFile("b.txt");
	dvec c = loadVectorFromFile("c.txt");
		
	cout << "\nVertices, a : " << endl;
	printVector(a);
	cout << "\nb : " << endl;
	printVector(b);
	cout << "\nc : " << endl;
	printVector(c);
	
	double A = areaoftriangle(a,b,c);	
	cout << "\nArea of triangle : " << A << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}