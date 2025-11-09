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
	Symbolic x("x"), y("y"), z("z"), f;
	f = 2*x-3*y+6*z + 1;
	dvec vecp0 = loadVectorFromFile("vectorx.txt");
		
	cout << "\nVector P0 : " << endl;
	printVector(vecp0);
	
	double dist = distance(f,x,y,z,vecp0);
	cout << "\nDistance between the plane 2x - 3y + 6z +1 = 0:\n " << dist << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}