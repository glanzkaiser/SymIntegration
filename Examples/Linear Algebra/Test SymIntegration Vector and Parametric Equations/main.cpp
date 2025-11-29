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
	f = x-y+2*z - 5;
		
	SymbolicMatrix ve = vectorequation(f,x,y,z);
	cout << "\nVector equation of the plane x - y + 2z - 5 = 0:\n" << ve << endl;

	vector<double> v1(3);
	vector<double> v2(3);
	vector<double> v3(3);
	vectorequationdecomp(f,x,y,z,v1,v2,v3);
	cout << "\nThe parametric equations in vector terms: \n(x,y,z) = v1 + v2*t + v2*s" << endl;
	cout << "v1 = " << endl;
	printVector(v1);
	cout << "\nv2 = " << endl;
	printVector(v2);
	cout << "\nv3 =" << endl;
	printVector(v3);
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}