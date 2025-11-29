// g++ -o result main.cpp -lsymintegration 
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include "symintegrationc++.h"

#include <chrono>

#define π 3.141592653589793238462643383279502884f

using namespace std::chrono;
using namespace std;

Symbolic divisionsym(Symbolic x, Symbolic y)
{
	return x/y;
}

double divisionint(int x, int y)
{
	return x/y;
}

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
		
	dvec y = loadVectorFromFile("Vectory.txt");
	dvec z = loadVectorFromFile("Vectorz.txt");
	//cout << "\nVector y : " << endl;
	//printVector(y);

	cout << endl;
	confidenceinterval_sigmaknown(y,z,6,8,0.04);
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}