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
		
	dvec x = loadVectorFromFile("Vectorx.txt");
	dvec y = loadVectorFromFile("Vectory.txt");
	//cout << "\nVector y : " << endl;
	//printVector(y);

	descriptivestatistics(x);
	descriptivestatistics(y);

	cout << endl;
	hypothesistest_equalunknownvariances_righttailed(x,y,2,0.05,1);
	hypothesistest_equalunknownvariances_lefttailed(x,y,2,0.05,1);
	hypothesistest_equalunknownvariances_twotailed(x,y,2,0.05,1);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}