// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

double computeCCELoss(const vector<double>& y_true, const vector<double>& y_pred) 
{
	double loss = 0.0;
	double epsilon = 1e-15; // To prevent log(0) which is undefined
	int n = y_true.size();
	for (int i = 0; i < n; ++i) 
	{
		// Clip predictions to avoid numerical instability
		double clipped_pred = std::max(epsilon, std::min(1.0 - epsilon, y_pred[i]));
		loss += y_true[i] * std::log(clipped_pred);
	}
	return -loss;
}
// Driver code
int main(int argc, char** argv)
{
		// Get starting timepoint
	auto start = high_resolution_clock::now();

	vector<double> y_true = {0, 1, 0};      // Target: Class 1
	vector<double> y_pred = {0.1, 0.8, 0.1}; // Model predictions
    
	cout << "CCE Loss: " << computeCCELoss(y_true, y_pred) << endl;
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}