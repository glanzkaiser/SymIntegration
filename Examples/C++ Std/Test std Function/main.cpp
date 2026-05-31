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

void calculate_offset(const double* input, std::vector<double>& output) 
{
    // Read from input, write to output
    for(size_t i = 0; i < output.size(); ++i) {
        output[i] = input[i] + 1.0; 
	cout << input[i] << endl;
    }
}

// Driver code
int main(int argc, char** argv)
{
		// Get starting timepoint
	auto start = high_resolution_clock::now();

	double data[] = {1.0, 2.0, 3.0};
	vector<double> results(3);

	// Assigning the function to the std::function wrapper / Wrapping the function
	std::function<void(const double*, std::vector<double>&)> func = calculate_offset;

	func(data, results);
	printVector(results);

	double data2[] = {5.0, 6.0, 7.0};
	func(data2,results);
	printVector(results);

	// Define the std::function variable
	std::function<vector<double>(vector<double>)> doubleElements;
	// Assign a lambda to it
	doubleElements = [](vector<double> input) 
	{
		for (double& val : input) 
		{
			val *= 2.0;
		}
		return input;
	};

	// Use the function
	vector<double> vec_data = {1.0, 2.5, 3.0};
	vector<double> result = doubleElements(vec_data);
	printVector(result);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}