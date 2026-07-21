/*
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#include <iterator>
#include <vector>
#include <map>
#include <algorithm> // For std::max_element,  std::sort
#include <numeric> // For std::accumulate
#include <iostream>
#include <string>
#include <fstream> // For file operations
#include <sstream> // Required for std::ostringstream
#include <string>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <climits> // for INT_MIN

#include <random> // For random number generation
#include <chrono>
#include <unordered_map>
#include <iomanip>
#include <numeric>
#include <thread>

// ANSI color codes
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define BOLD    "\033[1m"

// Debug macro
#define DEBUG(x) std::cout << YELLOW << "DEBUG: " << x << RESET << std::endl

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_MACHINELEARNING_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_MACHINELEARNING_DEFINE

using namespace std;


class Matrix2DToVector {
private:
	int height, width; // row, col
	vector<double> data;
public:
	Matrix2DToVector(int h, int w) : height(h), width(w), data(h * w , 0.0) {}

	// Set value using 2D indices
	void set(int i, int j, double value) 
	{
		data[(i * width) + j  ] = value;
	}

	// Get value using 2D indices
	double get(int i, int j) const 
	{
		return data[(i * width) + j ];
	}

};

vector<vector<double>> LinearKernelMatrix(vector<vector<double>> &input)
{
	int n = input.size();
	vector<vector<double>> result(n,vector<double>(n,0.0));

	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			double dotsum = 0.0;
			for (int k = 0; k < n ; ++k)
			{
				dotsum += input[i][k]*input[j][k];
			}
			
			result[i][j] = dotsum;
		}
	}
	return result;
}

#endif
#endif