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


class Matrix3DTo2D {
private:
	vector<double> data;
	size_t width, height, depth;

public:
	Matrix3DTo2D(size_t h, size_t w, size_t d) : height(h), width(w), depth(d), data(h * w * d, 0.0) {}

	// Set value using 3D indices
	void set(size_t i, size_t j, size_t k, double value) 
	{
		data[(i * width * depth) + (j * depth) + k] = value;
	}

	// Get value using 3D indices
	double get(size_t i, size_t j, size_t k) const 
	{
		return data[(i * width * depth) + (j * depth) + k];
	}

	// Sum over the depth at [i][j] into a 2D matrix
	vector<double> collapseTo2D() const 
	{
		vector<double> collapsed2D(height * width, 0.0);

		for (size_t i = 0; i < height; ++i) 
		{
			for (size_t j = 0; j < width; ++j) 
			{
				double sum = 0.0;
				for (size_t k = 0; k < depth; ++k) 
				{
					sum += get(i, j, k);
				}
				collapsed2D[(i * width) + j] = sum;
			}
	}
	return collapsed2D;
	}
};


// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	vector<double> inputVector = loadVectorFromFile("matrix.txt");

	
	vector<vector<vector<double>>> Matrix3D = Create3DMatrixfromVector(inputVector, 3, 3, 3);
	cout << "\nMatrix A : " << endl;
	print3DMatrix(Matrix3D);
	
	 // 3 rows, 3 columns, 3 depth layers
	Matrix3DTo2D mat(3, 3, 3);
    	
	for (int k= 0; k<3; ++k)
	{
		for (int i = 0; i < 3 ; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				mat.set(i, j, k, inputVector[9*k + 3*i + j]);
			}
		}
	}	

	/*
	// 2 rows, 3 columns, 2 depth layers
	Matrix3DTo2D mat(2, 3, 2);

	// Layer 0
	mat.set(0, 0, 0, 1.0); mat.set(0, 1, 0, 2.0); mat.set(0, 2, 0, 3.0);
	mat.set(1, 0, 0, 4.0); mat.set(1, 1, 0, 5.0); mat.set(1, 2, 0, 6.0);
    
	// Layer 1
	mat.set(0, 0, 1, 0.5); mat.set(0, 1, 1, 1.5); mat.set(0, 2, 1, 2.5);
	mat.set(1, 0, 1, 3.5); mat.set(1, 1, 1, 4.5); mat.set(1, 2, 1, 5.5);
	*/
	vector<double> sumMatrix = mat.collapseTo2D();

	// Display the resulting 2D values
	cout << "Collapsed 2D Matrix:\n";
	for (size_t i = 0; i < 3; ++i) 
	{
		for (size_t j = 0; j < 3; ++j) 
		{
			cout << sumMatrix[(i * 3) + j] << " ";
		}
		cout << "\n";
	}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}