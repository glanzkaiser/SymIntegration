// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..
// Compile with -DSTB_IMAGE_IMPLEMENTATION.

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>
#include <cstdint>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	 try {
        // Change to your actual file path
	vector<vector<vector<int>>> images = loadIDX3("train-images-idx3-ubyte");	
	vector<int> labels = loadIDX1("train-labels-idx1-ubyte");
        
	cout << "Successfully loaded " << images.size() << " images!" << endl;
	cout << "Dimension of Image 0: " << images[0].size() << "x" << images[0][0].size() << endl;

	save3DMatrixint(images,"images.txt");
	saveVectorint(labels,"labels.txt");

	}
	catch (const std::exception& e) 
	{
		std::cerr << "Error: " << e.what() << "\n";
	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}
