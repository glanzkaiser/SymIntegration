// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..
// Compile with -DSTB_IMAGE_IMPLEMENTATION.

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>
#include <fstream>
#include <vector>
#include <string>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;


int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Define your dataset specifications
	int target_width = 32;
	int target_height = 32;

	/*vector<string> images = 
	{
		"dataset/img_0.png",
		"dataset/img_1.png"
	};*/
	vector<string> images = loadStringVector("images.txt"); // more efficient
    
	// These must map 1:1 with the images vector above
	//vector<int> labels = { 5, 0 }; 
	vector<int> labels = loadIntVectorFromFile("labels.txt");
	convert_to_idx(images, labels, "train-images-idx3-ubyte", "train-labels-idx1-ubyte", target_width, target_height);


	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}
