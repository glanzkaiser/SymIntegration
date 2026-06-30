// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..
// Compile with -DSTB_IMAGE_IMPLEMENTATION.

// Makes life easier: bash createimageslist
#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// List of images to process
	//vector<string> files = {"img_0.png", "img_1.png","img_2.png", "img_3.png","img_4.png", "img_5.png","img_6.png", "img_7.png","img_8.png"};
	svec files = loadStringVector("images.txt"); // more efficient
	//printStringVector(files);
	int targetWidth = 32;
	int targetHeight = 32;

	// the for loop process to resize the image and convert it to grayscale one by one
	for (const auto& file : files) 
	{
		std::string outPath = "gray_" + file;
		ResizeAndGrayscale(file, outPath, targetWidth, targetHeight);
	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}
