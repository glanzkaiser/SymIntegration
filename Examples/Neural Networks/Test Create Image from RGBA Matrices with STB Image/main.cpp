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

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


int main() {
	int width = 434;
	int height = 679;
	int channels = 4; // RGBA

	// Get starting timepoint
	auto start = high_resolution_clock::now();
	// Allocate memory for the pixel array (width * height * 4 bytes)
	uint8_t* pixels = (uint8_t*)malloc(width * height * channels);	
	//int pixels[width * height * channels];
	if (pixels == NULL) 
	{
		printf("Memory allocation failed!\n");
		return 1;
	}

	vector<vector<int>> matrixR = loadIntMatrixFromFile("matrixR.txt");
	vector<vector<int>> matrixG = loadIntMatrixFromFile("matrixG.txt");
	vector<vector<int>> matrixB = loadIntMatrixFromFile("matrixB.txt");
	vector<vector<int>> matrixA = loadIntMatrixFromFile("matrixA.txt");
	
	// Generate the image from RGBA matrices
	for (int y = 0; y < height; y++) 
	{
		for (int x = 0; x < width; x++) 
		{
			int index = (y * width + x) * channels;
			pixels[index + 0] = (uint8_t)matrixR[y][x];        // Red channel
			pixels[index + 1] = (uint8_t)matrixG[y][x];        // Green channel
			pixels[index + 2] = (uint8_t)matrixB[y][x];	// Blue channel
			pixels[index + 3] = (uint8_t)matrixA[y][x];	// A channel
		}
	}

	// Calculate stride: distance in bytes between the start of consecutive rows
	int stride = width * channels; 
	stbi_flip_vertically_on_write(true); // tell stb_image.h to flip the saved image's on the y-axis.
	// Save the pixel data to a PNG file
	if (stbi_write_png("SineBamBam.png", width, height, channels, pixels, stride)) 
	{
		printf("Image saved successfully!\n");
	} 
	else 
	{
		printf("Failed to save image.\n");
	}
	 // stbi_write_jpg("output.jpg", width, height, channels, pixels, 100); // 100 is quality (1-100)
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	// Always clean up allocated memory
	free(pixels);
	return 0;
}