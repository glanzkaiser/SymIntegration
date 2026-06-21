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

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Helper function to convert Big-Endian to Little-Endian (for x86/ARM)
uint32_t reverseBytes(uint32_t i) 
{
		return (i >> 24) | 
		((i << 8) & 0x00FF0000) | 
		((i >> 8) & 0x0000FF00) | 
		(i << 24);
}

struct MNISTImages 
{
	uint32_t magic_number;
	uint32_t number_of_images;
	uint32_t rows;
	uint32_t columns;
	// Outer vector = list of images, Inner vector = raw pixel data (0-255)
	std::vector<std::vector<uint8_t>> data; 
};

MNISTImages parseIdx3Ubyte(const std::string& filepath) 
{
	MNISTImages dataset;

	// Open the file in binary mode
	std::ifstream file(filepath, std::ios::binary);
	if (!file.is_open()) 
	{
		throw std::runtime_error("Failed to open file: " + filepath);
	}

	// Read header information (all are 32-bit integers)
	file.read(reinterpret_cast<char*>(&dataset.magic_number), sizeof(dataset.magic_number));
	file.read(reinterpret_cast<char*>(&dataset.number_of_images), sizeof(dataset.number_of_images));
	file.read(reinterpret_cast<char*>(&dataset.rows), sizeof(dataset.rows));
	file.read(reinterpret_cast<char*>(&dataset.columns), sizeof(dataset.columns));

	// Convert from Big-Endian to the host system's endianness
	dataset.magic_number = reverseBytes(dataset.magic_number);
	dataset.number_of_images = reverseBytes(dataset.number_of_images);
	dataset.rows = reverseBytes(dataset.rows);
	dataset.columns = reverseBytes(dataset.columns);

	// Validate the magic number (2051 is the magic number for idx3)
	if (dataset.magic_number != 2051) 
	{
		throw std::runtime_error("Invalid magic number. Expected 2051 for idx3-ubyte.");
	}

	// Calculate the size of a single image in bytes
	size_t image_size = dataset.rows * dataset.columns;

	// Allocate memory and read the pixel data
	dataset.data.resize(dataset.number_of_images, std::vector<uint8_t>(image_size));
	for (uint32_t i = 0; i < dataset.number_of_images; ++i) 
	{
		file.read(reinterpret_cast<char*>(dataset.data[i].data()), image_size);
	}

	file.close();
	return dataset;
}

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	 try 
	{
		// Replace with the path to your dataset file
		MNISTImages mnistData = parseIdx3Ubyte("train-images-idx3-ubyte");
		
		cout << "Successfully loaded " << mnistData.number_of_images << " images.\n";
		cout << "Image dimensions: " << mnistData.rows << "x" << mnistData.columns << "\n";
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
