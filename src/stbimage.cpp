/*
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#include <iterator>
#include <vector>
#include <algorithm> // For std::max_element,  std::sort
#include <iostream>
#include <string>
#include <fstream> // For file operations
#include <sstream> // Required for std::ostringstream
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include <stb-master/stb_image_resize.h> // for stbir_resize_uint8
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb_image_write.h>

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_STBIMAGE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_STBIMAGE_DEFINE

using namespace std;

void ResizeAndGrayscale(const std::string& inputPath, const std::string& outputPath, int newWidth, int newHeight) 
{
	int width, height, channels;
    
	// 1. Load image
	unsigned char* img = stbi_load(inputPath.c_str(), &width, &height, &channels, 3);
	if (!img) 
	{
		std::cerr << "Failed to load image: " << inputPath << "\n";
		return;
	}

	// 2. Allocate buffer for resized RGB image / Prepare memory for resized image (retain original channels)
	int desired_channels = (channels >= 3) ? 3 : channels; 
	vector<unsigned char> resizedImg(newWidth * newHeight * desired_channels);
	 

	// 3. Resize image using STB (Default to 3 channels)
	stbir_resize_uint8(img, width, height, 0,
				resizedImg.data(), newWidth, newHeight, 0, 
				desired_channels);

	// 4. Convert resized RGB to Grayscale
	vector<unsigned char> grayImg(newWidth * newHeight);
	for (int y = 0; y < newHeight; ++y) 
	{
		for (int x = 0; x < newWidth; ++x) 
		{
		int srcIndex = (y * newWidth + x) * desired_channels;
		unsigned char r = resizedImg[srcIndex];
		unsigned char g = resizedImg[srcIndex + 1];
		unsigned char b = resizedImg[srcIndex + 2];
		//unsigned char g = (desired_channels >= 3) ? resizedImg[srcIndex + 1] : r;
		//unsigned char b = (desired_channels >= 3) ? resizedImg[srcIndex + 2] : r;
			// Luminosity formula (standard Rec. 709 weights)
			 // Standard ITU-R BT.709 weights for human perception of brightness -> 0.2126 * r + 0.7152 * g + 0.0722 * b
			// Alternative  Grayscale = 0.299*R + 0.587*G + 0.114*B
			grayImg[y * newWidth + x] = (unsigned char)(0.2126 * r + 0.7152 * g + 0.0722 * b);
		}
	}

	// 5. Save Grayscale image (using 1 for channel count)
	//int success = stbi_write_png(outputPath.c_str(), newWidth, newHeight, desired_channels, grayImg.data(), newWidth * desired_channels);
	int success = stbi_write_jpg(outputPath.c_str(), newWidth, newHeight, 1, grayImg.data(), 100); // quality = 100
	if (!success) 
	{
		cerr << "Error: Failed to save image " << outputPath << endl;
	} 
	else 
	{
		cout << "Successfully converted " << inputPath << " to " << outputPath << endl;
	}	
	
	// Free original image memory
	stbi_image_free(img);
	//cout << "Processed: " << outputPath << "\n";
}


// Helper function to convert little-endian to big-endian (required by IDX format)
uint32_t swap_endian(uint32_t val) 
{
	return ((val << 24) & 0xFF000000) |
		((val << 8)  & 0x00FF0000) |
		((val >> 8)  & 0x0000FF00) |
		((val >> 24) & 0x000000FF);
}

void convert_to_idx(const std::vector<std::string>& image_paths, const std::vector<uint8_t>& labels, 
                    const std::string& images_out_path, const std::string& labels_out_path, 
                    int width, int height) 
{
    
	// ==========================================
	// 1. Write Image Data (idx3-ubyte)
	// ==========================================
	std::ofstream img_file(images_out_path, std::ios::binary);
	if (!img_file.is_open()) 
	{
		cerr << "Failed to open output image file.\n";
		return;
	}
	
	// Write Magic Number (2051 for idx3-ubyte)
	uint32_t magic_image = swap_endian(2051); // 2051 corresponds to magic number 0x00000803
	uint32_t num_images = swap_endian(static_cast<uint32_t>(image_paths.size()));
	uint32_t img_rows = swap_endian(static_cast<uint32_t>(height));
	uint32_t img_cols = swap_endian(static_cast<uint32_t>(width));

	img_file.write(reinterpret_cast<char*>(&magic_image), sizeof(magic_image));
	img_file.write(reinterpret_cast<char*>(&num_images), sizeof(num_images));
	img_file.write(reinterpret_cast<char*>(&img_rows), sizeof(img_rows));
	img_file.write(reinterpret_cast<char*>(&img_cols), sizeof(img_cols));

	for (const auto& path : image_paths) 
	{
		int w, h, channels;
		// Force loading as grayscale (1 channel)
		unsigned char* data = stbi_load(path.c_str(), &w, &h, &channels, 1);
		
		if (data) 
		{
			if (w != width || h != height) 
			{
				cerr << "Image dimension mismatch: " << path << ". Expected " << width << "x" << height << "\n";
			}
			img_file.write(reinterpret_cast<char*>(data), width * height);
			stbi_image_free(data);
		} 
		else 
		{
			cerr << "Failed to load image: " << path << "\n";
			stbi_image_free(data);
		}
	}
	img_file.close();

	// ==========================================
	// 2. Write Label Data (idx1-ubyte)
	// ==========================================
	std::ofstream lbl_file(labels_out_path, std::ios::binary);
	if (!lbl_file.is_open()) 
	{
		cerr << "Failed to open output label file.\n";
		return;
	}

	// Write Magic Number (2049 for idx1-ubyte)
	uint32_t magic_label = swap_endian(2049); // 2049 corresponds to magic number 0x00000801
	uint32_t num_labels = swap_endian(static_cast<uint32_t>(labels.size()));

	lbl_file.write(reinterpret_cast<char*>(&magic_label), sizeof(magic_label));
	lbl_file.write(reinterpret_cast<char*>(&num_labels), sizeof(num_labels));
    
	for (uint8_t label : labels) 
	{
		lbl_file.write(reinterpret_cast<char*>(&label), sizeof(label));
	}
	lbl_file.close();

	cout << "Conversion completed successfully!\n";
}

#endif
#endif