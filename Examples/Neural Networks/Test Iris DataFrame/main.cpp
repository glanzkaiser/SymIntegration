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

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

struct IrisRecord 
{
	double sepal_length, sepal_width, petal_length, petal_width;
	std::string species;
};

std::vector<IrisRecord> loadIris(const std::string& filename) 
{
	std::vector<IrisRecord> data;
	std::ifstream file(filename);
	std::string line;

	while (std::getline(file, line)) 
	{
	if (line.empty()) continue;
	std::stringstream ss(line);
	std::string val;
	IrisRecord record;

	std::getline(ss, val, ','); record.sepal_length = std::stod(val);
	std::getline(ss, val, ','); record.sepal_width = std::stod(val);
	std::getline(ss, val, ','); record.petal_length = std::stod(val);
	std::getline(ss, val, ','); record.petal_width = std::stod(val);
	std::getline(ss, record.species, ',');

	data.push_back(record);
	}
	return data;
}

void showIrisVector(const std::vector<IrisRecord>& data) 
{
	cout << "Species \t | Sepal length  " << setw(10) << "| Sepal width  " << setw(10) << " | Petal length " << setw(10) << " | Petal width " << endl;
	for (const auto& record : data) 
	{
		cout << record.species 
		<< " \t | \t " << record.sepal_length << " \t | \t " << record.sepal_width 
		<< " \t | \t " << record.petal_length << " \t | \t" << record.petal_width 
		<< endl;
	}
}


// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	std::vector<IrisRecord> data = loadIris("iris.csv");

	int n = data.size();
	showIrisVector(data);
	cout << "\nNumber of data: " << n << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}