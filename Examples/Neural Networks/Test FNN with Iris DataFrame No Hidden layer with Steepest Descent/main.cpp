// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

//If you need to isolate the data into its own vector (e.g., getting all ages as a std::vector<int>), you can use std::transform from the <algorithm> library. 

#include <algorithm>
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
	cout << setw(15) << "Sepal Length" 
         << setw(15) << "Sepal Width" 
         << setw(15) << "Petal Length" 
         << setw(15) << "Petal Width" 
         << setw(15) << "Species" << endl;
	cout << string(75, '-') << endl; 
	
	for (const auto& record : data) 
	{
		cout << setw(15) << record.sepal_length 
                 << setw(15) << record.sepal_width 
                 << setw(15) << record.petal_length 
                 << setw(15) << record.petal_width 
                 << setw(15) << record.species << endl;
	}
}

void getColumnVector_IrisData(const std::vector<IrisRecord>& data, const string &inputString) 
{
	cout << setw(15) << "Sepal Length" 
         << setw(15) << "Sepal Width" 
         << setw(15) << "Petal Length" 
         << setw(15) << "Petal Width" 
         << setw(15) << "Species" << endl;
	cout << string(75, '-') << endl; 
	for (const auto& record : data) 
	{
		if(record.species ==inputString)
		{
			cout << setw(15) << record.sepal_length 
		         << setw(15) << record.sepal_width 
		         << setw(15) << record.petal_length 
		         << setw(15) << record.petal_width 
		         << setw(15) << record.species << endl;
		}
	}
}

//Memory Structure: A vector<struct> stores data in "Row-Major" order (all members of one struct are contiguous in memory). Extracting a column involves jumping across these row boundaries.

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	std::vector<IrisRecord> data_loaded = loadIris("iris.csv");

	int n = data_loaded.size();
	showIrisVector(data_loaded);
	cout << "\nNumber of data: " << n << endl;

	// to group by data by species name
	//getColumnVector_IrisData(data_loaded, "\"Setosa\"");

	// to take the column of sepal length, sepal_width, petal_length, petal_width
	vector<double> vec_x1;
	vector<double> vec_x2;
	vector<double> vec_x3;
	vector<double> vec_x4;
	std::transform(data_loaded.begin(), data_loaded.end(), std::back_inserter(vec_x1),
		[](const IrisRecord& p) { return p.sepal_length; });
	std::transform(data_loaded.begin(), data_loaded.end(), std::back_inserter(vec_x2),
		[](const IrisRecord& p) { return p.sepal_width; });
	std::transform(data_loaded.begin(), data_loaded.end(), std::back_inserter(vec_x3),
		[](const IrisRecord& p) { return p.petal_length; });
	std::transform(data_loaded.begin(), data_loaded.end(), std::back_inserter(vec_x4),
		[](const IrisRecord& p) { return p.petal_width; });
	
	vector<vector<double>> irisMatrix(n, vector<double>(4));
	
	irisMatrix = addColumn(irisMatrix,vec_x4,0);
	irisMatrix = addColumn(irisMatrix,vec_x3,0);
	irisMatrix = addColumn(irisMatrix,vec_x2,0);
	irisMatrix = addColumn(irisMatrix,vec_x1,0);
	for(int i = 0; i < 4; ++i)
	{
		irisMatrix = deleteColumn(irisMatrix,4);	
	}
	//printMatrix(newMatrix);
	vector<string> vec_x5;
	std::transform(data_loaded.begin(), data_loaded.end(), std::back_inserter(vec_x5),
		[](const IrisRecord& p) { return p.species; });
	//saveVectorstring(vec_x5,"output.txt");
	
	// We shuffle the order of the data for the training and testing
	std::vector<double> indices(n);
	std::iota(indices.begin(), indices.end(), 0);
	std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

	//cout << "\nIndex after shuffle:" << endl;
	//printVector(indices);

	//  The dataset is split into 80% training and 20% testing sets
	size_t train_size = n * 0.8;
	vector<vector<double>> X_train, X_test;
	vector<string> y_train, y_test;

	for (size_t i = 0; i < n; ++i) 
	{
		if (i < train_size) 
		{
			X_train.push_back(irisMatrix[indices[i]]);
			y_train.push_back(vec_x5[indices[i]]);
		} 
		else 
		{
			X_test.push_back(irisMatrix[indices[i]]);
			y_test.push_back(vec_x5[indices[i]]);
		}
	}
	//printStringVector(y_test);
	
	//cout << "\nTraining Feed-forward Neural Network with 0 hidden layer for Iris dataset:" << endl;
	//FNN_nohiddenlayer_irisdataset_training(X_train,y_train);
	
	cout << "\nTesting Feed-forward Neural Network with 0 hidden layer for Iris dataset:" << endl;
	//FNN_nohiddenlayer_irisdataset_testing(X_test,y_test);
	FNN_nohiddenlayer_irisdataset_testing(irisMatrix,vec_x5);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}