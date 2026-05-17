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

// Element-wise ReLU for a vector (layer activation)
void apply_relu(std::vector<double>& layer_output) 
{
	for (double& val : layer_output) 
	{
		val = std::max(0.0, val);
	}
}
//Basic FNN structures

void ANN_ReLU_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	
	std::vector<double> output(weights.size(), 0.0);
        
	for (int i = 0; i < weights.size(); ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < input.size(); ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
		// ReLU Activation: max(0, x)
		output[i] = std::max(0.0, output[i]);
	}
	
	printVector(output);
};


void ANN_Sigmoid_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	
	std::vector<double> output(weights.size(), 0.0);
        
	for (int i = 0; i < weights.size(); ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < input.size(); ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
		// Sigmoid Activation:  1 / (1 + exp(-x))
		output[i] = divisiond(1.0, 1.0 + std::exp(-output[i]));
	}
	
	printVector(output);
};


void ANN_Tanh_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	
	std::vector<double> output(weights.size(), 0.0);
        
	for (int i = 0; i < weights.size(); ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < input.size(); ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
		// Tanh Activation:  Tanh(x)
		output[i] = std::tanh(output[i]);
	}
	
	printVector(output);
};

void ANN_SoftMax_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	std::vector<double> output(weights.size(), 0.0);
        
	for (int i = 0; i < weights.size(); ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < input.size(); ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
	}
	// SoftMax Activation:  
	// 1. Find max for numerical stability
	double maxLogit = *std::max_element(output.begin(), output.end());
	//cout << "max element = " << maxLogit << endl;

	// 2. Compute exp(x_i - max) and sum
	for (int i = 0; i < input.size() ; ++i) 
	{
		output[i] = std::exp(output[i] - maxLogit);
	}
	double sumExp = accumulate(output.begin(),output.end(), 0.0);
	//cout << "sum Exp element = " << sumExp << endl;

	// 3. Normalize
	for (int i = 0; i < input.size() ; ++i) 
	{
		output[i] = divisiond(output[i],sumExp);
	}

	printVector(output);
};
// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	vector<vector<double>> weights = loadMatrixFromFile("weights.txt"); // [num_neurons][num_inputs]
	vector<double> biases = loadVectorFromFile("biases.txt");
	vector<double> input = loadVectorFromFile("input.txt");

	cout << "Input:" << endl;
	printVector(input);
	cout << "\nWeights :" << endl;
	printMatrix(weights);
	cout << "\nBiases:" << endl;
	printVector(biases);

	cout << "\nReLU Activation function:" << endl;
	ANN_ReLU_activationfunction(input,weights,biases);
	cout << "\nSigmoid Activation function:" << endl;
	ANN_Sigmoid_activationfunction(input,weights,biases);
	cout << "\nTanh Activation function:" << endl;
	ANN_Tanh_activationfunction(input,weights,biases);
	cout << "\nSoftMax Activation function:" << endl;
	ANN_SoftMax_activationfunction(input,weights,biases);

	//cout << "\nSoftMax Activation function:" << endl;
	//vector<double> sout = SoftMax_vectorresult_activationfunction(input);
	//printVector(sout);
	//cout << "\nSoftMax Activation function:" << endl;
	//SoftMax_logitsinput_activationfunction(input);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}