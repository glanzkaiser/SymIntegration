/*
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#include <iterator>
#include <vector>
#include <map>
#include <algorithm> // For std::max_element,  std::sort
#include <numeric> // For std::accumulate
#include <iostream>
#include <string>
#include <fstream> // For file operations
#include <sstream> // Required for std::ostringstream
#include <string>

#include <random> // For random number generation
#include <chrono>

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_DEFINE

void FNN_1hiddenlayer(vector<double>& input, vector<vector<double>>& weights, vector<vector<double>>& hiddenweights, vector<double>& biases, vector<double>& actual_output)
{
	int n_weights = weights.size();
	int n_hiddenweights = hiddenweights.size();
	int n_input = input.size();
	int n_output = actual_output.size();
	vector<double> hiddenlayer(weights.size(), 0.0);
        vector<double> predicted_output(actual_output.size(), 0.0);
        vector<double> error(actual_output.size(), 0.0);
        
	double learning_rate = 0.51;
	for (int i = 0; i < n_weights; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_input; ++j) 
		{
			hiddenlayer[i] += weights[i][j] * input[j];
		}
		hiddenlayer[i] += biases[i];
	}
	
	for (int i = 0; i < n_input; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-hiddenlayer[i]));
	}
	cout <<"\nHidden layer:" << endl;
	printVector(hiddenlayer);	
	for (int i = 0; i < n_hiddenweights; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_output; ++j) 
		{
			predicted_output[i] += hiddenweights[i][j] * hiddenlayer[j];
		}
		predicted_output[i] += biases[i];
	}
	for (int i = 0; i < n_output; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		predicted_output[i] = divisiond(1.0, 1.0 + std::exp(-predicted_output[i]));
		error[i] = actual_output[i] - predicted_output[i];
	}	
	cout <<"\nPredicted Output:" << endl;
	printVector(predicted_output);	
	cout <<"\nError vector:" << endl;
	printVector(error);	
	cout <<"\n||Error|| = "<< norm(error) << endl;

	vector<vector<double>> new_weights(n_weights, vector<double>(n_input));
	vector<vector<double>> new_hiddenweights(n_hiddenweights, vector<double>(n_output));
	vector<double> new_hiddenlayer(weights.size(), 0.0);
        vector<vector<double>> delta_hiddenweights(n_hiddenweights, vector<double>(n_output));
	vector<double> delta(actual_output.size(), 0.0);
	for (int i = 0; i < n_weights; ++i) 
	{
		for (int j = 0; j < n_input; ++j) 
		{
			new_weights[i][j] = weights[i][j] ;
		}
	}	
	for (int i = 0; i < n_hiddenweights; ++i) 
	{
		for (int j = 0; j < n_output; ++j) 
		{
			new_hiddenweights[i][j] = hiddenweights[i][j] ;
		}
	}	
	for (int i = 0; i < n_input; ++i) 
	{
		new_hiddenlayer[i] = hiddenlayer[i];
	}
	//printMatrix(new_weights);	
	//printMatrix(new_hiddenweights);	

	cout << "\n***************************************************************************************" << endl;
	cout <<"\n \t\t\t Backpropagation iteration" << endl;
	cout << "\n***************************************************************************************" << endl;
	cout <<"\nLearning rate: " << learning_rate << endl;
	// for loop till converge starts here
	for (int i = 0; i < 5; ++i) 
	{
		cout << "\nIteration = " << i+1 << endl;
		vector<double> new_predicted_output(actual_output.size(), 0.0);
		
		for (int i = 0; i < n_output; ++i) 
		{
			delta[i] = predicted_output[i]*(1-predicted_output[i])*(error[i]);
		}
		cout <<"\nDelta output vector:" << endl;
		printVector(delta);	
		for (int i = 0; i < n_hiddenweights; ++i) 
		{
			for (int j = 0; j < n_output; ++j) 
			{
				delta_hiddenweights[i][j] = new_hiddenlayer[i]*(1-new_hiddenlayer[i])*(new_hiddenweights[i][j]*delta[i]);
			}
		}
		cout <<"\nDelta hidden layer matrix:" << endl;
		printMatrix(delta_hiddenweights);	
		for (int i = 0; i < n_hiddenweights; ++i) 
		{
			for (int j = 0; j < n_output; ++j) 
			{
				new_hiddenweights[i][j] = new_hiddenweights[i][j] + (learning_rate*delta[i]*new_hiddenlayer[j]);
			}
		}	
		for (int i = 0; i < n_weights; ++i) 
		{
			for (int j = 0; j < n_input; ++j) 
			{
				new_weights[i][j] = new_weights[i][j] + (learning_rate*delta_hiddenweights[i][j]*input[i]);
			}
		}	
		cout <<"\nNew weights matrix:" << endl;
		printMatrix(new_weights);	
		cout <<"\nNew hiddenweights matrix:" << endl;
		printMatrix(new_hiddenweights);	

		for (int i = 0; i < n_input; ++i) 
		{
			new_hiddenlayer[i] = 0;
		}
		for (int i = 0; i < n_weights; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				new_hiddenlayer[i] += new_weights[i][j] * input[j];
			}
			new_hiddenlayer[i] += biases[i];
		}
		for (int i = 0; i < n_input; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-new_hiddenlayer[i]));
		}
		cout <<"\nHidden layer:" << endl;
		printVector(new_hiddenlayer);	
		for (int i = 0; i < n_hiddenweights; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_output; ++j) 
			{
				new_predicted_output[i] += new_hiddenweights[i][j] * new_hiddenlayer[j];
			}
			new_predicted_output[i] += biases[i];
		}
		for (int i = 0; i < n_output; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_predicted_output[i] = divisiond(1.0, 1.0 + std::exp(-new_predicted_output[i]));
			error[i] = actual_output[i] - new_predicted_output[i];
		}	
		cout <<"\nPredicted Output:" << endl;
		printVector(new_predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;
	}	
}

double random_double(double min, double max) 
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(min, max);
	return dis(gen);
}

void FNN_1hiddenlayer(vector<double>& input, vector<double>& actual_output, int n_hiddenlayer)
{
	int n_input = input.size();
	int n_output = actual_output.size();
	// n_hiddenlayer is the number of nodes in the hidden layer
	vector<vector<double>> weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> hiddenweights(n_output, vector<double>(n_hiddenlayer));
	int n_weights = weights.size();
	int n_hiddenweights = hiddenweights.size();

	for (int i = 0; i < n_hiddenlayer; ++i)
	{
		for(int j = 0; j < n_input; ++j)
		{
			weights[i][j] = random_double(0.1,1);	
		}	
	}
	for (int i = 0; i < n_output; ++i)
	{
		for(int j = 0; j < n_hiddenlayer; ++j)
		{
			hiddenweights[i][j] = random_double(0.1,1);	
		}	
	}
	cout <<"\nWeights from input to hidden layer:" << endl;
	printMatrix(weights);
	cout <<"\nWeights from hidden layer to output layer:" << endl;
	printMatrix(hiddenweights);
	vector<double> hiddenlayer(weights.size(), 0.0);
        vector<double> biases(weights.size(), 0.0); // biases set to 0
        vector<double> predicted_output(actual_output.size(), 0.0);
        vector<double> error(actual_output.size(), 0.0);
        
	double learning_rate = 0.51;
	for (int i = 0; i < n_weights; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_input; ++j) 
		{
			hiddenlayer[i] += weights[i][j] * input[j];
		}
		hiddenlayer[i] += biases[i];
	}
	for (int i = 0; i < n_input; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-hiddenlayer[i]));
	}
	cout <<"\nHidden layer:" << endl;
	printVector(hiddenlayer);	
	for (int i = 0; i < n_hiddenweights; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_output; ++j) 
		{
			predicted_output[i] += hiddenweights[i][j] * hiddenlayer[j];
		}
		predicted_output[i] += biases[i];
	}
	for (int i = 0; i < n_output; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		predicted_output[i] = divisiond(1.0, 1.0 + std::exp(-predicted_output[i]));
		error[i] = actual_output[i] - predicted_output[i];
	}	
	cout <<"\nPredicted Output:" << endl;
	printVector(predicted_output);	
	cout <<"\nError vector:" << endl;
	printVector(error);	
	cout <<"\n||Error|| = "<< norm(error) << endl;

	vector<vector<double>> new_weights(n_weights, vector<double>(n_input));
	vector<vector<double>> new_hiddenweights(n_hiddenweights, vector<double>(n_output));
	vector<double> new_hiddenlayer(weights.size(), 0.0);
        vector<vector<double>> delta_hiddenweights(n_hiddenweights, vector<double>(n_output));
	vector<double> delta(actual_output.size(), 0.0);
	for (int i = 0; i < n_weights; ++i) 
	{
		for (int j = 0; j < n_input; ++j) 
		{
			new_weights[i][j] = weights[i][j] ;
		}
	}	
	for (int i = 0; i < n_hiddenweights; ++i) 
	{
		for (int j = 0; j < n_output; ++j) 
		{
			new_hiddenweights[i][j] = hiddenweights[i][j] ;
		}
	}	
	for (int i = 0; i < n_input; ++i) 
	{
		new_hiddenlayer[i] = hiddenlayer[i];
	}
	//printMatrix(new_weights);	
	//printMatrix(new_hiddenweights);	

	cout << "\n***************************************************************************************" << endl;
	cout <<"\n \t\t\t Backpropagation iteration" << endl;
	cout << "\n***************************************************************************************" << endl;
	cout <<"\nLearning rate: " << learning_rate << endl;
	// for loop till converge starts here
	for (int i = 0; i < 5; ++i) 
	{
		cout << "\nIteration = " << i+1 << endl;
		vector<double> new_predicted_output(actual_output.size(), 0.0);
		
		for (int i = 0; i < n_output; ++i) 
		{
			delta[i] = predicted_output[i]*(1-predicted_output[i])*(error[i]);
		}
		cout <<"\nDelta output vector:" << endl;
		printVector(delta);	
		for (int i = 0; i < n_hiddenweights; ++i) 
		{
			for (int j = 0; j < n_output; ++j) 
			{
				delta_hiddenweights[i][j] = new_hiddenlayer[i]*(1-new_hiddenlayer[i])*(new_hiddenweights[i][j]*delta[i]);
			}
		}
		cout <<"\nDelta hidden layer matrix:" << endl;
		printMatrix(delta_hiddenweights);	
		for (int i = 0; i < n_hiddenweights; ++i) 
		{
			for (int j = 0; j < n_output; ++j) 
			{
				new_hiddenweights[i][j] = new_hiddenweights[i][j] + (learning_rate*delta[i]*new_hiddenlayer[j]);
			}
		}	
		for (int i = 0; i < n_weights; ++i) 
		{
			for (int j = 0; j < n_input; ++j) 
			{
				new_weights[i][j] = new_weights[i][j] + (learning_rate*delta_hiddenweights[i][j]*input[i]);
			}
		}	
		cout <<"\nNew weights matrix:" << endl;
		printMatrix(new_weights);	
		cout <<"\nNew hiddenweights matrix:" << endl;
		printMatrix(new_hiddenweights);	

		for (int i = 0; i < n_input; ++i) 
		{
			new_hiddenlayer[i] = 0;
		}
		for (int i = 0; i < n_weights; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				new_hiddenlayer[i] += new_weights[i][j] * input[j];
			}
			new_hiddenlayer[i] += biases[i];
		}
		for (int i = 0; i < n_input; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-new_hiddenlayer[i]));
		}
		cout <<"\nHidden layer:" << endl;
		printVector(new_hiddenlayer);	
		for (int i = 0; i < n_hiddenweights; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_output; ++j) 
			{
				new_predicted_output[i] += new_hiddenweights[i][j] * new_hiddenlayer[j];
			}
			new_predicted_output[i] += biases[i];
		}
		for (int i = 0; i < n_output; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_predicted_output[i] = divisiond(1.0, 1.0 + std::exp(-new_predicted_output[i]));
			error[i] = actual_output[i] - new_predicted_output[i];
		}	
		cout <<"\nPredicted Output:" << endl;
		printVector(new_predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;
	}	
}
#endif
#endif