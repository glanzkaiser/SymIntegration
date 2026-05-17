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

double LeakyReLUActivationfunction(double x, double alpha)
{
	return std::max(alpha, x);
}
double ReLUActivationfunction(double x)
{
	return std::max(0.0, x);
}

double SigmoidActivationfunction(double x)
{
	return divisiond(1.0, 1.0 + std::exp(-x));
}
double TanhActivationfunction(double x)
{
	return std::tanh(x);
}

double LeakyReLUDerivative(double x, double alpha) 
{
	return (x > 0) ? 1.0 : alpha;
}
double ReLUDerivative(double x) 
{
	return (x > 0) ? 1.0 : 0.0;
}
double SigmoidDerivative(double x) 
{
	return x*(1-x);
}
double TanhDerivative(double x) 
{
	double t = std::tanh(x);
	return 1.0 - t * t;
}

vector<double> SoftMax_vectorresult_activationfunction(vector<double>& input)
{
	int n = input.size();
	vector<double> output(n,0.0);// SoftMax Activation:  
	// 1. Find max for numerical stability
	auto maxLogit0 = std::max_element(input.begin(), input.end());
	double maxLogit = *maxLogit0;
	//cout << "max element = " << maxLogit << endl;
	vector<double> result(n,0.0);
	// 2. Compute exp(x_i - max) and sum
	for (int i = 0; i < n ; ++i) 
	{
		input[i] = std::exp(input[i] - maxLogit);
	}
	double sumExp = accumulate(input.begin(),input.end(), 0.0);
	//cout << "sum Exp element = " << sumExp << endl;

	// 3. Normalize
	for (int i = 0; i < n ; ++i)  
	{
		output[i] = divisiond(input[i],sumExp);
	}

	return output;
}

void SoftMax_logitsinput_activationfunction(vector<double>& input)
{
	int n = input.size();
	vector<double> output(n,0.0);
	// SoftMax Activation:  
	// 1. Find max for numerical stability
	auto maxLogit0 = std::max_element(input.begin(), input.end());
	double maxLogit = *maxLogit0;
	//cout << "max element = " << maxLogit << endl;
	// 2. Compute exp(x_i - max) and sum
	for (int i = 0; i < n ; ++i) 
	{
		input[i] = std::exp(input[i] - maxLogit);
	}
	double sumExp = accumulate(input.begin(),input.end(), 0.0);
	//cout << "sum Exp element = " << sumExp << endl;

	// 3. Normalize
	for (int i = 0; i < n ; ++i) 
	{
		output[i] = divisiond(input[i],sumExp);
	}

	printVector(output);
}

void ReLU_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	
	std::vector<double> output(weights.size(), 0.0);
        int w = weights.size();
	int n = input.size();
	for (int i = 0; i < w; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n; ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
		// ReLU Activation: max(0, x)
		output[i] = std::max(0.0, output[i]);
	}
	
	printVector(output);
}


void Sigmoid_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	
	std::vector<double> output(weights.size(), 0.0);
        int w = weights.size();
	int n = input.size();
	for (int i = 0; i < w; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n; ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
		// Sigmoid Activation:  1 / (1 + exp(-x))
		output[i] = divisiond(1.0, 1.0 + std::exp(-output[i]));
	}
	
	printVector(output);
}


void Tanh_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	
	std::vector<double> output(weights.size(), 0.0);
	int w = weights.size();
        int n = input.size();
	for (int i = 0; i < w; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n; ++j) 
		{
			output[i] += weights[i][j] * input[j];
		}
		output[i] += biases[i];
		// Tanh Activation:  Tanh(x)
		output[i] = std::tanh(output[i]);
	}
	
	printVector(output);
}

void SoftMax_activationfunction(vector<double>& input, vector<vector<double>>& weights, vector<double>& biases)
{
	std::vector<double> output(weights.size(), 0.0);
	int w = weights.size();
	int n = input.size();
	for (int i = 0; i < w; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n; ++j) 
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
	for (int i = 0; i < n ; ++i) 
	{
		output[i] = std::exp(output[i] - maxLogit);
	}
	double sumExp = accumulate(output.begin(),output.end(), 0.0);
	//cout << "sum Exp element = " << sumExp << endl;

	// 3. Normalize
	for (int i = 0; i < n ; ++i) 
	{
		output[i] = divisiond(output[i],sumExp);
	}

	printVector(output);
}


void FNN_1hiddenlayer(vector<double>& input, vector<vector<double>>& weights, vector<vector<double>>& hiddenweights, vector<double>& bias1, vector<double>& bias2, vector<double>& actual_output)
{
// Fixed on May 6th, 2026

	int n_hiddenlayer = weights.size();
	int n_input = input.size();
	int n_output = actual_output.size();
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
        vector<double> new_bias1(n_hiddenlayer, 0.0);
        vector<double> new_bias2(n_output, 0.0);
        vector<double> predicted_output(n_output, 0.0);
        vector<double> error(n_output, 0.0);
        
	double learning_rate = 1;
	for (int i = 0; i < n_hiddenlayer; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_input; ++j) 
		{
			hiddenlayer[i] += weights[i][j] * input[j];
		}
		hiddenlayer[i] += bias1[i];
	}
	
	for (int i = 0; i < n_input; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-hiddenlayer[i]));
	}
	cout <<"\nHidden layer:" << endl;
	printVector(hiddenlayer);	
	for (int i = 0; i < n_output; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			predicted_output[i] += hiddenweights[i][j] * hiddenlayer[j];
		}
		predicted_output[i] += bias2[i];
	}
	for (int i = 0; i < n_output; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		predicted_output[i] = SigmoidActivationfunction(predicted_output[i]);
		error[i] = actual_output[i] - predicted_output[i];
	}	
	cout <<"\nPredicted Output:" << endl;
	printVector(predicted_output);	
	cout <<"\nError vector:" << endl;
	printVector(error);	
	cout <<"\n||Error|| = "<< norm(error) << endl;

	vector<vector<double>> new_weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> new_hiddenweights(n_output, vector<double>(n_hiddenlayer));
	vector<double> new_hiddenlayer(n_hiddenlayer, 0.0);
        vector<double> delta_hiddenweights(n_hiddenlayer, 0.0);
	vector<double> delta(n_output, 0.0);
	for (int i = 0; i < n_hiddenlayer; ++i) 
	{
		new_bias1[i] = bias1[i];
		for (int j = 0; j < n_input; ++j) 
		{
			new_weights[i][j] = weights[i][j] ;
		}
	}	
	for (int i = 0; i < n_output; ++i) 
	{
		new_bias2[i] = bias2[i];
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			new_hiddenweights[i][j] = hiddenweights[i][j] ;
		}
	}	
	for (int i = 0; i < n_input; ++i) 
	{
		new_hiddenlayer[i] = hiddenlayer[i];
	}
	vector<double> new_predicted_output(actual_output.size(), 0.0);
	for (int i = 0; i < n_output; ++i) 
	{
		new_predicted_output[i] = predicted_output[i];
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
		
		for (int i = 0; i < n_output; ++i) 
		{
			delta[i] = SigmoidDerivative(new_predicted_output[i])*(error[i]);
		}
		cout <<"\nDelta output vector:" << endl;
		printVector(delta);	
			
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			double sum_output = 0;
			for (int j = 0; j < n_output; ++j) 
			{
				sum_output += delta[j]*new_hiddenweights[j][i];
			}
			delta_hiddenweights[i] = SigmoidDerivative(new_hiddenlayer[i])*(sum_output);
		}
		cout <<"\nDelta hidden layer vector:" << endl;
		printVector(delta_hiddenweights);	
		for (int i = 0; i < n_output; ++i) 
		{
			new_bias2[i] = new_bias2[i] + (learning_rate*delta[i]);	
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				new_hiddenweights[i][j] = new_hiddenweights[i][j] + (learning_rate*delta[i]*new_hiddenlayer[j]);
			}
		}	
			
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			new_bias1[i] = new_bias1[i] + (learning_rate*delta_hiddenweights[i]);
			for (int j = 0; j < n_input; ++j) 
			{
				new_weights[i][j] = new_weights[i][j] + (learning_rate*delta_hiddenweights[i]*input[j]);
			}
		}	
		cout <<"\nNew weights matrix:" << endl;
		printMatrix(new_weights);	
		cout <<"\nNew hiddenweights matrix:" << endl;
		printMatrix(new_hiddenweights);	

		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			new_hiddenlayer[i] = 0;
		}
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				new_hiddenlayer[i] += new_weights[i][j] * input[j];
			}
			new_hiddenlayer[i] += new_bias1[i];
		}
		for (int i = 0; i < n_input; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_hiddenlayer[i] = SigmoidActivationfunction(new_hiddenlayer[i]);
		}
		cout <<"\nHidden layer:" << endl;
		printVector(new_hiddenlayer);	
		for (int i = 0; i < n_output; ++i) 
		{
			new_predicted_output[i] = 0;
		}
		for (int i = 0; i < n_output; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				new_predicted_output[i] += new_hiddenweights[i][j] * new_hiddenlayer[j];
			}
			new_predicted_output[i] += new_bias2[i];
		}
		for (int i = 0; i < n_output; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_predicted_output[i] = SigmoidActivationfunction(new_predicted_output[i]);
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
// Fixed on May 6th, 2026
	int n_input = input.size();
	int n_output = actual_output.size();
	// n_hiddenlayer is the number of nodes in the hidden layer
	vector<vector<double>> weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> hiddenweights(n_output, vector<double>(n_hiddenlayer));

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
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
        vector<double> bias1(n_hiddenlayer, 0.0); // biases set to 0
        vector<double> bias2(n_output, 0.0); // biases set to 0
        vector<double> predicted_output(n_output, 0.0);
        vector<double> error(n_output, 0.0);
        
	double learning_rate = 0.51;
	for (int i = 0; i < n_hiddenlayer; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_input; ++j) 
		{
			hiddenlayer[i] += weights[i][j] * input[j];
		}
		hiddenlayer[i] += bias1[i];
	}
	for (int i = 0; i < n_input; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-hiddenlayer[i]));
	}
	cout <<"\nHidden layer:" << endl;
	printVector(hiddenlayer);	
	for (int i = 0; i < n_output; ++i) 
	{
		// Linear transformation: Wx + b
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			predicted_output[i] += hiddenweights[i][j] * hiddenlayer[j];
		}
		predicted_output[i] += bias2[i];
	}
	for (int i = 0; i < n_output; ++i) 
	{
		// Sigmoid Activation:  1 / (1 + exp(-x))
		predicted_output[i] = SigmoidActivationfunction(predicted_output[i]);
		error[i] = actual_output[i] - predicted_output[i];
	}	
	cout <<"\nPredicted Output:" << endl;
	printVector(predicted_output);	
	cout <<"\nError vector:" << endl;
	printVector(error);	
	cout <<"\n||Error|| = "<< norm(error) << endl;

	vector<vector<double>> new_weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> new_hiddenweights(n_output, vector<double>(n_hiddenlayer));
	vector<double> new_hiddenlayer(n_hiddenlayer, 0.0);
        vector<double> delta_hiddenweights(n_hiddenlayer, 0.0);
	vector<double> delta(n_output, 0.0);
	for (int i = 0; i < n_hiddenlayer; ++i) 
	{
		for (int j = 0; j < n_input; ++j) 
		{
			new_weights[i][j] = weights[i][j] ;
		}
	}	
	for (int i = 0; i < n_output; ++i) 
	{
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			new_hiddenweights[i][j] = hiddenweights[i][j] ;
		}
	}	
	for (int i = 0; i < n_hiddenlayer; ++i) 
	{
		new_hiddenlayer[i] = hiddenlayer[i];
	}
	vector<double> new_predicted_output(n_output, 0.0);	
	for (int i = 0; i < n_output; ++i) 
	{
		new_predicted_output[i] = predicted_output[i];
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
		
		for (int i = 0; i < n_output; ++i) 
		{
			delta[i] = SigmoidDerivative(new_predicted_output[i])*(error[i]);
		}
		cout <<"\nDelta output vector:" << endl;
		printVector(delta);	
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			double sum_output = 0;
			for (int j = 0; j < n_output; ++j) 
			{
				sum_output += delta[j]*new_hiddenweights[j][i];
			}
			delta_hiddenweights[i] = SigmoidDerivative(new_hiddenlayer[i])*(sum_output);
		}
		cout <<"\nDelta hidden layer vector:" << endl;
		printVector(delta_hiddenweights);	

		for (int i = 0; i < n_output; ++i) 
		{
			bias2[i] = bias2[i] + (learning_rate*delta[i]);
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				new_hiddenweights[i][j] = new_hiddenweights[i][j] + (learning_rate*delta[i]*new_hiddenlayer[j]);
			}
		}	
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			bias1[i] = bias1[i] + (learning_rate*delta_hiddenweights[i]);
			for (int j = 0; j < n_input; ++j) 
			{
				new_weights[i][j] = new_weights[i][j] + (learning_rate*delta_hiddenweights[i]*input[j]);
			}
		}	
		cout <<"\nNew weights matrix:" << endl;
		printMatrix(new_weights);	
		cout <<"\nNew hiddenweights matrix:" << endl;
		printMatrix(new_hiddenweights);	

		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			new_hiddenlayer[i] = 0;
		}
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				new_hiddenlayer[i] += new_weights[i][j] * input[j];
			}
			new_hiddenlayer[i] += bias1[i];
		}
		for (int i = 0; i < n_input; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_hiddenlayer[i] = divisiond(1.0, 1.0 + std::exp(-new_hiddenlayer[i]));
		}
		cout <<"\nHidden layer:" << endl;
		printVector(new_hiddenlayer);	
		for (int i = 0; i < n_output; ++i) 
		{
			new_predicted_output[i] = 0;
		}
		for (int i = 0; i < n_output; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				new_predicted_output[i] += new_hiddenweights[i][j] * new_hiddenlayer[j];
			}
			new_predicted_output[i] += bias2[i];
		}
		for (int i = 0; i < n_output; ++i) 
		{
			// Sigmoid Activation:  1 / (1 + exp(-x))
			new_predicted_output[i] = SigmoidActivationfunction(new_predicted_output[i]);
			error[i] = actual_output[i] - new_predicted_output[i];
		}	
		cout <<"\nPredicted Output:" << endl;
		printVector(new_predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;
	}	
}

void FNN_nohiddenlayer_irisdataset_training(vector<vector<double>>& input_matrix, vector<string>& vec_output )
{
// May 11th, 2026

	int n_classification = 3;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();
	// n_hiddenlayer is the number of nodes in the hidden layer
	vector<vector<double>> weights(n_input, vector<double>(n_classification));
        vector<double> bias_output(n_classification, 0.0); // biases set to 0
	vector<double> predicted_output(n_classification, 0.0);
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	vector<double> delta(n_classification, 0.0);

	double learning_rate = 0.1;
	int epochs = 1000;

	for (int i = 0; i < n_input; ++i)
	{
		bias_output[i] = random_double(-1,1);
		for(int j = 0; j < n_classification; ++j)
		{
			weights[i][j] = random_double(-1,1);	
		}	
	}
	
	cout <<"\nTrain data: " << n_data << endl;
	cout <<"\nInitial weights:" << endl;
	printMatrix(weights);
	cout <<"\nInitial bias in the output layer:" << endl;
	printVector(bias_output);

	/*

	For loop : from first shuffled sample / train data to the last train data

	*/
	for (int epoch = 0; epoch < epochs; ++epoch) 
	{
		double epoch_loss = 0.0;
		cout <<"\nEpoch:" << epoch+1 << endl;
		for (int sample = 0; sample < n_data; ++sample)
		{
			/*

				Predict each sample based on current weights
	
			*/
			//cout << "\n***************************************************************************************" << endl;
			//cout <<"\n \t\t\t\t\t Next data" << endl;
			//cout << "\n***************************************************************************************" << endl;
			//cout <<"\nSample:" << sample+1 << endl;
			if(vec_output[sample] == "\"Setosa\"")
			{
				for(int j = 0; j < n_classification; ++j)
				{
					actual_output[j] = 0;
				}
				actual_output[0] = 1;
			}
			else if(vec_output[sample] == "\"Versicolor\"")
			{
				for(int j = 0; j < n_classification; ++j)
				{
					actual_output[j] = 0;
				}
				actual_output[1] = 1;			
			}
			else if(vec_output[sample] == "\"Virginica\"")
			{
				for(int j = 0; j < n_classification; ++j)
				{
					actual_output[j] = 0;
				}
				actual_output[2] = 1;

			}
			vector<double> input = getRow(input_matrix,sample);
			//cout <<"\nInput layer:" << endl;
			//printVector(input);	

			//cout <<"\nOutput layer:" << endl;
			//printVector(actual_output);
			for (int j = 0; j < n_classification; ++j) 
			{
				predicted_output[j] = 0;
			}
			for (int j = 0; j < n_classification; ++j) 
			{
				// Linear transformation: Wx + b
				for (int i = 0; i < n_input; ++i) 
				{
					predicted_output[j] += weights[i][j] * input[i];
				}
				// Add bias vector
				predicted_output[j] += bias_output[j];
				// Activation function
				predicted_output[j] = SigmoidActivationfunction(predicted_output[j]);
			}
			/*

			 Calculate the error

			*/
			for (int j = 0; j < n_classification; ++j) 
			{
				error[j] = actual_output[j] - predicted_output[j];
			}	
			//cout <<"\nPredicted Output:" << endl;
			//printVector(predicted_output);	
			//cout <<"\nError vector:" << endl;
			//printVector(error);	
			//cout <<"\n||Error|| = "<< norm(error) << endl;
			epoch_loss += norm(error);

			/*

			 Backpropagation with conjugate gradient to update weights

			*/
			for (int j = 0; j < n_classification; ++j) 
			{
				delta[j] = SigmoidDerivative(predicted_output[j])*(error[j]);
			}
			//cout <<"\nDelta output vector:" << endl;
			//printVector(delta);	
			for (int j = 0; j < n_classification; ++j) 
			{
				bias_output[j] += (learning_rate*delta[j]);
				for (int i = 0; i < n_input; ++i) 
				{
					weights[i][j] +=  learning_rate*delta[j]*input[i];
				}
			}	
			
				//cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
				//printMatrix(weights);	
				
				//cout <<"\nNew bias in the output layer:" << endl;
				//printVector(bias_output);
				
		}	
		cout <<"\nEpoch loss:" << epoch_loss << endl;
	}
	cout <<"\nModel matrix (from input to output layer):" << endl;
	printMatrix(weights);	
	saveMatrixdouble(weights,"modelweights.txt");
	saveVectordouble(bias_output,"bias_output.txt");
}

void FNN_nohiddenlayer_irisdataset_training_conjugategradient(vector<vector<double>>& input_matrix, vector<string>& vec_output )
{
// May 11th, 2026

	int n_classification = 3;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();
	// n_hiddenlayer is the number of nodes in the hidden layer
	vector<vector<double>> weights(n_input, vector<double>(n_classification));
        vector<double> bias_output(n_classification, 0.0); // biases set to 0
	vector<double> predicted_output(n_classification, 0.0);
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	vector<double> delta(n_classification, 0.0);
	vector<double> vec_gk_input(n_input, 0.0);
	vector<double> vec_gk_input_current(n_input*n_classification, 0.0);
	vector<double> vec_gk_input_previous(n_input*n_classification, 0.0);

	vector<double> vec_dk_input(n_input, 0.0);
	vector<double> vec_dk_input_previous(n_input*n_classification, 0.0);
	vector<double> vec_dk_input_current(n_input*n_classification, 0.0);

	double learning_rate = 1;
	int epochs = 1000;
	for (int i = 0; i < n_input; ++i)
	{
		bias_output[i] = random_double(-1,1);
		for(int j = 0; j < n_classification; ++j)
		{
			weights[i][j] = random_double(-1,1);	
		}	
	}
	
	cout <<"\nTrain data: " << n_data << endl;
	cout <<"\nInitial weights:" << endl;
	printMatrix(weights);
	cout <<"\nInitial bias in the output layer:" << endl;
	printVector(bias_output);

	/*

	For loop : from first shuffled sample / train data to the last train data

	*/
	for (int epoch = 0; epoch < epochs; ++epoch) 
	{
		double epoch_loss = 0.0;
		cout <<"\nEpoch:" << epoch+1 << endl;
		for (int sample = 0; sample < n_data; ++sample)
		{
			/*

				Predict each sample based on current weights
	
			*/
			//cout << "\n***************************************************************************************" << endl;
			//cout <<"\n \t\t\t\t\t Next data" << endl;
			//cout << "\n***************************************************************************************" << endl;
			//cout <<"\nSample:" << sample+1 << endl;
			if(vec_output[sample] == "\"Setosa\"")
			{
				for(int j = 0; j < n_classification; ++j)
				{
					actual_output[j] = 0;
				}
				actual_output[0] = 1;
			}
			else if(vec_output[sample] == "\"Versicolor\"")
			{
				for(int j = 0; j < n_classification; ++j)
				{
					actual_output[j] = 0;
				}
				actual_output[1] = 1;			
			}
			else if(vec_output[sample] == "\"Virginica\"")
			{
				for(int j = 0; j < n_classification; ++j)
				{
					actual_output[j] = 0;
				}
				actual_output[2] = 1;

			}
			vector<double> input = getRow(input_matrix,sample);
			//cout <<"\nInput layer:" << endl;
			//printVector(input);	

			//cout <<"\nOutput layer:" << endl;
			//printVector(actual_output);
			for (int j = 0; j < n_classification; ++j) 
			{
				predicted_output[j] = 0;
			}
			for (int j = 0; j < n_classification; ++j) 
			{
				// Linear transformation: Wx + b
				for (int i = 0; i < n_input; ++i) 
				{
					predicted_output[j] += weights[i][j] * input[i];
				}
				// Add bias vector
				predicted_output[j] += bias_output[j];
				// Activation function
				predicted_output[j] = SigmoidActivationfunction(predicted_output[j]);
			}
			/*

			 Calculate the error

			*/
			for (int j = 0; j < n_classification; ++j) 
			{
				error[j] = actual_output[j] - predicted_output[j];
			}	
			//cout <<"\nPredicted Output:" << endl;
			//printVector(predicted_output);	
			//cout <<"\nError vector:" << endl;
			//printVector(error);	
			//cout <<"\n||Error|| = "<< norm(error) << endl;
			//epoch_loss += norm(error);

			//double beta_PR_input = 0;
			/*

			 Backpropagation with conjugate gradient to update weights

			*/
			for (int j = 0; j < n_classification; ++j) 
			{
				delta[j] = SigmoidDerivative(predicted_output[j])*(error[j]);
			}
			//cout <<"\nDelta output vector:" << endl;
			//printVector(delta);	
			for (int j = 0; j < n_classification; ++j) 
			{
				bias_output[j] += (learning_rate*delta[j]);
				for (int i = 0; i < n_input; ++i) 
				{
					weights[i][j] +=  learning_rate*delta[j]*input[i];
				}
			}	
			/*

			  To compute Beta with Polnak Ribiere formula for gradient of the error function with respect to weights from input layer to the output layer

			*/	
			/*
			if( epoch >=1) 
			{
				for (int i = 0; i < n_input*n_classification; ++i) 
				{
					vec_gk_input_previous[i] = 0;
				}			
				for (int i = 0; i < n_input*n_classification; ++i) 
				{
					vec_gk_input_previous[i] = vec_gk_input_current[i];
				}
				vec_gk_input_current.clear();
				for (int i = 0; i < n_classification; ++i) 
				{
					for (int j = 0; j < n_input; ++j) 
					{
						vec_gk_input_current.push_back(delta[i]*input[j]);
					}
				}	
				vector<double> vec_gk_diff;
					
				for (int i = 0; i < n_input*n_classification; ++i) 
				{
					vec_gk_diff.push_back(vec_gk_input_current[i] - vec_gk_input_previous[i]);
				}
				beta_PR_input = divisiond(dot(vec_gk_input_current,vec_gk_diff), dot(vec_gk_input_previous,vec_gk_input_previous));
			}
				
			if( epoch == 0) 
			{					
				vec_gk_input_current.clear();
				for (int i = 0; i < n_classification; ++i) 
				{
					bias_output[i] += (learning_rate*delta[i]);
					for (int j = 0; j < n_input; ++j) 
					{
						vec_gk_input_current.push_back(delta[i]*input[j]);
						vec_gk_input[j] = delta[i]*input[j];
						vec_dk_input[j] = vec_gk_input[j] ;
						//vec_dk_input[j] = -vec_gk_input[j] ;
						vec_dk_input_previous.push_back(-vec_gk_input[j]);

						weights[i][j] +=  (learning_rate*vec_dk_input[j]);
							
					}
				}	
			}
			if( epoch >=1 ) 
			{
				vec_gk_input_current.clear();
					//learning_rate = line_search;
					
				for (int i = 0; i < n_classification; ++i) 
				{
					bias_output[i] += (learning_rate*delta[i]);
					for (int j = 0; j < n_input; ++j) 
					{
						vec_gk_input_current.push_back(delta[i]*input[j]);
						vec_gk_input[j] = delta[i]*input[j];
						vec_dk_input[j] = vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j];
						//vec_dk_input[j] = -vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j];

						weights[i][j] +=  (learning_rate*vec_dk_input[j]);
							
					}
				}	
				vec_dk_input_current.clear();
			
				for (int i = 0; i < n_classification; ++i) 
				{
					for (int j = 0; j < n_input; ++j) 
					{
						vec_dk_input_current.push_back(vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j]);
						//vec_dk_input_current.push_back(-vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j]);
					}
				}	
				vec_dk_input_previous.clear();
						
				for (int i = 0; i < n_classification; ++i) 
				{
					for (int j = 0; j < n_input; ++j) 
					{
						vec_dk_input_previous.push_back(vec_dk_input_current[i+j]);
					}
				}	

			}
			*/
				//cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
				//printMatrix(weights);	
				
				//cout <<"\nNew bias in the output layer:" << endl;
				//printVector(bias_output);
				
		}	
		cout <<"\nEpoch loss:" << epoch_loss << endl;
	}
	cout <<"\nModel matrix (from input to output layer):" << endl;
	printMatrix(weights);	
	saveMatrixdouble(weights,"modelweights.txt");
	saveVectordouble(bias_output,"bias_output.txt");
}


void FNN_nohiddenlayer_irisdataset_testing(vector<vector<double>>& input_matrix, vector<string>& vec_output )
{
// May 11th, 2026

	int n_classification = 3;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();
	int correct = 0;
	int wrong = 0;	
	vector<vector<double>> weights = loadMatrixFromFile("modelweights.txt");
        vector<double> bias_output = loadVectorFromFile("bias_output.txt");
	vector<double> vec_error;
	vector<double> vec_CCEloss;
	vector<string> vec_correctorwrong;

	cout <<"\nTesting data: " << n_data << endl;
	cout <<"\nModel weights from input to output layer:" << endl;
	printMatrix(weights);
        
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	vector<vector<double>> mat_correct;	
	double correct_setosa=0, correct_virginica=0, correct_versicolor=0;
	for (int sample = 0; sample < n_data; ++sample)
	{
		vector<double> predicted_output(n_classification, 0.0);
		cout << "\n***************************************************************************************" << endl;
		cout <<"\n \t\t\t\t\t Next data" << endl;
		cout << "\n***************************************************************************************" << endl;
		cout <<"\nTesting sample:" << sample+1 << endl;
		if(vec_output[sample] == "\"Setosa\"")
		{
			for(int j = 0; j < n_classification; ++j)
			{
				actual_output[j] = 0;
			}
			actual_output[0] = 1;
		}
		else if(vec_output[sample] == "\"Versicolor\"")
		{
			for(int j = 0; j < n_classification; ++j)
			{
				actual_output[j] = 0;
			}
			actual_output[1] = 1;			
		}
		else if(vec_output[sample] == "\"Virginica\"")
		{
			for(int j = 0; j < n_classification; ++j)
			{
				actual_output[j] = 0;
			}
			actual_output[2] = 1;

		}
		vector<double> input = getRow(input_matrix,sample);
		cout <<"\nInput layer:" << endl;
		printVector(input);	

		cout <<"\nOutput layer:" << endl;
		printVector(actual_output);

		for (int j = 0; j < n_classification; ++j) 
		{
			// Linear transformation: Wx + b
			for (int i = 0; i < n_input; ++i) 
			{
				predicted_output[j] += weights[i][j] * input[i];
			}
			predicted_output[j] += bias_output[j];
			// Activation function
			predicted_output[j] = SigmoidActivationfunction(predicted_output[j]);
		}
	
		for (int j = 0; j < n_classification; ++j) 
		{
			error[j] = actual_output[j] - predicted_output[j];
		}	

		int max_index = MaxElementIndex(predicted_output);
		if(max_index==0)
		{
			if(vec_output[sample] == "\"Setosa\"")
			{
				vec_correctorwrong.push_back("Correct\n");
				correct +=1;
				correct_setosa +=1;
			}
			else 
			{
				vec_correctorwrong.push_back("Wrong\n");
				wrong +=1;
			}
		}
		else if(max_index==1)
		{
			if(vec_output[sample] == "\"Versicolor\"")
			{
				vec_correctorwrong.push_back("Correct\n");
				correct +=1;
				correct_versicolor +=1;
			}
			else 
			{
				vec_correctorwrong.push_back("Wrong\n");
				wrong +=1;
			}
		}
		else if(max_index==2)
		{
			if(vec_output[sample] == "\"Virginica\"")
			{
				vec_correctorwrong.push_back("Correct\n");
				correct +=1;
				correct_virginica+=1;
			}
			else 
			{
				vec_correctorwrong.push_back("Wrong\n");
				wrong +=1;
			}

		}

		cout <<"\nPredicted Output:" << endl;
		printVector(predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;
		vec_error.push_back(norm(error));
		for(int j = 0; j < n_classification; ++j)
		{
			if(actual_output[j] == 1)
			{
				vec_CCEloss.push_back(-log(predicted_output[j]));
			}
		}	
	}	
	mat_correct.push_back({double(correct), double(correct_setosa), double(correct_versicolor), double(correct_virginica), double(1)});
	//cout <<"Error vector : " << endl;
	//printVector(vec_error);
	//cout <<"\nCCE loss vector : " << endl;
	//printVector(vec_CCEloss);
	double CCE_accumulate = std::accumulate(vec_CCEloss.begin(), vec_CCEloss.end(), 0.0) ;
	double batch_loss = divisiond(CCE_accumulate,n_data);
	cout <<"\nAverage Loss over a batch of " << n_data << " samples : " << batch_loss << endl;
	cout <<"\nCorrect prediction : " << correct << endl;
	cout <<"Wrong prediction : " << wrong << endl;
	//printStringVector(vec_correctorwrong);
	cout <<"Correct matrix: " << endl;
	cout <<setw(23) <<"Total Correct" <<setw(23) << "Setosa" <<setw(23) << "Versicolor" << setw(23) <<"Virginica" << setw(23) << "Model" << endl;
	printMatrix(mat_correct);
}

void FNN_1hiddenlayer_irisdataset_training(vector<vector<double>>& input_matrix, vector<string>& vec_output )
{
// May 2nd, 2026 1 day after Full Moon
// Homework: weights matrix size has to be fixed -> input * next layer neuron
	int n_hiddenlayer = 5;
	int n_classification = 3;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();
	// n_hiddenlayer is the number of nodes in the hidden layer
	vector<vector<double>> model_weights(n_hiddenlayer, vector<double>(n_input, 0.0));
	vector<vector<double>> model_hiddenweights(n_classification, vector<double>(n_hiddenlayer, 0.0));
	vector<vector<double>> weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> hiddenweights(n_classification, vector<double>(n_hiddenlayer));
	vector<vector<double>> updated_weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> updated_hiddenweights(n_classification, vector<double>(n_hiddenlayer));
	vector<vector<double>> new_weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> new_hiddenweights(n_classification, vector<double>(n_hiddenlayer));
        vector<double> bias_hiddenlayer1(n_hiddenlayer, 0.0); // biases set to 0
        vector<double> bias_output(n_classification, 0.0); // biases set to 0
	double learning_rate = 0.5;

	for (int i = 0; i < n_hiddenlayer; ++i)
	{
		bias_hiddenlayer1[i] = random_double(-1,1);
		for(int j = 0; j < n_input; ++j)
		{
			weights[i][j] = random_double(-1,1);	
		}	
	}
	for (int i = 0; i < n_classification; ++i)
	{
		bias_output[i] = random_double(-1,1);
		for(int j = 0; j < n_hiddenlayer; ++j)
		{
			hiddenweights[i][j] = random_double(-1,1);	
		}	
	}
	cout <<"\nTrain data: " << n_data << endl;
	cout <<"\nInitial weights from input to hidden layer:" << endl;
	printMatrix(weights);
	cout <<"\nInitial weights from hidden layer to output layer:" << endl;
	printMatrix(hiddenweights);
	cout <<"\nInitial bias in the hidden layer:" << endl;
	printVector(bias_hiddenlayer1);
	cout <<"\nInitial bias in the output layer:" << endl;
	printVector(bias_output);
	
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	
	for (int i = 0; i < n_hiddenlayer; ++i)
	{
		for(int j = 0; j < n_input; ++j)
		{
			updated_weights[i][j] = weights[i][j];	
		}	
	}
	for (int i = 0; i < n_classification; ++i)
	{
		for(int j = 0; j < n_hiddenlayer; ++j)
		{
			updated_hiddenweights[i][j] = hiddenweights[i][j];	
		}	
	}	
	vector<double> new_hiddenlayer(n_hiddenlayer, 0.0);
	vector<double> delta_hiddenweights(n_hiddenlayer, 0.0);
	vector<double> delta(n_classification, 0.0);
	/*

	For loop : from first shuffled sample / train data to the last train data

	*/
	int epochs = 1;
	for (int i = 0; i < epochs ; ++i)
	{
	for (int sample = 0; sample < n_data; ++sample)
	{
		cout << "\n***************************************************************************************" << endl;
		cout <<"\n \t\t\t\t\t Next data" << endl;
		cout << "\n***************************************************************************************" << endl;
		cout <<"\nSample:" << sample+1 << endl;
		if(vec_output[sample] == "\"Setosa\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[0] = 1;
		}
		else if(vec_output[sample] == "\"Versicolor\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[1] = 1;			
		}
		else if(vec_output[sample] == "\"Virginica\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[2] = 1;

		}
		vector<double> input = getRow(input_matrix,sample);
		vector<double> hiddenlayer(n_hiddenlayer, 0.0);
		vector<double> predicted_output(n_classification, 0.0);

		cout <<"\nInput layer:" << endl;
		printVector(input);	

		cout <<"\nOutput layer:" << endl;
		printVector(actual_output);

		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				hiddenlayer[i] += updated_weights[i][j] * input[j];
			}
			hiddenlayer[i] += bias_hiddenlayer1[i];
		}
		cout <<"\nHidden layer:" << endl;
		printVector(hiddenlayer);
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Activation function
			hiddenlayer[i] = TanhActivationfunction(hiddenlayer[i]);
		}
		cout <<"\nHidden layer after activation function:" << endl;
		printVector(hiddenlayer);	
		
		for (int i = 0; i < n_classification; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				predicted_output[i] += updated_hiddenweights[i][j] * hiddenlayer[j];
			}
			predicted_output[i] += bias_output[i];
		}
		cout <<"\nLogits before softmax activation function:" << endl;
		printVector(predicted_output);	
		vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(predicted_output);
		/*

		 Softmax & Calculate the error

		*/
		for (int i = 0; i < n_classification; ++i) 
		{
			// Activation function
			predicted_output[i] = predicted_output_softmax[i];
			error[i] = actual_output[i] - predicted_output[i];
		}	
		cout <<"\nPredicted Output:" << endl;
		printVector(predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;

		
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			for (int j = 0; j < n_input; ++j) 
			{
				new_weights[i][j] = updated_weights[i][j] ;
			}
		}	
		for (int i = 0; i < n_classification; ++i) 
		{
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				new_hiddenweights[i][j] = updated_hiddenweights[i][j] ;
			}
		}	
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			new_hiddenlayer[i] = hiddenlayer[i];
		}
		vector<double> new_predicted_output(n_classification, 0.0);
		for (int i = 0; i < n_classification; ++i) 
		{
			new_predicted_output[i] = predicted_output[i];
		}
		//printMatrix(new_weights);	
		//printMatrix(new_hiddenweights);	

		cout << "\n***************************************************************************************" << endl;
		cout <<"\n \t\t\t Backpropagation iteration" << endl;
		cout << "\n***************************************************************************************" << endl;
		//cout <<"\nLearning rate: " << learning_rate << endl;
		// for loop till converge starts here
		int backpropagation_iter = 3;
		for (int i = 0; i < backpropagation_iter; ++i) 
		{
			/*

			 Backpropagation

			*/

			cout << "\nIteration = " << i+1 << endl;
			
			for (int i = 0; i < n_classification; ++i) 
			{
				delta[i] = TanhDerivative(new_predicted_output[i])*(error[i]);
			}
			//cout <<"\nDelta output vector:" << endl;
			//printVector(delta);	
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				double sum_output = 0;
				for (int j = 0; j < n_classification; ++j) 
				{
					sum_output += delta[j]*new_hiddenweights[j][i];
				}
				delta_hiddenweights[i] = TanhDerivative(new_hiddenlayer[i])*(sum_output);
			}
			//cout <<"\nDelta hidden layer vector:" << endl;
			//printVector(delta_hiddenweights);	

			for (int i = 0; i < n_classification; ++i) 
			{
				bias_output[i] += (learning_rate*delta[i]);
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					new_hiddenweights[i][j] += (learning_rate*delta[i]*new_hiddenlayer[j]);
				}
			}	
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				bias_hiddenlayer1[i] += (learning_rate*delta_hiddenweights[i]);
				for (int j = 0; j < n_input; ++j) 
				{
					new_weights[i][j] +=  (learning_rate*delta_hiddenweights[i]*input[j]);
					
				}
			}	
			
			//cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
			//printMatrix(new_weights);	
			//cout <<"\nNew W1 matrix (from hidden layer to output):" << endl;
			//printMatrix(new_hiddenweights);	
			//cout <<"\nNew bias in the hidden layer:" << endl;
			//printVector(bias_hiddenlayer1);
			//cout <<"\nNew bias in the output layer:" << endl;
			//printVector(bias_output);
			/*

			 Forward Propagation again..

			*/
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				new_hiddenlayer[i] = 0;
			}
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				// Linear transformation: Wx + b
				for (int j = 0; j < n_input; ++j) 
				{
					new_hiddenlayer[i] += new_weights[i][j] * input[j];
				}
				new_hiddenlayer[i] += bias_hiddenlayer1[i];
			}
			//cout <<"\nHidden layer before activation function:" << endl;
			//printVector(new_hiddenlayer);	
		
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				// Activation function
				new_hiddenlayer[i] = TanhActivationfunction(new_hiddenlayer[i]);
			}
			//cout <<"\nHidden layer after activation function:" << endl;
			//printVector(new_hiddenlayer);	
			for (int i = 0; i < n_classification; ++i) 
			{
				new_predicted_output[i] = 0;
			}

			for (int i = 0; i < n_classification; ++i) 
			{
				// Linear transformation: Wx + b
				for (int j = 0; j < n_classification; ++j) 
				{
					new_predicted_output[i] += new_hiddenweights[i][j] * new_hiddenlayer[j];
				}
				new_predicted_output[i] += bias_output[i];
			}
			//cout <<"\nLogits before softmax activation function:" << endl;
			//printVector(new_predicted_output);	
			vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(new_predicted_output);
			for (int i = 0; i < n_classification; ++i) 
			{
				// Activation function
				new_predicted_output[i] = predicted_output_softmax[i];
				error[i] = actual_output[i] - new_predicted_output[i];
				
			}	
			/*

			 Saving updated weights

			*/

			if(i == backpropagation_iter-1)
			{	
				for (int i = 0; i < n_hiddenlayer; ++i)
				{
					for(int j = 0; j < n_input; ++j)
					{
						updated_weights[i][j] = new_weights[i][j];	
						model_weights[i][j] = updated_weights[i][j];
					}	
				}
				for (int i = 0; i < n_classification; ++i)
				{
					for(int j = 0; j < n_hiddenlayer; ++j)
					{
						updated_hiddenweights[i][j] = new_hiddenweights[i][j];	
						model_hiddenweights[i][j] = updated_hiddenweights[i][j];
					}	
				}
			}
			cout <<"\nPredicted Output:" << endl;
			printVector(new_predicted_output);	
			cout <<"\nError vector:" << endl;
			printVector(error);	
			cout <<"\n||Error|| = "<< norm(error) << endl;
			/*if (i == backpropagation_iter-1)
			{
			cout << "\nIteration = " << i+1 << endl;
			cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
			printMatrix(new_weights);	
			cout <<"\nNew W1 matrix (from hidden layer to output):" << endl;
			printMatrix(new_hiddenweights);	

			cout <<"\nPredicted Output:" << endl;
			printVector(new_predicted_output);	
			cout <<"\nError vector:" << endl;
			printVector(error);	
			cout <<"\n||Error|| = "<< norm(error) << endl;

			}*/
		}// end of backpropagating for each sample
	}// end of iterating from sample 1 to last training sample
	}// end of epochs
	cout <<"\nModel W0 matrix (from input to hidden layer):" << endl;
	printMatrix(model_weights);	
	cout <<"\nModel W1 matrix (from hidden layer to output):" << endl;
	printMatrix(model_hiddenweights);	
	saveMatrixdouble(model_weights,"modelweights.txt");
	saveMatrixdouble(model_hiddenweights,"modelhiddenweights.txt");
	saveVectordouble(bias_hiddenlayer1,"bias_hiddenlayer1.txt");
	saveVectordouble(bias_output,"bias_output.txt");
}

void FNN_1hiddenlayer_irisdataset_training_conjugategradient(vector<vector<double>>& input_matrix, vector<string>& vec_output )
{
// May 9th, 2026
// Homework: weights matrix size has to be fixed -> input * next layer neuron
	int n_hiddenlayer = 5;
	int n_classification = 3;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();
	// n_hiddenlayer is the number of nodes in the hidden layer
	vector<vector<double>> model_weights(n_hiddenlayer, vector<double>(n_input, 0.0));
	vector<vector<double>> model_hiddenweights(n_classification, vector<double>(n_hiddenlayer, 0.0));
	vector<vector<double>> weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> hiddenweights(n_classification, vector<double>(n_hiddenlayer));
	vector<vector<double>> updated_weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> updated_hiddenweights(n_classification, vector<double>(n_hiddenlayer));
	vector<vector<double>> new_weights(n_hiddenlayer, vector<double>(n_input));
	vector<vector<double>> new_hiddenweights(n_classification, vector<double>(n_hiddenlayer));
        vector<double> bias_hiddenlayer1(n_hiddenlayer, 0.0); // biases set to 0
        vector<double> bias_output(n_classification, 0.0); // biases set to 0
	double learning_rate = 1;

	for (int i = 0; i < n_hiddenlayer; ++i)
	{
		bias_hiddenlayer1[i] = random_double(-1,1);
		for(int j = 0; j < n_input; ++j)
		{
			weights[i][j] = random_double(-1,1);	
		}	
	}
	for (int i = 0; i < n_classification; ++i)
	{
		bias_output[i] = random_double(-1,1);
		for(int j = 0; j < n_hiddenlayer; ++j)
		{
			hiddenweights[i][j] = random_double(-1,1);	
		}	
	}
	cout <<"\nTrain data: " << n_data << endl;
	cout <<"\nInitial weights from input to hidden layer:" << endl;
	printMatrix(weights);
	cout <<"\nInitial weights from hidden layer to output layer:" << endl;
	printMatrix(hiddenweights);
	cout <<"\nInitial bias in the hidden layer:" << endl;
	printVector(bias_hiddenlayer1);
	cout <<"\nInitial bias in the output layer:" << endl;
	printVector(bias_output);
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
        vector<double> predicted_output(n_classification, 0.0);
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	
	for (int i = 0; i < n_hiddenlayer; ++i)
	{
		for(int j = 0; j < n_input; ++j)
		{
			updated_weights[i][j] = weights[i][j];	
		}	
	}
	for (int i = 0; i < n_classification; ++i)
	{
		for(int j = 0; j < n_hiddenlayer; ++j)
		{
			updated_hiddenweights[i][j] = hiddenweights[i][j];	
		}	
	}	
	/*

	For loop : from first shuffled sample / train data to the last train data

	*/

	for (int sample = 0; sample < n_data; ++sample)
	{
		cout << "\n***************************************************************************************" << endl;
		cout <<"\n \t\t\t\t\t Next data" << endl;
		cout << "\n***************************************************************************************" << endl;
		cout <<"\nSample:" << sample+1 << endl;
		if(vec_output[sample] == "\"Setosa\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[0] = 1;
		}
		else if(vec_output[sample] == "\"Versicolor\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[1] = 1;			
		}
		else if(vec_output[sample] == "\"Virginica\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[2] = 1;

		}
		vector<double> input = getRow(input_matrix,sample);
		cout <<"\nInput layer:" << endl;
		printVector(input);	

		cout <<"\nOutput layer:" << endl;
		printVector(actual_output);

		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				hiddenlayer[i] += updated_weights[i][j] * input[j];
			}
			hiddenlayer[i] += bias_hiddenlayer1[i];
		}
		cout <<"\nHidden layer:" << endl;
		printVector(hiddenlayer);
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Activation function
			hiddenlayer[i] = TanhActivationfunction(hiddenlayer[i]);
		}
		cout <<"\nHidden layer after activation function:" << endl;
		printVector(hiddenlayer);	
		for (int i = 0; i < n_classification; ++i) 
		{
			predicted_output[i] = 0;
		}
		for (int i = 0; i < n_classification; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				predicted_output[i] += updated_hiddenweights[i][j] * hiddenlayer[j];
			}
			predicted_output[i] += bias_output[i];
		}
		cout <<"\nLogits before softmax activation function:" << endl;
		printVector(predicted_output);	
		vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(predicted_output);
		/*

		 Softmax & Calculate the error

		*/
		for (int i = 0; i < n_classification; ++i) 
		{
			// Activation function
			predicted_output[i] = predicted_output_softmax[i];
			error[i] = actual_output[i] - predicted_output[i];
		}	
		cout <<"\nPredicted Output:" << endl;
		printVector(predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;

		vector<vector<double>> new_weights(n_hiddenlayer, vector<double>(n_input));
		vector<vector<double>> new_hiddenweights(n_classification, vector<double>(n_hiddenlayer));
		vector<double> new_hiddenlayer(weights.size(), 0.0);
		vector<double> new_predicted_output(n_classification, 0.0);
		vector<double> delta_hiddenweights(n_hiddenlayer, 0.0);
		vector<double> delta(n_classification, 0.0);
		vector<double> vec_gk_hiddenlayer(n_hiddenlayer, 0.0);
		vector<double> vec_gk_input(n_input, 0.0);

		vector<double> vec_gk_hiddenlayer_current(n_hiddenlayer*n_classification, 0.0);
		vector<double> vec_gk_input_current(n_input*n_hiddenlayer, 0.0);
		vector<double> vec_gk_hiddenlayer_previous(n_hiddenlayer*n_classification, 0.0);
		vector<double> vec_gk_input_previous(n_input*n_hiddenlayer, 0.0);

		vector<double> vec_dk_hiddenlayer(n_hiddenlayer, 0.0);
		vector<double> vec_dk_input(n_input, 0.0);
		vector<double> vec_dk_hiddenlayer_previous(n_hiddenlayer*n_classification, 0.0);
		vector<double> vec_dk_input_previous(n_input*n_hiddenlayer, 0.0);
		vector<double> vec_dk_hiddenlayer_current(n_hiddenlayer*n_classification, 0.0);
		vector<double> vec_dk_input_current(n_input*n_hiddenlayer, 0.0);
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			for (int j = 0; j < n_input; ++j) 
			{
				new_weights[i][j] = updated_weights[i][j] ;
			}
		}	
		for (int i = 0; i < n_classification; ++i) 
		{
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				new_hiddenweights[i][j] = updated_hiddenweights[i][j] ;
			}
		}	
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			new_hiddenlayer[i] = hiddenlayer[i];
		}
		
		for (int i = 0; i < n_classification; ++i) 
		{
			new_predicted_output[i] = predicted_output[i];
		}
		//printMatrix(new_weights);	
		//printMatrix(new_hiddenweights);	

		cout << "\n***************************************************************************************" << endl;
		cout <<"\n \t\t\t Backpropagation iteration" << endl;
		cout << "\n***************************************************************************************" << endl;
		//cout <<"\nLearning rate: " << learning_rate << endl;
		// for loop till converge starts here
		int backpropagation_iter = 3;
		double beta_PR_input = 0;
		double beta_PR_hiddenlayer = 0;
		for (int i = 0; i < backpropagation_iter; ++i) 
		{
			/*

			 Backpropagation with conjugate gradient

			*/

			cout << "\nIteration = " << i+1 << endl;
			
			for (int i = 0; i < n_classification; ++i) 
			{
				delta[i] = TanhDerivative(new_predicted_output[i])*(error[i]);
			}
			//cout <<"\nDelta output vector:" << endl;
			//printVector(delta);	
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				double sum_output = 0;
				for (int j = 0; j < n_classification; ++j) 
				{
					sum_output += delta[j]*new_hiddenweights[j][i];
				}
				delta_hiddenweights[i] = TanhDerivative(new_hiddenlayer[i])*(sum_output);
			}
			//cout <<"\nDelta hidden layer vector:" << endl;
			//printVector(delta_hiddenweights);	
			
			if( i >=1) // to compute Beta with Polnak Ribiere formula for gradient of the error function with respect to weights from hidden layer to the output layer
			{
				
				for (int i = 0; i < n_classification*n_hiddenlayer; ++i) 
				{
					vec_gk_hiddenlayer_previous[i] = 0;
				}			
				for (int i = 0; i < n_classification*n_hiddenlayer; ++i) 
				{
					vec_gk_hiddenlayer_previous[i] = vec_gk_hiddenlayer_current[i];
				}
				vec_gk_hiddenlayer_current.clear();
				for (int i = 0; i < n_classification; ++i) 
				{
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						vec_gk_hiddenlayer_current.push_back(delta[i]*new_hiddenlayer[j]);
					}
				}	
				//printVector(vec_gk_hiddenlayer_current);
				vector<double> vec_gk_diff;
				
				for (int i = 0; i < n_classification*n_hiddenlayer; ++i) 
				{
					vec_gk_diff.push_back(vec_gk_hiddenlayer_current[i] - vec_gk_hiddenlayer_previous[i]);
				}
				//printVector(vec_gk_diff);
				beta_PR_hiddenlayer = divisiond(dot(vec_gk_hiddenlayer_current,vec_gk_diff), dot(vec_gk_hiddenlayer_previous,vec_gk_hiddenlayer_previous));
				
				
			}
			if( i >=1) // to compute Beta with Polnak Ribiere formula for gradient of the error function with respect to weights from input layer to the hidden layer 
			{
				for (int i = 0; i < n_input*n_hiddenlayer; ++i) 
				{
					vec_gk_input_previous[i] = 0;
				}			
				for (int i = 0; i < n_input*n_hiddenlayer; ++i) 
				{
					vec_gk_input_previous[i] = vec_gk_input_current[i];
				}
				vec_gk_input_current.clear();
				for (int i = 0; i < n_hiddenlayer; ++i) 
				{
					for (int j = 0; j < n_input; ++j) 
					{
						vec_gk_input_current.push_back(delta_hiddenweights[i]*input[j]);
					}
				}	
				vector<double> vec_gk_diff;
				
				for (int i = 0; i < n_input*n_hiddenlayer; ++i) 
				{
					vec_gk_diff.push_back(vec_gk_input_current[i] - vec_gk_input_previous[i]);
				}
				beta_PR_input = divisiond(dot(vec_gk_input_current,vec_gk_diff), dot(vec_gk_input_previous,vec_gk_input_previous));
			}
			
			if( i == 0) 
			{
				
				vec_gk_hiddenlayer_current.clear();
				vec_gk_input_current.clear();
				for (int i = 0; i < n_classification; ++i) 
				{
					bias_output[i] += (learning_rate*delta[i]);
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						vec_gk_hiddenlayer_current.push_back(delta[i]*new_hiddenlayer[j]);
						vec_gk_hiddenlayer[j] = delta[i]*new_hiddenlayer[j];
						vec_dk_hiddenlayer[j] = vec_gk_hiddenlayer[j] ;
						//vec_dk_hiddenlayer[j] = -vec_gk_hiddenlayer[j] ;
						vec_dk_hiddenlayer_previous.push_back(-vec_gk_hiddenlayer[j]);

						new_hiddenweights[i][j] += (learning_rate*vec_dk_hiddenlayer[j]);
					}
				}	
				for (int i = 0; i < n_hiddenlayer; ++i) 
				{
					bias_hiddenlayer1[i] += (learning_rate*delta_hiddenweights[i]);
					for (int j = 0; j < n_input; ++j) 
					{
						vec_gk_input_current.push_back(delta_hiddenweights[i]*input[j]);
						vec_gk_input[j] = delta_hiddenweights[i]*input[j];
						vec_dk_input[j] = vec_gk_input[j] ;
						//vec_dk_input[j] = -vec_gk_input[j] ;
						vec_dk_input_previous.push_back(-vec_gk_input[j]);

						new_weights[i][j] +=  (learning_rate*vec_dk_input[j]);
						
					}
				}	
			}
			if( i >=1 ) 
			{
				vec_gk_hiddenlayer_current.clear();
				vec_gk_input_current.clear();
				//learning_rate = line_search;
				for (int i = 0; i < n_classification; ++i) 
				{
					bias_output[i] += (learning_rate*delta[i]);
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						vec_gk_hiddenlayer_current.push_back(delta[i]*new_hiddenlayer[j]);
						vec_gk_hiddenlayer[j] = delta[i]*new_hiddenlayer[j];
						vec_dk_hiddenlayer[j] = vec_gk_hiddenlayer[j] + beta_PR_hiddenlayer*vec_dk_hiddenlayer_previous[i+j];
						//vec_dk_hiddenlayer[j] = -vec_gk_hiddenlayer[j] + beta_PR_hiddenlayer*vec_dk_hiddenlayer_previous[i+j];

						new_hiddenweights[i][j] += (learning_rate*vec_dk_hiddenlayer[j]);
					}
				}	
				
				for (int i = 0; i < n_hiddenlayer; ++i) 
				{
					bias_hiddenlayer1[i] += (learning_rate*delta_hiddenweights[i]);
					for (int j = 0; j < n_input; ++j) 
					{
						vec_gk_input_current.push_back(delta_hiddenweights[i]*input[j]);
						vec_gk_input[j] = delta_hiddenweights[i]*input[j];
						vec_dk_input[j] = vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j];
						//vec_dk_input[j] = -vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j];

						new_weights[i][j] +=  (learning_rate*vec_dk_input[j]);
						
					}
				}	
				vec_dk_input_current.clear();
				vec_dk_hiddenlayer_current.clear();
				for (int i = 0; i < n_classification; ++i) 
				{
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						vec_dk_hiddenlayer_current.push_back(vec_gk_hiddenlayer[j] + beta_PR_hiddenlayer*vec_dk_hiddenlayer_previous[i+j]);
						//vec_dk_hiddenlayer_current.push_back(-vec_gk_hiddenlayer[j] + beta_PR_hiddenlayer*vec_dk_hiddenlayer_previous[i+j]);
					}
				}	
				for (int i = 0; i < n_hiddenlayer; ++i) 
				{
					for (int j = 0; j < n_input; ++j) 
					{
						vec_dk_input_current.push_back(vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j]);
						//vec_dk_input_current.push_back(-vec_gk_input[j] + beta_PR_input*vec_dk_input_previous[i+j]);
					}
				}	
				vec_dk_input_previous.clear();
				vec_dk_hiddenlayer_previous.clear();
				for (int i = 0; i < n_classification; ++i) 
				{
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						vec_dk_hiddenlayer_previous.push_back(vec_dk_hiddenlayer_current[i+j]);
					}
				}	
				for (int i = 0; i < n_hiddenlayer; ++i) 
				{
					for (int j = 0; j < n_input; ++j) 
					{
						vec_dk_input_previous.push_back(vec_dk_input_current[i+j]);
					}
				}	

			}
			//cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
			//printMatrix(new_weights);	
			//cout <<"\nNew W1 matrix (from hidden layer to output):" << endl;
			//printMatrix(new_hiddenweights);	
			//cout <<"\nNew bias in the hidden layer:" << endl;
			//printVector(bias_hiddenlayer1);
			//cout <<"\nNew bias in the output layer:" << endl;
			//printVector(bias_output);
			/*

			 Forward Propagation again..

			*/
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				new_hiddenlayer[i] = 0;
			}
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				// Linear transformation: Wx + b
				for (int j = 0; j < n_input; ++j) 
				{
					new_hiddenlayer[i] += new_weights[i][j] * input[j];
				}
				new_hiddenlayer[i] += bias_hiddenlayer1[i];
			}
			//cout <<"\nHidden layer before activation function:" << endl;
			//printVector(new_hiddenlayer);	
		
			for (int i = 0; i < n_hiddenlayer; ++i) 
			{
				// Activation function
				new_hiddenlayer[i] = TanhActivationfunction(new_hiddenlayer[i]);
			}
			//cout <<"\nHidden layer after activation function:" << endl;
			//printVector(new_hiddenlayer);	
			for (int i = 0; i < n_classification; ++i) 
			{
				new_predicted_output[i] = 0;
			}

			for (int i = 0; i < n_classification; ++i) 
			{
				// Linear transformation: Wx + b
				for (int j = 0; j < n_classification; ++j) 
				{
					new_predicted_output[i] += new_hiddenweights[i][j] * new_hiddenlayer[j];
				}
				new_predicted_output[i] += bias_output[i];
			}
			//cout <<"\nLogits before softmax activation function:" << endl;
			//printVector(new_predicted_output);	
			vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(new_predicted_output);
			for (int i = 0; i < n_classification; ++i) 
			{
				// Activation function
				new_predicted_output[i] = predicted_output_softmax[i];
				error[i] = actual_output[i] - new_predicted_output[i];
				
			}	
			/*

			 Saving updated weights

			*/

			if(i == backpropagation_iter-1)
			{	
				for (int i = 0; i < n_hiddenlayer; ++i)
				{
					for(int j = 0; j < n_input; ++j)
					{
						updated_weights[i][j] = new_weights[i][j];	
						model_weights[i][j] = updated_weights[i][j];
					}	
				}
				for (int i = 0; i < n_classification; ++i)
				{
					for(int j = 0; j < n_hiddenlayer; ++j)
					{
						updated_hiddenweights[i][j] = new_hiddenweights[i][j];	
						model_hiddenweights[i][j] = updated_hiddenweights[i][j];
					}	
				}
			}
			cout <<"\nPredicted Output:" << endl;
			printVector(new_predicted_output);	
			cout <<"\nError vector:" << endl;
			printVector(error);	
			cout <<"\n||Error|| = "<< norm(error) << endl;
			/*if (i == backpropagation_iter-1)
			{
			cout << "\nIteration = " << i+1 << endl;
			cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
			printMatrix(new_weights);	
			cout <<"\nNew W1 matrix (from hidden layer to output):" << endl;
			printMatrix(new_hiddenweights);	

			cout <<"\nPredicted Output:" << endl;
			printVector(new_predicted_output);	
			cout <<"\nError vector:" << endl;
			printVector(error);	
			cout <<"\n||Error|| = "<< norm(error) << endl;

			}*/
		}
	}	
	cout <<"\nModel W0 matrix (from input to hidden layer):" << endl;
	printMatrix(model_weights);	
	cout <<"\nModel W1 matrix (from hidden layer to output):" << endl;
	printMatrix(model_hiddenweights);	
	saveMatrixdouble(model_weights,"modelweights.txt");
	saveMatrixdouble(model_hiddenweights,"modelhiddenweights.txt");
	saveVectordouble(bias_hiddenlayer1,"bias_hiddenlayer1.txt");
	saveVectordouble(bias_output,"bias_output.txt");
}

void FNN_1hiddenlayer_irisdataset_testing(vector<vector<double>>& input_matrix, vector<string>& vec_output )
{
// May 2nd, 2026 1 day after Full Moon
// Homework: weights matrix size has to be fixed -> input * next layer neuron
	int n_hiddenlayer = 5;
	int n_classification = 3;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();
	int correct = 0;
	int wrong = 0;	
	vector<vector<double>> weights = loadMatrixFromFile("modelweights.txt");
	vector<vector<double>> hiddenweights = loadMatrixFromFile("modelhiddenweights.txt");
        vector<double> bias_hiddenlayer1 = loadVectorFromFile("bias_hiddenlayer1.txt");
        vector<double> bias_output = loadVectorFromFile("bias_output.txt");
	vector<double> vec_error;
	vector<double> vec_CCEloss;
	vector<string> vec_correctorwrong;

	cout <<"\nTesting data: " << n_data << endl;
	cout <<"\nModel weights from input to hidden layer:" << endl;
	printMatrix(weights);
	cout <<"\nModel weights from hidden layer to output layer:" << endl;
	printMatrix(hiddenweights);
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
        
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
		
	for (int sample = 0; sample < n_data; ++sample)
	{
		vector<double> predicted_output(n_classification, 0.0);
		cout << "\n***************************************************************************************" << endl;
		cout <<"\n \t\t\t\t\t Next data" << endl;
		cout << "\n***************************************************************************************" << endl;
		cout <<"\nTesting sample:" << sample+1 << endl;
		if(vec_output[sample] == "\"Setosa\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[0] = 1;
		}
		else if(vec_output[sample] == "\"Versicolor\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[1] = 1;			
		}
		else if(vec_output[sample] == "\"Virginica\"")
		{
			for(int i = 0; i < n_classification; ++i)
			{
				actual_output[i] = 0;
			}
			actual_output[2] = 1;

		}
		vector<double> input = getRow(input_matrix,sample);
		cout <<"\nInput layer:" << endl;
		printVector(input);	

		cout <<"\nOutput layer:" << endl;
		printVector(actual_output);

		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_input; ++j) 
			{
				hiddenlayer[i] += weights[i][j] * input[j];
			}
			hiddenlayer[i] += bias_hiddenlayer1[i];
		}
		
		for (int i = 0; i < n_hiddenlayer; ++i) 
		{
			// Activation function
			hiddenlayer[i] = TanhActivationfunction(hiddenlayer[i]);
		}
		//cout <<"\nHidden layer:" << endl;
		//printVector(hiddenlayer);	
		for (int i = 0; i < n_classification; ++i) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				predicted_output[i] += hiddenweights[i][j] * hiddenlayer[j];
			}
			predicted_output[i] += bias_output[i];
		}
		vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(predicted_output);
		for (int i = 0; i < n_classification; ++i) 
		{
			// Activation function
			predicted_output[i] = predicted_output_softmax[i];
			error[i] = actual_output[i] - predicted_output[i];
		}	
	

		int max_index = MaxElementIndex(predicted_output);
		if(max_index==0)
		{
			if(vec_output[sample] == "\"Setosa\"")
			{
				vec_correctorwrong.push_back("Correct\n");
				correct +=1;
			}
			else 
			{
				vec_correctorwrong.push_back("Wrong\n");
				wrong +=1;
			}
		}
		else if(max_index==1)
		{
			if(vec_output[sample] == "\"Versicolor\"")
			{
				vec_correctorwrong.push_back("Correct\n");
				correct +=1;
			}
			else 
			{
				vec_correctorwrong.push_back("Wrong\n");
				wrong +=1;
			}
		}
		else if(max_index==2)
		{
			if(vec_output[sample] == "\"Virginica\"")
			{
				vec_correctorwrong.push_back("Correct\n");
				correct +=1;
			}
			else 
			{
				vec_correctorwrong.push_back("Wrong\n");
				wrong +=1;
			}

		}

		cout <<"\nPredicted Output:" << endl;
		printVector(predicted_output);	
		cout <<"\nError vector:" << endl;
		printVector(error);	
		cout <<"\n||Error|| = "<< norm(error) << endl;
		vec_error.push_back(norm(error));
		for(int i = 0; i < n_classification; ++i)
		{
			if(actual_output[i] == 1)
			{
				vec_CCEloss.push_back(-log(predicted_output[i]));
			}
		}
	}	
	//cout <<"Error vector : " << endl;
	//printVector(vec_error);
	//cout <<"\nCCE loss vector : " << endl;
	//printVector(vec_CCEloss);
	double CCE_accumulate = std::accumulate(vec_CCEloss.begin(), vec_CCEloss.end(), 0.0) ;
	double batch_loss = divisiond(CCE_accumulate,n_data);
	cout <<"\nAverage Loss over a batch of " << n_data << " samples : " << batch_loss << endl;
	cout <<"\nCorrect prediction : " << correct << endl;
	cout <<"Wrong prediction : " << wrong << endl;
	printStringVector(vec_correctorwrong);
	
}


#endif
#endif