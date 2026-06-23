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

#include <algorithm>
#include <iomanip>
#include <limits>
#include <climits> // for INT_MIN

#include <random> // For random number generation
#include <chrono>
#include <unordered_map>
#include <iomanip>
#include <numeric>
#include <thread>

// ANSI color codes
#define RESET   "\033[0m"
#define RED     "\033[31m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[33m"
#define BLUE    "\033[34m"
#define MAGENTA "\033[35m"
#define CYAN    "\033[36m"
#define BOLD    "\033[1m"

// Debug macro
#define DEBUG(x) std::cout << YELLOW << "DEBUG: " << x << RESET << std::endl

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_DEFINE

using namespace std;

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
	//return SigmoidActivationfunction(x)*(1-SigmoidActivationfunction(x));
}
double TanhDerivative(double x) 
{
	//double t = std::tanh(x);
	//return 1.0 - t * t;
	return 1 - (x*x);
	
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




void print_centered(const std::string& text, char fill = '=', int width = 60) 
{
	int padding = (width - text.length()) / 2;
	std::cout << std::string(padding, fill) << text << std::string(padding, fill) << std::endl;
}

double get_double_input(const string& prompt) 
{
	double value;
	string input;
	while (true) 
	{
		std::cout << prompt;
		std::getline(std::cin, input);
		std::istringstream iss(input);
		if (iss >> value && iss.eof()) 
		{
			return value;
		} 
		else 
		{
			cout << RED << "Invalid input. Please enter a valid number." << RESET << endl;
		}
	}
}

void print_progress_bar(double percentage) 
{
	int width = 30;
	int pos = width * percentage;
	cout << "[";
	for (int i = 0; i < width; ++i) 
	{
	if (i < pos) 
	{
		cout << "=";
	}
	else if (i == pos) 
	{
		cout << ">";
        }
	else 
	{
		cout << " ";
	}
	}
	cout << "] " << int(percentage * 100.0) << "%\r";
	cout.flush();
}

// Simple delay function
void delay() 
{
	for(int i = 0; i < 510000; ++i) 
	{
        // Do nothing, just waste some cycles
	}
}

vector<vector<int>> CNN_2DpadBorder(const vector<vector<int>>& matrix, int pad_size) 
{
	int pad_value = 0;
	if (matrix.empty())
	{ 
		return {};
	}

	int old_rows = matrix.size();
	int old_cols = matrix[0].size(); // Assumes a rectangular input matrix
    
	int new_rows = old_rows + 2 * pad_size;
	int new_cols = old_cols + 2 * pad_size;

	// Create a new matrix initialized entirely with the padding value
	vector<std::vector<int>> padded(new_rows, vector<int>(new_cols, pad_value));

	// Copy original data into the center region
	for (int i = 0; i < old_rows; ++i) 
	{
		std::copy(matrix[i].begin(), matrix[i].end(), padded[i + pad_size].begin() + pad_size);
	}

	return padded;
}

vector<double> CNN_1DConvolutionOperation(const vector<double>& u, const vector<double>& v) 
{
	int m = u.size();
	int n = v.size();
    
	// The output vector size is always M + N - 1
	vector<double> w(m + n - 1, 0.0);

	// Perform the convolution operation
	for (int i = 0; i < m; ++i) 
	{
		for (int j = 0; j < n; ++j) 
		{
			w[i + j] += u[i] * v[j];
		}
	}

	return w;
}

// Function to apply kernel convolution with stride
vector<vector<double>> CNN_2DConvolutionOperation(const vector<vector<int>>& input, const vector<vector<double>>& kernel, int stride) 
{
	int inputH = input.size();
	int inputW = input[0].size();
	int kernelH = kernel.size();
	int kernelW = kernel[0].size();

	// Calculate output dimensions
	int outputH = (inputH - kernelH) / stride + 1;
	int outputW = (inputW - kernelW) / stride + 1;

	// Initialize output matrix with zeros
	vector<vector<double>> output(outputH, vector<double>(outputW, 0.0));

	// Perform convolution
	for (int i = 0; i < outputH; ++i) 
	{
		for (int j = 0; j < outputW; ++j) 
		{
			double sum = 0.0;
		    
			// Map output position back to input position using stride
			int startRow = i * stride;
			int startCol = j * stride;

			// Multiply kernel with the matching input window
			for (int kh = 0; kh < kernelH; ++kh) 
			{
				for (int kw = 0; kw < kernelW; ++kw) 
				{
					sum += input[startRow + kh][startCol + kw] * kernel[kh][kw];
				}
			}
			output[i][j] = sum;
		}
	}

	return output;
}

void CNN_2DConvolutionOperation(vector<vector<int>>& input, vector<vector<double>>& kernel)
{
	int row_input = input.size();
	int col_input = input[0].size();
	int row_kernel = kernel.size();
	int col_kernel = kernel[0].size();
	if (row_kernel != col_kernel) 
	{
		throw std::invalid_argument("Kernel has to be square.");
	}
	if (row_kernel > row_input) 
	{
		throw std::invalid_argument("Kernel dimension has to be smaller than input matrix dimension.");
	}
	int n = row_input - row_kernel + 1;
	int m = col_input - col_kernel + 1;
	vector<vector<double>> outputmap(n, vector<double>(m,0.0));
   	
	vector<vector<double>> hadamardproduct(row_kernel, vector<double>(col_kernel,0.0));
	for (int i = 0; i < n; ++i) 
	{
		
		cout << "\ni: " << i << endl ;
		for (int j = 0; j < m; ++j) 
		{
			int R = i + row_kernel;
			int C = j + col_kernel;
			
			vector<vector<double>> subMatrix;
			subMatrix.assign(row_kernel, vector<double>(col_kernel, 0.0));
			
			int i_submatrix = 0;
			for (int i_index = i; i_index < R; ++i_index) 
			{
				int j_submatrix = 0;
				for (int j_index = j; j_index < C; ++j_index) 
				{
					subMatrix[i_submatrix][j_submatrix] = input[i_index][j_index];

					j_submatrix += 1;
				}
				i_submatrix += 1;
			}
			cout <<"\nSubmatrix at j: " << j << endl;
			printMatrix(subMatrix);
			for (int i1 = 0; i1 < row_kernel; ++i1) 
			{
				for (int j1 = 0; j1 < col_kernel; ++j1) 
				{
					hadamardproduct[i1][j1] = subMatrix[i1][j1] * kernel[i1][j1];
					outputmap[i][j] += hadamardproduct[i1][j1];
				}
			}
			cout <<"\nHadamard product: " << endl;
			printMatrix(hadamardproduct);
		}
		
	}
	cout <<"\nOutput map: " << endl;
	printMatrix(outputmap);
	
}

vector<vector<int>> CNN_2DmaxPooling(const vector<vector<int>>& input, int kernelSize, int stride)
{
	int inputHeight = input.size();
	int inputWidth = input[0].size();
    
	// Calculate output dimensions
	int outputHeight = (inputHeight - kernelSize) / stride + 1;
	int outputWidth = (inputWidth - kernelSize) / stride + 1;
    
	vector<vector<int>> output(outputHeight, vector<int>(outputWidth));
    
	for (int y = 0; y < outputHeight; ++y) 
	{
		for (int x = 0; x < outputWidth; ++x) 
		{
			int maxVal = INT_MIN;
            
			// Define the boundaries of the pooling window
			int startY = y * stride;
			int startX = x * stride;
			int endY = startY + kernelSize;
			int endX = startX + kernelSize;
            
			// Find the maximum value within the kernel window
			for (int j = startY; j < endY; ++j) 
			{
				for (int i = startX; i < endX; ++i) 
				{
					if (input[j][i] > maxVal) 
					{
						maxVal = input[j][i];
					}
				}
			}
			output[y][x] = maxVal;
		}
	}
	return output;
}

vector<vector<double>> CNN_2DaveragePooling(const vector<vector<int>>& input, int poolSize, int stride) 
{
	if (input.empty() || input[0].empty()) return {};

	int inputHeight = input.size();
	int inputWidth = input[0].size();

	// Calculate output dimensions
	int outputHeight = (inputHeight - poolSize) / stride + 1;
	int outputWidth = (inputWidth - poolSize) / stride + 1;

	vector<vector<double>> pooledOutput(outputHeight, std::vector<double>(outputWidth, 0.0));

	for (int i = 0; i < outputHeight; ++i) 
	{
		for (int j = 0; j < outputWidth; ++j) 
		{
			double sum = 0.0;
			int count = 0;

			// Define the window boundaries
			int startRow = i * stride;
			int startCol = j * stride;
			int endRow = std::min(startRow + poolSize, inputHeight);
			int endCol = std::min(startCol + poolSize, inputWidth);

			// Calculate average in the current window
			for (int r = startRow; r < endRow; ++r) 
			{
				for (int c = startCol; c < endCol; ++c) 
				{
					sum += input[r][c];
					count++;
				}
			}

			pooledOutput[i][j] = (count > 0) ? (sum / count) : 0.0;
		}
	}

	return pooledOutput;
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
double He_initialization(double n_in) 
{
	double mu = 0;
	double sigma = sqrt(divisiond(2,n_in));

	static std::random_device rd;
	static std::mt19937 gen(rd());
	std::normal_distribution<double> distribution(mu, sigma);
	return distribution(gen);
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

	vector<vector<double>> weights(n_input, vector<double>(n_classification));
        vector<double> bias_output(n_classification, 0.0); // biases set to 0
	vector<double> predicted_output(n_classification, 0.0);
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	vector<double> delta(n_classification, 0.0);

	double learning_rate = 0.5;
	int epochs = 1000;

	for (int i = 0; i < n_input; ++i)
	{
		for(int j = 0; j < n_classification; ++j)
		{
			weights[i][j] = random_double(-1,1);	
		}	
	}
	for(int j = 0; j < n_classification; ++j)
	{
		bias_output[j] = random_double(-1,1);
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

			 Backpropagation with steepest descent / gradient descent to update weights

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
// May 12th, 2026 Sweden Sexy
// Homework: weights matrix size has to be fixed -> input * next layer neuron

	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_data = input_matrix.size();
	int n_input = input_matrix[0].size();

	vector<vector<double>> weights(n_input, vector<double>(n_hiddenlayer));
	vector<vector<double>> hiddenweights(n_hiddenlayer, vector<double>(n_classification));
        vector<double> bias_output(n_classification, 0.0); // biases set to 0
	vector<double> bias_hiddenlayer1(n_hiddenlayer, 0.0); // biases set to 0
	vector<double> predicted_output(n_classification, 0.0);
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	vector<double> delta(n_classification, 0.0);

	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
	vector<double> delta_hiddenweights(n_hiddenlayer, 0.0);

	double learning_rate = 0.5;
	int epochs = 1000;

	for (int i = 0; i < n_input; ++i)
	{
		for(int j = 0; j < n_hiddenlayer; ++j)
		{
			weights[i][j] = random_double(-1,1);	
		}	
	}
	for(int j = 0; j < n_hiddenlayer; ++j)
	{
		bias_hiddenlayer1[j] = random_double(-1,1);
	}
	for(int k = 0; k < n_classification; ++k)
	{
		bias_output[k] = random_double(-1,1);
	}
	for (int j = 0; j < n_hiddenlayer; ++j)
	{
		for(int k = 0; k < n_classification; ++k)
		{
			hiddenweights[j][k] = random_double(-1,1);	
		}	
	}
	
	cout <<"\nTrain data: " << n_data << endl;
	cout <<"\nInitial weights from input layer to hidden layer:" << endl;
	printMatrix(weights); 
	cout <<"\nInitial weights from hidden layer to output layer:" << endl;
	printMatrix(hiddenweights); 
	cout <<"\nInitial bias in the hidden layer:" << endl;
	printVector(bias_hiddenlayer1);
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
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				hiddenlayer[j] = 0;
			}
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				// Linear transformation: Wx + b
				for (int i = 0; i < n_input; ++i) 
				{
					hiddenlayer[j] += weights[i][j] * input[i];
				}
				// Add bias vector
				hiddenlayer[j] += bias_hiddenlayer1[j];
				// Activation function
				hiddenlayer[j] = SigmoidActivationfunction(hiddenlayer[j]);
			}
			for (int k = 0; k < n_classification; ++k) 
			{
				predicted_output[k] = 0;
			}
			for (int k = 0; k < n_classification; ++k) 
			{
				// Linear transformation: Wx + b
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					predicted_output[k] += hiddenweights[j][k] * hiddenlayer[j];
				}
				// Add bias vector
				predicted_output[k] += bias_output[k];
			}
			// Softmax activation function
			vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(predicted_output);
			
			for (int k = 0; k < n_classification; ++k) 
			{
				// Save the activation function
				predicted_output[k] = predicted_output_softmax[k];
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

			 Backpropagation to update weights

			*/
			for (int k = 0; k < n_classification; ++k) 
			{
				delta[k] = SigmoidDerivative(predicted_output[k])*(error[k]);
			}
			//cout <<"\nDelta output vector:" << endl;
			//printVector(delta);	
			for (int k = 0; k < n_classification; ++k) 
			{
				bias_output[k] += (learning_rate*delta[k]);
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					hiddenweights[j][k] +=  learning_rate*delta[k]*hiddenlayer[j];
				}
			}	
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				double sum_output = 0;
				for (int k = 0; k < n_classification; ++k) 
				{
					sum_output += delta[k]*hiddenweights[j][k];
				}
				delta_hiddenweights[j] = SigmoidDerivative(hiddenlayer[j])*(sum_output);
			}
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				bias_hiddenlayer1[j] += (learning_rate*delta_hiddenweights[j]);
				for (int i = 0; i < n_input; ++i) 
				{
					weights[i][j] +=  learning_rate*delta_hiddenweights[j]*input[i];
				}
			}	
				//cout <<"\nNew W0 matrix (from input to hidden layer):" << endl;
				//printMatrix(weights);	
				
				//cout <<"\nNew bias in the output layer:" << endl;
				//printVector(bias_output);
				
		}	
		cout <<"\nEpoch loss:" << epoch_loss << endl;
	}
	cout <<"\nModel W0 matrix (from input to hidden layer):" << endl;
	printMatrix(weights);	
	cout <<"\nModel W1 matrix (from hidden layer to output):" << endl;
	printMatrix(hiddenweights);	
	saveMatrixdouble(weights,"modelweights.txt");
	saveMatrixdouble(hiddenweights,"modelhiddenweights.txt");
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
	int correct = 0, correct_setosa = 0, correct_versicolor = 0, correct_virginica = 0;
	int wrong = 0;	
	vector<vector<double>> weights = loadMatrixFromFile("modelweights.txt");
	vector<vector<double>> hiddenweights = loadMatrixFromFile("modelhiddenweights.txt");
        vector<double> bias_hiddenlayer1 = loadVectorFromFile("bias_hiddenlayer1.txt");
        vector<double> bias_output = loadVectorFromFile("bias_output.txt");
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);        
        vector<double> error(n_classification, 0.0);
	vector<double> actual_output(n_classification,0.0);
	vector<vector<double>> mat_correct;	
	vector<double> vec_error;
	vector<double> vec_CCEloss;
	vector<string> vec_correctorwrong;

	cout <<"\nTesting data: " << n_data << endl;
	cout <<"\nModel weights from input to hidden layer:" << endl;
	printMatrix(weights);
	cout <<"\nModel weights from hidden layer to output layer:" << endl;
	printMatrix(hiddenweights);
	

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


		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			hiddenlayer[j] = 0;
		}
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			// Linear transformation: Wx + b
			for (int i = 0; i < n_input; ++i) 
			{
				hiddenlayer[j] += weights[i][j] * input[i];
			}
			hiddenlayer[j] += bias_hiddenlayer1[j];
			// Activation function
			hiddenlayer[j] = SigmoidActivationfunction(hiddenlayer[j]);
		}
		
		//cout <<"\nHidden layer:" << endl;
		//printVector(hiddenlayer);	
		for (int k = 0; k < n_classification; ++k) 
		{
			predicted_output[k] = 0;
		}
		for (int k = 0; k < n_classification; ++k) 
		{
			// Linear transformation: Wx + b
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				predicted_output[k] += hiddenweights[j][k] * hiddenlayer[j];
			}
			predicted_output[k] += bias_output[k];
		}
		
		vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(predicted_output);
		for (int k = 0; k < n_classification; ++k) 
		{
			// Activation function
			predicted_output[k] = predicted_output_softmax[k];
			error[k] = actual_output[k] - predicted_output[k];
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
	

	mat_correct.push_back({double(correct), double(correct_setosa), double(correct_versicolor), double(correct_virginica), double(1)});
	cout <<"Correct matrix: " << endl;
	cout <<setw(23) <<"Total Correct" <<setw(23) << "Setosa" <<setw(23) << "Versicolor" << setw(23) <<"Virginica" << setw(23) << "Model" << endl;
	printMatrix(mat_correct);
	
}


// Utility function to read the dataset from a CSV file
vector<vector<double>> read_dataset_iris(const string &filename, vector<int> &labels, vector<string> &label_names) 
{
	std::ifstream file(filename);
	if (!file.is_open()) 
	{
		throw std::runtime_error("Unable to open file: " + filename);
	}

	vector<vector<double>> dataset;
	string line;
	unordered_map<string, int> label_map;
	int label_counter = 0;

	std::getline(file, line); // Skip header

	while (getline(file, line)) 
	{
		stringstream ss(line);
		vector<double> data_point;
		string value;

		// Skip the ID column
		std::getline(ss, value, ',');

		// Read the 4 feature values
		for (int i = 0; i < 4; ++i) 
		{
			if (!std::getline(ss, value, ',')) 
			{
				throw std::runtime_error("Error reading CSV: unexpected end of line");
			}
			try 
			{
				data_point.push_back(std::stod(value));
			} catch (const std::exception& e) 
			{
				throw std::runtime_error("Error converting to double: " + value);
			}
		}

		// Read the species label
		if (!std::getline(ss, value)) 
		{
			throw std::runtime_error("Error reading CSV: missing species label");
		}
		// this will give integer labeling 
		// 0 corresponding to Setosa species, 1 corresponding to Versicolor species, and 2 corresponding to Virginica species.
		if (label_map.find(value) == label_map.end()) 
		{
			label_map[value] = label_counter++;
			label_names.push_back(value);
		}
		labels.push_back(label_map[value]);

		dataset.push_back(data_point);
	} // end while

	file.close();
	return dataset;
}


	
// Constructor initializing the array with random numbers
FNN_NoHiddenLayer_Iris::FNN_NoHiddenLayer_Iris() 
{
	int n_classification = 3;
	int n_input = 4;
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_classification; ++j) 
		{
			weights[i][j] = random_double(-0.5, 0.5);
		}
	}
	for (int j = 0; j < n_classification; ++j) 
	{
		bias[j] = random_double(-0.5, 0.5);
        }
	learning_rate = 0.01;
}

vector<double> FNN_NoHiddenLayer_Iris::predict(const vector<double> &inputs) 
{
	int n_classification = 3;
	int n_input = 4;
	vector<double> outputs(n_classification, 0.0);
	for (int j = 0; j < n_classification; ++j) 
	{
		for (int i = 0; i < n_input; ++i) 
		{
			outputs[j] += inputs[i] * weights[i][j];
		}
		outputs[j] += bias[j];
		outputs[j] = SigmoidActivationfunction(outputs[j]);
	}
	return outputs;
}

void FNN_NoHiddenLayer_Iris::train(vector<vector<double>> &X, vector<int> &y, int epochs) 
{
	cout << BLUE << BOLD << "\nTraining Progress:\n" << RESET << endl;

	int bar_width = 50;
	int n_classification = 3;
	int n_input = 4;
	vector<vector<vector<double>>> Matrix3D_weights;
	vector<vector<double>> mat_w(n_input, vector<double>(n_classification, 0));
	vector<vector<double>> mat_bias;
	vector<double> vec_bias(n_classification,0.0);
	for (int epoch = 0; epoch < epochs; ++epoch) 
	{
		double epoch_loss = 0.0;
		for (size_t i = 0; i < X.size(); ++i) 
		{
		vector<double> outputs = predict(X[i]);
		vector<double> target(n_classification, 0.0);
		target[y[i]] = 1.0;
				
 		std::vector<double> errors(n_classification, 0.0);
			for (int j = 0; j < n_classification; ++j) 
			{
				errors[j] = target[j] - outputs[j];
				epoch_loss += errors[j] * errors[j];
			}

			for (int j = 0; j < n_classification; ++j) 
			{
				for (int k = 0; k < n_input; ++k) 
				{
					weights[k][j] += learning_rate * errors[j] * SigmoidDerivative(outputs[j]) * X[i][k];
					mat_w[k][j] = weights[k][j];
				}
			bias[j] += learning_rate * errors[j] * SigmoidDerivative(outputs[j]);
			vec_bias[j] = bias[j];
			}
		}

		if ((epoch + 1) % (epochs / 100) == 0 || epoch == epochs - 1) 
		{
			float progress = static_cast<float>(epoch + 1) / epochs;
			int pos = static_cast<int>(bar_width * progress);

			std::cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) 
				{
					std::cout << "=";
				}
				else if (i == pos)
				{
					std::cout << ">";
				}
				else
				{
					std::cout << " ";
				}
			}
			std::cout << "] " << int(progress * 100.0) << "% ";
			std::cout << "Epoch " << epoch + 1 << "/" << epochs << " - Loss: " << std::fixed << std::setprecision(4) << epoch_loss / X.size() << "\r";
			std::cout.flush();
		}
		Matrix3D_weights.push_back(mat_w);
		mat_bias.push_back(vec_bias);
	} 
	save3DMatrixdouble(Matrix3D_weights,"mat3d_weights.txt");
	saveMatrixdouble(mat_bias,"matrixbias.txt");
	cout << endl;
}

void FNN_NoHiddenLayer_Iris::save_model(const std::string &filename) 
{
	int n_classification = 3;
	int n_input = 4;
	std::ofstream file(filename);
	file << std::fixed << std::setprecision(6);
	vector<vector<double>> mat_w_dummy(n_input, vector<double>(n_classification, 0));
	vector<double> vec_bias(n_classification,0.0);
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_classification; ++j) 
		{
			file << weights[i][j] << " ";
			mat_w_dummy[i][j] = weights[i][j];
		}
	}
	for (int j = 0; j < n_classification; ++j) 
	{
		file << bias[j] << " ";
		vec_bias[j]= bias[j];
	}
	file.close();
	cout << "\nThe model weights:" << endl;
	printMatrix(mat_w_dummy);
	cout << "\nThe model bias:" << endl;
	printVector(vec_bias);
}

void FNN_NoHiddenLayer_Iris::load_model(const std::string& filename) 
{
	int n_classification = 3;
	int n_input = 4;
	std::ifstream file(filename);
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_classification; ++j) 
		{
			file >> weights[i][j];
		}
	}
	for (int j = 0; j < n_classification; ++j) 
	{
		file >> bias[j];
	}
	file.close();
}

// Constructor initializing the array with random numbers
FNN_NoHiddenLayer_Iris_ConjugateGradient::FNN_NoHiddenLayer_Iris_ConjugateGradient() 
{
	int n_classification = 3;
	int n_input = 4;
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_classification; ++j) 
		{
			weights[i][j] = random_double(-0.5, 0.5);
		}
	}
	for (int j = 0; j < n_classification; ++j) 
	{
		bias[j] = random_double(-0.5, 0.5);
        }
	learning_rate = 0.51;
}

vector<double> FNN_NoHiddenLayer_Iris_ConjugateGradient::predict(const vector<double> &inputs) 
{
	int n_classification = 3;
	int n_input = 4;
	vector<double> outputs(n_classification, 0.0);
	for (int j = 0; j < n_classification; ++j) 
	{
		for (int i = 0; i < n_input; ++i) 
		{
			outputs[j] += inputs[i] * weights[i][j];
		}
		outputs[j] += bias[j];
		outputs[j] = SigmoidActivationfunction(outputs[j]);
	}
	return outputs;
}


void FNN_NoHiddenLayer_Iris_ConjugateGradient::train(vector<vector<double>> &X, vector<int> &y, int epochs) 
{
	// May 27th, 2026: Computation with beta cannot converge, so we change beta into fix 0.01 and it converges well.
	// the learning rate is adjusted at every step with Strong Wolfe condition modified by us
	cout << BLUE << BOLD << "\nTraining Progress:\n" << RESET << endl;

	int bar_width = 50;
	int n_classification = 3;
	int n_input = 4;
	vector<vector<vector<double>>> Matrix3D_weights;
	vector<vector<double>> mat_w(n_input, vector<double>(n_classification, 0));
	vector<vector<double>> mat_bias;
	vector<double> vec_bias(n_classification,0.0);
	double beta = 0.01;
	double new_learning_rate = learning_rate;
	vector<vector<double>> mat_target;
	vector<vector<double>> mat_gk;
	vector<double> learningrate;
	//vector<double> beta;
	vector<double> norm_gk;
	vector<double> gk;
	vector<double> gk_previous;
	vector<double> dk;
	vector<double> dk_previous;
	int reset_iter = 0;
	int n_dk= n_input*n_classification + n_classification;

	for (int epoch = 0; epoch < epochs; ++epoch) 
	{
		double epoch_loss = 0.0;
		for (size_t i = 0; i < X.size(); ++i) 
		{
			vector<double> outputs = predict(X[i]);
			vector<double> target(n_classification, 0.0);
			target[y[i]] = 1.0;
			mat_target.push_back(outputs);
			mat_target.push_back(target);
			
	 		vector<double> errors(n_classification, 0.0);
			for (int j = 0; j < n_classification; ++j) 
			{
				errors[j] = target[j] - outputs[j];
				epoch_loss += errors[j] * errors[j];
			}
			double errors_currentweights =norm(errors);	
			
			
			if(epoch ==0 && i ==0)
			{
				// compute g_{0}/ create vector g_{0} it is the default steepest descent or gradient descent algorithm
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						gk.push_back(errors[j] *  SigmoidDerivative(outputs[j]) * X[i][k]);
					}
				}
				for (int j = 0; j < n_classification; ++j) 
				{
					gk.push_back(errors[j] * SigmoidDerivative(outputs[j]));
				}
				// d_{0} = -g_{0}
				for (int j = 0; j < n_dk; ++j) 
				{
					dk.push_back(gk[j]);
				}	
				// Compute the first weights w_{1} = w_{0} + learning_rate * g_{0}
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						weights[k][j] += new_learning_rate  * ( gk[k + n_input*j]);
						mat_w[k][j] = weights[k][j];
					}
				bias[j] += new_learning_rate * ( gk[n_input*n_classification + j] );
				vec_bias[j] = bias[j];
				}
				i = i+1;
			}

			if(epoch !=0 && i !=0)
			{
			gk_previous.clear();
			dk_previous.clear();
			for (int j = 0; j < n_dk; ++j) 
			{
				gk_previous.push_back(gk[j]);
				dk_previous.push_back(dk[j]);
			}	
	
			gk.clear();
			dk.clear();
			// compute g_{k}/ starts from g_{1}
			for (int j = 0; j < n_classification; ++j) 
			{
				for (int k = 0; k < n_input; ++k) 
				{
					gk.push_back(errors[j] *  SigmoidDerivative(outputs[j]) * X[i][k]);
				}
			}
			for (int j = 0; j < n_classification; ++j) 
			{
				gk.push_back(errors[j] * SigmoidDerivative(outputs[j]));
			}
			vector<double> gk_diff = subtract(gk,gk_previous);
					
			// Compute beta Polak-Ribiere, Fletcher-Reeves, OFR
			//double mu = 1.1;
			//double beta_PR = divisiond(dot(gk,gk_diff), dot(gk_previous,gk_previous));
			//double beta_FR = divisiond(dot(gk,gk), dot(gk_previous,gk_previous));		
			//double beta_OFR = divisiond(dot(gk,gk), mu*abs(dot(gk,dk_previous)) + norm(gk_previous)*norm(dk_previous));		
			//if (norm(gk_previous) ==0 || norm(dk_previous) ==0)
			//{
			//	beta_PR = 0;
			//}
			//beta.push_back(beta_PR);
			// Compute d_{k+1} = -g_{k+1} + beta * d_{k}
			if(beta >= norm(gk)/norm(gk_previous))
			{
				beta = norm(gk)/norm(gk_previous); // GLANZ beta
			}	
			else if(beta < norm(gk)/norm(gk_previous))
			{
				beta = 0.01; // GLANZ beta
			}
			for (int j = 0; j < n_dk; ++j) 
			{
				//dk.push_back(gk[j]); // the steepest descent algorithm
				dk.push_back(gk[j] + beta*dk_previous[j]); // the modified conjugate gradient algorithm
			}

						
			/*
			 Compute the learning rate with Strong Wolfe line search, iterative method
			*/
			// Strong wolfe line search
			double delta = 0.1, sigma = 0.9;
			double low = 0, high = 1e10;
			int sw_iter = 10;
			for (int i = 0; i<sw_iter ; ++i)
			{
				vector<double> gk_strongwolfe_newweights;
				// Compute new weights based on current learning rate / step size w_{k} + learning_rate_{k} * d_{k}
				
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						weights[k][j] += new_learning_rate * dk[k + n_input*j];	
					}
					bias[j] += new_learning_rate *  dk[n_input*n_classification + j] ;
				}

				// Compute the new error based on the new weights
				vector<double> new_outputs = predict(X[i]);
	 			vector<double> new_errors(n_classification, 0.0);
				for (int j = 0; j < n_classification; ++j) 
				{
					new_errors[j] = target[j] - new_outputs[j]; 
				} 
				double errors_newweights = norm(new_errors);
				// Neutralize weights and bias
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						weights[k][j] -= new_learning_rate * dk[k + n_input*j];	
					}
					bias[j] -= new_learning_rate *  dk[n_input*n_classification + j] ;
				}
			
				// Creating vector g_{k} with push_back, for the strong wolfe conditions iteration
				// for weights
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						gk_strongwolfe_newweights.push_back(new_errors[j] *  SigmoidDerivative(outputs[j]) * X[i][k]);
					}
				}
				for (int j = 0; j < n_classification; ++j) 
				{
					gk_strongwolfe_newweights.push_back(new_errors[j] * SigmoidDerivative(outputs[j]));
				}

				double strongwolfeeq1_lhs = errors_newweights;
				double strongwolfeeq1_rhs = errors_currentweights + delta*new_learning_rate*dot(gk,dk_previous);
					
				double strongwolfeeq2_lhs = abs(dot(gk_strongwolfe_newweights,dk_previous));
				double strongwolfeeq2_rhs = sigma*abs(dot(gk,dk_previous));

				if (strongwolfeeq2_lhs <= strongwolfeeq2_rhs && strongwolfeeq1_lhs <= strongwolfeeq1_rhs)
				{
					new_learning_rate = new_learning_rate;
					i = sw_iter;
				}
						
				// Check Armijo Condition
				if (strongwolfeeq1_lhs > strongwolfeeq1_rhs ) 
				{
					high = new_learning_rate;
					new_learning_rate = 0.5*(low + high);
				} 
				// Check Curvature Condition
				else if (strongwolfeeq2_lhs < strongwolfeeq2_rhs ) 
				{
					/*low = new_learning_rate;
					if (high > 1e9) 
					{
						new_learning_rate *= 2.0;
					}
					else if (high <= 1e9) 
					{
						new_learning_rate = 0.5*(low + high) ;
					}*/
					new_learning_rate = new_learning_rate;
					i = sw_iter;
				} 

				if (new_learning_rate < 1e-5) 
				{
					// for iris dataset the best is 1e-1
					new_learning_rate =  1e-1; // return small number if the learning rate becoming too small.
					i = sw_iter; // prevent infinite loop
				}
					
			}// End of strong wolfe conditions iteration for inexact line search			
			learningrate.push_back(new_learning_rate);


			if ((reset_iter) % (n_dk) != 0) 
			{
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						weights[k][j] += new_learning_rate * ( dk[k + n_input*j]);
						mat_w[k][j] = weights[k][j];
					}
					bias[j] += new_learning_rate * ( dk[n_input*n_classification + j] );
					vec_bias[j] = bias[j];
				}	
			}
			if ((reset_iter) % (n_dk) == 0) // Periodically reset d_{k} = -g_{k}
			{
				for (int j = 0; j < n_classification; ++j) 
				{
					for (int k = 0; k < n_input; ++k) 
					{
						weights[k][j] += new_learning_rate  * ( gk[k + n_input*j]);
						mat_w[k][j] = weights[k][j];
					}
				bias[j] += new_learning_rate * ( gk[n_input*n_classification + j] );
				vec_bias[j] = bias[j];
				}
			}
			mat_gk.push_back(gk);
			norm_gk.push_back(norm(gk));

			//if (norm(gk) ==0)
			{
				//break;
			}	
			}
			reset_iter += 1;	
		}// end iterating all training samples

		if ((epoch + 1) % (epochs / 100) == 0 || epoch == epochs - 1) 
		{
			float progress = static_cast<float>(epoch + 1) / epochs;
			int pos = static_cast<int>(bar_width * progress);

			std::cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) 
				{
					std::cout << "=";
				}
				else if (i == pos)
				{
					std::cout << ">";
				}
				else
				{
					std::cout << " ";
				}
			}
			std::cout << "] " << int(progress * 100.0) << "% ";
			std::cout << "Epoch " << epoch + 1 << "/" << epochs << " - Loss: " << std::fixed << std::setprecision(4) << epoch_loss / X.size() << "\r";
			std::cout.flush();
		}
		Matrix3D_weights.push_back(mat_w);
		mat_bias.push_back(vec_bias);
	} 
	//save3DMatrixdouble(Matrix3D_weights,"mat3d_weights.txt");
	//saveMatrixdouble(mat_bias,"matrixbias.txt");
	saveMatrixdouble(mat_target,"matrixoutputstarget.txt");
	saveMatrixdouble(mat_gk,"matrixgk.txt");
	//saveVectordouble(norm_gk,"normgk.txt");
	//saveVectordouble(beta,"beta.txt");
	saveVectordouble(learningrate,"strongwolfelearningrate.txt");
	cout << endl;
}

void FNN_NoHiddenLayer_Iris_ConjugateGradient::save_model(const std::string &filename) 
{
	int n_classification = 3;
	int n_input = 4;
	std::ofstream file(filename);
	file << std::fixed << std::setprecision(6);
	vector<vector<double>> mat_w_dummy(n_input, vector<double>(n_classification, 0));
	vector<double> vec_bias(n_classification,0.0);
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_classification; ++j) 
		{
			file << weights[i][j] << " ";
			mat_w_dummy[i][j] = weights[i][j];
		}
	}
	for (int j = 0; j < n_classification; ++j) 
	{
		file << bias[j] << " ";
		vec_bias[j]= bias[j];
	}
	file.close();
	cout << "\nThe model weights:" << endl;
	printMatrix(mat_w_dummy);
	cout << "\nThe model bias:" << endl;
	printVector(vec_bias);
}

void FNN_NoHiddenLayer_Iris_ConjugateGradient::load_model(const std::string& filename) 
{
	int n_classification = 3;
	int n_input = 4;
	std::ifstream file(filename);
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_classification; ++j) 
		{
			file >> weights[i][j];
		}
	}
	for (int j = 0; j < n_classification; ++j) 
	{
		file >> bias[j];
	}
	file.close();
}

/*

Iris class with 1 hidden layer

*/	

// Constructor initializing the array with random numbers
FNN_1HiddenLayer_Iris::FNN_1HiddenLayer_Iris() 
{
	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_input = 4;
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			weights[i][j] = random_double(-0.5, 0.5);
		}
	}
	for (int j = 0; j < n_hiddenlayer; ++j) 
	{
		for (int k = 0; k < n_classification; ++k) 
		{
			hiddenweights[j][k] = random_double(-0.5, 0.5);
		}
	}
	for (int j = 0; j < n_hiddenlayer; ++j) 
	{
		biashidden[j] = random_double(-0.5, 0.5);
        }
	for (int k = 0; k < n_classification; ++k) 
	{
		biasoutput[k] = random_double(-0.5, 0.5);
        }
	learning_rate = 1;
}
	

vector<double> FNN_1HiddenLayer_Iris::predict(const vector<double> &inputs) 
{
	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_input = 4;
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
	vector<double> outputs(n_classification, 0.0);
	for (int j = 0; j < n_hiddenlayer; ++j) 
	{
		for (int i = 0; i < n_input; ++i) 
		{
			hiddenlayer[j] += inputs[i] * weights[i][j];
		}
		hiddenlayer[j] += biashidden[j];
		hiddenlayer[j] = SigmoidActivationfunction(hiddenlayer[j]);
		//hiddenlayer[j] =  TanhActivationfunction(hiddenlayer[j]);
	}
	for (int k = 0; k < n_classification; ++k) 
	{
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			outputs[k] += hiddenlayer[j] * hiddenweights[j][k];
		}
		outputs[k] += biasoutput[k];
	}
	vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(outputs);
	return predicted_output_softmax;
}

void FNN_1HiddenLayer_Iris::train(vector<vector<double>> &X, vector<int> &y, int epochs) 
{
	cout << BLUE << BOLD << "\nTraining Progress:\n" << RESET << endl;

	int bar_width = 50;
	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_input = 4;

	//vector<vector<vector<double>>> Matrix3D_weights;
	//vector<vector<vector<double>>> Matrix3D_hiddenweights;
	
	//vector<vector<double>> mat_w(n_input, vector<double>(n_hiddenlayer, 0.0));
	//vector<vector<double>> mat_whidden(n_hiddenlayer, vector<double>(n_classification, 0.0));
	//vector<vector<double>> mat_biashidden;
	//vector<vector<double>> mat_biasoutput;
	//vector<vector<double>> mat_target;
	//vector<double> vec_biashidden(n_hiddenlayer,0.0);
	//vector<double> vec_biasoutput(n_classification,0.0);
	
	for (int epoch = 0; epoch < epochs; ++epoch) 
	{
		double epoch_loss = 0.0;
		for (size_t iter = 0; iter < X.size(); ++iter) 
		{
			vector<vector<double>> mat_jacobian(n_classification, vector<double>(n_classification, 0.0)); // Softmax derivative
			vector<double> delta_outputlayer(n_classification,0.0);
			vector<double> outputs = predict(X[iter]);
			vector<double> target(n_classification, 0.0);
			target[y[iter]] = 1.0;
			//mat_target.push_back(outputs);	
			//mat_target.push_back(target);		
	 		std::vector<double> errors(n_classification, 0.0);
			for (int k = 0; k < n_classification; ++k) 
			{
				//errors[k] = target[k] - outputs[k];
				errors[k] = -target[k] * log(outputs[k]); // Categorical Cross-Entropy (CCE) Loss because we are using softmax activation function at the output layer 
				epoch_loss += errors[k];
			}
			vector<double> hiddenlayer(n_hiddenlayer, 0.0);
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				for (int i = 0; i < n_input; ++i) 
				{
					hiddenlayer[j] += X[iter][i] * weights[i][j];
				}
				hiddenlayer[j] += biashidden[j];
				hiddenlayer[j] = SigmoidActivationfunction(hiddenlayer[j]);
				//hiddenlayer[j] = TanhActivationfunction(hiddenlayer[j]);
			}
			for (int k_row = 0; k_row < n_classification; ++k_row) 
			{
				for (int k_col = 0; k_col < n_classification; ++k_col) 
				{
					if (k_row == k_col)
					{
						mat_jacobian[k_row][k_col] = outputs[k_row]*(1 - outputs[k_col]);
					}
					else if (k_row != k_col)
					{
						mat_jacobian[k_row][k_col] = outputs[k_row]*(- outputs[k_col]);
					}
				}
			}
			for (int k_row = 0; k_row < n_classification; ++k_row) 
			{
				for (int k_col = 0; k_col < n_classification; ++k_col) 
				{
					delta_outputlayer[k_row] += mat_jacobian[k_row][k_col]*errors[k_col];
				}
			}
			// this delta code makes the hiddenweights to grow big and make the epoch loss > 1
			/*for (int k = 0; k < n_classification; ++k) 
			{
				delta_outputlayer[k] = outputs[k] - target[k];
			}*/
			for (int k = 0; k < n_classification; ++k) 
			{
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					hiddenweights[j][k] += learning_rate * delta_outputlayer[k] * hiddenlayer[j]; 
					//mat_whidden[j][k] = hiddenweights[j][k];
				}
				biasoutput[k] += learning_rate * delta_outputlayer[k];
				//vec_biasoutput[k] = biasoutput[k];
			}

			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				double sum_delta = 0;			
				for (int k = 0; k < n_classification; ++k) 
				{
					sum_delta += hiddenweights[j][k] *  delta_outputlayer[k] ; 
				}
				for (int i = 0; i < n_input; ++i) 
				{
					weights[i][j] += learning_rate * sum_delta * SigmoidDerivative(hiddenlayer[j]) * X[iter][i];					
					//weights[i][j] += learning_rate * sum_delta * TanhDerivative(hiddenlayer[j]) * X[iter][i];
					//mat_w[i][j] = weights[i][j];
				}
			biashidden[j] += learning_rate *  sum_delta * SigmoidDerivative(hiddenlayer[j]) ;
			//biashidden[j] += learning_rate *  sum_delta * TanhDerivative(hiddenlayer[j]) ;
			//vec_biashidden[j] = biashidden[j];
			}
		}

		if ((epoch + 1) % (epochs / 100) == 0 || epoch == epochs - 1) 
		{
			float progress = static_cast<float>(epoch + 1) / epochs;
			int pos = static_cast<int>(bar_width * progress);

			std::cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) 
				{
					std::cout << "=";
				}
				else if (i == pos)
				{
					std::cout << ">";
				}
				else
				{
					std::cout << " ";
				}
			}
			std::cout << "] " << int(progress * 100.0) << "% ";
			std::cout << "Epoch " << epoch + 1 << "/" << epochs << " - Loss: " << std::fixed << std::setprecision(4) << epoch_loss / X.size() << "\r";
			std::cout.flush();
		}
		//Matrix3D_weights.push_back(mat_w);
		//Matrix3D_hiddenweights.push_back(mat_whidden);
		//mat_biashidden.push_back(vec_biashidden);
		//mat_biasoutput.push_back(vec_biasoutput);
		
	} 
	//save3DMatrixdouble(Matrix3D_weights,"mat3d_weights.txt");
	//save3DMatrixdouble(Matrix3D_hiddenweights,"mat3d_hiddenweights.txt");
	//saveMatrixdouble(mat_biashidden,"matrixbiashidden.txt");
	//saveMatrixdouble(mat_biasoutput,"matrixbiasoutput.txt");
	//saveMatrixdouble(mat_target,"matrixoutputstarget.txt");
	cout << endl;
}

void FNN_1HiddenLayer_Iris::save_model(const std::string &filename) 
{
	std::ofstream file(filename);
	file << std::fixed << std::setprecision(6);
	vector<vector<double>> mat_w_dummy(4, vector<double>(5, 0));
	vector<vector<double>> mat_whidden_dummy(5, vector<double>(3, 0));
	vector<double> vec_biashidden(5,0.0);
	vector<double> vec_biasoutput(3,0.0);
	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 5; ++j) 
		{
			file << weights[i][j] << " ";
			mat_w_dummy[i][j] = weights[i][j];
		}
	}
	for (int j = 0; j < 5; ++j) 
	{
		file << biashidden[j] << " ";
		vec_biashidden[j]= biashidden[j];
	}
	for (int j = 0; j < 5; ++j) 
	{
		for (int k = 0; k < 3; ++k) 
		{
			file << hiddenweights[j][k] << " ";
			mat_whidden_dummy[j][k] = hiddenweights[j][k];
		}
	}
	for (int k = 0; k < 3; ++k) 
	{
		file << biasoutput[k] << " ";
		vec_biasoutput[k]= biasoutput[k];
	}
	file.close();
	cout << "\nThe model weights:" << endl;
	printMatrix(mat_w_dummy);
	cout << "\nThe model hiddenweights:" << endl;
	printMatrix(mat_whidden_dummy);
	cout << "\nThe model bias in the hidden layer:" << endl;
	printVector(vec_biashidden);
	cout << "\nThe model bias in the output layer:" << endl;
	printVector(vec_biasoutput);

}

void FNN_1HiddenLayer_Iris::load_model(const std::string& filename) 
{
	std::ifstream file(filename);
	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 5; ++j) 
		{
			file >> weights[i][j];
		}
	}
	for (int j = 0; j < 5; ++j) 
	{
		file >> biashidden[j];
	}
	for (int j = 0; j < 5; ++j) 
	{
		for (int k = 0; k < 3; ++k) 
		{
			file >> hiddenweights[j][k];
		}
	}
	for (int k = 0; k < 3; ++k) 
	{
		file >> biasoutput[k];
	}
	file.close();
}

// Constructor initializing the array with random numbers
FNN_1HiddenLayer_Iris_ConjugateGradient::FNN_1HiddenLayer_Iris_ConjugateGradient() 
{
	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_input = 4;
	for (int i = 0; i < n_input; ++i) 
	{
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			weights[i][j] = He_initialization(n_input);
		}
	}
	for (int j = 0; j < n_hiddenlayer; ++j) 
	{
		for (int k = 0; k < n_classification; ++k) 
		{
			hiddenweights[j][k] = He_initialization(n_input);
		}
	}
	for (int j = 0; j < n_hiddenlayer; ++j) 
	{
		biashidden[j] = He_initialization(n_input);
        }
	for (int k = 0; k < n_classification; ++k) 
	{
		biasoutput[k] = He_initialization(n_input);
        }
	learning_rate = 0.51;
}
	

vector<double> FNN_1HiddenLayer_Iris_ConjugateGradient::predict(const vector<double> &inputs) 
{
	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_input = 4;
	vector<double> hiddenlayer(n_hiddenlayer, 0.0);
	vector<double> outputs(n_classification, 0.0);
	for (int j = 0; j < n_hiddenlayer; ++j) 
	{
		for (int i = 0; i < n_input; ++i) 
		{
			hiddenlayer[j] += inputs[i] * weights[i][j];
		}
		hiddenlayer[j] += biashidden[j];
		hiddenlayer[j] = ReLUActivationfunction(hiddenlayer[j]);
		//hiddenlayer[j] =  TanhActivationfunction(hiddenlayer[j]);
	}
	for (int k = 0; k < n_classification; ++k) 
	{
		for (int j = 0; j < n_hiddenlayer; ++j) 
		{
			outputs[k] += hiddenlayer[j] * hiddenweights[j][k];
		}
		outputs[k] += biasoutput[k];
	}
	vector<double>predicted_output_softmax = SoftMax_vectorresult_activationfunction(outputs);
	return predicted_output_softmax;
}

void FNN_1HiddenLayer_Iris_ConjugateGradient::train(vector<vector<double>> &X, vector<int> &y, int epochs) 
{
// May 25th, 2026
	cout << BLUE << BOLD << "\nTraining Progress:\n" << RESET << endl;

	int bar_width = 50;
	int n_classification = 3;
	int n_hiddenlayer = 5;
	int n_input = 4;
	int n_dk = n_input*n_hiddenlayer + n_hiddenlayer*n_classification + n_hiddenlayer + n_classification;

	vector<double> gk;
	vector<double> gk_previous(n_classification*n_hiddenlayer + n_input*n_hiddenlayer + n_classification + n_hiddenlayer, 0.0);
	vector<double> dk;
	vector<double> dk_previous(n_classification*n_hiddenlayer + n_input*n_hiddenlayer + n_classification + n_hiddenlayer, 0.0);

	double beta = 0.01;
	double new_learning_rate = learning_rate;
	vector<vector<double>> mat_target;
	vector<vector<double>> mat_gk;
	vector<vector<double>> mat_error1;
	vector<vector<double>> mat_error2;
	vector<double> learningrate;
	vector<double> vec_beta;
	vector<double> norm_gk;
	vector<double> gk_strongwolfe;
	vector<double> vec_error;

	//vector<vector<vector<double>>> Matrix3D_weights;
	//vector<vector<vector<double>>> Matrix3D_hiddenweights;
	
	//vector<vector<double>> mat_w(n_input, vector<double>(n_hiddenlayer, 0.0));
	//vector<vector<double>> mat_whidden(n_hiddenlayer, vector<double>(n_classification, 0.0));
	//vector<vector<double>> mat_biashidden;
	//vector<vector<double>> mat_biasoutput;
	//vector<vector<double>> mat_target;
	//vector<double> vec_biashidden(n_hiddenlayer,0.0);
	//vector<double> vec_biasoutput(n_classification,0.0);
	int reset_iter=0;
	for (int epoch = 0; epoch < epochs; ++epoch) 
	{
		double epoch_loss = 0.0;
		for (size_t iter = 0; iter < X.size(); ++iter) 
		{
			vector<vector<double>> mat_jacobian(n_classification, vector<double>(n_classification, 0.0)); // Softmax derivative
			vector<double> delta_outputlayer(n_classification,0.0);
			vector<double> outputs = predict(X[iter]);
			vector<double> target(n_classification, 0.0);
			target[y[iter]] = 1.0;
			//mat_target.push_back(outputs);	
			//mat_target.push_back(target);		
	 		vector<double> errors(n_classification, 0.0);
			//vector<double> errors1(n_classification, 0.0);
			for (int k = 0; k < n_classification; ++k) 
			{
				//errors[k] = target[k] - outputs[k];
				errors[k] = -target[k] * log(outputs[k]); // Categorical Cross-Entropy (CCE) Loss because we are using softmax activation function at the output layer 
				epoch_loss += errors[k];
			} 
			double errors_currentweights = std::accumulate(errors.begin(), errors.end(), 0.0);
			vec_error.push_back(std::accumulate(errors.begin(), errors.end(), 0.0));
			vector<double> hiddenlayer(n_hiddenlayer, 0.0);
			for (int j = 0; j < n_hiddenlayer; ++j) 
			{
				for (int i = 0; i < n_input; ++i) 
				{
					hiddenlayer[j] += X[iter][i] * weights[i][j];
				}
				hiddenlayer[j] += biashidden[j];
				//hiddenlayer[j] = SigmoidActivationfunction(hiddenlayer[j]);
				hiddenlayer[j] = ReLUActivationfunction(hiddenlayer[j]);
			}
			for (int k_row = 0; k_row < n_classification; ++k_row) 
			{
				for (int k_col = 0; k_col < n_classification; ++k_col) 
				{
					if (k_row == k_col)
					{
						mat_jacobian[k_row][k_col] = outputs[k_row]*(1 - outputs[k_col]);
					}
					else if (k_row != k_col)
					{
						mat_jacobian[k_row][k_col] = outputs[k_row]*(- outputs[k_col]);
					}
				}
			}
			for (int k_row = 0; k_row < n_classification; ++k_row) 
			{
				for (int k_col = 0; k_col < n_classification; ++k_col) 
				{
					delta_outputlayer[k_row] += mat_jacobian[k_row][k_col]*errors[k_col];
				}
			}
			// this delta code makes the hiddenweights to grow big and make the epoch loss > 1
			/*for (int k = 0; k < n_classification; ++k) 
			{
				delta_outputlayer[k] = outputs[k] - target[k];
			}*/
			
			if(epoch  == 0 &&  iter==0)
			{
				gk.clear();
				// Creating vector g_{k} with push_back, the first one here is g_{0}
				// compute g_{k}/ create vector g_{k} it is the default steepest descent or gradient descent algorithm
				// for weights
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					double sum_delta = 0;			
					for (int k = 0; k < n_classification; ++k) 
					{
						sum_delta += hiddenweights[j][k] *  delta_outputlayer[k] ; 
					}
					for (int i = 0; i < n_input; ++i) 
					{
						gk.push_back(sum_delta * ReLUDerivative(hiddenlayer[j]) * X[iter][i]);
					}
				}
				// for bias in hidden layer
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					double sum_delta = 0;			
					for (int k = 0; k < n_classification; ++k) 
					{
						sum_delta += hiddenweights[j][k] *  delta_outputlayer[k] ; 
					}
					gk.push_back(sum_delta * ReLUDerivative(hiddenlayer[j]));
				}
				// for hiddenweights			
				for (int k = 0; k < n_classification; ++k) 
				{
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						gk.push_back(delta_outputlayer[k] * hiddenlayer[j]);
					}
				}
				// for bias output
				for (int k = 0; k < n_classification; ++k) 
				{
					gk.push_back(delta_outputlayer[k]);
					
				}
				dk.clear(); 
				// d_{0} = -g_{0}
				for (int j = 0; j < n_dk ; ++j) 
				{
					dk.push_back(gk[j]);
				}
				// Compute the first weights w_{1} = w_{0} + learning_rate * g_{0}
				for (int k = 0; k < n_classification; ++k) 
				{
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						hiddenweights[j][k] += new_learning_rate * delta_outputlayer[k] * hiddenlayer[j];
						//mat_whidden[j][k] = hiddenweights[j][k];
					}
				}
				for (int k = 0; k < n_classification; ++k) 
				{
					biasoutput[k] += new_learning_rate * delta_outputlayer[k];
					//vec_biasoutput[k] = biasoutput[k];
				}
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					double sum_delta = 0;			
					for (int k = 0; k < n_classification; ++k) 
					{
						sum_delta += hiddenweights[j][k] *  delta_outputlayer[k] ; 
					}
					for (int i = 0; i < n_input; ++i) 
					{
						weights[i][j] += new_learning_rate * sum_delta * ReLUDerivative(hiddenlayer[j]) * X[iter][i];
						//mat_w[i][j] = weights[i][j];
					}
					biashidden[j] += new_learning_rate *  sum_delta * ReLUDerivative(hiddenlayer[j]) ;
					//vec_biashidden[j] = biashidden[j];
				}
				
			}
			
			if(epoch  != 0 &&  iter!=0)
			{
				gk_previous.clear();
				dk_previous.clear();
				for (int j = 0; j < n_dk ; ++j) 
				{
					gk_previous.push_back(gk[j]);
					dk_previous.push_back(dk[j]);
				}

				gk.clear();
				dk.clear();
				// Creating vector g_{k} with push_back, the first one here is g_{1}
				// compute g_{k}/ create vector g_{k} it is the default steepest descent or gradient descent algorithm
				// for weights
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					double sum_delta = 0;			
					for (int k = 0; k < n_classification; ++k) 
					{
						sum_delta += hiddenweights[j][k] *  delta_outputlayer[k] ; 
					}
					for (int i = 0; i < n_input; ++i) 
					{
						gk.push_back(sum_delta * ReLUDerivative(hiddenlayer[j]) * X[iter][i]);
					}
				}
				// for bias in hidden layer
				for (int j = 0; j < n_hiddenlayer; ++j) 
				{
					double sum_delta = 0;			
					for (int k = 0; k < n_classification; ++k) 
					{
						sum_delta += hiddenweights[j][k] *  delta_outputlayer[k] ; 
					}
					gk.push_back(sum_delta * ReLUDerivative(hiddenlayer[j]));
				}
				// for hiddenweights			
				for (int k = 0; k < n_classification; ++k) 
				{
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						gk.push_back(delta_outputlayer[k] * hiddenlayer[j]);
					}
				}
				// for bias output
				for (int k = 0; k < n_classification; ++k) 
				{
					gk.push_back(delta_outputlayer[k]);
					
				}
				//reverse(gk.begin(),gk.end()); // reverse the vector
				mat_gk.push_back(gk);

				vector<double> gk_diff = subtract(gk,gk_previous);
						
				// Compute beta Polak-Ribiere, Fletcher-Reeves, OFR
				//double mu = 1.1;
				//double betavalue = divisiond(dot(gk,gk_diff), dot(gk_previous,gk_previous)); // PR
				//double betavalue = divisiond(dot(gk,gk), dot(gk_previous,gk_previous));	// FR	
				//double betavalue = divisiond(dot(gk,gk), mu*abs(dot(gk,dk_previous)) + norm(gk_previous)*norm(dk_previous));		// OFR
				//if (norm(gk_previous) ==0 || norm(dk_previous) ==0 || betavalue > 100 || betavalue <1e-4) // probably beta problem
				//{
				//	betavalue = 0.1;
				//}
				
				// Compute d_{k+1} = -g_{k+1} + beta * d_{k}
				if(beta >= divisiond(norm(gk),norm(gk_previous)) && divisiond(norm(gk),norm(gk_previous)) < 0.01)
				{
					beta = divisiond(norm(gk),norm(gk_previous)); // GLANZ beta
				}	
				else if(beta < norm(gk)/norm(gk_previous))
				{
					beta = 0.01; // GLANZ beta
				}
				for (int j = 0; j < n_dk; ++j) 
				{
					//dk.push_back(gk[j]); // this is the normal steepest descent / gradient descent
					dk.push_back(gk[j] + beta*dk_previous[j]); //We have problem on this, we cannot follow the algorithm and get the fast convergence, there is a leakage
				}
				vec_beta.push_back(beta);
				/*

				 Compute the learning rate with Strong Wolfe line search, iterative method

				*/
				double delta = 0.1, sigma = 0.9;
				double low = 0, high = 1e10;
				int sw_iter = 10;
				for (int i = 0; i<sw_iter ; ++i)
				{
					vector<double> gk_strongwolfe_newweights;
					// Compute new weights based on current learning rate / step size w_{k} + learning_rate_{k} * d_{k}
					for (int k = 0; k < n_classification; ++k) 
					{
						for (int j = 0; j < n_hiddenlayer; ++j) 
						{
							hiddenweights[j][k] += new_learning_rate * dk[ n_input*n_hiddenlayer + n_hiddenlayer + j + n_hiddenlayer*k]; 
						}
					}
					for (int k = 0; k < n_classification; ++k) 
					{
						biasoutput[k] += new_learning_rate * dk[  n_input*n_hiddenlayer + n_hiddenlayer + n_classification*n_hiddenlayer + k];
					}
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						for (int i = 0; i < n_input; ++i) 
						{
							weights[i][j] += new_learning_rate * dk[i + n_input*j];	
						}
					}

					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						biashidden[j] += new_learning_rate *  dk[n_input*n_hiddenlayer + j] ;
					}

					// Compute delta of each layer
					vector<vector<double>> mat_jacobian_sw(n_classification, vector<double>(n_classification, 0.0)); // Softmax derivative
					vector<double> delta_outputlayer_sw(n_classification,0.0);
				
					vector<double> hiddenlayer_sw(n_hiddenlayer, 0.0);
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						for (int i = 0; i < n_input; ++i) 
						{
							hiddenlayer_sw[j] += X[iter][i] * weights[i][j];
						}
						hiddenlayer_sw[j] += biashidden[j];
						hiddenlayer_sw[j] = ReLUActivationfunction(hiddenlayer[j]);
					}
					for (int k_row = 0; k_row < n_classification; ++k_row) 
					{
						for (int k_col = 0; k_col < n_classification; ++k_col) 
						{
							if (k_row == k_col)
							{
								mat_jacobian_sw[k_row][k_col] = outputs[k_row]*(1 - outputs[k_col]);
							}
							else if (k_row != k_col)
							{
								mat_jacobian_sw[k_row][k_col] = outputs[k_row]*(- outputs[k_col]);
							}
						}
					}
					for (int k_row = 0; k_row < n_classification; ++k_row) 
					{
						for (int k_col = 0; k_col < n_classification; ++k_col) 
						{
							delta_outputlayer_sw[k_row] += mat_jacobian_sw[k_row][k_col]*errors[k_col];
						}
					}

					gk_strongwolfe_newweights.clear();
				
					// Creating vector g_{k} with push_back, for the strong wolfe conditions iteration
					// for weights
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						double sum_delta = 0;			
						for (int k = 0; k < n_classification; ++k) 
						{
							sum_delta += hiddenweights[j][k] *  delta_outputlayer_sw[k] ; 
						}
						for (int i = 0; i < n_input; ++i) 
						{
							gk_strongwolfe_newweights.push_back(sum_delta * ReLUDerivative(hiddenlayer[j]) * X[iter][i]);
						}
					}
					// for bias in hidden layer
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						double sum_delta = 0;			
						for (int k = 0; k < n_classification; ++k) 
						{
							sum_delta += hiddenweights[j][k] *  delta_outputlayer_sw[k] ; 
						}
						gk_strongwolfe_newweights.push_back(sum_delta * ReLUDerivative(hiddenlayer_sw[j]));
					}
					// for hiddenweights			
					for (int k = 0; k < n_classification; ++k) 
					{
						for (int j = 0; j < n_hiddenlayer; ++j) 
						{
							gk_strongwolfe_newweights.push_back(delta_outputlayer_sw[k] * hiddenlayer_sw[j]);
						}
					}
					// for bias output
					for (int k = 0; k < n_classification; ++k) 
					{
						gk_strongwolfe_newweights.push_back(delta_outputlayer_sw[k]);
						
					}
					// Compute the new error based on the new weights
					vector<double> new_outputs = predict(X[iter]);
	 				vector<double> new_errors(n_classification, 0.0);
					for (int k = 0; k < n_classification; ++k) 
					{
						new_errors[k] = -target[k] * log(new_outputs[k]); // Categorical Cross-Entropy (CCE) Loss because we are using softmax activation function at the output layer 
					} 
					double errors_newweights = std::accumulate(new_errors.begin(), new_errors.end(), 0.0);

					// Neutralize weights and bias
					for (int k = 0; k < n_classification; ++k) 
					{
						for (int j = 0; j < n_hiddenlayer; ++j) 
						{
							hiddenweights[j][k] -= new_learning_rate * dk[ n_input*n_hiddenlayer + n_hiddenlayer + j + n_hiddenlayer*k]; 
						}
					}
					for (int k = 0; k < n_classification; ++k) 
					{
						biasoutput[k] -= new_learning_rate * dk[  n_input*n_hiddenlayer + n_hiddenlayer + n_classification*n_hiddenlayer + k];
					}
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						for (int i = 0; i < n_input; ++i) 
						{
							weights[i][j] -= new_learning_rate * dk[i + n_input*j];	
						}
					}

					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						biashidden[j] -= new_learning_rate *  dk[n_input*n_hiddenlayer + j] ;
					}

					double strongwolfeeq1_lhs = errors_newweights;
					double strongwolfeeq1_rhs = errors_currentweights + delta*new_learning_rate*dot(gk,dk_previous);
					
					double strongwolfeeq2_lhs = abs(dot(gk_strongwolfe_newweights,dk_previous));
					double strongwolfeeq2_rhs = sigma*abs(dot(gk,dk_previous));

					if (strongwolfeeq2_lhs <= strongwolfeeq2_rhs && strongwolfeeq1_lhs <= strongwolfeeq1_rhs)
					{
						new_learning_rate = new_learning_rate;
						i = sw_iter;
					}
						
					// Check Armijo Condition
					if (strongwolfeeq1_lhs > strongwolfeeq1_rhs ) 
					{
						high = new_learning_rate;
						new_learning_rate = 0.5*(low + high);
					} 
					// Check Curvature Condition
					else if (strongwolfeeq2_lhs < strongwolfeeq2_rhs ) 
					{
						/*low = new_learning_rate;
						if (high > 1e9) 
						{
							new_learning_rate *= 2.0;
						}
						else if (high <= 1e9) 
						{
							new_learning_rate = 0.5*(low + high) ;
						}*/
						new_learning_rate = new_learning_rate;
						i= sw_iter;
					} 

					if (new_learning_rate < 1e-6) 
					{
						new_learning_rate =  1e-1; // return small number if the learning rate becoming too small.
						i = sw_iter; // prevent infinite loop
					}
					/*if(epoch_loss / X.size() > 5)
					{
						high = new_learning_rate;
						new_learning_rate = (low + high);
						i = sw_iter;
					}*/
					
				}// End of strong wolfe conditions iteration for inexact line search

				learningrate.push_back(new_learning_rate);


				if ((reset_iter) % (n_dk) != 0) 
				{
					for (int k = 0; k < n_classification; ++k) 
					{
						for (int j = 0; j < n_hiddenlayer; ++j) 
						{
							hiddenweights[j][k] += new_learning_rate * dk[ n_input*n_hiddenlayer + n_hiddenlayer + j + n_hiddenlayer*k]; 
							//mat_whidden[j][k] = hiddenweights[j][k];
						}
					}
					for (int k = 0; k < n_classification; ++k) 
					{
						biasoutput[k] += new_learning_rate * dk[  n_input*n_hiddenlayer + n_hiddenlayer + n_classification*n_hiddenlayer + k];
						//vec_biasoutput[k] = biasoutput[k];
					}
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						for (int i = 0; i < n_input; ++i) 
						{
							weights[i][j] += new_learning_rate * dk[i + n_input*j];	
							//mat_w[i][j] = weights[i][j];
						}
					}

					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						biashidden[j] += new_learning_rate *  dk[n_input*n_hiddenlayer + j] ;
						//vec_biashidden[j] = biashidden[j];
					}
				}
				if ((reset_iter) % (n_dk) == 0) // Periodically reset d_{k} = -g_{k}
				{
					for (int k = 0; k < n_classification; ++k) 
					{
						for (int j = 0; j < n_hiddenlayer; ++j) 
						{
							hiddenweights[j][k] += new_learning_rate *  gk[ n_input*n_hiddenlayer + n_hiddenlayer + j + n_hiddenlayer*k]; 
							//mat_whidden[j][k] = hiddenweights[j][k];
						}
					}
					for (int k = 0; k < n_classification; ++k) 
					{
						biasoutput[k] += new_learning_rate * gk[  n_input*n_hiddenlayer + n_hiddenlayer + n_classification*n_hiddenlayer + k];
						//vec_biasoutput[k] = biasoutput[k];
					}
					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						for (int i = 0; i < n_input; ++i) 
						{
							weights[i][j] += new_learning_rate * gk[i + n_input*j];	
							//mat_w[i][j] = weights[i][j];
						}
					}

					for (int j = 0; j < n_hiddenlayer; ++j) 
					{
						biashidden[j] += new_learning_rate *  gk[n_input*n_hiddenlayer + j] ;
						//vec_biashidden[j] = biashidden[j];
					}
				}
			}
	
			reset_iter += 1;	
		}// end iterating all training samples

		if ((epoch + 1) % (epochs / 100) == 0 || epoch == epochs - 1) 
		{
			float progress = static_cast<float>(epoch + 1) / epochs;
			int pos = static_cast<int>(bar_width * progress);

			std::cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) 
				{
					std::cout << "=";
				}
				else if (i == pos)
				{
					std::cout << ">";
				}
				else
				{
					std::cout << " ";
				}
			}
			std::cout << "] " << int(progress * 100.0) << "% ";
			std::cout << "Epoch " << epoch + 1 << "/" << epochs << " - Loss: " << std::fixed << std::setprecision(4) << epoch_loss / X.size() << "\r";
			std::cout.flush();
		}
		//Matrix3D_weights.push_back(mat_w);
		//Matrix3D_hiddenweights.push_back(mat_whidden);
		//mat_biashidden.push_back(vec_biashidden);
		//mat_biasoutput.push_back(vec_biasoutput);
		
	} 
	//save3DMatrixdouble(Matrix3D_weights,"mat3d_weights.txt");
	//save3DMatrixdouble(Matrix3D_hiddenweights,"mat3d_hiddenweights.txt");
	//saveMatrixdouble(mat_biashidden,"matrixbiashidden.txt");
	//saveMatrixdouble(mat_biasoutput,"matrixbiasoutput.txt");
	//saveMatrixdouble(mat_target,"matrixoutputstarget.txt");
	//saveMatrixdouble(mat_gk,"matrixgk.txt");
	saveVectordouble(vec_beta,"beta.txt");
	saveVectordouble(learningrate,"strongwolfelearningrate.txt");
	saveVectordouble(vec_error,"error.txt");
	cout << endl;
}

void FNN_1HiddenLayer_Iris_ConjugateGradient::save_model(const std::string &filename) 
{
	std::ofstream file(filename);
	file << std::fixed << std::setprecision(6);
	vector<vector<double>> mat_w_dummy(4, vector<double>(5, 0));
	vector<vector<double>> mat_whidden_dummy(5, vector<double>(3, 0));
	vector<double> vec_biashidden(5,0.0);
	vector<double> vec_biasoutput(3,0.0);
	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 5; ++j) 
		{
			file << weights[i][j] << " ";
			mat_w_dummy[i][j] = weights[i][j];
		}
	}
	for (int j = 0; j < 5; ++j) 
	{
		file << biashidden[j] << " ";
		vec_biashidden[j]= biashidden[j];
	}
	for (int j = 0; j < 5; ++j) 
	{
		for (int k = 0; k < 3; ++k) 
		{
			file << hiddenweights[j][k] << " ";
			mat_whidden_dummy[j][k] = hiddenweights[j][k];
		}
	}
	for (int k = 0; k < 3; ++k) 
	{
		file << biasoutput[k] << " ";
		vec_biasoutput[k]= biasoutput[k];
	}
	file.close();
	cout << "\nThe model weights:" << endl;
	printMatrix(mat_w_dummy);
	cout << "\nThe model hiddenweights:" << endl;
	printMatrix(mat_whidden_dummy);
	cout << "\nThe model bias in the hidden layer:" << endl;
	printVector(vec_biashidden);
	cout << "\nThe model bias in the output layer:" << endl;
	printVector(vec_biasoutput);

}

void FNN_1HiddenLayer_Iris_ConjugateGradient::load_model(const std::string& filename) 
{
	std::ifstream file(filename);
	for (int i = 0; i < 4; ++i) 
	{
		for (int j = 0; j < 5; ++j) 
		{
			file >> weights[i][j];
		}
	}
	for (int j = 0; j < 5; ++j) 
	{
		file >> biashidden[j];
	}
	for (int j = 0; j < 5; ++j) 
	{
		for (int k = 0; k < 3; ++k) 
		{
			file >> hiddenweights[j][k];
		}
	}
	for (int k = 0; k < 3; ++k) 
	{
		file >> biasoutput[k];
	}
	file.close();
}

void print_confusion_matrix(const vector<vector<int>>& confusion_matrix, const vector<string>& label_names) 
{
	cout << CYAN << BOLD << "\nConfusion Matrix:\n" << RESET << std::endl;
    
	// Calculate max width for labels
	size_t max_width = std::max_element(label_names.begin(), label_names.end(),
		[](const std::string& a, const std::string& b) { return a.length() < b.length(); })->length();
	max_width = std::max(max_width, size_t(15));  // Minimum width of 15

	// Print header
	cout << std::setw(max_width + 2) << "Predicted >";
	for (const auto &name : label_names) 
	{
		cout << std::setw(max_width + 2) << name;
	}
	cout << endl;

	// Print rows
	for (size_t i = 0; i < confusion_matrix.size(); ++i) 
	{
		cout << std::setw(max_width + 2) << label_names[i];
		for (size_t j = 0; j < confusion_matrix[i].size(); ++j) 
		{
			cout << std::setw(max_width + 2) << confusion_matrix[i][j];
		}
		cout << endl;
	}
}


void print_metrics(const vector<vector<int>>& confusion_matrix, const vector<string>& label_names) 
{
	std::cout << MAGENTA << BOLD << "\nPrecision, Recall, and F1-score:\n" << RESET << std::endl;
    
	for (size_t i = 0; i < label_names.size(); ++i) 
	{
		int true_positive = confusion_matrix[i][i];
		int false_positive = 0;
		int false_negative = 0;

		for (size_t j = 0; j < confusion_matrix.size(); ++j) 
		{
			if (i != j) 
			{
				false_positive += confusion_matrix[j][i];
				false_negative += confusion_matrix[i][j];
			}
		}

		double precision = true_positive / static_cast<double>(true_positive + false_positive);
		double recall = true_positive / static_cast<double>(true_positive + false_negative);
		double f1_score = 2 * (precision * recall) / (precision + recall);

		cout << BOLD << label_names[i] << ":" << RESET << endl;
		cout << "  Precision: " << std::fixed << std::setprecision(2) << precision << endl;
		cout << "  Recall:    " << std::fixed << std::setprecision(2) << recall << endl;
		cout << "  F1-score:  " << std::fixed << std::setprecision(2) << f1_score << endl;
		cout << std::endl;
	}
}

void evaluate_model(FNN_NoHiddenLayer_Iris &nn, vector<vector<double>> &X, vector<int> &y, const vector<string> &label_names) 
{
	int correct = 0;
	vector<vector<int>> confusion_matrix(3, vector<int>(3, 0));

	cout << BLUE << BOLD << "\nEvaluating the model on the test set:" << RESET << endl;

	int progress = 0;
	int bar_width = 50;
	for (size_t i = 0; i < X.size(); ++i) 
	{
		vector<double> prediction = nn.predict(X[i]);
		int predicted_class = std::distance(prediction.begin(), std::max_element(prediction.begin(), prediction.end()));
        
		if (predicted_class == y[i]) 
		{
			correct++;
		}
        
		confusion_matrix[y[i]][predicted_class]++;

		// Update progress bar
		int new_progress = static_cast<int>((i + 1) * 100 / X.size());
		if (new_progress > progress) 
		{
			progress = new_progress;
			int pos = bar_width * progress / 100;
			cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) cout << "=";
				else if (i == pos) cout << ">";
				else cout << " ";
			}
			cout << "] " << progress << "%\r";
			cout.flush();
		}
	}
	cout << endl;

	double accuracy = static_cast<double>(correct) / X.size();
	cout << GREEN << BOLD << "\nAccuracy: " << std::fixed << std::setprecision(2) << accuracy * 100 << "%" << RESET << endl;

	print_confusion_matrix(confusion_matrix, label_names);
	print_metrics(confusion_matrix, label_names);
}

void evaluate_model(FNN_NoHiddenLayer_Iris_ConjugateGradient &nn, vector<vector<double>> &X, vector<int> &y, const vector<string> &label_names) 
{
	int correct = 0;
	vector<vector<int>> confusion_matrix(3, vector<int>(3, 0));

	cout << BLUE << BOLD << "\nEvaluating the model on the test set:" << RESET << endl;

	int progress = 0;
	int bar_width = 50;
	for (size_t i = 0; i < X.size(); ++i) 
	{
		vector<double> prediction = nn.predict(X[i]);
		int predicted_class = std::distance(prediction.begin(), std::max_element(prediction.begin(), prediction.end()));
        
		if (predicted_class == y[i]) 
		{
			correct++;
		}
        
		confusion_matrix[y[i]][predicted_class]++;

		// Update progress bar
		int new_progress = static_cast<int>((i + 1) * 100 / X.size());
		if (new_progress > progress) 
		{
			progress = new_progress;
			int pos = bar_width * progress / 100;
			cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) cout << "=";
				else if (i == pos) cout << ">";
				else cout << " ";
			}
			cout << "] " << progress << "%\r";
			cout.flush();
		}
	}
	cout << endl;

	double accuracy = static_cast<double>(correct) / X.size();
	cout << GREEN << BOLD << "\nAccuracy: " << std::fixed << std::setprecision(2) << accuracy * 100 << "%" << RESET << endl;

	print_confusion_matrix(confusion_matrix, label_names);
	print_metrics(confusion_matrix, label_names);
}

void evaluate_model(FNN_1HiddenLayer_Iris &nn, vector<vector<double>> &X, vector<int> &y, const vector<string> &label_names) 
{
	int correct = 0;
	vector<vector<int>> confusion_matrix(3, vector<int>(3, 0));

	cout << BLUE << BOLD << "\nEvaluating the model on the test set:" << RESET << endl;

	int progress = 0;
	int bar_width = 50;
	for (size_t i = 0; i < X.size(); ++i) 
	{
		vector<double> prediction = nn.predict(X[i]);
		int predicted_class = std::distance(prediction.begin(), std::max_element(prediction.begin(), prediction.end()));
        
		if (predicted_class == y[i]) 
		{
			correct++;
		}
        
		confusion_matrix[y[i]][predicted_class]++;

		// Update progress bar
		int new_progress = static_cast<int>((i + 1) * 100 / X.size());
		if (new_progress > progress) 
		{
			progress = new_progress;
			int pos = bar_width * progress / 100;
			cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) cout << "=";
				else if (i == pos) cout << ">";
				else cout << " ";
			}
			cout << "] " << progress << "%\r";
			cout.flush();
		}
	}
	cout << endl;

	double accuracy = static_cast<double>(correct) / X.size();
	cout << GREEN << BOLD << "\nAccuracy: " << std::fixed << std::setprecision(2) << accuracy * 100 << "%" << RESET << endl;

	print_confusion_matrix(confusion_matrix, label_names);
	print_metrics(confusion_matrix, label_names);
}

void evaluate_model(FNN_1HiddenLayer_Iris_ConjugateGradient &nn, vector<vector<double>> &X, vector<int> &y, const vector<string> &label_names) 
{
	int correct = 0;
	vector<vector<int>> confusion_matrix(3, vector<int>(3, 0));

	cout << BLUE << BOLD << "\nEvaluating the model on the test set:" << RESET << endl;

	int progress = 0;
	int bar_width = 50;
	for (size_t i = 0; i < X.size(); ++i) 
	{
		vector<double> prediction = nn.predict(X[i]);
		int predicted_class = std::distance(prediction.begin(), std::max_element(prediction.begin(), prediction.end()));
        
		if (predicted_class == y[i]) 
		{
			correct++;
		}
        
		confusion_matrix[y[i]][predicted_class]++;

		// Update progress bar
		int new_progress = static_cast<int>((i + 1) * 100 / X.size());
		if (new_progress > progress) 
		{
			progress = new_progress;
			int pos = bar_width * progress / 100;
			cout << "[";
			for (int i = 0; i < bar_width; ++i) 
			{
				if (i < pos) cout << "=";
				else if (i == pos) cout << ">";
				else cout << " ";
			}
			cout << "] " << progress << "%\r";
			cout.flush();
		}
	}
	cout << endl;

	double accuracy = static_cast<double>(correct) / X.size();
	cout << GREEN << BOLD << "\nAccuracy: " << std::fixed << std::setprecision(2) << accuracy * 100 << "%" << RESET << endl;

	print_confusion_matrix(confusion_matrix, label_names);
	print_metrics(confusion_matrix, label_names);
}

#endif
#endif