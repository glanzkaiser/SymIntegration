/*
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This library is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include <random>
#include <iostream>
#include <vector>

#ifndef SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_FORWARD

#ifndef FNN_NOHIDDENLAYER_IRIS_H
#define FNN_NOHIDDENLAYER_IRIS_H

class FNN_NoHiddenLayer_Iris {
private:
	double weights[4][3];
	double bias[3]; // Fixed-size array member
	double learning_rate;

public:
	FNN_NoHiddenLayer_Iris(); // Constructor
	
	vector<double> predict(const vector<double> &inputs);

	void train(vector<vector<double>> &X, vector<int> &y, int epochs);

	void save_model(const string &filename) ; // Void function declaration
	void load_model(const string& filename); 
};

#endif

#ifndef FNN_NOHIDDENLAYER_IRIS_CONJUGATEGRADIENT_H
#define FNN_NOHIDDENLAYER_IRIS_CONJUGATEGRADIENT_H

class FNN_NoHiddenLayer_Iris_ConjugateGradient {
private:
	double weights[4][3];
	double bias[3]; // Fixed-size array member
	double learning_rate;

public:
	FNN_NoHiddenLayer_Iris_ConjugateGradient(); // Constructor
	
	vector<double> predict(const vector<double> &inputs);
	
	void train(vector<vector<double>> &X, vector<int> &y, int epochs);

	void save_model(const string &filename) ; // Void function declaration
	void load_model(const string& filename); 
};

#endif

#ifndef FNN_1HIDDENLAYER_IRIS_H
#define FNN_1HIDDENLAYER_IRIS_H

class FNN_1HiddenLayer_Iris {
private:
	double weights[4][5];
	double hiddenweights[5][3];
	double biashidden[5]; // Fixed-size array member
	double biasoutput[3];
	double learning_rate;

public:
	FNN_1HiddenLayer_Iris(); // Constructor
	
	vector<double> predict(const vector<double> &inputs);

	void train(vector<vector<double>> &X, vector<int> &y, int epochs);

	void save_model(const string &filename) ; // Void function declaration
	void load_model(const string& filename); 
};

#endif

#ifndef FNN_1HIDDENLAYER_IRIS_CONJUGATEGRADIENT_H
#define FNN_1HIDDENLAYER_IRIS_CONJUGATEGRADIENT_H

class FNN_1HiddenLayer_Iris_ConjugateGradient {
private:
	double weights[4][5];
	double hiddenweights[5][3];
	double biashidden[5]; // Fixed-size array member
	double biasoutput[3];
	double learning_rate;

public:
	FNN_1HiddenLayer_Iris_ConjugateGradient(); // Constructor
	
	vector<double> predict(const vector<double> &inputs);
	
	void train(vector<vector<double>> &X, vector<int> &y, int epochs);

	void save_model(const string &filename) ; // Void function declaration
	void load_model(const string& filename); 
};

#endif

#ifndef CNN_LENET5_H
#define CNN_LENET5_H

class CNN_LeNet5 {
private:
	
	double c1_weights_kernel[6][5][5]; // first convolutional layer, 6 feature maps
	double s2_weights_kernel[6];
	double c3_weights_kernel_first6[6][3][5][5]; // second convolutional layer, 16 feature maps
	double c3_weights_kernel_next6[6][4][5][5]; // second convolutional layer, 16 feature maps
	double c3_weights_kernel_next3[3][4][5][5]; // second convolutional layer, 16 feature maps
	double c3_weights_kernel_last1[6][5][5]; // second convolutional layer, 16 feature maps

	double s4_weights_kernel[16];

	double c5_weights_kernel[120][16][5][5]; // third convolutional layer, 120 feature maps, kernel size of 5x5
	
	double c1_bias[6];
	double s2_bias[6];
	double c3_bias[16];
	double s4_bias[16];
	double c5_bias[120];

	vector<vector<double>> myMatrix;

	double biashidden[5]; // Fixed-size array member
	double biasoutput[3];
	double learning_rate;

public:
	CNN_LeNet5(); // Constructor	
	// 6 kernels for the first convolution layer 
	//vector<vector<vector<double>>> c1_weights(6, vector<vector<double>>(5,vector<double>(5,0.0)));
	// 16 kernels for the second convolution layer
	//vector<vector<vector<double>>> c2_weights(16, vector<vector<double>>(5,vector<double>(5,0.0)));

	// Constructor to initialize the 2D vector (optional)
	CNN_LeNet5(int rows, int cols) 
        : myMatrix(rows, vector<double>(cols, 0.0)) {}

	vector<double> predict(const vector<vector<int>> &inputs);

	void train(vector<vector<double>> &X, vector<int> &y, int epochs);

	void save_model(const string &filename) ; // Void function declaration
	void load_model(const string& filename); 
};

#endif

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS
#ifndef SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_ARTIFICIALNEURALNETWORKS_DECLARE

double LeakyReLUActivationfunction(double, double);
double ReLUActivationfunction(double);
double SigmoidActivationfunction(double);
double TanhActivationfunction(double);

double LeakyReluDerivative(double, double);
double ReLUDerivative(double);
double SigmoidDerivative(double); 
double TanhDerivative(double);

void ReLU_activationfunction(vector<double>&, vector<vector<double>>&, vector<double>& );
void Sigmoid_activationfunction(vector<double>&, vector<vector<double>>&, vector<double>& );
void Tanh_activationfunction(vector<double>&, vector<vector<double>>&, vector<double>& );
void Softmax_activationfunction(vector<double>&, vector<vector<double>>&, vector<double>& );
void SoftMax_logitsinput_activationfunction(vector<double>&);

vector<double> SoftMax_vectorresult_activationfunction(vector<double>& );

void print_centered(const std::string&, char , int);
double get_double_input(const string&);
void print_progress_bar(double);
void delay() ;

vector<vector<int>> CNN_2DpadBorder(const vector<vector<int>>&, int);
vector<double> CNN_1DConvolutionOperation(const vector<double>&, const vector<double>&);
void CNN_2DConvolutionOperation(vector<vector<int>>&, vector<vector<double>>&);
vector<vector<double>> CNN_2DConvolutionOperation(const vector<vector<int>>& , const vector<vector<double>>& , int);
vector<vector<double>> CNN_2DConvolutionOperation(const vector<vector<double>>& , const vector<vector<double>>& , int);
vector<vector<double>> CNN_2DmaxPooling(const vector<vector<double>>& , int, int);
vector<vector<int>> CNN_2DmaxPooling(const vector<vector<int>>& , int, int);
vector<vector<double>> CNN_2DaveragePooling(const vector<vector<int>>&, int, int);

void FNN_1hiddenlayer(vector<double>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>& );
void FNN_1hiddenlayer(vector<double>&, vector<double>&, int );

void FNN_NoHiddenLayer_Iris_irisdataset_training(vector<vector<double>>&, vector<string>&);
void FNN_NoHiddenLayer_Iris_irisdataset_testing(vector<vector<double>>&, vector<string>&);

void FNN_1hiddenlayer_irisdataset_training(vector<vector<double>>&, vector<string>&);
void FNN_1hiddenlayer_irisdataset_testing(vector<vector<double>>&, vector<string>&);

vector<vector<double>> read_dataset_iris(const string &, vector<int> &, vector<string> &); 

void print_confusion_matrix(const vector<vector<int>>& , const vector<string>& );
void print_metrics(const vector<vector<int>>& , const vector<string>& );
void evaluate_model(FNN_NoHiddenLayer_Iris &, vector<vector<double>> &, vector<int> &, const vector<string> &);
void evaluate_model(FNN_NoHiddenLayer_Iris_ConjugateGradient &, vector<vector<double>> &, vector<int> &, const vector<string> &);
void evaluate_model(FNN_1HiddenLayer_Iris &, vector<vector<double>> &, vector<int> &, const vector<string> &);
void evaluate_model(FNN_1HiddenLayer_Iris_ConjugateGradient &, vector<vector<double>> &, vector<int> &, const vector<string> &);
#endif
#endif


#endif
