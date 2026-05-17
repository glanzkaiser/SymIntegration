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

void FNN_1hiddenlayer(vector<double>&, vector<vector<double>>&, vector<vector<double>>&, vector<double>&, vector<double>&, vector<double>& );
void FNN_1hiddenlayer(vector<double>&, vector<double>&, int );
void FNN_nohiddenlayer_irisdataset_training(vector<vector<double>>&, vector<string>&);
void FNN_nohiddenlayer_irisdataset_training_conjugategradient(vector<vector<double>>&, vector<string>&);
void FNN_nohiddenlayer_irisdataset_testing(vector<vector<double>>&, vector<string>&);

void FNN_1hiddenlayer_irisdataset_training(vector<vector<double>>&, vector<string>&);
void FNN_1hiddenlayer_irisdataset_training_conjugategradient(vector<vector<double>>&, vector<string>&);
void FNN_1hiddenlayer_irisdataset_testing(vector<vector<double>>&, vector<string>&);


#endif
#endif


#endif
