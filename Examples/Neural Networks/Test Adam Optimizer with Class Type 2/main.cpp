// g++ -o result main.cpp -lsymintegration
// Merci beaucoup Freya et Sentinel

#include<bits/stdc++.h>
#include<iostream>
#include "symintegrationc++.h"
#include<vector>
#include <chrono>
#include <algorithm> // For std::next_permutation
#include <string>
using namespace std::chrono;
using namespace std;


class AdamOptimizer1 {
private:
public:

	AdamOptimizer1(){}; // Constructor	
	
vector<double> output(vector<double>& m, vector<double>& v, double beta1, double beta2, const vector<double>& gradients, int timestep) 
{
	if (m.size() != gradients.size() ||  gradients.size() != v.size() || m.size() != v.size()) 
	{
		throw std::invalid_argument("Vector size mismatch.");
	}
	int n = gradients.size();
	vector<double> adamresult(n,0.0);
	double epsilon = 1e-8f;

	// Prefetch bias corrections to avoid calculating them inside the loop
	double bias_correction1 = 1.0f - std::pow(beta1, timestep);
	float bias_correction2 = 1.0f - std::pow(beta2, timestep);
	for (int i = 0; i < n; ++i) 
	{
		// Update moment estimates
		m[i] = beta1 * m[i] + (1.0f - beta1) * gradients[i];
		v[i] = beta2 * v[i] + (1.0f - beta2) * (gradients[i] * gradients[i]);

		// Compute corrected estimates
		double m_hat = m[i] / bias_correction1;
		double v_hat = v[i] / bias_correction2;

		// Update step
		adamresult[i] = ( m_hat) / (std::sqrt(v_hat) + epsilon);
	}
	return adamresult;
}

void update(vector<double>& m, vector<double>& v, double beta1, double beta2, const vector<double>& gradients) 
{
	int n = gradients.size();

	for (int i = 0; i < n; ++i) 
	{
		// Update moment estimates
		m[i] = beta1 * m[i] + (1.0f - beta1) * gradients[i];
		v[i] = beta2 * v[i] + (1.0f - beta2) * (gradients[i] * gradients[i]);
	}
}
};


// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Track 3 model weights
    	vector<double> adam_m_bias_output(10,0.0);
	vector<double> adam_v_bias_output(10,0.0);
	vector<double> delta_outputlayer(10,0.12); 
	// Instantiate optimizer for 3 parameters

	AdamOptimizer1 optimizer_output_bias;

	double beta1 = 0.9, beta2 = 0.999;
	// Simulate 3 training steps
	for (int step = 1; step <= 3; ++step) 
	{
		// Dummy mock gradients (usually computed via backpropagation)
		
		vector<double> adamoptimizer_bias_output = optimizer_output_bias.output(adam_m_bias_output, adam_v_bias_output, beta1, beta2, delta_outputlayer, step ); 
		optimizer_output_bias.update(adam_m_bias_output, adam_v_bias_output, beta1, beta2, delta_outputlayer);
		printVector(adamoptimizer_bias_output); 
		beta1 *= 0.99;
		beta2 *= 0.99;
	}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}