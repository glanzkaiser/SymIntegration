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


class AdamOptimizer {
private:
	double learning_rate;
	double beta1;
	double beta2;
	double epsilon;
	int timestep;

	vector<double> m; // First moments
	vector<double> v; // Second moments

public:
	AdamOptimizer(size_t num_params, double learning_rate = 0.001f, double beta1 = 0.9f, double beta2 = 0.999f, double epsilon = 1e-8f)
	: learning_rate(learning_rate), beta1(beta1), beta2(beta2), epsilon(epsilon), timestep(0) 
	{
		m.resize(num_params, 0.0f);
		v.resize(num_params, 0.0f);
	}

	void update(std::vector<double>& weights, const std::vector<double>& gradients) 
	{
		if (weights.size() != gradients.size() || weights.size() != m.size()) 
		{
			throw std::invalid_argument("Vector size mismatch.");
		}

		// Increment the overall training step index
		timestep++; 

		// Prefetch bias corrections to avoid calculating them inside the loop
		double bias_correction1 = 1.0f - std::pow(beta1, timestep);
		float bias_correction2 = 1.0f - std::pow(beta2, timestep);

		for (size_t i = 0; i < weights.size(); ++i) 
		{
		// Update moment estimates
		m[i] = beta1 * m[i] + (1.0f - beta1) * gradients[i];
		v[i] = beta2 * v[i] + (1.0f - beta2) * (gradients[i] * gradients[i]);

		// Compute corrected estimates
		double m_hat = m[i] / bias_correction1;
		double v_hat = v[i] / bias_correction2;

		// Update step
		weights[i] -= (learning_rate * m_hat) / (std::sqrt(v_hat) + epsilon);
		}
	}
};


// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	// Track 3 model weights
	vector<double> weights = {0.5f, -0.2f, 0.1f};
    
	// Instantiate optimizer for 3 parameters
	AdamOptimizer optimizer(weights.size(), 0.01f); 

	// Simulate 3 training steps
	for (int step = 1; step <= 3; ++step) 
	{
		// Dummy mock gradients (usually computed via backpropagation)
		vector<double> gradients = {0.1f, -0.05f, 0.02f}; 

		optimizer.update(weights, gradients);

		cout << "Weights after step " << step << ": ";
		for (double w : weights) std::cout << w << " ";
		cout << "\n";
	}
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}