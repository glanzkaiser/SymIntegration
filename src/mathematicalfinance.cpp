/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_DEFINE

double presentvalue(double A, double r, int n)
{
	return A*((1-(1)/(pow(1+r,n)))/r);
}

double bondpricing(double C, double M, double r, int n)
{
	return (C*((1-(1)/(pow(1+r,n)))/r) ) + M/(pow(1+r,n));
}
double amortization(double L, double r, int n)
{
	double M = L*((r*(pow(1+r,n)))/(pow(1+r,n) - 1));
	double interest, principal, endingbalance, v;
	//endingbalance  = L;
	cout << "\nPeriod" << setw(23) << "Payment" << setw(23)<< "Interest" << setw(23) << "Principal" << setw(30) << "Outstanding Loan Balance" << endl;
	for (int i = 1; i <= n ; i++ )
	{
		v = pow(1+r,-1);
		interest = (1 - (pow(v,n-i+1)))*M; // (1-v^n) * M
		principal = M - interest;
		//endingbalance  -= principal ;
		endingbalance = ( ((1-pow(v,n-i+1))/(r)) - pow(v,n-i+1) )*M; // an - v^n
		if(i==n)
		{
			endingbalance  = 0;		
		}
		cout << i << setw(28) << M << setw(23) << interest << setw(23)  << principal << setw(23) << endingbalance << endl; 
	}
	cout << "\nTotal payment:" << endl;
	return M*n;
}

// 1. ARMA Residual Filter (Helper function)
vector<double> computeARMAResiduals(const std::vector<double>& data, double constant, const std::vector<double>& phi, const std::vector<double>& theta) 
{
	size_t N = data.size();
	size_t p = phi.size();
	size_t q = theta.size();
    
	vector<double> fitted(N, 0.0);
  	vector<double> epsilon(N, 0.0);

	// Calculate long-run mean for cold-start initialization
	double sum_phi = std::accumulate(phi.begin(), phi.end(), 0.0);
	double unconditional_mean = 0.0;
	if (std::abs(1.0 - sum_phi) > 1e-5) 
	{
		unconditional_mean = constant / (1.0 - sum_phi);
	} 
	else 
	{
		double sum_data = std::accumulate(data.begin(), data.end(), 0.0);
		unconditional_mean = sum_data / static_cast<double>(N);
	}

	for (size_t t = 0; t < N; ++t) 
	{
		double current_fit = constant;

		// AR component
		for (size_t i = 1; i <= p; ++i) 
		{
			current_fit += phi[i - 1] * ((t >= i) ? data[t - i] : unconditional_mean);
		}

		// MA component
		for (size_t j = 1; j <= q; ++j) 
		{
			if (t >= j) 
			{
				current_fit += theta[j - 1] * epsilon[t - j];
			}
		}

		fitted[t] = current_fit;
		epsilon[t] = data[t] - current_fit;
	}
	return epsilon;
}

// 2. Negative Log-Likelihood Function for ARMA
double computeARMANegativeLogLikelihood(const std::vector<double>& params, const std::vector<double>& data, size_t p, size_t q) 
{
	double constant = params[0];
	vector<double> phi(p);
	vector<double> theta(q);
    
	for (size_t i = 0; i < p; ++i) phi[i] = params[1 + i];
	for (size_t j = 0; j < q; ++j) theta[j] = params[1 + p + j];

	// Penalty for stationarity boundaries (sum of AR absolute roots should be < 1 roughly)
	double sum_ar = 0.0;
	for (double phi_val : phi) sum_ar += std::abs(phi_val);
	if (sum_ar >= 0.995) return 1e10;

	// Penalty for invertibility boundaries (sum of MA absolute roots should be < 1 roughly)
	double sum_ma = 0.0;
	for (double theta_val : theta) sum_ma += std::abs(theta_val);
	if (sum_ma >= 0.995) return 1e10;

	// Filter out residual arrays
	vector<double> epsilon = computeARMAResiduals(data, constant, phi, theta);
    
	// Compute residual variance (sigma^2)
	double sum_sq_errors = 0.0;
	for (double eps : epsilon) 
	{
		sum_sq_errors += eps * eps;
	}
	double N = static_cast<double>(data.size());
	double sigma2 = sum_sq_errors / N;

	if (sigma2 <= 1e-8) return 1e10; // Avoid log(0) errors

	// Calculate conditional log-likelihood
	double nll = (N / 2.0) * (std::log(2.0 * M_PI) + std::log(sigma2) + 1.0);
    
	return std::isnan(nll) ? 1e10 : nll;
}

// 3. Nelder-Mead Optimization implementation for ARMA
vector<double> optimizeARMAMLE(const std::vector<double>& data, size_t p, size_t q, int max_iter) 
{
	size_t dim = 1 + p + q; // constant + AR parameters + MA parameters
	size_t num_vertices = dim + 1;
    
	vector<std::vector<double>> simplex(num_vertices, std::vector<double>(dim));
	vector<double> f_values(num_vertices);

	// Initial Guess Construction
	double sum_data = std::accumulate(data.begin(), data.end(), 0.0);
	double initial_mean = sum_data / static_cast<double>(data.size());
    
	simplex[0][0] = initial_mean * 0.1; // Baseline intercept guess
	for (size_t i = 1; i <= p; ++i) simplex[0][i] = 0.1 / p; // AR guesses
	for (size_t j = 1; j <= q; ++j) simplex[0][1 + p + j - 1] = 0.1 / q; // MA guesses

	// Expand initial guess across vertices
	for (size_t i = 1; i < num_vertices; ++i) 
	{
		simplex[i] = simplex[0];
		simplex[i][i - 1] += (simplex[0][i - 1] == 0.0) ? 0.01 : 0.05 * simplex[0][i - 1];
	}

	// Evaluate current simplex points
	for (size_t i = 0; i < num_vertices; ++i) 
	{
		f_values[i] = computeARMANegativeLogLikelihood(simplex[i], data, p, q);
	}

	// Simplex shrinking loop
	for (int iter = 0; iter < max_iter; ++iter) 
	{
	 	vector<size_t> idx(num_vertices);
		std::iota(idx.begin(), idx.end(), 0);
		std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return f_values[a] < f_values[b]; });

		size_t best = idx[0];
		size_t worst = idx[dim];
		size_t next_worst = idx[dim - 1];

		if (std::abs(f_values[worst] - f_values[best]) < 1e-6) break;

		// Centroid computation
		vector<double> centroid(dim, 0.0);
		for (size_t i = 0; i < dim; ++i) 
		{
			for (size_t d = 0; d < dim; ++d) centroid[d] += simplex[idx[i]][d];
		}
		for (double& val : centroid) val /= static_cast<double>(dim);

		// Reflection
		vector<double> reflected(dim);
		for (size_t d = 0; d < dim; ++d) reflected[d] = centroid[d] + 1.0 * (centroid[d] - simplex[worst][d]);
		double f_reflected = computeARMANegativeLogLikelihood(reflected, data, p, q);

		if (f_reflected >= f_values[best] && f_reflected < f_values[next_worst]) 
		{
			simplex[worst] = reflected;
			f_values[worst] = f_reflected;
			continue;
		}

		// Expansion
		if (f_reflected < f_values[best]) 
		{
			vector<double> expanded(dim);
			for (size_t d = 0; d < dim; ++d) expanded[d] = centroid[d] + 2.0 * (reflected[d] - centroid[d]);
			double f_expanded = computeARMANegativeLogLikelihood(expanded, data, p, q);
            
			if (f_expanded < f_reflected) 
			{
				simplex[worst] = expanded;
				f_values[worst] = f_expanded;
			} 
			else 
			{
				simplex[worst] = reflected;
				f_values[worst] = f_reflected;
			}
			continue;
		}

		// Contraction
		vector<double> contracted(dim);
		for (size_t d = 0; d < dim; ++d) contracted[d] = centroid[d] + 0.5 * (simplex[worst][d] - centroid[d]);
		double f_contracted = computeARMANegativeLogLikelihood(contracted, data, p, q);

		if (f_contracted < f_values[worst]) 
		{
			simplex[worst] = contracted;
			f_values[worst] = f_contracted;
			continue;
		}

		// Shrink step fallback
		for (size_t i = 1; i < num_vertices; ++i) 
		{
			for (size_t d = 0; d < dim; ++d) 
			{
				simplex[idx[i]][d] = simplex[best][d] + 0.5 * (simplex[idx[i]][d] - simplex[best][d]);
			}
			f_values[idx[i]] = computeARMANegativeLogLikelihood(simplex[idx[i]], data, p, q);
		}
	}

	size_t absolute_best_idx = std::distance(f_values.begin(), std::min_element(f_values.begin(), f_values.end()));
	return simplex[absolute_best_idx];
}

// 4. Primary Public Interface Function for fitting ARMA
ARMAParameters fitARMA(const vector<double>& raw_data, size_t p, size_t q) 
{
	if (raw_data.size() <= (p + q + 1)) 
	{
		throw std::invalid_argument("Dataset sample size is insufficient for the targeted parameters.");
	}

	vector<double> optimal_params = optimizeARMAMLE(raw_data, p, q);

	ARMAParameters results;
	results.constant = optimal_params[0];
	results.phi.assign(optimal_params.begin() + 1, optimal_params.begin() + 1 + p);
	results.theta.assign(optimal_params.begin() + 1 + p, optimal_params.end());
    
	// Extract optimized statistical values
	vector<double> epsilon = computeARMAResiduals(raw_data, results.constant, results.phi, results.theta);
	double sum_sq_errors = 0.0;
	for (double eps : epsilon) sum_sq_errors += eps * eps;
    
	results.sigma2 = sum_sq_errors / static_cast<double>(raw_data.size());
	double min_nll = computeARMANegativeLogLikelihood(optimal_params, raw_data, p, q);
	results.log_likelihood = -min_nll;

	return results;
}

/**
 * Generates out-of-sample forecasts for a stationary ARMA(p,q) series.
 * 
 * Formula: y_hat_{T+h} = c + sum_{i=1}^p (phi_i * y_{T+h-i}) + sum_{j=1}^q (theta_j * epsilon_{T+h-j})
 * 
 * @param historical_data The observed data vector up to time T (y_t)
 * @param constant The constant intercept term (c)
 * @param phi Autoregressive (AR) coefficients vector, size p
 * @param theta Moving Average (MA) coefficients vector, size q
 * @param steps_to_forecast Number of forward periods to predict (h)
 * @return std::vector<double> Vector containing the forward predicted mean track points
 */
vector<double> forecastARMA(const std::vector<double>& historical_data, double constant, const vector<double>& phi, const vector<double>& theta,
                                 size_t steps_to_forecast) 
{
	if (steps_to_forecast == 0) return std::vector<double>();

	size_t N = historical_data.size();
	size_t p = phi.size();
	size_t q = theta.size();

	// 1. Reconstruct historical residual errors conditional on the dataset
	// We instantiate vectors to compute past fits and track historical innovations
	vector<double> fitted(N, 0.0);
	vector<double> epsilon(N, 0.0);

	double sum_phi = std::accumulate(phi.begin(), phi.end(), 0.0);
	double unconditional_mean = (std::abs(1.0 - sum_phi) > 1e-6) ? (constant / (1.0 - sum_phi)) : 0.0;

	for (size_t t = 0; t < N; ++t) 
	{
		double current_fit = constant;
		for (size_t i = 1; i <= p; ++i) 
		{
			current_fit += phi[i - 1] * ((t >= i) ? historical_data[t - i] : unconditional_mean);
		}
		for (size_t j = 1; j <= q; ++j) 
		{
			if (t >= j) current_fit += theta[j - 1] * epsilon[t - j];
		}
		fitted[t] = current_fit;
		epsilon[t] = historical_data[t] - current_fit;
	}

	// 2. Setup out-of-sample tracking structures
	vector<double> y_extended = historical_data;
	vector<double> eps_extended = epsilon;

	// Pad arrays to hold future values
	y_extended.resize(N + steps_to_forecast, 0.0);
	eps_extended.resize(N + steps_to_forecast, 0.0); // Future unobserved shocks default to 0.0

	// 3. Recursive out-of-sample forecasting loop
	for (size_t h = 0; h < steps_to_forecast; ++h) 
	{
		size_t t = N + h;
		double future_y = constant;

		// AR components: dynamic feedback from history or prior forecasted steps
		for (size_t i = 1; i <= p; ++i) 
		{
			future_y += phi[i - 1] * y_extended[t - i];
		}

		// MA components: rely only on known historical errors (decays to zero when h > q)
		for (size_t j = 1; j <= q; ++j) 
		{
			future_y += theta[j - 1] * eps_extended[t - j];
		}

		y_extended[t] = future_y;
		eps_extended[t] = 0.0; // The expected value of future white noise is 0
	}

	// Isolate and return only the out-of-sample forecast segment
	vector<double> forecast(y_extended.begin() + N, y_extended.end());
	return forecast;
}

// 1. GARCH Variance Engine (Helper function)
std::vector<double> computeGARCHVariances(const std::vector<double>& epsilon, double omega, 
	const std::vector<double>& alpha, const std::vector<double>& beta) 
{
	size_t N = epsilon.size();
	size_t q = alpha.size();
	size_t p = beta.size();
	vector<double> sigma2(N, 0.0);

	// Compute sample variance as initial fallback
	double sample_var = 0.0;
	for (double eps : epsilon) sample_var += eps * eps;
	sample_var /= static_cast<double>(N);

	double sum_alpha = std::accumulate(alpha.begin(), alpha.end(), 0.0);
	double sum_beta = std::accumulate(beta.begin(), beta.end(), 0.0);
    
	double unconditional_var = (1.0 - sum_alpha - sum_beta > 0.0) ? omega / (1.0 - sum_alpha - sum_beta) : sample_var;

	for (size_t t = 0; t < N; ++t) 
	{
		double current_variance = omega;
		for (size_t i = 1; i <= q; ++i) 
		{
			current_variance += alpha[i - 1] * ((t >= i) ? epsilon[t - i] * epsilon[t - i] : unconditional_var);
		}
		for (size_t j = 1; j <= p; ++j) 
		{
			current_variance += beta[j - 1] * ((t >= j) ? sigma2[t - j] : unconditional_var);
		}
		sigma2[t] = std::max(current_variance, 1e-6); // Prevent division by zero or negative variance
	}
	return sigma2;
}

// 2. Negative Log-Likelihood Function (The objective function to minimize)
double computeNegativeLogLikelihood(const std::vector<double>& params, const std::vector<double>& epsilon, 
                                    size_t p, size_t q) 
{
	double omega = params[0];
	vector<double> alpha(q);
	vector<double> beta(p);
    
 	for (size_t i = 0; i < q; ++i) alpha[i] = params[1 + i];
	for (size_t j = 0; j < p; ++j) beta[j] = params[1 + q + j];

	// Impose strict penalties for violating GARCH constraints
	if (omega <= 1e-4) return 1e10;
	double sum_coef = 0.0;
	for (double a : alpha) { if (a < 0) return 1e10; sum_coef += a; }
	for (double b : beta)  { if (b < 0) return 1e10; sum_coef += b; }
	if (sum_coef >= 0.999) return 1e10; // Ensure covariance stationarity

	vector<double> sigma2 = computeGARCHVariances(epsilon, omega, alpha, beta);
    
	double nll = 0.0;
	const double log_2pi = std::log(2.0 * M_PI);
    
	for (size_t t = 0; t < epsilon.size(); ++t) 
	{
		nll += 0.5 * (log_2pi + std::log(sigma2[t]) + (epsilon[t] * epsilon[t]) / sigma2[t]);
	}
	return std::isnan(nll) ? 1e10 : nll;
}

// 3. Nelder-Mead Optimization Algorithm implementation
vector<double> optimizeNelderMead(const std::vector<double>& epsilon, size_t p, size_t q, int max_iter) 
{
	size_t dim = 1 + q + p; // omega + alpha elements + beta elements
	size_t num_vertices = dim + 1;
    
	// Initialize Simplex matrix
	vector<std::vector<double>> simplex(num_vertices, std::vector<double>(dim));
	vector<double> f_values(num_vertices);

	// Initial guess vertex
	simplex[0][0] = 0.1; // Baseline omega guess
	for (size_t i = 1; i <= q; ++i) simplex[0][i] = 0.05; // alpha guess
	for (size_t j = 1; j <= p; ++j) simplex[0][1 + q + j - 1] = 0.7 / p; // beta guess

	// Perturb initial guess to build the rest of the simplex matrix
	for (size_t i = 1; i < num_vertices; ++i) 
	{
		simplex[i] = simplex[0];
		simplex[i][i - 1] += (simplex[0][i - 1] == 0.0) ? 0.00025 : 0.05 * simplex[0][i - 1];
	}

	// Evaluate the starting vertices
	for (size_t i = 0; i < num_vertices; ++i) 
	{
		f_values[i] = computeNegativeLogLikelihood(simplex[i], epsilon, p, q);
	}

	// Main optimization loop
	for (int iter = 0; iter < max_iter; ++iter) 
	{
		// Sort indices based on function evaluation values (ascending)
		vector<size_t> idx(num_vertices);
		std::iota(idx.begin(), idx.end(), 0);
		std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return f_values[a] < f_values[b]; });

		size_t best = idx[0];
		size_t worst = idx[dim];
		size_t next_worst = idx[dim - 1];

		// Break if variance between best and worst is microscopic
		if (std::abs(f_values[worst] - f_values[best]) < 1e-7) break;

		// Compute Centroid of all points except the worst one
		std::vector<double> centroid(dim, 0.0);
		for (size_t i = 0; i < dim; ++i) 
		{
			for (size_t d = 0; d < dim; ++d) 
			{
				centroid[d] += simplex[idx[i]][d];
			}
		}
		for (double& val : centroid) val /= static_cast<double>(dim);

		// Reflection step
		vector<double> reflected(dim);
		for (size_t d = 0; d < dim; ++d) reflected[d] = centroid[d] + 1.0 * (centroid[d] - simplex[worst][d]);
		double f_reflected = computeNegativeLogLikelihood(reflected, epsilon, p, q);

		if (f_reflected >= f_values[best] && f_reflected < f_values[next_worst]) 
		{
			simplex[worst] = reflected;
			f_values[worst] = f_reflected;
			continue;
		}	

		// Expansion step
		if (f_reflected < f_values[best]) 
		{
			vector<double> expanded(dim);
			for (size_t d = 0; d < dim; ++d) expanded[d] = centroid[d] + 2.0 * (reflected[d] - centroid[d]);
			double f_expanded = computeNegativeLogLikelihood(expanded, epsilon, p, q);
		    
			if (f_expanded < f_reflected) 
			{
				simplex[worst] = expanded;
				f_values[worst] = f_expanded;
			} 
			else 
			{
				simplex[worst] = reflected;
				f_values[worst] = f_reflected;
			}
			continue;
		}

		// Contraction step
		vector<double> contracted(dim);
		for (size_t d = 0; d < dim; ++d) contracted[d] = centroid[d] + 0.5 * (simplex[worst][d] - centroid[d]);
		double f_contracted = computeNegativeLogLikelihood(contracted, epsilon, p, q);

		if (f_contracted < f_values[worst]) 
		{
			simplex[worst] = contracted;
			f_values[worst] = f_contracted;
			continue;
		}

		// Shrink step (fallback if all other vectors fail)
		for (size_t i = 1; i < num_vertices; ++i) 
		{
			for (size_t d = 0; d < dim; ++d) 
			{
				simplex[idx[i]][d] = simplex[best][d] + 0.5 * (simplex[idx[i]][d] - simplex[best][d]);
			}
			f_values[idx[i]] = computeNegativeLogLikelihood(simplex[idx[i]], epsilon, p, q);
		}
	}

	// Return the absolute best parameters located
	size_t absolute_best_idx = std::distance(f_values.begin(), std::min_element(f_values.begin(), f_values.end()));
	return simplex[absolute_best_idx];
}

// 4. Primary Public Interface Function
GARCHParameters fitGARCH(const std::vector<double>& raw_data, size_t p, size_t q) 
{
	if (raw_data.size() <= (p + q)) 
	{
		throw std::invalid_argument("Dataset sample size is too small for the specified GARCH dimensions.");
	}

	// Demean raw return data to yield clean residuals
	double sum = std::accumulate(raw_data.begin(), raw_data.end(), 0.0);
	double mu = sum / static_cast<double>(raw_data.size());
	vector<double> epsilon(raw_data.size());
	for (size_t i = 0; i < raw_data.size(); ++i) 
	{
		epsilon[i] = raw_data[i] - mu;
	}

	// Run Optimization
	vector<double> optimal_params = optimizeNelderMead(epsilon, p, q);

	// Build return structure
	GARCHParameters results;
	results.omega = optimal_params[0];
	results.alpha.assign(optimal_params.begin() + 1, optimal_params.begin() + 1 + q);
	results.beta.assign(optimal_params.begin() + 1 + q, optimal_params.end());
    
	double min_nll = computeNegativeLogLikelihood(optimal_params, epsilon, p, q);
	results.log_likelihood = -min_nll; // Recover actual log-likelihood

	return results;
}

/**
 * Generates out-of-sample conditional variance forecasts for a GARCH(p,q) process.
 * 
 * @param raw_data The historical input series (e.g., returns or asset data used during fitting)
 * @param fit The structural results containing optimized MLE parameters from fitGARCH
 * @param steps_to_forecast Number of forward out-of-sample steps to predict (h)
 * @return GARCHForecast Parallel vectors of forecasted variances and volatilities
 */
GARCHForecast forecastGARCH(const std::vector<double>& raw_data, const GARCHParameters& fit, size_t steps_to_forecast) 
{
	if (steps_to_forecast == 0) return GARCHForecast();

	size_t N = raw_data.size();
	size_t q = fit.alpha.size(); // ARCH order
	size_t p = fit.beta.size();  // GARCH order

	// 1. Reconstruct historical residuals (epsilon_t = y_t - mu)
	double sum_data = std::accumulate(raw_data.begin(), raw_data.end(), 0.0);
	double mu = sum_data / static_cast<double>(N);
    
	vector<double> epsilon(N, 0.0);
	for (size_t t = 0; t < N; ++t) 
	{
		epsilon[t] = raw_data[t] - mu;
	}

	// 2. Re-filter the history to get conditional variances up to time T
	vector<double> sigma2_hist = computeGARCHVariances(epsilon, fit.omega, fit.alpha, fit.beta);

	// 3. Setup out-of-sample forecast tracking vectors
	vector<double> var_forecast(steps_to_forecast, 0.0);
	vector<double> vol_forecast(steps_to_forecast, 0.0);

	// Compute the long-run unconditional variance to use as an algorithmic fallback
	double sum_alpha = std::accumulate(fit.alpha.begin(), fit.alpha.end(), 0.0);
	double sum_beta = std::accumulate(fit.beta.begin(), fit.beta.end(), 0.0);
	double unconditional_var = fit.omega;
	if (1.0 - sum_alpha - sum_beta > 0.0) 
	{
		unconditional_var = fit.omega / (1.0 - sum_alpha - sum_beta);
	}

	// 4. Recursive out-of-sample projection loop
	for (size_t h = 0; h < steps_to_forecast; ++h) 
	{
		double current_forecast_var = fit.omega;

		// ARCH Lags (alpha_i * epsilon_{T+h-i}^2)
		for (size_t i = 1; i <= q; ++i) 
		{
			// Case A: The lag points inside our known historical timeline (t = N - i + h < N)
			if (h < i) 
			{
				size_t hist_idx = N - i + h;
				current_forecast_var += fit.alpha[i - 1] * epsilon[hist_idx] * epsilon[hist_idx];
			} 
			// Case B: The lag falls out-of-sample (E[epsilon_{T+h-i}^2] = sigma_{T+h-i}^2)
			else 
			{
				size_t forecast_idx = h - i;
				current_forecast_var += fit.alpha[i - 1] * var_forecast[forecast_idx];
			}
		}

		// GARCH Lags (beta_j * sigma_{T+h-j}^2)
		for (size_t j = 1; j <= p; ++j) 	
		{
			// Case A: The lag points inside known history
			if (h < j) 
			{
				size_t hist_idx = N - j + h;
				current_forecast_var += fit.beta[j - 1] * sigma2_hist[hist_idx];
			} 
			// Case B: The lag falls out-of-sample
			else 
			{
				size_t forecast_idx = h - j;
				current_forecast_var += fit.beta[j - 1] * var_forecast[forecast_idx];
			}
		}

		// Save variance and map to volatility scale
		var_forecast[h] = std::max(current_forecast_var, 1e-6);
		vol_forecast[h] = std::sqrt(var_forecast[h]);
	}

	GARCHForecast results;
	results.forecast_variance = std::move(var_forecast);
	results.forecast_volatility = std::move(vol_forecast);
	return results;
}

#endif
#endif