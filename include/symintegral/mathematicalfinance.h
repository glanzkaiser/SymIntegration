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


// dsolve.h

#ifndef SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE
#ifndef SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_DECLARE

double presentvalue(double, double, int);
double bondpricing(double, double, double, int);
double amortization(double, double, int);

#ifndef STRUCT_ARMAPARAMETERS
#define STRUCT_ARMAPARAMETERS

// Structure to hold the estimated ARMA parameters
struct ARMAParameters {
	double constant;
	vector<double> phi;
  	vector<double> theta;
	double sigma2; // Residual variance
	double log_likelihood;
};

#endif

vector<double> computeARMAResiduals(const std::vector<double>& data, double constant, const std::vector<double>& phi, const std::vector<double>& theta); 
double computeARMANegativeLogLikelihood(const std::vector<double>& params, const std::vector<double>& data, size_t p, size_t q);
vector<double> optimizeARMAMLE(const std::vector<double>& data, size_t p, size_t q, int max_iter = 1200);
ARMAParameters fitARMA(const vector<double>& raw_data, size_t p, size_t q);
vector<double> forecastARMA(const std::vector<double>& historical_data, double constant, const vector<double>& phi, const vector<double>& theta, size_t steps_to_forecast);

#ifndef STRUCT_GARCHPARAMETERS
#define STRUCT_GARCHPARAMETERS

// Structure to hold the estimated GARCH parameters
struct GARCHParameters {
	double omega;
	vector<double> alpha;
	vector<double> beta;
	double log_likelihood;
};

#endif

#ifndef STRUCT_GARCHFORECAST
#define STRUCT_GARCHFORECAST

// Structure to wrap out-of-sample GARCH forecast vectors
struct GARCHForecast {
	vector<double> forecast_variance;   // Projected sigma^2
	vector<double> forecast_volatility; // Projected sigma (standard deviation)
};
#endif

vector<double> computeGARCHVariances(const std::vector<double>& epsilon, double omega, const std::vector<double>& alpha, const std::vector<double>& beta);
double computeNegativeLogLikelihood(const std::vector<double>& params, const std::vector<double>& epsilon, size_t p, size_t q);
vector<double> optimizeNelderMead(const std::vector<double>& epsilon, size_t p, size_t q, int max_iter = 1000);
GARCHParameters fitGARCH(const std::vector<double>& raw_data, size_t p, size_t q) ;
GARCHForecast forecastGARCH(const std::vector<double>& raw_data, const GARCHParameters& fit, size_t steps_to_forecast);


#endif
#endif


#endif
