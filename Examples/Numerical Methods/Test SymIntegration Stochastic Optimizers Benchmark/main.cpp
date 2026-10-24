/*
    SymIntegration is branching from SymbolicC++ 3.35
    SymbolicC++ : An object oriented computer algebra system written in C++

    Copyright (C) 2008 Yorick Hardy and Willi-Hans Steeb

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/
// THANKS SENTINEL!!! and Freya too

// g++ -o result main.cpp -lsymintegration

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;



// ============================================================================
// MATHEMATICAL BENCHMARK OBJECTIVE FUNCTIONS
// ============================================================================

// 1. Rosenbrock Function: Global minimum at f(1,1,...,1) = 0
double rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size() - 1; ++i) {
        double t1 = x[i + 1] - x[i] * x[i];
        double t2 = 1.0 - x[i];
        sum += 100.0 * t1 * t1 + t2 * t2;
    }
    return sum;
}

// 2. Rastrigin Function: Global minimum at f(0,0,...,0) = 0
double rastrigin(const std::vector<double>& x) {
    double sum = 10.0 * x.size();
    for (double val : x) {
        sum += (val * val) - 10.0 * std::cos(2.0 * M_PI * val);
    }
    return sum;
}

// 3. Ackley Function: Global minimum at f(0,0,...,0) = 0
double ackley(const std::vector<double>& x) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (double val : x) {
        sum1 += val * val;
        sum2 += std::cos(2.0 * M_PI * val);
    }
    double n = static_cast<double>(x.size());
    return -20.0 * std::exp(-0.2 * std::sqrt(sum1 / n)) - std::exp(sum2 / n) + 20.0 + M_E;
}

// 4. Griewank Function: Global minimum at f(0,0,...,0) = 0
double griewank(const std::vector<double>& x) {
    double sum = 0.0;
    double prod = 1.0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += (x[i] * x[i]) / 4000.0;
        prod *= std::cos(x[i] / std::sqrt(static_cast<double>(i + 1)));
    }
    return sum - prod + 1.0;
}

// ============================================================================
// HELPER UTILITY FOR RESULTS LOGGER & CSV GENERATION
// ============================================================================
void process_benchmark_records(const std::string& benchmark_name, 
                               const std::vector<BenchmarkRecord>& records, 
                               int dimensions, 
                               int budget, 
                               std::ofstream& csv_file) {
    // --- Phase A: Terminal Output Presentation ---
    std::cout << "\n============================================================\n";
    std::cout << " BENCHMARK: " << benchmark_name << "\n";
    std::cout << "============================================================\n";
    std::left;
    std::cout << std::setw(23) << "Solver" 
              << std::setw(20) << "Best Energy" 
              << "Execution Time (ms)\n";
    std::cout << "------------------------------------------------------------\n";
    
    for (const auto& rec : records) {
        std::cout << std::left << std::setw(23) << rec.solver_name 
                  << std::setw(20) << std::setprecision(6) << std::scientific << rec.best_energy 
                  << std::fixed << std::setprecision(2) << rec.execution_time_ms << " ms\n";
        
        // --- Phase B: Dynamic Appending to CSV Stream Layer ---
        if (csv_file.is_open()) {
            csv_file << benchmark_name << ","
                     << rec.solver_name << ","
                     << dimensions << ","
                     << budget << ","
                     << std::scientific << std::setprecision(8) << rec.best_energy << ","
                     << std::fixed << std::setprecision(4) << rec.execution_time_ms << "\n";
        }
    }
}

// ============================================================================
// MAIN SYSTEM ORCHESTRATOR
// ============================================================================
int main() {
    // Structural Testing Framework Target Metrics
    const int dimensions = 5;
    const int evaluation_budget = 20000;
    const std::string csv_filename = "benchmark_results.csv";
    
    // Open dynamic file writer stream and generate structural table columns header
    std::ofstream csv_file(csv_filename);
    if (csv_file.is_open()) {
        csv_file << "Benchmark,Solver,Dimensions,EvaluationBudget,BestEnergy,ExecutionTimeMs\n";
    } else {
        std::cerr << "Warning: Unable to create or unlock targeted logging channel: " << csv_filename << "\n";
    }
    
    // Define an initial guessing state offset from global solutions
    std::vector<double> initial_base(dimensions, 2.0);
    
    // ------------------------------------------------------------------------
    // 1. Run Rosenbrock Suite
    // ------------------------------------------------------------------------
    std::vector<double> rosen_low(dimensions, -5.0);
    std::vector<double> rosen_high(dimensions, 5.0);
    auto rosen_results = run_unified_stochastic_benchmark(
        dimensions, rosen_low, rosen_high, initial_base, evaluation_budget, rosenbrock
    );
    process_benchmark_records("Rosenbrock", rosen_results, dimensions, evaluation_budget, csv_file);

    // ------------------------------------------------------------------------
    // 2. Run Rastrigin Suite
    // ------------------------------------------------------------------------
    std::vector<double> rast_low(dimensions, -5.12);
    std::vector<double> rast_high(dimensions, 5.12);
    auto rast_results = run_unified_stochastic_benchmark(
        dimensions, rast_low, rast_high, initial_base, evaluation_budget, rastrigin
    );
    process_benchmark_records("Rastrigin", rast_results, dimensions, evaluation_budget, csv_file);

    // ------------------------------------------------------------------------
    // 3. Run Ackley Suite
    // ------------------------------------------------------------------------
    std::vector<double> ackley_low(dimensions, -32.768);
    std::vector<double> ackley_high(dimensions, 32.768);
    auto ackley_results = run_unified_stochastic_benchmark(
        dimensions, ackley_low, ackley_high, initial_base, evaluation_budget, ackley
    );
    process_benchmark_records("Ackley", ackley_results, dimensions, evaluation_budget, csv_file);

    // ------------------------------------------------------------------------
    // 4. Run Griewank Suite
    // ------------------------------------------------------------------------
    std::vector<double> griewank_low(dimensions, -600.0);
    std::vector<double> griewank_high(dimensions, 600.0);
    auto griewank_results = run_unified_stochastic_benchmark(
        dimensions, griewank_low, griewank_high, initial_base, evaluation_budget, griewank
    );
    process_benchmark_records("Griewank", griewank_results, dimensions, evaluation_budget, csv_file);

    // Explicitly safely secure streams allocations closure processes
    if (csv_file.is_open()) {
        csv_file.close();
        std::cout << "\n>>> Dynamic log export process complete. Data saved to: " << csv_filename << "\n";
    }

    std::cout << "\n============================================================\n";
    std::cout << " Benchmark processing execution complete.\n";
    std::cout << "============================================================\n";

    return 0;
}
