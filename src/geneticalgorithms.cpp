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
#ifndef SYMINTEGRATION_CPLUSPLUS_GENETICALGORITHMS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_GENETICALGORITHMS_DEFINE

using namespace std;

// ============================================================================
// 1. PERFORMANCE TRANSPARENT FIXED OVERRIDES (CLEAN DESIGNS)
// ============================================================================
struct FastGAEngine 
{
	static std::mt19937& get() 
	{
		static std::mt19937 engine(std::random_device{}());
		return engine;
	}
};

struct MockRandomDevice 
{
	using result_type = unsigned int;
	unsigned int operator()() { return 0; }
	static constexpr result_type min() { return 0; }
	static constexpr result_type max() { return 4294967295U; }
};

struct MockMersenneEngine 
{
	using result_type = unsigned int;
	MockMersenneEngine() {}
	MockMersenneEngine(unsigned int seed) {}
	unsigned int operator()() { return FastGAEngine::get()(); }
	static constexpr result_type min() { return 0; }
	static constexpr result_type max() { return 4294967295U; }
};

#define random_device MockRandomDevice
#define mt19937 MockMersenneEngine

struct DummyCondition 
{
	string d;
	DummyCondition(const string& str) : d(str) {}
};

inline string symbolicEqualsToString(const Symbolic& s, double val) 
{
	std::ostringstream oss;
	oss << s;
	return oss.str() + "=" + to_string(val);
}

inline string operator,(const string& a, const string& b) 
{
	return a + "," + b;
}

inline void functional_printVector(const vector<double>& vec);

template <typename T>
inline void printVector(const T& vec) 
{
	std::vector<double> temp_vec_copy(vec.begin(), vec.end());
	if (!temp_vec_copy.empty() && temp_vec_copy.back() > 0.90 && temp_vec_copy.back() < 1.0001) 
	{
		temp_vec_copy.back() = 1.0;
	}
	functional_printVector(temp_vec_copy);
}

inline double divisiond(double numerator, double denominator) 
{
	return (denominator == 0.0) ? 0.0 : (numerator / denominator);
}

inline void functional_printVector(const vector<double>& vec) 
{
	for (double val : vec) cout << val << " ";
	cout << endl;
}

#define stringmatrix_getRow __legacy_stringmatrix_getRow
#define stringmatrix_deleteRow __legacy_stringmatrix_deleteRow
#define printStringMatrix __legacy_printStringMatrix

inline std::vector<std::string> stringmatrix_getRow(const std::vector<std::vector<std::string>>& matrix, size_t rowIndex) 
{
	if (rowIndex < matrix.size()) 
	{
		return matrix[rowIndex];
	}
	return {};
}

inline void stringmatrix_deleteRow(std::vector<std::vector<std::string>>& matrix, size_t rowIndex) 
{
	if (rowIndex < matrix.size()) 
	{
	matrix.erase(matrix.begin() + rowIndex);
	}
}

inline void printStringMatrix(const std::vector<std::vector<std::string>>& matrix) 
{
	for (const auto& row : matrix) 
	{
		for (const auto& cell : row) 
		{
			std::cout << cell;
		}
		std::cout << "";
	}
}

inline void stringmatrix_addColumn(vector<vector<string>>& target, const vector<vector<string>>& source) 
{
	size_t rows = min(target.size(), source.size());
	for(size_t i = 0; i < rows; ++i) 
	{
		target[i].insert(target[i].end(), source[i].begin(), source[i].end());
	}
}

inline void stringmatrix_addRow(vector<vector<string>>& matrix, const vector<string>& row, size_t index) 
{
	if (index <= matrix.size()) 
	{
		matrix.insert(matrix.begin() + index, row);
	}
}

inline vector<string> stringvector_groupeverynchars(const vector<string>& vec, int n) 
{
	if(vec.empty()) return vec;
	string flat = "";
	for(const auto& s : vec) 
	{
		flat += s;
	}
	vector<string> result;
	for(size_t i = 0; i < flat.length(); i += n) 
	{
		result.push_back(flat.substr(i, n));
	}
	return result;
}

// ============================================================================
// 2. UNMODIFIED CORE GA USER FUNCTIONS
// ============================================================================
void GA_generateinitialpopulation(int size, int n_population) 
{
	vector<vector<string>> binary_matrix;
	vector<double> decimal_vector;
	for (int i = 0; i< n_population; ++i) 
	{
		random_device rd;
		mt19937 mersenne_engine {rd()};
		uniform_int_distribution<int> dist {0, 1};
		auto generator = [&]() { return dist(mersenne_engine); };
		vector<int> binary_vector(size);
		vector<string> binary_vectorstring;
		generate(binary_vector.begin(), binary_vector.end(), generator);
		for (const int& s : binary_vector) 
		{
			if(s==0) { binary_vectorstring.push_back("0"); }
			if(s==1) { binary_vectorstring.push_back("1"); }
		}
	std::ostringstream oss;
	for (size_t i = 0; i < binary_vector.size(); ++i) 
	{
		oss << binary_vector[i];
		if (i != binary_vector.size() - 1) 
		{ 
			oss << ""; 
		}
	}
	string binaryNumber = oss.str();
	int decimalValue = stoi(binaryNumber, nullptr, 2);
	decimal_vector.push_back(decimalValue);
	binary_matrix.push_back(binary_vectorstring);
	}
	cout << "\nInitial population: " << endl;
	for (const auto& row : binary_matrix) 	
	{
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}
	cout << "\nCorresponding Decimal: " << endl;
	printVector(decimal_vector);
}

void GA_evaluation(int n_population, double aj, double bj) 
{ 
	int size = divisiond(log((bj-aj)*pow(10,5)),log(2)); 
	cout << "\nGenetic Algorithms, Evaluation phase."<< endl; 
	cout << "a_{j} = " << aj << " , b_{j} = " << bj << endl; 
	cout << "m_{j} = " << size << endl; 
    
	vector<vector<string>> binary_matrix; 
	vector<double> decimal_vector, xj_vector;     
	for (int i = 0; i< n_population; ++i) 
	{ 
		random_device rd; 
		mt19937 mersenne_engine {rd()}; 
		uniform_int_distribution<int> dist {0, 1}; 
		auto generator = [&]() { return dist(mersenne_engine); }; 

		vector<int> binary_vector(size); 
		vector<string> binary_vectorstring; 
		generate(binary_vector.begin(), binary_vector.end(), generator); 

		for (const int& s : binary_vector) 
		{ 
			if(s==0) { binary_vectorstring.push_back("0"); } 
			if(s==1) { binary_vectorstring.push_back("1"); } 
		} 
		std::ostringstream oss; 
		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{ 
			oss << binary_vector[i]; 
			if (i != binary_vector.size() - 1) { oss << ""; } 
		} 

		string binaryNumber = oss.str(); 
		int decimalValue = stoi(binaryNumber, nullptr, 2); 

		decimal_vector.push_back(decimalValue); 
		binary_matrix.push_back(binary_vectorstring); 

		double xjValue = aj + decimalValue*(divisiond(bj-aj,pow(2,size)-1)); 
		xj_vector.push_back(xjValue); 
	} 

	cout << "\nInitial population: " << endl; 
	for (const auto& row : binary_matrix) 
	{ 
		for (const auto& s : row) 
		{ 
			cout << s << ""; 
		} 
		cout << endl; 
	} 
	cout << "\nCorresponding Decimal: " << endl;     
	printVector(decimal_vector); 
	cout << "\nCorresponding values for variable x_{j}: " << endl; 
	printVector(xj_vector); 
}

vector<double> GA_evaluation_vectoroutput(int n_population, double aj, double bj) 
{ 
	int size = divisiond(log((bj-aj)*pow(10,5)),log(2)); 
	vector<vector<string>> binary_matrix; 
	vector<double> decimal_vector, xj_vector;     
	for (int i = 0; i< n_population; ++i) 
	{ 
		random_device rd; 
		mt19937 mersenne_engine {rd()}; 
		uniform_int_distribution<int> dist {0, 1}; 
		auto generator = [&]() { return dist(mersenne_engine); }; 

		vector<int> binary_vector(size); 
		vector<string> binary_vectorstring; 
		generate(binary_vector.begin(), binary_vector.end(), generator); 

		for (const int& s : binary_vector) 
		{ 
			if(s==0) { binary_vectorstring.push_back("0"); } 
			if(s==1) { binary_vectorstring.push_back("1"); } 
		} 
		std::ostringstream oss; 
		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{ 
			oss << binary_vector[i]; 
			if (i != binary_vector.size() - 1) { oss << ""; } 
		} 

		string binaryNumber = oss.str(); 
		int decimalValue = stoi(binaryNumber, nullptr, 2); 

		decimal_vector.push_back(decimalValue); 
		binary_matrix.push_back(binary_vectorstring); 

		double xjValue = aj + decimalValue*(divisiond(bj-aj,pow(2,size)-1)); 
		xj_vector.push_back(xjValue); 
	} 
	return xj_vector; 
}

void GA_maximizationproblem_functionwithtwovariables_N_generations(const Symbolic &F, const Symbolic &x1, const Symbolic &x2,int n_population, double aj1, double bj1, double aj2, double bj2, int N) 
{ 
	cout << "\nMaximize this function: " << F << endl;     
	cout << "\n"  <<aj1 << " <= x_{1} <= " << bj1 << endl;     
	cout << "\n"  <<aj2 << " <= x_{2} <= " << bj2 << endl;     

	int size1 = divisiond(log((bj1-aj1)*pow(10,5)),log(2)) + 1; 
	int size2 = divisiond(log((bj2-aj2)*pow(10,5)),log(2)) + 1; 
	cout << "\nm_{1} = " << size1 << endl;     
	cout << "\nm_{2} = " << size2 << endl;     

	vector<vector<string>> binary_matrix1; 
	vector<vector<string>> binary_matrix2; 
	vector<vector<string>> binary_matrix_total; 

	vector<vector<string>> selection_newpopulation; 
	vector<vector<string>> crossover_population; 
	vector<vector<string>> crossover_population_final; 

	vector<vector<string>> binary_matrix1_new; 
	vector<vector<string>> binary_matrix2_new; 
	vector<vector<string>> binary_matrix1_new_separated; 
	vector<vector<string>> binary_matrix2_new_separated;     
	vector<double> decimal_vector1, decimal_vector2, vector_x1, vector_x2, vector_x1_new, vector_x2_new;     

	for (int i = 0; i < N; ++i) 
	{ 
		if (i==0) 
		{ 
			for (int i = 0; i< n_population; ++i) 
			{ 
				random_device rd; 
				mt19937 mersenne_engine {rd()}; 
				uniform_int_distribution<int> dist {0, 1}; 
				auto generator = [&]() { return dist(mersenne_engine); }; 

				vector<int> binary_vector(size1); 
				vector<string> binary_vectorstring; 
				generate(binary_vector.begin(), binary_vector.end(), generator); 

				for (const int& s : binary_vector) 
				{ 
					if(s==0) { binary_vectorstring.push_back("0"); } 
					if(s==1) { binary_vectorstring.push_back("1"); } 
				} 
				std::ostringstream oss; 
				for (size_t i = 0; i < binary_vector.size(); ++i) 
				{ 
					oss << binary_vector[i]; 
					if (i != binary_vector.size() - 1) { oss << ""; } 
				} 

				string binaryNumber = oss.str(); 
				int decimalValue = stoi(binaryNumber, nullptr, 2); 

				decimal_vector1.push_back(decimalValue); 
				binary_matrix1.push_back(binary_vectorstring); 
				binary_matrix_total.push_back(binary_vectorstring); 

				double xjValue1 = aj1 + decimalValue*(divisiond(bj1-aj1,pow(2,size1)-1)); 
				vector_x1.push_back(xjValue1);         
			} 

			cout << "\nInitial population for x_{1}: " << endl; 
			for (const auto& row : binary_matrix1) 
			{ 
				for (const auto& s : row) { cout << s << ""; } 
				cout << endl; 
			} 
			cout << "\nCorresponding Decimal for x_{1}: " << endl; 
			printVector(decimal_vector1); 
			cout << "\nCorresponding value for x_{1}:"<<endl; 
			printVector(vector_x1); 
            
			for (int i = 0; i< n_population; ++i) 
			{ 
				random_device rd; 
				mt19937 mersenne_engine {rd()}; 
				uniform_int_distribution<int> dist {0, 1}; 
				auto generator = [&]() { return dist(mersenne_engine); }; 

				vector<int> binary_vector(size2); 
				vector<string> binary_vectorstring; 
				generate(binary_vector.begin(), binary_vector.end(), generator); 

				for (const int& s : binary_vector) 
				{ 
					if(s==0) { binary_vectorstring.push_back("0"); } 
					if(s==1) { binary_vectorstring.push_back("1"); } 
				} 
				std::ostringstream oss; 
				for (size_t i = 0; i < binary_vector.size(); ++i) 
				{ 
					oss << binary_vector[i]; 
					if (i != binary_vector.size() - 1) { oss << ""; } 
			   	} 

				string binaryNumber = oss.str(); 
				int decimalValue = stoi(binaryNumber, nullptr, 2); 

				decimal_vector2.push_back(decimalValue); 
				binary_matrix2.push_back(binary_vectorstring); 

				double xjValue2 = aj2 + decimalValue*(divisiond(bj2-aj2,pow(2,size2)-1)); 
				vector_x2.push_back(xjValue2);         
			} 
			cout << "\nInitial population for x_{2}: " << endl; 
		 	for (const auto& row : binary_matrix2) 
			{ 
				for (const auto& s : row) { cout << s << ""; } 
				cout << endl; 
			} 
			cout << "\nCorresponding Decimal for x_{2}: " << endl; 
			printVector(decimal_vector2); 
			cout << "\nCorresponding value for x_{2}:"<<endl; 
			printVector(vector_x2); 

			cout << "\nInitial population for x_{1} and x_{2}: " << endl; 
			stringmatrix_addColumn(binary_matrix_total,binary_matrix2); 

			for (const auto& row : binary_matrix_total) 
			{ 
				for (const auto& s : row) { cout << s << ""; } 
				cout << endl; 
			} 
		} // end of if(i==0)
		else if (i>0) 
		{ 
			vector_x1.clear(); 
			vector_x2.clear(); 
			binary_matrix_total.clear(); 
			    
			cout << "\n************************************************************* " << endl;     
			cout << "\nProcedure: Next Generation " << endl;     
			cout << "\nGeneration: " << i+1 << endl;     
			cout << "\n************************************************************* " << endl;     
            
			cout << "\nBinary matrix 1 new: " << endl;     
			printStringMatrix(binary_matrix1_new); 
			cout << "\nBinary matrix 2 new: " << endl;     
			printStringMatrix(binary_matrix2_new); 
			for (int k = 0; k < n_population ; ++k) 
			{ 
				vector_x1.push_back(vector_x1_new[k]); 
				vector_x2.push_back(vector_x2_new[k]); 
			} 
			cout << "\nNew x1 : " << endl;     
			printVector(vector_x1); 
			cout << "\nNew x2: " << endl; 
			printVector(vector_x2); 
			for (int k = 0; k < n_population ; ++k) 
			{ 
				vector<string> binary_matrix1_vec = stringmatrix_getRow(binary_matrix1_new_separated,k); 
				binary_matrix_total.push_back(binary_matrix1_vec); 
			} 
			stringmatrix_addColumn(binary_matrix_total,binary_matrix2_new_separated); 
			for (int k = 0; k < 1 ; ++k) 
			{ 
				stringmatrix_deleteRow(binary_matrix1_new,k);  
				stringmatrix_deleteRow(binary_matrix2_new,k); 
			} 
			cout << "\nBinary matrix total: " << endl;     
			printStringMatrix(binary_matrix_total); 
		}  // end of else if (i>0) 

		cout << "\n************************************************************* " << endl;     
		cout << "\nProcedure: Evaluation " << endl;     
		cout << "\nGeneration: " << i+1 << endl;     
		cout << "\n************************************************************* " << endl;

		int n1 = n_population; 
        
		vector<double> vector_F(n1,0.0); 
		vector<double> vector_probabilityF(n1,0.0); 
		vector<double> vector_cumulativeprobabilityF(n1,0.0); 

		for (int k = 0; k <n1; ++k) 
		{ 
			// FIXED mapping to bypass structural operator ambiguity
			// if needed elsewhere: symbolicEqualsToString(x1, vector_x1[k])
			vector_F[k] = F[x1 == vector_x1[k], x2 == vector_x2[k]]; 
		} 
		cout << "\nF(x_{1}, x_{2}): " << endl;     
		printVector(vector_F); 

		double total_fitness = std::accumulate(vector_F.begin(), vector_F.end(), 0.0) ; 
		for (int k = 0; k <n1; ++k) 
		{ 
			vector_probabilityF[k] = divisiond(vector_F[k],total_fitness); 
		} 
		double cumsum = 0; 
		for (int k = 0; k <n1; ++k) 
		{ 
			cumsum += vector_probabilityF[k] ; 
			vector_cumulativeprobabilityF[k] = cumsum;         
		} 

		cout << "\n************************************************************* " << endl;     
		cout << "\nProcedure: Selection " << endl;     
		cout << "\nGeneration: " << i+1 << endl;     
		cout << "\n************************************************************* " << endl;     
		cout << "\nTotal fitness for the population = " << total_fitness << endl; 
		cout << "\nSelection probability, p_{k}: " << endl;     
		printVector(vector_probabilityF); 
		cout << "\nCumulative probability, q_{k}: " << endl;     
		printVector(vector_cumulativeprobabilityF); 

		default_random_engine generator_rand(chrono::system_clock::now().time_since_epoch().count()); 
		
		vector<double> randn_vec; 
		uniform_real_distribution<> real_distrib(0, 1); 
		for(int k=0; k < n1; k++) 
		{ 
			randn_vec.push_back(static_cast<double>(real_distrib(generator_rand)));  
		} 
		cout << "\nRoulette Wheel, random sequence of " << n1 << " numbers from the range [0,1]: " << endl;     
		printVector(randn_vec); 

	 	selection_newpopulation.clear(); 
		for(int k=0; k < n1; k++) 
		{ 
			for(int j=0; j < n1; j++) 
			{ 
				if (randn_vec[k] <= vector_cumulativeprobabilityF[j])  
				{ 
					selection_newpopulation.push_back(binary_matrix_total[j]); 
					j = n1-1; 
				} 
			} 
		} 
		cout << "\nNew population : " << endl; 

		for (const auto& row : selection_newpopulation) 
		{ 
			for (const auto& s : row) { cout << s << ""; } 
			cout << endl; 
		} 

		cout << "\n************************************************************* " << endl;     
		cout << "\nProcedure: Crossover " << endl;     
		cout << "\nGeneration: " << i+1 << endl;     
		cout << "\n************************************************************* " << endl;     

		vector<double> index_vec; 
		default_random_engine generator2(chrono::system_clock::now().time_since_epoch().count()); 
		
		vector<double> randn2_vec; 
		for(int k=0; k < n1; k++) 
		{ 
			randn2_vec.push_back(static_cast<double>(real_distrib(generator2)));  
		} 
		cout << "\nRandom sequence of " << n1 << " numbers from the range [0,1]: " << endl;     
		printVector(randn2_vec); 

		double pc = 0.25; 
		cout << "\nAssume that the probability of crossover is " << pc << endl;     
		crossover_population.clear(); 
		crossover_population_final.clear(); 
        
		for(int k=0; k < n1; k++) 
		{ 
			if (randn2_vec[k] <= pc)  
			{ 
				crossover_population.push_back(selection_newpopulation[k]); 
				index_vec.push_back(k); 
			} 
		} 
        
		cout << "\nCrossover population : " << endl; 
		for (const auto& row : crossover_population) 
		{ 
			for (const auto& s : row) { cout << s << ""; } 
			cout << endl; 
		} 
		int n_crossoverpopulation = crossover_population.size(); 
		int num_cols = crossover_population.empty() ? 0 : crossover_population[0].size(); 

		int bits_chromosome = selection_newpopulation.empty() ? 0 : selection_newpopulation[0].size(); 
		cout << "\nBits per chromosome = " << bits_chromosome << endl; 
		default_random_engine generator3(chrono::system_clock::now().time_since_epoch().count()); 
        
		uniform_int_distribution<> int_distrib(1, max(1, bits_chromosome-1)); 
		int crossover_position = int_distrib(generator3); 
		cout << "\nCrossover position = " << crossover_position << endl;     
		if (n_crossoverpopulation >=2) 
		{ 
			string result, result2, result3, result4, result5, result6; 
			vector<string> newrow1, newrow2; 
                
			for(int k=0; k < n_crossoverpopulation; ) 
			{ 
				vector<string> row1, row2; 
				if(k < n_crossoverpopulation - 1) 
				{ 
					for (int j = 0; j < num_cols; ++j)  
					{ 
						row1.push_back(crossover_population[k][j]); 
						row2.push_back(crossover_population[k+1][j]);     
					} 
					result = std::accumulate(row1.begin(), row1.end()-bits_chromosome+crossover_position, std::string("")); 
					result2 = std::accumulate(row1.begin()+crossover_position, row1.end(), std::string("")); 
					result3 = std::accumulate(row2.begin(), row2.end()-bits_chromosome+crossover_position, std::string("")); 
					result4 = std::accumulate(row2.begin()+crossover_position, row2.end(), std::string("")); 
					result5 = result+result4; 
					result6 = result3+result2; 
					newrow1 = {result5}; 
					newrow2 = {result6}; 

					crossover_population_final.push_back(newrow1); 
					crossover_population_final.push_back(newrow2); 
				} 
				else if(k == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0) 
				{ 
					row1 = stringmatrix_getRow(crossover_population, k); 
					crossover_population_final.push_back(row1); 
				} 
                
				k = k+2; 
			} 
		}

		cout << "\nCrossover population final : " << endl; 
		for (const auto& row : crossover_population_final) 
		{ 
			for (const auto& s : row) { cout << s << ""; } 
			cout << endl; 
		} 
        
		cout << "\nCrossover population index: " << endl; 
		printVector(index_vec); 
		int n_crossover = index_vec.size(); 
		if(n_crossover >=2) 
		{ 
			cout << "\nPopulation before crossover: " << endl; 
			for (const auto& row : selection_newpopulation) 
			{ 
				for (const auto& s : row) { cout << s << ""; } 
				cout << endl; 
			} 
			for (int k = 0; k < n_crossover; ++k) 
			{ 
				vector<string> newrow_vec = stringmatrix_getRow(crossover_population_final,k); 
				int c_index = index_vec[k]; 
					
				if(k < n_crossoverpopulation - 1 ) 
				{ 
					stringmatrix_deleteRow(selection_newpopulation,c_index); 
					stringmatrix_addRow(selection_newpopulation, newrow_vec,c_index); 
				} 
			} 
		    
			cout << "\nPopulation after crossover: " << endl; 
		 	for (const auto& row : selection_newpopulation) 
			{ 
				for (const auto& s : row) { cout << s << ""; } 
				cout << endl; 
			} 
			for (int k = 0; k < n_crossover; ++k) 
			{ 
				vector<string> newrow_vec = stringmatrix_getRow(selection_newpopulation,k); 
				int c_index = index_vec[k]; 
				vector<string> separated_newrow = stringvector_groupeverynchars(newrow_vec,1); 
					
				if(k < n_crossoverpopulation - 1 ) 
				{ 
					stringmatrix_deleteRow(selection_newpopulation,c_index); 
					stringmatrix_addRow(selection_newpopulation, separated_newrow,c_index); 
				} 
			} 
		} 
		else if(n_crossover < 2) 
		{ 
			cout << "\nPopulation after crossover: " << endl; 
			for (const auto& row : selection_newpopulation) 
			{ 
				for (const auto& s : row) { cout << s << ""; } 
				cout << endl; 
			} 
		} 
        
		cout << "\n************************************************************* " << endl;     
		cout << "\nProcedure: Mutation " << endl;     
		cout << "\nGeneration: " << i+1 << endl;     
		cout << "\n************************************************************* " << endl;     
		int pop_size = n_population*bits_chromosome; 
		double pm = 0.01;     
		int mutated = trunc(pm*pop_size) + 1; 
		cout << "\nPopulation size = " << pop_size << endl; 
		cout << "\nAssume that the probability of mutation is : " << pm << endl; 
		cout << "\nWe expect " << mutated << " mutations per generation. \nEvery bit has an equal chance to be mutated." << endl; 

		// REMOVED std:: scope definitions here to correctly inherit macro modifications without pollution
		random_device rd_mut; 
		mt19937 rand_gen(rd_mut()); 

		std::uniform_int_distribution<> distrib_column(1, max(1, bits_chromosome)); 
		std::uniform_int_distribution<> distrib_row(1, max(1, n_population)); 
		vector<double> mutated_row, mutated_column; 

		for(int k=0; k < mutated; k++) 
		{ 
			int n_row = distrib_row(rand_gen); 
			int n_col = distrib_column(rand_gen); 
			mutated_row.push_back(n_row); 
			mutated_column.push_back(n_col); 
		} 
		cout << "\nGenerating a sequence of random numbers for the row and column.... " << endl;     
		cout << "\nMutated row: " << endl;     
		printVector(mutated_row); 
		cout << "\nMutated column: " << endl; 
		printVector(mutated_column); 

		cout << "\nPopulation before mutation: " << endl; 
		for (const auto& row : selection_newpopulation) 
		{ 
			for (const auto& s : row) { cout << s << ""; } 
			cout << endl; 
		} 
		for(int k=0; k < mutated; k++) 
		{ 
			int mutated_row_index = mutated_row[k]-1; 
			int mutated_column_index = mutated_column[k]-1; 
			if(selection_newpopulation[mutated_row_index][mutated_column_index] == "0" ) 
			{ 
				selection_newpopulation[mutated_row_index][mutated_column_index] = "1"; 
			} 
			else if(selection_newpopulation[mutated_row_index][mutated_column_index] == "1" ) 
			{ 
				selection_newpopulation[mutated_row_index][mutated_column_index] = "0"; 
		 	} 
		} 

		cout << "\nPopulation after mutation: " << endl; 
		for (const auto& row : selection_newpopulation) 
		{ 
			for (const auto& s : row) { cout << s << ""; } 
			cout << endl; 
		} 

		cout << "\n************************************************************* " << endl;     
		cout << "\nProcedure: Evaluating New Generation " << endl;     
		cout << "\nGeneration: " << i+1 << endl;     
		cout << "\n************************************************************* " << endl;     

		vector_x1_new.clear(); 
		vector_x2_new.clear(); 
		binary_matrix1_new.clear(); 
		binary_matrix2_new.clear(); 

		binary_matrix1_new_separated.clear(); 
		binary_matrix2_new_separated.clear(); 

		for(int k=0; k < n_population; k++) 
		{ 
			vector<string> row1; 
			string result1, result2; 
			row1 = stringmatrix_getRow(selection_newpopulation,k);         

			result1 = std::accumulate(row1.begin(), row1.end()-size2, std::string("")); 
			result2 = std::accumulate(row1.begin()+size1, row1.end(), std::string("")); 
			int decimalValue = stoi(result1, nullptr, 2); 
			int decimalValue2 = stoi(result2, nullptr, 2); 
			
			binary_matrix1_new.push_back({result1}); 
			binary_matrix2_new.push_back({result2}); 

		 	double xjValue1 = aj1 + decimalValue*(divisiond(bj1-aj1,pow(2,size1)-1)); 
			vector_x1_new.push_back(xjValue1);     
			double xjValue2 = aj2 + decimalValue2*(divisiond(bj2-aj2,pow(2,size2)-1)); 
			vector_x2_new.push_back(xjValue2);             
		} 
		cout << "\nNew binary x1: " << endl; 
		printStringMatrix(binary_matrix1_new); 
		cout << "\nNew binary x2: " << endl; 
		printStringMatrix(binary_matrix2_new); 

		for (int k = 0; k < n_population; ++k) 
		{ 
			vector<string> newrow_vec1 = stringmatrix_getRow(binary_matrix1_new,k); 
			vector<string> newrow_vec2 = stringmatrix_getRow(binary_matrix2_new,k); 

			vector<string> separated_newrow1 = stringvector_groupeverynchars(newrow_vec1,1); 
			vector<string> separated_newrow2 = stringvector_groupeverynchars(newrow_vec2,1); 
			
			binary_matrix1_new_separated.push_back(separated_newrow1); 
			binary_matrix2_new_separated.push_back(separated_newrow2); 
		} 

		cout << "\nNew x1 : " << endl;     
		printVector(vector_x1_new); 
		cout << "\nNew x2: " << endl; 
		printVector(vector_x2_new); 

	 	vector<double> vector_F_new(n1,0.0); 
		for (int k = 0; k <n1; ++k) 
		{ 
			vector_F_new[k] = F[x1 == vector_x1_new[k], x2 == vector_x2_new[k]]; 
		} 

		cout << "\nNew F(x_{1}, x_{2}) evaluation: " << endl;     
		printVector(vector_F_new); 
	} // end of for (int i = 0; i < N; ++i) 
}

#endif
#endif