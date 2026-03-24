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

void GA_generateinitialpopulation(int size, int n_population)
{
	vector<vector<string>> binary_matrix;
	vector<double> decimal_vector;	
	for (int i = 0; i< n_population; ++i)
	{
		// 1. Seed the random number engine
		random_device rd; 
		// Use a Mersenne Twister engine for high quality pseudo-random numbers
		mt19937 mersenne_engine {rd()}; 

		// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
		uniform_int_distribution<int> dist {0, 1};

		// 3. Create a generator function object (lambda function) that binds the distribution and engine
		auto generator = [&]() { return dist(mersenne_engine); };

		// 4. Create the vector and fill it using std::generate
		vector<int> binary_vector(size);
		vector<string> binary_vectorstring;
		generate(binary_vector.begin(), binary_vector.end(), generator);

		// 5. Print the generated vector (optional)
		/*
		cout << "Generated binary vector: ";
		for (int bit : binary_vector) 
		{
			cout << bit;
		}
		cout << endl;*/
		
		for (const int& s : binary_vector) 
		{
			if(s==0)
			{
				binary_vectorstring.push_back("0");
			}
			if(s==1)
			{
				binary_vectorstring.push_back("1");
			}
		}
		// Convert vector<int> into string
		// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
		std::ostringstream oss;

		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{
		oss << binary_vector[i];
		if (i != binary_vector.size() - 1) 
		{ // Add a space if it's not the last element
			oss << "";
		}
		}

		// Accumulate messenger_RNA as one string with no space
		string binaryNumber = oss.str();
		// string binaryNumber = binary_matrix[i][0]; // Example binary input
		int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

		//cout << "Binary: " << binaryNumber << endl;
		//cout << "Decimal: " << decimalValue << endl;

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
	// for precision of five places after the decimal point
	int size = divisiond(log((bj-aj)*pow(10,5)),log(2));
	cout << "\nGenetic Algorithms, Evaluation phase."<< endl;
	cout << "a_{j} = " << aj << " , b_{j} = " << bj << endl;
	cout << "m_{j} = " << size << endl;
	
	vector<vector<string>> binary_matrix;
	vector<double> decimal_vector, xj_vector;	
	for (int i = 0; i< n_population; ++i)
	{
		// 1. Seed the random number engine
		random_device rd; 
		// Use a Mersenne Twister engine for high quality pseudo-random numbers
		mt19937 mersenne_engine {rd()}; 

		// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
		uniform_int_distribution<int> dist {0, 1};

		// 3. Create a generator function object (lambda function) that binds the distribution and engine
		auto generator = [&]() { return dist(mersenne_engine); };

		// 4. Create the vector and fill it using std::generate
		vector<int> binary_vector(size);
		vector<string> binary_vectorstring;
		generate(binary_vector.begin(), binary_vector.end(), generator);

		// 5. Print the generated vector (optional)
		/*
		cout << "Generated binary vector: ";
		for (int bit : binary_vector) 
		{
			cout << bit;
		}
		cout << endl;*/
		
		for (const int& s : binary_vector) 
		{
			if(s==0)
			{
				binary_vectorstring.push_back("0");
			}
			if(s==1)
			{
				binary_vectorstring.push_back("1");
			}
		}
		// Convert vector<int> into string
		// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
		std::ostringstream oss;

		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{
		oss << binary_vector[i];
		if (i != binary_vector.size() - 1) 
		{ // Add a space if it's not the last element
			oss << "";
		}
		}

		// Accumulate messenger_RNA as one string with no space
		string binaryNumber = oss.str();
		// string binaryNumber = binary_matrix[i][0]; // Example binary input
		int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

		//cout << "Binary: " << binaryNumber << endl;
		//cout << "Decimal: " << decimalValue << endl;

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

void GA_maximizationproblem_functionwithtwovariables(const Symbolic &F, const Symbolic &x1, const Symbolic &x2,int n_population, double aj1, double bj1, double aj2, double bj2)
{
	cout << "\nMaximize this function: " << F << endl;	
	cout << "\n"  <<aj1 << " <= x_{1} <= " << bj1 << endl;	
	cout << "\n"  <<aj2 << " <= x_{2} <= " << bj2 << endl;	

	int size1 = divisiond(log((bj1-aj1)*pow(10,5)),log(2));
	int size2 = divisiond(log((bj2-aj2)*pow(10,5)),log(2));
	
	vector<vector<string>> binary_matrix1;
	vector<vector<string>> binary_matrix2;
	vector<vector<string>> binary_matrix_total;
	vector<vector<string>> selection_newpopulation;
	vector<vector<string>> crossover_population;
	vector<vector<string>> crossover_population_final;
	
	vector<double> decimal_vector1, decimal_vector2, vector_x1, vector_x2, vector_x1_new, vector_x2_new;	
	for (int i = 0; i< n_population; ++i)
	{
		// 1. Seed the random number engine
		random_device rd; 
		// Use a Mersenne Twister engine for high quality pseudo-random numbers
		mt19937 mersenne_engine {rd()}; 

		// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
		uniform_int_distribution<int> dist {0, 1};

		// 3. Create a generator function object (lambda function) that binds the distribution and engine
		auto generator = [&]() { return dist(mersenne_engine); };

		// 4. Create the vector and fill it using std::generate
		vector<int> binary_vector(size1);
		vector<string> binary_vectorstring;
		generate(binary_vector.begin(), binary_vector.end(), generator);

		// 5. Print the generated vector (optional)
		/*
		cout << "Generated binary vector: ";
		for (int bit : binary_vector) 
		{
			cout << bit;
		}
		cout << endl;*/
		
		for (const int& s : binary_vector) 
		{
			if(s==0)
			{
				binary_vectorstring.push_back("0");
			}
			if(s==1)
			{
				binary_vectorstring.push_back("1");
			}
		}
		// Convert vector<int> into string
		// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
		std::ostringstream oss;

		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{
		oss << binary_vector[i];
		if (i != binary_vector.size() - 1) 
		{ // Add a space if it's not the last element
			oss << "";
		}
		}

		string binaryNumber = oss.str();
		// string binaryNumber = binary_matrix[i][0]; // Example binary input
		int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

		//cout << "Binary: " << binaryNumber << endl;
		//cout << "Decimal: " << decimalValue << endl;

		decimal_vector1.push_back(decimalValue);
		binary_matrix1.push_back(binary_vectorstring);
		binary_matrix_total.push_back(binary_vectorstring);

		double xjValue1 = aj1 + decimalValue*(divisiond(bj1-aj1,pow(2,size1)-1));
		vector_x1.push_back(xjValue1);		
	}

	cout << "\nInitial population for x_{1}: " << endl;
	
	for (const auto& row : binary_matrix1) 
	{
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}
	cout << "\nCorresponding Decimal for x_{1}: " << endl;
	printVector(decimal_vector1);
	cout << "\nCorresponding value for x_{1}:"<<endl;
	printVector(vector_x1);
	for (int i = 0; i< n_population; ++i)
	{
		// 1. Seed the random number engine
		random_device rd; 
		// Use a Mersenne Twister engine for high quality pseudo-random numbers
		mt19937 mersenne_engine {rd()}; 

		// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
		uniform_int_distribution<int> dist {0, 1};

		// 3. Create a generator function object (lambda function) that binds the distribution and engine
		auto generator = [&]() { return dist(mersenne_engine); };

		// 4. Create the vector and fill it using std::generate
		vector<int> binary_vector(size2);
		vector<string> binary_vectorstring;
		generate(binary_vector.begin(), binary_vector.end(), generator);

		// 5. Print the generated vector (optional)
		/*
		cout << "Generated binary vector: ";
		for (int bit : binary_vector) 
		{
			cout << bit;
		}
		cout << endl;*/
		
		for (const int& s : binary_vector) 
		{
			if(s==0)
			{
				binary_vectorstring.push_back("0");
			}
			if(s==1)
			{
				binary_vectorstring.push_back("1");
			}
		}
		// Convert vector<int> into string
		// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
		std::ostringstream oss;

		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{
		oss << binary_vector[i];
		if (i != binary_vector.size() - 1) 
		{ // Add a space if it's not the last element
			oss << "";
		}
		}

		string binaryNumber = oss.str();
		// string binaryNumber = binary_matrix[i][0]; // Example binary input
		int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

		//cout << "Binary: " << binaryNumber << endl;
		//cout << "Decimal: " << decimalValue << endl;

		decimal_vector2.push_back(decimalValue);
		binary_matrix2.push_back(binary_vectorstring);

		double xjValue2 = aj2 + decimalValue*(divisiond(bj2-aj2,pow(2,size2)-1));
		vector_x2.push_back(xjValue2);		
	}
	cout << "\nInitial population for x_{2}: " << endl;
	for (const auto& row : binary_matrix2) 
	{
		for (const auto& s : row) 
		{
			cout << s << "";
		}
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
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}

	cout << "\n************************************************************* " << endl;	
	cout << "\nProcedure: Evaluation " << endl;	
	cout << "\n************************************************************* " << endl;	
	
	int n1 = vector_x1.size();
	
	vector<double> vector_F(n1,0.0);
	vector<double> vector_probabilityF(n1,0.0);
	vector<double> vector_cumulativeprobabilityF(n1,0.0);

	for (int i = 0; i <n1; ++i)
	{
		vector_F[i] = F[x1 == vector_x1[i], x2 == vector_x2[i]];
	}
	cout << "\nF(x_{1}, x_{2}): " << endl;	
	printVector(vector_F);

	double total_fitness = std::accumulate(vector_F.begin(), vector_F.end(), 0.0) ;
	for (int i = 0; i <n1; ++i)
	{
		vector_probabilityF[i] = divisiond(vector_F[i],total_fitness);
	}
	double cumsum = 0;
	for (int i = 0; i <n1; ++i)
	{
		cumsum += vector_probabilityF[i] ;
		vector_cumulativeprobabilityF[i] = cumsum;		
	}
	cout << "\n************************************************************* " << endl;	
	cout << "\nProcedure: Selection " << endl;	
	cout << "\n************************************************************* " << endl;	
	cout << "\nTotal fitness for the population = " << total_fitness << endl;
	cout << "\nSelection probability, p_{k}: " << endl;	
	printVector(vector_probabilityF);
	cout << "\nCumulative probability, q_{k}: " << endl;	
	printVector(vector_cumulativeprobabilityF);

	// Create random number from the range [0,1] 
	default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
	
	vector<double> randn_vec;
	uniform_real_distribution<> real_distrib(0, 1);
	for(int i=0; i < n1; i++)
	{
		randn_vec.push_back(static_cast<double>(real_distrib(generator))); 
	}
	cout << "\nRoulette Wheel, random sequence of " << n1 << " numbers from the range [0,1]: " << endl;	
	printVector(randn_vec);

	for(int i=0; i < n1; i++)
	{
		for(int j=0; j < n1; j++)
		{
			if  (randn_vec[i] <= vector_cumulativeprobabilityF[j]) 
			{
				selection_newpopulation.push_back(binary_matrix_total[j]);
				j = n1-1;
			}
			else if  (randn_vec[i] > vector_cumulativeprobabilityF[j]) 
			{
				
			}
		}
	}
	cout << "\nNew population : " << endl;

	for (const auto& row : selection_newpopulation) 
	{
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}
	cout << "\n************************************************************* " << endl;	
	cout << "\nProcedure: Crossover " << endl;	
	cout << "\n************************************************************* " << endl;	

	vector<double> index_vec;
	// Create random number from the range [0,1] for crossover
	default_random_engine generator2(chrono::system_clock::now().time_since_epoch().count());
	
	vector<double> randn2_vec;
	for(int i=0; i < n1; i++)
	{
		randn2_vec.push_back(static_cast<double>(real_distrib(generator2))); 
	}
	cout << "\nRandom sequence of " << n1 << " numbers from the range [0,1]: " << endl;	
	printVector(randn2_vec);

	double pc = 0.25;
	cout << "\nAssume that the probability of crossover is " << pc << endl;	
	
	for(int i=0; i < n1; i++)
	{
		if  (randn2_vec[i] <= pc) 
		{
			crossover_population.push_back(selection_newpopulation[i]);
			index_vec.push_back(i);
		}
		else if  (randn2_vec[i] > pc) 
		{

		}
	}
	
	cout << "\nCrossover population : " << endl;

	for (const auto& row : crossover_population) 
	{
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}
	int n_crossoverpopulation = crossover_population.size();
	int num_cols = crossover_population[0].size();

	int bits_chromosome = selection_newpopulation[0].size();
	cout << "\nBits per chromosome = " << bits_chromosome << endl;
	// Create random int from the range [1,bits_chromosome-1] 
	default_random_engine generator3(chrono::system_clock::now().time_since_epoch().count());
	
	uniform_int_distribution<> int_distrib(1, bits_chromosome-1);
	int crossover_position = int_distrib(generator3); 
	cout << "\nCrossover position = " << crossover_position << endl;	
	if (n_crossoverpopulation >=2)
	{
		string result, result2, result3, result4, result5, result6;
		vector<string> newrow1, newrow2;
			
		for(int i=0; i < n_crossoverpopulation; )
		{
			//cout << "\ni = " << i << endl;	
	
			vector<string> row1, row2;
			if(i < n_crossoverpopulation - 1)
			{
				for (int j = 0; j < num_cols; ++j) 
				{
					row1.push_back(crossover_population[i][j]);
					row2.push_back(crossover_population[i+1][j]);	
				}
				result = std::accumulate(row1.begin(), row1.end()-bits_chromosome+crossover_position, std::string(""));
				//cout << "\ncutted result row 1 left= " << endl;
				//cout << result << endl;
				result2 = std::accumulate(row1.begin()+crossover_position, row1.end(), std::string(""));
				//cout << "\ncutted result row 1 right= " << endl;
				//cout << result2 << endl;
				result3 = std::accumulate(row2.begin(), row2.end()-bits_chromosome+crossover_position, std::string(""));
				//cout << "\ncutted result row 2 left = " << endl;
				//cout << result3 << endl;
				result4 = std::accumulate(row2.begin()+crossover_position, row2.end(), std::string(""));
				//cout << "\ncutted result row 2 right= " << endl;
				//cout << result4 << endl;
				result5 = result+result4;
				//cout << "\ncombined result row 1 left +  row 2 right= " << endl;
				//cout << result5 << endl;
				result6 = result3+result2;
				//cout << "\ncombined result row 2 left +  row 1 right= " << endl;
				//cout << result6 << endl;
				newrow1 = {result5};
				newrow2 = {result6};

				crossover_population_final.push_back(newrow1);
				crossover_population_final.push_back(newrow2);
			}
			else if(i == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0)
			{
				row1 = stringmatrix_getRow(crossover_population, i);
				crossover_population_final.push_back(row1);
			}
			
			i = i+2;
		}
	}
	else if(n_crossoverpopulation < 2)
	{

	}
	
	cout << "\nCrossover population final : " << endl;

	for (const auto& row : crossover_population_final) 
	{
		for (const auto& s : row) 
		{
			cout << s << "";
		}
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
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
		vector<string> newrow_vec;
		for (int i = 0; i < n_crossover; ++i)
		{
			newrow_vec = stringmatrix_getRow(crossover_population_final,i);
			int c_index = index_vec[i];
			
			if(i < n_crossoverpopulation - 1 )
			{
				stringmatrix_deleteRow(selection_newpopulation,c_index);
				stringmatrix_addRow(selection_newpopulation, newrow_vec,c_index);
			}
			else if(i == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0)
			{

			}

		}
		
		cout << "\nPopulation after crossover: " << endl;
		for (const auto& row : selection_newpopulation) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
		// The rows that have undergone crossover become 1 string that represent 37 bits (the chromosome bits), while it should be 37 strings instead of 1.
		// This is the separation of 1 string of 101010100101010...101001 into 37 strings to represent the chromosome bits.
		for (int i = 0; i < n_crossover; ++i)
		{
			newrow_vec = stringmatrix_getRow(selection_newpopulation,i);
			int c_index = index_vec[i];
			vector<string> separated_newrow = stringvector_groupeverynchars(newrow_vec,1);
			
			if(i < n_crossoverpopulation - 1 )
			{
				stringmatrix_deleteRow(selection_newpopulation,c_index);
				stringmatrix_addRow(selection_newpopulation, separated_newrow,c_index);
			}
			else if(i == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0)
			{

			}

		}

		//cout << "\nPopulation after crossover and refashioned: " << endl;
		//printStringMatrix(selection_newpopulation);

	}
	else if(n_crossover < 2)
	{
		cout << "\nPopulation after crossover: " << endl;

		for (const auto& row : selection_newpopulation) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
	}
	cout << "\n************************************************************* " << endl;	
	cout << "\nProcedure: Mutation " << endl;	
	cout << "\n************************************************************* " << endl;	
	int pop_size = n_population*bits_chromosome;
	double pm = 0.01;	
	int mutated = trunc(pm*pop_size) + 1;
	cout << "\nPopulation size = " << pop_size << endl;
	cout << "\nAssume that the probability of mutation is : " << pm << endl;
	cout << "\nWe expect " << mutated << " mutations per generation. \nEvery bit has an equal chance to be mutated." << endl;

	std::random_device rd;
	std::mt19937 generate(rd()); // Mersenne Twister engine seeded with rd

	std::uniform_int_distribution<> distrib_column(1, bits_chromosome);
	std::uniform_int_distribution<> distrib_row(1, n_population);
	vector<double> mutated_row, mutated_column;

	for(int i=0; i < mutated; i++)
	{
		int n_row = distrib_row(generate);
		int n_col = distrib_column(generate);
		
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
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}
	for(int i=0; i < mutated; i++)
	{
		int mutated_row_index = mutated_row[i]-1;
		int mutated_column_index = mutated_column[i]-1;
		//cout << "mutated entry= " << selection_newpopulation[mutated_row_index][mutated_column_index] << endl;
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
		for (const auto& s : row) 
		{
			cout << s << "";
		}
		cout << endl;
	}

	cout << "\n************************************************************* " << endl;	
	cout << "\nProcedure: Evaluating New Generation " << endl;	
	cout << "\n************************************************************* " << endl;	
	for(int i=0; i < n_population; i++)
	{
		//cout << "\ni = " << i << endl;
		vector<string> row1;
		string result1, result2;
		row1 = stringmatrix_getRow(selection_newpopulation,i);		

		result1 = std::accumulate(row1.begin(), row1.end()-size1, std::string(""));
		//cout << "\nsplit left= " << endl;
		//cout << result1 << endl;
		result2 = std::accumulate(row1.begin()+size2, row1.end(), std::string(""));
		//cout << "\nsplit right= " << endl;
		//cout << result2 << endl;
		int decimalValue = stoi(result1, nullptr, 2); // Convert string to int with base 2
		int decimalValue2 = stoi(result2, nullptr, 2); // Convert string to int with base 2

		//cout << "Decimal 1: " << decimalValue << endl;
		//cout << "Decimal 2: " << decimalValue2 << endl;
	
		double xjValue1 = aj1 + decimalValue*(divisiond(bj1-aj1,pow(2,size1)-1));
		vector_x1_new.push_back(xjValue1);	
		double xjValue2 = aj2 + decimalValue2*(divisiond(bj2-aj2,pow(2,size2)-1));
		vector_x2_new.push_back(xjValue2);			
	}
	cout << "\nNew x1 : " << endl;	
	printVector(vector_x1_new);
	cout << "\nNew x2: " << endl;
	printVector(vector_x2_new);

	vector<double> vector_F_new(n1,0.0);

	for (int i = 0; i <n1; ++i)
	{
		vector_F_new[i] = F[x1 == vector_x1_new[i], x2 == vector_x2_new[i]];
	}
	cout << "\nF(x_{1}, x_{2}) new generation: " << endl;	
	printVector(vector_F_new);

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

	for (int i = 0; i< N; ++i)
	{
	if (i==0)
	{
		for (int i = 0; i< n_population; ++i)
		{
			// 1. Seed the random number engine
			random_device rd; 
			// Use a Mersenne Twister engine for high quality pseudo-random numbers
			mt19937 mersenne_engine {rd()}; 

			// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
			uniform_int_distribution<int> dist {0, 1};

			// 3. Create a generator function object (lambda function) that binds the distribution and engine
			auto generator = [&]() { return dist(mersenne_engine); };

			// 4. Create the vector and fill it using std::generate
			vector<int> binary_vector(size1);
			vector<string> binary_vectorstring;
			generate(binary_vector.begin(), binary_vector.end(), generator);

			// 5. Print the generated vector (optional)
			/*
			cout << "Generated binary vector: ";
			for (int bit : binary_vector) 
			{
				cout << bit;
			}
			cout << endl;*/
			
			for (const int& s : binary_vector) 
			{
				if(s==0)
				{
					binary_vectorstring.push_back("0");
				}
				if(s==1)
				{
					binary_vectorstring.push_back("1");
				}
			}
			// Convert vector<int> into string
			// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
			std::ostringstream oss;

			for (size_t i = 0; i < binary_vector.size(); ++i) 
			{
			oss << binary_vector[i];
			if (i != binary_vector.size() - 1) 
			{ // Add a space if it's not the last element
				oss << "";
			}
			}

			string binaryNumber = oss.str();
			// string binaryNumber = binary_matrix[i][0]; // Example binary input
			int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

			//cout << "Binary: " << binaryNumber << endl;
			//cout << "Decimal: " << decimalValue << endl;

			decimal_vector1.push_back(decimalValue);
			binary_matrix1.push_back(binary_vectorstring);
			binary_matrix_total.push_back(binary_vectorstring);

			double xjValue1 = aj1 + decimalValue*(divisiond(bj1-aj1,pow(2,size1)-1));
			vector_x1.push_back(xjValue1);		
		}

		cout << "\nInitial population for x_{1}: " << endl;
		
		for (const auto& row : binary_matrix1) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
		cout << "\nCorresponding Decimal for x_{1}: " << endl;
		printVector(decimal_vector1);
		cout << "\nCorresponding value for x_{1}:"<<endl;
		printVector(vector_x1);
		for (int i = 0; i< n_population; ++i)
		{
			// 1. Seed the random number engine
			random_device rd; 
			// Use a Mersenne Twister engine for high quality pseudo-random numbers
			mt19937 mersenne_engine {rd()}; 

			// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
			uniform_int_distribution<int> dist {0, 1};

			// 3. Create a generator function object (lambda function) that binds the distribution and engine
			auto generator = [&]() { return dist(mersenne_engine); };

			// 4. Create the vector and fill it using std::generate
			vector<int> binary_vector(size2);
			vector<string> binary_vectorstring;
			generate(binary_vector.begin(), binary_vector.end(), generator);

			// 5. Print the generated vector (optional)
			/*
			cout << "Generated binary vector: ";
			for (int bit : binary_vector) 
			{
				cout << bit;
			}
			cout << endl;*/
			
			for (const int& s : binary_vector) 
			{
				if(s==0)
				{
					binary_vectorstring.push_back("0");
				}
				if(s==1)
				{
					binary_vectorstring.push_back("1");
				}
			}
			// Convert vector<int> into string
			// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
			std::ostringstream oss;

			for (size_t i = 0; i < binary_vector.size(); ++i) 
			{
			oss << binary_vector[i];
			if (i != binary_vector.size() - 1) 
			{ // Add a space if it's not the last element
				oss << "";
			}
			}

			string binaryNumber = oss.str();
			// string binaryNumber = binary_matrix[i][0]; // Example binary input
			int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

			//cout << "Binary: " << binaryNumber << endl;
			//cout << "Decimal: " << decimalValue << endl;

			decimal_vector2.push_back(decimalValue);
			binary_matrix2.push_back(binary_vectorstring);

			double xjValue2 = aj2 + decimalValue*(divisiond(bj2-aj2,pow(2,size2)-1));
			vector_x2.push_back(xjValue2);		
		}
		cout << "\nInitial population for x_{2}: " << endl;
		for (const auto& row : binary_matrix2) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
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
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
	}
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
		for (int i = 0; i < n_population ; ++i)
		{
			vector_x1.push_back(vector_x1_new[i]);
			vector_x2.push_back(vector_x2_new[i]);
		}
		cout << "\nNew x1 : " << endl;	
		printVector(vector_x1);
		cout << "\nNew x2: " << endl;
		printVector(vector_x2);
		for (int i = 0; i < n_population ; ++i)
		{
			vector<string> binary_matrix1_vec = stringmatrix_getRow(binary_matrix1_new_separated,i);
			binary_matrix_total.push_back(binary_matrix1_vec);
		}
		stringmatrix_addColumn(binary_matrix_total,binary_matrix2_new_separated);
		for (int i = 0; i < 1 ; ++i)
		{
			stringmatrix_deleteRow(binary_matrix1_new,i); // clear the new binary matrix obtained at the end of loop number i
			stringmatrix_deleteRow(binary_matrix2_new,i);
			//stringmatrix_deleteRow(binary_matrix1_new_separated,i); // clear the new binary matrix obtained at the end of loop number i
			//stringmatrix_deleteRow(binary_matrix2_new_separated,i);
			
		}
		cout << "\nBinary matrix total: " << endl;	
		printStringMatrix(binary_matrix_total);
	}
		cout << "\n************************************************************* " << endl;	
		cout << "\nProcedure: Evaluation " << endl;	
		cout << "\nGeneration: " << i+1 << endl;	
		cout << "\n************************************************************* " << endl;	
		
		int n1 = n_population;
		
		vector<double> vector_F(n1,0.0);
		vector<double> vector_probabilityF(n1,0.0);
		vector<double> vector_cumulativeprobabilityF(n1,0.0);

		for (int i = 0; i <n1; ++i)
		{
			vector_F[i] = F[x1 == vector_x1[i], x2 == vector_x2[i]];
		}
		cout << "\nF(x_{1}, x_{2}): " << endl;	
		printVector(vector_F);

		double total_fitness = std::accumulate(vector_F.begin(), vector_F.end(), 0.0) ;
		for (int i = 0; i <n1; ++i)
		{
			vector_probabilityF[i] = divisiond(vector_F[i],total_fitness);
		}
		double cumsum = 0;
		for (int i = 0; i <n1; ++i)
		{
			cumsum += vector_probabilityF[i] ;
			vector_cumulativeprobabilityF[i] = cumsum;		
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

		// Create random number from the range [0,1] 
		default_random_engine generator(chrono::system_clock::now().time_since_epoch().count());
		
		vector<double> randn_vec;
		uniform_real_distribution<> real_distrib(0, 1);
		for(int i=0; i < n1; i++)
		{
			randn_vec.push_back(static_cast<double>(real_distrib(generator))); 
		}
		cout << "\nRoulette Wheel, random sequence of " << n1 << " numbers from the range [0,1]: " << endl;	
		printVector(randn_vec);

		selection_newpopulation.clear();
		for(int i=0; i < n1; i++)
		{
			for(int j=0; j < n1; j++)
			{
				if  (randn_vec[i] <= vector_cumulativeprobabilityF[j]) 
				{
					selection_newpopulation.push_back(binary_matrix_total[j]);
					j = n1-1;
				}
				else if  (randn_vec[i] > vector_cumulativeprobabilityF[j]) 
				{
					
				}
			}
		}
		cout << "\nNew population : " << endl;

		for (const auto& row : selection_newpopulation) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
		cout << "\n************************************************************* " << endl;	
		cout << "\nProcedure: Crossover " << endl;	
		cout << "\nGeneration: " << i+1 << endl;	
		cout << "\n************************************************************* " << endl;	

		vector<double> index_vec;
		// Create random number from the range [0,1] for crossover
		default_random_engine generator2(chrono::system_clock::now().time_since_epoch().count());
		
		vector<double> randn2_vec;
		for(int i=0; i < n1; i++)
		{
			randn2_vec.push_back(static_cast<double>(real_distrib(generator2))); 
		}
		cout << "\nRandom sequence of " << n1 << " numbers from the range [0,1]: " << endl;	
		printVector(randn2_vec);

		double pc = 0.25;
		cout << "\nAssume that the probability of crossover is " << pc << endl;	
		crossover_population.clear();
		crossover_population_final.clear();
		
		for(int i=0; i < n1; i++)
		{
			if  (randn2_vec[i] <= pc) 
			{
				crossover_population.push_back(selection_newpopulation[i]);
				index_vec.push_back(i);
			}
			else if  (randn2_vec[i] > pc) 
			{

			}
		}
		
		cout << "\nCrossover population : " << endl;

		for (const auto& row : crossover_population) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
		int n_crossoverpopulation = crossover_population.size();
		int num_cols = crossover_population[0].size();

		int bits_chromosome = selection_newpopulation[0].size();
		cout << "\nBits per chromosome = " << bits_chromosome << endl;
		// Create random int from the range [1,bits_chromosome-1] 
		default_random_engine generator3(chrono::system_clock::now().time_since_epoch().count());
		
		uniform_int_distribution<> int_distrib(1, bits_chromosome-1);
		int crossover_position = int_distrib(generator3); 
		cout << "\nCrossover position = " << crossover_position << endl;	
		if (n_crossoverpopulation >=2)
		{
			string result, result2, result3, result4, result5, result6;
			vector<string> newrow1, newrow2;
				
			for(int i=0; i < n_crossoverpopulation; )
			{
				//cout << "\ni = " << i << endl;	
		
				vector<string> row1, row2;
				if(i < n_crossoverpopulation - 1)
				{
					for (int j = 0; j < num_cols; ++j) 
					{
						row1.push_back(crossover_population[i][j]);
						row2.push_back(crossover_population[i+1][j]);	
					}
					result = std::accumulate(row1.begin(), row1.end()-bits_chromosome+crossover_position, std::string(""));
					//cout << "\ncutted result row 1 left= " << endl;
					//cout << result << endl;
					result2 = std::accumulate(row1.begin()+crossover_position, row1.end(), std::string(""));
					//cout << "\ncutted result row 1 right= " << endl;
					//cout << result2 << endl;
					result3 = std::accumulate(row2.begin(), row2.end()-bits_chromosome+crossover_position, std::string(""));
					//cout << "\ncutted result row 2 left = " << endl;
					//cout << result3 << endl;
					result4 = std::accumulate(row2.begin()+crossover_position, row2.end(), std::string(""));
					//cout << "\ncutted result row 2 right= " << endl;
					//cout << result4 << endl;
					result5 = result+result4;
					//cout << "\ncombined result row 1 left +  row 2 right= " << endl;
					//cout << result5 << endl;
					result6 = result3+result2;
					//cout << "\ncombined result row 2 left +  row 1 right= " << endl;
					//cout << result6 << endl;
					newrow1 = {result5};
					newrow2 = {result6};

					crossover_population_final.push_back(newrow1);
					crossover_population_final.push_back(newrow2);
				}
				else if(i == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0)
				{
					row1 = stringmatrix_getRow(crossover_population, i);
					crossover_population_final.push_back(row1);
				}
				
				i = i+2;
			}
		}
		else if(n_crossoverpopulation < 2)
		{

		}
		
		cout << "\nCrossover population final : " << endl;

		for (const auto& row : crossover_population_final) 
		{
			for (const auto& s : row) 
			{
				cout << s << "";
			}
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
				for (const auto& s : row) 
				{
					cout << s << "";
				}
				cout << endl;
			}
			for (int i = 0; i < n_crossover; ++i)
			{
				vector<string> newrow_vec = stringmatrix_getRow(crossover_population_final,i);
				int c_index = index_vec[i];
				
				if(i < n_crossoverpopulation - 1 )
				{
					stringmatrix_deleteRow(selection_newpopulation,c_index);
					stringmatrix_addRow(selection_newpopulation, newrow_vec,c_index);
				}
				else if(i == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0)
				{

				}

			}
			
			cout << "\nPopulation after crossover: " << endl;
			for (const auto& row : selection_newpopulation) 
			{
				for (const auto& s : row) 
				{
					cout << s << "";
				}
				cout << endl;
			}
			// The rows that have undergone crossover become 1 string that represent 37 bits (the chromosome bits), while it should be 37 strings instead of 1.
			// This is the separation of 1 string of 101010100101010...101001 into 37 strings to represent the chromosome bits.
			for (int i = 0; i < n_crossover; ++i)
			{
				vector<string> newrow_vec = stringmatrix_getRow(selection_newpopulation,i);
				int c_index = index_vec[i];
				vector<string> separated_newrow = stringvector_groupeverynchars(newrow_vec,1);
				
				if(i < n_crossoverpopulation - 1 )
				{
					stringmatrix_deleteRow(selection_newpopulation,c_index);
					stringmatrix_addRow(selection_newpopulation, separated_newrow,c_index);
				}
				else if(i == n_crossoverpopulation - 1 && n_crossoverpopulation %2 != 0)
				{

				}

			}

			//cout << "\nPopulation after crossover and refashioned: " << endl;
			//printStringMatrix(selection_newpopulation);

		}
		else if(n_crossover < 2)
		{
			cout << "\nPopulation after crossover: " << endl;

			for (const auto& row : selection_newpopulation) 
			{
				for (const auto& s : row) 
				{
					cout << s << "";
				}
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

		std::random_device rd;
		std::mt19937 generate(rd()); // Mersenne Twister engine seeded with rd

		std::uniform_int_distribution<> distrib_column(1, bits_chromosome);
		std::uniform_int_distribution<> distrib_row(1, n_population);
		vector<double> mutated_row, mutated_column;

		for(int i=0; i < mutated; i++)
		{
			int n_row = distrib_row(generate);
			int n_col = distrib_column(generate);
			
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
			for (const auto& s : row) 
			{
				cout << s << "";
			}
			cout << endl;
		}
		for(int i=0; i < mutated; i++)
		{
			int mutated_row_index = mutated_row[i]-1;
			int mutated_column_index = mutated_column[i]-1;
			//cout << "mutated entry= " << selection_newpopulation[mutated_row_index][mutated_column_index] << endl;
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
			for (const auto& s : row) 
			{
				cout << s << "";
			}
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

		for(int i=0; i < n_population; i++)
		{
			//cout << "\ni = " << i << endl;

			vector<string> row1;
			string result1, result2;
			row1 = stringmatrix_getRow(selection_newpopulation,i);		

			result1 = std::accumulate(row1.begin(), row1.end()-size2, std::string(""));
			//cout << "\nsplit left= " << endl;
			//cout << result1 << endl;
			result2 = std::accumulate(row1.begin()+size1, row1.end(), std::string(""));
			//cout << "\nsplit right= " << endl;
			//cout << result2 << endl;
			int decimalValue = stoi(result1, nullptr, 2); // Convert string to int with base 2
			int decimalValue2 = stoi(result2, nullptr, 2); // Convert string to int with base 2

			//cout << "Decimal 1: " << decimalValue << endl;
			//cout << "Decimal 2: " << decimalValue2 << endl;
		
			binary_matrix1_new.push_back({result1});
			binary_matrix2_new.push_back({result2});

			double xjValue1 = aj1 + decimalValue*(divisiond(bj1-aj1,pow(2,size1)-1));
			vector_x1_new.push_back(xjValue1);	
			double xjValue2 = aj2 + decimalValue2*(divisiond(bj2-aj2,pow(2,size2)-1));
			vector_x2_new.push_back(xjValue2);			

		}
		//printStringMatrix(selection_newpopulation);
		cout << "\nNew binary x1: " << endl;
		printStringMatrix(binary_matrix1_new);
		cout << "\nNew binary x2: " << endl;
		printStringMatrix(binary_matrix2_new);

		for (int i = 0; i < n_population; ++i)
		{
			vector<string> newrow_vec1 = stringmatrix_getRow(binary_matrix1_new,i);
			vector<string> newrow_vec2 = stringmatrix_getRow(binary_matrix2_new,i);

			vector<string> separated_newrow1 = stringvector_groupeverynchars(newrow_vec1,1);
			vector<string> separated_newrow2 = stringvector_groupeverynchars(newrow_vec2,1);
		
			binary_matrix1_new_separated.push_back(separated_newrow1);
			binary_matrix2_new_separated.push_back(separated_newrow2);
	
		}

		//cout << "\nNew binary x1: " << endl;
		//printStringMatrix(binary_matrix1_new_separated);
		//cout << "\nNew binary x2: " << endl;
		//printStringMatrix(binary_matrix2_new_separated);

		cout << "\nNew x1 : " << endl;	
		printVector(vector_x1_new);
		cout << "\nNew x2: " << endl;
		printVector(vector_x2_new);

		vector<double> vector_F_new(n1,0.0);

		for (int i = 0; i <n1; ++i)
		{
			vector_F_new[i] = F[x1 == vector_x1_new[i], x2 == vector_x2_new[i]];
		}

		cout << "\nNew F(x_{1}, x_{2}) evaluation: " << endl;	
		printVector(vector_F_new);
	}
}

vector<double> GA_evaluation_vectoroutput(int n_population, double aj, double bj)
{
	// for precision of five places after the decimal point
	int size = divisiond(log((bj-aj)*pow(10,5)),log(2));
	
	vector<vector<string>> binary_matrix;
	vector<double> decimal_vector, xj_vector;	
	for (int i = 0; i< n_population; ++i)
	{
		// 1. Seed the random number engine
		random_device rd; 
		// Use a Mersenne Twister engine for high quality pseudo-random numbers
		mt19937 mersenne_engine {rd()}; 

		// 2. Define the distribution (uniform integer distribution between 0 and 1, inclusive)
		uniform_int_distribution<int> dist {0, 1};

		// 3. Create a generator function object (lambda function) that binds the distribution and engine
		auto generator = [&]() { return dist(mersenne_engine); };

		// 4. Create the vector and fill it using std::generate
		vector<int> binary_vector(size);
		vector<string> binary_vectorstring;
		generate(binary_vector.begin(), binary_vector.end(), generator);

		// 5. Print the generated vector (optional)
		/*
		cout << "Generated binary vector: ";
		for (int bit : binary_vector) 
		{
			cout << bit;
		}
		cout << endl;*/
		
		for (const int& s : binary_vector) 
		{
			if(s==0)
			{
				binary_vectorstring.push_back("0");
			}
			if(s==1)
			{
				binary_vectorstring.push_back("1");
			}
		}
		// Convert vector<int> into string
		// Using std::ostringstream, an output stream that writes to a string. It's efficient and flexible, especially if you need specific formatting or delimiters
		std::ostringstream oss;

		for (size_t i = 0; i < binary_vector.size(); ++i) 
		{
		oss << binary_vector[i];
		if (i != binary_vector.size() - 1) 
		{ // Add a space if it's not the last element
			oss << "";
		}
		}

		// Accumulate messenger_RNA as one string with no space
		string binaryNumber = oss.str();
		// string binaryNumber = binary_matrix[i][0]; // Example binary input
		int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

		//cout << "Binary: " << binaryNumber << endl;
		//cout << "Decimal: " << decimalValue << endl;

		decimal_vector.push_back(decimalValue);
		binary_matrix.push_back(binary_vectorstring);

		double xjValue = aj + decimalValue*(divisiond(bj-aj,pow(2,size)-1));
		xj_vector.push_back(xjValue);
		
	}

	return xj_vector;
}

#endif
#endif