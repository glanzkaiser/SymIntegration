/*
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#include <vector>
#include <map>
#include <algorithm> // For std::max_element,  std::sort
#include <numeric> // For std::accumulate
#include <iostream>
#include <string>
#include <fstream> // For file operations

#include <random> // For random number generation
#include <chrono>

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_BIOLOGYDNA_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_BIOLOGYDNA_DEFINE

void printStringMatrix(vector<vector<string>> matrix) 
{
	 for (const auto& row : matrix) 
	{
		for (const auto& s : row) 
		{
			cout << s << " ";
		}
		cout << endl;
	}
}

vector<vector<string>> transposeStringMatrix(vector<vector<string>> &matrix) // Done beautifully
{
	if (matrix.empty() || matrix[0].empty() ) 
	{
        	// Handle empty matrices or invalid dimensions
		return {};
	}	
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();

	vector<vector<string>> transposed_matrix(cols, vector<string>(rows, ""));
	for (size_t i = 0; i < rows; ++i) 
	{
		for (size_t j = 0; j < cols; ++j) 
		{
			transposed_matrix[j][i] = matrix[i][j];	
		}
	}
	return transposed_matrix;
}

vector<vector<string>> DNA_full(vector<string> DNA)
{
	vector<string> DNA_pair;
	vector<string> messenger_RNA;

	vector<vector<string>> DNAfull_matrix;

	for (const std::string& s : DNA) 
	{
		if(s=="T")
		{
			DNA_pair.push_back("A");
		}
		if(s=="C")
		{
			DNA_pair.push_back("G");
		}
		if(s=="G")
		{
			DNA_pair.push_back("C");
		}
		if(s=="A")
		{
			DNA_pair.push_back("T");
		}
	}
	
	for (const std::string& s : DNA_pair) 
	{
		if(s=="T")
		{
			messenger_RNA.push_back("A");
		}
		if(s=="C")
		{
			messenger_RNA.push_back("G");
		}
		if(s=="G")
		{
			messenger_RNA.push_back("C");
		}
		if(s=="A")
		{
			messenger_RNA.push_back("U");
		}
	}

	DNAfull_matrix.push_back(DNA);
	DNAfull_matrix.push_back(DNA_pair);
	DNAfull_matrix.push_back(messenger_RNA);
		
	return DNAfull_matrix;

}

vector<vector<string>> polypeptide(vector<string> DNA)
{
	vector<string> DNA_pair;
	vector<string> messenger_RNA;

	vector<vector<string>> polypeptide_matrix;

	for (const std::string& s : DNA) 
	{
		if(s=="T")
		{
			DNA_pair.push_back("A");
		}
		if(s=="C")
		{
			DNA_pair.push_back("G");
		}
		if(s=="G")
		{
			DNA_pair.push_back("C");
		}
		if(s=="A")
		{
			DNA_pair.push_back("T");
		}
	}
	
	for (const std::string& s : DNA_pair) 
	{
		if(s=="T")
		{
			messenger_RNA.push_back("A");
		}
		if(s=="C")
		{
			messenger_RNA.push_back("G");
		}
		if(s=="G")
		{
			messenger_RNA.push_back("C");
		}
		if(s=="A")
		{
			messenger_RNA.push_back("U");
		}
	}

	// Accumulate messenger_RNA as one string with no space
	string polypeptide_accumulate = accumulate(messenger_RNA.begin(), messenger_RNA.end(), std::string(""));
	
	vector<string> first_three_letters;
	int N = polypeptide_accumulate.length()/3;
	
	int j = 0;
	for (int i =0; i < N; i++)
	{
		string amino = polypeptide_accumulate.substr(j,3);
		first_three_letters.push_back(amino);
		j = j+3;
		
	}
	polypeptide_matrix.push_back(first_three_letters);
		
	// translation of mRNA code to polypeptide / amino acids
	size_t rows = polypeptide_matrix.size();
	size_t cols = polypeptide_matrix[0].size();

	for (size_t i = 0; i < rows; ++i) 
	{
		for (size_t j = 0; j < cols; ++j) 
		{
			if(polypeptide_matrix[i][j] == "UUU" || polypeptide_matrix[i][j] == "UUC")
			{
				polypeptide_matrix[i][j] = "Phe";
			}
			if(polypeptide_matrix[i][j] == "UUA" || polypeptide_matrix[i][j] == "UUG")
			{
				polypeptide_matrix[i][j] = "Leu";
			}
			if(polypeptide_matrix[i][j] == "CUU" || polypeptide_matrix[i][j] == "CUC" || polypeptide_matrix[i][j] == "CUA" || polypeptide_matrix[i][j] == "CUG")
			{
				polypeptide_matrix[i][j] = "Leu";
			}
			if(polypeptide_matrix[i][j] == "AUU" || polypeptide_matrix[i][j] == "AUC" || polypeptide_matrix[i][j] == "AUA" )
			{
				polypeptide_matrix[i][j] = "Ile";
			}	
			if(polypeptide_matrix[i][j] == "AUG" )
			{
				polypeptide_matrix[i][j] = "Met";
			}	
			if(polypeptide_matrix[i][j] == "GUU" || polypeptide_matrix[i][j] == "GUC" || polypeptide_matrix[i][j] == "GUA" || polypeptide_matrix[i][j] == "GUG")
			{
				polypeptide_matrix[i][j] = "Val";
			}
			if(polypeptide_matrix[i][j] == "UCU" || polypeptide_matrix[i][j] == "UCC" || polypeptide_matrix[i][j] == "UCA" || polypeptide_matrix[i][j] == "UCG")
			{
				polypeptide_matrix[i][j] = "Ser";
			}			
			if(polypeptide_matrix[i][j] == "CCU" || polypeptide_matrix[i][j] == "CCC" || polypeptide_matrix[i][j] == "CCA" || polypeptide_matrix[i][j] == "CCG")
			{
				polypeptide_matrix[i][j] = "Pro";
			}
			if(polypeptide_matrix[i][j] == "ACU" || polypeptide_matrix[i][j] == "ACC" || polypeptide_matrix[i][j] == "ACA" || polypeptide_matrix[i][j] == "ACG")
			{
				polypeptide_matrix[i][j] = "Thr";
			}
			if(polypeptide_matrix[i][j] == "GCU" || polypeptide_matrix[i][j] == "GCC" || polypeptide_matrix[i][j] == "GCA" || polypeptide_matrix[i][j] == "GCG")
			{
				polypeptide_matrix[i][j] = "Ala";
			}	
			if(polypeptide_matrix[i][j] == "UAU" || polypeptide_matrix[i][j] == "UAC" )
			{
				polypeptide_matrix[i][j] = "Tyr";
			}
			if(polypeptide_matrix[i][j] == "UAA" || polypeptide_matrix[i][j] == "UAG")
			{
				polypeptide_matrix[i][j] = "Stop";
			}
			if(polypeptide_matrix[i][j] == "CAU" || polypeptide_matrix[i][j] == "CAC")
			{
				polypeptide_matrix[i][j] = "His";
			}
			if(polypeptide_matrix[i][j] == "CAA" || polypeptide_matrix[i][j] == "CAG")
			{
				polypeptide_matrix[i][j] = "Gln";
			}
			if(polypeptide_matrix[i][j] == "AAU" || polypeptide_matrix[i][j] == "AAC")
			{
				polypeptide_matrix[i][j] = "Asn";
			}
			if(polypeptide_matrix[i][j] == "AAA" || polypeptide_matrix[i][j] == "AAG")
			{
				polypeptide_matrix[i][j] = "Lys";
			}
			if(polypeptide_matrix[i][j] == "GAU" || polypeptide_matrix[i][j] == "GAC")
			{
				polypeptide_matrix[i][j] = "Asp";
			}
			if(polypeptide_matrix[i][j] == "GAA" || polypeptide_matrix[i][j] == "GAG")
			{
				polypeptide_matrix[i][j] = "Glu";
			}
			if(polypeptide_matrix[i][j] == "UGU" || polypeptide_matrix[i][j] == "UGC")
			{
				polypeptide_matrix[i][j] = "Cys";
			}
			if(polypeptide_matrix[i][j] == "UGA" )
			{
				polypeptide_matrix[i][j] = "Stop";
			}
			if(polypeptide_matrix[i][j] == "UGG")
			{
				polypeptide_matrix[i][j] = "Trp";
			}
			if(polypeptide_matrix[i][j] == "CGU" || polypeptide_matrix[i][j] == "CGC" || polypeptide_matrix[i][j] == "CGA" || polypeptide_matrix[i][j] == "CGG")
			{
				polypeptide_matrix[i][j] = "Arg";
			}
			if(polypeptide_matrix[i][j] == "AGU" || polypeptide_matrix[i][j] == "AGC")
			{
				polypeptide_matrix[i][j] = "Ser";
			}
			if(polypeptide_matrix[i][j] == "AGA" || polypeptide_matrix[i][j] == "AGG")
			{
				polypeptide_matrix[i][j] = "Arg";
			}
			if(polypeptide_matrix[i][j] == "GGU" || polypeptide_matrix[i][j] == "GGC" || polypeptide_matrix[i][j] == "GGA" || polypeptide_matrix[i][j] == "GGG")
			{
				polypeptide_matrix[i][j] = "Gly";
			}
		}
	}

	// Insert the new row at the beginning (before the first element)
	polypeptide_matrix.insert(polypeptide_matrix.begin(), first_three_letters);

	return polypeptide_matrix;
}

void mRNA(vector<string> DNA)
{	
	vector<vector<string>> DNA_full;
	vector<string> DNA_pair;
	vector<string> messenger_RNA;
	for (const std::string& s : DNA) 
	{
		if(s=="T")
		{
			DNA_pair.push_back("A");
		}
		if(s=="C")
		{
			DNA_pair.push_back("G");
		}
		if(s=="G")
		{
			DNA_pair.push_back("C");
		}
		if(s=="A")
		{
			DNA_pair.push_back("T");
		}
	}
	
	for (const std::string& s : DNA_pair) 
	{
		if(s=="T")
		{
			messenger_RNA.push_back("A");
		}
		if(s=="C")
		{
			messenger_RNA.push_back("G");
		}
		if(s=="G")
		{
			messenger_RNA.push_back("C");
		}
		if(s=="A")
		{
			messenger_RNA.push_back("U");
		}
	}
	

	cout << "\nDNA pair:" << endl;
	/*for (const std::string& s : DNA_pair) 
	{
		cout << s << endl;
	}*/
	 // Iterating and printing
	DNA_full.push_back(DNA);
	DNA_full.push_back(DNA_pair);

	for (const auto& row : DNA_full) 
	{
		for (const auto& s : row) 
		{
			cout << s << " ";
		}
		cout << endl;
	}
	//cout << DNA_full[0].size() << endl;
	cout << "\nmRNA:" << endl;
	for (const std::string& s : messenger_RNA) 
	{
		cout << s << endl;
	}
}



#endif
#endif