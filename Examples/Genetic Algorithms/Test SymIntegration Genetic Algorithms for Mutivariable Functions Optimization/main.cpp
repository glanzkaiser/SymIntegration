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
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

double division(double x, double y)
{
	return x/y;
}

int main(void)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic x1("x1"), x2("x2");
	Symbolic F  = 21.5 +x1*sin(4*π*x1) + x2*sin(20**x2);
	/*svec DNA = loadStringVector("initialpopulation.txt");
	
	printStringVector(DNA);
	int n = DNA.size();
	for (int i = 0; i < n; ++i)
	{
		string binaryNumber = DNA[i]; // Example binary input
		int decimalValue = stoi(binaryNumber, nullptr, 2); // Convert string to int with base 2

		cout << "Binary: " << binaryNumber << endl;
		cout << "Decimal: " << decimalValue << endl;
	}*/

	//GA_generateinitialpopulation(18,  10);
	//GA_evaluation(10, -3,12.1);
	//GA_evaluation(10, 4.1,5.8);

	//vector<double> vector_x1 = GA_evaluation_vectoroutput(10, -3,12.1);
	//vector<double> vector_x2 = GA_evaluation_vectoroutput(10, 4.1,5.8);
	//cout << "\nx_{1}:"<<endl;
	//printVector(vector_x1);
	//cout << "\nx_{2}:"<<endl;
	//printVector(vector_x2);

	GA_maximizationproblem_functionwithtwovariables_N_generations(F,x1, x2, 10, -3,12.1, 4.1, 5.8,2);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	

	return 0; 
}
