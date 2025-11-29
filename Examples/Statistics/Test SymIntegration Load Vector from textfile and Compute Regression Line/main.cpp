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

#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#define π 3.1415926535897f
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

int main() {
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	vector<int> vectorx;
	vector<int> vectory;
	ifstream inputFile("vectorx.txt"); // your file name
	ifstream inputFile2("vectory.txt"); // your file name

	if (inputFile.is_open()) 
	{
		int num;
		while (inputFile >> num) 
		{ // Reads numbers separated by whitespace
			vectorx.push_back(num);
		}
		inputFile.close();
		//cout << "File loaded successfully into vector<int>." << endl;
		} 
		else 
		{
			cerr << "Error: Unable to open file." << endl;
			return 1;
	}
	if (inputFile2.is_open()) 
	{
		int num2;
		while (inputFile2 >> num2) 
		{ // Reads numbers separated by whitespace
			vectory.push_back(num2);
		}
		inputFile2.close();
		//cout << "File loaded successfully into vector<int>." << endl;
		} 
		else 
		{
			cerr << "Error: Unable to open file." << endl;
			return 1;
	}

	int N = vectorx.size() ;

	// Copy matrix A from Armadillo into W Symbolic matrix for SymIntegration
	Matrix<Symbolic> W(N,2);
	for(int i=0; i < N; ++i)
	{
		W[i][0]= vectorx[i];
		W[i][1]= vectory[i];
	}
	
	cout << "W:\n" << W <<endl;
	
	cout << rpearson(W,N) << endl;
	cout << "\nRegression line, y = " << regressionline(W,N) << endl;


	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}
