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

#define π 3.1415926535897f

using namespace std;
using namespace SymbolicConstant;

#include <fstream>
#include <iostream>
#include <string>
#include <vector>


int main() {
	std::vector<std::string> lines;
	std::ifstream inputFile("data.txt"); // Replace "data.txt" with your file name

	if (inputFile.is_open()) 
	{
        std::string line;
	while (std::getline(inputFile, line)) 
	{
		lines.push_back(line);
        }
	inputFile.close();
	std::cout << "File loaded successfully into vector<string>." << std::endl;
	} 
	else 
	{
	std::cerr << "Error: Unable to open file." << std::endl;
	return 1;
    }

	// Optional: Print the vector contents to verify
	for (const std::string& s : lines) 
	{
		std::cout << s << std::endl;
	}

    return 0;
}
