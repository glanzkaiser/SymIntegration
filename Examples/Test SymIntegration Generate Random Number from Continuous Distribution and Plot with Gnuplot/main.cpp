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

// g++ -o main main.cpp -lsymintegration -lboost_iostreams
// make
// ./main > gamma.dat

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>

#include "gnuplot-iostream.h"

#define π 3.1415926535897f
#include <chrono>

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
	
	cout << randomnumbergamma(2,1.5,50) << endl;
	
	Gnuplot gp;

	// Don't forget to put "\n" at the end of each line!
	gp << "unset key\n";
	gp << "bin_width = 1.0\n";
	gp << "bin_number(x) = floor(x/bin_width)\n";
	gp << "rounded(x) = bin_width * ( bin_number(x) + 0.5 )\n";

	gp << "set xrange [0:8]\n";
	gp << "set yrange [0:14]\n";
	gp << "set xtics 1\n";
	gp << "set title 'Gamma Distribution' \n";
	gp << "set key off\n"; // supress the key / no legend
	gp << "plot 'gamma.dat' using (rounded($1)):(1) smooth frequency with histeps title 'Data Frequency' lc rgb 'blue'\n";

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	//cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	
	return 0; 
}
