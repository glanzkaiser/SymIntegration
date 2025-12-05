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

using namespace std;
using namespace SymbolicConstant;

int main() {
	dmat M = loadMatrixFromFile("Matrix.txt");
	dvec x = loadVectorFromFile("Vectorx.txt");
	dvec y = loadVectorFromFile("Vectory.txt");

	printMatrix(M);

	cout << "\nAdd column" << endl;
	dmat M1 = addColumn(M,x,0);
	printMatrix(M1);
	cout << "\nAdd row" << endl;
	dmat M2 = addRow(M,y,1);
	printMatrix(M2);
	
	cout << "\nCreate matrix" << endl;
	dmat M3 = createMatrix(3,5,82);
	printMatrix(M3);
	cout << "\nCreate vector" << endl;
	dvec V1 = createVector(10,108.0);
	printVector(V1);

	return 0;
}
