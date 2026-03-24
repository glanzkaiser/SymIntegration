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

#include <random>
#include <iostream>
#include <vector>

#ifndef SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS

#ifdef  SYMBOLIC_FORWARD
#ifndef SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS_FORWARD
#define SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS_FORWARD

#endif
#endif

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS
#ifndef SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_VECTORMATRIXOPERATIONS_DECLARE

using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<std::vector<Complex>>;

vector<string> loadStringVector(const string&);
vector<vector<string>> loadStringMatrix(const string&);
void printStringVector(vector<string>);
void printStringMatrix(vector<vector<string>>); 
//vector<vector<string>> transposeStringMatrix(vector<vector<string>>&); 

void stringmatrix_addColumn(vector<vector<string>>&, const vector<vector<string>>&);
vector<string> stringmatrix_getColumn(vector<vector<string>>&, int);
void stringmatrix_deleteColumn(vector<vector<string>>&, int);

void stringmatrix_addRow(vector<vector<string>>&, vector<string>&, int);
vector<string> stringmatrix_getRow(vector<vector<string>>&, int);
void stringmatrix_deleteRow(vector<vector<string>>&, int);

vector<vector<string>> transpose_stringmatrix(const vector<vector<string>>& ); 
void stringmatrix_removeduplicaterows_sort(vector<vector<string>>& );
vector<string> stringvector_groupeverynchars(const vector<string>&, int); 

#endif
#endif


#endif
