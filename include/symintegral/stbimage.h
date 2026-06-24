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

#ifndef SYMINTEGRATION_CPLUSPLUS_STBIMAGE

#ifdef  SYMBOLIC_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_STBIMAGE
#ifndef SYMINTEGRATION_CPLUSPLUS_STBIMAGE_DECLARE
#define SYMINTEGRATION_CPLUSPLUS_STBIMAGE_DECLARE

void ResizeAndGrayscale(const std::string &, const std::string &, int, int); 

uint32_t swap_endian(uint32_t);
void convert_to_idx(const vector<string> & , const vector<int> & , const std::string &, const std::string & , int , int ); 

int reverseInt(int);
vector<vector<vector<int>>> loadIDX3(const std::string&);
vector<int> loadIDX1(const std::string&);
#endif
#endif


#endif
