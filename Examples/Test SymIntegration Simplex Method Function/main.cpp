// g++ -o result main.cpp -lsymbolicc++
// Merci beaucoup Freya et Sentinel

#include<bits/stdc++.h>
#include<iostream>
#include "symintegrationc++.h"
#include<vector>
#include <chrono>
#include <algorithm> // For std::next_permutation
#include <string>
using namespace std::chrono;
using namespace std;

#define R 2 // number of rows
#define C 7 // number of columns

void comb(int N, int K)
{
	string bitmask(K, 1); // K leading 1's
	bitmask.resize(N, 0); // N-K trailing 0's
 
	// print integers and permute bitmask
	do 
	{
	for (int i = 0; i < N; ++i) // [0..N-1] integers
	{
		if (bitmask[i]) 
		{
			cout << " " << i;
		}
	}
	cout << endl;
	} 
	while (prev_permutation(bitmask.begin(), bitmask.end()));
}

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

 	// Construct a symbolic matrix of size 2 X 7
	Matrix<Symbolic> B_mat(2,7);
	B_mat[0][0] = 2500;	B_mat[0][1] = 3000;	B_mat[0][2] = 600;  B_mat[0][3] = 100;	B_mat[0][4] = 600;  B_mat[0][5] = 0;	B_mat[0][6] = 3000;
	B_mat[1][0] = 80;		B_mat[1][1] = 150;	B_mat[1][2] = 20;	B_mat[1][3] = 10;	B_mat[1][4] = 40;	B_mat[1][5] = 0;	B_mat[1][6] = 100;
	//cout << "The matrix A:\n" << B_mat <<endl;
	
	vector<double> obj_function = {3,10,1,2,3}; // 3x1 + 10x2+ x3 + 2x4 + 3x5

	cout << "Simplex Method :\n" << simplexmethod(B_mat, obj_function, C, R, 5) <<endl;
	//comb(5, 2) ;	

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}