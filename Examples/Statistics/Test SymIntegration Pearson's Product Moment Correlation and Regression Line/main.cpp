// g++ main.cpp -o result -lsymintegration

#include <iostream>
#include "symintegrationc++.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	dmat doubleMatrix = loadMatrixFromFile("MatrixA.txt");
	cout << "\nMatrix A : " << endl;
	printMatrix(doubleMatrix);
	int N = doubleMatrix.size() ;

	// Copy doubleMatrix into W Symbolic matrix for SymIntegration
	Matrix<Symbolic> W(N,2);
	for(int i=0; i < N; ++i)
	{
		W[i][0]= doubleMatrix[i][0];
		W[i][1]= doubleMatrix[i][1];
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


