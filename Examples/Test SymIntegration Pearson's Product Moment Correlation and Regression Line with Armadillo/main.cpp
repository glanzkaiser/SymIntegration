// g++ main.cpp -o result -larmadillo

#include <iostream>
#include <armadillo>
#include "symintegrationc++.h"
#include <chrono>

using namespace std::chrono;

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	fmat A;
	
	A.load("MatrixA.txt");	
	cout << "Matrix A:" << endl;
	A.print();
	int N = A.n_rows ;

	// Copy matrix A from Armadillo into W Symbolic matrix for SymIntegration
	Matrix<Symbolic> W(N,2);
	for(int i=0; i < N; ++i)
	{
		W[i][0]= A.col(0)[i];
		W[i][1]= A.col(1)[i];
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


