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

	dmat X0 = loadMatrixFromFile("MatrixX.txt");
	int n = X0.size();
	// Create a column vector of one
	dvec ones = createVector(n,1);
	// add the vector of ones to the first column 
	dmat X = addColumn(X0,ones,0);

	cout << "\nMatrix X : " << endl;
	printMatrix(X);
	dmat y = loadMatrixFromFile("Vectory.txt");
	cout << "\nVector y : " << endl;
	printMatrix(y);
	
	cout << "Multiple regression with matrix algebra:" << endl;
	dmat multiplereg = multipleregressionfull(X,y,0.05);
	//cout << "\n(X^{T} X)^{-1} X^{T} y: " << endl;
	//printMatrix(multiplereg);
	/*int N = X.size() ;
	int k = X[0].size();
	cout << "\ny" << "\t\t\t\t" << "y estimated" << endl;
	for(int i = 0; i<N; ++i)
	{
		double y_estimated = multiplereg[0][0]*X[i][0] + multiplereg[1][0]*X[i][1] + multiplereg[2][0]*X[i][2] + multiplereg[3][0]*X[i][3] ;
		cout << y[i][0] << "\t\t\t" << y_estimated << endl;
	} */
	//cout << "\nFor 50% humidity, 76^{0} F and 29.30 barometric pressure, the estimated multiple regression is:" << endl;
	//double y_estimated2 = multiplereg[0][0]*1 + multiplereg[1][0]*50 + multiplereg[2][0]*76+ multiplereg[3][0]*29.30;
	//cout << y_estimated2<<endl;
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	
	return 0;
}


