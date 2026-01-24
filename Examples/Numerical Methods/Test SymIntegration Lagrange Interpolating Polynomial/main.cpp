// g++ -o result main.cpp -larmadillo -lsymintegration
// Merci beaucoup Freya..
// C++ program to find the solution of system of linear equations

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector>
#include "symintegrationc++.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

double division(double x, double y)
{
	return x/y;
}
Symbolic divisionsym(Symbolic x, Symbolic y)
{
	return x/y;
}

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic x("x");
	Symbolic f = x^(Symbolic(-1));
	dvec X = loadVectorFromFile("vectorX.txt");

	int N = X.size() ;
	dvec F(N,0.0);
	
	for (int i = 0; i<N; i++)
	{
		F[i] = f[x==X[i]];
	}
	
	cout <<"Vector x:" << endl;
	printVector(X);
	cout <<"\nf(x) = " << f <<endl;
	cout <<"\nVector f:" <<endl;
	printVector(F);

	Matrix<Symbolic> L(N,1);
	L[0][0] = 1;	
	L[1][0] = 1;
	L[2][0] = 1;
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			if (i != j)
			{
				L[i][0] *= (x-X[j])/(X[i]-X[j]);
			}			
		}
	}
	cout << "\nL :\n" << L << endl;

	Symbolic P ;
	for (int i = 0; i<N; i++)
	{
		P += F[i]*L[i][0];
	}
	cout <<"\nP(x) = " << P <<endl;
	cout <<"P(3) = " << P[x==3] <<endl;
	cout <<"f(3) = " << f[x==3] <<endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;


}