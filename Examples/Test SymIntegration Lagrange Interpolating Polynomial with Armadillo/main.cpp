// g++ -o result main.cpp -larmadillo -lsymintegration
// Merci beaucoup Freya..
// C++ program to find the solution of system of linear equations

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector>
#include <armadillo>
#include "symintegrationc++.h"
#include <chrono>

using namespace std::chrono;
using namespace std;
using namespace arma;

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
	mat X;
	X.load("vectorX.txt");
	mat F;

	int N = X.n_rows ;
	F.set_size(N,1);
	
	for (int i = 0; i<N; i++)
	{
		F(i) = f[x==X(i)];
	}
	
	cout <<"Vector x:" << "\n" << X <<endl;
	cout <<"f(x) = " << f <<endl;
	cout <<"\nVector f:" << "\n" << F <<endl;
	
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
				L[i][0] *= (x-X(j))/(X(i)-X(j));
			}			
		}
	}
	cout << "L :\n" << L << endl;

	Symbolic P ;
	for (int i = 0; i<N; i++)
	{
		P += F(i)*L[i][0];
	}
	cout <<"\nP(x) = " << P <<endl;
	cout <<"P(3) = " << P[x==3] <<endl;
	cout <<"f(3) = " << f[x==3] <<endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	/*Symbolic L0, L1, L2;
	L0 = 1, L1= 1, L2 = 1;
	for (int i = 0; i<N; i++)
	{
		if( i!=0)
		{
			L0 *= (x-X(i))/(X(0)-X(i));
		}
		if( i!=1)
		{
			L1 *= (x-X(i))/(X(1)-X(i));
		}
		if( i!=2)
		{
			L2 *= (x-X(i))/(X(2)-X(i));
		}
	}
	cout <<"L_{0}(x) = " << L0 <<endl;
	cout <<"L_{1}(x) = " << L1 <<endl;
	cout <<"L_{2}(x) = " << L2 <<endl;*/

}