// g++ -o result main.cpp -lsymintegration
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

// Driver program
int main()
{	
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	dvec vecu = loadVectorFromFile("vectorx.txt");
	dvec veca = loadVectorFromFile("vectory.txt");
		
	cout << "\nVector u : " << endl;
	printVector(vecu);
	
	cout << "\nVector a : " << endl;
	printVector(veca);
	
	cout << "\nnorm(a)^2 : " << norm(veca)*norm(veca) << endl;
	
	cout << "\ndot(u,a) : " << dot(vecu,veca) << endl;
	
	dvec w1 = orthogonalprojection(vecu,veca);
	
	cout << "\nOrthogonal projection of u on a: " << endl;
	printVector(w1);
	
	cout << "\nVector component of u orthogonal to a: " << endl;
	printVector(subtract(vecu,w1));
	
	cout << "\nu + a : " << endl;
	printVector(add(vecu,veca));
	
	cout << "\nu - a : " << endl;
	printVector(subtract(vecu,veca));
	
	cout << "\n5*u  : " << endl;
	printVector(scalarmultiplication(vecu,5));
	
	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}