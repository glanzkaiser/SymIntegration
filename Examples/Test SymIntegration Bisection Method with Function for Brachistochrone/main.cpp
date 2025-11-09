// Merci beaucoup Freya et Sentinel
// g++ main.cpp -o result -lsymintegration
#include <iostream>
#include "symintegrationc++.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic θ0("θ0"), g, f,h;
	
	g = 1 - cos(θ0);
	h = 2*θ0 - 2*sin(θ0);
	f = g-h;
	cout << "f = " << f << endl;

	float a = 1;
	float b = 2;
	int N = 17;

	bisectionmethod(f,θ0,a,b,N);

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}