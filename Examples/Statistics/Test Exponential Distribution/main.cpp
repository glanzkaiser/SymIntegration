// g++ -o result main.cpp -lsymintegration 
// Merci beaucoup Freya..

#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
 
#include <chrono>

using namespace std::chrono;
using namespace std;

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	random_device rd;
	mt19937 gen(rd());
 
	// if particles decay once per second on average,
	// how much time, in seconds, until the next one?
	exponential_distribution<> d(1);
 
	map<int, int> hist;
	for (int n = 0; n != 10000; ++n)
	{
		++hist[2 * d(gen)];
 	}
	for (auto const& [x, y] : hist)
	{		
		if (y / 100.0 > 0.5)
		{            
		cout << fixed << setprecision(1) << x / 2.0 << '-' << (x + 1) / 2.0 << ' ' << string(y / 100, '*') << '\n';
		}
	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

}