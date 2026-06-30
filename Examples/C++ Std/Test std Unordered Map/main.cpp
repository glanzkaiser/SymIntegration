// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include "symintegrationc++.h"
#include <bits/stdc++.h>
#include <cmath>
#include <chrono>

#define π 3.1415926535897f

using namespace std::chrono;
using namespace std;
using namespace SymbolicConstant;

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

//The std::unordered_map<std::string, int> is a C++ Standard Library container that stores key-value pairs (where the key is a string and the value is an integer) using a hash table for fast access

// To use it, you must include the <unordered_map> and <string> headers.

// Driver code
int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();
	std::unordered_map<string, int> ages = { {"Alice", 25}, {"Bob", 30} };

	// Insertion
	ages["Charlie"] = 35;                      // Using operator[]
	ages.insert({"Diana", 28});               // Using insert()

	// Access
	cout << "Alice's age: " << ages["Alice"] << endl;

	// Check if key exists (C++20)
	//if (ages.contains("Bob")) 
	//{ 
	//	cout << "Bob is in the map." << endl;
	//}

	// Iteration
	for (const auto& [name, age] : ages) 
	{    // Structured bindings (C++17)
		cout << name << ": " << age << endl;
	}

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}

/*


#include <iostream>
#include <unordered_set> // Required header

int main() {
    // 1. Initialization
    std::unordered_set<int> mySet = {10, 20, 30};

    // 2. Insertion
    mySet.insert(40);
    mySet.insert(20); // Duplicate: ignored silently

    // 3. Search Lookup
    if (mySet.find(20) != mySet.end()) {
        std::cout << "20 is in the set\n";
    }

    // C++20 approach for lookup
    if (mySet.contains(40)) {
        std::cout << "40 is in the set\n";
    }

    // 4. Deletion
    mySet.erase(10);

    // 5. Iterating through the set
    for (int num : mySet) {
        std::cout << num << " "; // Order is unpredictable
    }
    
    return 0;
}
*/