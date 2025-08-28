// g++ -o result main.cpp -lsymintegration  
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector> // For std::vector (example container)
#include "symintegrationc++.h"
#include <algorithm> // For std::sort

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f

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
	int n, t1 = 0, t2 = 1, nextTerm = 0;

    	cout << "Enter the number of terms: ";
    	cin >> n;

    	cout << "Fibonacci Series: ";

    	for (int i = 1; i <= n; ++i) {
        // Prints the first two terms.
        if(i == 1) {
            cout << t1 << ", ";
            continue;
     	   }
        if(i == 2) {
            cout << t2 << ", ";
            continue;
        }
        nextTerm = t1 + t2;
        t1 = t2;
        t2 = nextTerm;
        
        cout << nextTerm << ", " ;
    	}
	cout << "\n" << endl;

	cout << fibonacciseries(n) << endl ;
	return 0;
}