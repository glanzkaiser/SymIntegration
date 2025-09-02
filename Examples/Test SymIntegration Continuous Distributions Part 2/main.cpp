// g++ -o result main.cpp -lsymintegration 
// Merci beaucoup Freya..

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
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

	cout << "\ncauchypdf(2) = " << cauchypdf(2) << endl;
	
	cout << "\nchisquaredpdf(2.167,7) = " << chisquaredpdf(2.167,7) << endl;
	cout << "\nchisquaredcdf(2.167,7) = " << chisquaredcdf(2.167,7) << endl;
	cout << "\nchisquaredmgf(2.167,7) = " << chisquaredmgf(2.167,7) << endl;
	cout << "\nchisquaredmean(2.167,7) = " << chisquaredmean(2.167,7) << endl;
	cout << "\nchisquaredvar(2.167,7) = " << chisquaredvar(2.167,7) << endl;
	
	cout << "\nFpdf(3.22, 6, 10) = " << Fpdf(3.22,6,10) << endl;
	cout << "\nFcdf(3.22, 6, 10) = " << Fcdf(3.22,6,10) << endl;
	cout << "\nFmean(3.22, 6, 10) = " << Fmean(3.22,6,10) << endl;
	cout << "\nFvar(3.22, 6, 10) = " << Fvar(3.22,6,10) << endl;

	cout << "\ntpdf(3, 6) = " << tpdf(3,6) << endl;
	cout << "\ntcdf(3, 6) = " << tcdf(3,6) << endl;
	cout << "\ntmean(3, 6) = " << tmean(3,6) << endl;
	cout << "\ntvar(3, 6) = " << tvar(3,6) << endl;
	
	cout << "\nLaplacepdf(2,3.3) = " << laplacepdf(2,3.3) << endl;
	cout << "\nLaplacecdf(2,3.3) = " << laplacecdf(2,3.3) << endl;
	cout << "\nLaplacemgf(2,3.3) = " << laplacemgf(2,3.3) << endl;
	cout << "\nLaplacemean(2,3.3) = " << laplacemean(2,3.3) << endl;
	cout << "\nLaplacevar(2,3.3) = " << laplacevar(2,3.3) << endl;

	cout << "\nLogisticpdf(2,3.3) = " << logisticpdf(2,3.3) << endl;
	cout << "\nLogisticcdf(2,3.3) = " << logisticcdf(2,3.3) << endl;
	cout << "\nLogisticmean(2,3.3) = " << logisticmean(2,3.3) << endl;
	cout << "\nLogisticvar(2,3.3) = " << logisticvar(2,3.3) << endl;
	

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}