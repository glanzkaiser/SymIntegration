#include <cmath>
#include <hamzstplot/hamzstplot.h>
#include <random>
#include "symintegrationc++.h"

using namespace std::chrono;
using namespace std;

double division(double x, double y)
{
	return x/y;
}

int main() {
	using namespace hamzstplot;

	// Get starting timepoint
	//auto start = high_resolution_clock::now();
	dvec alloyA =loadVectorFromFile("alloyA.txt");
 	dvec alloyB =loadVectorFromFile("alloyB.txt");
	dvec quantileA = normalquantile(alloyA, 0.3);
	dvec quantileB = normalquantile(alloyB, 0.3);

	// sort the vector of alloy A and alloy B
	sort(alloyA.begin(), alloyA.end()); 
	sort(alloyB.begin(), alloyB.end()); 

	//printVector(quantileA);
	auto x =quantileB;
	auto y = alloyB;
	double sz = 6; // size for the scatter plot

	auto l = scatter(x, y, sz);
	l->marker_color({0.f, .5f, .5f});
	l->marker_face_color({0.f, .7f, .7f});

	xlabel("Normal Quantile");
	ylabel("Quantile");
	//yrange({73, 93});
	title("Normal Quantile-Quantile Plot for Alloy B");

	show();

	// Get ending timepoint
	//auto stop = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(stop - start);

	//cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}