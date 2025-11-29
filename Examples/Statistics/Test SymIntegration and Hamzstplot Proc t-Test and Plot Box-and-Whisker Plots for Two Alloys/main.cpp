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
	dvec alloy =loadVectorFromFile("alloy.txt");
 	dvec alloyA =loadVectorFromFile("alloyA.txt");
 	dvec alloyB =loadVectorFromFile("alloyB.txt");
 	svec alloyname = loadStringVector("alloyname.txt");
	
	hypothesistest_unequalunknownvariances_righttailed(alloyA, alloyB, 0, 0.05,1);
	hypothesistest_twosamplesvariances_righttailed(alloyA,alloyB, 0.05);
	hypothesistest_twosamplesvariances_twotailed(alloyA,alloyB, 0.05);

	//hold(on);
	boxplot(alloy, alloyname);
	//hold(off);
	xlabel("Alloys");
	ylabel("Deflection");
	yrange({73, 93});
	title("Breaking Strength for Alloys");

	show();

	// Get ending timepoint
	//auto stop = high_resolution_clock::now();
	//auto duration = duration_cast<microseconds>(stop - start);

	//cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}