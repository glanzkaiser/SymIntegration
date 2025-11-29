// g++ main.cpp -o result -larmadillo

#include <iostream>
#include <armadillo>
#include <chrono>

using namespace std::chrono;

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	fmat A;
	
	A.load("MatrixA.txt");	
	cout << "Matrix A:" << endl;
	A.print();
	int N = A.n_rows ;
	int C = A.n_cols;
		
	double sum_xy, sum_x, sum_y, x_bar, y_bar, sum_xsquared, sum_ysquared, Sxy, Sx, Sy, r_pearson;
	cout << "Column vector from A: " << endl;


	for(int i=0; i < C; ++i)
	{
		cout << "column " << i << " :" <<endl;
		cout << A.col(i) << endl;
		
	}
	for(int i=0; i < N; ++i)
	{
		//cout << A.col(0)[i]*A.col(1)[i] << endl;
		sum_xy += A.col(0)[i]*A.col(1)[i] ;
	}
	x_bar = sum(A.col(0))/N;
	y_bar = sum(A.col(1))/N;
	sum_xsquared = dot(A.col(0),A.col(0))	;
	sum_ysquared = dot(A.col(1),A.col(1))	;

	Sxy = (sum_xy/N) - x_bar*y_bar;
	Sx = sqrt(sum_xsquared/N - (x_bar*x_bar));
	Sy = sqrt(sum_ysquared/N - (y_bar*y_bar));

	cout << "sum xy = " << sum_xy << endl;
	cout << "sum x^2 = " << sum_xsquared << endl;
	cout << "sum y^2 = " << sum_ysquared << endl;
	cout << "x bar = " << x_bar << endl;
	cout << "y bar = " << y_bar << endl;

	cout << "\nS_{xy} = " << Sxy << endl;
	cout << "S_{x} = " << Sx << endl;
	cout << "S_{y} = " << Sy << endl;
	
	r_pearson = Sxy/(Sx*Sy);
	cout << "r = " <<  r_pearson << endl;

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;

	return 0;
}


