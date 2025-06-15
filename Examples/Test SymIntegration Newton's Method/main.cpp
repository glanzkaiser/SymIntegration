// Merci beaucoup Freya et Sentinel
// g++ main.cpp -o result -lsymintegration
#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include "symintegrationc++.h"

#define DEGTORAD 0.0174532925199432957f
#define RADTODEG 57.295779513082320876f
#define pi  3.1415926535897

using namespace std;


int main()
{
	Symbolic x("x"), f, fd, fp, fpd, pn, p0, radtodeg;
	p0 = (pi/4);
	int N = 5;
	radtodeg = 57.295779513082320876;

	f = cos(x) - x; // the input x has to be in Radian
	fd = df(f,x);	
	fp = f[x==p0] ;
	fpd = fd[x==p0] ;
	cout << "f(x) = " << f <<endl;
	cout << "f(x)[x==pi/4] = " << f[x==p0] <<endl;
	cout << endl;
	cout << "f'(x) = " << fd <<endl;
	cout << "f'(x)[x==pi/4] = " << fd[x==p0] <<endl;

	cout << endl;
	cout << setw(6) << "n" << "\t\t" << "p_{n}"  << "\t\t\t\t" << "p_{n} in Degree" << "\n";
	cout << setprecision(14) << setw(6) << "0" << "\t\t" << p0 <<  "\t\t" << p0*radtodeg << "\n";	
	for (int i = 1; i <=N; i++)
	{
		fp = f[x==p0] ;
		fpd = fd[x==p0] ;
		pn = p0 - (fp/fpd);

		cout << setprecision(14) << setw(6) << i << "\t\t" << pn << "\t\t" << pn*radtodeg << "\n";
		double err = p0-pn;
		if (abs(err) < pow(10,-5))
		{
			cout << "The procedure was successful." << endl;			
			break;
		}
		p0 = pn;
	}

	return 0;
}