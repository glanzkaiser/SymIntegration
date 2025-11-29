// Merci beaucoup Freya et Sentinel
// g++ main.cpp -o result -lsymintegration
#include <iostream>
#include "symintegrationc++.h"

using namespace std;

int main()
{
	Symbolic x("x"), xa, xb, xp, z1, z2, fa, fp, f;
	float a = 1;
	float b = 2;
	int N = 17;

	float p = a + (b-a)/2 ;

	z1 = 3;
	z2 = 2;
	f = pow(x,z1) + 4*pow(x,z2) - 10;
	
	cout << "f(x) = " << f <<endl;
	cout << endl;
	cout << "a = " << a <<endl;
	cout << "b = " << b <<endl;
	cout << "p = " << p <<endl;

	xp = p;
	cout << "f(p) = " << pow(xp,z1) + 4*pow(xp,z2) - 10 <<endl;
	xb = b;
	cout << "f(b) = " << pow(xb,z1) + 4*pow(xb,z2) - 10 <<endl;
	xa = a;
	cout << "f(a) = " << pow(xa,z1) + 4*pow(xa,z2) - 10 << endl;
	cout << "f(a)*f(p) =" << (pow(xa,z1) + 4*pow(xa,z2) - 10) * (pow(xp,z1) + 4*pow(xp,z2) - 10) << endl;
	
	cout << setw(6) << "iteration" << "\t\t" << "a" << "\t\t\t" << "b" << "\t\t\t" << "p" << "\t\t\t" << "f(p)" << "\n";
	for (int i = 1; i <=N; i++)
	{
		p = a + (b-a)/2 ;
		xa = a;
		xp = p;
		fa = pow(xa,z1) + 4*pow(xa,z2) - 10;
		fp = pow(xp,z1) + 4*pow(xp,z2) - 10;	
		double fafp = fa*fp;
	
		if (fafp> 0)
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t" << b << "\t\t\t" << p << "\t\t\t" << pow(xp,z1) + 4*pow(xp,z2) - 10 << "\n";
			a = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}	
		else
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t"<< b << "\t\t\t" << p << "\t\t\t" << pow(xp,z1) + 4*pow(xp,z2) - 10 << "\n";
			b = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}
		
		cout << endl;
	}	
	return 0;
}