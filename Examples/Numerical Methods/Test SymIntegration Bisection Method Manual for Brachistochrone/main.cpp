// Merci beaucoup Freya et Sentinel
// g++ main.cpp -o main -lsymintegration
#include <iostream>
#include "symintegrationc++.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

int main()
{
	// Get starting timepoint
	auto start = high_resolution_clock::now();

	Symbolic x("x"), xa, xb, xp, fa, fp;
	Symbolic θ0("θ0"), g, f,h;
	
	g = 1 - cos(θ0);
	h = 2*θ0 - 2*sin(θ0);
	f = g-h;
	cout << "f = " << f << endl;

	float a = 1;
	float b = 2;
	int N = 17;

	float p = a + (b-a)/2 ;

	cout << endl;
	cout << "a = " << a <<endl;
	cout << "b = " << b <<endl;
	cout << "p = " << p <<endl;

	xp = p;
	cout << "f(p) = " << f[θ0==xp] <<endl;
	xb = b;
	cout << "f(b) = " <<  f[θ0==xb] <<endl;
	xa = a;
	cout << "f(a) = " <<  f[θ0==xa] << endl;
	cout << "f(a)*f(p) =" <<  f[θ0==xa] * f[θ0==xp]  << endl;
	
	cout << setw(6) << "iteration" << "\t\t" << "a" << "\t\t\t" << "b" << "\t\t\t" << "p" << "\t\t\t" << "f(p)" << "\n";
	for (int i = 1; i <=N; i++)
	{
		p = a + (b-a)/2 ;
		xa = a;
		xp = p;
		fa = f[θ0==xa] ;
		fp = f[θ0==xp] ;
		double fafp = fa*fp;
	
		if (fafp> 0)
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t" << b << "\t\t\t" << p << "\t\t\t" <<  f[θ0==xp]  << "\n";
			a = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}	
		else
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t"<< b << "\t\t\t" << p << "\t\t\t" <<  f[θ0==xp]  << "\n";
			b = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}
		
		cout << endl;
	}	

	// Get ending timepoint
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "\nTime taken by function: " << duration.count() << " microseconds" << endl;
	return 0;
}