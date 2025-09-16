/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_NUMERICALMETHOD_DEFINE

Symbolic bisectionmethod(const Symbolic &f, const Symbolic &x, double a, double b, int N)
{
 	Symbolic xa, xb, xp, fa, fp;
	
	float p = a + (b-a)/2 ;

	cout << setw(6) << "iteration" << "\t\t" << "a" << "\t\t\t" << "b" << "\t\t\t" << "p" << "\t\t\t" << "f(p)" << "\n";
	for (int i = 1; i <=N; i++)
	{
		p = a + (b-a)/2 ;
		xa = a;
		xp = p;
		fa = f[x==xa] ;
		fp = f[x==xp] ;
		double fafp = fa*fp;
	
		if (fafp> 0)
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t" << b << "\t\t\t" << p << "\t\t\t" <<  f[x==xp]  << "\n";
			a = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}	
		else
		{
			cout << setw(6) << i << "\t\t\t" << a << "\t\t\t"<< b << "\t\t\t" << p << "\t\t\t" <<  f[x==xp]  << "\n";
			b = p;
			if ((b-a)/2 < pow(10,-4))
			{
				cout << setw(6) << "Procedure completed successfully" << "\n";
				break;
			}
		}
		
		cout << endl;
	}	
	cout << "solution = " << p << endl;
	return 0;
}

Symbolic newtonmethod(const Symbolic &f, const Symbolic &x, const Symbolic &x0, int N)
{
 	Symbolic fd, fp, fpd, pn, p0;

	p0 = x0;
	fd = df(f,x);	
	fp = f[x==p0] ;
	fpd = fd[x==p0] ;
	cout << "\nf(x) = " << f <<endl;
	cout << "f'(x) = " << fd <<endl;

	cout << endl;
	cout << setw(6) << "n" << "\t\t" << "p_{n}"  << "\n";
	cout << setprecision(14) << setw(6) << "0" << "\t\t" << p0 << "\n";	
	for (int i = 1; i <=N; i++)
	{
		fp = f[x==p0] ;
		fpd = fd[x==p0] ;
		pn = p0 - (fp/fpd);

		cout << setprecision(14) << setw(6) << i << "\t\t" << pn << "\n";
		double err = p0-pn;
		if (abs(err) < pow(10,-5))
		{
			cout << "The procedure was successful." << endl;			
			break;
		}
		p0 = pn;
	}
	cout << "solution = "<< endl;
	return pn;
}

Symbolic eulermethod(const Symbolic &f, const Symbolic &y, const Symbolic &x, const Symbolic &y0, const Symbolic &x0, const Symbolic &x1, double h)
{
 	Symbolic tangent, t("t");
	double t_now, y_now, f_now;
	double N1 = (x1-x0)/h;
	int N = N1;
	y_now = y0;
	t_now = x0;

	f_now = f[x==t_now, y==y_now] ;
	tangent = y_now + f_now*(t-t_now);

	cout << "\nf(x) = " << f <<endl;
	
	cout << endl;
	cout << setw(6) << "t" << "\t\t" << "Euler approximation y_{i}" << "\t\t\t" << "Tangent line"  << "\n";
	cout << setprecision(6) << setw(6) << x0 << "\t\t" << y_now << "\t\t\t\t\t" << tangent << "\n";	
	for (int i = x0; i < N;  i++)
	{
		t_now = t_now + h;
		y_now = y_now + f_now*h;
		f_now = f[x==t_now, y==y_now] ;
		tangent = y_now + f_now*(t-t_now);

		cout << setprecision(6) << setw(6) << t_now << "\t\t" << y_now << "\t\t\t\t\t" << tangent << "\n";
		
		
	}
	cout << "\nsolution = "<< endl;
	return y_now;
}



#endif
#endif