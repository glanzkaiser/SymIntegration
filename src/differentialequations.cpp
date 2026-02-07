/*
   
*/
#include "symintegral/symintegrationc++.h"
#include <cmath> // For erfc and M_SQRT1_2 (or define M_SQRT1_2 if not available)

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DIFFERENTIALEQUATIONS_DEFINE
#define π 3.1415926535897f

Symbolic picarditeration(const Symbolic &F, double p0, Symbolic &t, Symbolic &y, int N) 
{
//	Symbolic t("t"), y("y");
	Symbolic p = F[y==p0];
	cout<< "\nFor dy/dt =" << F << endl;
	cout << "\nΦ_{0}(t) =" << p << endl;	
		
	for(int i=1; i<=N; i++)
	{
		p = integrate(p,t);
		cout << "Φ_{" << i <<"}(t) =" << p << endl;	
		p = F[y==p];
	}
	return p;
}

Symbolic picarditeration_mathematicalinduction(const Symbolic &F, double p0, Symbolic &t, Symbolic &y) 
{
//	Symbolic t("t"), y("y");
	Symbolic p = F[y==p0];
	cout<< "\nFor dy/dt =" << F << endl;
	cout << "\nΦ_{0}(t) =" << p << endl;	
	list<Symbolic> eq1;	
	list<Symbolic> eq2;
	list<Symbolic> eq3;
	list<Symbolic> eq4;
	list<Symbolic>::iterator i;

	for(int i=1; i<2; i++)
	{
		p = integrate(p,t);
		cout << "Φ_{" << i <<"}(t) =" << p << endl;	
		eq1.push_back(p);
		p = F[y==p];
	}
	for(int i=2; i<3; i++)
	{
		p = integrate(p,t);
		cout << "Φ_{" << i <<"}(t) =" << p << endl;	
		eq2.push_back(p);
		p = F[y==p];
	}
	for(int i=3; i<4; i++)
	{
		p = integrate(p,t);
		cout << "Φ_{" << i <<"}(t) =" << p << endl;	
		eq3.push_back(p);
		p = F[y==p];
	}
	for(int i=4; i<5; i++)
	{
		p = integrate(p,t);
		cout << "Φ_{" << i <<"}(t) =" << p << endl;	
		eq4.push_back(p);
		p = F[y==p];
	}
	
	/*for(i=eq1.begin(); i!=eq1.end(); ++i)
	{
		cout << eq1 << endl;
	}
	for(i=eq2.begin(); i!=eq2.end(); ++i)
	{
		cout << eq2 << endl;
	}
	for(i=eq3.begin(); i!=eq3.end(); ++i)
	{
		cout << eq3 << endl;
	}*/

	Symbolic a0 = eq1 ;
	Symbolic a1 = eq2 - eq1 ;
	Symbolic a2 = eq3 - eq2 ;
	Symbolic a3 = eq4 - eq3 ;
	
	cout <<  "\nSequence 1: " << a0 << endl;
	cout <<  "Sequence 2: " << a1 << endl;
	cout <<  "Sequence 3: " << a2 << endl;
	cout <<  "Sequence 4: " << a3 << endl;

	cout <<  "\nSequence 4 / sequence 3: " << a3/a2 << endl;
	cout <<  "Sequence 3 / sequence 2: " << a2/a1 << endl;
	cout <<  "Sequence 2 / sequence 1: " << a1/a0 << endl;

	Symbolic n("n"), s("t");
	Symbolic dsol;
	int npower;
	list<Equations> eq_final;
	list<Equations> eq_final2;
	list<Equations> eq_final3;
	list<Equations>::iterator i2;
	list<Equations>::iterator i3;
	list<Equations>::iterator i4;

	UniqueSymbol a, b;
	eq_final = (a*(t^(b))).match(a3, (a,b));
	for(i2=eq_final.begin(); i2!=eq_final.end(); ++i2)
	{	
		try {
		Symbolic ap = rhs(*i2, a), bp = rhs(*i2, b);
		dsol =ap;
		//cout << ap << endl;
		//cout << bp << endl;

		UniqueSymbol c, d;
		eq_final2 = (c*(t^(d))).match(a2, (c,d));

		for(i3=eq_final2.begin(); i3!=eq_final2.end(); ++i3)
		{
			try {
			Symbolic cp = rhs(*i3, c), dp = rhs(*i3, d);
			dsol =cp;
			//cout << cp << endl;
			//cout << dp << endl;

			UniqueSymbol f, g;
			eq_final3 = (f*(t^(g))).match(a1, (f,g));

			for(i4=eq_final3.begin(); i4!=eq_final3.end(); ++i4)
			{
				try {
				Symbolic fp = rhs(*i4, f), gp = rhs(*i4, g);
				dsol =fp;
				//cout << fp << endl;
				//cout << gp << endl;

				// Create the pattern finding rule
				if (bp - dp == 2 && dp - gp == 2 && cp/ap == 4 && fp/cp == 3 && 1/fp == 2)
				{
					npower = 2;
					cout << "\nΦ_{n}(t) = " << divisions(pow(s,(npower*n)),SymbolicConstant::n_factorial) << endl;
				}
				} catch(const SymbolicError &se) {}
			}	

			} catch(const SymbolicError &se) {}
		}	
		
		return dsol;
		
		} catch(const SymbolicError &se) {}
	}
	
	return p;
}

void adjointequation(const Symbolic &c1, const Symbolic &c2, const Symbolic &c3, Symbolic &x) 
{
	Symbolic μx("μ(x)"), dμx("μ'(x)"), ddμx("μ''(x)"), y ("y"), dy("y'"), ddy("y''") ;
	Symbolic fx ("f(x)"), dfx("f'(x)");
	Symbolic F = c1*μx*ddy + c2*μx*dy + c3*μx*y;
	cout << "\nExact equation :\n" << F << " = 0" << endl;
	
	Symbolic F2 = df(c1,x)*μx*dy + c1*dμx*dy + c1*μx*ddy + dfx*y + fx*dy;
	cout << F2 << " = 0" << endl;

	Symbolic coeff_dy = F2.coeff(dy,1);
	cout << "\nEquate the coefficients :\n" << coeff_dy << " = " << c2 << " * μ(x)" << endl;
	cout << "\nf'(x) = " << c3 << " * μ(x)" << endl;
	
	cout << "\nDifferentiate both sides of the first equation with respect to : " << x << endl;
	Symbolic a = coeff_dy.coeff(dμx,1);
	Symbolic b = coeff_dy.coeff(μx,1);
	Symbolic c = coeff_dy.coeff(fx,1);
	Symbolic d = c2*dμx + df(c2,x)*μx;

	Symbolic F3 = (df(a,x)*dμx) + (a*ddμx) + (df(b,x)*μx) + (b*dμx) + dfx - d;
	cout << F3 << endl;
	
	cout << "\nSubstitute " << c3 << " * μ(x) for f'(x)\n" << endl;
	cout << F3[dfx==c3*μx] << " = 0" << endl;

}

#endif
#endif