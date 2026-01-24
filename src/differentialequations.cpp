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

#endif
#endif