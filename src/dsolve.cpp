/*
   
*/

#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_DSOLVE_DEFINE

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, mu, C("C");
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, r;
		// Case 1 : ay' + ty = b
		eq = (a*x*y + d).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = exp(integrate(-ap*x,x));
		dsol = (integrate(mu*dp,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : aty' + by = ct^2	/ y' = -(a/t)y + b*t	
		eq = ((a/x)*y + d*x).match(fx, (a,d));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), dp = rhs(*i, d);
		mu = x^(-ap);
		dsol = (integrate(mu*dp*x,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : ay' + by = c*exp(d*t)	/ y' = -ay + b*exp(c*t) 	
		eq = (a*y + b*exp(c*x)).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(-ap*x);
		dsol = (integrate(mu*bp*exp(cp*x),x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : ay' + by = c*t+d	/ y' = -ay + b*t + c 	
		eq = (a*y + b*x+c).match(fx, (a,b,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b), cp = rhs(*i,c);
		mu = exp(-ap*x);
		dsol = (integrate(mu*(bp*x+cp),x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 5 :  ay' + by = c / y' = ay +by +c
		eq = (a*y+b).match(fx, (a,b)); // pow(u,r)[r==4]
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
		mu = exp(integrate(-ap,x));
		dsol = (integrate(bp*mu,x))/(mu) + (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 6 :  y' + ay = 0
		eq = (a*y).match(fx, (a,b)); // pow(u,r)[r==4]
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		mu = exp(integrate(-ap,x));
		dsol = (C)/(mu);
		if(df(rhs(*i, a), x) == 0) 
		{
			return dsol;
		}
		} catch(const SymbolicError &se) {}
		}
	}
	return dsol;
}

Symbolic dsolve(const Symbolic &fx, const Symbolic &y, const Symbolic &x, const Symbolic &z) // for 1st order ODE with rate r or k
{
	Symbolic dsol, mu, gs,  C("C"), k("k");
 	
	if(fx != 0)
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d;
		// Will work for this model: y' + ry/b = r/c	/ y' = r/b - ry/c 
		// y' = ry - k	
		// k = (fx.coeff(z,0)).coeff(y,0)	
		mu = exp(-fx.coeff(y*z,1)*x*z);
		gs = z*(fx.coeff(z,1)).coeff(y,0);
		dsol = (integrate(mu*(gs + (fx.coeff(z,0)).coeff(y,0)),x))/(mu) + (C)/(mu);
		
		/*// Case 1 : y' - ry = -k	/ y' = ry - k 	Useless this cannot be matched. Still keep it for now maybe will be of use one day.
		eq = (a*y*z - b*k).match(fx, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i,a), bp = rhs(*i, b);
		mu = exp(ap*z*x);
		dsol = (integrate(mu*(bp*z),x))/(mu) + (C)/(mu);
		
		return 8;
		
		} catch(const SymbolicError &se) {}
		}*/
	}
	return dsol;
}

Symbolic dsolvelogistic(const Symbolic &fx, const Symbolic &y, const Symbolic &y0, const Symbolic &t, const Symbolic &r, const Symbolic &K, const Symbolic &T) // for 1st order ODE logistic growth
{
	Symbolic dsol;
 	
	if(fx.coeff(y,1) == r && fx.coeff(y,2) == -r*(K^(-1)))
 	{
		dsol = (y0*K)/(y0 + (K - y0)*exp(-r*t));
	}
	if(fx.coeff(y,1) == -r && fx.coeff(y,2) == r*(T^(-1)))
 	{
		dsol = (y0*T)/(y0 + (T - y0)*exp(r*t));
	}
	if(fx.coeff(y,1) == -r && fx.coeff(y,2) == r*(T^(-1)) + r*(K^(-1)))
 	{
		cout << "For T < y < K and y > K"<< endl;
		cout << "y(t) = " << (y0*K)/(y0 + (K - y0)*exp(-r*t)) << endl;
		cout << "For 0 < y < T, \ny(t) ="<< endl;
		dsol = (y0*T)/(y0 + (T - y0)*exp(r*t));
	}
	return dsol;
}

Symbolic ivp(const Symbolic &fx, const Symbolic &x, const Symbolic &c, const Symbolic &so)
{
	Symbolic ivpsol, f0,  C("C");
 
	if(fx != 0)
 	{
		// Will work for this model: dS/dt = rS-k	
		f0 = fx[x==0];
		C = solve(f0-so,c).front().rhs ;
		ivpsol = fx[c==C];	
	}
	return ivpsol;
}

Symbolic dsolveseparable(const Symbolic &fdy, const Symbolic &fdx, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, dsol_y, dsol_x, C("C"), c_y;
 	
	if(fdx == 0 && fdy !=0)
 	{
		dsol = C;
	}
	else if(df(fdx,x) == 0 && fdy !=0)
 	{
		c_y = integrate(fdy,y).coeff(y);
		dsol = integrate(fdy,y) - integrate(fdx,x)  - C;
	}
	else if(fdx != 0 && fdy != 0 && fdy.coeff(x,1) == 0 && fdx.coeff(y,1) == 0 && df(fdx,x) != 0 )
 	{
		dsol = integrate(fdy,y) - integrate(fdx,x)  - C; 		
		
	}
	else if(fdy.coeff(x,1) !=0 && fdx.coeff(y,1) !=0 && fdy.coeff(x,2) ==0 && fdx.coeff(y,2) ==0 )
	{
		dsol_y = fdy/x;
		dsol_x = fdx/x ;		
		dsol = dsol_x / dsol_y;
		//dsol = dsol[y*(x^-1)==x];
		//dsol = dsol -(x*dsol_y[y*(x^-1)==x] / dsol_y[y*(x^-1)==x] );
		//dsol = 1/dsol;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x] - x*dsol_y[y*(x^-1)==x];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
		
	}
	else if(fdy.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 )
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}
	else if(fdx.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 && fdy.coeff(x,1) !=0 ) // is this necessary?
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}
	return dsol;
}

void dsolvesecondorderlinear(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			y1 = exp((-b/(2*a))*x) *(cos(D_real*x) + SymbolicConstant::i*sin(D_real*x));
			y2 = exp((-b/(2*a))*x) *(cos(D_real*x) - SymbolicConstant::i*sin(D_real*x));

			cout <<"\nThe general solution is:" << endl;
			cout << "\ny_{1} (t) = " << y1 << endl;
			cout << "\ny_{2} (t) = " << y2 << endl;
		}
	}
}

void ivpsecondorderlinear(double a, double b, double c, const Symbolic &y, const Symbolic &x, double t0, double y0, double dy0)
{
	Symbolic dsol, yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

			double c1_ans = y0;
			double c2_ans = dy0 - (r1*c1_ans);
			cout <<"\nThe solution for the initial value problem is:" << endl;
			yt = yt[c1 == c1_ans, c2 == c2_ans] ;
			cout << yt << endl;

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			cout <<"\nThe general solution is:" << endl;
			cout << yt << endl;

			double c1_ans = divisiond(dy0-(r2*y0),r1-r2)*exp(-r1*t0);
			double c2_ans = divisiond((y0*r1)-dy0,r1-r2)*exp(-r2*t0);
			cout <<"\nThe solution for the initial value problem is:" << endl;
			yt = yt[c1 == c1_ans, c2 == c2_ans] ;
			cout << yt << endl;

			Symbolic dyt = df(yt,x);
			cout << "\ny' = " << dyt << endl;

			double tm = NewtonRaphson(dyt,x,0);

			Equations rules = (  SymbolicConstant::e == exp(1), SymbolicConstant::i == sqrt(-1));
			yt = yt.subst_all(rules);
			double ym = yt[x==tm];

			cout << "\nCritical value t_{m} = " << tm << endl;
			cout << "\nMaximum value y_{m} = " << ym << endl;
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);
			double D_real = divisiond(imag(D_sqrt),2*a); 			
		
			y1 = exp((-b/(2*a))*x) *(cos(D_real*x) + SymbolicConstant::i*sin(D_real*x));
			y2 = exp((-b/(2*a))*x) *(cos(D_real*x) - SymbolicConstant::i*sin(D_real*x));

			cout <<"\nThe complex-valued general solution is:" << endl;
			cout << "\ny_{1} (t) = " << y1 << endl;
			cout << "\ny_{2} (t) = " << y2 << endl;

			Symbolic ut = exp((-b/(2*a))*x) *(cos(D_real*x));
			Symbolic vt = exp((-b/(2*a))*x) *(sin(D_real*x));
			Symbolic d_ut = df(ut,x);
			Symbolic d_vt = df(vt,x);

			cout << "\nu(t) = " << ut << endl;
			cout << "\nv(t) = " << vt << endl;

			double a11 = ut[x==t0];
			double a12 = vt[x==t0];
			double a21 = d_ut[x==t0];
			double a22 = d_vt[x==t0];
			
			vector<vector<double>> A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			A[0][0] = a11;
			A[0][1] = a12;
			A[1][0] = a21;
			A[1][1] = a22;
			vec_b[0][0] = y0;
			vec_b[1][0] = dy0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(A,vec_b,c_solution);
			//printVector(c_solution);

			Symbolic y_solution = exp((-b/(2*a))*x) * (c_solution[0]*(cos(D_real*x)) + c_solution[1]*(sin(D_real*x)));
			cout <<"\nThe real-valued initial value problem solution is:" << endl;
			cout << "\ny (t) = " << y_solution << endl;
			
		}
	}
	 
}

void wronskian(double a, double b, double c, const Symbolic &y, const Symbolic &x, double t0)
{
	Symbolic dsol, yt, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		r1 = divisiond(-b + sqrt((b*b) - (4*a*c)),2*a );
		r2 = divisiond(-b - sqrt((b*b) - (4*a*c)),2*a );
		yt = c1*exp(r1*x) + c2*exp(r2*x);
	}
	cout <<"\nThe general solution is:" << endl;
	cout << yt << endl;
	Symbolic w11 = exp(r1*x);
	Symbolic w12 = exp(r2*x);
	Symbolic w21 = df(w11,x);
	Symbolic w22 = df(w12,x);
	Symbolic W = (w11*w22) - (w12*w21);
	cout <<"\nThe Wronskian is:" << W << endl;

	cout << "\nW(t_{0}) = " << W[x==t0] << endl;
}

void wronskian_fundamentalsetofsolutions(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic dsol, yt, c1("c1"), c2("c2"), k1("k1"), k2("k2");
 	double r1, r2;
	if(a != 0 )
 	{
		r1 = divisiond(-b + sqrt((b*b) - (4*a*c)),2*a );
		r2 = divisiond(-b - sqrt((b*b) - (4*a*c)),2*a );
		yt = c1*exp(r1*x) + c2*exp(r2*x);
	}
	double t0 = 0;
	double y0 = 1;
	double dy0 = 0;
	double c1_ans = divisiond(dy0-(r2*y0),r1-r2)*exp(-r1*t0);
	double c2_ans = divisiond((y0*r1)-dy0,r1-r2)*exp(-r2*t0);
	cout <<"\nThe solution that satisfies the initial value problem y(0)=1 and y'(0)=0:" << endl;
	Symbolic yt_1 = yt[c1 == c1_ans, c2 == c2_ans] ;
	cout << yt_1 << endl;


	double y0_2= 0;
	double dy0_2 = 1;
	double c1_ans_2 = divisiond(dy0_2-(r2*y0_2),r1-r2)*exp(-r1*t0);
	double c2_ans_2 = divisiond((y0_2*r1)-dy0_2,r1-r2)*exp(-r2*t0);

	cout <<"\nThe solution that satisfies the initial value problem y(0)=0 and y'(0)=1:" << endl;
	Symbolic yt_2 = yt[c1 == c1_ans_2, c2 == c2_ans_2] ;
	cout << yt_2 << endl;

	cout <<"\nThe general solution :" << endl;
	Symbolic yt_general = k1*yt_1 + k2*yt_2;
	cout << yt_general << endl;

}


#endif
#endif