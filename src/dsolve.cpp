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
	else if((fdy)/(y^(-1)) == 1 && x.coeff(x,-1) != 0) // 2tv' - v = 0 
	{
		double coeff_x = x.coeff(x,-1);
		dsol = C*(x^(coeff_x));		
	}
	/*else if(fdx.coeff(x,2) !=0 && fdx.coeff(y,2) !=0 && fdy.coeff(x,1) !=0 ) // is this necessary?
	{
		dsol_y = fdy/(x*x);
		dsol_x = fdx/(x*x) ;		
		dsol = dsol_x / dsol_y;
		dsol_y = dsol_y[y*(x^-1)==x];
		dsol_x = dsol_x[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)] - x*dsol_y[y*(x^-1)==x, (y^2)*(x^-2)==(x^2)];
		dsol = fractionintegrate(dsol_y,dsol_x,x)[x==y*(x^-1)] - ln(x) - C;
		
	}*/
	return dsol;
}

void dsolvesecondorderlinear(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic yt, y1, y2, c1("c1"), c2("c2");
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
	Symbolic yt, y1, y2, c1("c1"), c2("c2");
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

Symbolic wronskian_resultonly(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic yt, c1("c1"), c2("c2");
 	double r1, r2;
	Symbolic W;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = x*exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));	

			W = D_real*exp(2*(-b/(2*a))*x);		
		}
	}

	return W;
}

void wronskian(double a, double b, double c, const Symbolic &y, const Symbolic &x, double t0)
{
	Symbolic W, yt, c1("c1"), c2("c2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = x*exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			Symbolic w11 = exp(r1*x);
			Symbolic w12 = exp(r2*x);
			Symbolic w21 = df(w11,x);
			Symbolic w22 = df(w12,x);
			W = (w11*w22) - (w12*w21);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));			
			W = D_real*exp(2*(-b/(2*a))*x);
		}
	}
	cout <<"\nThe general solution is:" << endl;
	cout << yt << endl;

	cout <<"\nThe Wronskian is: " << W << endl;

	cout << "\nW(t_{0}) = " << W[x==t0] << endl;
}

void wronskian_fundamentalsetofsolutions(double a, double b, double c, const Symbolic &y, const Symbolic &x)
{
	Symbolic yt, c1("c1"), c2("c2"), k1("k1"), k2("k2");
 	double r1, r2;
	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);

		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));
		}
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

#include "polynomial.h"
#include "rational.h"

void nonhomogeneousequationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, Polynomial<double> &rhs_function)
{
	Symbolic A("A"), t("t");
	Symbolic lhs_final, Yt_final;

	Polynomial<double> p = rhs_function;
	cout << "p(t) = " << p << endl;
	int n = rhs_function.terms.front().second; // degree of the polynomial
	if (n >= 0 ) // this is working for Polynomial of degree n, n = positive integer
	{
		Symbolic Yt;
		for (int i = 0; i <= n; ++i)
		{
			Yt += A[i]*(t^(n-i)); // A[i] is symbolic A_{i} so we don't need to declare infinite number of Symbolic A1("A1"),... ,An("An")
		}
		Symbolic dy = df(Yt,t);
		Symbolic ddy = df(dy,t);
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;

		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		
		for (int i = 0; i <= n; ++i)
		{
			Symbolic coeff_ti = lhs_final.coeff(t,n-i); 
			for (int j = 0; j <= n; ++j)
			{
				mat_A[i][j] = coeff_ti.coeff(A[j],1);
			}
		}

		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		for (int i = 0 ;  i < n+1 ; ++i)
		{
			Yt_final += c_solution[i]*(t^(n-i));
		}
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}
	/*if (n == 3 )
	{
		Symbolic Yt = A*t*t*t + B*t*t + C*t + D;
		Symbolic dy = 3*A*t*t + 2*B*t +C;
		Symbolic ddy = 6*A*t + 2*B;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t3 =  lhs_final.coeff(t,3) ;
		Symbolic coeff_t2 =  lhs_final.coeff(t,2) ;
		Symbolic coeff_t1 =  lhs_final.coeff(t,1) ;
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t3.coeff(A,1);
		mat_A[0][1] = coeff_t3.coeff(B,1);
		mat_A[0][2] = coeff_t3.coeff(C,1);
		mat_A[0][3] = coeff_t3.coeff(D,1);
		mat_A[1][0] = coeff_t2.coeff(A,1);
		mat_A[1][1] = coeff_t2.coeff(B,1);
		mat_A[1][2] = coeff_t2.coeff(C,1);
		mat_A[1][3] = coeff_t2.coeff(D,1);
		mat_A[2][0] = coeff_t1.coeff(A,1);
		mat_A[2][1] = coeff_t1.coeff(B,1);
		mat_A[2][2] = coeff_t1.coeff(C,1);
		mat_A[2][3] = coeff_t1.coeff(D,1);
		mat_A[3][0] = coeff_t0.coeff(A,1);
		mat_A[3][1] = coeff_t0.coeff(B,1);
		mat_A[3][2] = coeff_t0.coeff(C,1);
		mat_A[3][3] = coeff_t0.coeff(D,1);

		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0]*(t^3) + c_solution[1]*(t^2)+ c_solution[2]*t + c_solution[3];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}
	if (n == 2 )
	{
		Symbolic Yt = A*t*t + B*t + C;
		Symbolic dy = 2*A*t + B;
		Symbolic ddy = 2*A;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t2 =  lhs_final.coeff(t,2) ;
		Symbolic coeff_t1 =  lhs_final.coeff(t,1) ;
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t2.coeff(A,1);
		mat_A[0][1] = coeff_t2.coeff(B,1);
		mat_A[0][2] = coeff_t2.coeff(C,1);
		mat_A[1][0] = coeff_t1.coeff(A,1);
		mat_A[1][1] = coeff_t1.coeff(B,1);
		mat_A[1][2] = coeff_t1.coeff(C,1);
		mat_A[2][0] = coeff_t0.coeff(A,1);
		mat_A[2][1] = coeff_t0.coeff(B,1);
		mat_A[2][2] = coeff_t0.coeff(C,1);

		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0]*(t^2) + c_solution[1]*t+ c_solution[2];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}
	if (n == 1 )
	{
		Symbolic Yt = A*t + B;
		Symbolic dy = A;
		Symbolic ddy = 0;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t1 =  lhs_final.coeff(t,1) ;
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t1.coeff(A,1);
		mat_A[0][1] = coeff_t1.coeff(B,1);
		mat_A[1][0] = coeff_t0.coeff(A,1);
		mat_A[1][1] = coeff_t0.coeff(B,1);
	
		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0]*(t) + c_solution[1];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}

	if (n == 0 )
	{
		Symbolic Yt = A;
		Symbolic dy = 0;
		Symbolic ddy = 0;
		
		lhs_final = ddy*lhs_a + dy* lhs_b + Yt*lhs_c;
		//cout << lhs_final << endl;
		
		Symbolic coeff_t0 =  lhs_final.coeff(t,0) ;
		vector<vector<double>> mat_A(n+1, vector<double>(n+1));
		vector<vector<double>> vec_b(n+1, vector<double>(1));
		mat_A[0][0] = coeff_t0.coeff(A,1);
	
		int n1 = n;
		for (int i = 0 ;  i < n+1 ; ++i)
		{
			if(p.terms.front().second == n1) // p.terms.front().second = get the degree of the polynomial
			{
				vec_b[i][0] = p.terms.front().first; //  p.terms.front().first = get the coefficient  
				p.terms.pop_front(); // remove value from front of list
			}
			else if(p.terms.front().second != n1)
			{
				vec_b[i][0] = 0;  
			}
			n1 = n1-1;
		} // this is for all my stressful day

		//printMatrix(mat_A);
		//printMatrix(vec_b);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		Yt_final = c_solution[0];
		//printVector(c_solution);
		cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
	}*/
}

void nonhomogeneousequationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const Symbolic &rhs_function, const Symbolic &y, const Symbolic &x)
{
	Symbolic A("A"), B("B");
	
	double c_final, c1, c2, c3, lhs_final, rhs_final;
	Symbolic Yt_final;

	if(rhs_function != 0 )
 	{
		list<Equations> eq;
		list<Equations>::iterator i;
		UniqueSymbol a, b, c, d, f;
		// Case 1 : g(t) = a*exp(b*t)
		eq = (a*exp(b*x)).match(rhs_function, (a,b));
		
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
		
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);
		if(c_final != INFINITY) 
		{
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		if(c_final ==  INFINITY) 
		{
			Symbolic Yt = x*exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			Symbolic c1 = ddy.coeff(exp(bp*x),1);
			Symbolic c2 = dy.coeff(exp(bp*x),1);
			Symbolic c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c ;
			Symbolic subtract = lhs_final;
			lhs_final = lhs_final - subtract.coeff(x*exp(bp*x),1)*x*exp(bp*x) ;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
	
			Yt_final = c_final*Yt;
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 2 : g(t) = exp(b*t)
		eq = (exp(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic bp = rhs(*i, b);
		Symbolic Yt = exp(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(bp*x),1);
		c2 = dy.coeff(exp(bp*x),1);
		c3 = Yt.coeff(exp(bp*x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(bp*x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;
		if(df(rhs(*i, b), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 3 : g(t) = a*exp(t)
		eq = (a*exp(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = exp(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);

		c1 = ddy.coeff(exp(x),1);
		c2 = dy.coeff(exp(x),1);
		c3 = Yt.coeff(exp(x),1);
		lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
		rhs_final = rhs_function.coeff(exp(x),1);
		c_final = divisiond(rhs_final,lhs_final);

		Yt_final = c_final*Yt;

		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 4 : g(t) = a*sin(t)
		eq = (a*sin(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 5 : g(t) = a*cos(t)
		eq = (a*cos(x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a);
		Symbolic Yt = A*sin(x) + B*cos(x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(x),1);
		Symbolic coeff_cos = Ly.coeff(cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 6 : g(t) = a*sin(b*t)
		eq = (a*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 7 : g(t) = a*cos(b*t)
		eq = (a*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
		Symbolic dy = df(Yt,x);
	 	Symbolic ddy = df(dy,x);
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 8 : g(t) = a*exp(t)*sin(b*t)
		eq = (a*exp(x)*sin(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 9 : g(t) = a*exp(t)*cos(b*t)
		eq = (a*exp(x)*cos(b*x)).match(rhs_function, (a,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
		Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
		Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
	 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 10 : g(t) = a*exp(c*t)*sin(b*t)
		eq = (a*exp(c*x)*sin(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(bp*x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 11 : g(t) = a*exp(c*t)*cos(b*t)
		eq = (a*exp(c*x)*cos(b*x)).match(rhs_function, (a,c,b));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp= rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
		Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(bp*x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 12 : g(t) = a*exp(c*t)*sin(t)
		eq = (a*exp(c*x)*sin(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(x),1);
		vec_b[1][0] = 0;
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}

		// Case 13 : g(t) = a*exp(c*t)*cos(t)
		eq = (a*exp(c*x)*cos(x)).match(rhs_function, (a,c));
		for(i=eq.begin(); i!=eq.end(); ++i)
		{
		try {
		Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
		Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
		Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
	 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
		//cout << "\nY(t) = " << Yt << endl;
		//cout << "\nY'(t) = " << dy << endl;
		//cout << "\nY''(t) = " << ddy << endl;
		
		Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
		Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
		Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

		//cout << Ly << endl;
		//cout << coeff_sin << endl;
		//cout << coeff_cos << endl;

		// We use Gaussian elimination here to obtain A and B
		vector<vector<double>> mat_A(2, vector<double>(2));
		vector<vector<double>> vec_b(2, vector<double>(1));
		mat_A[0][0] = coeff_sin.coeff(A,1);
		mat_A[0][1] = coeff_sin.coeff(B,1);
		mat_A[1][0] = coeff_cos.coeff(A,1);
		mat_A[1][1] = coeff_cos.coeff(B,1);
		vec_b[0][0] = 0;
		vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(x),1);
		vector<double> c_solution;
		solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

		//printVector(c_solution);
		Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
		if(df(rhs(*i, a), x) == 0) 
		{
			cout << "\nThe particular solution is\nY(t) = " << Yt_final << endl;
		}
		} catch(const SymbolicError &se) {}
		}
	}
}

void nonhomogeneousequationssolution(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const SymbolicMatrix &Matrix_A, const Symbolic &y, const Symbolic &x)
{ //Code it in 51 minutes, a silly bug on redeclaring Symbolic Yt_final occurs, it should be done in 30 minutes, on February 14th, 2026
	Symbolic A("A"), B("B");
	
	Symbolic Yt_solution;
 	int n_row = Matrix_A.rows();

	for (int i = 0; i < n_row ; ++i)
	{
		Symbolic Yt_final;	
		double c_final, c1, c2, c3, lhs_final, rhs_final;
		Symbolic rhs_function = Matrix_A[i][0];
		cout << "\n***********************************************************"<< endl;
		cout << "\nFor rhs = " << rhs_function << endl;
		if(rhs_function != 0)
	 	{
			list<Equations> eq;
			list<Equations>::iterator i;
			UniqueSymbol a, b, c, d, f;
			// Case 1 : g(t) = a*exp(b*t)
			eq = (a*exp(b*x)).match(rhs_function, (a,b));
			
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i, b);
			Symbolic Yt = exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			c1 = ddy.coeff(exp(bp*x),1);
			c2 = dy.coeff(exp(bp*x),1);
			c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
			Yt_final = c_final*Yt ;

			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n " << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 2 : g(t) = exp(b*t)
			eq = (exp(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic bp = rhs(*i, b);
			Symbolic Yt = exp(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			c1 = ddy.coeff(exp(bp*x),1);
			c2 = dy.coeff(exp(bp*x),1);
			c3 = Yt.coeff(exp(bp*x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
			rhs_final = rhs_function.coeff(exp(bp*x),1);
			c_final = divisiond(rhs_final,lhs_final);
			Yt_final = c_final*Yt ;

			if(df(rhs(*i, b), x) == 0) 
			{
				cout << "\n " << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 3 : g(t) = a*exp(t)
			eq = (a*exp(x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a);
			Symbolic Yt = exp(x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);

			c1 = ddy.coeff(exp(x),1);
			c2 = dy.coeff(exp(x),1);
			c3 = Yt.coeff(exp(x),1);
			lhs_final = c1*lhs_a + c2* lhs_b + c3*lhs_c;
			rhs_final = rhs_function.coeff(exp(x),1);
			c_final = divisiond(rhs_final,lhs_final);
			Yt_final = c_final*Yt ;

			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 4 : g(t) = a*sin(t)
			eq = (a*sin(x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a);
			Symbolic Yt = A*sin(x) + B*cos(x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(x),1);
			Symbolic coeff_cos = Ly.coeff(cos(x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(sin(x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 5 : g(t) = a*cos(t)
			eq = (a*cos(x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a);
			Symbolic Yt = A*sin(x) + B*cos(x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(x),1);
			Symbolic coeff_cos = Ly.coeff(cos(x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(cos(x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(x) + c_solution[1]*cos(x) ;
			
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n " << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 6 : g(t) = a*sin(b*t)
			eq = (a*sin(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(sin(bp*x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
			cout << Yt_final << endl;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 7 : g(t) = a*cos(b*t)
			eq = (a*cos(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*sin(bp*x) + B*cos(bp*x);
			Symbolic dy = df(Yt,x);
		 	Symbolic ddy = df(dy,x);
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(cos(bp*x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(cos(bp*x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*sin(bp*x) + c_solution[1]*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 8 : g(t) = a*exp(t)*sin(b*t)
			eq = (a*exp(x)*sin(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
			Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
		 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(exp(x)*sin(bp*x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 9 : g(t) = a*exp(t)*cos(b*t)
			eq = (a*exp(x)*cos(b*x)).match(rhs_function, (a,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b);
			Symbolic Yt = A*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x);
			Symbolic dy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) ;
		 	Symbolic ddy = A*exp(x)*sin(bp*x) + A*bp*exp(x)*cos(bp*x) + A*bp*exp(x)*cos(bp*x) - A*bp*bp*exp(x)*sin(bp*x) + B*exp(x)*cos(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*exp(x)*sin(bp*x) - B*bp*bp*exp(x)*cos(bp*x) ;
		
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(x)*cos(bp*x),1);

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(exp(x)*cos(bp*x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(x)*sin(bp*x) + c_solution[1]*exp(x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n" << endl;
			}
			} catch(const SymbolicError &se) {}
			} 

			// Case 10 : g(t) = a*exp(c*t)*sin(b*t)
			eq = (a*exp(c*x)*sin(b*x)).match(rhs_function, (a,c,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
			Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(bp*x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 11 : g(t) = a*exp(c*t)*cos(b*t)
			eq = (a*exp(c*x)*cos(b*x)).match(rhs_function, (a,c,b));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), bp = rhs(*i,b), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(bp*x) + B*exp(cp*x)*cos(bp*x);
			Symbolic dy = A*cp*exp(cp*x)*sin(bp*x) + A*bp*exp(cp*x)*cos(bp*x) + B*cp*exp(cp*x)*cos(bp*x) - B*bp*exp(cp*x)*sin(bp*x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(bp*x) + A*cp*bp*exp(cp*x)*cos(bp*x) + A*bp*cp*exp(cp*x)*cos(bp*x) - A*bp*bp*exp(cp*x)*sin(bp*x) + B*cp*cp*exp(cp*x)*cos(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*cp*bp*exp(cp*x)*sin(bp*x) - B*bp*bp*exp(cp*x)*cos(bp*x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(bp*x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(bp*x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(bp*x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(bp*x) + c_solution[1]*exp(cp*x)*cos(bp*x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 12 : g(t) = a*exp(c*t)*sin(t)
			eq = (a*exp(c*x)*sin(x)).match(rhs_function, (a,c));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
			Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = rhs_function.coeff(exp(cp*x)*sin(x),1);
			vec_b[1][0] = 0;
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}

			// Case 13 : g(t) = a*exp(c*t)*cos(t)
			eq = (a*exp(c*x)*cos(x)).match(rhs_function, (a,c));
			for(i=eq.begin(); i!=eq.end(); ++i)
			{
			try {
			Symbolic ap = rhs(*i, a), cp = rhs(*i,c);
			Symbolic Yt = A*exp(cp*x)*sin(x) + B*exp(cp*x)*cos(x);
			Symbolic dy = A*cp*exp(cp*x)*sin(x) + A*exp(cp*x)*cos(x) + B*cp*exp(cp*x)*cos(x) - B*exp(cp*x)*sin(x) ;
		 	Symbolic ddy = A*cp*cp*exp(cp*x)*sin(x) + A*cp*exp(cp*x)*cos(x) + A*cp*exp(cp*x)*cos(x) - A*exp(cp*x)*sin(x) + B*cp*cp*exp(cp*x)*cos(x) - B*cp*exp(cp*x)*sin(x) - B*cp*exp(cp*x)*sin(x) - B*exp(cp*x)*cos(x) ;
			//cout << "\nY(t) = " << Yt << endl;
			//cout << "\nY'(t) = " << dy << endl;
			//cout << "\nY''(t) = " << ddy << endl;
			
			Symbolic Ly = lhs_a*ddy + lhs_b*dy + lhs_c*Yt;
			Symbolic coeff_sin = Ly.coeff(exp(cp*x)*sin(x),1);
			Symbolic coeff_cos = Ly.coeff(exp(cp*x)*cos(x),1);

			//cout << Ly << endl;
			//cout << coeff_sin << endl;
			//cout << coeff_cos << endl;

			// We use Gaussian elimination here to obtain A and B
			vector<vector<double>> mat_A(2, vector<double>(2));
			vector<vector<double>> vec_b(2, vector<double>(1));
			mat_A[0][0] = coeff_sin.coeff(A,1);
			mat_A[0][1] = coeff_sin.coeff(B,1);
			mat_A[1][0] = coeff_cos.coeff(A,1);
			mat_A[1][1] = coeff_cos.coeff(B,1);
			vec_b[0][0] = 0;
			vec_b[1][0] = rhs_function.coeff(exp(cp*x)*cos(x),1);
			vector<double> c_solution;
			solve_nhsystem_resultsonly(mat_A,vec_b,c_solution);

			//printVector(c_solution);
			Yt_final = c_solution[0]*exp(cp*x)*sin(x) + c_solution[1]*exp(cp*x)*cos(x) ;
			if(df(rhs(*i, a), x) == 0) 
			{
				cout << "\n"<< endl;
			}
			} catch(const SymbolicError &se) {}
			}
		}
			cout << "\nThe current solution is\nY(t) = " << Yt_final ;
			Yt_solution += Yt_final;
		}
		cout << "\n***********************************************************"<< endl;
		cout << "\n***********************************************************"<< endl;

		cout << "\nThe particular solution is\nY(t) = " << Yt_solution ;

}

void nonhomogeneousequationssolution_variationofparameters(const Symbolic &lhs_a, const Symbolic &lhs_b, const Symbolic &lhs_c, const Symbolic &rhs_function, const Symbolic &y, const Symbolic &x)
{

	Symbolic Yt, yt, y1, y2, c1("c1"), c2("c2");
 	double r1, r2;
	double a  = lhs_a, b = lhs_b, c = lhs_c;
	Symbolic W = wronskian_resultonly(lhs_a, lhs_b, lhs_c,y,x);

	if(a != 0 )
 	{
		double D = (b*b) - (4*a*c);
		if (D == 0)
		{
			r1 = divisiond(-b, 2*a );
			r2 = divisiond(-b ,2*a );
			yt = c1*exp(r1*x) + c2*x*exp(r2*x);
			y1 = c1*exp(r1*x);
			y2 = c2*x*exp(r2*x);
		}
		if (D > 0)
		{
			r1 = divisiond(-b + sqrt(D),2*a );
			r2 = divisiond(-b - sqrt(D),2*a );
			yt = c1*exp(r1*x) + c2*exp(r2*x);
			y1 = c1*exp(r1*x);
			y2 = c2*exp(r2*x);
		}
		if (D < 0)
		{
			complex<double> Dc(D,0);
			complex<double> D_sqrt = sqrt(Dc);		
			double D_real = divisiond(imag(D_sqrt),2); 	
		
			yt = exp((-b/(2*a))*x) * (c1*(cos(D_real*x)) + c2*(sin(D_real*x)));			
			y1 = exp((-b/(2*a))*x) * (c1*(cos(D_real*x))) ;
			y2 = exp((-b/(2*a))*x) * (c2*(sin(D_real*x))) ;
		}
	}	
	cout << "y_{1} (t) = " << y1/c1 << endl;
	cout << "y_{2} (t) = " << y2/c2 << endl;
	//cout << "W = " << W << endl;
	Symbolic Y1 = (y2/c2)*rhs_function/(W);
	Y1 = simplifybeforeintegrate(Y1,x);
	Y1 = -(y1/c1)*integrate(Y1,x);	
	Symbolic Y2 = (y1/c1)*rhs_function/(W);
	Y2 = simplifybeforeintegrate(Y2,x);
	Y2 = (y2/c2)*integrate(Y2,x);	
	//Symbolic du1 = -((y2/c2)*rhs_function)/W;	
	//Symbolic du2 = ((y1/c1)*rhs_function)/W;	
	//cout << "u_{1}' = " << du1 << endl;
	//cout << "u_{2}' = " << du2 << endl;


	Yt = Y1 + Y2;
	cout << "Y_{t} = " << Yt << endl;
	yt = Yt + y1 + y2;
	cout << "\nThe general solution is\ny_{t} = " << yt << endl;

}
#endif
#endif