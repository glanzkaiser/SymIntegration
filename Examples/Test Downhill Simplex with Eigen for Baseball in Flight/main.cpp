// g++ -o result main.cpp -lsymintegration
// Merci beaucoup Freya..
// C++ program to plot the particular and general solution of linear system Ax=b and Ax=0

#include <iostream>
#include <iomanip> // to declare the manipulator of setprecision()
#include <fstream>
#include <bits/stdc++.h> //for setw(6) at display() function
#include <vector> // For std::vector (example container)
#include "symintegrationc++.h"
#include <algorithm> // For std::sort
#include <eigen3/Eigen/Dense> // header file

using namespace std;

double division(double x, double y)
{
	return x/y;
}
Symbolic divisionsym(Symbolic x, Symbolic y)
{
	return x/y;
}
// Driver code
int main(int argc, char** argv)
{
	float μ = 0.05;
	int num_insidecontraction = 0, num_outsidecontraction = 0, num_expansion = 0, num_reflection = 0, num_shrink = 0;
	Symbolic x("x"), y("y");
	Symbolic f_x0, f_x1, f_x2, f_xr, f_xe, f_xc;
	double fx0, fx1, fx2, fxr, fxe, fxc, fx0_final, fx1_final, fx2_final;
	
	int N = 20; // Number of iteration
	int n = 2; // number of  simplex points

	Eigen::Vector2d x0b, x1b, x2b, x0, x1, x2;
	x0b << 145.3, 0.5;
	x1b << x0b[0]+μ, x0b[1]+0 ;
	x2b << x0b[0]+0, x0b[1]+μ ;

	//cout <<"Vector x_{0}: \n" << x0b <<endl;
	//cout <<"Vector x_{1}: \n" << x1b <<endl;
	//cout <<"Vector x_{2}: \n" << x2b <<endl;

	Eigen::Vector2d c, xr, xe, xdb, xd, xc, s0, s1;

	Equations rules_x0b = (x == x0b[0],
                        		y == x0b[1]);
	Equations rules_x1b = (x == x1b[0],
                        		y == x1b[1]);
	Equations rules_x2b = (x == x2b[0],
                        		y == x2b[1]);
	Symbolic f = -(800*( ln(1-70/(x*cos(y))) ) + 56000/(x*cos(y)) + 350*sin(y)*(cos(y)^(-1)) - 7);
	Symbolic f_x0b = evalf(f.subst_all(rules_x0b),1,1);
	Symbolic f_x1b = evalf(f.subst_all(rules_x1b),1,1);
	Symbolic f_x2b = evalf(f.subst_all(rules_x2b),1,1);
	double fx0b = f_x0b;
	double fx1b = f_x1b;
	double fx2b = f_x2b;
	// Sort in ascending order
	vector<double> v_fxyb = {fx0b, fx1b, fx2b};
	sort(v_fxyb.begin(), v_fxyb.end()); 
	
		if(v_fxyb[0] == fx0b) // sorting x0, x1, x2
		{
			x0 = x0b;
			if (v_fxyb[1] == fx1b)
			{
				x1 = x1b;
				x2 = x2b;
			}
			else if (v_fxyb[1] == fx2b)
			{
				xdb = x1b;
				x1 = x2b;
				x2 = xdb;
			}
		}
		else if(v_fxyb[0] == fx1b)
		{
			xdb = x0b;
			x0 = x1b;
			if (v_fxyb[1] == fx0b)
			{
				x1 = xdb;
				x2 = x2b;
			}
			else if (v_fxyb[1] == fx2b)
			{
				x1 = x2b;
				x2 = xdb;
			}			
		}
		else if(v_fxyb[0] == fx2b)
		{
			xdb = x0b;
			x0 = x2b;
			if (v_fxyb[1] == fx0b)
			{
				x2 = x1b;
				x1 = xdb;
			}
			else if (v_fxyb[1] == fx1b)
			{
				x1 = x1b;
				x2 = xdb;
			}			
		}
		//cout << "After sorting, \nx_{0}: \n" << x0 <<endl;
		//cout << "x_{1}: \n" << x1 <<endl;
		//cout << "x_{2}: \n" << x2 <<endl;
	// end of sorting in ascending order

	float α = 1, β = 0.5, γ = 2, δ = 0.5;
	
	cout << "f(x,y) = " << f << endl;
	// Start the downhill simplex computation
	for (int i = 1; i <= N; i++) 
	{
		c = division(1,n)*(x0+x1);
		
		xr = c+ α*(c-x2);

		Equations rules_xr = (x== xr[0],
                        		y == xr[1]);
		Equations rules_x0 = (x == x0[0],
                        		y == x0[1]);
		Equations rules_x1 = (x == x1[0],
                        		y == x1[1]);
		Equations rules_x2 = (x == x2[0],
                        		y == x2[1]);
		f_x0 = evalf(f.subst_all(rules_x0),1,1);
		f_x1 = evalf(f.subst_all(rules_x1),1,1);
		f_x2 = evalf(f.subst_all(rules_x2),1,1);
		f_xr = evalf(f.subst_all(rules_xr),1,1);
		
		fx0 = f_x0;
		fx1 = f_x1;
		fx2 = f_x2;
		fxr = f_xr;
		if (fxr < fx1 && fx0 <= fxr) // I. Reflection
		{
			//cout << "*** Reflection *** " << endl;
			x2 = xr;
			num_reflection += 1 ; 
		}
		else if (fxr < fx0) // continue to explore and do expansion
		{
			//cout << "*** Expansion *** " << endl;
			
			xe = c + γ*(xr - c);
			Equations rules_xe = (x == xe[0],
                        		y == xe[1]);
			f_xe = evalf(f.subst_all(rules_xe),1,1);
			fxe = f_xe;
			
			// Greedy minimization
			if (fxe < fxr )
			{
				//cout << "*** Greedy minimization *** " << endl;
				x2 = xe;
				num_expansion += 1 ;
			}
			else if (fxr <= fxe)
			{
				//cout << "*** Greedy minimization *** " << endl;
				x2 = xr;
				num_reflection += 1 ; 
			}
			// end of Greedy minimization
		}
		else if (fxr >= fx1 && fxr <fx2 ) // Do outside contraction
		{
			
			xc = c + β*(xr-c);
			Equations rules_xc = (x == xc[0],
                       		y == xc[1]);
			f_xc = evalf(f.subst_all(rules_xc),1,1);
			fxc = f_xc;
			if (fxc <= fxr)
			{
				//cout << "*** Outside Contraction *** " << endl;
				x2 = xc; 
				num_outsidecontraction += 1 ;
			}
			else if (fxc > fxr)//  Peform shrink operation
			{
				//cout << "*** Shrink from outside contraction *** " << endl;
					
				s0 = x0 + δ*(x1-x0);
				s1 = x0 + δ*(x2-x0);	

				x1 = s0;
				x2 = s1;
					num_shrink += 1;
				
			}
		}
		else if (fxr >= fx2) // Do Inside contraction
		{
			xc = c + β*(x2-c);
			Equations rules_xc = (x == xc[0],
                       		y == xc[1]);
			f_xc = evalf(f.subst_all(rules_xc),1,1);
			fxc = f_xc;

			if (fxc <= fxr)
			{
				//cout << "*** Inside contraction *** " << endl;
				x2 =xc;
				num_insidecontraction += 1;
			}
			else if (fxc > fxr) //  Peform shrink operation
			{
				//cout << "*** Shrink from inside contraction *** " << endl;
					
				s0 = x0 + δ*(x1-x0);
				s1 = x0 + δ*(x2-x0);	
					
				x1 = s0;
				x2 = s1;
				num_shrink += 1;
			}
			
		} // end of downhill simplex

		/*cout << "\ni = " << i <<endl;
		cout <<"c:" << "\n" << c <<endl;
		cout << "\nx_{0}: \n" << x0 <<endl;
		cout << "x_{1}: \n" << x1 <<endl;
		cout << "x_{2}: \n" << x2 <<endl;
		cout << "x_{r}: \n" << xr <<endl;
		cout << "x_{e}: \n" << xe <<endl;
		cout << "x_{c}: \n" << xc <<endl;

		cout << "\n Total iteration: \t " << N <<endl;
		cout << "Reflection: \t" << num_reflection <<endl;
		cout << "Expansion: \t" << num_expansion <<endl;
		cout << "Outside contraction: \t" << num_outsidecontraction <<endl;
		cout << "Inside contraction: \t" << num_insidecontraction <<endl;
		cout << "Shrink: \t" << num_shrink <<endl;
		
		cout << "\nf(x,y) from x_{0}= " << fx0 << endl;
		cout << "f(x,y) from x_{1}= " << fx1 << endl;
		cout << "f(x,y) from x_{2}= " << fx2 << endl;
		cout << "f(x,y) from x_{r}= " << fxr << endl;
		cout << "f(x,y) from x_{e}= " << fxe << endl;
		cout << "f(x,y) from x_{c}= " << fxc << endl;*/
		
		f_x0 = evalf(f.subst_all(rules_x0),1,1);
		f_x1 = evalf(f.subst_all(rules_x1),1,1);
		f_x2 = evalf(f.subst_all(rules_x2),1,1);
		
		fx0 = f_x0;
		fx1 = f_x1;
		fx2 = f_x2;

		// Sort in ascending order
		vector<double> v_fxy = {fx0, fx1, fx2};
		sort(v_fxy.begin(), v_fxy.end()); 
		// Sort in descending order using a custom comparator
        	// sort(numbers.begin(), numbers.end(), greater<double>());

		if(v_fxy[0] == fx0) // sorting x0, x1, x2
		{
			x0 = x0;
			if (v_fxy[1] == fx1)
			{
				x1 = x1;
				x2 = x2;
			}
			else if (v_fxy[1] == fx2)
			{
				xd = x1;
				x1 = x2;
				x2 = xd;
			}
		}
		else if(v_fxy[0] == fx1)
		{
			xd = x0;
			x0 = x1;
			if (v_fxy[1] == fx0)
			{
				x1 = xd;
				x2 = x2;
			}
			else if (v_fxy[1] == fx2)
			{
				x1 = x2;
				x2 = xd;
			}			
		}
		else if(v_fxy[0] == fx2)
		{
			xd = x0;
			x0 = x2;
			if (v_fxy[1] == fx0)
			{
				x2 = x1;
				x1 = xd;
			}
			else if (v_fxy[1] == fx1)
			{
				x1 = x1;
				x2 = xd;
			}			
		}
		Symbolic f_x0_final = evalf(f.subst_all(rules_x0),1,1);
		Symbolic f_x1_final = evalf(f.subst_all(rules_x1),1,1);
		Symbolic f_x2_final = evalf(f.subst_all(rules_x2),1,1);
		
		fx0_final = f_x0_final;
		fx1_final = f_x1_final;
		fx2_final = f_x2_final;
		//cout << "After sorting, \nx_{0}: \n" << x0 <<endl;
		//cout << "x_{1}: \n" << x1 <<endl;
		//cout << "x_{2}: \n" << x2 <<endl;
	}	// end of for looping
	cout << "\nTotal iteration: \t " << N <<endl;
	cout << "Reflection: \t\t" << num_reflection <<endl;
	cout << "Expansion: \t\t" << num_expansion <<endl;
	cout << "Outside contraction: \t" << num_outsidecontraction <<endl;
	cout << "Inside contraction: \t" << num_insidecontraction <<endl;
	cout << "Shrink: \t\t" << num_shrink <<endl;
	cout << "\nx_{0}: \n" << x0 <<endl;
	cout << "x_{1}: \n" << x1 <<endl;
	cout << "x_{2}: \n" << x2 <<endl;	
	cout << "\nf(x,y) from x_{0}= " << setprecision(12) << fx0_final << endl;
	cout << "f(x,y) from x_{1}= " << setprecision(12) << fx1_final << endl;
	cout << "f(x,y) from x_{2}= " << setprecision(12) << fx2_final << endl;
	return 0;
}