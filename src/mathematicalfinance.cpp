/*
   
*/
#include "symintegral/symintegrationc++.h"

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_MATHEMATICALFINANCE_DEFINE

double presentvalue(double A, double r, int n)
{
	return A*((1-(1)/(pow(1+r,n)))/r);
}

double bondpricing(double C, double M, double r, int n)
{
	return (C*((1-(1)/(pow(1+r,n)))/r) ) + M/(pow(1+r,n));
}
double amortization(double L, double r, int n)
{
	double M = L*((r*(pow(1+r,n)))/(pow(1+r,n) - 1));
	double interest, principal, endingbalance, v;
	//endingbalance  = L;
	cout << "\nPeriod" << setw(23) << "Payment" << setw(23)<< "Interest" << setw(23) << "Principal" << setw(30) << "Outstanding Loan Balance" << endl;
	for (int i = 1; i <= n ; i++ )
	{
		v = pow(1+r,-1);
		interest = (1 - (pow(v,n-i+1)))*M; // (1-v^n) * M
		principal = M - interest;
		//endingbalance  -= principal ;
		endingbalance = ( ((1-pow(v,n-i+1))/(r)) - pow(v,n-i+1) )*M; // an - v^n
		if(i==n)
		{
			endingbalance  = 0;		
		}
		cout << i << setw(28) << M << setw(23) << interest << setw(23)  << principal << setw(23) << endingbalance << endl; 
	}
	cout << "\nTotal payment:" << endl;
	return M*n;
}

#endif
#endif