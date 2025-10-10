/*
   Thanks Freya the Goddess, RK, Sentinel, Berlin, Mother Mary.
*/
#include "symintegral/symintegrationc++.h"

using namespace std;

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_LINEARALGEBRA_DEFINE

Symbolic gaussianelimination(const SymbolicMatrix &A, int n)
{
	Matrix<Symbolic> B_mat(n,n+1);
	vector<vector<double>> augmentedMatrix(n, vector<double>(n + 1));

	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			augmentedMatrix[i][j] = A[i][j];
			B_mat[i][j] = A[i][j];
		}
	}
	cout << "A :\n" << B_mat <<endl;
	
	// Forward Elimination
	for (int i = 0; i < n; ++i) 
	{
		// Partial Pivoting (optional but recommended for stability)
		int pivotRow = i;
		for (int k = i + 1; k < n; ++k) 
		{
			if (abs(augmentedMatrix[k][i]) > abs(augmentedMatrix[pivotRow][i])) 
			{
				pivotRow = k;
			}
		}
		swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);

		// Check for singular matrix (no unique solution)
		if (abs(augmentedMatrix[i][i]) < 1e-9) // Using a small epsilon
		{ 
			cout << "No unique solution or infinite solutions exist." << endl;
	 		return 0;
		}

		// Eliminate elements below the pivot
		for (int k = i + 1; k < n; ++k) 
		{
			double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
			for (int j = i; j <= n; ++j) 
			{ // Iterate up to n for the constant term
				augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
			}
		}
	}

	// Back Substitution
	vector<double> solution(n);
	for (int i = n - 1; i >= 0; --i) 
	{
		double sum = 0.0;
		for (int j = i + 1; j < n; ++j) 
		{
			sum += augmentedMatrix[i][j] * solution[j];
		}
		solution[i] = (augmentedMatrix[i][n] - sum) / augmentedMatrix[i][i];
	}	

	for(int i = 0; i<n; i++)
	{
		for(int j=0; j<=n; j++)
		{
			B_mat[i][j] = augmentedMatrix[i][j];
		}
	}
	
	cout << "A (in row reduced echelon form) :\n" << B_mat << endl;

	// Print Solution
	cout << "Solution:" << endl;
	for (int i = 0; i < n; ++i) 
	{
	cout << "x" << i + 1 << " = " << fixed << setprecision(5) << solution[i] << endl;
	}
	cout << endl;
	return 0;
}

Symbolic gaussianeliminationtest(const SymbolicMatrix &A, int C, int R)
{
	Matrix<Symbolic> B_mat(R,C);
	
	for(int i = 0; i<R; i++)
	{
		for(int j=0; j<C; j++)
		{
			B_mat[i][j] = A[i][j];
		}
	}
	cout << "A :\n" << B_mat <<endl;
	
	for (int k=0; k<R; k++)
	{
		int i_max = k;
	
		if (i_max != k)
		{
			for (int m=0; m<C; k++)
			{
				Symbolic temp = B_mat[k][m];
				B_mat[k][m] = B_mat[i_max][m];
				B_mat[i_max][m] = temp;
			}
		}
		for (int i=k+1; i<R; i++)
		{
			 // factor f to set current row kth element to 0 and subsequently remaining kth column to 0 
			Symbolic f = B_mat[i][k]/B_mat[k][k];
			// subtract fth multiple of corresponding kth row element
			for (int j=k+1; j<C; j++)
			{
				B_mat[i][j] -= B_mat[k][j]*f;
			}
			B_mat[i][k] = 0;
		}
		
	}
	
	cout << "A (in reduced row form) :\n" << B_mat <<endl;
	for (int i = 0; i < R; i++)
		{
			double pivot = B_mat[i][i];	
			if (B_mat[i][i] != 0)
			{
				pivot = B_mat[i][i];					
			}
			else if (B_mat[i][i] == 0)
			{
				int j = 0;
				while (j < C)
				{
					if (B_mat[i][j] !=0)
					{
						pivot = B_mat[i][j];
					}
					j++;
				}					
			}
			for (int j = 0; j < C; j++) 
			{	
				B_mat[i][j] = B_mat[i][j]/pivot;
			}
		}

	cout << "A (in row reduced echelon form) :\n" << B_mat << endl;

	return 0;
}

#endif
#endif