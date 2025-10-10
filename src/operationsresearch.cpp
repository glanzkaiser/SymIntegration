/*
   Thanks Freya the Goddess, RK, Sentinel, Berlin, Mother Mary.
*/
#include "symintegral/symintegrationc++.h"

using namespace std;

#ifdef  SYMBOLIC_DEFINE
#ifndef SYMINTEGRATION_CPLUSPLUS_OPERATIONSRESEARCH_DEFINE
#define SYMINTEGRATION_CPLUSPLUS_OPERATIONSRESEARCH_DEFINE

// Factorial and combinations
int factorialint(int n) {
	if (n <= 1) 
	{
		return 1;
	}
	int result = 1;
	for (int i = 2; i <= n; ++i) 
	{
		result *= i;
	}
	return result;
}

int combinationsint(int n, int r) {
	if (r < 0 || r > n) 
	{
		return 0; // Invalid input
	}
	return factorial(n) / (factorial(r) * factorial(n - r));
}

Symbolic simplexmethod(const SymbolicMatrix &A, vector<double> F, int C, int R, int n)
{
	//int nbv = n-R;
	//int combi =combinationsint(n,nbv);
	int m = R; // number of equality constraints = number of absic variables
	
	vector<double> v_objfunction; // Create a vector to store the value of objective function

	Matrix<Symbolic> B_mat(R,C);
	Matrix<Symbolic> C_mat(R,C);
	Matrix<Symbolic> D_mat(R,R);
	Matrix<Symbolic> E_mat(R,R+1);
	
	for(int i = 0; i<R; i++)
	{
		for(int j=0; j<C; j++)
		{
			B_mat[i][j] = A[i][j];
		}
	}

	cout << "A :\n" << B_mat <<endl;
	
	for(int i = 0; i<R; i++)
	{
		C_mat[i][C-2] = B_mat[i][C-2];
		C_mat[i][C-1] = B_mat[i][C-1];	
	}
	
	// Generate all possible combinations of (n,m)
	string bitmask(m, 1); // m leading 1's
	bitmask.resize(n, 0); // n-m trailing 0's

	// print integers and permute bitmask
	do 
	{
		vector<int> v;
		vector<double> x;
		cout << "\n*********************Basic Variables***********************" << endl;
		cout << "***********************************************************" << endl;
		for (int i = 0; i < n; ++i) // [0..N-1] integers
		{
			if (bitmask[i]) 
			{	
				cout << " " << i;
				v.push_back(i);
				for (int j2 = 0; j2 < C-2; j2++) 
				{
					if(i == j2)
					{
						for (int j1 = 0; j1 < R; j1++)
						{					
							C_mat[j1][j2] = B_mat[j1][j2];
						}
					}
				}
			}
		}
		cout << "\n***********************************************************" << endl;
		//cout << "\nMy wife is the most beautiful Goddess, she teaches me this !" << endl;
		cout << "A :\n" << C_mat << endl;
		
		// Gaussian elimination on play
		for (int k=0; k<R; k++)
		{
			int i_max = k;
		
			if (i_max != k)
			{
				for (int m=0; m<C; k++)
				{
					double temp = C_mat[k][m];
					C_mat[k][m] = C_mat[i_max][m];
					C_mat[i_max][m] = temp;
				}
			}
			for (int i=k+1; i<R; i++)
			{
				 // factor f to set current row kth element to 0 and subsequently remaining kth column to 0 
				double f = C_mat[i][k]/C_mat[k][k];
				// subtract fth multiple of corresponding kth row element
				for (int j=k+1; j<C; j++)
				{
					C_mat[i][j] -= C_mat[k][j]*f;
				}
				C_mat[i][k] = 0;
			}
			
		}
		
		cout << "A (in reduced row form) :\n" << C_mat <<endl;
		for (int i = 0; i < R; i++)
		{
			double pivot = C_mat[i][i];	
			if (C_mat[i][i] != 0)
			{
				pivot = C_mat[i][i];					
			}
			else if (C_mat[i][i] == 0)
			{
				int j = 0;
				while (j < C)
				{
					if (C_mat[i][j] !=0)
					{
						pivot = C_mat[i][j];
					}
					j++;
				}					
			}
			for (int j = 0; j < C; j++) 
			{	
				C_mat[i][j] = C_mat[i][j]/pivot;
			}			
		}

		cout << "A (in row reduced echelon form) :\n" << C_mat << endl;

		// Converting A row reduced echelon form into simplified final matrix.
		for(int i = 0; i<R; i++)
		{
			for(int j=0; j<R; j++)
			{
				D_mat[i][j] = C_mat.transpose()[v[i]][j];
			}
		}
		for(int i = 0; i<R; i++)
		{
			for(int j=0; j<R; j++)
			{
				E_mat[i][j] = D_mat.transpose()[i][j];
			}
			E_mat[i][R] = C_mat.transpose()[C-1][i];
		}
		cout << "Simplified final matrix :\n" << E_mat <<endl;

		// Gaussian elimination
		vector<vector<double>> augmentedMatrix(R, vector<double>(R + 1));

		for(int i = 0; i<R; i++)
		{
			for(int j=0; j<=R; j++)
			{
				augmentedMatrix[i][j] = E_mat[i][j];
			}
		}		
	
		for (int i = 0; i < R; ++i) 
		{
			// Partial Pivoting (optional but recommended for stability)
			int pivotRow = i;
			for (int k = i + 1; k < R; ++k) 
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
			for (int k = i + 1; k < R; ++k) 
			{
				double factor = augmentedMatrix[k][i] / augmentedMatrix[i][i];
				for (int j = i; j <= R; ++j) 
				{ // Iterate up to n for the constant term
					augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
				}
			}
		}

		// Back Substitution
		vector<double> solution(n);
		for (int i = R - 1; i >= 0; --i) 
		{
			double sum = 0.0;
			for (int j = i + 1; j < R; ++j) 
			{
				sum += augmentedMatrix[i][j] * solution[j];
			}
			solution[i] = (augmentedMatrix[i][R] - sum) / augmentedMatrix[i][i];
		}	
		// End of Gaussian elimination
		cout <<"Solution:" <<endl;
		for(int i=0; i<R; i++)
		{
			cout << "x_{" << v[i]+1 << "} = " << solution[i] << endl ;
		
		}	
		double value_objectivefunction=0;
		for(int i=0; i<R; i++)
		{
			value_objectivefunction += F[v[i]]*solution[i] ;
		
		}
		cout << "\nValue of the objective function : " << value_objectivefunction <<endl;	
		v_objfunction.push_back(value_objectivefunction);
		// Refresh the matrix to 0s except the last two columns, they are reverted to the original value again
		
		for(int i = 0; i<R; i++)
		{
			for(int j=0; j<C-2; j++)
			{
				C_mat[i][j] = 0;
			}
		}
		for(int i = 0; i<R; i++)
		{
			C_mat[i][C-2] = B_mat[i][C-2];
			C_mat[i][C-1] = B_mat[i][C-1];	
		}
	}// end of do re mi fa sol la ti do
		while (prev_permutation(bitmask.begin(), bitmask.end()));

	// print integers and permute bitmask
	int i2 = 0;
	cout << "\n***********************************************************" << endl;
	cout << "********************End of Simplex Method*******************" << endl;
	cout << "***********************************************************" << endl;
	cout << "\nObjective Function Value" << setw(33) << "Decision Variables" << endl;
	
	do 
	{
		cout << setprecision(2) << v_objfunction[i2] <<  setw(37) ;	
	
	for (int i = 0; i < n; ++i) 
	{
		if (bitmask[i]) 
		{
			cout << "   " << i;
		}
	}
	
	cout << endl;
	i2 = i2+1;
	} 
	while (prev_permutation(bitmask.begin(), bitmask.end()));

	return 0;
}

#endif
#endif