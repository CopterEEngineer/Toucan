#include "stdafx.h"
#include "Algorithm.h"

#include "MatrixTemplate.h"

void Msolver(double *A, double *b, const int n, const int nrhs) {
	// A will be overwritten by L U.
	// b will be overwritten by inv(A)b
	int info;
	int *ipiv;
	ipiv = new int[n];

	info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);
	if (info != 0) {
		cout << "Wrong LU factorization." << endl;
		system("pause");
	}
	info = LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', n, nrhs, A, n, ipiv, b, n);
	if (info != 0) {
		cout << "Wrong solving Ax = b." << endl;
		system("pause");
	}
	
	if (ipiv != NULL) { delete[] ipiv; }
}


void Msolver(float *A, float *b, const int n, const int nrhs) {
	// A will be overwritten by L U.
	// b will be overwritten by inv(A)b
	int info, *ipiv;
	ipiv = new int[n];

	info = LAPACKE_sgetrf(LAPACK_COL_MAJOR, n, n, A, n, ipiv);
	if (info != 0) {
		cout << "Wrong LU factorization." << endl;
		system("pause");
	}
	info = LAPACKE_sgetrs(LAPACK_COL_MAJOR, 'N', n, nrhs, A, n, ipiv, b, n);
	if (info != 0) {
		cout << "Wrong solving Ax = b." << endl;
		system("pause");
	}
	delete[] ipiv;
}
