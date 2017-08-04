#include "stdafx.h"
#include "Algorithm.h"

#include "MatrixTemplate.h"

void Msolver(double *A, double *b, const int n, const int nrhs) {
	// A will be overwritten by L U.
	// b will be overwritten by inv(A)b
	int info;
	int *ipiv;
	ipiv = new int[n];

	//info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
	dgetrf(&n, &n, A, &n, ipiv, &info);
	if (info != 0) {
		cout << "Wrong LU factorization." << endl;
		system("pause");
	}
	dgetrs("N", &n, &nrhs, A, &n, ipiv, b, &n, &info);
	//info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, nrhs, A, n, ipiv, b, n);
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

	info = LAPACKE_sgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
	if (info != 0) {
		cout << "Wrong LU factorization." << endl;
		system("pause");
	}
	info = LAPACKE_sgetrs(LAPACK_ROW_MAJOR, 'N', n, nrhs, A, n, ipiv, b, n);
	if (info != 0) {
		cout << "Wrong solving Ax = b." << endl;
		system("pause");
	}
	delete[] ipiv;
}
