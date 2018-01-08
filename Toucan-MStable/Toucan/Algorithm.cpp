#include "stdafx.h"
#include "Algorithm.h"



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


void SolveSpEig(SpMtrx<double> &K, SpMtrx<double> &M, Matrix1<double> &lam, Matrix2<double> &Vn, int Nwmax, double emin, double emax, int &Nw)
{
	//只能计算对称矩阵
	int N = K.N, m0 = Nwmax;
	//input
	int  fpm[128];
	int  info;
	char uplo;

	//output
	double epsout;
	int loop;
	double *e; e = new double[m0]; //the first m entries of e are eigenvalues foundi in the interval
								   //int m; //The total number of eigenvalues found in the interval[emin, emax]: 0≤m≤m0.
	double *res; res = new double[m0];
	int nn, mm0;
	double coef;
	double *x;
	x = new double[N*m0];
	feastinit(fpm);
	fpm[0] = 1;
	fpm[1] = 20;
	fpm[2] = 10;//误差小于1e-10停止计算
	fpm[3] = 6;
	//fpm[4] = 1;
	if (K.SpMtxType == 0 && M.SpMtxType == 0) {
		uplo = 'U';
	}
	else if (K.SpMtxType == 1 && M.SpMtxType == 1) {
		cout << " Only full Sym Mtrx can be solved!" << endl;
		uplo = 'F';
	}
	else {
		cout << "Wrong in SolveSpEig: The SpMtxType of K and M is not the same!" << endl;
		system("pause");
	}
	mm0 = m0;
	Nw = m0;
	nn = N;
	dfeast_scsrgv(&uplo, &nn, K.v, K.ia, K.ja, M.v, M.ia, M.ja, fpm, &epsout, &loop, &emin, &emax, &mm0, e, x, &Nw, res, &info);

	if (info >= 100 || info < 0) {
		cout << "EIG Error" << endl;
		getchar();
	}

	switch (info)
	{
	case 1:
		cout << "EIG warning : No eigenvalue found in the search interval.";
		//getchar();
		break;
	case 2:
		cout << "EIG warning : No Convergence (number of iteration loops >fpm[3])";
		//getchar();
		break;
	case 3:
		cout << "EIG warning : Size of the subspace m0 is too small (m0<m)";
		//getchar();
		break;
	case 4:
		cout << "EIG warning : Successful return of only the computed subspace after call with fpm[13] = 1";
		//getchar();
		break;


	default:
		break;
	}
	if (lam.Nv == 0)	lam.allocate(Nw);
	if (Vn.Nv == 0) Vn.allocate(N, Nw);

	for (int i = 0; i < Nw; i++) {
		lam(i) = e[i];
	}

	for (int i = 0; i < Nw; i++) {
		for (int j = 0; j < K.N; j++) {
			Vn(j, i) = x[j + i*N];
		}
	}

	for (int i = 0; i < Nw; i++) {
		coef = Abs(Vn(1, i));
		for (int j = 1; j < K.N; j++) {
			coef = Max(coef, Abs(Vn(j, i)));
		}

		coef = Max(coef, 1e-12);

		for (int j = 0; j < K.N; j++) {
			Vn(j, i) /= coef;
		}


	}
	delete[] e, res;
}
