#pragma once
#ifndef Algorithm_h
#define Algorithm_h


#ifndef PI
#define PI 3.141592653589793238462643383279502884197169
#endif // !PI


#include "MatrixTemplate.h"


template <class Type> 
bool Matrix_LU(Type *L, Type *U, const Type *K, int n) {
	//对方阵K进行LU分解.分解失败返回False.成功返回True以及分解得到的L与U
	int i, j, a, b, c, d;
	Type temp;
	for (i = 0, a = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			L[a + j] = U[a + j] = 0;
		}
		U[a + i] = 1;
		a += n;
	}
	for (j = 0, d = 0; j<n; j++)
	{
		for (i = j, b = d; i<n; i++)
		{
			temp = 0;
			a = 0, c = j;
			while (a<j)
			{
				temp += L[b + a] * U[c];
				c += n;
				a++;
			}
			L[b + j] = K[b + j] - temp;
			b += n;
		}
		i = j + 1;
		while (i<n)
		{
			temp = 0;
			a = 0, c = i;
			while (a<j)
			{
				temp += L[d + a] * U[c];
				a++;
				c += n;
			}
			if (L[d + j] == 0)
			{
				cout << "L[d + j] == 0 fail in Matrix_LU()" << endl;
				system("pause");
				return false;
			}
			U[d + i] = (K[d + i] - temp) / L[d + j];
			i++;
		}
		d += n;
	}
	return true;
}


template <class Type> 
bool Matrix_Inv(Type *InvK, Type *K, int n, Type *buf) {
	//采用LU分解方法求方阵K的逆InvK,K[n][n]
	if (1 == n)
	{
		if (K[0] == 0)
		{
			cout << "K[0]==0 in Matrix_Inv()" << endl;
			system("pause");
			return false;
		}
		else
		{

			InvK[0] = 1 / K[0];
		}
	}
	else if (n<1)
	{
		cout << "n<1 in Matrix_Inv()" << endl;
		system("pause");
		return false;
	}
	else
	{
		int i, j, a, b;
		Type *d, *x, *e, temp, *L, *U;
		a = n*n;
		L = buf;
		U = &buf[n*n];
		if (Matrix_LU(L, U, K, n))
		{
			d = &buf[2 * n*n];
			x = &buf[2 * n*n + n];
			e = &buf[2 * n*n + 2 * n];

			for (i = 0; i<n; i++)
			{
				x[i] = d[i] = 0;
			}
			for (i = 0; i<n; i++)
			{
				for (j = 0; j<n; j++)
				{
					e[j] = 0;
				}
				e[i] = 1;
				j = 0;
				b = 0;
				while (j<n)
				{
					temp = 0;
					a = 0;
					while (a<j)
					{
						temp += d[a] * L[b + a];
						a++;
					}
					d[j] = e[j] - temp;
					d[j] /= L[b + j];
					j++;
					b += n;
				}
				j = n - 1;
				b -= n;
				while (j>-1)
				{
					temp = 0;
					a = j + 1;
					while (a<n)
					{
						temp += U[b + a] * x[a];
						a++;
					}
					x[j] = d[j] - temp;
					x[j] /= U[b + j];
					j--;
					b -= n;
				}
				for (j = 0, b = i; j<n; j++)
				{
					InvK[b] = x[j];
					b += n;
				}
			}

		}
		else
		{
			cout << "Matrix_LU() fail in Matrix_Inv()" << endl;
			system("pause");
			return false;
		}

	}
	return true;
}


template <class Type> 
void InvDenFast(Type *InvA, Type &err, int &Niter, Type *InvA0, Type *A, Type *R0, int *ipiv, int N, int NiterMax, const Type errMax) {
	//利用迭代法计算密矩阵的逆
	Type aa, bb;

	//计算初始残差矩阵R0
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			aa = 0;
			for (int k = 0; k < N; ++k) {
				aa += A[i*N + k] * InvA0[k*N + j];
			}
			R0[i*N + j] = -aa;
		}
		R0[i*N + i] += 1;
	}

	//计算R0无穷范数
	aa = 0;
	for (int i = 0; i < N; i++) {
		bb = 0;
		for (int j = 0; j < N; j++) {
			bb += Abs(R0[i*N + j]);
		}
		aa = Max(aa, bb);//无穷范数
	}

	//根据R0无穷范数判断
	//cout << aa;
	if (aa >= 1 || isnan(aa)) {
		//如果范数过大就直接用直接解法
		for (int i = 0; i < N*N; i++) InvA[i] = A[i];
		//Inverse(InvA, ipiv, N);
		Matrix_Inv(InvA, A, N, R0);
		Niter = 0;
		err = 0;
		return;
	}
	else {
		//如果范数比较小，用迭代法求解
		err = 0;
		for (Niter = 1; Niter <= NiterMax; Niter++) {
			for (int i = 0; i < N; i++) R0[i*N + i] += 1;

			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					aa = 0;
					for (int k = 0; k < N; k++) {
						aa += InvA0[i*N + k] * R0[k*N + j];
					}
					InvA[i*N + j] = aa;
				}
			}


			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					aa = 0;
					for (int k = 0; k < N; k++) {
						aa += A[i*N + k] * InvA[k*N + j];
					}
					R0[i*N + j] = -aa;
				}
				R0[i*N + i] += 1;
			}

			err = 0;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					err += Abs(R0[i*N + j]);
				}
			}
			err /= N;
			if (err < errMax) break;

			for (int i = 0; i < N*N; i++) InvA0[i] = InvA[i];
		}
	}
}


template <class Type> 
Type Atan2(const Type &m, const Type &n) {
	Type aoa = 0;
	if (((m > 0) && (n > 0)) || ((m < 0) && (n > 0))) { aoa = atan(m / n); }
	else if ((m > 0) && (n < 0)) { aoa = atan(m / n) + PI; }
	else if ((m < 0) && (n < 0)) { aoa = atan(m / n) - PI; }
	else if ((m == 0) && (n > 0)) { aoa = 0.0; }
	else if ((m == 0) && (n < 0)) { aoa = PI; }
	else if ((m > 0) && (n == 0)) { aoa = PI / 2; }
	else if ((m < 0) && (n == 0)) { aoa = -PI / 2; }
	else { aoa = 0.0; }
	return aoa;
}


template<class Type> 
void Cross(Type A[3], const Type &x1, const Type &y1, const Type &z1, const Type &x2, const Type &y2, const Type &z2) {
	A[0] = y1*z2 - y2*z1;
	A[1] = z1*x2 - x1*z2;
	A[2] = x1*y2 - y1*x2;
}


template<class Type> 
void Cross(Type A[3], const Type x[3], const Type y[3]) {
	A[0] = x[1] * y[2] - y[1] * x[2]; //y1*z2 - y2*z1;
	A[1] = x[2] * y[0] - x[0] * y[2]; //z1*x2 - x1*z2;
	A[2] = x[0] * y[1] - x[1] * y[0]; //x1*y2 - y1*x2;

}


template<class Type> 
Type Norm(const Type &x, const Type &y, const Type &z) {
	return sqrt(x*x + y*y + z*z);
}


template<class Type> 
Type Dot(const Type &x1, const Type &y1, const Type &z1, const Type &x2, const Type &y2, const Type &z2) {
	return x1*x2 + y1*y2 + z1*z2;
}


template<class Type> 
Type Dot(const Type x[3], const Type y[3]) {
	return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}




template<class _Ty>
void KK2SpM(SpMtrx<_Ty> &Kp, Matrix3<_Ty> &KK, Matrix2<int> &t, int SpMtrixType) {
	int NI = KK.NI;
	int NE = KK.NK;
	int ep = t.NJ;
	int D = NI / ep;
	Matrix1<int> nk(NI);
	Matrix3<int> ikm, jkm;

	/*time_t t0, t1;
	t0 = clock();*/

	ikm.allocate(NI, NI, NE);
	jkm.allocate(NI, NI, NE);

	for (int kk = 1; kk <= NE; kk++) {
		for (int ii = 1; ii <= ep; ii++) {
			for (int jj = 1; jj <= D; jj++) {
				nk((ii - 1)*D + jj - 1) = D*(t(kk - 1, ii - 1) - 1) + jj;
			}
		}

		for (int ii = 0; ii < NI; ii++) {
			for (int jj = 0; jj < NI; jj++) {
				ikm(ii, jj, kk - 1) = nk(jj);
				jkm(ii, jj, kk - 1) = nk(ii);
			}
		}
	}

	//nk.output();
	//ikm.output();
	//jkm.output();

	//t1 = clock();	cout << t1 - t0 << endl; t0 = t1;

	int N = 0;
	for (int i = 0; i < t.Nv; i++) {
		N = Max(N, t.v_p[i]);
	}
	//t1 = clock();	cout << t1 - t0 << endl; t0 = t1;
	xij2SpM(Kp, KK.v_p, ikm.v_p, jkm.v_p, N*D, NI*NI*NE, SpMtrixType);
	//t1 = clock();	cout << t1 - t0 << endl; t0 = t1;
	//Kp.refine();//不需要
	//t1 = clock();	cout << t1 - t0 << endl; t0 = t1;
}


template<class _Ty>
void xij2SpM(SpMtrx<_Ty> &Kp, _Ty *X, int *ix, int *jx, int N, int Nv, int SpMtrixType)
{
	//三元数矩阵转CSR转*对称*稀疏矩阵
	//输入的参数X, ix, jx都将被破坏


	//首先把ix和jx和X相同元素进行缩并

	/*time_t t0, t1;
	t0 = clock();*/

	int i1, i2, i, j;

	for (i = 0; i < Nv; i++) {
		if (X[i] == 0 || ix[i]>jx[i]) {
			ix[i] = 0;
			jx[i] = 0;
		}
	}
	//去除0元素
	int Nv1 = 0;
	if (SpMtrixType == 0) {
		//对称矩阵只保留上半部
		for (int i = 0; i < Nv; i++) {
			if (ix[i] <= jx[i] && ix[i] != 0) {
				ix[Nv1] = ix[i];
				jx[Nv1] = jx[i];
				X[Nv1] = X[i];
				Nv1++;
			}

		}
	}
	else if (SpMtrixType == 1) {
		//非对称矩阵
		for (int i = 0; i < Nv; i++) {
			if (ix[i] != 0) {
				ix[Nv1] = ix[i];
				jx[Nv1] = jx[i];
				X[Nv1] = X[i];
				Nv1++;
			}

		}
	}
	Nv = Nv1;
	//int N1;
	//	DuplicateSum(N1, X, ix, jx, Nv);
	//cout << Nv<<endl;
	for (i = 0; i < Nv; i++) {
		if (ix[i] != 0) {

			i1 = ix[i];
			i2 = jx[i];

			for (j = i + 1; j < Nv; j++) {
				if (i1 == ix[j] && i2 == jx[j]) {
					X[i] += X[j];
					X[j] = 0;
					ix[j] = 0;
					jx[j] = 0;
				}
			}
		}
		/*if (i%( Nv / 1000)==0){
		cout << i / (Nv / 1000)<<" / 1000" << endl;
		}*/
	}

	//t1 = clock();	cout << t1 - t0 << endl; t0 = t1;



	Kp.allocate(N, Nv1);


	//统计每一行的元素个数，暂时存到Kp.ia
	for (int i = 0; i < N + 1; i++) {
		Kp.ia[i] = 0;
	}

	for (int i = 0; i < Nv; i++) {
		if (ix[i] <= jx[i] && ix[i] != 0) {
			Kp.ia[ix[i]]++;
		}
	}

	//生成ia
	Kp.ia[0] = 1;
	for (int i = 1; i < N + 1; i++) {
		Kp.ia[i] += Kp.ia[i - 1];
	}

	//生成ja和v
	for (int i = 0; i < Nv1; i++) {
		Kp.ja[i] = 0;
		Kp.v[i] = 0;
	}
	//Kp.Output("_Kp.output");

	if (SpMtrixType == 0) {
		//对称矩阵只保留上半部
		for (int i = 0; i < Nv; i++) {
			if (ix[i] <= jx[i] && ix[i] != 0) {
				for (int ii = Kp.ia[ix[i] - 1]; ii < Kp.ia[ix[i]]; ii++) {
					//cout << endl;
					//cout << i << endl;
					//cout << ix[i] << endl;
					//cout << Kp.ia[ix[i] - 1] << endl;
					//cout << endl;
					if (Kp.ja[ii - 1] == 0) {
						Kp.ja[ii - 1] = jx[i];
						Kp.v[ii - 1] = X[i];
						break;
					}
				}
			}
		}
	}
	else if (SpMtrixType == 1) {
		//非对称矩阵
		for (int i = 0; i < Nv; i++) {
			if (ix[i] != 0) {
				for (int ii = Kp.ia[ix[i] - 1]; ii < Kp.ia[ix[i]]; ii++) {
					if (Kp.ja[ii - 1] == 0) {
						Kp.ja[ii - 1] = jx[i];
						Kp.v[ii - 1] = X[i];
						break;
					}
				}
			}
		}
	}
	//对ja和v重排序
	int jb;
	double vb;
	for (int i = 1; i < N; i++) {
		for (int ii = Kp.ia[i - 1]; ii < Kp.ia[i]; ii++) {
			for (int jj = Kp.ia[i - 1]; jj < ii; jj++) {
				if (Kp.ja[jj - 1] > Kp.ja[ii - 1]) {
					jb = Kp.ja[jj - 1];
					Kp.ja[jj - 1] = Kp.ja[ii - 1];
					Kp.ja[ii - 1] = jb;

					vb = Kp.v[jj - 1];
					Kp.v[jj - 1] = Kp.v[ii - 1];
					Kp.v[ii - 1] = vb;
				}
			}

		}

	}
	Kp.SpMtxType = SpMtrixType;
}


template<class _Ty>
void SpM2Mtr(Matrix2<_Ty> &A, SpMtrx<_Ty> &Kp, const int ni, const int nj)
{
	A.allocate(ni, nj);
	for (int j = 0; j < nj; ++j)
		for (int i = 0; i < ni; ++i)
			A(i, j) = Kp(i + 1, j + 1);
}


template<class _Ty>
bool Msolver(SpMtrx<_Ty> &K, Matrix1<_Ty> &F, Matrix1<_Ty> &X)
{
	__int64 pt[64];
	int maxfct = 1, mnum = 1, mtype, phase, error, msglvl;
	int iparm[64];
	int perm;
	int nrhs = 1;

	//Initiliaze the internal solver memory pointer.This is only necessary for the FIRST call of the PARDISO solver.
	for (int i = 0; i < 64; i++) {
		iparm[i] = 0;
		pt[i] = 0;
	}

	//ompnum = omp_get_thread_num();
	iparm[0] = 1;//no solver default
	iparm[1] = 2;//fill - in reordering from METIS
	iparm[2] = 0;//numbers of processors, value of OMP_NUM_THREADS
	iparm[3] = 0;//no iterative - direct algorithm
	iparm[4] = 0;//no user fill - in reducing permutation
	iparm[5] = 0;//= 0 solution on the first n compoments of x
	iparm[6] = 16;//default logical fortran unit number for output
	iparm[7] = 9;//numbers of iterative refinement steps
	iparm[8] = 0;//not in use
	iparm[9] = 13;//perturbe the pivot elements with 1E-13
	iparm[10] = 1;//use nonsymmetric permutation and scaling MPS
	iparm[11] = 0;//not in use
	iparm[12] = 0;//not in use
	iparm[13] = 0;//Output: number of perturbed pivots
	iparm[14] = 0;//not in use
	iparm[15] = 0;//not in use
	iparm[16] = 0;//not in use
	iparm[17] = -1;//Output : number of nonzeros in the factor LU
	iparm[18] = -1;//Output : Mflops for LU factorization
	iparm[19] = 0;//Output : Numbers of CG Iterations
	error = 0;//initialize error flag
	msglvl = 0;//don't print statistical information
	if (K.SpMtxType == 0) {
		mtype = -2;//real and symmetric indefinite
	}
	else if (K.SpMtxType == 1) {
		mtype = 11;//real and nonsymmetric
	}
	phase = 13;// Analysis, numerical factorization, solve, iterative refinement
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &K.N, K.v, K.ia, K.ja, &perm, &nrhs, iparm, &msglvl, F.v_p, X.v_p, &error);
	if (error != 0) {
		K.Output("K.output");
		F.output("F.output", 4);
		cout << "ERROR1 IN SolverMKL was detected: " << error;
		cout << endl;
		getchar();
	}

	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &K.N, K.v, K.ia, K.ja, &perm, &nrhs, iparm, &msglvl, F.v_p, X.v_p, &error);

	if (error != 0) {
		cout << "ERROR2 IN SolverMKL was detected: " << error;
		cout << endl;
		getchar();
	}
	return (error == 0);
}


template<class _Ty>
bool Msolver(SpMtrx<_Ty> &K, Matrix1<_Ty> &X)
{
	return Msolver(SpMtrx<_Ty> &K, Matrix1<_Ty> &X, Matrix1<_Ty> &X);
}


template<class _Ty>
void Msolver(SpMtrx<_Ty> &K, Matrix1<_Ty> &F, Matrix1<_Ty> &X, int solveid)
{
	//用OpenMP的方式单核(solveid = 0)求解,然后广播给所有进程
	//int solveid = 0;//求解的核
	//int myid, NumCore;
	//MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	//MPI_Comm_size(MPI_COMM_WORLD, &NumCore);

	//NumCore = min(NumCore, 4);//一般不超过4个核

	//if (myid == solveid) {
	//	int NumCore0 = omp_get_thread_num();//原始核数
	//	NumCore0 = max(NumCore0, 1);
	//	//cout << "Solve MPI : " << myid << " " << NumCore << NumCore0 << endl;
	//	omp_set_num_threads(NumCore);//新的核数
	//	Msolver(K, F, X);
	//	omp_set_num_threads(NumCore0);//复原
	//}
}


void Msolver(double *A, double *b, const int n, const int nrhs = 1);


void Msolver(float *A, float *b, const int n, const int nrhs = 1);


void SolveSpEig(SpMtrx<double> &K, SpMtrx<double> &M, Matrix1<double> &lam, Matrix2<double> &Vn, int Nwmax, double emin, double emax, int &Nw);



#endif // !Algorithm_h
