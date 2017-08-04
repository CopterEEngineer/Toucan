#pragma once
#ifndef Algorithm_h
#define Algorithm_h


#ifndef PI
#define PI 3.141592653589793238462643383279502884197169
#endif // !PI


template <class Type> bool Matrix_LU(Type *L, Type *U, const Type *K, int n) {
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


template <class Type> bool Matrix_Inv(Type *InvK, Type *K, int n, Type *buf) {
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


template <class Type> void InvDenFast(Type *InvA, Type &err, int &Niter, Type *InvA0, Type *A, Type *R0, int *ipiv, int N, int NiterMax, const Type errMax) {
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


template <class Type> Type Atan2(const Type &m, const Type &n) {
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


template<class Type> void Cross(Type A[3], const Type &x1, const Type &y1, const Type &z1, const Type &x2, const Type &y2, const Type &z2) {
	A[0] = y1*z2 - y2*z1;
	A[1] = z1*x2 - x1*z2;
	A[2] = x1*y2 - y1*x2;
}


template<class Type> void Cross(Type A[3], const Type x[3], const Type y[3]) {
	A[0] = x[1] * y[2] - y[1] * x[2]; //y1*z2 - y2*z1;
	A[1] = x[2] * y[0] - x[0] * y[2]; //z1*x2 - x1*z2;
	A[2] = x[0] * y[1] - x[1] * y[0]; //x1*y2 - y1*x2;

}


template<class Type> Type Norm(const Type &x, const Type &y, const Type &z) {
	return sqrt(x*x + y*y + z*z);
}


template<class Type> Type Dot(const Type &x1, const Type &y1, const Type &z1, const Type &x2, const Type &y2, const Type &z2) {
	return x1*x2 + y1*y2 + z1*z2;
}


template<class Type> Type Dot(const Type x[3], const Type y[3]) {
	return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}


void Msolver(double *A, double *b, const int n, const int nrhs = 1);


void Msolver(float *A, float *b, const int n, const int nrhs = 1);


#endif // !Algorithm_h
