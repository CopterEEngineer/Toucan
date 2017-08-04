#include "stdafx.h"
#include "UnitTest.h"
#include <time.h>


void test_random(myTYPE &x) {

	unsigned int time_ui = unsigned int(time(NULL));

	srand(time_ui);
	cout << "Seed: " << time_ui << endl;
	x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	printf("x = %.4f \n", x);
}


void test_random(Matrix1<myTYPE> &A) {
	myTYPE r;
	int isize = A.Nv;
	unsigned int time_ui = unsigned int(time(NULL));

	srand(time_ui);
	cout << time_ui << endl;
	r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	for (int i = 0; i < isize; ++i) {
		r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		A.v_p[i] = r;
	}

#ifdef _DEBUG
	printf("A: \n");
	A.output(4);
#endif // _DEBUG	
}


void test_random(Matrix2<myTYPE> &A) {
	myTYPE r;
	int isize = A.Nv;
	unsigned int time_ui = unsigned int(time(NULL));

	srand(time_ui);
	cout << time_ui << endl;
	r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	for (int i = 0; i < isize; ++i) {
		r = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		A.v_p[i] = r;
	}
#ifdef _DEBUG
	printf("A: \n");
	A.output(4);
#endif // _DEBUG	
}


void test_findmaxmin(Matrix1<myTYPE> &A) {
	clock_t tStart;
	myTYPE maxvalue, minvalue;

	printf("Max: %.2f, min: %.2f \n", A.findmax(), A.findmin());
	printf("\n");
	tStart = clock();
	for (int i = 0; i<100000; ++i)
	{
		A.findmaxmin(maxvalue, minvalue);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("Max: %.2f, min: %.2f \n", maxvalue, minvalue);
}


void test_findmaxmin(Matrix2<myTYPE> &A) {
	clock_t tStart;
	myTYPE maxvalue, minvalue;

	printf("Max: %.2f, min: %.2f \n", A.findmax(), A.findmin());
	printf("\n");
	tStart = clock();
	for (int i = 0; i<100000; ++i)
	{
		A.findmaxmin(maxvalue, minvalue);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("Max: %.2f, min: %.2f \n\n", maxvalue, minvalue);

	printf("row max: \n");
	A.findmax(0).output(4);
	printf("row min: \n");
	A.findmin(0).output(4);

	printf("col max: \n");
	A.findmax(1).output(4);
	printf("col min: \n");
	A.findmin(1).output(4);

}


void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, myTYPE xs) {
	myTYPE v = 0;
	clock_t tStart;

	printf("X: ");
	X.output(4);
	printf("Y: ");
	Y.output(4);

	tStart = clock();
	for (int i = 0; i < 100000; ++i)
	{
		v = Y.interplinear_fast(X, xs);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("x = %.4f, y = %.4f \n", xs, v);

	tStart = clock();
	for (int i = 0; i < 100000; ++i)
	{
		v = Y.interplinear(X, xs);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("x = %.4f, y = %.4f \n", xs, v);
}


void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, Matrix1<myTYPE> &XN) {
	Matrix1<myTYPE> temp;
	clock_t tStart;

	printf("X= ");
	X.output(4);
	printf("Y= ");
	Y.output(4);

	tStart = clock();
	for (int i = 0; i < 10000; ++i)
	{
		temp = Y.interplinear_fast(X, XN);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

	printf("XN= ");
	XN.output(4);
	printf("YN= ");
	temp.output(4);

	tStart = clock();
	for (int i = 0; i < 10000; ++i)
	{
		temp = Y.interplinear(X, XN);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);


	printf("XN= ");
	XN.output(4);
	printf("YN= ");
	temp.output(4);
}


void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, Matrix2<myTYPE> &Z, myTYPE xs, myTYPE ys) {
	myTYPE v = 0;
	clock_t tStart;

	printf("X= ");
	X.output(4);
	printf("Y= ");
	Y.output(4);
	printf("Z= ");
	Z.output(4);

	tStart = clock();
	for (int i = 0; i < 10000; ++i)
	{
		v = Z.interplinear(X, Y, xs, ys);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("x = %.4f, y = %.4f, z = %.4f \n", xs, ys, v);

	tStart = clock();
	for (int i = 0; i < 10000; ++i)
	{
		v = Z.interplinear_fast(X, Y, xs, ys);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("x = %.4f, y = %.4f, z = %.4f \n", xs, ys, v);
}


void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, Matrix2<myTYPE> &Z, Matrix2<myTYPE> XN, Matrix2<myTYPE> YN)
{
	Matrix2<myTYPE> temp;
	clock_t tStart;

	printf("X= ");
	X.output(4);
	printf("Y= ");
	Y.output(4);
	printf("Z= ");
	Z.output(4);

	printf("XN= ");
	XN.output(4);
	printf("YN= ");
	YN.output(4);

	tStart = clock();
	for (int i = 0; i < 1000; ++i)
	{
		temp = Z.interplinear(X, Y, XN, YN);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("ZN= ");
	temp.output(4);


	tStart = clock();
	for (int i = 0; i < 1000; ++i)
	{
		temp = Z.interplinear_fast(X, Y, XN, YN);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("ZN= ");
	temp.output(4);
}


void test_matrixmultip(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y) {
	myTYPE temp(X.Nv);
	clock_t tStart;
	printf("Matrix1 dot Matrix1: \n");
	int icount = 0;
	const int COUNT = 10000;

	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		temp = X.dot(Y);
		++icount;
	}
	printf("Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("COUNT: %d\n", icount);
	printf("X dot Y = %.4f \n", temp);

	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		temp = X.dotP(Y);
		++icount;
	}
	printf("Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("COUNT: %d\n", icount);
	printf("X dot Y = %.4f \n", temp);
}


void test_matrixmultip(Matrix2<myTYPE> &A, Matrix1<myTYPE> &b) {
	Matrix1<myTYPE> temp(b.Nv);
	clock_t tStart;

	tStart = clock();
	for (int i = 0; i < 10000; ++i)
	{
		temp = A.matrixmultiply(b);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("A * b = \n");
	temp.output(2);

	tStart = clock();
	for (int i = 0; i < 10000; ++i)
	{
		temp = A.matrixmultiplyP(b);
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("A * b = \n");
	temp.output(2);

}


void test_matrixmultip(Matrix2<myTYPE> &X, Matrix2<myTYPE> &Y) {
	Matrix2<myTYPE> temp(X.NI, Y.NJ);
	clock_t tStart;
	int icount = 0;

	tStart = clock();
	for (int i = 0; i < 10000; ++i) {
		temp = X.matrixmultiply(Y);
		++icount;
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("X * Y = \n");
	//temp.output(4);

	icount = 0;
	tStart = clock();
	for (int i = 0; i < 10000; ++i) {
		temp = X.matrixmultiplyT(Y);
		++icount;
	}
	printf("Time taken: %fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	printf("X * Y = \n");
	temp.output("test.txt");
}


void test_vectormultip(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y) {
#define JSIZE 3000
	Matrix1<myTYPE> temp(X.Nv);
	clock_t tStart;
	int icount = 0;
	int jsize = 0;
	const int COUNT = 70000;
	myTYPE XX[JSIZE], YY[JSIZE], ZZ[JSIZE];
	myTYPE *XXX, *YYY, *ZZZ;
	XXX = new myTYPE[JSIZE];
	YYY = new myTYPE[JSIZE];
	ZZZ = new myTYPE[JSIZE];


	jsize = X.Nv;
	/************************************************************************************/
	tStart = clock();
	for (int i = 0; i < X.Nv; ++i) {
		XX[i] = X.v_p[i];
		YY[i] = Y.v_p[i];
		ZZ[i] = 0;
		XXX[i] = XX[i];
		YYY[i] = YY[i];
		ZZZ[i] = 0;
	}
	for (int i = 0; i < COUNT; ++i) {
		//for (int j = 0; j < jsize; ++j) { ZZ[j] = XX[j] * YY[j]; }
		for (int j = 0; j < jsize; ++j) { ZZ[j] = (*(XX + j))*(*(YY + j)); }

	}
	printf("Time taken (Static Vector): %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y: \n");
	cout.precision(2);
	cout << "[";
	for (int i = 0; i < 10; ++i) { cout << ZZ[i] << " "; }
	cout << "]";
	cout << endl << endl;
	/************************************************************************************/
	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		//vdMul(jsize, XX, YY, ZZ);
		//vsMul(jsize, XX, YY, ZZ);
	}
	printf("Time taken (Static Vector MKL): %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y: \n");
	cout.precision(2);
	cout << "[";
	for (int i = 0; i < 10; ++i) { cout << ZZ[i] << " "; }
	cout << "]";
	cout << endl << endl;
	/************************************************************************************/
	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		for (int j = 0; j < jsize; ++j) { ZZZ[j] = XXX[j] * YYY[j]; }
	}
	printf("Time taken (New Vector): %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y: \n");
	cout.precision(2);
	cout << "[";
	for (int i = 0; i < 10; ++i) { cout << ZZZ[i] << " "; }
	cout << "]";
	cout << endl << endl;
	delete[] XXX;
	delete[] YYY;
	delete[] ZZZ;
	/************************************************************************************/
	XXX = new myTYPE[JSIZE];
	YYY = new myTYPE[JSIZE];
	ZZZ = new myTYPE[JSIZE];
	for (int i = 0; i < X.Nv; ++i) {
		XXX[i] = X.v_p[i];
		YYY[i] = Y.v_p[i];
		ZZZ[i] = 0;
	}
	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		//vdMul(jsize, XXX, YYY, ZZZ);
		//vsMul(jsize, XXX, YYY, ZZZ);
	}
	printf("Time taken (New Vector MKL): %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y: \n");
	cout.precision(2);
	cout << "[";
	for (int i = 0; i < 10; ++i) { cout << ZZZ[i] << " "; }
	cout << "]";
	cout << endl << endl;
	delete[] XXX;
	delete[] YYY;
	delete[] ZZZ;
	/************************************************************************************/
	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		//vdMul(jsize, X.v_p, Y.v_p, temp.v_p);
		//vsMul(jsize, X.v_p, Y.v_p, temp.v_p);
	}
	printf("Time taken (Matrix MKL): %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y: \n");
	temp.output(0, 9);
	/************************************************************************************/
	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		for (int j = 0; j < jsize; ++j) {
			temp.v_p[j] = X.v_p[j] * Y.v_p[j];
		}
	}
	printf("Time taken (Matrix expansion): %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y: \n");
	temp.output(0, 9);
	/************************************************************************************/
	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		temp = X * Y;
		++icount;
	}
	printf("Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("Z = X .* Y (Matrix direct): \n");
	temp.output(0, 9);

}


void test_vectormultip(Matrix2<myTYPE> &X, Matrix2<myTYPE> &Y) {
	Matrix2<myTYPE> temp(X.NI, X.NJ);
	clock_t tStart;
	int icount = 0;
	const int COUNT = 10000;

	tStart = clock();
	for (int i = 0; i < COUNT; ++i) {
		temp = X*Y;
		++icount;
	}
	printf("Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
	printf("X .* Y = \n");
	temp.output();
}


void test_coordinate(Coordinate &X, Coordinate &Y) {
	myTYPE *a[3], t12[3][3], r12[3];
	Matrix2<myTYPE> tt(3, 3);
	X.Transfer(t12, r12, X, Y);
	a[0] = t12[0];
	a[1] = t12[1];
	a[2] = t12[2];
	tt = Arr2Mat(a, 3, 3);
	matrix_output_on_screen(tt, 4);
}