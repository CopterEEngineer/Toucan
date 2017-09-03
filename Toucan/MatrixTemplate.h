#pragma once
#ifndef MatrixTemplate_h
#define MatrixTemplate_h

//define matrix template and simple operator
#ifndef PI
#define PI 3.141592653589793238462643383279502884197169
#endif // !PI
#ifndef Sign
#define Sign(a) ((a) >= 0 ? 1 : -1)
#endif // !sign
#ifndef Max
#define Max(a,b) ((a)>(b)?(a):(b))
#endif // !max
#ifndef Min
#define Min(a,b) ((a)<(b)?(a):(b))
#endif // !min
#ifndef Abs
#define Abs(a) ((a)>=0?(a):(-(a)))
#endif // !abs

#define INI 0

#ifdef _DEBUG
#define BOUNDS_CHECK
#endif // bound check
#ifdef USE_MKL
#define _USE_MKL
//#undef _USE_MKL
#endif // USE_MKL

#define USE_DOUBLE
typedef double myTYPE; //float //


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "mkl.h"
#include "mkl_vml.h"
#include <omp.h>


using std::cout;
using std::endl;
using std::cerr;
using std::exception;
using std::string;
using std::ifstream;
using std::ofstream;


template<class _Ty>
_Ty _Atan2(const _Ty &m, const _Ty &n)
{
	_Ty aoa = 0;
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


//define matrix template
template<class Type> class Matrix1
{
public:
	Type *v_p;
	int I0, I1, Nv;

	Matrix1(int i) {
		Nv = i;
		I0 = 0;
		I1 = Nv - 1;
		v_p = new Type[Nv];
#if(INI)
		for (int i = I0; i <= I1; ++i) { v_p[i] = NAN; }
#else
		for (int i = I0; i <= I1; ++i) { v_p[i] = Type(0); }
#endif
	}


	Matrix1(double v, const Matrix1<Type> &x);


	Matrix1(int v, const Matrix1<Type> &x);


	Matrix1() {
		I0 = 0;
		I1 = 0;
		Nv = 0;
		v_p = new Type[Nv];
	}


	Matrix1(const Matrix1<Type> &A) {
		//cout << "Copy constructor." << endl;
		Nv = A.Nv;
		I0 = A.I0;
		I1 = A.I1;
		v_p = new Type[Nv];
		*v_p = *A.v_p;
		for (int i = I0; i <= I1; ++i) { v_p[i] = A.v_p[i]; }
	}


	~Matrix1() {
		I0 = 0;
		I1 = 0;
		Nv = 0;
		if (v_p != NULL) delete[] v_p;
		//cout << "[] is Deleted." << endl;
	}


	inline void check_bound(int i) {
#ifdef BOUNDS_CHECK
		try {
			if ((i < I0) || (i > I1)) throw i;
		}
		catch (int i) {
			printf("Bound error: i = %d is out of [%d, %d]. \n", i, I0, I1);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_bound_v_sorted(int x) {
#ifdef BOUNDS_CHECK
		try {
			if ((x < v_p[0]) || (x > v_[I1])) throw 99;
		}
		catch (int i) {
			printf("Bound error: input value %d is out of [%f, %f]. \n", x, v_p[0], v_p[Nv - 1]);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_bound_v_sorted(float x) {
#ifdef BOUNDS_CHECK
		try {
			if ((x < v_p[0]) || (x > v_[I1])) throw 99;
		}
		catch (int i) {
			printf("Bound error: input value %f is out of [%f, %f]. \n", x, v_p[0], v_p[Nv - 1]);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_bound_v_sorted(double x) {
#ifdef BOUNDS_CHECK
		try {
			if ((x < v_p[0]) || (x > v_[I1])) throw 99;
		}
		catch (int i) {
			printf("Bound error: input value %f is out of [%f, %f]. \n", x, v_p[0], v_p[Nv - 1]);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_bound_v_sorted(const Matrix1<float> &A) {
#ifdef BOUNDS_CHECK
		try {
			if ((A.v_p[0] < v_p[0]) || (A.v_p[Nv - 1] > v_p[I1])) throw 99;
		}
		catch (int i) {
			printf("Bound error: input Matrix1 value bound [%f, %f] is out of [%f, %f]. \n", A.v_p[0], A.v_p[Nv - 1], v_p[0], v_p[Nv - 1]);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_bound_v_sorted(const Matrix1<double> &A) {
#ifdef BOUNDS_CHECK
		try {
			if ((A.v_p[0] < v_p[0]) || (A.v_p[Nv - 1] > v_p[I1])) throw 99;
		}
		catch (int i) {
			printf("Bound error: input Matrix1 value bound [%f, %f] is out of [%f, %f]. \n", A.v_p[0], A.v_p[Nv - 1], v_p[0], v_p[Nv - 1]);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_bound_v_sorted(const Matrix1<int> &A) {
#ifdef BOUNDS_CHECK
		try {
			if ((A.v_p[0] < v_p[0]) || (A.v_p[Nv - 1] > v_p[I1])) throw 99;
		}
		catch (int i) {
			printf("Bound error: input Matrix1 value bound [%d, %d] is out of [%f, %f]. \n", A.v_p[0], A.v_p[Nv - 1], v_p[0], v_p[Nv - 1]);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_shape(const Matrix1<Type> &A) {
#ifdef BOUNDS_CHECK
		try {
			if (this->Nv != A.Nv) throw 99;
		}
		catch (int i) {
			printf("Shape error %d: shape = %d is assigned to %d. \n", i, A.Nv, this->Nv);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif
	}


	inline void check_shape(const Matrix1<Type> &A, const Matrix1<Type> &B) {
#ifdef BOUNDS_CHECK
		try {
			if (A.Nv != B.Nv) throw 99;
		}
		catch (int i) {
			printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
#endif // BOUNDS_CHECK
	}


	inline void allocate(int i) {
		if (!Nv) {
			delete[] v_p;
			Nv = i;
			I0 = 0;
			I1 = Nv - 1;
			v_p = new Type[Nv];
		}
		else {
			cout << "Matrix has been allocated." << endl;
			//system("pause");
			cout << "Go" << endl;
		}
#if(INI)
		for (int i = I0; i <= I1; ++i) v_p[i] = NAN;
#else
		for (int i = I0; i <= I1; ++i) v_p[i] = Type(0);
#endif
	}


	inline void deallocate() {
		I0 = 0;
		I1 = 0;
		Nv = 0;
		delete[] v_p;
		v_p = new Type[0];
	}


	inline void setvalue(Type x) {
		for (int i = I0; i <= I1; ++i) { v_p[i] = x; }
	}


	inline Type sum() {
		Type s = 0.0;
		for (int i = I0; i <= I1; ++i) { s += v_p[i]; }
		return s;
	}


	inline Type sum(int i0, int i1) {
		Type s = 0.0;
		try { if (i0 > i1) throw 99; }
		catch (int i) {
			cerr << "Order error: " << i0 << " > " << i1 << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		check_bound(i0);
		check_bound(i1);
		for (int i = i0; i <= i1; ++i) { s += v_p[i - i0]; }
		return s;
	}


	inline Matrix1<Type> msqrt(void) {
		Matrix1<Type> temp(Nv);
		for (int i = 0; i < Nv; ++i) {
#if _DEBUG
			if (v_p[i] < 0) {
				cout << "Negitive number got in sqrt." << endl;
				system("pause");
			}
#endif // _DEBUG
			temp.v_p[i] = sqrt(v_p[i]);
		}
		return temp;
	}


	inline Matrix1<int> step(int i0, int i1) {
		try { if (i0 > i1) throw 99; }
		catch (int i) {
			cerr << "Order error: " << i0 << " > " << i1 << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		Matrix1<int> temp(i1 - i0 + 1);
		for (int i = 0; i <= i1 - i0; ++i) {
			temp.v_p[i] = i + i0;
		}
		return temp;
	}


	inline Type findmax(void) {
		Type maxvalue = *(v_p + Nv - 1);

		if (Nv > 1) {
			for (int i = Nv - 2; i >= 0; --i) {
				if (*(v_p + i) > maxvalue) maxvalue = *(v_p + i);
			}
		}
		return maxvalue;
	}


	inline Type findmin(void) {
		Type minvalue = *(v_p + Nv - 1);

		if (Nv > 1) {
			for (int i = Nv - 2; i >= 0; --i) {
				if (*(v_p + i) < minvalue) minvalue = *(v_p + i);
			}
		}
		return minvalue;
	}


	inline void findmaxmin2(Type &maxvalue, Type &minvalue) {
		if (Nv == 1) {
			maxvalue = *v_p;
			minvalue = *v_p;
		}
		else if (Nv == 2) {
			maxvalue = Max(*(v_p + Nv - 1), *(v_p + Nv - 2));
			minvalue = Min(*(v_p + Nv - 1), *(v_p + Nv - 2));
		}
		else {
			int istart = Nv - 2;
			if (Nv % 2 == 0) {
				maxvalue = Max(*(v_p + Nv - 1), *(v_p + Nv - 2));
				minvalue = Min(*(v_p + Nv - 1), *(v_p + Nv - 2));
				istart = Nv - 3;
			}
			else {
				maxvalue = *(v_p + Nv - 1);
				minvalue = *(v_p + Nv - 1);
			}

			for (int i = istart; i >= 0; i -= 2) {
				if (*(v_p + i) > *(v_p + i - 1)) {
					if (*(v_p + i) > maxvalue) maxvalue = *(v_p + i);
					if (*(v_p + i - 1) < minvalue) minvalue = *(v_p + i - 1);
				}
				else {
					if (*(v_p + i) < minvalue) minvalue = *(v_p + i);
					if (*(v_p + i - 1) > maxvalue) maxvalue = *(v_p + i - 1);
				}
			}
		}

	}


	inline void findmaxmin(Type &maxvalue, Type &minvalue) {
		maxvalue = *(v_p + Nv - 1);
		minvalue = *(v_p + Nv - 1);

		if (Nv > 1) {
			for (int i = Nv - 2; i >= 0; --i) {
				if (*(v_p + i) < minvalue) minvalue = *(v_p + i);
				if (*(v_p + i) > maxvalue) maxvalue = *(v_p + i);
			}
		}
	}


	inline Type interplinear(Matrix1<Type> &X, Type xs) {
		Type v = 0;
		Type drhs, dlhs, maxvalue, minvalue;
		int  icount = 0;

#ifdef BOUNDS_CHECK
		check_shape(X);
		X.findmaxmin(maxvalue, minvalue);
		if ((xs > maxvalue) || (xs < minvalue)) {
			cout << "Bound error: " << xs << " is out of [" << minvalue << ", " << maxvalue << "]." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		for (int i = Nv - 1; i >= 1; --i) {
			icount++;
			if (X.v_p[i] == X.v_p[i - 1]) {
				cout << "Warning: linear may be unsteady near " << X.v_p[i] << endl;
				system("pause");
			}
			drhs = X.v_p[i] - xs;
			dlhs = X.v_p[i - 1] - xs;
			if (drhs * dlhs <= 0) {
				v = (v_p[i] - v_p[i - 1]) / (X.v_p[i] - X.v_p[i - 1]) * (xs - X.v_p[i - 1]) + v_p[i - 1];
				break;
			}
		}
		if (drhs * dlhs > 0) {
			X.findmaxmin(maxvalue, minvalue);
			cout << "Bound error: " << xs << " is out of [" << minvalue << ", " << maxvalue << "]." << endl;
			system("pause");
		}
		//printf("%d\n", icount);
		return v;
	}


	inline Type interplinear_fast(Matrix1<Type> &X, Type xs) {
		Type v = 0;
		Type drhs, dlhs, maxvalue, minvalue;
		int li = 0;
		int ri = X.Nv - 1;
		int m, index;
		int icount = 0;

#ifdef BOUNDS_CHECK
		check_shape(X);
		X.findmaxmin(maxvalue, minvalue);
		if ((xs > maxvalue) || (xs < minvalue)) {
			cout << "Bound error: " << xs << " is out of [" << minvalue << ", " << maxvalue << "]." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		for (int i = Nv - 1; i >= 0; --i) {
			icount++;
			if (li > ri) {
				v = (v_p[li] - v_p[ri]) / (X.v_p[li] - X.v_p[ri]) * (xs - X.v_p[ri]) + v_p[ri];
				break;
			}
			else {
				m = (li + ri) / 2;
				if (*(X.v_p + m) < xs) li = m + 1;
				else if (*(X.v_p + m) > xs) ri = m - 1;
				else {
					v = v_p[m];
					break;
				}
			}
		}
		if ((ri < 0) || (li > X.Nv)) {
			X.findmaxmin(maxvalue, minvalue);
			cout << "Bound error: " << xs << " is out of [" << minvalue << ", " << maxvalue << "]." << endl;
			system("pause");
		}
		//printf("%d\n", icount);

		return v;
	}


	inline Matrix1<Type> interplinear(Matrix1<Type> &X, Matrix1<Type> &XN) {
		Matrix1<Type> temp(XN.Nv);
		Type x_max, x_min, xn_max, xn_min;
		Type xs, drhs, dlhs;

#ifdef BOUNDS_CHECK
		check_shape(X);
		X.findmaxmin(x_max, x_min);
		XN.findmaxmin(xn_max, xn_min);
		if ((xn_max > x_max) || (xn_min < x_min)) {
			cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		for (int j = XN.Nv - 1; j >= 0; --j) {
			//temp.v_p[i] = this->interplinear(X, XN.v_p[i]);
			xs = XN.v_p[j];
			for (int i = Nv - 1; i >= 1; --i) {
				if (X.v_p[i] == X.v_p[i - 1]) {
					cout << "Warning: linear may be unsteady near " << X.v_p[i] << endl;
					system("pause");
				}
				drhs = X.v_p[i] - xs;
				dlhs = X.v_p[i - 1] - xs;
				if (drhs * dlhs <= 0) {
					temp.v_p[j] = (v_p[i] - v_p[i - 1]) / (X.v_p[i] - X.v_p[i - 1]) * (xs - X.v_p[i - 1]) + v_p[i - 1];
					break;
				}
			}

			if (drhs * dlhs > 0) {
				X.findmaxmin(x_max, x_min);
				XN.findmaxmin(xn_max, xn_min);
				cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
				system("pause");
			}

		}
		return temp;
	}


	inline Matrix1<Type> interplinear_fast(Matrix1<Type> &X, Matrix1<Type> &XN) {
		Matrix1<Type> temp(XN.Nv);
		Type x_max, x_min, xn_max, xn_min, xs;
		int li, ri, m, jsize;

#ifdef BOUNDS_CHECK
		check_shape(X);
		X.findmaxmin(x_max, x_min);
		XN.findmaxmin(xn_max, xn_min);
		if ((xn_max > x_max) || (xn_min < x_min)) {
			cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		jsize = XN.Nv;
		for (int j = jsize - 1; j >= 0; --j) {
			xs = XN.v_p[j];
			li = 0;
			ri = X.Nv - 1;
			for (int i = Nv - 1; i >= 0; --i) {
				if (li <= ri) {
					m = (li + ri) / 2;
					if (*(X.v_p + m) < xs) li = m + 1;
					else if (*(X.v_p + m) > xs) ri = m - 1;
					else {
						temp.v_p[j] = v_p[m];
						break;
					}
				}
				else {
					temp.v_p[j] = (v_p[li] - v_p[ri]) / (X.v_p[li] - X.v_p[ri]) * (xs - X.v_p[ri]) + v_p[ri];
					break;
				}
			}
			if ((ri < 0) || (li > X.Nv)) {
				X.findmaxmin(x_max, x_min);
				XN.findmaxmin(xn_max, xn_min);
				cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
				system("pause");
			}
		}
		return temp;
	}

	inline void input(string filename) {
		int m;
		ifstream InFile(filename);
		InFile >> m;
#ifdef BOUNDS_CHECK
		try {
			if (m != Nv) throw 99;
		}
		catch (int i) {
			printf("Shape error: input %d is not consistent with %d. \n", m, Nv);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
		if (InFile) {
			for (int i = 0; i < m; ++i) {
				InFile >> this->v_p[i];
			}
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}
	}

	inline void output(int presc = 2) {
		cout.precision(presc);
		cout << "[";
		for (int i = I0; i <= I1; ++i) { cout << v_p[i] << " "; }
		cout << "]";
		cout << endl << endl;
	}


	inline void output(int i0, int i1, int presc = 2) {
		cout.precision(presc);
		cout << "[";
		for (int i = i0; i <= i1; ++i) { cout << v_p[i] << " "; }
		cout << "]";
		cout << endl << endl;
	}


	inline void output(string filename, int presc = 2) {
		ofstream OutFile(filename);
		if (OutFile) {
			OutFile.precision(presc);
			for (int i = I0; i <= I1; ++i) {
				OutFile << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
				OutFile << v_p[i] << '\t';
			}
			OutFile << endl;
			OutFile.close();
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}


		OutFile.close();
	}


	inline void outputs(string filename, int presc = 2) {
		ofstream OutFile(filename, std::ios::app);
		if (OutFile) {
			OutFile.precision(presc);
			for (int i = I0; i <= I1; ++i) {
				OutFile << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
				OutFile << v_p[i] << '\t';
			}
			OutFile << endl;
			OutFile.close();
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}


		OutFile.close();
	}

	inline void output(ofstream &outfile) {
		for (int i = I0; i <= I1; ++i) { outfile << v_p[i] << "\t"; }
		outfile.close();
	}


	inline void output(int i0, int i1, string filename) {
		check_bound(i0);
		check_bound(i1);
		ofstream OutFile(filename);
		OutFile.precision(10);
		for (int i = i0; i <= i1; ++i) { OutFile << v_p[i - I0] << "\t"; }
		OutFile.close();
	}


	inline void output(int i0, int i1, ofstream &outfile) {
		check_bound(i0);
		check_bound(i1);
		for (int i = i0; i <= i1; ++i) { outfile << v_p[i - I0] << "\t"; }
	}


	inline Type dot(const Matrix1<Type> &Y) {
		check_shape(Y);
		Type temp = 0;
		int isize = Y.Nv;
		for (int i = 0; i < isize; ++i) {
			temp += (*this)(i) * Y.v_p[i];
		}
		return temp;
	}


	inline Type dotP(const Matrix1<Type> &Y) {
		check_shape(Y);
		Type temp = 0;
		Type *X_p, *Y_p;
		int isize = Y.Nv;
		X_p = (this->v_p);
		Y_p = Y.v_p;
		for (int i = 0; i < isize; ++i) {
			temp += (*(X_p++))*(*(Y_p++));
		}
		return temp;
	}




	// operator overloading
	inline Matrix1<Type> & operator=(const Matrix1<Type> &A) {
		if (this->Nv != A.Nv) {
			if (this->Nv != 0) {
				cout << "Warning: inconsistency of matrix size." << endl;
				//system("pause");
			}

			this->deallocate();
			this->allocate(A.Nv);
		}

		int isize = A.Nv;
		for (int i = isize - 1; i >= 0; --i) {
			this->v_p[i] = A.v_p[i];
		}

		return (*this);
	}


	inline Type & operator()(int i) {
		check_bound(i);
		return *(v_p + (i - I0));
	}


	inline Matrix1<Type> operator()(const Matrix1<int> &A);


	inline Matrix1<Type> & operator+=(const Matrix1<Type> &A) {
		check_shape(A);
		int isize = Nv;
		for (int i = isize - 1; i >= 0; --i) {
			this->v_p[i] += A.v_p[i];
		}
		return (*this);
	}


	inline Matrix1<Type> & operator-=(const Matrix1<Type> &A) {
		check_shape(A);
		int isize = Nv;
		for (int i = isize - 1; i >= 0; --i) {
			this->v_p[i] -= A.v_p[i];
		}
		return (*this);
	}


	inline Matrix1<Type> & operator*=(const Matrix1<Type> &A) {
		check_shape(A);
		int isize = Nv;
		for (int i = isize - 1; i >= 0; --i) {
			this->v_p[i] *= A.v_p[i];
		}
		return (*this);
	}





	inline Matrix1<Type> & operator/=(const Matrix1<Type> &A) {
		check_shape(A);
		int isize = Nv;
		try {
			for (int i = isize - 1; i >= 0; --i) {
				this->v_p[i] /= A.v_p[i];
			}
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return (*this);
	}


	inline Matrix1<Type> operator+(const Matrix1<Type> &A);


	inline Matrix1<Type> operator-(const Matrix1<Type> &A);


	inline Matrix1<Type> operator*(const Matrix1<Type> &A);


	inline Matrix1<Type> operator/(const Matrix1<Type> &A);


	inline Matrix1<Type> operator+(const Type x) {
		Matrix1<Type> temp(Nv);
		for (int i = I0; i <= I1; ++i) { temp.v_p[i] = v_p[i] + x; }
		return temp;
	}


	inline Matrix1<Type> operator-(const Type x) {
		Matrix1<Type> temp(Nv);
		for (int i = I0; i <= I1; ++i) { temp.v_p[i] = v_p[i] - x; }
		return temp;
	}


	inline Matrix1<Type> operator*(const Type x) {
		Matrix1<Type> temp(Nv);
		for (int i = I0; i <= I1; ++i) { temp.v_p[i] = v_p[i] * x; }
		return temp;
	}


	inline Matrix1<Type> operator/(const Type x) {
		Matrix1<Type> temp(Nv);
		try {
			for (int i = I0; i <= I1; ++i) { temp.v_p[i] = v_p[i] / x; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return temp;
	}

};


template<class Type> class Matrix2
{
public:
	int I0, I1, J0, J1, NI, NJ, Nv;
	Type *v_p;

	Matrix2(int i, int j) {
		NI = i;
		NJ = j;
		I0 = 0;
		J0 = 0;
		I1 = NI - 1;
		J1 = NJ - 1;
		Nv = NI * NJ;
		v_p = new Type[Nv];
#if(INI)
		for (int i = I0; i < Nv; ++i) { v_p[i] = NAN; }
#else
		for (int i = I0; i < Nv; ++i) { v_p[i] = Type(0); }
#endif
	}


	Matrix2(double v, const Matrix2<Type> &x);


	Matrix2(int v, const Matrix2<Type> &x);


	Matrix2() {
		NI = 0;
		NJ = 0;
		I0 = 0;
		J0 = 0;
		I1 = 0;
		J1 = 0;
		Nv = 0;
		v_p = new Type[Nv];
	}


	Matrix2(const Matrix2<Type> &A) {
		//cout << "Copy constructor." << endl;
		NI = A.NI;
		NJ = A.NJ;
		I0 = A.I0;
		I1 = A.I1;
		J0 = A.J0;
		J1 = A.J1;
		Nv = A.Nv;
		v_p = new Type[Nv];
		*v_p = *A.v_p;
		for (int i = I0; i < Nv; ++i) { v_p[i] = A.v_p[i]; }
	}


	~Matrix2() {
		NI = 0;
		NJ = 0;
		I0 = 0;
		J0 = 0;
		I1 = 0;
		J1 = 0;
		Nv = 0;
		if (v_p != NULL) delete[] v_p;
		v_p = NULL;
		//cout << "[] is deleted." << endl;
	}


	inline void check_bound(int i, int j) {
#ifdef BOUNDS_CHECK
		try {
			if ((i < I0) || (i > I1) || (j < J0) || (j > J1)) throw 99;
		}
		catch (int id) {
			printf("Bound error: i, j = %d, %d is out of [%d,%d, %d,%d]. \n", i, j, I0, I1, J0, J1);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
	}


	inline void check_shape(const Matrix2<Type> &A) {
#ifdef BOUNDS_CHECK
		try {
			if (this->Nv != A.Nv) throw 99;
			else if (this->NI != A.NI) throw 99;
			else if (this->NJ != A.NJ) throw 99;
		}
		catch (int i) {
			printf("Shape error: %d %d is not consisten with %d %d. \n", NI, NJ, A.NI, A.NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif	
	}


	inline void check_shape(const Matrix2<Type> &A, const Matrix2<Type> &B) {
#ifdef BOUNDS_CHECK
		try {
			if (B.Nv != A.Nv) throw 99;
			else if (B.NI != A.NI) throw 99;
			else if (B.NJ != A.NJ) throw 99;
		}
		catch (int i) {
			printf("Shape error: %d %d is not consisten with %d %d. \n", NI, NJ, A.NI, A.NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
	}


	inline void allocate(int i, int j) {
		if (!Nv) {
			delete[] v_p;
			NI = i;
			NJ = j;
			I0 = 0;
			J0 = 0;
			I1 = NI - 1;
			J1 = NJ - 1;
			Nv = NI * NJ;
			v_p = new Type[Nv];
		}
		else {
			//cout << "Matrix has been allocated." << endl;
			//system("pause");
			//cout << "Go" << endl;
		}
#if(INI)
		for (int i = I0; i < Nv; ++i) { v_p[i] = NAN; }
#else
		for (int i = I0; i < Nv; ++i) { v_p[i] = Type(0); }
#endif
	}


	inline void deallocate() {
		I0 = 0;
		I1 = 0;
		I1 = 0;
		J1 = 0;
		NI = 0;
		NJ = 0;
		Nv = 0;
		if (v_p != NULL) delete[] v_p;
		v_p = new Type[0];
	}


	inline void setvalue(Type x) {
		int isize = Nv;
		for (int i = isize - 1; i >= 0; --i) { v_p[i] = x; }
	}


	inline Type sum(void) {
		Type s = 0;
		int isize = Nv;
		for (int i = isize - 1; i >= 0; --i) { s += v_p[i]; }
		return s;
	}


	inline Matrix1<Type> sum(int index0, int index1, int axis);


	inline Matrix2<Type> transpose(void) {
		int isize, jsize;
		isize = NI;
		jsize = NJ;
		Matrix2<Type> temp(jsize, isize);
		for (int j = jsize - 1; j >= 0; --j) {
			for (int i = isize - 1; i >= 0; --i) {
				temp(j, i) = (*this)(i, j);
			}
		}
		return temp;
	}


	inline Matrix2<Type> transpose(const Matrix2<Type> &A) {
		int isize, jsize;
		isize = A.NI;
		jsize = A.NJ;
		Matrix2<Type> temp(jsize, isize);
		for (int j = jsize - 1; j >= 0; --j) {
			for (int i = isize - 1; i >= 0; --i) {
				temp(j, i) = A.v_p(i, j);
			}
		}
		return temp;
	}


	inline Matrix2<Type> msqrt(void) {
		Matrix2<Type> temp(NI, NJ);
		for (int i = Nv - 1; i >= 0; --i) {
#if _DEBUG
			if (v_p[i] < 0) {
				cout << "Negitive number got in sqrt." << endl;
				system("pause");
			}
#endif // _DEBUG
			temp.v_p[i] = sqrt(v_p[i]);
		}
		return temp;
	}


	inline Type findmax(void) {
		Type maxvalue = *(v_p + Nv - 1);

		if (Nv > 1) {
			for (int i = Nv - 2; i >= 0; --i) {
				if (*(v_p + i) > maxvalue) maxvalue = *(v_p + i);
			}
		}
		return maxvalue;
	}


	inline Matrix1<Type> findmax(int axis) {
		Matrix1<Type> temp;
		Matrix1<int> id;

		if (axis == 0) {
			// row max
			id = id.step(0, NJ - 1);
			temp.allocate(NI);
			for (int i = NI - 1; i >= 0; --i) {
				temp(i) = (*this)(i, id).findmax();
			}
		}
		else if (axis == 1) {
			// column max
			id = id.step(0, NI - 1);
			temp.allocate(NJ);
			for (int j = NJ - 1; j >= 0; --j) {
				temp(j) = (*this)(id, j).findmax();
			}
		}
		else {
			printf("Error axis: %d", axis);
			system("pause");
		}
		return temp;
	}


	inline Type findmin(void) {
		Type minvalue = *(v_p + Nv - 1);

		if (Nv > 1) {
			for (int i = Nv - 2; i >= 0; --i) {
				if (*(v_p + i) < minvalue) minvalue = *(v_p + i);
			}
		}
		return minvalue;
	}


	inline Matrix1<Type> findmin(int axis) {
		Matrix1<Type> temp;
		Matrix1<int> id;

		if (axis == 0) {
			// row max
			id = id.step(0, NJ - 1);
			temp.allocate(NI);
			for (int i = NI - 1; i >= 0; --i) {
				temp(i) = (*this)(i, id).findmin();
			}
		}
		else if (axis == 1) {
			// column max
			id = id.step(0, NI - 1);
			temp.allocate(NJ);
			for (int j = NJ - 1; j >= 0; --j) {
				temp(j) = (*this)(id, j).findmin();
			}
		}
		else {
			printf("Error axis: %d", axis);
			system("pause");
		}
		return temp;
	}


	inline void findmaxmin(Type &maxvalue, Type &minvalue) {
		maxvalue = *(v_p + Nv - 1);
		minvalue = *(v_p + Nv - 1);

		if (Nv > 1) {
			for (int i = Nv - 2; i >= 0; --i) {
				if (*(v_p + i) < minvalue) minvalue = *(v_p + i);
				if (*(v_p + i) > maxvalue) maxvalue = *(v_p + i);
			}
		}
	}


	inline Type interplinear(Matrix1<Type> &X, Matrix1<Type> &Y, Type xs, Type ys) {
		//:param X: row
		//:param Y: col
		Type v, x1, x0, y1, y0;
		Type f00, f01, f10, f11, drhs_x, dlhs_x, drhs_y, dlhs_y;
		Type x_max, x_min, y_max, y_min;
		int ix, iy;
#ifdef BOUNDS_CHECK
		if ((NI != X.Nv) || (NJ != Y.Nv)) {
			printf("Shape error: Z(%d, %d) is inconsistent with X %d, Y %d. \n", NI, NJ, X.Nv, Y.Nv);
			system("pause");
		}

		X.findmaxmin(x_max, x_min);
		Y.findmaxmin(y_max, y_min);

		if ((xs > x_max) || (xs < x_min)) {
			cout << "Bound error: " << xs << " is out of [" << x_min << ", " << x_max << "]." << endl;
			system("pause");
		}
		if ((ys > y_max) || (ys < y_min)) {
			cout << "Bound error: " << ys << " is out of [" << y_min << ", " << y_max << "]." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		for (int i = NI - 1; i >= 1; --i) {
			if (X.v_p[i] == X.v_p[i - 1]) {
				cout << "Warning: linear may be unsteady near " << X.v_p[i] << endl;
				system("pause");
			}
			drhs_x = X.v_p[i] - xs;
			dlhs_x = X.v_p[i - 1] - xs;
			if (drhs_x * dlhs_x <= 0) {
				x1 = X.v_p[i];
				x0 = X.v_p[i - 1];
				ix = i;
				break;
			}
		}


		for (int j = NJ - 1; j >= 1; --j) {
			if (Y.v_p[j] == Y.v_p[j - 1]) {
				cout << "Warning: linear may be unsteady near " << Y.v_p[j] << endl;
				system("pause");
			}
			drhs_y = Y.v_p[j] - ys;
			dlhs_y = Y.v_p[j - 1] - ys;
			if (drhs_y * dlhs_y <= 0) {
				y1 = Y.v_p[j];
				y0 = Y.v_p[j - 1];
				iy = j;
				break;
			}
		}

		if ((drhs_x * dlhs_x > 0) || (drhs_y * dlhs_y > 0)) {
			X.findmaxmin(x_max, x_min);
			Y.findmaxmin(y_max, y_min);
			cout << "Bound error: " << xs << " is out of [" << x_min << ", " << x_max << "]." << endl;
			cout << "Bound error: " << ys << " is out of [" << y_min << ", " << y_max << "]." << endl;
			system("pause");
		}

		f00 = (*this)(ix - 1, iy - 1);
		f01 = (*this)(ix - 1, iy);
		f10 = (*this)(ix, iy - 1);
		f11 = (*this)(ix, iy);

		v = ((x1 - xs) * f00 + (xs - x0) * f10) * (y1 - ys);
		v += (f01*(x1 - xs) + f11 * (xs - x0)) * (ys - y0);
		v /= (x1 - x0) * (y1 - y0);

		return v;
	}


	inline Type interplinear_fast(Matrix1<Type> &X, Matrix1<Type> &Y, Type xs, Type ys) {
		//:param X: row
		//:param Y: col
		Type v, x1, x0, y1, y0;
		Type f00, f01, f10, f11;
		Type x_max, x_min, y_max, y_min;
		int li, ri, m, index_r, index_c;
#ifdef BOUNDS_CHECK
		if ((NI != X.Nv) || (NJ != Y.Nv)) {
			printf("Shape error: Z(%d, %d) is inconsistent with X %d, Y %d. \n", NI, NJ, X.Nv, Y.Nv);
			system("pause");
		}

		X.findmaxmin(x_max, x_min);
		Y.findmaxmin(y_max, y_min);

		if ((xs > x_max) || (xs < x_min)) {
			cout << "Bound error: " << xs << " is out of [" << x_min << ", " << x_max << "]." << endl;
			system("pause");
		}
		if ((ys > y_max) || (ys < y_min)) {
			cout << "Bound error: " << ys << " is out of [" << y_min << ", " << y_max << "]." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		li = 0;
		ri = X.Nv - 1;
		index_r = X.Nv;
		for (int i = NI - 1; i >= 0; --i) {
			if (li > ri) {
				index_r = ri;
				break;
			}
			else {
				m = (li + ri) / 2;
				if (*(X.v_p + m) < xs) li = m + 1;
				else if (*(X.v_p + m) > xs) ri = m - 1;
				else {
					index_r = m;
					break;
				}
			}
		}
		li = 0;
		ri = Y.Nv - 1;
		index_c = Y.Nv;
		for (int j = NJ - 1; j >= 0; --j) {
			if (li > ri) {
				index_c = ri;
				break;
			}
			else {
				m = (li + ri) / 2;
				if (*(Y.v_p + m) < ys) li = m + 1;
				else if (*(Y.v_p + m) > ys) ri = m - 1;
				else {
					index_c = m;
					break;
				}
			}
		}

		if ((index_r > X.Nv - 1) || (index_c > Y.Nv - 1)) {
			X.findmaxmin(x_max, x_min);
			Y.findmaxmin(y_max, y_min);
			cout << "Bound error: " << xs << " is out of [" << x_min << ", " << x_max << "]." << endl;
			cout << "Bound error: " << ys << " is out of [" << y_min << ", " << y_max << "]." << endl;
			system("pause");
		}

		f00 = (*this)(index_r, index_c);
		f01 = (*this)(index_r, index_c + 1);
		f10 = (*this)(index_r + 1, index_c);
		f11 = (*this)(index_r + 1, index_c + 1);
		x0 = X.v_p[index_r];
		x1 = X.v_p[index_r + 1];
		y0 = Y.v_p[index_c];
		y1 = Y.v_p[index_c + 1];

		v = ((x1 - xs) * f00 + (xs - x0) * f10) * (y1 - ys);
		v += (f01*(x1 - xs) + f11 * (xs - x0)) * (ys - y0);
		v /= (x1 - x0) * (y1 - y0);

		return v;
	}


	inline Matrix2<Type> interplinear(Matrix1<Type> &X, Matrix1<Type> &Y, Matrix2<Type> &XN, Matrix2<Type> &YN) {
		Matrix2<Type> temp(XN.NI, YN.NJ);
		Type v, x1, x0, y1, y0;
		Type xs, ys, drhs_x, dlhs_x, drhs_y, dlhs_y;
		Type f00, f01, f10, f11, x_max, x_min, y_max, y_min, xn_max, xn_min, yn_max, yn_min;
		int ix, iy, ni, nj;
#ifdef BOUNDS_CHECK
		if ((NI != X.Nv) || (NJ != Y.Nv)) {
			printf("Shape error: Z(%d, %d) is inconsistent with X %d, Y %d. \n", NI, NJ, X.Nv, Y.Nv);
			system("pause");
		}

		check_shape(XN, YN);

		X.findmaxmin(x_max, x_min);
		Y.findmaxmin(y_max, y_min);
		XN.findmaxmin(xn_max, xn_min);
		YN.findmaxmin(yn_max, yn_min);

		if ((xn_max > x_max) || (xn_min < x_min)) {
			cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
			system("pause");
		}
		if ((yn_max > y_max) || (yn_min < y_min)) {
			cout << "Bound error: " << yn_max << ", " << yn_min << " is out of " << y_max << ", " << y_min << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		ni = NI;
		nj = NJ;
		for (int k = XN.Nv - 1; k >= 0; --k) {
			xs = XN.v_p[k];
			ys = YN.v_p[k];
			for (int i = ni - 1; i >= 1; --i) {
				if (X.v_p[i] == X.v_p[i - 1]) {
					cout << "Warning: linear may be unsteady near " << X.v_p[i] << endl;
					system("pause");
				}
				drhs_x = X.v_p[i] - xs;
				dlhs_x = X.v_p[i - 1] - xs;
				if (drhs_x * dlhs_x <= 0) {
					x1 = X.v_p[i];
					x0 = X.v_p[i - 1];
					ix = i;
					break;
				}
			}

			for (int j = nj - 1; j >= 1; --j) {
				if (Y.v_p[j] == Y.v_p[j - 1]) {
					cout << "Warning: linear may be unsteady near " << Y.v_p[j] << endl;
					system("pause");
				}
				drhs_y = Y.v_p[j] - ys;
				dlhs_y = Y.v_p[j - 1] - ys;
				if (drhs_y * dlhs_y <= 0) {
					y1 = Y.v_p[j];
					y0 = Y.v_p[j - 1];
					iy = j;
					break;
				}
			}

			if ((drhs_x * dlhs_x > 0) || (drhs_y * dlhs_y > 0)) {
				X.findmaxmin(x_max, x_min);
				Y.findmaxmin(y_max, y_min);
				XN.findmaxmin(xn_max, xn_min);
				YN.findmaxmin(yn_max, yn_min);
				cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
				cout << "Bound error: " << yn_max << ", " << yn_min << " is out of " << y_max << ", " << y_min << endl;
				system("pause");
			}

			f00 = (*this)(ix - 1, iy - 1);
			f01 = (*this)(ix - 1, iy);
			f10 = (*this)(ix, iy - 1);
			f11 = (*this)(ix, iy);

			v = ((x1 - xs) * f00 + (xs - x0) * f10) * (y1 - ys);
			v += (f01*(x1 - xs) + f11 * (xs - x0)) * (ys - y0);
			v /= (x1 - x0) * (y1 - y0);
			temp.v_p[k] = v;
		}


		return temp;
	}


	inline Matrix2<Type> interplinear_fast(Matrix1<Type> &X, Matrix1<Type> &Y, Matrix2<Type> &XN, Matrix2<Type> &YN) {
		Matrix2<Type> temp(XN.NI, YN.NJ);
		Type v, x1, x0, y1, y0;
		Type xs, ys;
		Type f00, f01, f10, f11, x_max, x_min, y_max, y_min, xn_max, xn_min, yn_max, yn_min;
		int li, ri, m, index_r, index_c, ni, nj;
#ifdef BOUNDS_CHECK
		if ((NI != X.Nv) || (NJ != Y.Nv)) {
			printf("Shape error: Z(%d, %d) is inconsistent with X %d, Y %d. \n", NI, NJ, X.Nv, Y.Nv);
			system("pause");
		}

		check_shape(XN, YN);

		X.findmaxmin(x_max, x_min);
		Y.findmaxmin(y_max, y_min);
		XN.findmaxmin(xn_max, xn_min);
		YN.findmaxmin(yn_max, yn_min);

		if ((xn_max > x_max) || (xn_min < x_min)) {
			cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
			system("pause");
		}
		if ((yn_max > y_max) || (yn_min < y_min)) {
			cout << "Bound error: " << yn_max << ", " << yn_min << " is out of " << y_max << ", " << y_min << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		ni = NI;
		nj = NJ;
		for (int k = XN.Nv - 1; k >= 0; --k) {
			xs = XN.v_p[k];
			ys = YN.v_p[k];

			li = 0;
			ri = X.Nv - 1;
			index_r = X.Nv;
			for (int i = ni - 1; i >= 0; --i) {
				if (li > ri) {
					index_r = ri;
					break;
				}
				else {
					m = (li + ri) / 2;
					if (*(X.v_p + m) < xs) li = m + 1;
					else if (*(X.v_p + m) > xs) ri = m - 1;
					else {
						index_r = m;
						break;
					}
				}
			}
			li = 0;
			ri = Y.Nv - 1;
			index_c = Y.Nv;
			for (int j = nj - 1; j >= 0; --j) {
				if (li > ri) {
					index_c = ri;
					break;
				}
				else {
					m = (li + ri) / 2;
					if (*(Y.v_p + m) < ys) li = m + 1;
					else if (*(Y.v_p + m) > ys) ri = m - 1;
					else {
						index_c = m;
						break;
					}
				}
			}

			if ((index_r > X.Nv - 1) || (index_c > Y.Nv - 1)) {
				X.findmaxmin(x_max, x_min);
				Y.findmaxmin(y_max, y_min);
				cout << "Bound error: " << xs << " is out of [" << x_min << ", " << x_max << "]." << endl;
				cout << "Bound error: " << ys << " is out of [" << y_min << ", " << y_max << "]." << endl;
				system("pause");
			}

			f00 = (*this)(index_r, index_c);
			f01 = (*this)(index_r, index_c + 1);
			f10 = (*this)(index_r + 1, index_c);
			f11 = (*this)(index_r + 1, index_c + 1);
			x0 = X.v_p[index_r];
			x1 = X.v_p[index_r + 1];
			y0 = Y.v_p[index_c];
			y1 = Y.v_p[index_c + 1];

			v = ((x1 - xs) * f00 + (xs - x0) * f10) * (y1 - ys);
			v += (f01*(x1 - xs) + f11 * (xs - x0)) * (ys - y0);
			v /= (x1 - x0) * (y1 - y0);
			temp.v_p[k] = v;
		}
		return temp;
	}


	inline Matrix1<Type> interplinear_fast(Matrix1<Type> &X, Matrix1<Type> &Y, Matrix1<Type> &XN, Matrix1<Type> &YN) {
		Matrix1<Type> temp(XN.Nv);
		Type v, x1, x0, y1, y0;
		Type xs, ys;
		Type f00, f01, f10, f11, x_max, x_min, y_max, y_min, xn_max, xn_min, yn_max, yn_min;
		int li, ri, m, index_r, index_c, ni, nj;
		Type *x_ptr = X.v_p;
		Type *y_ptr = Y.v_p;
#ifdef BOUNDS_CHECK
		if ((NI != X.Nv) || (NJ != Y.Nv)) {
			printf("Shape error: Z(%d, %d) is inconsistent with X %d, Y %d. \n", NI, NJ, X.Nv, Y.Nv);
			system("pause");
		}

		XN.check_shape(XN, YN);

		X.findmaxmin(x_max, x_min);
		Y.findmaxmin(y_max, y_min);
		XN.findmaxmin(xn_max, xn_min);
		YN.findmaxmin(yn_max, yn_min);

		if ((xn_max > x_max) || (xn_min < x_min)) {
			cout << "Bound error: " << xn_max << ", " << xn_min << " is out of " << x_max << ", " << x_min << endl;
			system("pause");
		}
		if ((yn_max > y_max) || (yn_min < y_min)) {
			cout << "Bound error: " << yn_max << ", " << yn_min << " is out of " << y_max << ", " << y_min << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		ni = NI;
		nj = NJ;
		for (int k = XN.Nv - 1; k >= 0; --k) {
			xs = XN.v_p[k];
			ys = YN.v_p[k];

			li = 0;
			//ri = X.Nv - 1;
			//index_r = X.Nv;
			ri = ni - 1;
			index_r = ni;
			for (int i = ni - 1; i >= 0; --i) {
				if (li > ri) {
					index_r = ri;
					break;
				}
				else {
					m = (li + ri) / 2;
					//if (*(X.v_p + m) < xs) li = m + 1;
					//else if (*(X.v_p + m) > xs) ri = m - 1;
					//else {
					//	index_r = m;
					//	break;
					//}
					if (*(x_ptr + m) < xs) { li = m + 1; }
					else if (*(x_ptr + m) > xs) { ri = m - 1; }
					else {
						index_r = m;
						break;
					}
				}
			}
			li = 0;
			/*ri = Y.Nv - 1;
			index_c = Y.Nv;*/
			ri = nj - 1;
			index_c = nj;
			for (int j = nj - 1; j >= 0; --j) {
				if (li > ri) {
					index_c = ri;
					break;
				}
				else {
					m = (li + ri) / 2;
					/*if (*(Y.v_p + m) < ys) li = m + 1;
					else if (*(Y.v_p + m) > ys) ri = m - 1;
					else {
						index_c = m;
						break;
					}*/
					if (*(y_ptr + m) < ys) { li = m + 1; }
					else if (*(y_ptr + m) > ys) { ri = m - 1; }
					else {
						index_c = m;
						break;
					}
				}
			}

			if ((index_r > ni - 1) || (index_c > nj - 1)) {
				X.findmaxmin(x_max, x_min);
				Y.findmaxmin(y_max, y_min);
				cout << "Bound error: " << xs << " is out of [" << x_min << ", " << x_max << "]." << endl;
				cout << "Bound error: " << ys << " is out of [" << y_min << ", " << y_max << "]." << endl;
				system("pause");
			}

			// deal with boundary
			if (index_c == nj - 1) { --index_c; }
			if (index_c == -1) { ++index_c; }
			if (index_r == ni - 1) { --index_r; }
			if (index_r == -1) { ++index_r; }

			f00 = (*this)(index_r, index_c);
			f01 = (*this)(index_r, index_c + 1);
			f10 = (*this)(index_r + 1, index_c);
			f11 = (*this)(index_r + 1, index_c + 1);
			x0 = *(x_ptr + index_r);// X.v_p[index_r];
			x1 = *(x_ptr + index_r + 1);// X.v_p[index_r + 1];
			y0 = *(y_ptr + index_c);// Y.v_p[index_c];
			y1 = *(y_ptr + index_c + 1);// Y.v_p[index_c + 1];

			v = ((x1 - xs) * f00 + (xs - x0) * f10) * (y1 - ys);
			v += (f01*(x1 - xs) + f11 * (xs - x0)) * (ys - y0);
			v /= (x1 - x0) * (y1 - y0);
			temp.v_p[k] = v;
		}
		return temp;
	}


	inline void output(int presc = 2) {
		int ni = NI;
		int nj = NJ;
		cout.precision(presc);
		cout << "[";
		for (int i = 0; i < ni; ++i) {
			for (int j = 0; j < nj; ++j) {
				cout << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
				cout << (*this)(i, j) << " ";
			}
			cout << endl;
		}
		cout << "]";
		cout << endl << endl;
	}


	inline void output(string filename, int presc = 2) {
		ofstream OutFile(filename, std::ios::out | std::ios::trunc);
		int ni = NI;
		int nj = NJ;
		if (OutFile) {
			OutFile.precision(presc);
			for (int i = 0; i < ni; ++i) {
				for (int j = 0; j < nj; ++j) {
					OutFile << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
					OutFile << (*this)(i, j) << "\t";
				}
				OutFile << endl;
			}
			OutFile.close();
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}
	}


	inline void output(ofstream &outfile) {
		int ni = NI;
		int nj = NJ;
		if (outfile) {
			for (int i = 0; i < ni; ++i) {
				for (int j = 0; j < nj; ++j) {
					outfile << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
					outfile << (*this)(i, j) << "\t";
				}
				outfile << endl;
			}
			outfile.close();
		}
		else {
			cout << "File open failed." << endl;
			system("pause");
		}
	}


	inline void output(int i0, int i1, int j0, int j1, string filename);


	inline void output(int i0, int i1, int j0, int j1, ofstream &outfile);


	inline void input(string filename) {
		int m, n;
		ifstream InFile(filename);
		InFile >> m;
		InFile >> n;
#ifdef BOUNDS_CHECK
		try {
			if (!((m == NI) && (n == NJ))) throw 99;
		}
		catch (int i) {
			printf("Shape error: input %d %d is not consistent with %d %d. \n", m, n, NI, NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
		if (InFile) {
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					InFile >> this->v_p[j*NI + i];
				}
			}
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}
	}


	inline void input(string filename, Matrix2<Type> &A) {
		int m, n;
		ifstream InFile(filename);
		InFile >> m;
		InFile >> n;
#ifdef BOUNDS_CHECK
		try {
			if (!((m == A.NI) && (n == A.NJ))) throw 99;
		}
		catch (int i) {
			printf("Shape error: input %d %d is not consistent with %d %d. \n", m, n, A.NI, A.NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
		if (InFile) {
			A.deallocate();
			A.allocate(m, n);
			for (int i = 0; i < m; ++i) {
				for (int j = 0; j < n; ++j) {
					InFile >> A(i, j);
				}
			}
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}

	}


	inline Matrix1<Type> matrixmultiply(const Matrix1<Type> &b) {
#ifdef BOUNDS_CHECK
		try {
			if (NJ != b.Nv) throw 99;
		}
		catch (int i) {
			printf("shape error: shape = (%d, %d) multiplies shape %d. \n", NI, NJ, b.Nv);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		Matrix1<Type> temp(NI);
		Type btemp;
		for (int j = NJ - 1; j >= 0; --j) {
			btemp = b.v_p[j];
			for (int i = NI - 1; i >= 0; --i) {
				temp(i) += (*this)(i, j) * btemp;
			}
		}
		return temp;
	}


	inline Matrix1<Type> matrixmultiplyP(const Matrix1<Type> &b) {
#ifdef BOUNDS_CHECK
		try {
			if (NJ != b.Nv) throw 99;
		}
		catch (int i) {
			printf("shape error: shape = (%d, %d) multiplies shape %d. \n", NI, NJ, b.Nv);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");


		}
#endif // BOUNDS_CHECK
		Matrix1<Type> temp(NI);
		Type *b_p = &b.v_p[0];
		Type *A_p = &(this->v_p[0]);
		int ni = NI;
		int nj = NJ;
		for (int j = 0; j < nj; ++j) {
			for (int i = 0; i < ni; ++i) {
				temp(i) += (*(A_p++)) * (*b_p);
			}
			++b_p;
		}
		return temp;
	}


	inline Matrix1<Type> matrixmultiplyP2(const Matrix1<Type> &b) {
#ifdef BOUNDS_CHECK
		try {
			if (NJ != b.Nv) throw 99;
		}
		catch (int i) {
			printf("shape error: shape = (%d, %d) multiplies shape %d. \n", NI, NJ, b.Nv);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		int ni = NI;
		int nj = NJ;
		Matrix1<Type> temp(ni);
		Matrix2<Type> AT(nj, ni);

		AT = this->transpose();
		Type *A_p = &(AT.v_p[0]);
		for (int i = 0; i < ni; ++i) {
			Type temptemp = 0;
			Type *b_p = &b.v_p[0];
			for (int j = 0; j < nj; ++j) {
				temptemp += (*A_p++) * (*b_p);
				++b_p;
			}

			temp(i) = temptemp;
		}
		return temp;
	}


	inline Matrix1<Type> matrixmultiplyPB(const Matrix1<Type> &b);


	inline Matrix2<Type> matrixmultiply(const Matrix2<Type> &B) {
#ifdef BOUNDS_CHECK
		try {
			if (NJ != B.NI) throw 99;
		}
		catch (int i) {
			printf("shape error: shape = (%d, %d) multiplies shape (%d, %d). \n", NI, NJ, B.NI, B.NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		int ni = NI;
		int nj = NJ;
		int bnj = B.NJ;
		Matrix2<Type> temp(ni, bnj);
		Type temptemp = 0;
		for (int j = bnj - 1; j >= 0; --j) {
			for (int i = ni - 1; i >= 0; --i) {
				temptemp = 0;
				for (int k = nj - 1; k >= 0; --k) {
					temptemp += (*this)(i, k) * B.v_p[j*bnj + k]; // B(k, j)
				}
				temp(i, j) = temptemp;
			}
		}
		return temp;
	}


	inline Matrix2<Type> matrixmultiplyT(const Matrix2<Type> &B) {
#ifdef BOUNDS_CHECK
		try {
			if (NJ != B.NI) throw 99;
		}
		catch (int i) {
			printf("shape error: shape = (%d, %d) multiplies shape (%d, %d). \n", NI, NJ, B.NI, B.NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		int ni = NI;
		int nj = NJ;
		int bnj = B.NJ;
		int bni = B.NI;
		Matrix2<Type> temp(ni, bnj), AT(nj, ni);
		Type temptemp = 0;
		AT = this->transpose();
		for (int j = bnj - 1; j >= 0; --j) {
			for (int i = ni - 1; i >= 0; --i) {
				temptemp = 0;
				for (int k = nj - 1; k >= 0; --k) {
					temptemp += AT(k, i) * B.v_p[j*bni + k];
				}
				temp(i, j) = temptemp;
			}
		}
		return temp;
	}


	inline Matrix2<Type> matrixmultiplyTP(const Matrix2<Type> &B) {
#ifdef BOUNDS_CHECK
		try {
			if (NJ != B.NI) throw 99;
		}
		catch (int i) {
			printf("shape error: shape = (%d, %d) multiplies shape (%d, %d). \n", NI, NJ, B.NI, B.NJ);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK
		int ni = NI;
		int nj = NJ;
		int bnj = B.NJ;
		int bni = B.NI;
		Type *A_p, *B_p;
		Matrix2<Type> temp(ni, bnj), AT(nj, ni);
		Type temptemp = 0;
		AT = this->transpose();
		for (int j = 0; j < bnj; ++j) {
			A_p = &AT.v_p[0];
			for (int i = 0; i < ni; ++i) {
				temptemp = 0;
				B_p = &B.v_p[j*bni];
				for (int k = 0; k < nj; ++k) {
					temptemp += (*(A_p++)) * (*(B_p++));
				}
				temp(i, j) = temptemp; //左边是地址递增或递减
			}
		}
		return temp;
	}



	//Matrix3<Type> reshape(int ni, int nj, int nk)
	//	{
	//		Matrix3<Type> temp(ni, nj, nk);
	//#ifdef BOUNDS_CHECK
	//		if (Nv != ni*nj*nk) {
	//			cout << "Shape error when reshape Matrix2 to Matrix3." << endl;
	//			system("pause");
	//		}
	//#endif // BOUNDS_CHECK
	//		for (int k = 0; k < nk; ++k) {
	//			for (int j = 0; j < ni*nj; ++j) {
	//				temp(j / nk, j % nk, k) = (*this)(j, k);
	//			}
	//		}
	//		return temp;
	//	}


	// operator overloading
	inline Type & operator()(int i, int j) {
		check_bound(i, j);
		return *(v_p + j*NI + i);
	}


	inline Matrix1<Type> operator()(int i, const Matrix1<int> &A);


	inline Matrix1<Type> operator()(const Matrix1<int> &A, int j);


	inline Matrix2<Type> operator()(const Matrix2<int> &A);


	inline Matrix2<Type> & operator=(const Matrix2<Type> &A) {
		//check_shape(A);	
		if (~((this->Nv == A.Nv) && (this->NI == A.NI))) {
			this->deallocate();
			this->allocate(A.NI, A.NJ);
		}
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] = A.v_p[i]; }
		return (*this);
	}


	inline Matrix2<Type> & operator+=(const Matrix2<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] += A.v_p[i]; }
		return (*this);
	}


	inline Matrix2<Type> & operator-=(const Matrix2<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] -= A.v_p[i]; }
		return (*this);
	}


	inline Matrix2<Type> & operator*=(const Matrix2<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] *= A.v_p[i]; }
		return (*this);
	}


	inline Matrix2<Type> & operator/=(const Matrix2<Type> &A) {
		check_shape(A);
		try {
			for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] /= A.v_p[i]; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return (*this);
	}


	inline Matrix2<Type> operator+(const Matrix2<Type> &A);


	inline Matrix2<Type> operator-(const Matrix2<Type> &A);


	inline Matrix2<Type> operator*(const Matrix2<Type> &A);


	inline Matrix2<Type> operator/(const Matrix2<Type> &A);


	inline Matrix2<Type> operator+(const Type x) {
		Matrix2<Type> temp(NI, NJ);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] + x; }
		return temp;
	}


	inline Matrix2<Type> operator-(const Type x) {
		Matrix2<Type> temp(NI, NJ);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] - x; }
		return temp;
	}


	inline Matrix2<Type> operator*(const Type x) {
		Matrix2<Type> temp(NI, NJ);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] * x; }
		return temp;
	}


	inline Matrix2<Type> operator/(const Type x) {
		Matrix2<Type> temp(NI, NJ);
		try {
			for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] / x; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return temp;
	}


};


template<class Type> class Matrix3
{
public:
	int I0, I1, J0, J1, K0, K1, NI, NJ, NK, Nv;
	Type *v_p;

	Matrix3(int i, int j, int k) {
		NI = i;
		NJ = j;
		NK = k;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		I1 = NI - 1;
		J1 = NJ - 1;
		K1 = NK - 1;
		Nv = NI * NJ * NK;
		v_p = new Type[Nv];
#if(INI)
		for (int i = I0; i < Nv; ++i) { v_p[i] = NAN; }
#else
		for (int i = I0; i < Nv; ++i) { v_p[i] = Type(0); }
#endif
	}


	Matrix3() {
		NI = 0;
		NJ = 0;
		NK = 0;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		I1 = 0;
		J1 = 0;
		K1 = 0;
		Nv = 0;
		v_p = new Type[Nv];
	}


	Matrix3(const Matrix3<Type> &A) {
		cout << "Copy constructor." << endl;
		NI = A.NI;
		NJ = A.NJ;
		NK = A.NK;
		I0 = A.I0;
		I1 = A.I1;
		J0 = A.J0;
		J1 = A.J1;
		K0 = A.K0;
		K1 = A.K1;
		Nv = A.Nv;
		v_p = new Type[Nv];
		*v_p = *A.v_p;
		for (int i = I0; i < Nv; ++i) { v_p[i] = A.v_p[i]; }
	}


	~Matrix3() {
		NI = 0;
		NJ = 0;
		NK = 0;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		I1 = 0;
		J1 = 0;
		K1 = 0;
		Nv = 0;
		if (v_p != NULL) delete[] v_p;
		//cout << "[] is deleted." << endl;
	}


	inline void check_bound(int i, int j, int k) {
#ifdef BOUNDS_CHECK
		try {
			if ((i < I0) || (i > I1) || (j < J0) || (j > J1) || (k < K0) || (k > K1)) throw 99;
		}
		catch (int id) {
			printf("Bound error: i, j, k = %d, %d, %d is out of [%d,%d, %d,%d, %d,%d]. \n", i, j, k, I0, I1, J0, J1, K0, K1);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
	}


	inline void check_shape(const Matrix3<Type> &A) {
#ifdef BOUNDS_CHECK
		try {
			if (this->Nv != A.Nv) throw 99;
			else if (this->NI != A.NI) throw 99;
			else if (this->NJ != A.NJ) throw 99;
			else if (this->NK != A.NK) throw 99;
		}
		catch (int i) {
			printf("Shape error: %d %d %d is not consisten with %d %d %d. \n", NI, NJ, NK, A.NI, A.NJ, A.NK);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif	
	}


	inline void allocate(int i, int j, int k) {
		if (!Nv) {
			delete[] v_p;
			NI = i;
			NJ = j;
			NK = k;
			I0 = 0;
			J0 = 0;
			K0 = 0;
			I1 = NI - 1;
			J1 = NJ - 1;
			K1 = NK - 1;
			Nv = NI * NJ * NK;
			v_p = new Type[Nv];
		}
		else {
			cout << "Matrix has been allocated." << endl;
			//system("pause");
			//cout << "Go" << endl;
		}
#if(INI)
		for (int i = I0; i < Nv; ++i) { v_p[i] = NAN; }
#else
		for (int i = I0; i < Nv; ++i) { v_p[i] = Type(0); }
#endif
	}


	inline void deallocate() {
		NI = 0;
		NJ = 0;
		NK = 0;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		I1 = 0;
		J1 = 0;
		K1 = 0;
		Nv = 0;
		delete[] v_p;
		v_p = new Type[0];
	}


	inline void setvalue(Type x) {
		for (int i = Nv - 1; i >= 0; --i) { v_p[i] = x; }
	}


	inline Type sum(void) {
		Type s = 0;
		for (int i = Nv - 1; i >= 0; --i) { s += v_p[i]; }
		return s;
	}


	inline Matrix3<Type> sum(int index0, int index1, int index3, int axis);


	inline Matrix3<Type> transpose(int index0, int index1, int index2);


	inline Matrix3<Type> msqrt(void) {
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) {
#if _DEBUG
			if (v_p[i] < 0) {
				cout << "Negitive number got in sqrt." << endl;
				system("pause");
			}
#endif // _DEBUG
			temp.v_p[i] = sqrt(v_p[i]);
		}
		return temp;
	}


	inline void output(int presc = 2) {
		int ni, nj, nk;
		ni = NI;
		nj = NJ;
		nk = NK;
		cout.precision(presc);
		for (int k = 0; k < nk; ++k) {
			printf("k = %d \n", k);
			for (int i = 0; i < ni; ++i) {
				for (int j = 0; j < nj; ++j) {
					cout << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
					cout << (*this)(i, j, k) << " ";
				}
				cout << endl;
			}
			cout << endl << endl;
		}
		cout << endl;
	}


	inline void output(string filename, int presc = 2) {
		int ni, nj, nk;
		//ofstream OutFile(filename);
		ni = NI;
		nj = NJ;
		nk = NK;
		string _file;
		for (int k = 0; k < nk; ++k)
		{
			//printf("k = %d \n", k);
			//OutFile << "k = " << k << endl;
			_file = "id" + std::to_string(k) + "-" + filename;
			ofstream OutFile(_file, std::ios::out | std::ios::trunc);
			if (OutFile)
			{
				OutFile.precision(presc);
				for (int i = 0; i < ni; ++i) {
					for (int j = 0; j < nj; ++j) {
						OutFile << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
						OutFile << (*this)(i, j, k) << "\t";
					}
					OutFile << endl;
				}
				OutFile << endl << endl;
				OutFile.close();
			}
			else {
				cout << filename << " open failed." << endl;
				system("pause");
			}
		}
	}


	inline void output(ofstream &outfile) {
		int ni, nj, nk;
		ni = NI;
		nj = NJ;
		nk = NK;
		if (outfile) {
			for (int k = 0; k < nk; ++k) {
				printf("k = %d \n", k);
				for (int i = 0; i < ni; ++i) {
					for (int j = 0; j < nj; ++j) {
						outfile << std::right << std::setw(presc + 2) << std::setprecision(presc) << std::fixed << std::showpoint;
						outfile << (*this)(i, j, k) << "\t";
					}
					outfile << endl;
				}
				outfile << endl << endl;
			}
			outfile.close();
		}
		else {
			cout << "File open failed." << endl;
			system("pause");
		}

	}


	inline void output(int i0, int i1, int j0, int j1, int k0, int k1, string filename);


	inline void output(int i0, int i1, int j0, int j1, int k0, int k1, ofstream &outfile);


	inline void input(string filename) {
		int m, n, k;
		ifstream InFile(filename);
		InFile >> m;
		InFile >> n;
		InFile >> k;
#ifdef BOUNDS_CHECK
		try {
			if (!((m == NI) && (n == NJ) && (k == NK))) throw 99;
		}
		catch (int id) {
			printf("Shape error: input %d %d %d is not consistent with %d %d %d. \n", m, n, k, NI, NJ, NK);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
		if (InFile) {
			for (int ik = 0; ik < k; ++ik) {
				for (int ii = 0; ii < m; ++ii) {
					for (int ij = 0; ij < n; ++ij) {
						InFile >> (*this)(ii, ij, ik);
					}
				}
			}
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}
	}


	inline void input(string filename, Matrix3<Type> &A) {
		int m, n, k;
		ifstream InFile(filename);
		InFile >> m;
		InFile >> n;
		InFile >> k;
#ifdef BOUNDS_CHECK
		try {
			if (!((m == A.NI) && (n == A.NJ) && (k == A.NK))) throw 99;
		}
		catch (int id) {
			printf("Shape error: input %d %d %d is not consistent with %d %d %d. \n", m, n, k, A.NI, A.NJ, A.NK);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
		if (InFile) {
			A.deallocate();
			A.allocate(m, n, k);
			for (int ik = 0; ik < k; ++i) {
				for (int ii = 0; ii < m; ++ii) {
					for (int ij = 0; ij < n; ++ij) {
						InFile >> A(ii, ij, ik);
					}
				}
			}
		}
		else {
			cout << filename << " open failed." << endl;
			system("pause");
		}
	}


	inline Matrix2<Type> reshape(int aixs, int ni, int nj) {
		Matrix2<Type> temp(ni, nj);
#ifdef BOUNDS_CHECK
		try {
			if (ni*nj != Nv) throw 99;
		}
		catch (int id) {
			printf("Shape error: reshape (%d, %d, %d) to (%d, %d).", NI, NJ, NK, ni, nj);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif // BOUNDS_CHECK

		if (aixs == 3) {
			for (int i = ni - 1; i >= 0; --i) {
				for (int j = nj - 1; j >= 0; --j) {
					temp(i, j) = *(v_p + i*nj + j);
				}
			}
		}
		return temp;
	}


	// operator overloading
	inline Type & operator()(int i, int j, int k) {
		check_bound(i, j, k);
		return *(v_p + k*NI*NJ + j*NI + i);
	}


	inline Matrix1<Type> operator()(int axis, const Matrix1<int> &A);


	inline Matrix2<Type> operator()(int axis0, int axis1, const Matrix2<int> &A);


	inline Matrix2<Type> operator()(int axis, int id) {
		Matrix2<Type> temp;
		if (axis == 0) {
			temp.allocate(NJ, NK);
#ifdef BOUNDS_CHECK
			if (id > NI || id < 0) {
				printf("Bound error");
				system("pause");
			}
#endif // BOUND_CHECK

			for (int i = NJ - 1; i >= 0; --i) {
				for (int j = NK - 1; j >= 0; --j) {
					temp(i, j) = (*this)(id, i, j);
				}
			}
		}
		else if (axis == 1) {
			temp.allocate(NI, NK);
#ifdef BOUNDS_CHECK
			if (id > NJ || id < 0) {
				printf("Bound error");
				system("pause");
			}
#endif // BOUND_CHECK
			for (int i = NI - 1; i >= 0; --i) {
				for (int j = NK - 1; j >= 0; --j) {
					temp(i, j) = (*this)(i, id, j);
				}
			}
		}
		else if (axis == 2) {
			temp.allocate(NI, NJ);
#ifdef BOUNDS_CHECK
			if (id > NK || id < 0) {
				printf("Bound error");
				system("pause");
			}
#endif // BOUND_CHECK
			for (int i = NI - 1; i >= 0; --i) {
				for (int j = NJ - 1; j >= 0; --j) {
					temp(i, j) = (*this)(i, j, id);
				}
			}
		}
		else {
			//Matrix2<Type> temp(1, 1);
			cout << "Wrong axis gotten in Matrix3 data drawn." << endl;
			system("pause");
		}
		return temp;
	}


	inline Matrix3<Type> operator()(const Matrix3<int> &A);


	inline Matrix3<Type> & operator=(const Matrix3<Type> &A) {
		//check_shape(A);	
		if (~((this->Nv == A.Nv) && (this->NI == A.NI) && (this->NJ == A.NJ))) {
			this->deallocate();
			this->allocate(A.NI, A.NJ, A.NK);
		}
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] = A.v_p[i]; }
		return (*this);
	}


	inline Matrix3<Type> & operator+=(const Matrix3<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] += A.v_p[i]; }
		return (*this);
	}


	inline Matrix3<Type> & operator-=(const Matrix3<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] -= A.v_p[i]; }
		return (*this);
	}


	inline Matrix3<Type> & operator*=(const Matrix3<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] *= A.v_p[i]; }
		return (*this);
	}


	inline Matrix3<Type> & operator/=(const Matrix3<Type> &A) {
		check_shape(A);
		try {
			for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] /= A.v_p[i]; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return (*this);
	}


	inline Matrix3<Type> operator+(const Matrix3<Type> &A) {
		check_shape(A);
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] + A.v_p[i];
		}
		return temp;
	}


	inline Matrix3<Type> operator-(const Matrix3<Type> &A) {
		check_shape(A);
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] - A.v_p[i];
		}
		return temp;
	}


	inline Matrix3<Type> operator*(const Matrix3<Type> &A) {
		check_shape(A);
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] * A.v_p[i];
		}
		return temp;
	}


	inline Matrix3<Type> operator/(const Matrix3<Type> &A) {
		check_shape(A);
		Matrix3<Type> temp(NI, NJ, NK);
		try {
			for (int i = Nv - 1; i >= 0; --i) {
				temp.v_p[i] = v_p[i] / A.v_p[i];
			}
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return temp;
	}


	inline Matrix3<Type> operator+(const Type x) {
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] + x; }
		return temp;
	}


	inline Matrix3<Type> operator-(const Type x) {
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] - x; }
		return temp;
	}


	inline Matrix3<Type> operator*(const Type x) {
		Matrix3<Type> temp(NI, NJ, NK);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] * x; }
		return temp;
	}


	inline Matrix3<Type> operator/(const Type x) {
		Matrix3<Type> temp(NI, NJ, NK);
		try {
			for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] / x; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return temp;
	}
};


template<class Type> class Matrix4
{
public:
	int I0, I1, J0, J1, K0, K1, W0, W1, NI, NJ, NK, NW, Nv;
	Type *v_p;

	Matrix4(int i, int j, int k, int w) {
		NI = i;
		NJ = j;
		NK = k;
		NW = w;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		W0 = 0;
		I1 = NI - 1;
		J1 = NJ - 1;
		K1 = NK - 1;
		W1 = NW - 1;
		Nv = NI * NJ * NK * NW;
		v_p = new Type[Nv];
#if(INI)
		for (int i = I0; i < Nv; ++i) { v_p[i] = NAN; }
#else
		for (int i = I0; i < Nv; ++i) { v_p[i] = Type(0); }
#endif
	}


	Matrix4() {
		NI = 0;
		NJ = 0;
		NK = 0;
		NW = 0;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		W0 = 0;
		I1 = 0;
		J1 = 0;
		K1 = 0;
		W1 = 0;
		Nv = 0;
		v_p = new Type[Nv];
	}


	Matrix4(const Matrix4<Type> &A) {
		cout << "Copy constructor." << endl;
		NI = A.NI;
		NJ = A.NJ;
		NK = A.NK;
		NW = A.NW;
		I0 = A.I0;
		I1 = A.I1;
		J0 = A.J0;
		J1 = A.J1;
		K0 = A.K0;
		K1 = A.K1;
		W0 = A.W0;
		W1 = A.W1;
		Nv = A.Nv;
		v_p = new Type[Nv];
		*v_p = *A.v_p;
		for (int i = I0; i < Nv; ++i) { v_p[i] = A.v_p[i]; }
	}


	~Matrix4() {
		NI = 0;
		NJ = 0;
		NK = 0;
		NW = 0;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		W0 = 0;
		I1 = 0;
		J1 = 0;
		K1 = 0;
		W1 = 0;
		Nv = 0;
		if (v_p != NULL) delete[] v_p;
		//cout << "[] is deleted." << endl;
	}


	inline void check_bound(int i, int j, int k, int w) {
#ifdef BOUNDS_CHECK
		try {
			if ((i < I0) || (i > I1) || (j < J0) || (j > J1) || (k < K0) || (k > K1) || (w < W0) || (w > W1)) throw 99;
		}
		catch (int id) {
			printf("Bound error: i, j, k, w = %d, %d, %d, %d is out of [%d,%d, %d,%d, %d,%d, %d,%d]. \n", i, j, k, w, I0, I1, J0, J1, K0, K1, W0, W1);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif
	}


	inline void check_shape(const Matrix4<Type> &A) {
#ifdef BOUNDS_CHECK
		try {
			if (this->Nv != A.Nv) throw 99;
			else if (this->NI != A.NI) throw 99;
			else if (this->NJ != A.NJ) throw 99;
			else if (this->NK != A.NK) throw 99;
			else if (this->NW != A.NW) throw 99;
		}
		catch (int i) {
			printf("Shape error: %d %d %d %d is not consisten with %d %d %d %d. \n", NI, NJ, NK, NW, A.NI, A.NJ, A.NK, A.NW);
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause");
		}
#endif	
	}


	inline void allocate(int i, int j, int k, int w) {
		if (!Nv) {
			delete[] v_p;
			NI = i;
			NJ = j;
			NK = k;
			NW = w;
			I0 = 0;
			J0 = 0;
			K0 = 0;
			W0 = 0;
			I1 = NI - 1;
			J1 = NJ - 1;
			K1 = NK - 1;
			W1 = NW - 1;
			Nv = NI * NJ * NK * NW;
			v_p = new Type[Nv];
		}
		else {
			cout << "Matrix has been allocated." << endl;
			//system("pause");
			//cout << "Go" << endl;
		}
#if(INI)
		for (int i = Nv - 1; i >= 0; --i) { v_p[i] = NAN; }
#else
		for (int i = Nv - 1; i >= 0; --i) { v_p[i] = Type(0); }
#endif
	}


	inline void deallocate() {
		NI = 0;
		NJ = 0;
		NK = 0;
		NW = 0;
		I0 = 0;
		J0 = 0;
		K0 = 0;
		W0 = 0;
		I1 = 0;
		J1 = 0;
		K1 = 0;
		W1 = 0;
		Nv = 0;
		delete[] v_p;
		v_p = new Type[0];
	}


	inline void setvalue(Type x) {
		for (int i = Nv - 1; i >= 0; --i) { v_p[i] = x; }
	}


	inline Type sum(void) {
		Type s = 0;
		for (int i = Nv - 1; i >= 0; --i) { s += v_p[i]; }
		return s;
	}



	inline Matrix4<Type> msqrt(void) {
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) {
#if _DEBUG
			if (v_p[i] < 0) {
				cout << "Negitive number got in sqrt." << endl;
				system("pause");
			}
#endif // _DEBUG
			temp.v_p[i] = sqrt(v_p[i]);
		}
		return temp;
	}


	inline void output(int presc = 2) {
		int ni, nj, nk, nw;
		ni = NI;
		nj = NJ;
		nk = NK;
		nw = NW;
		cout.precision(presc);
		for (int w = 0; w < nw; ++w) {
			for (int k = 0; k < nk; ++k) {
				printf("w = %d, k = %d \n", w, k);
				for (int i = 0; i < ni; ++i) {
					for (int j = 0; j < nj; ++j) {
						cout << (*this)(i, j, k, w) << " ";
					}
					cout << endl;
				}
				cout << endl << endl;
			}
		}

	}


	inline void output(string filename) {
		int ni, nj, nk, nw;
		ofstream OutFile(filename);
		OutFile.precision(10);
		ni = NI;
		nj = NJ;
		nk = NK;
		nw = NW;
		for (int w = 0; w < nw; ++w) {
			for (int k = 0; k < nk; ++k) {
				printf("w = %d, k = %d \n", w, k);
				for (int i = 0; i < ni; ++i) {
					for (int j = 0; j < nj; ++j) {
						OutFile << (*this)(i, j, k, w) << "\t";
					}
					OutFile << endl;
				}
				OutFile << endl << endl;
			}
		}


		OutFile.close();
	}


	inline void output(ofstream &outfile) {
		int ni, nj, nk, nw;
		ni = NI;
		nj = NJ;
		nk = NK;
		nw = NW;
		for (int w = 0; w < nw; ++w) {
			for (int k = 0; k < nk; ++k) {
				printf("w = %d, k = %d \n", w, k);
				for (int i = 0; i < ni; ++i) {
					for (int j = 0; j < nj; ++j) {
						OutFile << (*this)(i, j, k, w) << "\t";
					}
					outfile << endl;
				}
				outfile << endl << endl;
			}
		}

		outfile.close();
	}


	inline void output(int i0, int i1, int j0, int j1, int k0, int k1, string filename);


	inline void output(int i0, int i1, int j0, int j1, int k0, int k1, ofstream &outfile);


	// operator overloading
	inline Type & operator()(int i, int j, int k, int w) {
		check_bound(i, j, k, w);
		return *(v_p + w*NI*NJ*NK + k*NI*NJ + j*NI + i);
	}


	inline Matrix1<Type> operator()(int axis, const Matrix1<int> &A);


	inline Matrix2<Type> operator()(int axis0, int axis1, const Matrix2<int> &A);


	inline Matrix4<Type> operator()(const Matrix4<int> &A);


	inline Matrix4<Type> & operator=(const Matrix4<Type> &A) {
		//check_shape(A);	
		if (~((this->Nv == A.Nv) && (this->NI == A.NI) && (this->NJ == A.NJ) && (this->NK == A.NK))) {
			this->deallocate();
			this->allocate(A.NI, A.NJ);
		}
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] = A.v_p[i]; }
		return (*this);
	}


	inline Matrix4<Type> & operator+=(const Matrix4<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] += A.v_p[i]; }
		return (*this);
	}


	inline Matrix4<Type> & operator-=(const Matrix4<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] -= A.v_p[i]; }
		return (*this);
	}


	inline Matrix4<Type> & operator*=(const Matrix4<Type> &A) {
		check_shape(A);
		for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] *= A.v_p[i]; }
		return (*this);
	}


	inline Matrix4<Type> & operator/=(const Matrix4<Type> &A) {
		check_shape(A);
		try {
			for (int i = Nv - 1; i >= 0; --i) { this->v_p[i] /= A.v_p[i]; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return (*this);
	}


	inline Matrix4<Type> operator+(const Matrix4<Type> &A) {
		check_shape(A);
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] + A.v_p[i];
		}
		return temp;
	}


	inline Matrix4<Type> operator-(const Matrix4<Type> &A) {
		check_shape(A);
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] - A.v_p[i];
		}
		return temp;
	}


	inline Matrix4<Type> operator*(const Matrix4<Type> &A) {
		check_shape(A);
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] * A.v_p[i];
		}
		return temp;
	}


	inline Matrix4<Type> operator/(const Matrix4<Type> &A) {
		check_shape(A);
		Matrix4<Type> temp(NI, NJ, NK, NW);
		try {
			for (int i = Nv - 1; i >= 0; --i) {
				temp.v_p[i] = v_p[i] / A.v_p[i];
			}
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return temp;
	}


	inline Matrix4<Type> operator+(const Type x) {
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] + x; }
		return temp;
	}


	inline Matrix4<Type> operator-(const Type x) {
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] - x; }
		return temp;
	}


	inline Matrix4<Type> operator*(const Type x) {
		Matrix4<Type> temp(NI, NJ, NK, NW);
		for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] * x; }
		return temp;
	}


	inline Matrix4<Type> operator/(const Type x) {
		Matrix4<Type> temp(NI, NJ, NK, NW);
		try {
			for (int i = Nv - 1; i >= 0; --i) { temp.v_p[i] = v_p[i] / x; }
		}
		catch (exception &e) {
			cerr << "Exception caught: " << e.what() << endl;
			system("pause");
		}
		catch (...) {
			cerr << "Exception caught." << endl;
			system("pause"); //exit(EXIT_FAILURE);
		}
		return temp;
	}
};



template<class Type> class SpMtrx {
	Type aa = 0;//用于()重载的时候返回零元素地址，防止野指针
public:
	Type *v;
	int *ia, *ja; // 行偏移, 列编号，均从1开始计数

	//稀疏矩阵的类型：
	//0（默认）对称稀疏矩阵，只存储右上半
	//1 非对称矩阵存储全部
	//在做矩阵乘以向量和求解线性方程组的时候是要区别对待
	//其余时候，特别是（）取值操作时一样对待，
	//也就是说，对于SpMtxType = 0的对称矩阵，当取(I,J)时，若I>J则返回0，并不会取（J,I）的值
	int SpMtxType = 0;


	int Nv, N; // 元素总数，行数


	SpMtrx() {
		Nv = 0;
		N = 0;
		ia = new int[0];
		ja = new int[0];
		v = new Type[0];
	};

	SpMtrx(int n, int nv) {
		Nv = nv;
		N = n;
		ia = new int[N + 1];
		ja = new int[Nv];
		v = new Type[Nv];
	}

	SpMtrx(const SpMtrx<Type> &A) {
		Nv = A.Nv;
		N = A.N;
		ia = new int[N + 1];
		ja = new int[Nv];
		v = new Type[Nv];
		for (int i = 0; i < Nv; i++) {
			v[i] = A.v[i];
			ja[i] = A.ja[i];
		}
		for (int i = 0; i < N + 1; i++) {
			ia[i] = A.ia[i];
		}
	}

	~SpMtrx() {
		Nv = 0;
		N = 0;
		delete[] ia, ja, v;
	};

	void allocate(int n, int nv) {
		delete[] ia, ja, v;
		Nv = nv;
		N = n;
		ja = new int[Nv];
		v = new Type[Nv];
		ia = new int[N + 1];

	}

	void allocate(int n, int nv, Type *vin, int *iain, int *jain) {
		delete[] ia, ja, v;
		Nv = nv;
		N = n;
		ia = new int[N + 1];
		ja = new int[Nv];
		v = new Type[Nv];
		for (int i = 0; i < Nv; i++) {
			v[i] = vin[i];
			ja[i] = jain[i];
		}
		for (int i = 0; i < N + 1; i++) {
			ia[i] = iain[i];
		}

	}

	void deallocate() {
		Nv = 0;
		N = 0;
		delete[] ia, ja, v;
		v = new Type[0];
		ia = new int[0];
		ja = new int[0];
	}

	void refine() {
		//剔除稀疏矩阵中的零元
		int Nv2 = Nv;
		for (int i = 0; i < Nv; i++) {
			if (v[i] == 0) {
				ja[i] = 0;
				Nv2--;
			}
		}
		double *v2;
		int k = 0;
		v2 = new double[Nv2];
		for (int i = 0; i < Nv; i++) {
			if (ja[i] != 0) {
				v2[k] = v[i];
				k++;
			}
		}

		delete[] v;
		v = v2;//指针操作

		int *ia2;
		ia2 = new int[N + 1];
		ia2[0] = 1;
		for (int i = 0; i < N; i++) {
			ia2[i + 1] = ia[i + 1];
			for (int j = ja[0]; j < ia[i + 1]; j++) {
				if (ja[j - 1] == 0) {
					ia2[i + 1]--;
				}
			}
		}
		delete[] ia;
		ia = ia2;


		int *ja2;
		k = 0;
		ja2 = new int[Nv2];
		for (int i = 0; i < Nv; i++) {
			if (ja[i] != 0) {
				ja2[k] = ja[i];
				k++;
			}
		}

		delete[] ja;
		ja = ja2;

		Nv = Nv2;
	}

	inline Type & operator()(int I) {
#ifdef BOUNDS_CHECK
		if (I<1) {
			cout << "Error: I = " << I << " < 1 " << " in SpMatrix operator()" << endl;
			system("pause");
			I = 1;
		}
		if (I>Nv) {
			cout << "Error: I = " << I << " > Nv = " << Nv << " in SpMatrix operator()" << endl;
			system("pause");
			I = Nv;
		}
#endif
		return v[I - 1];
	}

	inline Type & operator()(int I, int J) {
#ifdef BOUNDS_CHECK

		if (I < 1) {
			cout << "Error: I = " << I << " < 1 " << " in SpMtx" << endl;
			system("pause");
			I = 1;
		}
		if (I > N) {
			cout << "Error: I = " << I << " > N = " << N << " in SpMtx" << endl;
			system("pause");
			I = N;
		}
		if (J < 1) {
			cout << "Error: J = " << J << " < 1 = " << " in SpMtx" << endl;
			system("pause");
			J = 1;
		}
		if (J > N) {
			cout << "Error: J = " << J << " > N = " << N << " in SpMtx" << endl;
			system("pause");
			J = N;
		}
#endif

		for (int ii = ia[I - 1]; ii < ia[I]; ii++) {
			if (J == ja[ii - 1]) {
				return v[ii - 1];
			}
		}

#ifdef BOUNDS_CHECK
		if (aa != 0) {
			//上次对零元强行赋值，并没有赋上
			cout << "Warning: Sparse Element wasn't set: ";
			cout << "  I = " << I << "  J = " << J;
			getchar();
		}
#endif
		aa = 0;//一定要将aa清零
		return aa;
	}

	inline int & Ia(int I) {
#ifdef BOUNDS_CHECK
		if (I<1) {
			cout << "Error: I = " << I << " < 1 " << " in SpMatrix Ia(int I)" << endl;
			system("pause");
			I = 1;
		}
		if (I>N + 1) {
			cout << "Error: I = " << I << " > N+1 = " << N + 1 << " in SpMatrix Ia(int I)" << endl;
			system("pause");
			I = N + 1;
		}

#endif
		return ia[I - 1];
	}

	inline int & Ja(int I) {
#ifdef BOUNDS_CHECK
		if (I<1) {
			cout << "Error: I = " << I << " < 1 " << " in SpMatrix Ja(int I)" << endl;
			system("pause");
			I = 1;
		}
		if (I>Nv) {
			cout << "Error: I = " << I << " > Nv = " << Nv << " in SpMatrix Ja(int I)" << endl;
			system("pause");
			I = Nv;
		}
#endif
		return ja[I - 1];
	}

	inline void Output(string filename) {
		ofstream Outfile(filename);
		Outfile.precision(10);
		for (int i = 0; i < Nv; i++) {
			Outfile << v[i] << endl;
		}
		for (int i = 0; i < N + 1; i++) {
			Outfile << ia[i] << endl;
		}
		for (int i = 0; i < Nv; i++) {
			Outfile << ja[i] << endl;
		}
		Outfile.close();
	}

	inline SpMtrx<Type> & operator+=(const Type x) {
		for (int i = 0; i < Nv; i++) {
			this->v[i] += x;
		}

		return *this;
	};

	inline SpMtrx<Type> & operator-=(const Type x) {
		for (int i = 0; i < Nv; i++) {
			this->v[i] -= x;
		}

		return *this;
	};

	inline SpMtrx<Type> & operator*=(const Type x) {
		for (int i = 0; i < Nv; i++) {
			this->v[i] *= x;
		}

		return *this;
	};
	inline SpMtrx<Type> & operator/=(const Type x) {
		for (int i = 0; i < Nv; i++) {
			this->v[i] /= x;
		}

		return *this;
	};

	inline SpMtrx<Type> & operator=(const SpMtrx<Type> & A) {
		if (this->Nv != A.Nv) {
			if (this->Nv != 0) {
				cout << "Warning: \"operator=\" in Matrix1 may be wrong, because Nv is not equal" << endl;
				system("pause");
			}
			this->deallocate();
			this->allocate(A.N, A.Nv);
		}
		this->SpMtxType = A.SpMtxType;
		for (int i = 0; i <= A.N; i++) {
			this->ia[i] = A.ia[i];
		}

		for (int i = 0; i < A.Nv; i++) {
			this->v[i] = A.v[i];
			this->ja[i] = A.ja[i];
		}

		return (*this);
	};
};

inline SpMtrx<double> operator+(const SpMtrx<double> &A, const  SpMtrx<double> &B) {

	SpMtrx<double> Temp, c;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N && A.N*B.N != 0) {
		cout << "Wrong in SpMtrx operator+: A.N != B.N && A.N*B.N != 0" << endl;
		system("pause");
	};
	if (A.N == 0 && A.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: A.N == 0 && A.Nv != 0" << endl;
		system("pause");
	}
	if (B.N == 0 && B.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: B.N == 0 && B.Nv != 0" << endl;
		system("pause");
	}
#endif
	if (A.Nv == 0) {
		c = B;
		return c;
	}
	if (B.Nv == 0) {
		c = A;
		return c;
	}
	N = A.N;
	Nv = A.Nv + B.Nv + 1;
	Temp.allocate(N, Nv);
	double c0 = 1.0;
	char trans = 'N';
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	if (ierr != 0) {
		cout << "Wrong in SpMtrx<double> operator+: ierr = " << ierr << endl;
		system("pause");
	}
	Nv = Temp.Ia(N + 1) - 1;
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
	return c;
};
inline void Mplus(SpMtrx<double> &c, const SpMtrx<double> &A, const SpMtrx<double> &B) {
	SpMtrx<double> Temp;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N && A.N*B.N != 0) {
		cout << "Wrong in SpMtrx operator+: A.N != B.N && A.N*B.N != 0" << endl;
		system("pause");
	};
	if (A.N == 0 && A.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: A.N == 0 && A.Nv != 0" << endl;
		system("pause");
	}
	if (B.N == 0 && B.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: B.N == 0 && B.Nv != 0" << endl;
		system("pause");
	}
#endif
	if (A.Nv == 0) {
		c = B;
	}
	if (B.Nv == 0) {
		c = A;
	}
	N = A.N;
	Nv = A.Nv + B.Nv + 1;
	Temp.allocate(N, Nv);
	double c0 = 1.0;
	char trans = 'N';
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	if (ierr != 0) {
		cout << "Wrong in SpMtrx<double> operator+: ierr = " << ierr << endl;
		system("pause");
	}
	Nv = Temp.Ia(N + 1) - 1;
	c.deallocate();
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
};
inline void Mplus(SpMtrx<double> &c, const SpMtrx<double> &A, const SpMtrx<double> &B, double c0) {
	//c=A+c0*B;
	SpMtrx<double> Temp;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N && A.N*B.N != 0) {
		cout << "Wrong in SpMtrx operator+: A.N != B.N && A.N*B.N != 0" << endl;
		system("pause");
	};
	if (A.N == 0 && A.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: A.N == 0 && A.Nv != 0" << endl;
		system("pause");
	}
	if (B.N == 0 && B.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: B.N == 0 && B.Nv != 0" << endl;
		system("pause");
	}
#endif
	if (A.Nv == 0) {
		c = B;
	}
	if (B.Nv == 0) {
		c = A;
	}
	N = A.N;
	Nv = A.Nv + B.Nv + 1;
	Temp.allocate(N, Nv);
	char trans = 'N';
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	if (ierr != 0) {
		cout << "Wrong in SpMtrx<double> operator+: ierr = " << ierr << endl;
		system("pause");
	}
	Nv = Temp.Ia(N + 1) - 1;
	c.deallocate();
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
};
inline void Mplus(SpMtrx<double> &c, const SpMtrx<double> &A, const SpMtrx<double> &B, double c0, double d0) {
	//c=c0*A+d0*B;
	if (c0 != 0) {
		double d1 = d0 / c0;
		//c = (A + d1*B)*c0
		Mplus(c, A, B, d1);
		c *= c0;
	}
	else if (d0 != 0) {
		double c1 = c0 / d0;
		//c = (c1*A + B)*d0
		Mplus(c, B, A, c1);
		c *= d0;
	}
	else {
		c *= 0;
	}
}

inline SpMtrx<double> MplusT(const SpMtrx<double> &A, const  SpMtrx<double> &B) {
	//计算 c=A+B^(T)
	SpMtrx<double> Temp, c;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N && A.N*B.N != 0) {
		cout << "Wrong in SpMtrx operator+: A.N != B.N && A.N*B.N != 0" << endl;
		system("pause");
	};
	if (A.N == 0 && A.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: A.N == 0 && A.Nv != 0" << endl;
		system("pause");
	}
	if (B.N == 0 && B.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: B.N == 0 && B.Nv != 0" << endl;
		system("pause");
	}
#endif
	if (A.Nv == 0) {
		c = B;
		return c;
	}
	if (B.Nv == 0) {
		c = A;
		return c;
	}
	N = A.N;
	Nv = A.Nv + B.Nv + 1;
	Temp.allocate(N, Nv);
	double c0 = 1.0;
	char trans = 'T';//
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	if (ierr != 0) {
		cout << "Wrong in SpMtrxAddT: ierr = " << ierr << endl;
		system("pause");
	}
	Nv = Temp.Ia(N + 1) - 1;
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
	return c;
};
inline void MplusT(SpMtrx<double> &c, const SpMtrx<double> &A, const  SpMtrx<double> &B) {
	//计算 c=A+B^(T)
	SpMtrx<double> Temp;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N && A.N*B.N != 0) {
		cout << "Wrong in SpMtrx operator+: A.N != B.N && A.N*B.N != 0" << endl;
		system("pause");
	};
	if (A.N == 0 && A.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: A.N == 0 && A.Nv != 0" << endl;
		system("pause");
	}
	if (B.N == 0 && B.Nv != 0) {
		cout << "Wrong in SpMtrx operator+: B.N == 0 && B.Nv != 0" << endl;
		system("pause");
	}
#endif
	if (A.Nv == 0) {
		c = B;
	}
	if (B.Nv == 0) {
		c = A;
	}
	N = A.N;
	Nv = A.Nv + B.Nv + 1;
	Temp.allocate(N, Nv);
	double c0 = 1.0;
	char trans = 'T';//
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	if (ierr != 0) {
		cout << "Wrong in SpMtrxAddT: ierr = " << ierr << endl;
		system("pause");
	}
	Nv = Temp.Ia(N + 1) - 1;
	c.deallocate();
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
};

inline SpMtrx<double> operator-(const SpMtrx<double> &A, const  SpMtrx<double> &B) {

	SpMtrx<double> Temp, c;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N) {
		cout << "Error: A.N = " << A.N << " != B.N = " << B.N << " in SpMtrx<double> operator- (const SpMtrx<double> &A, const  SpMtrx<double> &B)" << endl;
		system("pause");
		return A;
	};
#endif
	N = A.N;
	Nv = A.Nv + B.Nv;
	Temp.allocate(N, Nv);
	double c0 = -1.0;
	char trans = 'N';
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	Nv = Temp.Ia(N + 1) - 1;
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
	return c;
};
inline void Msubs(SpMtrx<double> &c, const SpMtrx<double> &A, const  SpMtrx<double> &B) {

	SpMtrx<double> Temp;
	int N, Nv;
#ifdef BOUNDS_CHECK
	if (A.N != B.N) {
		cout << "Error: A.N = " << A.N << " != B.N = " << B.N << " in SpMtrx<double> operator- (const SpMtrx<double> &A, const  SpMtrx<double> &B)" << endl;
		system("pause");
	};
#endif
	N = A.N;
	Nv = A.Nv + B.Nv;
	Temp.allocate(N, Nv);
	double c0 = -1.0;
	char trans = 'N';
	int ierr;
	int request = 0, sort = 0;
	mkl_dcsradd(&trans, &request, &sort, &N, &N, A.v, A.ja, A.ia, &c0, B.v, B.ja, B.ia, Temp.v, Temp.ja, Temp.ia, &Nv, &ierr);
	Nv = Temp.Ia(N + 1) - 1;
	c.deallocate();
	c.allocate(N, Nv);

	for (int i = 0; i <= c.N; i++) {
		c.ia[i] = Temp.ia[i];
	}

	for (int i = 0; i < c.Nv; i++) {
		c.v[i] = Temp.v[i];
		c.ja[i] = Temp.ja[i];
	}
	if (A.SpMtxType == 1 || B.SpMtxType == 1) {
		c.SpMtxType = 1;
	}
	else {
		c.SpMtxType = 0;
	}
};

inline SpMtrx<double> operator*(const double x, const SpMtrx<double> & B) {

	SpMtrx<double> Temp;

	Temp.allocate(B.N, B.Nv);

	for (int i = 0; i <= B.N; i++) {
		Temp.ia[i] = B.ia[i];
		//cout << Temp.ia[i] << endl;
	}

	for (int i = 0; i < B.Nv; i++) {
		Temp.v[i] = x*B.v[i];
		Temp.ja[i] = B.ja[i];
	}
	Temp.SpMtxType = B.SpMtxType;

	return Temp;

};


inline Matrix1<double> operator*(const SpMtrx<double> & A, const Matrix1<double> & X) {
#ifdef BOUNDS_CHECK
	if (A.N != X.Nv) {
		cout << "Error: A.N = " << A.N << " != X.Nv = " << X.Nv << " in Matrix1<double> operator*(const SpMtrx<double> & A, const Matrix1<double> & X)" << endl;
		system("pause");
		return X;
	};
#endif
	Matrix1<double> b(X.I0, X.I1);

	if (A.SpMtxType == 0) {
		//对称稀疏矩阵
		const char uplo = 'U';
		int NN = A.N;
		mkl_dcsrsymv(&uplo, &NN, A.v, A.ia, A.ja, X.v_p, b.v_p);
	}
	else if (A.SpMtxType == 1) {
		//非对称稀疏矩阵
		const char transa = 'N';
		int NN = A.N;
		mkl_dcsrgemv(&transa, &NN, A.v, A.ia, A.ja, X.v_p, b.v_p);
	}
	return b;
};
inline void Mmul(Matrix1<double> &b, const SpMtrx<double> & A, const Matrix1<double> & X) {
#ifdef BOUNDS_CHECK
	if (A.N != X.Nv) {
		cout << "Error: A.N = " << A.N << " != X.Nv = " << X.Nv << " in Matrix1<double> operator*(const SpMtrx<double> & A, const Matrix1<double> & X)" << endl;
		system("pause");
		return;
	};
#endif
	if (b.Nv == 0) {
		b.allocate(X.Nv);
	}
	else if (b.Nv != X.Nv) {
		cout << "Error: b.Nv != X.Nv in SpMtrxMtpM1!" << endl;
		system("pause");
		return;
	}

	if (A.SpMtxType == 0) {
		//对称稀疏矩阵
		const char uplo = 'U';
		int NN = A.N;
		mkl_dcsrsymv(&uplo, &NN, A.v, A.ia, A.ja, X.v_p, b.v_p);
	}
	else if (A.SpMtxType == 1) {
		//非对称稀疏矩阵
		const char transa = 'N';
		int NN = A.N;
		mkl_dcsrgemv(&transa, &NN, A.v, A.ia, A.ja, X.v_p, b.v_p);
	}
}

inline void Msolve(SpMtrx<double> &K, Matrix1<double> &F, Matrix1<double> &X) {

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
		K.Output("K.dat");
		F.output("F.dat");
		cout << "ERROR1 IN SolverMKL was detected: " << error;
		getchar();
	}

	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &K.N, K.v, K.ia, K.ja, &perm, &nrhs, iparm, &msglvl, F.v_p, X.v_p, &error);

	if (error != 0) {
		cout << "ERROR2 IN SolverMKL was detected: " << error;
		getchar();
	}
}
inline void Msolve(SpMtrx<double> &K, Matrix1<double> &F) {
	//solve X = K\X;
	__int64 pt[64];
	int maxfct = 1, mnum = 1, mtype, phase, error, msglvl;
	int iparm[64];
	int perm;
	int nrhs = 1;
	Matrix1<double> X;
	X = F;
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
	iparm[5] = 1;//= 0 solution on the first n compoments of x///**************************************
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
		K.Output("K.dat");
		F.output("F.dat");
		cout << "ERROR1 IN SolverMKL was detected: " << error;
		getchar();
	}

	phase = -1;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &K.N, K.v, K.ia, K.ja, &perm, &nrhs, iparm, &msglvl, F.v_p, X.v_p, &error);

	if (error != 0) {
		cout << "ERROR2 IN SolverMKL was detected: " << error;
		getchar();
	}
}

inline void Msolve(SpMtrx<double> &K, Matrix1<double> &F, Matrix1<double> &X, int solveid);



// operator overloading

inline Matrix1<double> Matrix1<double>::operator+(const Matrix1<double> &A) {
	check_shape(A);
	Matrix1<double> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] + A.v_p[i];
	}
	return temp;
}


inline Matrix1<double> Matrix1<double>::operator-(const Matrix1<double> &A) {
	check_shape(A);
	Matrix1<double> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] - A.v_p[i];
	}
	return temp;
}


inline Matrix1<double> Matrix1<double>::operator*(const Matrix1<double> &A) {
	check_shape(A);
	Matrix1<double> temp(A.Nv);

#ifdef _USE_MKL
	vdMul(A.Nv, v_p, A.v_p, temp.v_p);
#else
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] * A.v_p[i];
		//temp.v_p[i] = v_p[i] * A_p[i];
	}
#endif // _USE_MKL

	return temp;

}


inline Matrix1<double> Matrix1<double>::operator/(const Matrix1<double> &A) {
	check_shape(A);
	Matrix1<double> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	try {
		for (int i = i1; i >= i0; --i) {
			temp.v_p[i] = v_p[i] / A.v_p[i];
		}
	}
	catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	return temp;
}


inline Matrix1<int> Matrix1<int>::operator+(const Matrix1<int> &A) {
	check_shape(A);
	Matrix1<int> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] + A.v_p[i];
	}
	return temp;
}


inline Matrix1<int> Matrix1<int>::operator-(const Matrix1<int> &A) {
	check_shape(A);
	Matrix1<int> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] - A.v_p[i];
	}
	return temp;
}


inline Matrix1<int> Matrix1<int>::operator*(const Matrix1<int> &A) {
	check_shape(A);
	Matrix1<int> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] * A.v_p[i];
	}
	return temp;
}


inline Matrix1<int> Matrix1<int>::operator/(const Matrix1<int> &A) {
	check_shape(A);
	Matrix1<int> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	try {
		for (int i = i1; i >= i0; --i) {
			temp.v_p[i] = v_p[i] / A.v_p[i];
		}
	}
	catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	return temp;
}


inline Matrix1<float> Matrix1<float>::operator+(const Matrix1<float> &A) {
	check_shape(A);
	Matrix1<float> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] + A.v_p[i];
	}
	return temp;
}


inline Matrix1<float> Matrix1<float>::operator-(const Matrix1<float> &A) {
	check_shape(A);
	Matrix1<float> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] - A.v_p[i];
	}
	return temp;
}


inline Matrix1<float> Matrix1<float>::operator*(const Matrix1<float> &A) {
	check_shape(A);
	Matrix1<float> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	for (int i = i1; i >= i0; --i) {
		temp.v_p[i] = v_p[i] * A.v_p[i];
	}
	return temp;
}


inline Matrix1<float> Matrix1<float>::operator/(const Matrix1<float> &A) {
	check_shape(A);
	Matrix1<float> temp(A.Nv);
	int i0 = I0;
	int i1 = I1;
	try {
		for (int i = i1; i >= i0; --i) {
			temp.v_p[i] = v_p[i] / A.v_p[i];
		}
	}
	catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	return temp;
}


inline Matrix2<double> Matrix2<double>::operator+(const Matrix2<double> &A) {
	check_shape(A);
	Matrix2<double> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] + A.v_p[i];
	}
	return temp;
}


inline Matrix2<double> Matrix2<double>::operator-(const Matrix2<double> &A) {
	check_shape(A);
	Matrix2<double> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] - A.v_p[i];
	}
	return temp;
}


inline Matrix2<double> Matrix2<double>::operator*(const Matrix2<double> &A) {
	check_shape(A);
	Matrix2<double> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] * A.v_p[i];
	}
	return temp;
}


inline Matrix2<double> Matrix2<double>::operator/(const Matrix2<double> &A) {
	check_shape(A);
	Matrix2<double> temp(NI, NJ);
	int nv = Nv;
	try {
		for (int i = nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] / A.v_p[i];
		}
	}
	catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
		system("pause");
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	return temp;
}


inline Matrix2<int> Matrix2<int>::operator+(const Matrix2<int> &A) {
	check_shape(A);
	Matrix2<int> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] + A.v_p[i];
	}
	return temp;
}


inline Matrix2<int> Matrix2<int>::operator-(const Matrix2<int> &A) {
	check_shape(A);
	Matrix2<int> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] - A.v_p[i];
	}
	return temp;
}


inline Matrix2<int> Matrix2<int>::operator*(const Matrix2<int> &A) {
	check_shape(A);
	Matrix2<int> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] * A.v_p[i];
	}
	return temp;
}


inline Matrix2<int> Matrix2<int>::operator/(const Matrix2<int> &A) {
	check_shape(A);
	Matrix2<int> temp(NI, NJ);
	int nv = Nv;
	try {
		for (int i = nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] / A.v_p[i];
		}
	}
	catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
		system("pause");
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	return temp;
}


inline Matrix2<float> Matrix2<float>::operator+(const Matrix2<float> &A) {
	check_shape(A);
	Matrix2<float> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] + A.v_p[i];
	}
	return temp;
}


inline Matrix2<float> Matrix2<float>::operator-(const Matrix2<float> &A) {
	check_shape(A);
	Matrix2<float> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] - A.v_p[i];
	}
	return temp;
}


inline Matrix2<float> Matrix2<float>::operator*(const Matrix2<float> &A) {
	check_shape(A);
	Matrix2<float> temp(NI, NJ);
	int nv = Nv;
	for (int i = nv - 1; i >= 0; --i) {
		temp.v_p[i] = v_p[i] * A.v_p[i];
	}
	return temp;
}


inline Matrix2<float> Matrix2<float>::operator/(const Matrix2<float> &A) {
	check_shape(A);
	Matrix2<float> temp(NI, NJ);
	int nv = Nv;
	try {
		for (int i = nv - 1; i >= 0; --i) {
			temp.v_p[i] = v_p[i] / A.v_p[i];
		}
	}
	catch (exception &e) {
		cerr << "Exception caught: " << e.what() << endl;
		system("pause");
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	return temp;
}


inline Matrix1<float> Matrix1<float>::operator()(const Matrix1<int> &A) {
	check_bound(*(A.v_p));
	check_bound(*(A.v_p + A.Nv - 1));
	Matrix1<float> temp(A.Nv);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = *(v_p + (*(A.v_p + i))); }
	return temp;
}


inline Matrix1<double> Matrix1<double>::operator()(const Matrix1<int> &A) {
	check_bound(*(A.v_p));
	check_bound(*(A.v_p + A.Nv - 1));
	Matrix1<double> temp(A.Nv);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = *(v_p + (*(A.v_p + i))); }
	return temp;
}


inline Matrix1<int> Matrix1<int>::operator()(const Matrix1<int> &A) {
	check_bound(*(A.v_p));
	check_bound(*(A.v_p + A.Nv - 1));
	Matrix1<int> temp(A.Nv);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = *(v_p + (*(A.v_p + i))); }
	return temp;
}


inline Matrix1<float> Matrix2<float>::operator()(int i, const Matrix1<int> &A) {
	check_bound(i, *(A.v_p));
	check_bound(i, *(A.v_p + A.Nv - 1));
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int k = anv - 1; k >= 0; --k) { temp.v_p[k] = (*this)(i, A.v_p[k]); }
	return temp;
}


inline Matrix1<double> Matrix2<double>::operator()(int i, const Matrix1<int> &A) {
	check_bound(i, *(A.v_p));
	check_bound(i, *(A.v_p + A.Nv - 1));
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int k = anv - 1; k >= 0; --k) { temp.v_p[k] = (*this)(i, A.v_p[k]); }
	return temp;
}


inline Matrix1<int> Matrix2<int>::operator()(int i, const Matrix1<int> &A) {
	check_bound(i, *(A.v_p));
	check_bound(i, *(A.v_p + A.Nv - 1));
	int anv = A.Nv;
	Matrix1<int> temp(anv);
	for (int k = anv - 1; k >= 0; --k) { temp.v_p[k] = (*this)(i, A.v_p[k]); }
	return temp;
}


inline Matrix1<float> Matrix2<float>::operator()(const Matrix1<int> &A, int j) {
	check_bound(*(A.v_p), j);
	check_bound(*(A.v_p + A.Nv - 1), j);
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int k = anv - 1; k >= 0; --k) { temp.v_p[k] = (*this) (A.v_p[k], j); }
	return temp;
}


inline Matrix1<double> Matrix2<double>::operator()(const Matrix1<int> &A, int j) {
	check_bound(*(A.v_p), j);
	check_bound(*(A.v_p + A.Nv - 1), j);
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int k = anv - 1; k >= 0; --k) { temp.v_p[k] = (*this) (A.v_p[k], j); }
	return temp;
}

//template<class Type>
//inline Matrix2<Type> Matrix2<Type>::operator()(const Matrix2<int>& A)
//{
//	Matrix2<Type> temp(A.NI, A.NJ);
//#ifdef BOUNDS_CHECK
//	try {
//		if (A(0) < 0) throw 99;
//		if (A(0) > A(A.Nv-1))　throw 99;
//	}
//	catch (int i) {
//		printf_wrong_msg("Index beyond.");
//	}
//#endif // BOUNDS_CHECK
//	for (int j = A.NJ - 1; j >= 0; --j) {
//		for (int i = A.NI - 1; i >= 0; --i) {
//			temp(i, j) = this->v_p[A(i, j)];
//		}
//	}
//	return temp;
//}


inline Matrix1<int> Matrix2<int>::operator()(const Matrix1<int> &A, int j) {
	check_bound(*(A.v_p), j);
	check_bound(*(A.v_p + A.Nv - 1), j);
	int anv = A.Nv;
	Matrix1<int> temp(anv);
	for (int k = anv - 1; k >= 0; --k) { temp.v_p[k] = (*this) (A.v_p[k], j); }
	return temp;
}


// function definition
inline void check_shape(const Matrix1<double> &A, const Matrix1<double> &B) {
#ifdef BOUNDS_CHECK
	try {
		if (A.Nv != B.Nv) throw 99;
	}
	catch (int i) {
		printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
#endif // BOUNDS_CHECK
}


inline void check_shape(const Matrix1<int> &A, const Matrix1<int> &B) {
#ifdef BOUNDS_CHECK
	try {
		if (A.Nv != B.Nv) throw 99;
	}
	catch (int i) {
		printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
#endif // BOUNDS_CHECK
}


inline void check_shape(const Matrix1<float> &A, const Matrix1<float> &B) {
#ifdef BOUNDS_CHECK
	try {
		if (A.Nv != B.Nv) throw 99;
	}
	catch (int i) {
		printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
#endif // BOUNDS_CHECK
}


inline void check_shape(const Matrix2<double> &A, const Matrix2<double> &B) {
#ifdef BOUNDS_CHECK
	try {
		if (A.Nv != B.Nv) throw 99;
	}
	catch (int i) {
		printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
#endif // BOUNDS_CHECK
}


inline void check_shape(const Matrix2<int> &A, const Matrix2<int> &B) {
#ifdef BOUNDS_CHECK
	try {
		if (A.Nv != B.Nv) throw 99;
	}
	catch (int i) {
		printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
#endif // BOUNDS_CHECK
}


inline void check_shape(const Matrix2<float> &A, const Matrix2<float> &B) {
#ifdef BOUNDS_CHECK
	try {
		if (A.Nv != B.Nv) throw 99;
	}
	catch (int i) {
		printf("Shape error %d: shape = %d is assigned to %d. \n", i, B.Nv, A.Nv);
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
#endif // BOUNDS_CHECK
}


inline Matrix1<int> step(int i0, int i1) {
	try { if (i0 > i1) throw 99; }
	catch (int i) {
		cerr << "Order error: " << i0 << " > " << i1 << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	Matrix1<int> temp(i1 - i0 + 1);
	for (int i = 0; i <= i1 - i0; ++i) {
		temp.v_p[i] = i + i0;
	}
	return temp;
}


inline void step(int i0, int i1, Matrix1<float> &A) {
	try { if (i0 > i1) throw 99; }
	catch (int i) {
		cerr << "Order error: " << i0 << " > " << i1 << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	for (int i = 0; i <= i1 - i0; ++i) {
		A.v_p[i] = i + i0;
	}
}


inline void step(int i0, int i1, Matrix1<double> &A) {
	try { if (i0 > i1) throw 99; }
	catch (int i) {
		cerr << "Order error: " << i0 << " > " << i1 << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	for (int i = 0; i <= i1 - i0; ++i) {
		A.v_p[i] = i + i0;
	}
}


inline void step(int i0, int i1, Matrix1<int> &A) {
	try { if (i0 > i1) throw 99; }
	catch (int i) {
		cerr << "Order error: " << i0 << " > " << i1 << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	catch (...) {
		cerr << "Exception caught." << endl;
		system("pause"); //exit(EXIT_FAILURE);
	}
	for (int i = 0; i <= i1 - i0; ++i) {
		A.v_p[i] = i + i0;
	}
}


inline Matrix1<float> msin(const Matrix1<float> &A) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = sin(A.v_p[i]); }
	return temp;
}


inline Matrix1<float> mcos(const Matrix1<float> &A) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = cos(A.v_p[i]); }
	return temp;
}


inline Matrix1<float> mtan(const Matrix1<float> &A) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = tan(A.v_p[i]); }
	return temp;
}


inline Matrix1<double> msin(const Matrix1<double> &A) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = sin(A.v_p[i]); }
	return temp;
}


inline Matrix1<double> mcos(const Matrix1<double> &A) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = cos(A.v_p[i]); }
	return temp;
}


inline Matrix1<double> mtan(const Matrix1<double> &A) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = tan(A.v_p[i]); }
	return temp;
}


inline Matrix1<float> atan2(const Matrix1<float> &A, const Matrix1<float> &B) {
	// -PI to +PI
	check_shape(A, B);
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	float m, n;
	for (int i = anv - 1; i >= 0; --i) {
		m = A.v_p[i];
		n = B.v_p[i];
		temp.v_p[i] = _Atan2(m, n);
#ifdef _DEBUG
		if ((temp.v_p[i] > PI) || (temp.v_p[i] < -PI)) {
			cerr << "Wrong atan." << endl;
			system("pause");
		}
		//try {
		//	if (~isdigit(temp.v_p[i])) throw 99;
		//}
		//catch (int i) {
		//	cerr << "NON-digit." << endl;
		//	system("pause");
		//}
#endif // _DEBUG

	}
	return temp;
}


inline Matrix1<double> atan2(const Matrix1<double> &A, const Matrix1<double> &B) {
	// -PI to +PI
	check_shape(A, B);
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	float m, n;
	for (int i = anv - 1; i >= 0; --i) {
		m = A.v_p[i];
		n = B.v_p[i];
		temp.v_p[i] = _Atan2(m, n);
#ifdef _DEBUG
		if ((temp.v_p[i] > PI) || (temp.v_p[i] < -PI)) {
			cerr << "Wrong atan." << endl;
			system("pause");
		}
		//try {
		//	if (~isdigit(temp.v_p[i])) throw 99;
		//}
		//catch (int i) {
		//	cerr << "NON-digit." << endl;
		//	system("pause");
		//}
#endif // _DEBUG

	}
	return temp;
}


inline Matrix1<double> msqrt(const Matrix1<double> &A) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) {
#if _DEBUG
		if (A.v_p[i] < 0) {
			cout << "Negitive number got in sqrt." << endl;
			system("pause");
		}
#endif // _DEBUG
		temp.v_p[i] = sqrt(A.v_p[i]);
	}
	return temp;
}


inline Matrix1<float> msqrt(const Matrix1<float> &A) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) {
#if _DEBUG
		if (A.v_p[i] < 0) {
			cout << "Negitive number got in sqrt." << endl;
			system("pause");
		}
#endif // _DEBUG
		temp.v_p[i] = sqrt(A.v_p[i]);
	}
	return temp;
}


inline Matrix1<double> mmax(const Matrix1<double> &A, double x) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix1<float> mmax(const Matrix1<float> &A, float x) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix1<int> mmax(const Matrix1<int> &A, int x) {
	int anv = A.Nv;
	Matrix1<int> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix1<double> mmin(const Matrix1<double> &A, double x) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix1<float> mmin(const Matrix1<float> &A, float x) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix1<int> mmin(const Matrix1<int> &A, int x) {
	int anv = A.Nv;
	Matrix1<int> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix1<double> mmax(double x, const Matrix1<double> &A) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix1<float> mmax(float x, const Matrix1<float> &A) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix1<int> mmax(int x, const Matrix1<int> &A) {
	int anv = A.Nv;
	Matrix1<int> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix1<double> mmin(double x, const Matrix1<double> &A) {
	int anv = A.Nv;
	Matrix1<double> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix1<float> mmin(float x, const Matrix1<float> &A) {
	int anv = A.Nv;
	Matrix1<float> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix1<int> mmin(int x, const Matrix1<int> &A) {
	int anv = A.Nv;
	Matrix1<int> temp(anv);
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<float> msin(const Matrix2<float> &A) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = sin(A.v_p[i]); }
	return temp;
}


inline Matrix2<float> mcos(const Matrix2<float> &A) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = cos(A.v_p[i]); }
	return temp;
}


inline Matrix2<float> mtan(const Matrix2<float> &A) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = tan(A.v_p[i]); }
	return temp;
}


inline Matrix2<double> msin(const Matrix2<double> &A) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = sin(A.v_p[i]); }
	return temp;
}


inline Matrix2<double> mcos(const Matrix2<double> &A) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = cos(A.v_p[i]); }
	return temp;
}


inline Matrix2<double> mtan(const Matrix2<double> &A) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = tan(A.v_p[i]); }
	return temp;
}


inline Matrix2<float> atan2(const Matrix2<float> &A, const Matrix2<float> &B) {
	// -PI to +PI
	check_shape(A, B);
	Matrix2<float> temp(A.NI, A.NJ);
	float m, n;
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) {
		m = A.v_p[i];
		n = B.v_p[i];
		temp.v_p[i] = _Atan2(m, n);
#ifdef _DEBUG
		if ((temp.v_p[i] > PI) || (temp.v_p[i] < -PI)) {
			cerr << "Wrong atan." << endl;
			system("pause");
		}
		//try {
		//	if (~isdigit(temp.v_p[i])) throw 99;
		//}
		//catch (int i) {
		//	cerr << "NON-digit." << endl;
		//	system("pause");
		//}
#endif // _DEBUG

	}
	return temp;
}


inline Matrix2<double> atan2(const Matrix2<double> &A, const Matrix2<double> &B) {
	// -PI to +PI
	check_shape(A, B);
	Matrix2<double> temp(A.NI, A.NJ);
	float m, n;
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) {
		m = A.v_p[i];
		n = B.v_p[i];
		temp.v_p[i] = _Atan2(m, n);
#ifdef _DEBUG
		if ((temp.v_p[i] > PI) || (temp.v_p[i] < -PI)) {
			cerr << "Wrong atan." << endl;
			system("pause");
		}
		//try {
		//	if (~isdigit(temp.v_p[i])) throw 99;
		//}
		//catch (int i) {
		//	cerr << "NON-digit." << endl;
		//	system("pause");
		//}
#endif // _DEBUG

	}
	return temp;
}


inline Matrix2<double> msqrt(const Matrix2<double> &A) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) {
#if _DEBUG
		if (A.v_p[i] < 0) {
			cout << "Negitive number got in sqrt." << endl;
			system("pause");
		}
#endif // _DEBUG
		temp.v_p[i] = sqrt(A.v_p[i]);
	}
	return temp;
}


inline Matrix2<float> msqrt(const Matrix2<float> &A) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) {
#if _DEBUG
		if (A.v_p[i] < 0) {
			cout << "Negitive number got in sqrt." << endl;
			system("pause");
		}
#endif // _DEBUG
		temp.v_p[i] = sqrt(A.v_p[i]);
	}
	return temp;
}


inline Matrix2<double> mmax(const Matrix2<double> &A, double x) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix2<float> mmax(const Matrix2<float> &A, float x) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix2<int> mmax(const Matrix2<int> &A, int x) {
	Matrix2<int> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix2<double> mmin(const Matrix2<double> &A, double x) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<float> mmin(const Matrix2<float> &A, float x) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<int> mmin(const Matrix2<int> &A, int x) {
	Matrix2<int> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<double> mmax(double x, const Matrix2<double> &A) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix2<float> mmax(float x, const Matrix2<float> &A) {
	Matrix2<float> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix2<int> mmax(int x, const Matrix2<int> &A) {
	Matrix2<int> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Max(A.v_p[i], x); }
	return temp;
}


inline Matrix2<double> mmin(double x, const Matrix2<double> &A) {
	Matrix2<double> temp(A.NI, A.NJ);
	int anv = A.Nv;
	for (int i = anv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<float> mmin(float x, const Matrix2<float> &A) {
	Matrix2<float> temp(A.NI, A.NJ);
	for (int i = A.Nv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<int> mmin(int x, const Matrix2<int> &A) {
	Matrix2<int> temp(A.NI, A.NJ);
	for (int i = A.Nv - 1; i >= 0; --i) { temp.v_p[i] = Min(A.v_p[i], x); }
	return temp;
}


inline Matrix2<double> Arr2Mat(double **B, int m, int n) {
	Matrix2<double> temp(m, n);
	for (int j = m - 1; j >= 0; --j)
	{
		for (int i = n - 1; i >= 0; --i) {
			temp(i, j) = B[i][j];
		}
	}
	return temp;
}


inline Matrix1<double> Arr2Mat(double *B, int n) {
	Matrix1<double> temp(n);
	for (int i = n - 1; i >= 0; --i) { temp.v_p[i] = B[i]; }
	return temp;
}


inline Matrix2<int> Arr2Mat(int **B, int m, int n) {
	Matrix2<int> temp(m, n);
	for (int j = m - 1; j >= 0; --j)
	{
		for (int i = n - 1; i >= 0; --i) {
			temp(i, j) = B[i][j];
		}
	}
	return temp;
}


inline Matrix1<int> Arr2Mat(int *B, int n) {
	Matrix1<int> temp(n);
	for (int i = n - 1; i >= 0; --i) { temp.v_p[i] = B[i]; }
	return temp;
}


inline Matrix2<float> Arr2Mat(float **B, int m, int n) {
	Matrix2<float> temp(m, n);
	for (int j = m - 1; j >= 0; --j)
	{
		for (int i = n - 1; i >= 0; --i) {
			temp(i, j) = B[i][j];
		}
	}
	return temp;
}


inline Matrix1<float> Arr2Mat(float *B, int n) {
	Matrix1<float> temp(n);
	for (int i = n - 1; i >= 0; --i) { temp.v_p[i] = B[i]; }
	return temp;
}

// Matrix computation




//template<class Type> Type glinear1(const Matrix1<double> &A, const Matrix1<double> &B, Type x) {
//	check_shape(A, B);
//	check_bound(A, x);
//	
//}


inline void cross(double A[3], const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2) {
	A[0] = y1*z2 - y2*z1;
	A[1] = z1*x2 - x1*z2;
	A[2] = x1*y2 - y1*x2;
}


inline void cross(float A[3], const float &x1, const float &y1, const float &z1, const float &x2, const float &y2, const float &z2) {
	A[0] = y1*z2 - y2*z1;
	A[1] = z1*x2 - x1*z2;
	A[2] = x1*y2 - y1*x2;
}


inline double norm(const double &x, const double &y, const double &z) {
	return sqrt(x*x + y*y + z*z);
}


inline float norm(const float &x, const float &y, const float &z) {
	return sqrt(x*x + y*y + z*z);
}


inline double dot(const double &x1, const double &y1, const double &z1, const double &x2, const double &y2, const double &z2) {

	return x1*x2 + y1*y2 + z1*z2;
}


inline float dot(const float &x1, const float &y1, const float &z1, const float &x2, const float &y2, const float &z2) {

	return x1*x2 + y1*y2 + z1*z2;
}




#endif


