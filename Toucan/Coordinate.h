#pragma once
#ifndef Coordinate_h
#define Coordinate_h


#include "MatrixTemplate.h"
#include "Algorithm.h"
#include "DebugHelper.h"


class Coordinate {

public:
	myTYPE origin[3], euler[3];
	myTYPE Ttransf[3][3], Etransf[3][3];
	const Coordinate* base;

	Coordinate() {
		//print_cons_on_screen("Coordinate");

		for (int i = 2; i >= 0; --i) {
			origin[i] = 0;
			euler[i] = 0;
			for (int j = 2; j >= 0; --j) {
				if (i == j) {
					Ttransf[i][j] = 1;
					Etransf[i][j] = 1;
				}
				else {
					Ttransf[i][j] = 0;
					Etransf[i][j] = 0;

				}
			}
		}
		base = NULL;
	}

	Coordinate(const Coordinate &A) {
		base = A.base;
		for (int i = 2; i >= 0; --i) {
			origin[i] = A.origin[i];
			euler[i] = A.euler[i];
			for (int j = 2; j >= 0; --j) {
				Ttransf[i][j] = A.Ttransf[i][j];
				Etransf[i][j] = A.Etransf[i][j];
			}
		}
	}

	~Coordinate() {
		print_dscons_on_screen("Coordinate");
		for (int i = 2; i >= 0; --i) {
			origin[i] = 0;
			euler[i] = 0;
			for (int j = 2; j >= 0; --j) {
				if (i == j) {
					Ttransf[i][j] = 1;
					Etransf[i][j] = 1;
				}
				else {
					Ttransf[i][j] = 0;
					Etransf[i][j] = 0;
				}
			}
		}
		base = NULL;
	}

	
	template <class Type, class Coordinate> void SetCoordinate(const Type *og, const Type *el, const Coordinate *C) {
		Type rx[3][3] = { 0 };
		Type ry[3][3] = { 0 };
		Type rz[3][3] = { 0 };

		base = C;

		for (int i = 2; i >= 0; --i) {
			origin[i] = og[i];
			euler[i] = el[i];
		}
		// roll
		rx[0][0] = 1;
		rx[1][1] = cos(euler[0]);
		rx[2][2] = rx[1][1];
		rx[1][2] = sin(euler[0]);
		rx[2][1] = -sin(euler[0]);
		// pitch
		ry[0][0] = cos(euler[1]);
		ry[1][1] = 1;
		ry[2][2] = ry[0][0];
		ry[0][2] = -sin(euler[1]);
		ry[2][0] = sin(euler[1]);
		// yaw
		rz[0][0] = cos(euler[2]);
		rz[1][1] = rz[0][0];
		rz[2][2] = 1;
		rz[0][1] = sin(euler[2]);
		rz[1][0] = -sin(euler[2]);

		// Ttransf
		Type temp[3][3] = { 0 };
#ifdef USE_DOUBLE
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, *rx, 3, *ry, 3, 0, *temp, 3);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, *temp, 3, *rz, 3, 0, *Ttransf, 3);
#else
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, *rx, 3, *ry, 3, 0, *temp, 3);
		cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1, *temp, 3, *rz, 3, 0, *Ttransf, 3);
#endif // USE_DOUBLE

		// Etransf
		Etransf[0][0] = 1;
		Etransf[0][1] = Etransf[1][0] = Etransf[2][0] = 0;
		Etransf[0][2] = -sin(euler[1]);
		Etransf[1][1] = cos(euler[0]);
		Etransf[1][2] = sin(euler[0])*cos(euler[1]);
		Etransf[2][1] = -sin(euler[0]);
		Etransf[2][2] = cos(euler[0])*cos(euler[1]);
	}

	template <class Type> void SetCoordinate(const Type *x, const char *c) {
		if (!strcmp(c, "origin")) { SetCoordinate(x, euler, base); }
		else if (!strcmp(c, "euler")) { SetCoordinate(origin, x, base); }
		else { wrong_coordinate_setting(c); }
	}

	template <class Type, class Coordinate> void Transfer(Type t12[3][3], Type r12[3], const Coordinate &coord1, const Coordinate &coord2) {
		const Type *origin1_ptr = coord1.origin;
		const Type *origin2_ptr = coord2.origin;
		Type temp[3][3];

		if (coord1.base == coord2.base) {
			for (int i = 2; i >= 0; --i) {
				r12[i] = *(origin2_ptr + i) - *(origin1_ptr + i);
				for (int j = 2; j >= 0; --j) {
					t12[i][j] = coord2.Ttransf[i][j];
					temp[i][j] = coord1.Ttransf[i][j];
				}

			}
			//cout << typeid(t12).name() << endl;
			//cout << typeid(temp).name() << endl;

			Msolver(*temp, *t12, 3, 3);
		}
		else { wrong_comp_coordinate_diff_base(); }

	}

	
	
	
	void SetBase(const Coordinate *C) {
		base = C;
	}
};



#endif // !Coordinate_h

