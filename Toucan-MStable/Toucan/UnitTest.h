#pragma once
#ifndef UnitTest_h
#define UnitTest_h

#include "stdafx.h"
#include "MatrixTemplate.h"
#include "Coordinate.h"


void test_random(myTYPE &x);

void test_random(Matrix1<myTYPE> &A);

void test_random(Matrix2<myTYPE> &A);

void test_random(myTYPE &x);

void test_findmaxmin(Matrix1<myTYPE> &A);

void test_findmaxmin(Matrix2<myTYPE> &A);

void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, myTYPE xs);

void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, Matrix1<myTYPE> &XN);

void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, Matrix2<myTYPE> &Z, myTYPE xs, myTYPE ys);

void test_interplinear(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y, Matrix2<myTYPE> &Z, Matrix2<myTYPE> XN, Matrix2<myTYPE> YN);

void test_matrixmultip(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y);

void test_matrixmultip(Matrix2<myTYPE> &A, Matrix1<myTYPE> &b);

void test_matrixmultip(Matrix2<myTYPE> &X, Matrix2<myTYPE> &Y);

void test_vectormultip(Matrix1<myTYPE> &X, Matrix1<myTYPE> &Y);

void test_vectormultip(Matrix2<myTYPE> &X, Matrix2<myTYPE> &Y);

void test_coordinate(Coordinate &X, Coordinate &Y);

#endif // !UnitTest_h
