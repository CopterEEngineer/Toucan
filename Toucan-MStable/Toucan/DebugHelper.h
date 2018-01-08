#pragma once
#ifndef DebugHelper_h
#define DebugHelper_h

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using std::cout;
using std::endl;
using std::cerr;
using std::exception;
using std::string;
using std::ifstream;
using std::ofstream;


inline void print_cons_on_screen(char *s) {
#ifdef _DEBUG
	//cout << s << " constructor." << endl;

#endif // _DEBUG

}

inline void print_cpcons_on_screen(char *s) {
#ifdef _DEBUG
	//cout << s << " copy constructor." << endl;

#endif // _DEBUG

}

inline void print_dscons_on_screen(char *s) {
#ifdef _DEBUG
	//cout << s << " destructor." << endl;

#endif // _DEBUG

}

inline void wrong_coordinate_setting(const char *s) {
#ifdef _DEBUG

	cout << "Wrong coordinate member " << s << " gotten." << endl;
	system("pause");
#endif // _DEBUG

}

inline void wrong_comp_coordinate_diff_base(void) {
#ifdef _DEBUG

	cout << "Wrong coordinate compare with different bases." << endl;
	system("pause");
#endif // _DEBUG

}

template <class Matrix> inline void matrix_output_on_screen(Matrix &A, int pres) {
#ifdef _DEBUG
	A.output(pres);
#endif // _DEBUG

}


inline void print_wrong_msg(char *s) {
#ifdef _DEBUG
	cout << s << endl;
	system("pause");
#endif // _DEBUG

}

#endif // !DebugHelper_h
//#pragma once
