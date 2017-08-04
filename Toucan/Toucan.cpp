// Toucan.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MatrixTemplate.h"
#include "Coordinate.h"
#include "Component.h"
#include "Solver.h"
#include "UnitTest.h"

#include <vector>
#include <memory>
#include <time.h>
#include <omp.h>


int main()
{
	Solver trimsolver;
	Copter helicopter;
	Component component;
	Coordinate *BASE_COORD = &helicopter.refcoord;
	std::vector<std::unique_ptr<Component>> comps;
	clock_t tStart;

	// order fixed
	comps.emplace_back(new Rotor("main"));
	comps.emplace_back(new Fuselage);
	comps.emplace_back(new Wing("hor"));
	comps.emplace_back(new Wing("ver", 1));
	comps.emplace_back(new Wing("ver", 2));
	comps.emplace_back(new Rotor("tail"));

	component.refcoord.SetBase(BASE_COORD);
	for (auto it = comps.cbegin(); it != comps.cend(); ++it) {
		(*it)->refcoord.SetBase(BASE_COORD);
	}
	cout << comps.size() << endl;
#ifdef FLIGHT_TRIM
	comps[0]->SetCoordBase();
#ifndef TROTOR
	comps.back()->SetCoordBase();
#endif // TROTOR

#endif // FLIGHT_TRIM

	trimsolver.TrimSolver(helicopter, component, comps);
	//comps[TROTOR]->SetAirfm();
	//myTYPE jacob[6][6], t12[3][3], temp[3][3];
	//myTYPE a0[6], b[3];
	//Matrix2<myTYPE> jacob_M(6, 6), temp_M(3, 3), t12_M(3, 3);
	//Matrix1<myTYPE> a0_M(6);
	//int icount = 1;

	//temp[0][0] = 1;
	//temp[0][1] = 2;
	//temp[0][2] = 3;
	//temp[1][0] = 4;
	//temp[1][1] = 5;
	//temp[1][2] = 6;
	//temp[2][0] = 8;
	//temp[2][1] = 8;
	//temp[2][2] = 9;

	//t12[0][0] = 30;
	//t12[0][1] = 20;
	//t12[0][2] = 10;
	//t12[1][0] = 99;
	//t12[1][1] = 84;
	//t12[1][2] = 69;
	//t12[2][0] = 162;
	//t12[2][1] = 138;
	//t12[2][2] = 114;

	//b[0] = t12[0][0];
	//b[1] = t12[1][0];
	//b[2] = t12[2][0];

	//for (int i = 0; i < 3; ++i) {
	//	for (int j = 0; j < 3; ++j) {
	//		temp_M(i, j) = temp[i][j];
	//	}
	//}	
	//for (int i = 0; i < 3; ++i) {
	//	for (int j = 0; j < 3; ++j) {
	//		t12_M(i, j) = t12[i][j];
	//	}
	//}
	
	//Msolver(*temp, *t12, 3);
	//Msolver(*temp, b, 3);

	/*Msolver(temp_M.v_p, t12_M.v_p, 3, 3);

	cout << b[0] << "\t" << b[1] << "\t" << b[2] << endl;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << t12[i][j] << "\t";
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << t12_M(i,j) << "\t";
		}
		cout << endl;
	}
	
	
	*/
	
	
	/*jacob_M.input("jacob_temp.txt");
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 6; ++j) {
			jacob[i][j] = jacob_M(i, j);
		}
	}
	a0_M(0) = 3.2234;
	a0_M(1) = 2.8996;
	a0_M(2) = -15.6208;
	a0_M(3) = -0.4372;
	a0_M(4) = -0.2065;
	a0_M(5) = -0.24;
	for (int i = 0; i < 6; ++i) {
		a0[i] = a0_M(i);
	}

	Msolver(*jacob, a0, 6);
	Msolver(jacob_M.v_p, a0_M.v_p, 6);

	t12[2][2] = t12[3][3];
	Msolver(*jacob, *t12, 6, 36);
	for (int i = 0; i < FREEDOM; ++i) {
		for (int j = 0; j < FREEDOM; ++j) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << t12[i][j] << "\t";
		}
		cout << endl;
	}*/
	/*for (int i = 0; i < 6; ++i) {
		cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
		cout << a0[i] << "\t" << a0_M(i) << endl;
	}*/


	//const int size = 256;
	//double sinTable[size];
	//int count = 0;
	//tStart = clock();
	/*for (int n = 0; n<size; ++n)
		sinTable[n] = std::sin(2 * PI * n / size);
	printf("%f s", (double)(clock() - tStart) / CLOCKS_PER_SEC);*/

	printf("\n");
	printf("Completed.\n");
	system("pause");
	return 0;
}

