// Toucan.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MatrixTemplate.h"
#include "Coordinate.h"
#include "Components.h"
#include "Solver.h"
#include "SimCase.h"
#include "Optimization.h"

#include <vector>
#include <memory>
#include <time.h>
#include <omp.h>

int main()
{	
	clock_t tStart;
	//Optimization paraSweep;
	LBDynamicStall lbdtest;

	tStart = clock();

	//LevelFlight();
	//LevelFlightMP(1); //54 s, 195s
	//RPMSweepMP(8); // 1043 s
	//RPMSweep(6, 0);
	//OPT_RPMSweep(); // 560s; 100+ s; 1219s; 2514s; 1531s

	LinearModel(1);

	//lbdtest.FuncTest();

	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);

	system("pause");
	return 0;
}

