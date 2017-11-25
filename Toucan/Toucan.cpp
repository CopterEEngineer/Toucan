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
	Optimization paraSweep;

	tStart = clock();

	//LevelFlight();
	//LevelFlightMP(); //54 s
	//RPMSweepMP(); // 1043 s
	//rpmsweep(); // 3483 s
	//RPMSweep(6, 0);
	OPT_RPMSweep(); // 560s; 100+ s; 1219s; 2514s; 1531s

	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);

	system("pause");
	return 0;
}

