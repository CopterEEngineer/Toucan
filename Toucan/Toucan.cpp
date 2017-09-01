// Toucan.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MatrixTemplate.h"
#include "Coordinate.h"
//#include "Component.h"
#include "Components.h"
#include "Solver.h"

#include "SimCase.h"
//#include "Solver.h"
#include "UnitTest.h"

#include <vector>
#include <memory>
#include <time.h>
#include <omp.h>


int main()
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e;
	clock_t tStart;

	ul496.GetProb();
	ul496.GetModel();
	copter.InitRotorCraft(ul496);
	
	jobs.InitProject();
	s = 0, e = 2; // jobs.nCase;
	tStart = clock();
	for (int i = s; i < 1; ++i)
	{
		jobs.SetSimCond(copter, i);
		solver.CopterSimulation(copter);
		jobs.PostProcess(copter, i, s, e);
	}
	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);

	system("pause");
	return 0;
}

