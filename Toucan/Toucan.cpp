// Toucan.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MatrixTemplate.h"
#include "Coordinate.h"
#include "Components.h"
#include "Solver.h"
#include "SimCase.h"

#include <vector>
#include <memory>
#include <time.h>
#include <omp.h>

void levelflight(void)
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

	s = 0, e = jobs.nCase;
	tStart = clock();
	for (int i = s; i < e; ++i)
	{
		cout << i << endl;
		jobs.SetSimCond(copter, i);
		solver.CopterSimulation(copter);
		jobs.PostProcess(copter, i, s, e);
	}
	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
}

void rpmsweep(void)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e, np;
	clock_t tStart;

	ul496.GetProb();
	ul496.GetModel();

	copter.InitRotorCraft(ul496);
	jobs.InitProject();

	jobs.ParamSweep(copter);
	s = 0, e = jobs.nCase, np = jobs.nParams;
	tStart = clock();
	for (int i = s; i < e; ++i)
	{
		//jobs.SetSimCond(copter, i);
		for (int j = 0; j < np; ++j)
		{
			printf("\n\n");
			printf("Case ID: (%d, %d) \n", i, j);
			jobs.SetSimCond(copter, i);
			jobs.UpdateParam(copter, i, j);
			solver.CopterSimulation(copter);
			jobs.PostProcess(copter, i, j, s, e);
		}
	}
	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);
}


int main()
{
	//levelflight();
	//rpmsweep();
	clock_t tStart;

	tStart = clock();

	//LevelFlight();
	LevelFlightMP();

	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);

//	double a[100], b[100], c[100];
//	double *d, *ap, *bp;
//	d = new double[100];
//	ap = a;
//	bp = b;
//	for (int i = 0; i < 100; i++)
//	{
//		a[i] = 1.3*i;
//		b[i] = 0.5 + i;
//		++ap = a + i;
//		++bp = b + i;
//	}
//
//	bp = b; ap = a;
//	//for (int i = 0; i < 100; i++)
//	//	cout << *bp++ << "\t";
//	int i = 0;
//#pragma omp parallel shared(a,b,d) private(i)
//	{
//#pragma omp for
//		for (i = 0; i < 100; i++)
//		{
//			c[i] = a[i] + b[i];
//			*(d+i) = *(a + i) + *(b + i);
//		}
//	}
//
//	for (i = 0; i < 100; i++)
//		cout << c[i] << "\t";
//
//	for (i = 0; i < 100; i++)
//		cout << d[i] << "\t";
	//delete[]d;


	system("pause");
	return 0;
}

