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
	comps.back()->SetCoordBase();
#endif // FLIGHT_TRIM

	tStart = clock();
	for (int ic = 0; ic < CASE_NUM; ++ic) {
		helicopter.SetInit(ic);
		trimsolver.TrimSolver(helicopter, component, comps);
	}
	printf("\n Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);

	printf("\n");
	printf("Completed.\n");
	system("pause");
	return 0;
}

