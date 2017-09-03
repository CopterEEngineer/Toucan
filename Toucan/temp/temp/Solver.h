#pragma once
#ifndef Solver_h
#define Solver_h

// header files
#include "Component.h"
#include <vector>
#include <memory>

// macro definition

class Solver// :private Copter, Component
{
public:
	Solver();

	Solver(const Solver &S);

	~Solver();

	void TrimSolver(Copter &H, Component &CP, const std::vector<std::unique_ptr<Component>> &C);

	void TrimSolver_Ad(Copter &H, Component &CP, const std::vector<std::unique_ptr<Component>> &C);

	void TransientSolver(const Copter &H, const std::vector<Component> &C);
};

#endif // !Solver_h

