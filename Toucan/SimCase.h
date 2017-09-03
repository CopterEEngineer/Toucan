#pragma once
#define SimCase_h

#include "Components.h"
#include "Solver.h"

class Model_UL496 :public ModelCase
{
public:
	Model_UL496() : ModelCase() { ; }
	//Model_UL496(const Model_UL496 &M) :ModelCase(M) { ; }
	~Model_UL496() { ; }

	void GetModel(void);
	void GetProb(void);
	void InitMainRotor(Rotor &R);
	void InitTailRotor(Rotor &R, double w);
};

class Jobs
{
public:
	Jobs() { ; }
	~Jobs() { ; }

	void InitProject(void);
	void SetSimCond(Copter &C, const int ic);
	void PostProcess(Copter &C, const int ic, const int s, const int e);

private:
	Matrix1<double> Mus, Pits;
	Matrix1<int> Kwtips;
public:
	int nCase;
};

