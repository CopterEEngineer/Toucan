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
	void ParamSweep(void);

public:
	Matrix2<double> RPMs;
};

class Jobs
{
public:
	Jobs() { ; }
	~Jobs() { ; }

	void InitProject(void);
	void SetSimCond(Copter &C, const int ic);
	void PostProcess(Copter &C, const int ic, const int s, const int e);
	void PostProcess(Copter &C, const int ic, const int ip, const int s, const int e);

	void ParamSweep(const Copter &C);
	void UpdateParam(Copter &C, const int ic, const int ip);

private:
	Matrix1<double> Mus, Pits, Vfs;
	Matrix1<int> Kwtips;
	Matrix2<double> RPMs;
	Matrix1<double> param0;
public:
	int nCase, nParams;
};

