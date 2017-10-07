#pragma once
#define SimCase_h

#include "Components.h"
#include "Solver.h"
#include "Optimization.h"


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
	void InitProject(Jobs &J, const int ic);
	void SetSimCond(Copter &C, const int ic);
	void SetSimCond(Copter &C);
	void PostProcess(Copter &C);
	void PostProcess(Copter &C, const int ic, const int s, const int e);
	void PostProcess(Copter &C, const int ic, const int ip, const int s, const int e);
	void PostProcessMP(Copter &C, const int ic, const int s, const int e);
	void PostProcessMP(Copter &C, const int ic, const int ip, const int s, const int e, const int np);

	void ParamSweep(const Copter &C);
	void UpdateParam(Copter &C, const int ic, const int ip);

	//void LevelFlight(void);
	//void RPMSweep(void);

private:
	Matrix1<double> Mus, Pits, Vfs;
	Matrix2<double> Consini;
	Matrix1<int> Kwtips;
	Matrix2<double> RPMs, Chds, Rads;
	Matrix1<double> param0;
	Matrix2<double> uctrl, beta;
	Matrix2<double> _power, _torque;
	Matrix3<double> uctrl3, beta3, _power3, _torque3, err3;
	Matrix2<double> flightspeed2, wmega2;
	int flg;
public:
	int nCase, nParams;
	JobsType jtype;
	Matrix1<JobsType> jtype_M;
};

void LevelFlight(void);
void LevelFlightMP(void);
void RPMSweepMP(void);
void ChordSweepMP(void);
void MultiSweepMP(void);
void RPMSweep(const int ic, const int ip);
void OPT_RPMSweep(void);

