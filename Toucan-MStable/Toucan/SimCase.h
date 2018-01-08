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
	void GetProb(double);
	void GetProb(int, double, SimType);
	void InitMainRotor(Rotor &R);
	void InitTailRotor(Rotor &R, double w);
	void ParamSweep(void);

public:
	Matrix2<double> RPMs;
};

class Model_BO105 :public ModelCase
{
public:
	Model_BO105() :ModelCase() { ; }
	~Model_BO105() { ; }

	void GetModel(void);
	//void GetProb(void);
	//void GetProb(double);
	void GetProb(int, double, SimType);
	void InitFuselage(void);
	void InitWing(Wing &W);
	void InitFin(Wing &W);
	void InitMainRotor(Rotor &R);
	void InitTailRotor(Rotor &R, double w);

};

class Model_Puma330 :public ModelCase
{
public:
	Model_Puma330() :ModelCase() { ; }
	~Model_Puma330() { ; }

	void GetModel(void);
	//void GetProb(int, double, SimType);
	void GetProb(int, double, SimType, bool);
	void InitFuselage(void);
	void InitWing(Wing &W);
	void InitFin(Wing &W);
	void InitMainRotor(Rotor &R);
	void InitTailRotor(Rotor &R, double w);
};

class Jobs
{
public:
	Jobs() { ; }
	~Jobs() { ; }

	void InitProject(JobsType);
	void InitProject(JobsType, string);
	void InitProject(Jobs &J, const int ic);
	void SetSimCond(Copter &C, const int ic);
	void SetSimCond(Copter &C);
	void PostProcess(Copter &C);
	void PostProcess(Copter &, string);
	void PostProcess(Copter &, string, const int, const int, const int);
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
	Matrix2<double> uctrl, beta, err;
	Matrix2<double> _power, _torque;
	Matrix3<double> uctrl3, beta3, _power3, _torque3, err3;

	Matrix2<double> uctrl_tr, beta_tr;
	Matrix2<double> _power_tr, _torque_tr;
	Matrix3<double> uctrl3_tr, beta3_tr, _power3_tr, _torque3_tr, err3_tr;

	Matrix3<double> sol3_mr, azmth_mr;

	Matrix2<int> niter;

	Matrix2<double> flightspeed2, wmega2;
	Matrix2<double> dXdvel, dYdvel, dZdvel, dXdomg, dYdomg, dZdomg, dXdeul, dYdeul, dZdeul;
	Matrix2<double> dMdvel, dNdvel, dLdvel, dMdomg, dNdomg, dLdomg, dMdeul, dNdeul, dLdeul;
	Matrix2<double> dXdctrl, dYdctrl, dZdctrl, dMdctrl, dNdctrl, dLdctrl;
	Matrix3<double> AMatrix, BMatrix, AMatrixRe;
	Matrix3<double> ALongM, ALateM, ALongMRe, ALateMRe;

	Matrix2<double> StableInitPara;

	int flg;
	string path;
public:
	int nCase, nParams;
	JobsType jtype;
	Matrix1<JobsType> jtype_M;
};

void LevelFlight(void);
void LevelFlightMP(int);
void LinearModel(int);
void RPMSweepMP(int);
void ChordSweepMP(void);
void MultiSweepMP(void);
void RPMSweep(const int ic, const int ip);
void OPT_RPMSweep(void);

