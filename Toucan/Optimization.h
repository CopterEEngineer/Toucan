#pragma once
#define OPTIMIZATION_h

#include "MatrixTemplate.h"
#include "Components.h"
#include "Solver.h"
#include <vector>
#include <map>

struct Fitness
{
	bool converge;
	Matrix1<double> uctrl, beta, power, torque, err;
	int niter, niter_r0, niter_r1;
};

struct DesignParas
{
	double fspd, wme;
};

class Optimization
{
public:
	Optimization() {; }
	~Optimization() {; }

	void InitOptimization(void);
	void InitOptimization(std::map<std::vector<int>, Fitness> &fitMap, std::map<std::vector<int>, DesignParas> &deMap);
	bool haveKey(std::vector<int> k);
	void RPMSweepMP(void);
	void SetSimCond(Copter &C, std::vector<int> k, std::vector<int> ku, std::vector<int> kd);
	void SetSimCond(Copter &C, std::vector<int> k);
	void UpdateParam(Copter &C, DesignParas &dM);
	void ResetParam(Copter &C, DesignParas &dM);
	void RecordFitness(Copter &C, CopterSolver &S, Fitness &F);
	void DPMonitor(DesignParas &dM, bool c);
	void DPMonitor(DesignParas &dM, bool c, const int id);
	void DPMonitor(DesignParas &dM, bool c, const int id, const int niter);
	void PostProcess(void);
	void PostProcess(std::map<std::vector<int>, Fitness> _fitMap);

	std::vector<int> findNearKey(std::vector<int> k);
	void findNearKey(std::vector<int> k0, std::vector<int> &kup, std::vector<int> &kdn);
	void findNearKey(std::vector<int> k0, std::vector<int> &kup, std::vector<int> &kdn, std::vector<int> &ki);
	void UpdateDesigns(std::vector<int> k, std::vector<int> knearest);
	void UpdateDesigns(std::vector<int> k);

private:
	void _setsimcond(Copter &C, double *uctrl, const int ic);

public:
	int numPara;
	std::vector<int> range;
	std::vector<int> key0, key1;
	Fitness fitness;
	DesignParas designs;
	std::map<std::vector<int>, Fitness> fitMap;
	std::map<std::vector<int>, DesignParas> deMap;
	Matrix2<double> Consini;
	Matrix1<int> Kwtips;
};

