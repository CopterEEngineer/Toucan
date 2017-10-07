#include "stdafx.h"
#include "Optimization.h"
#include "Components.h"
#include "SimCase.h"

void Optimization::InitOptimization(void) {
	Matrix1<double> _vfs;
	Matrix2<double> _rpms;
	int _numPara1, _numPara2;
	
	numPara = 2, _numPara1 = 13, _numPara2 = 20;
	_vfs.allocate(_numPara1), _rpms.allocate(2, _numPara1);
	Consini.allocate(_numPara1, 6);
	Kwtips.allocate(_numPara1);
	fitness.converge = false;
	fitness.uctrl.allocate(6), fitness.beta.allocate(3);
	fitness.power.allocate(6), fitness.torque.allocate(6);
	fitness.err.allocate(3);
	
	Consini.input("Controls.in");
	Kwtips.input("Kwtips.in");
	_vfs.input("Vfs.in");
	_rpms.input("RPMs.in");
	
	range.emplace_back(_numPara1), range.emplace_back(_numPara2);
	key0.emplace_back(0), key0.emplace_back(0);
	
	for (int i = _numPara1 - 1; i >= 0; i--)
	{
		designs.fspd = _vfs(i);
		for (int j = _numPara2 - 1; j >= 0; j--)
		{
			key0[0] = i, key0[1] = j;
			designs.wme = _rpms(1, i) + (_rpms(0, i) - _rpms(1, i)) / (_numPara2 - 1)*j;
			
			deMap[key0] = designs;
			fitMap[key0] = fitness;
		}
	}
	key1 = key0;
}

void Optimization::InitOptimization(std::map<std::vector<int>, Fitness> &_fitMap, std::map<std::vector<int>, DesignParas> &_deMap)
{
	Matrix1<double> _vfs;
	Matrix2<double> _rpms;
	int _numPara1, _numPara2;

	numPara = 2, _numPara1 = 13, _numPara2 = 21;
	_vfs.allocate(_numPara1), _rpms.allocate(2, _numPara1);
	Consini.allocate(_numPara1, 6);
	Kwtips.allocate(_numPara1);
	fitness.converge = false;
	fitness.uctrl.allocate(6), fitness.beta.allocate(3);
	fitness.power.allocate(6), fitness.torque.allocate(6);
	fitness.err.allocate(3);

	Consini.input("Controls.in");
	Kwtips.input("Kwtips.in");
	_vfs.input("Vfs.in");
	_rpms.input("RPMs.in");

	range.emplace_back(_numPara1), range.emplace_back(_numPara2);
	key0.emplace_back(0), key0.emplace_back(0);

	for (int i = _numPara1 - 1; i >= 0; i--)
	{
		designs.fspd = _vfs(i);
		for (int j = _numPara2 - 1; j >= 0; j--)
		{
			key0[0] = i, key0[1] = j;
			designs.wme = _rpms(1, i) + (_rpms(0, i) - _rpms(1, i)) / (_numPara2 - 1)*j;

			_deMap[key0] = designs, deMap[key0] = designs;
			_fitMap[key0] = fitness, fitMap[key0] = fitness;
		}
	}
	key1 = key0;
}


bool Optimization::haveKey(std::vector<int> k)
{
	bool havekey;
	if (deMap.find(k) != deMap.end() && fitMap.find(k) != fitMap.end())
		return true;
	else
		return false;
}

std::vector<int> Optimization::findNearKey(std::vector<int> k)
{
	std::vector<int> ki(k.size());
	ki[0] = k[0];

	if (k[1] > range[1] >> 1)
	{
		if (k[1] == range[1] - 1)
			ki[1] = k[1] - 1;
		else
			ki[1] = k[1] + 1;
	}
	else
	{
		if (k[1] == 0)
			ki[1] = 1;
		else
			ki[1] = k[1] - 1;
	}
	return ki;
}

void Optimization::findNearKey(std::vector<int> k0, std::vector<int> &kup, std::vector<int> &kdn)
{
	kup[0] = kdn[0] = k0[0];
	if (range[1] == 1)
		kup[1] = kdn[1] = k0[1];
	else
	{
		if (k0[1] == range[1] - 1)
		{
			kup[1] = k0[1];
			kdn[1] = k0[1] - 1;
		}
		else if (k0[1] == 0)
		{
			kdn[1] = 0;
			kup[1] = 1;
		}
		else if (k0[1] >= range[1] >> 1)
		{
			kup[1] = k0[1] + 1;
			kdn[1] = k0[1] - 1;
		}
		else
		{
			kdn[1] = k0[1] - 1;
			kup[1] = k0[1] + 1;
		}
	}
}

void Optimization::findNearKey(std::vector<int> k0, std::vector<int> &kup, std::vector<int> &kdn, std::vector<int> &ki)
{
	ki[0] = kup[0] = kdn[0] = k0[0];
	if (range[1] == 1)
		ki[1] = kup[1] = kdn[1] = k0[1];
	else
	{
		kdn = deMap.upper_bound(k0)->first;
		/*if (k0[1] == range[1] - 1)
		{
			kup[1] = k0[1];
			kdn[1] = k0[1] - 1;
			ki[1] = kdn[1];
		}
		else if (k0[1] == 0)
		{
			kdn[1] = 0;
			kup[1] = 1;
			ki[1] = 1;
		}
		else if (k0[1] >= range[1] >> 1)
		{
			kup[1] = k0[1] + 1;
			kdn[1] = k0[1] - 1;
			ki[1] = kdn[1];
		}
		else
		{
			kdn[1] = k0[1] - 1;
			kup[1] = k0[1] + 1;
			ki[1] = kup[1];
		}*/
	}
}

void Optimization::SetSimCond(Copter &C, std::vector<int> k, std::vector<int> ku, std::vector<int> kd)
{
	double _uctrl[6] = { 0,0,0,0,0,0 };
	//std::vector<int> _k(k.size());

	C.vel_g[0] = deMap[k].fspd;

	if (fitMap[k].converge)
		_setsimcond(C, fitMap[k].uctrl.v_p, k[0]);
	else
	{
		//_k = findNearKey(k);
		if (fitMap[ku].converge)
		{
			_setsimcond(C, fitMap[ku].uctrl.v_p, k[0]);
		}
		else if (fitMap[kd].converge)
		{
			_setsimcond(C, fitMap[kd].uctrl.v_p, k[0]);
		}
		else
		{
			for (int i = 0; i < 6; i++)
				_uctrl[i] = RAD(Consini(k[0], i));
			_setsimcond(C, _uctrl, k[0]);
		}
	}
}

void Optimization::SetSimCond(Copter &C, std::vector<int> k)
{
	double _uctrl[6] = { 0,0,0,0,0,0 };
	std::vector<int> k_temp(numPara);

	for (int i = 0; i < 6; i++)
		_uctrl[i] = RAD(Consini(k[0], i));
	_setsimcond(C, _uctrl, k[0]);

	C.vel_g[0] = deMap[k].fspd;
	/*if (fitMap[k].converge)
		_setsimcond(C, fitMap[k].uctrl.v_p, k[0]);
	else
	{
		auto it = --fitMap.lower_bound(k);
		if (it != fitMap.cend())
			k_temp = it->first;
		else
			k_temp = k;

		if (fitMap[k_temp].converge)
			_setsimcond(C, fitMap[k_temp].uctrl.v_p, k[0]);
		else 
		{
			auto it = fitMap.upper_bound(k);
			if (it != fitMap.cend())
				k_temp = it->first;
			else
				k_temp = k;
			if (fitMap[k_temp].converge)
				_setsimcond(C, fitMap[k_temp].uctrl.v_p, k[0]);
			else
			{
				for (int i = 0; i < 6; i++)
					_uctrl[i] = RAD(Consini(k[0], i));
				_setsimcond(C, _uctrl, k[0]);
			}
		}
	}*/
}

void Optimization::_setsimcond(Copter &C, double *_uctrl, const int ic)
{
	double euler[3] = { 0,0,0 };

	euler[0] = _uctrl[4], euler[1] = _uctrl[5];
	
	C.refcoord.SetCoordinate(euler, "euler");
	C.SetCtrl(_uctrl, C.RotorV[0]);
	C.SetCtrl(_uctrl + 3, C.RotorV[1]);

	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		C.RotorV[i].InitVariables();
		if (C.RotorV[i].adyna > 0)
			C.RotorV[i].WakeModelPrams(Kwtips(ic));
	}
}

void Optimization::UpdateDesigns(std::vector<int> k, std::vector<int> knearest)
{
	deMap[k].wme += deMap[knearest].wme;
	deMap[k].wme *= 0.5;
}

void Optimization::UpdateDesigns(std::vector<int> k)
{
	std::vector<int> k_temp(numPara);
	if (range[1] < 3)
	{
		printf("Parameter separate is not enough in UpdateDesigns() \n");
		return;
	}
	if (deMap[k].wme >= 1.0)
	{
		auto it = fitMap.upper_bound(k);
		if (it != fitMap.cend())
			k_temp = it->first;
		else
			k_temp = k;
	}
	else
	{
		auto it = --fitMap.lower_bound(k);
		if (it != fitMap.cend())
			k_temp = it->first;
		else
			k_temp = k;
	}
	if ((deMap[k].wme - 1.0)*(deMap[k_temp].wme - 1.0) <= 0)
	{
		deMap[k].wme += 1.0;
		deMap[k].wme *= 0.5;
	}
	else 
	{
		deMap[k].wme += deMap[k_temp].wme;
		deMap[k].wme *= 0.5;
	}

}

void Optimization::UpdateParam(Copter &C, DesignParas &dM)
{
	C.RotorV[0].omega *= dM.wme;
	C.RotorV[0].vtipa *= dM.wme;
}

void Optimization::ResetParam(Copter &C, DesignParas &dM)
{
	C.RotorV[0].omega /= dM.wme;
	C.RotorV[0].vtipa /= dM.wme;
}

void Optimization::RecordFitness(Copter &C, CopterSolver &S, Fitness &_ft)
{
	C.GetCtrl(_ft.uctrl.v_p);
	C.RotorV[0].GetBeta(_ft.beta.v_p);
	C.RotorV[0].GetPower(_ft.power.v_p, _ft.torque.v_p);
	C.GetErr(_ft.err.v_p);
	_ft.niter = S.niter;
	_ft.niter_r0 = S.niter_r(0), _ft.niter_r1 = S.niter_r(1);
}

void Optimization::DPMonitor(DesignParas &dM, bool c)
{
	printf("(Vf, Omg, Converge) = (%f, %f, %d) \n", dM.fspd, dM.wme, c);
}

void Optimization::DPMonitor(DesignParas &dM, bool c, const int id)
{
	printf("(threadID, Vf, Omg, Converge) = (%d, %f, %f, %d) \n", id, dM.fspd, dM.wme, c);
}

void Optimization::DPMonitor(DesignParas &dM, bool c, const int id, const int niter)
{
	printf("(threadID, Vf, Omg, Converge, Count) = (%d, %f, %f, %d, %d) \n", id, dM.fspd, dM.wme, c, niter);
}

void Optimization::PostProcess(void)
{
	std::vector<int> _k(numPara);
	Matrix3<double> uctrl3(range[0], 6, range[1]), beta3(range[0], 3, range[1]);
	Matrix3<double> power3(range[0], 6, range[1]), torque3(range[0], 6, range[1]);
	Matrix3<double> err3(range[0], 3, range[1]);
	Matrix2<int> niter2(range[0], range[1]);
	Matrix2<double> flightspeed(range[0], range[1]), wmega(range[0], range[1]);
	for (auto &p : fitMap)
	{
		_k = p.first;
		for (int i = 0; i < 6; i++)
			uctrl3(_k[0], i, _k[1]) = DEG(p.second.uctrl(i));
		
		for (int i = 0; i < 6; i++)
		{
			power3(_k[0], i, _k[1]) = p.second.power(i);
			torque3(_k[0], i, _k[1]) = p.second.torque(i);
		}

		for (int i = 0; i < 3; i++)
			beta3(_k[0], i, _k[1]) = DEG(p.second.beta(i));

		for (int i = 0; i < 3; i++)
			err3(_k[0], i, _k[1]) = p.second.err(i);

		niter2(_k[0], _k[1]) = p.second.niter;
	}
	for (auto &p : deMap)
	{
		_k = p.first;
		flightspeed(_k[0], _k[1]) = p.second.fspd;
		wmega(_k[0], _k[1]) = p.second.wme;
	}

	uctrl3.output2("uctrl_rpm.output", 6);
	niter2.output("niter_rpm.output", 0);
	power3.output2("power_rpm.output", 6);
	torque3.output2("torque_rpm.output", 6);
	beta3.output2("beta_rpm.output", 4);
	err3.output2("err_rpm.output", 10);

	flightspeed.output("flightspeed.output", 5);
	wmega.output("wmega.output", 4);
}

void Optimization::PostProcess(std::map<std::vector<int>, Fitness> _fitMap)
{
	std::vector<int> _k(numPara);
	Matrix3<double> uctrl3(range[0], 6, range[1]);
	for (auto &p : _fitMap)
	{
		_k = p.first;
		for (int i = 0; i < 6; i++)
			uctrl3(_k[0], i, _k[1]) = DEG(p.second.uctrl(i));
	}
	uctrl3.output2("uctrl_rpm.output", 6);
}

void Optimization::RPMSweepMP(void)
{
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e, i, j, k, np, allcase, icount;
	int maxRepeat = 2;
	std::vector<int> _key(2), _keyup(2), _keydn(2);
	std::map<std::vector<int>, Fitness> _fitMap;
	std::map<std::vector<int>, DesignParas> _deMap;

	ul496.GetProb();
	ul496.GetModel();
	copter.InitRotorCraft(ul496);

	//InitOptimization(_fitMap, _deMap);
	InitOptimization();

//	s = i = j = k = icount = 0, e = range[0], np = range[1];
//	e = 1;
//	allcase = np*(e - s);
//
////#pragma omp parallel num_threads(8) shared(s, e, fitMap, deMap) firstprivate(i, j, k, icount, solver, copter)
////#pragma omp parallel num_threads(7) shared(s, e, _fitMap, _deMap) firstprivate(i, j, icount, _key, _keyup, _keydn, solver, copter)
//	{
//		int k, kstart, kend;
//		int Nthrds = omp_get_num_threads(), id = omp_get_thread_num();
//		if (allcase / Nthrds == 0)
//		{
//			kstart = id*(allcase / Nthrds);
//			kend = (id + 1)*(allcase / Nthrds);
//		}
//		else
//		{
//			kstart = id*(allcase / Nthrds + 1);
//			kend = (id + 1)*(allcase / Nthrds + 1);
//		}
//		kend = Min(kend, allcase);
////#pragma omp for
//		for (k = 0; k < allcase; k++)
//		{
//			i = k / np; // flight speed index
//			j = k % np; // parameter index
//			_key[0] = i, _key[1] = j;
//			findNearKey(_key, _keyup, _keydn);
//			if (haveKey(_key))
//			{
//				icount = 0;
//				while (!fitMap[_key].converge && icount < maxRepeat)
//				{
//					SetSimCond(copter, _key, _keyup, _keydn);
//					UpdateParam(copter, deMap[_key]);
//					solver.CopterSimulation(copter);
//					ResetParam(copter, deMap[_key]);
//
//					if (solver.converge)
//					{
//						RecordFitness(copter, solver, fitMap[_key]);
//						fitMap[_key].converge = true;
//					}
//					else
//					{
//						icount++;
//						UpdateDesigns(_key, _key == _keydn ? _keyup : _keydn);
//						fitMap[_key].converge = false;
//					}
//					DPMonitor(deMap[_key], fitMap[_key].converge, id);
//				}
//				if (icount == maxRepeat)
//				{
//					fitMap.erase(_key);
//					deMap.erase(_key);
//				}
//			}
//		}
//#pragma omp barrier
//	}
//	PostProcess(fitMap);
}