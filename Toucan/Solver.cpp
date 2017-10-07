#include "stdafx.h"
#include "Components.h"
#include "Solver.h"


void CopterSolver::CopterSimulation(Copter &C)
{
	InitCopterSolver();
	simtype = C.simtype;
	switch (simtype)
	{
	case FreeTrim0:
	case FreeTrim1:
	case FreeTrim2:
		_FreeTrimSolver(C);
		break;
	case WindTrim:
		_WindTrimSolver(C);
		break;
	case Transient:
		_TransientSolver(C);
		break;
	default:
		break;
	}
}


void CopterSolver::_FreeTrimSolver(Copter &C)
{
	
	double uctrl[6] = { 0.0 };
	double uctrl_temp[6] = { 1.0,1.0,1.0,1.0,1.0,1.0 };
	double _vc[3], _wc[3], _dvc[3], _dwc[3];
	double deltt[6], _dv_p[3], _dv_n[3], _dw_p[3], _dw_n[3];
	int iter = 0;
	bool doWake = false;

	//for (int i = C.RotorV.size() - 1; i >= 0; --i)
	//	if (C.RotorV[i].adyna > 0)
	//		doWake |= true;

	jacobM.allocate(C.nfree, C.nfree);

	C.GetCtrl(uctrl);
#ifdef OUTPUT_MODE
	cout << endl;
	printf("Initial Controls:\n");
	for (int iq = 0; iq < 6; ++iq) {
		cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
		cout << DEG(uctrl[iq]) << "\t";
	}
	cout << endl << endl;
#endif // OUTPUT_MODE

	C.SetStates();
	C.GetStates(_vc, _wc, _dvc, _dwc);

	_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
	_CompsSetAirFM(C);
	_Assemble(C);
	_SetDerivs(dv, dw, SigmaF, SigmaM, C);
	for (int i = 0; i < 3; ++i)
	{
		deltt[i] = dv[i];
		deltt[i + 3] = dw[i];
	}

	for (iter = 0; iter < nitermax; ++iter)
		{
			// tail rotor first
			for (int i = C.RotorV.size() - 1; i >= 0; --i)
			{
				switch (C.RotorV[i].type)
				{
				case Trotor:
					uctrl[3] += epsilon;
					C.SetCtrl((uctrl + 3), C.RotorV[i]);
					_UpdateDerivs(_dv_p, _dw_p, C.RotorV[i], C);

					uctrl[3] -= 2 * epsilon;
					C.SetCtrl((uctrl + 3), C.RotorV[i]);
					_UpdateDerivs(_dv_n, _dw_n, C.RotorV[i], C);

					_SetJacob(_dv_p, _dv_n, _dw_p, _dw_n, 3);
					uctrl[3] += epsilon;
					C.SetCtrl((uctrl + 3), C.RotorV[i]);
					C.RotorV[i].SetAirfm();
					break;
				case Mrotor:
					for (int iq = 0; iq < 3; ++iq)
					{
						uctrl[iq] += epsilon;
						C.SetCtrl(uctrl, C.RotorV[i]);
						_UpdateDerivs(_dv_p, _dw_p, C.RotorV[i], C);

						uctrl[iq] -= 2 * epsilon;
						C.SetCtrl(uctrl, C.RotorV[i]);
						_UpdateDerivs(_dv_n, _dw_n, C.RotorV[i], C);

						_SetJacob(_dv_p, _dv_n, _dw_p, _dw_n, iq);
						uctrl[iq] += epsilon;
						C.SetCtrl(uctrl, C.RotorV[i]);
					}
					break;
				default:
					break;
				}
			}

			// euler
			for (int i = 4; i < 6; ++i)
			{
				uctrl[i] += epsilon;
				_UpdateEuler(uctrl, C);
				C.SetStates(_vc, _wc, _dvc, _dwc);
				_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
				_CompsSetAirFM(C);
				_Assemble(C);
				_SetDerivs(_dv_p, _dw_p, SigmaF, SigmaM, C);

				uctrl[i] -= 2 * epsilon;
				_UpdateEuler(uctrl, C);
				C.SetStates(_vc, _wc, _dvc, _dwc);
				_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
				_CompsSetAirFM(C);
				_Assemble(C);
				_SetDerivs(_dv_n, _dw_n, SigmaF, SigmaM, C);

				_SetJacob(_dv_p, _dv_n, _dw_p, _dw_n, i);
				uctrl[i] += epsilon;
			}
			_UpdateEuler(uctrl, C);

#ifdef OUTPUT_MODE
			cout << endl;
			printf("Accelaration:\n");
			for (int iq = 0; iq < 6; ++iq) {
				cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
				cout << deltt[iq] << "\t";
			}
			cout << endl << endl;

			printf("Jacob:\n");
			for (int i = 0; i < 6; ++i) {
				for (int j = 0; j < 6; ++j) {
					cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
					cout << jacobM(i, j) << "\t";
				}
				cout << endl;
			}
#endif // !OUTPUT_MODE_1

			Msolver(jacobM.v_p, deltt, 6);
			if (_UpdateCtrls(uctrl, deltt) == 6)
			{
				printf("Controls reach boundary.");
				break;
			}

			_UpdateEuler(uctrl, C);
			_UpdateRotorCtrls(uctrl, C);

			C.SetStates();
			C.GetStates(_vc, _wc, _dvc, _dwc);

			_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
			//if (iter > 0 && Max(sum_a1_del(iter - 1), sum_a2_del(iter - 1)) < err_a*err_a * 4)
			//	doWake = true;
			//if (doWake)
			//	_EnableWake(C);

			_EnableWake(C);
			_CompsSetAirFM(C);
			_Assemble(C);
			_SetDerivs(dv, dw, SigmaF, SigmaM, C);
			for (int i = 0; i < 3; ++i)
			{
				deltt[i] = dv[i];
				deltt[i + 3] = dw[i];
			}

			if (isExit(uctrl, deltt, iter))
				break;
		}


	for (int i = 0; i < C.nfree; ++i)
		C.controls[i] = uctrl[i];
	niter = iter == nitermax ? iter - 1 : iter;
	C.sum_a1_del = sum_a1_del(niter);
	C.sum_a2_del = sum_a2_del(niter);
	C.max_c_del = max_c_del(niter);
	C.Niter = niter;
	if (Max(C.sum_a1_del, C.sum_a2_del) < 100 * err_a*err_a && C.max_c_del < 10 * err_c)
		converge = true;
	else
		converge = false;

}

void CopterSolver::_WindTrimSolver(Copter &C)
{
	;
}

void CopterSolver::_TransientSolver(Copter &C)
{
	;
}

void CopterSolver::_CompsSetAirFM(Copter &C)
{
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		C.RotorV[i].SetAirfm();
		niter_r(i) = C.RotorV[i].niter_w;
	}
	for (int i = C.WingV.size() - 1; i >= 0; --i)
		C.WingV[i].SetAirfm();
	C.fuselage.SetAirfm();

#ifdef OUTPUT_MODE
	double ftemp[3], mtemp[3];
	printf("Component aerodynamics.\n");
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		C.RotorV[i].GetAirfm(ftemp, mtemp);
		printf("Rotor Comps: %d \n", C.RotorV[i].type);
		printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
		printf("M: %f, %f, %f \n\n", mtemp[0], mtemp[1], mtemp[2]);
	}
	for (int i = C.WingV.size() - 1; i >= 0; --i)
	{
		C.WingV[i].GetAirfm(ftemp, mtemp);
		printf("Rotor Comps: %d \n", C.WingV[i].type);
		printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
		printf("M: %f, %f, %f \n\n", mtemp[0], mtemp[1], mtemp[2]);
	}	
	C.fuselage.GetAirfm(ftemp, mtemp);
	printf("Fuselage \n");
	printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
	printf("M: %f, %f, %f \n\n", mtemp[0], mtemp[1], mtemp[2]);
#endif // OUTPUT_MODE
}

void CopterSolver::_Assemble(Copter &C)
{
	myTYPE ftemp[3] = { 0,0,0 };
	myTYPE mtemp[3] = { 0,0,0 };
	SigmaF[0] = SigmaF[1] = SigmaF[2] = 0;
	SigmaM[0] = SigmaM[1] = SigmaM[2] = 0;

	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		C.RotorV[i].SetAirfm_cg(C.refcoord.base);
		C.RotorV[i].GetAirfm_cg(ftemp, mtemp);

		for (int j = 0; j < 3; ++j)
		{
			SigmaF[j] += ftemp[j];
			SigmaM[j] += mtemp[j];
		}
	}
		
	for (int i = C.WingV.size() - 1; i >= 0; --i)
	{
		C.WingV[i].SetAirfm_cg(C.refcoord.base);
		C.WingV[i].GetAirfm_cg(ftemp, mtemp);

		for (int j = 0; j < 3; ++j)
		{
			SigmaF[j] += ftemp[j];
			SigmaM[j] += mtemp[j];
		}
	}

	C.fuselage.SetAirfm_cg(C.refcoord.base);
	C.fuselage.GetAirfm_cg(ftemp, mtemp);
	for (int j = 0; j < 3; ++j)
	{
		SigmaF[j] += ftemp[j];
		SigmaM[j] += mtemp[j];
	}
}

void CopterSolver::_EnableWake(Copter &C)
{
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
		if (C.RotorV[i].adyna > 0)
			C.RotorV[i].WakeInducedVel();
}