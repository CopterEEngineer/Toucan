#include "stdafx.h"
#include "Components.h"
#include "Solver.h"
#include <vector>

void CopterSolver::CopterSimulation(Copter &C)
{
	InitCopterSolver();
	simtype = C.simtype;
	
	std::vector<AeroDynType> ad;
	for (auto it = C.RotorV.begin(); it != C.RotorV.end(); ++it)
		ad.emplace_back(it->adyna);

	if (C.RotorV[0].adyna>100 || C.RotorV[1].adyna>100)
	{
		C.RotorV[0].adyna = Averaged;
		C.RotorV[1].adyna = Min(Averaged, ad[1]);
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
			//_TransientSolver(C);
			break;
		default:
			break;
		}
	}

	C.RotorV[0].adyna = ad[0];
	C.RotorV[1].adyna = ad[1];
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
		//C.simtype = FreeTrim1;
		//_FreeTrimSolver(C);
		//C.simtype = Transient;
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
	if (Max(C.sum_a1_del, C.sum_a2_del) < 1 * err_a*err_a && C.max_c_del < 1 * err_c)
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
	double _vc[3], _wc[3], _dvc[3], _dwc[3];
	double ftempp[3], mtempp[3], ftempn[3], mtempn[3];
	double v2;
	double uctrl[6] = { 0.0 };

	// setting trim controls, compute trim forces and moments
	C.SetStates();
	C.GetStates(_vc, _wc, _dvc, _dwc); // body coordinate states
	_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
	C.RotorV[0].AvrgInducedVel(); // compute cirlb
	_EnableWake(C);
	_CompsSetAirFM(C);

	// check if trimmed
	_Assemble(C);
	_SetDerivs(dv, dw, SigmaF, SigmaM, C);
	for (int i = 0; i < 3; i++)
		cout << dv[i] << " ";
	cout << endl;
	for (int i = 0; i < 3; i++)
		cout << dw[i] << " ";
	cout << endl;

#ifdef _DEBUG
	//system("pause");
#endif // _DEBUG

	// velocity perturbance
	for (int iq = 0; iq < 3; iq++)
	{
		double vtemp = _vc[iq];
		double delta;

		if (C.fuselage.si_unit)
			delta = 1; // m/s
		else
			delta = 0.3048; // ft/s

		_vc[iq] = vtemp + delta; // vtemp*(1 + delta);
		_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		//C.RotorV[0].DiskOutput("DiskVelP");
		_Assemble(C);
		_GetSigmaFM(C, ftempp, mtempp);

		_vc[iq] = vtemp - delta; // vtemp*(1 - delta);
		_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		//C.RotorV[0].DiskOutput("DiskVelN");
		_Assemble(C);
		_GetSigmaFM(C, ftempn, mtempn);

		if (C.fuselage.si_unit)
		{
			dXdvel(iq) = (ftempp[0] - ftempn[0]) / 2 / delta / C.mass; // unit s-1
			dYdvel(iq) = (ftempp[1] - ftempn[1]) / 2 / delta / C.mass;
			dZdvel(iq) = (ftempp[2] - ftempn[2]) / 2 / delta / C.mass;
			dLdvel(iq) = (mtempp[0] - mtempn[0]) / 2 / delta;
			dMdvel(iq) = (mtempp[1] - mtempn[1]) / 2 / delta; // unit kg*m/s
			dNdvel(iq) = (mtempp[2] - mtempn[2]) / 2 / delta;
		}
		else
		{
			dXdvel(iq) = (ftempp[0] - ftempn[0]) / 2 / delta / (C.mass / SLUG_CONST); // unit s-1
			dYdvel(iq) = (ftempp[1] - ftempn[1]) / 2 / delta / (C.mass / SLUG_CONST);
			dZdvel(iq) = (ftempp[2] - ftempn[2]) / 2 / delta / (C.mass / SLUG_CONST);
			dLdvel(iq) = 0.3048 * (mtempp[0] - mtempn[0]) / 2 / delta;
			dMdvel(iq) = 0.3048 * (mtempp[1] - mtempn[1]) / 2 / delta; // unit slug*m/s
			dNdvel(iq) = 0.3048 * (mtempp[2] - mtempn[2]) / 2 / delta;
		}
		_vc[iq] = vtemp;
	}
	// omega perturbance
	for (int iq = 0; iq < 3; iq++)
	{
		double otemp = _wc[iq];
		double delta = 0.1; // rad/s
		_wc[iq] = otemp + delta; // otemp*(1 + delta);
		_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		_Assemble(C);
		_GetSigmaFM(C, ftempp, mtempp);

		_wc[iq] = otemp - delta; // otemp*(1 - delta);
		_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		_Assemble(C);
		_GetSigmaFM(C, ftempn, mtempn);

		if (C.fuselage.si_unit)
		{
			dXdomg(iq) = (ftempp[0] - ftempn[0]) / 2 / delta / C.mass; // m/s
			dYdomg(iq) = (ftempp[1] - ftempn[1]) / 2 / delta / C.mass;
			dZdomg(iq) = (ftempp[2] - ftempn[2]) / 2 / delta / C.mass;
			dLdomg(iq) = (mtempp[0] - mtempn[0]) / 2 / delta; // kg*m2/s
			dMdomg(iq) = (mtempp[1] - mtempn[1]) / 2 / delta;
			dNdomg(iq) = (mtempp[2] - mtempn[2]) / 2 / delta;
		}
		else
		{
			dXdomg(iq) = 0.3048 * (ftempp[0] - ftempn[0]) / 2 / delta / (C.mass / SLUG_CONST); // m/s
			dYdomg(iq) = 0.3048 * (ftempp[1] - ftempn[1]) / 2 / delta / (C.mass / SLUG_CONST);
			dZdomg(iq) = 0.3048 * (ftempp[2] - ftempn[2]) / 2 / delta / (C.mass / SLUG_CONST);
			dLdomg(iq) = 0.3048 * 0.3048 * (mtempp[0] - mtempn[0]) / 2 / delta; // slug*m2/s
			dMdomg(iq) = 0.3048 * 0.3048 * (mtempp[1] - mtempn[1]) / 2 / delta;
			dNdomg(iq) = 0.3048 * 0.3048 * (mtempp[2] - mtempn[2]) / 2 / delta;
		}
		_wc[iq] = otemp;
	}
	// euler perturbance
	C.GetCtrl(uctrl);
	for (int i = 4; i < 6; i++)
	{
		double etemp = uctrl[i];
		double delta = 0.1; // rad
		uctrl[i] = etemp + delta; // etemp*(1 + delta);
		_UpdateEuler(uctrl, C);
		C.SetStates(_vc, _wc, _dvc, _dwc);
		_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		_Assemble(C);
		_GetSigmaFM(C, ftempp, mtempp);

		uctrl[i] = etemp - delta; // etemp*(1 - delta);
		_UpdateEuler(uctrl, C);
		C.SetStates(_vc, _wc, _dvc, _dwc);
		_CompsSetStates(C, _vc, _wc, _dvc, _dwc);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		_Assemble(C);
		_GetSigmaFM(C, ftempn, mtempn);

		if (C.fuselage.si_unit)
		{
			dXdeul(i - 4) = (ftempp[0] - ftempn[0]) / 2 / delta / C.mass; // unit m/s2
			dYdeul(i - 4) = (ftempp[1] - ftempn[1]) / 2 / delta / C.mass;
			dZdeul(i - 4) = (ftempp[2] - ftempn[2]) / 2 / delta / C.mass;
			dLdeul(i - 4) = (mtempp[0] - mtempn[0]) / 2 / delta; // kg*m2/s2
			dMdeul(i - 4) = (mtempp[1] - mtempn[1]) / 2 / delta;
			dNdeul(i - 4) = (mtempp[2] - mtempn[2]) / 2 / delta;
		}
		else
		{
			dXdeul(i - 4) = 0.3048 * (ftempp[0] - ftempn[0]) / 2 / delta / (C.mass / SLUG_CONST); // unit m/s2
			dYdeul(i - 4) = 0.3048 * (ftempp[1] - ftempn[1]) / 2 / delta / (C.mass / SLUG_CONST);
			dZdeul(i - 4) = 0.3048 * (ftempp[2] - ftempn[2]) / 2 / delta / (C.mass / SLUG_CONST);
			dLdeul(i - 4) = 0.3048 * 0.3048 * (mtempp[0] - mtempn[0]) / 2 / delta; // slug*m2/s2
			dMdeul(i - 4) = 0.3048 * 0.3048 * (mtempp[1] - mtempn[1]) / 2 / delta;
			dNdeul(i - 4) = 0.3048 * 0.3048 * (mtempp[2] - mtempn[2]) / 2 / delta;
		}

		uctrl[i] = etemp;
	}
	_UpdateEuler(uctrl, C);
	C.SetStates(_vc, _wc, _dvc, _dwc);
	_CompsSetStates(C, _vc, _wc, _dvc, _dwc);

	// controls perturbance
	//double _delta[4] = { 0.1,0.001,0.001,0.1 };
	for (int i = 0; i < 4; i++)
	{
		double ctemp = uctrl[i];
		double delta = 0.1;// _delta[i]; // rad
		uctrl[i] = ctemp + delta; // ctemp*(1 + delta);
		_UpdateRotorCtrls(uctrl, C);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		//C.RotorV[0].DiskOutput("DiskCtrlP");
		_Assemble(C);
		_GetSigmaFM(C, ftempp, mtempp);

		uctrl[i] = ctemp - delta; // ctemp*(1 - delta);
		_UpdateRotorCtrls(uctrl, C);
		//C.RotorV[0].AvrgInducedVel(); // should be dynamic inflow model
		_EnableWake(C);
		_CompsSetAirFM(C);
		//C.RotorV[0].DiskOutput("DiskCtrlN");
		_Assemble(C);
		_GetSigmaFM(C, ftempn, mtempn);

		if (C.fuselage.si_unit)
		{
			dXdctrl(i) = (ftempp[0] - ftempn[0]) / 2 / delta / C.mass;// unit m/s2
			dYdctrl(i) = (ftempp[1] - ftempn[1]) / 2 / delta / C.mass;
			dZdctrl(i) = (ftempp[2] - ftempn[2]) / 2 / delta / C.mass;
			dLdctrl(i) = (mtempp[0] - mtempn[0]) / 2 / delta; // kg*m2/s2
			dMdctrl(i) = (mtempp[1] - mtempn[1]) / 2 / delta;
			dNdctrl(i) = (mtempp[2] - mtempn[2]) / 2 / delta;
		}
		else
		{
			dXdctrl(i) = 0.3048 * (ftempp[0] - ftempn[0]) / 2 / delta / (C.mass / SLUG_CONST);// unit m/s2
			dYdctrl(i) = 0.3048 * (ftempp[1] - ftempn[1]) / 2 / delta / (C.mass / SLUG_CONST);
			dZdctrl(i) = 0.3048 * (ftempp[2] - ftempn[2]) / 2 / delta / (C.mass / SLUG_CONST);
			dLdctrl(i) = 0.3048 * 0.3048 * (mtempp[0] - mtempn[0]) / 2 / delta; // slug*m2/s2
			dMdctrl(i) = 0.3048 * 0.3048 * (mtempp[1] - mtempn[1]) / 2 / delta;
			dNdctrl(i) = 0.3048 * 0.3048 * (mtempp[2] - mtempn[2]) / 2 / delta;
		}

		uctrl[i] = ctemp;
	}
	_UpdateRotorCtrls(uctrl, C);
	for (int i = 0; i < C.nfree; ++i)
		C.controls[i] = uctrl[i];

	// save derivates to Copter
	//C.dXdvel = dXdvel, C.dYdvel = dYdvel, C.dZdvel = dZdvel;
	//C.dXdomg = dXdomg, C.dYdomg = dYdomg, C.dZdomg = dZdomg;
	//C.dXdeul = dXdeul, C.dYdeul = dYdeul, C.dZdeul = dZdeul;

	//C.dLdvel = dLdvel, C.dMdvel = dMdvel, C.dNdvel = dNdvel;
	//C.dLdomg = dLdomg, C.dMdomg = dMdomg, C.dNdomg = dNdomg;
	//C.dLdeul = dLdeul, C.dMdeul = dMdeul, C.dNdeul = dNdeul;

	//C.dXdctrl = dXdctrl, C.dYdctrl = dYdctrl, C.dZdctrl = dZdctrl;
	//C.dLdctrl = dLdctrl, C.dMdctrl = dMdctrl, C.dNdctrl = dNdctrl;

	// compute A matrix
	_ComputeAMatrix(AMatrix, C);
	C.AMatrix = AMatrix;

	// compute B matrix
	_ComputeBMatrix(BMatrix, C);
	C.BMatrix = BMatrix;
}

void CopterSolver::_CompsSetAirFM(Copter &C)
{
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		C.RotorV[i].SetAirfm();
		niter_r(i) = C.RotorV[i].niter_w;
	}
	for (int i = C.WingV.size() - 1; i >= 0; --i)
	{
		switch (C.WingV[i].type)
		{
		case Hwing:
			C.WingV[i].SetAirfm(0.0);
			break;
		case Vwing:
			C.WingV[i].SetAirfm(0.0);
			break;
		default:
			break;
		}
	}
		
	C.fuselage.SetAirfm(C.RotorV[0].GetLambdi());

#ifdef OUTPUT_MODE
	double ftemp[3], mtemp[3];
	//printf("Component aerodynamics.\n");
	//for (int i = C.RotorV.size() - 1; i >= 0; --i)
	//{
	//	C.RotorV[i].GetAirfm(ftemp, mtemp);
	//	printf("Rotor Comps: %d \n", C.RotorV[i].type);
	//	printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
	//	printf("M: %f, %f, %f \n\n", mtemp[0], mtemp[1], mtemp[2]);
	//}
	//for (int i = C.WingV.size() - 1; i >= 0; --i)
	//{
	//	C.WingV[i].GetAirfm(ftemp, mtemp);
	//	printf("WingV Comps: %d \n", C.WingV[i].type);
	//	printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
	//	printf("M: %f, %f, %f \n\n", mtemp[0], mtemp[1], mtemp[2]);
	//}	
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
	
	//for (int i = C.WingV.size() - 1; i >= 0; --i)
	//{
	//	C.WingV[i].SetAirfm_cg(C.refcoord.base);
	//	C.WingV[i].GetAirfm_cg(ftemp, mtemp);

	//	for (int j = 0; j < 3; ++j)
	//	{
	//		SigmaF[j] += ftemp[j];
	//		SigmaM[j] += mtemp[j];
	//	}
	//}

	//C.fuselage.SetAirfm_cg(C.refcoord.base);
	//C.fuselage.GetAirfm_cg(ftemp, mtemp);
	//for (int j = 0; j < 3; ++j)
	//{
	//	SigmaF[j] += ftemp[j];
	//	SigmaM[j] += mtemp[j];
	//}

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
}

void CopterSolver::_EnableWake(Copter &C)
{
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
		if (C.RotorV[i].adyna > 100)
			C.RotorV[i].WakeInducedVel();
}

void CopterSolver::_GetSigmaFM(Copter & C, myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = SigmaF[i];
		m[i] = SigmaM[i];
	}
}

void CopterSolver::_ComputeAMatrix(Matrix2<double> &A, Copter &C)
{

	double uctrl[6] = { 0.0 };
	double Gg = 9.8;
	double Muu, Mvv, Mww, Mpp, Mqq, Mrr, Lpp, Lqq, Lrr, Nrr, Luu, Lww, Lvv, Nuu, Nww, Nqq, Nvv, Npp;
	double Lsita, Lfi, Msita, Mfi, Nsita, Nfi;
	double f1, f2, f3;
	double vector1, vector2, vector3, vector4, vector5;
	double unit_factor = 1.0;

	A.allocate(xdof, xdof);
	C.GetCtrl(uctrl);

	if (!C.fuselage.si_unit)
		unit_factor = 0.3048;

	vector1 = 1 / unit_factor / unit_factor * C.inmatx_M(2, 2) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector2 = 1 / unit_factor / unit_factor * C.inmatx_M(0, 2) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector3 = 1 / unit_factor / unit_factor * C.inmatx_M(0, 2) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector4 = 1 / unit_factor / unit_factor * C.inmatx_M(0, 0) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));

	Luu = dot(vector1, vector2, 0, dLdvel(0), dNdomg(0), 0);
	Lvv = dot(vector1, vector2, 0, dLdvel(1), dNdomg(1), 0);
	Lww = dot(vector1, vector2, 0, dLdvel(2), dNdomg(2), 0);
	Lpp = dot(vector1, vector2, 0, dLdomg(0), dNdomg(0), 0);
	Lqq = dot(vector1, vector2, 0, dLdomg(1), dNdomg(1), 0);
	Lrr = dot(vector1, vector2, 0, dLdomg(2), dNdomg[2], 0);
	Lfi = dot(vector1, vector2, 0, dLdeul(0), dNdeul(1), 0);
	Lsita = dot(vector1, vector2, 0, dLdeul(1), dNdeul(1), 0);

	Nuu = dot(vector3, vector4, 0, dLdvel[0], dNdomg[0], 0);
	Nvv = dot(vector3, vector4, 0, dLdvel[1], dNdomg[1], 0);
	Nww = dot(vector3, vector4, 0, dLdvel[2], dNdomg[2], 0);
	Npp = dot(vector3, vector4, 0, dLdomg[0], dNdomg[0], 0);
	Nqq = dot(vector3, vector4, 0, dLdomg[1], dNdomg[1], 0);
	Nrr = dot(vector3, vector4, 0, dLdomg[2], dNdomg[2], 0);
	Nfi = dot(vector3, vector4, 0, dLdeul(0), dNdeul(0), 0);
	Nsita = dot(vector3, vector4, 0, dLdeul(1), dNdeul(1), 0);

	vector5 = (unit_factor * unit_factor * C.inmatx_M(1, 1));
	Muu = dMdvel[0] / vector5;
	Mvv = dMdvel[1] / vector5;
	Mww = dMdvel[2] / vector5;
	Mpp = dMdomg[0] / vector5;
	Mqq = dMdomg[1] / vector5;
	Mrr = dMdomg[2] / vector5;
	Mfi = dMdeul[0] / vector5;
	Msita = dMdeul[1] / vector5;

	f1 = C.inmatx_M(0, 2)*(C.inmatx_M(2, 2) + C.inmatx_M(0, 0) - C.inmatx_M(1, 1)) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	f2 = (C.inmatx_M(2, 2)*(C.inmatx_M(2, 2) - C.inmatx_M(1, 1)) + C.inmatx_M(0, 2)*C.inmatx_M(0, 2)) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	f3 = (C.inmatx_M(0, 0)*(C.inmatx_M(1, 1) - C.inmatx_M(0, 0)) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2)) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));


	// x vector as [u, w, q, sita, v, p, fi, r, Omega]
	C.dXdvel[0] = A(0, 0) = dXdvel[0];
	C.dXdvel[2] = A(0, 1) = dXdvel[2] - C.omg_c[1];
	C.dXdomg[1] = A(0, 2) = dXdomg[1] - C.vel_c[2];
	A(0, 3) = -Gg * cos(-uctrl[5]);
	C.dXdvel[1] = A(0, 4) = dXdvel[1] + C.omg_c[2];
	C.dXdomg[0] = A(0, 5) = dXdomg[0];
	A(0, 6) = 0;
	C.dXdomg[2] = A(0, 7) = dXdomg[2] + C.vel_c[1];

	//for (int i = 0; i < 8; i++)
	//	cout << A(0, i) << '\t';
	//cout << endl;

	C.dZdvel[0] = A(1, 0) = dZdvel[0] + C.omg_c[1];
	C.dZdvel[2] = A(1, 1) = dZdvel[2];
	C.dZdomg[1] = A(1, 2) = dZdomg[1] + C.vel_c[0];
	A(1, 3) = -Gg*cos(uctrl[4])*sin(-uctrl[5]);
	C.dZdvel[1] = A(1, 4) = dZdvel[1] - C.omg_c[0];
	C.dZdomg[0] = A(1, 5) = dZdomg[0] - C.vel_c[1];
	A(1, 6) = -Gg*sin(uctrl[4])*cos(-uctrl[5]);
	C.dZdomg[2] = A(1, 7) = dZdomg[2];

	//for (int i = 0; i < 8; i++)
	//	cout << A(1, i) << '\t';
	//cout << endl;

	C.dMdvel[0] = A(2, 0) = Muu;
	C.dMdvel[2] = A(2, 1) = Mww;
	C.dMdomg[1] = A(2, 2) = Mqq;
	A(2, 3) = 0;
	C.dMdvel[1] = A(2, 4) = Mvv;
	C.dMdomg[0] = A(2, 5) = Mpp - 2 * C.omg_c[0] * C.inmatx_M(0, 2) / C.inmatx_M(1, 1) - C.omg_c[2] * (C.inmatx_M(0, 0) - C.inmatx_M(2, 2)) / C.inmatx_M(1, 1);
	A(2, 6) = 0;
	C.dMdomg[2] = A(2, 7) = Mrr - 2 * C.omg_c[2] * C.inmatx_M(0, 2) / C.inmatx_M(1, 1) - C.omg_c[0] * (C.inmatx_M(0, 0) - C.inmatx_M(2, 2)) / C.inmatx_M(1, 1);

	//for (int i = 0; i < 8; i++)
	//	cout << A(2, i) << '\t';
	//cout << endl;

	A(3, 2) = cos(-uctrl[5]);
	A(3, 7) = -sin(-uctrl[5]);

	//for (int i = 0; i < 8; i++)
	//	cout << A(3, i) << ' ';
	//cout << endl;

	C.dYdvel[0] = A(4, 0) = dYdvel[0] - C.omg_c[2];
	C.dYdvel[2] = A(4, 1) = dYdvel[2] - C.omg_c[0];
	C.dYdomg[1] = A(4, 2) = dYdomg[1];
	A(4, 3) = -Gg*sin(uctrl[4])*sin(-uctrl[5]);
	C.dYdvel[1] = A(4, 4) = dYdvel[1];
	C.dYdomg[0] = A(4, 5) = dYdomg[0] + C.vel_c[2];
	A(4, 6) = Gg*cos(uctrl[4])*cos(-uctrl[5]);
	C.dYdomg[2] = A(4, 7) = dYdomg[2] -C.vel_c[0];

	//for (int i = 0; i < 8; i++)
	//	cout << A(4, i) << '\t';
	//cout << endl;

	C.dLdvel[0] = A(5, 0) = Luu;
	C.dLdvel[2] = A(5, 1) = Lww;
	C.dLdomg[1] = A(5, 2) = Lqq + f1*C.omg_c[0] - f2*C.omg_c[2];
	A(5, 3) = 0;
	C.dLdvel[1] = A(5, 4) = Lvv;
	C.dLdomg[0] = A(5, 5) = Lpp +f1*C.omg_c[1];
	A(5, 6) = 0;
	C.dLdomg[2] = A(5, 7) = Lrr - f2*C.omg_c[1];

	//for (int i = 0; i < 8; i++)
	//	cout << A(5, i) << '\t';
	//cout << endl;

	A(6, 0) = A(6, 1) = 0;
	A(6, 2) = sin(uctrl[4])*tan(-uctrl[5]);
	A(6, 3) = 0;
	A(6, 4) = 0;
	A(6, 5) = 1;
	A(6, 6) = 0;
	A(6, 7) = cos(uctrl[4])*tan(-uctrl[5]);

	//for (int i = 0; i < 8; i++)
	//	cout << A(6, i) << '\t';
	//cout << endl;

	C.dNdvel[0] = A(7, 0) = Nuu;
	C.dNdvel[2] = A(7, 1) = Nww;
	C.dNdomg[1] = A(7, 2) = Nqq - f1*C.omg_c[2] - f3*C.omg_c[0];
	A(7, 3) = 0;
	C.dNdvel[1] = A(7, 4) = Nvv;
	C.dNdomg[0] = A(7, 5) = Npp - f3*C.omg_c[1];
	A(7, 6) = 0;
	C.dNdomg[2] = A(7, 7) = Nrr - f1*C.omg_c[1];

	//for (int i = 0; i < 8; i++)
	//	cout << A(7, i) << '\t';
	//cout << endl;

	//C.dXdvel = dXdvel, C.dYdvel = dYdvel, C.dZdvel = dZdvel;
	//C.dXdomg = dXdomg, C.dYdomg = dYdomg, C.dZdomg = dZdomg;
	C.dXdeul = dXdeul, C.dYdeul = dYdeul, C.dZdeul = dZdeul;

	//C.dLdvel = dLdvel, C.dMdvel = dMdvel, C.dNdvel = dNdvel;
	//C.dLdomg = dLdomg, C.dMdomg = dMdomg, C.dNdomg = dNdomg;
	//C.dLdeul = dLdeul, C.dMdeul = dMdeul, C.dNdeul = dNdeul;
	C.dLdeul[0] = Lfi, C.dLdeul[1] = Lsita;
	C.dMdeul[0] = Mfi, C.dMdeul[1] = Msita;
	C.dNdeul[0] = Nfi, C.dNdeul[1] = Nsita;

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			C.ALongM(i, j) = A(i, j);
			C.ALongMRe(i, j) = A(i, j);
			C.ALateM(i, j) = A(i + 4, j + 4);
			C.ALateMRe(i, j) = A(i + 4, j + 4);
		}
	}

	for (int i = 0; i < 8; i++)
		for (int j = 0; j < 8; j++)
			C.AMatrixRe(i, j) = A(i, j);

	C.ALongMRe(0, 3) += dXdeul[1];
	C.ALongMRe(1, 3) += dZdeul[1];
	C.ALongMRe(2, 3) += Msita;

	C.ALateMRe(0, 2) += dYdeul[0];
	C.ALateMRe(1, 2) += Lfi;
	C.ALateMRe(3, 2) += Nfi;

	C.AMatrixRe(0, 3) += dXdeul[1];
	C.AMatrixRe(1, 3) += dZdeul[1];
	C.AMatrixRe(2, 3) += Msita;
	C.AMatrixRe(4, 3) += dYdeul[1];
	C.AMatrixRe(5, 3) += Lsita;
	C.AMatrixRe(7, 3) += Nsita;

	C.AMatrixRe(0, 6) += dXdeul[0];
	C.AMatrixRe(1, 6) += dZdeul[0];
	C.AMatrixRe(2, 6) += Mfi;
	C.AMatrixRe(4, 6) += dYdeul[0];
	C.AMatrixRe(5, 6) += Lfi;
	C.AMatrixRe(7, 6) += Nfi;
}

void CopterSolver::_ComputeBMatrix(Matrix2<double> &B, Copter &C)
{
	double uctrl[6] = { 0.0 };
	double L00, Lcc, Lss, Ltt, N00, Ncc, Nss, Ntt, M00, Mcc, Mss, Mtt;
	double vector1, vector2, vector3, vector4, vector5;
	double unit_factor = 1.0;

	if (!C.fuselage.si_unit)
		unit_factor = 0.3048;

	B.allocate(xdof, udof); 
	C.GetCtrl(uctrl);

	vector1 = 1 / unit_factor / unit_factor * C.inmatx_M(2, 2) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector2 = 1 / unit_factor / unit_factor * C.inmatx_M(0, 2) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector3 = 1 / unit_factor / unit_factor * C.inmatx_M(0, 2) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector4 = 1 / unit_factor / unit_factor * C.inmatx_M(0, 0) / (C.inmatx_M(0, 0)*C.inmatx_M(2, 2) - C.inmatx_M(0, 2)*C.inmatx_M(0, 2));
	vector5 = (unit_factor * unit_factor * C.inmatx_M(1, 1));

	L00 = dot(vector1, vector2, 0, dLdctrl[0], dNdctrl[0], 0);
	Lcc = dot(vector1, vector2, 0, dLdctrl[1], dNdctrl[1], 0);
	Lss = dot(vector1, vector2, 0, dLdctrl[2], dNdctrl[2], 0);
	Ltt = dot(vector1, vector2, 0, dLdctrl[3], dNdctrl[3], 0);
	N00 = dot(vector3, vector4, 0, dLdctrl[0], dNdctrl[0], 0);
	Ncc = dot(vector3, vector4, 0, dLdctrl[1], dNdctrl[1], 0);
	Nss = dot(vector3, vector4, 0, dLdctrl[2], dNdctrl[2], 0);
	Ntt = dot(vector3, vector4, 0, dLdctrl[3], dNdctrl[3], 0);
	M00 = dMdctrl[0] / vector5;
	Mcc = dMdctrl[1] / vector5;
	Mss = dMdctrl[2] / vector5;
	Mtt = dMdctrl[3] / vector5;

	B(0, 0) = dXdctrl[0], B(0, 1) = dXdctrl[1], B(0, 2) = dXdctrl[2], B(0, 3) = dXdctrl[3];
	B(1, 0) = dZdctrl[0], B(1, 1) = dZdctrl[1], B(1, 2) = dZdctrl[2], B(1, 3) = dZdctrl[3];
	B(2, 0) = M00, B(2, 1) = Mcc, B(2, 2) = Mss, B(2, 3) = Mtt;
	B(4, 0) = dYdctrl[0], B(4, 1) = dYdctrl[1], B(4, 2) = dYdctrl[2], B(4, 3) = dYdctrl[3];
	B(5, 0) = L00, B(5, 1) = Lcc, B(5, 2) = Lss, B(5, 3) = Ltt;
	B(7, 0) = N00, B(7, 1) = Ncc, B(7, 2) = Nss, B(7, 3) = Ntt;

	C.dXdctrl = dXdctrl, C.dYdctrl = dYdctrl, C.dZdctrl = dZdctrl;
	C.dLdctrl[0] = L00, C.dLdctrl[1] = Lcc, C.dLdctrl[2] = Lss, C.dLdctrl[3] = Ltt;
	C.dMdctrl[0] = M00, C.dMdctrl[1] = Mcc, C.dMdctrl[2] = Mss, C.dMdctrl[3] = Mtt;
	C.dNdctrl[0] = N00, C.dNdctrl[1] = Ncc, C.dNdctrl[2] = Nss, C.dNdctrl[3] = Ntt;
}

