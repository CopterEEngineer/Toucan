#pragma once
#define Solver_h

#include "Components.h"

class CopterSolver
{
public:
	CopterSolver() { ; }
	~CopterSolver() { ; }

	void InitCopterSolver(void);
	void CopterSimulation(Copter &C);

private:

	void _FreeTrimSolver(Copter &C);
	void _WindTrimSolver(Copter &C);
	void _TransientSolver(Copter &C);
	void _CompsSetAirFM(Copter &C);
	void _Assemble(Copter &C);
	void _EnableWake(Copter &C);

	template <class _Ty>
	void _CompsSetStates(Copter &C, const _Ty *vc, const _Ty *wc, const _Ty *dvc, const _Ty *dwc);

	template <class _Ty>
	void _SetDerivs(_Ty dv[3], _Ty dw[3], const _Ty *f, const _Ty *m, Copter &C);

	template <class _Ty>
	void _UpdateDerivs(_Ty dv[3], _Ty dw[3], Rotor &R, Copter &C);

	template <class _Ty>
	void _UpdateEuler(_Ty *xctrl, Copter &C);

	template <class _Ty>
	void _SetJacob(_Ty *dvp, _Ty *dvn, _Ty *dwp, _Ty *dwn, const int dimid);

	template <class _Ty>
	int _UpdateCtrls(_Ty xctrl[6], _Ty *deltt);

	template <class _Ty>
	void _UpdateRotorCtrls(_Ty xctrl[6], Copter &C);

	template <class _Ty>
	bool isExit(_Ty *xctrl, _Ty *deltt, const int iter);

public:
	SimType simtype;
	double err_a, err_c, epsilon;
	Matrix1<double> sum_a1_del, sum_a2_del, max_c_del;
	int nitermax;
	double SigmaF[3], SigmaM[3];
	double dv[3], dw[3];
	double jacob66[6][6], jacob33[3][3];
	Matrix2<double> jacobM;
	double sita_coll_max, sita_cycl_max, euler_max;
	double xctrl_temp[6];
};



template <class _Ty>
void CopterSolver::_CompsSetStates(Copter &C, const _Ty *vc, const _Ty *wc, const _Ty *dvc, const _Ty *dwc)
{
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
		C.RotorV[i].SetStates(vc, wc, dvc, dwc);
	for (int i = C.WingV.size() - 1; i >= 0; --i)
		C.WingV[i].SetStates(vc, wc, dvc, dwc);
	C.fuselage.SetStates(vc, wc, dvc, dwc);
}

template <class _Ty>
void CopterSolver::_SetDerivs(_Ty dv[3], _Ty dw[3], const _Ty *f, const _Ty *m, Copter &C)
{
	_Ty temp_cross_v[3] = { 0 };
	_Ty temp_cross_o[3] = { 0 };
	_Ty temp_omg[3] = { 0 };

	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			temp_omg[i] += C.inmatx[i][j] * C.omg_c[j];
		}
	}

	Cross(temp_cross_o, C.omg_c, temp_omg);
	Cross(temp_cross_v, C.omg_c, C.vel_c);

	for (int i = 2; i >= 0; --i) {
		dv[i] = f[i] / (C.mass / SLUG_CONST) + C.refcoord.Ttransf[i][2] * UNIT_CONST * SLUG_CONST - temp_cross_v[i];
		dw[i] = m[i] - temp_cross_o[i];
	}

	Msolver(C.inmatx_M.v_p, dw, 3, 1);
}

template <class _Ty>
void CopterSolver::_UpdateDerivs(_Ty dv[3], _Ty dw[3], Rotor &R, Copter &C)
{
	R.SetAirfm();
	_Assemble(C);
	_SetDerivs(dv, dw, SigmaF, SigmaM, C);
}

template <class _Ty>
void CopterSolver::_UpdateEuler(_Ty *xctrl, Copter &C)
{
	_Ty euler[3];
	switch (simtype)
	{
	case FreeTrim0:
		printf("Undefined. \n");
		return;
	case FreeTrim1:
		euler[0] = xctrl[4];
		euler[1] = xctrl[5];
		euler[2] = 0;
		break;
	case FreeTrim2:
		euler[0] = 0;
		euler[1] = xctrl[5];
		euler[2] = xctrl[4];
		break;
	default:
		break;
	}
	C.refcoord.SetCoordinate(euler, "euler");
}

template <class _Ty>
void CopterSolver::_SetJacob(_Ty *dvp, _Ty *dvn, _Ty *dwp, _Ty *dwn, const int dimid)
{
	switch (simtype)
	{
	case FreeTrim0:
	case FreeTrim1:
	case FreeTrim2:
		for (int i = 0; i < 3; ++i)
		{
			jacobM(i, dimid) = jacob66[i][dimid] = (dvp[i] - dvn[i]) / 2.0 / epsilon;
			jacobM(i + 3, dimid) = jacob66[i + 3][dimid] = (dwp[i] - dwn[i]) / 2.0 / epsilon;		
		}
		break;
	case GeneralTrim:
		break;
	case WindTrim:
		break;
	case Transient:
		break;
	default:
		break;
	}
}

template <class _Ty>
int CopterSolver::_UpdateCtrls(_Ty xctrl[6], _Ty *deltt)
{
	int ic_flg = 0;
	int nfree = 0;
	for (int i = 0; i < 6; ++i)
		xctrl[i] -= deltt[i];
	if (Abs(xctrl[0]) > sita_coll_max)
	{
		xctrl[0] = 0.5 * sita_coll_max*Sign(xctrl[0]) + 0.5 * (xctrl[0] + deltt[0]);
		++ic_flg;
	}
	switch (simtype)
	{
	case FreeTrim0:
	case FreeTrim1:
	case FreeTrim2:
	case GeneralTrim:
		nfree = 6;
		if (Abs(xctrl[3]) > sita_coll_max)
		{
			xctrl[3] = 0.5 * sita_coll_max*Sign(xctrl[3]) + 0.5 * (xctrl[3] + deltt[3]);
			++ic_flg;
		}
		if (Abs(xctrl[1]) > sita_cycl_max)
		{
			xctrl[1] = 0.5 * sita_cycl_max*Sign(xctrl[1]) + 0.5*(xctrl[1] + deltt[1]);
			++ic_flg;
		}
		if (Abs(xctrl[2]) > sita_cycl_max)
		{
			xctrl[2] = 0.5 * sita_cycl_max*Sign(xctrl[2]) + 0.5*(xctrl[2] + deltt[2]);
			++ic_flg;
		}
		if (Abs(xctrl[4]) > euler_max)
		{
			xctrl[4] = 0.5 * euler_max*Sign(xctrl[4]) + 0.5*(xctrl[4] + deltt[4]);
			++ic_flg;
		}
		if (Abs(xctrl[5]) > euler_max)
		{
			xctrl[5] = 0.5 * euler_max*Sign(xctrl[5]) + 0.5*(xctrl[5] + deltt[5]);
			++ic_flg;
		}
		break;
	case WindTrim:
		break;
	default:
		break;
	}

#ifdef OUTPUT_MODE
	cout << endl;
	printf("Delta controls (degs).\n");
	for (int iq = 0; iq < nfree; ++iq) {
		cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
		cout << DEG(deltt[iq]) << "\t";
	}
	cout << endl << endl;

	printf("New controls (degs).\n");
	for (int iq = 0; iq < nfree; ++iq) {
		cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
		cout << DEG(xctrl[iq]) << "\t";
	}
	cout << endl << endl;
#endif // !OUTPUT_MODE_1

	return ic_flg;
}

template <class _Ty>
void CopterSolver::_UpdateRotorCtrls(_Ty xctrl[6], Copter &C)
{
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		switch (C.RotorV[i].type)
		{
		case Trotor:
			C.SetCtrl((xctrl + 3), C.RotorV[i]);
			break;
		case Mrotor:
			C.SetCtrl(xctrl, C.RotorV[i]);
			break;
		default:
			break;
		}
	}
}

template <class _Ty>
bool CopterSolver::isExit(_Ty *xctrl, _Ty *deltt, const int iter)
{
	//static Matrix1<double> sum_a1_del(nitermax), sum_a2_del(nitermax), max_c_del(nitermax);
	//static double xctrl_temp[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
	//if (iter == 0)
	//{
	//	//sum_a1_del.allocate(nitermax), sum_a2_del.allocate(nitermax), max_c_del.allocate(nitermax);
	//	xctrl_temp[0] = xctrl_temp[1] = xctrl_temp[2] = xctrl_temp[3] = xctrl_temp[4] = xctrl_temp[5] = 1.0;
	//}
	if (iter == 0)
	{
		xctrl_temp[0] = xctrl_temp[1] = xctrl_temp[2] = xctrl_temp[3] = xctrl_temp[4] = xctrl_temp[5] = 1.0;
	}

	max_c_del(iter) = Abs(xctrl[0] / xctrl_temp[0] - 1.0);
	xctrl_temp[0] = xctrl[0];
	for (int i = 1; i < 6; ++i)
	{
		max_c_del(iter) = Max(max_c_del(iter), Abs(xctrl[i] / xctrl_temp[i] - 1.0));
		xctrl_temp[i] = xctrl[i];
	}
	sum_a1_del.v_p[iter] = deltt[0] * deltt[0] + deltt[1] * deltt[1] + deltt[2] * deltt[2];
	sum_a2_del.v_p[iter] = deltt[3] * deltt[3] + deltt[4] * deltt[4] + deltt[5] * deltt[5];

	if (Max(sum_a1_del.v_p[iter], sum_a2_del.v_p[iter]) >= err_a*err_a || max_c_del.v_p[iter] > err_c)
	{
		if (iter < nitermax - 1)
			return false;
		else
		{

			cout << endl;
			printf("**********************************************\n");
			printf("*        Warning: may not be convergent.     *\n");
			printf("**********************************************\n");

			printf("Iter: %d, Total aerodynamics at CG: \n", iter);

			printf("Finally controls: \n");
			for (int iq = 0; iq < 6; ++iq) {
				cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
				cout << DEG(xctrl[iq]) << "\t";
			}
			cout << endl;
			printf("Translation Accelation: %f \n", sum_a1_del(iter));
			printf("Angular Accelation: %f \n", sum_a2_del(iter));
			printf("Max Controls Delt: %f \n", max_c_del(iter));

			sum_a1_del.outputs("sum_a1_del.output", 8);
			sum_a2_del.outputs("sum_a2_del.output", 8);
			max_c_del.outputs("max_c_del.output", 8);

			return false;
		}
	}
	else
	{
		cout << endl;
		printf("Convergent. \n");
		printf("Iter: %d, Total aerodynamics at CG: \n", iter);
#ifdef OUTPUT_MODE_1
		printf("Accelaration error: %f, %f\n", sum_a1_del.v_p[iter], sum_a2_del.v_p[iter]);
		printf("Control converge: %f \n", max_c_del.v_p[iter]);
#endif // OUTPUT_MODE_1
		printf("Finally controls: \n");
		for (int iq = 0; iq < 6; ++iq) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << DEG(xctrl[iq]) << "\t";
		}

#ifdef _DEBUG
		sum_a1_del.outputs("sum_a1_del.output", 8);
		sum_a2_del.outputs("sum_a2_del.output", 8);
		max_c_del.outputs("max_c_del.output", 8);
#endif // _DEBUG
		return true;
	}
}