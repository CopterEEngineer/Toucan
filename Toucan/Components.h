#pragma once
#define Components_h

// header files
#include "MatrixTemplate.h"
#include "Coordinate.h"
#include "Algorithm.h"
#include <vector>


// macro definition
#define SLUG_CONST 32.174
#define UNIT_CONST (32.1522/32.174)
#define RAD(x) (PI*x)/180.0
#define DEG(x) (x*180.0)/PI

#define OUTPUT_MODE
//#define OUTPUT_MODE_1

// enum
enum CompType {
	Fuse = 0, Mrotor = 1, Trotor = 2, Hwing = 3, Vwing = 4
};

enum BladeSolType {
	Rotation = 0, HubFixed = 1
};

enum SimType {
	FreeTrim0 = 0, FreeTrim1 = 1, FreeTrim2 = 2, GeneralTrim = 3,
	WindTrim = 4, Transient = 5
};

enum AeroDynType {
	Averaged = 0, PWake = 1, FWake = 2
};

class Ambience
{
public:
	Ambience() { ; }
	~Ambience() { ; }
	void SetAmbience(void);

public:
	myTYPE rho, vsound;

};

class Fuselage
{
public:
	Fuselage() { ; }
	~Fuselage() { ; }
	void SetAirfm(void);
	void GetAirfm(myTYPE f[3], myTYPE m[3]);
	void GetAirfm_cg(myTYPE f[3], myTYPE m[3]);
	void SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc);
	void SetAirfm_cg(const Coordinate *base);
	void InitVariables(void);

private:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3];
public:
	Coordinate refcoord;
	myTYPE dragA;
	Ambience amb;
};

class Wing
{
public:
	Wing() { ; }
	~Wing() { ; }
	void SetAirfm(void);
	void GetAirfm(myTYPE f[3], myTYPE m[3]);
	void SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc);
	void SetAirfm_cg(const Coordinate *base);
	void GetAirfm_cg(myTYPE f[3], myTYPE m[3]);
	void InitVariables(void);

private:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3];
public:
	Coordinate refcoord;
	CompType type;
	myTYPE a0, cd0, cd1, cd2;
	myTYPE span, chord, taper;
	Ambience amb;
};

class BladeSolver
{
public:
	BladeSolver() { ; }
	~BladeSolver() { ; }

	template<class _Ty> void _GenArfPrepare(_Ty dt);
	template<class _Ty> bool _GenArfStarter(_Ty Qt, _Ty q, _Ty dq);
	template<class _Ty> void _GenArfTimeMarch(_Ty Qt, _Ty dt, int niter);

	bool isGenArfExit(int niter);
	void InitBladeSolver(bool teeter);

public:
	myTYPE Marf, Carf, Karf, Ka, err_b;
	Matrix2<myTYPE> Marf22, Carf22, Karf22, Ka22;
	myTYPE Qt, q0, dq0, ddq0, q1, q, dq, ddq;
	Matrix1<myTYPE> Qt22, q02, dq02, ddq02, q12, dq2, ddq2;
	Matrix1<myTYPE> sol;

	int nitermax, niter, dff, nperiod;
	myTYPE beta[3];
	GenArf GAf;
	BladeSolType soltype;
};

//class RotorWake
//{
//public:
//	RotorWake() { ; }
//	~RotorWake() { ; }
//	void InitWakePram(const int naz, const int nrs);
//	void Allocate(void);
//	//void BladePosition(myTYPE b[3]);
//	void WakeInducedVel(myTYPE *b);
//
//private:
//	void _BladePosition(myTYPE *b);
//	void _Bd_WakeGeo(void);
//	void _Bd_WakeGeo_v0(void);
//	void _Bd_WakeGeo_v1(void);
//	void _Bd_WakeGeo_v2(void);
//
//private:
//	Matrix2<myTYPE> tipstr, rotstr, shdstr, trlstr, cirlb;
//	Matrix3<myTYPE> bladedeform, tipgeometry;
//	Matrix2<myTYPE> lambdi, lambdh, lambdt, lambdx, lambdy;
//	myTYPE beta[3];
//	bool haveGeo;
//public:
//	int kwtip, kwrot, nk, nf, ns, ni, nbn, naf, nnr;
//	myTYPE rtip, rc0, outboard;
//};


class Rotor
{
public:
	Rotor() { ; }
	~Rotor() { ; }
	void SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc);
	void SetAirfm_cg(const Coordinate *base);
	void Allocate(void);
	void InitVariables(void);
	void SetCoordBase(const Coordinate *base);
	void SetAirfm(void);
	void GetAirfm(myTYPE f[3], myTYPE m[3]);
	void GetAirfm_cg(myTYPE f[3], myTYPE m[3]);
	void AvrgInducedVel(void);
	void BladeDynamics(void);
	void InitWakePram(void);
private:
	void _allocate(void);
	void _initvariables(void);
	void _avrgInducedVel_MR(void);
	void _avrgInducedVel_TR(void);
	void _avrgInducedVel_Rt(void);	
	void _avrgInducedVel_Fx(void);
	void _avrgInducedVel_TRFx(void);
	bool _avrgIndVelStarter(void);

	void _wakeStarter(void);
	void _wakeInducedVel(void);
	void _bonVortexStr(void);
	bool _tipVortexStr(void);
	void _bladePosition(void);
	bool _wakeGeoBd(void);
	void _wakeIndVelCalc(void);

	void _teeterdynamics_rt(void);
	void _hingedynamics_rt(void);
	void _teeterdynamics_fx(void);
	void _hingedynamics_fx(void);

	void _teeterflap_rt(void);
	void _hingeflap_rt(void);
	void _teeterflap_fx(void);
	void _hingeflap_fx(void);

	void _setairfm(Matrix1<double> &_dfx, Matrix1<double> &_dfz, Matrix1<double> &_dfr, const double &it, const double &b, const double &db, const int ia);
	void _setairfm(Matrix1<float> &_dfx, Matrix1<float> &_dfz, Matrix1<float> &_dfr, const float &it, const float &b, const float &db);
	void _setbladevelc(Matrix1<double> &_ut, Matrix1<double> &_up, double &ia, const double &it, const double &b, const double &db);
	void _setbladevelc(Matrix1<float> &_ut, Matrix1<float> &_up, float &ia, const float &it, const float &b, const float &db);
	void _setairfm_sp(double f[3], double m[3]);
	void _setairfm_sp(float f[3], double m[3]);


	template <class _Ty>
	void _aerodynacoef(Matrix1<_Ty> &_cl, Matrix1<_Ty> &_cd, Matrix1<_Ty> &incidn, Matrix1<_Ty> &ma_n);

	template <class _Ty>
	void _aerodynacoef(Matrix2<_Ty> &_cl, Matrix2<_Ty> &_cd, Matrix2<_Ty> &incidn, Matrix2<_Ty> &ma_n);

private:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3];
	int niter_w, niter_a;
	myTYPE mul;
	myTYPE beta[3];
	myTYPE lambdi_ag, lambdh_ag, lambdt_ag; // NOTE: these three variations defined at different coordinates
	myTYPE power, torque;
	Matrix2<myTYPE> bflap, dbflap, sfth;
	Matrix2<myTYPE> ut, un, up, ua, ma_n;
	Matrix1<myTYPE> _cl, _ua;
	Matrix2<myTYPE> incidn, cl, cd;
	Matrix2<myTYPE> lambdi, lambdh, lambdt, lambdx, lambdy;

	Matrix2<myTYPE> tipstr, rotstr, shdstr, trlstr, cirlb;
	Matrix3<myTYPE> bladedeform, tipgeometry;
	//Matrix2<myTYPE> lambdi, lambdh, lambdt, lambdx, lambdy;
	//myTYPE beta[3];
	bool haveGeo, haveStr;
	
public:
	CompType type;
	bool teeter;
	Coordinate refcoord;
	int nb, nf, ns, ni;
	myTYPE eflap, khub, del, pitchroot, radius, bt, rroot, disk_A;
	myTYPE precone, omega, vtipa;
	myTYPE sigma, gama, a0;
	Matrix1<myTYPE> chord, twist, sweep;
	myTYPE iflap, m1;
	Coordinate hubfxcoord, hubrtcoord, bladecoord, tppcoord;
	Matrix2<myTYPE> cltc, cdtc, cmtc; //
	Matrix2<myTYPE> rastation, ristation, azstation; //, chord , twist, sweep;
	myTYPE t0;
	//GenArf Gaf;
	BladeSolver bld;
	Ambience amb;

	myTYPE sita[3];
	AeroDynType adyna;
	//RotorWake wake;

	int kwtip, kwrot, nk, nbn, naf, nnr;
	//int nf, ns, ni;
	myTYPE rtip, rc0, outboard;

};

class ModelCase
{
public:
	ModelCase() { ; }
	~ModelCase() { ; }
	virtual void GetModel();
	virtual void GetProb();

public:
	std::vector<Wing> WingV;
	std::vector<Rotor> RotorV;
	Fuselage fuselage;
	myTYPE mass, inmatx[3][3];
	Matrix2<myTYPE> inmatx_M;
	Ambience amb;
	Coordinate refcoord;
	SimType simtype;
	int nfree;
};

class Copter
{
public:
	Copter() { ; }
	~Copter() { ; }
	void InitProject(void);
	void InitRotorCraft(ModelCase &M);
	void InitCoord(ModelCase &M);
	void SetSimCond(const int ic);
	void CopterSimulation(void);
	void SetStates(void);
	void PostProcess(void);
	void functest(void);

	template <class _Ty> 
	void GetCtrl(_Ty *xctrl);

	template <class _Ty>
	void SetCtrl(_Ty *xctrl, Rotor &R);

	template <class _Ty>
	void GetStates(_Ty v[3], _Ty w[3], _Ty dv[3], _Ty dw[3]);

	template <class _Ty>
	void SetStates(_Ty v[3], _Ty w[3], _Ty dv[3], _Ty dw[3]);

private:
	myTYPE airforce_sigma[3], airmoment_sigma[3]; // synthesized airdynamics												  

	Ambience amb;
	ModelCase model;
	//CopterSolver CSolver;

public:
	SimType simtype;
	std::vector<Wing> WingV;
	std::vector<Rotor> RotorV;
	Fuselage fuselage;
	Coordinate refcoord;

	myTYPE mass, inmatx[3][3];
	Matrix2<myTYPE> inmatx_M;

	myTYPE vel_c[3], omg_c[3], dvel_c[3], domg_c[3];
	myTYPE vel_g[3], omg_g[3], dvel_g[3], domg_g[3];

	int nfree;
	int Ncase;
	myTYPE controls[6];

private:
	void _Allocate(void);
	void _InitVariables(void);
};


template <class _Ty>
void _setstates(_Ty _vel[3], _Ty _omg[3], _Ty _dvel[3], _Ty _domg[3], const _Ty *vc, const _Ty *wc, const _Ty *dvc, const _Ty *dwc, const Coordinate &refcoord)
{
	// transfer needed component base coordinate to copter.

	//Type temp_a1[3] = { 0,0,0 };
	_Ty temp_v[3], temp_w[3], temp_dv[3], temp_dw[3];
	for (int i = 2; i >= 0; --i) {
		_vel[i] = _omg[i] = _dvel[i] = _domg[i] = 0;
		temp_v[i] = temp_w[i] = temp_dv[i] = temp_dw[i] = 0;
	}
	Cross(_vel, wc, refcoord.origin);
	Cross(temp_dv, wc, _vel); // w x (w x r)
	Cross(_dvel, dwc, refcoord.origin); // e x r

	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			_vel[i] += refcoord.Ttransf[i][j] * (_vel[j] + vc[j]);
			_omg[i] += refcoord.Etransf[i][j] * wc[j];
			_dvel[i] += refcoord.Ttransf[i][j] * (temp_dv[j] + dvc[j] + _dvel[j]);
			_domg[i] += refcoord.Etransf[i][j] * dwc[j];
		}
	}
}

template <class _Ty>
bool _setairfm_cg(_Ty fcg[3], _Ty mcg[3], const _Ty *f, const _Ty *m, const Coordinate &refcoord, const Coordinate *base)
{
	_Ty rcb[3];
	const _Ty *origin1_ptr = refcoord.origin;
	const _Ty *origin2_ptr = base->origin;

	for (int j = 2; j >= 0; --j) { fcg[j] = mcg[j] = 0; }

	//refcoord.Transfer(tcb, rcb, refcoord, *base);
	if (refcoord.base == base) {
		for (int i = 2; i >= 0; --i) {
			rcb[i] = *(origin2_ptr + i) - *(origin1_ptr + i);
		}
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				// give a transpose by j i swap of Ttransf[][]
				//airforce_cg[i] += refcoord.Ttransf[j][i] * airforce[j];
				// index friendly
				fcg[j] += refcoord.Ttransf[i][j] * f[i];
			}
		}
		Cross(mcg, fcg, rcb);
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				// give a transpose by j i swap of Ttransf[][]
				//airmoment_cg[i] += refcoord.Ttransf[j][i] * airmoment[j];
				// index friendly
				mcg[j] += refcoord.Ttransf[i][j] * m[i];
			}
		}
	}
	else { 
		printf("Wrong Coordinate in _setairfm_cg() Func. \n");
		return false;
	}
	return true;
}

template<class _Ty> 
void BladeSolver::_GenArfPrepare(_Ty dt)
{
	GAf.InitCoef(dt);
	Karf = Ka*GAf.ck + Marf*GAf.c0 + Carf*GAf.c1;
}

template<class _Ty>
bool BladeSolver::_GenArfStarter(_Ty _Qt, _Ty _q, _Ty _dq)
{
	ddq = Marf / (_Qt - Carf*_dq - Karf*_q) ;
	dq = _dq;
	q = _q;
	sol(0) = q;
	return true;
}

template<class _Ty>
void BladeSolver::_GenArfTimeMarch(_Ty _Qt, _Ty dt, int niter)
{
	_Qt += Marf*(q*GAf.c0 + dq*GAf.c2 + ddq*GAf.c3) - Ka*q*GAf.af;
	_Qt += Carf*(q*GAf.c1 + dq*GAf.c4 + ddq*GAf.c5);
	q1 = _Qt / Karf;
	dq += ddq*(1 - GAf.r)*dt;
	ddq = (q1 - q)*(GAf.bt / dt / dt) - dq*GAf.bt / dt + ddq*(GAf.bt / 2.0 - GAf.r*GAf.bt + 1.0);
	dq += ddq*GAf.r*dt;
	q = q1;
	sol(niter) = q;
}

template <class _Ty>
void Rotor::_aerodynacoef(Matrix1<_Ty> &_cl, Matrix1<_Ty> &_cd, Matrix1<_Ty> &incidn, Matrix1<_Ty> &ma_n)
{
	_cl = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn / PI*180.0, ma_n);
	_cd = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn / PI*180.0, ma_n);
}

template <class _Ty>
void Rotor::_aerodynacoef(Matrix2<_Ty> &_cl, Matrix2<_Ty> &_cd, Matrix2<_Ty> &incidn, Matrix2<_Ty> &ma_n)
{
	_cl = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn / PI*180.0, ma_n);
	_cd = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn / PI*180.0, ma_n);
}

template <class _Ty> 
void Copter::GetCtrl(_Ty *xctrl)
{
	for (auto it = RotorV.cbegin(); it != RotorV.cend(); ++it)
	{
		switch (it->type)
		{
		case Mrotor:
			xctrl[0] = it->sita[0];
			xctrl[1] = it->sita[1];
			xctrl[2] = it->sita[2];
			break;
		case Trotor:
			xctrl[3] = it->sita[0];
			break;
		default:
			break;
		}
	}
	switch (simtype)
	{
	case FreeTrim0:
		printf("Undefined. \n");
		return;
	case FreeTrim1:
		xctrl[4] = refcoord.euler[0];
		xctrl[5] = refcoord.euler[1];
		break;
	case FreeTrim2:
		xctrl[4] = refcoord.euler[2];
		xctrl[5] = refcoord.euler[1];
		break;
	default:
		break;
	}
}

template <class _Ty>
void Copter::SetCtrl(_Ty *xctrl, Rotor &R)
{
	for (int i = 0; i < 3; ++i)
		R.sita[i] = xctrl[i];
}

template <class _Ty>
void Copter::GetStates(_Ty v[3], _Ty w[3], _Ty dv[3], _Ty dw[3])
{
	for (int i = 0; i < 3; ++i)
	{
		v[i] = vel_c[i];
		w[i] = omg_c[i];
		dv[i] = dvel_c[i];
		dw[i] = domg_c[i];
	}
}

template <class _Ty>
void Copter::SetStates(_Ty v[3], _Ty w[3], _Ty dv[3], _Ty dw[3])
{
	for (int i = 2; i >= 0; --i)
		v[i] = w[i] = dv[i] = dw[i] = 0;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			v[i] += refcoord.Ttransf[i][j] * vel_g[j];
			w[i] += refcoord.Etransf[i][j] * omg_g[j];
		}
	}
	switch (simtype)
	{
	case FreeTrim0:
	case FreeTrim1:
	case FreeTrim2:
	case WindTrim:
		dv[0] = dv[1] = dv[2] = 0;
		dw[0] = dw[1] = dw[2] = 0;
		break;
	case GeneralTrim:
		break;
	case Transient:
		break;
	default:
		break;
	}
}
