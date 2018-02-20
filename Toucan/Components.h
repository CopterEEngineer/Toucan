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
#define HORSEPOWER (0.0018182)
#define RAD(x) (PI*(x))/180.0
#define DEG(x) ((x)*180.0)/PI

#ifdef _DEBUG
#define OUTPUT_MODE
#define OUTPUT_MODE_1
#endif // _DEBUG

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
	Simple = -1, Averaged = 0, LinearInflow = 1, PWake = 101, FWake = 102,
};

enum JobsType
{
	SimTrim = 0, RPMSwp = 1, ChordSwp = 2, RadiusSwp = 3, Stable = 4,
};

enum HingeType {
	Teeter = 0, Cantilever = 1, Hinged = 2, Hingeless = 3
};

enum FuseModeling {
	Parasite = 0, Fitting = 1
};

enum WingModeling {
	WSimple = 0, WFitting = 1
};

enum AirfoilAero {
	Liftslope = 0, C81Table = 1, Padfield = 2, LBStallMethod = 3
};

class Monitor
{
public:
	Monitor() { ; }
	~Monitor() { ; }

public:
	double TEMP[100];
	double AOA, CL, KA, Alpha, Beta, CM, CX, CY, CZ, CML;
	double af[3], mf[3], vel[3], omg[3];
	double afcg[3], mfcg[3];
	double afp[3], afn[3], mfp[3], mfn[3];
	double afpB[9][3], afnB[9][3], mfpB[9][3], mfnB[9][3];
	double afcgd[2][3], mfcgd[2][3];
	double AfcgPB[9][3], AfcgNB[9][3], MfcgPB[9][3], MfcgNB[9][3];
	double Alphad[2], Betad[2];
	double KAd[2], AOAd[2];
	double Nvel[3], Npvel[3];
	double veld[2][3], omgd[2][3];
	double niterw[2], errw[2];
	double flapc[3];
	double flapcd[2][3];
	double velw[3], omgw[3];
	double velwd[2][3], omgwd[2][3];
	double errdb2[2], errdw2[2];
	double VelPB[9][3], OmgPB[9][3], VelNB[9][3], OmgNB[9][3];
	double VelwPB[9][3], OmgwPB[9][3], VelwNB[9][3], OmgwNB[9][3];
	double errw2, errb2;
	double LBerrSum;
	int Countsb, Countsw;
	double Ct_s;
};


class Ambience
{
public:
	Ambience() { ; }
	~Ambience() { ; }
	void SetAmbience(double);

public:
	myTYPE rho, vsound, height;

};

class Fuselage
{
public:
	Fuselage() { ; }
	~Fuselage() { ; }
	void SetAirfm(void);
	void SetAirfm(double);
	void SetAirfm(double, double, double);
	void GetAirfm(myTYPE f[3], myTYPE m[3]);
	void GetAirfm_cg(myTYPE f[3], myTYPE m[3]);
	void SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc);
	void SetAirfm_cg(const Coordinate *base);
	void InitVariables(void);
	void GetAngle(double &, double&);
	void GetStates(double v[3], double w[3]);

private:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3];
	double Alpha, Beta;

public:
	Coordinate refcoord;
	myTYPE dragA, Sp, Ss, Lf;
	Ambience amb;
	bool si_unit;
	FuseModeling fmdling;
	double Vtest;
	double xf0, xf1, xf2, xf3;
	double yf0, yf1, yf2, yf3;
	double zf0, zf1, zf2, zf3;
	double mf0, mf1, mf2, mf3;
	double nf0, nf1, nf2, nf3;
	Matrix2<double> cxtc, cytc, cztc, cmtc, cntc;
	double KLT, Inter0, Inter1;
	double Krho;
	Monitor monitor;
};

class Wing
{
public:
	Wing() { ; }
	~Wing() { ; }
	void SetAirfm(void);
	void SetAirfm(double);
	void GetAirfm(myTYPE f[3], myTYPE m[3]);
	void SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc);
	void SetAirfm_cg(const Coordinate *base);
	void GetAirfm_cg(myTYPE f[3], myTYPE m[3]);
	void InitVariables(void);
	void GetAngle(double &);
	void GetStates(double v[3], double w[3]);
private:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3];
	double AOA;
public:
	Coordinate refcoord;
	CompType type;
	myTYPE span, chord, taper;
	Ambience amb;
	bool si_unit;
	WingModeling wmdling;
	double a0, a1, a3, a5, alpha0, cd0, cd1, cd2;
	double KLT, KLTail, Inter0, Inter1;

	Monitor monitor;
};

class BladeSolver
{
public:
	BladeSolver() { ; }
	~BladeSolver() { ; }

	template<class _Ty> void _GenArfPrepare(_Ty dt);
	template<class _Ty> bool _GenArfStarter(_Ty Qt, _Ty q, _Ty dq);
	template<class _Ty> bool _GenArfStarter(_Ty Qt, _Ty q, _Ty dq, int n);
	template<class _Ty> void _GenArfTimeMarch(_Ty Qt, _Ty dt, int niter);

	bool isGenArfExit(int niter);
	void InitBladeSolver(bool teeter);

public:
	myTYPE Marf, Carf, Karf, Ka, epsb, err_b;
	//Matrix2<myTYPE> Marf22, Carf22, Karf22, Ka22;
	myTYPE Qt, q0, dq0, ddq0, q1, q, dq, ddq;
	//Matrix1<myTYPE> Qt22, q02, dq02, ddq02, q12, dq2, ddq2;
	Matrix1<myTYPE> sol, azmuth;

	int nitermax, niter, dff, nperiod;
	myTYPE beta[3];
	GenArf GAf;
	BladeSolType soltype;
};

class Airfoil
{
public:
	Airfoil() { ; }
	~Airfoil() { ; }

	void SetAirfoil(string, int, int);

public:
	string nid;
	Matrix1<double> Ma, CLa, alpha1, dalpha1, alpha0;
	Matrix1<double> s1, s2, K0, K1, K2, CD0, Df, CNl;
	Matrix1<double> Tp, Tf, Tv, Tvl;
	double r0, alphads0, alphass, alphacr, Talpha;
	double alphamin0, alphamin, Tr;
};

class LBStall
{
public:
	enum FlowState
	{
		Attach = 0, Separate = 1, ReAttach = 2
	};
	enum VortexState
	{
		NoVortex = 0, Conv = 1, Nearby = 2, Far = 3, Secd = 4
	};
	typedef std::pair<FlowState, VortexState> StatePair;


	LBStall() { ; }
	~LBStall() { ; }

	void Allocate(int, int, int);
	void SetConstants(void);
	void Prepare(Airfoil af, double c, double a, double dt, int nf, int ns, int nk, Matrix2<double> &aoa, Matrix2<double> &Ma, Matrix2<double> &q0, double r);
	void Prepare(Airfoil af, double c, double a, double dt, int nf, int ns, int nk, Matrix2<double> &aoa, Matrix2<double> &Ma, Matrix2<double> &q0);
	void Prepare(double c, double a, double dt, int nf, int ns, int nk, Matrix2<double> &aoa, Matrix2<double> &Ma, Matrix2<double> &q0);
	void Starter(double c, double a, double dt);
	void Complete(void);
	void Complete(Matrix1<double> &cl, Matrix1<double> &cd, Matrix1<double> &cn, Matrix1<double> &cc);
	void AttachFlow(int , int);
	void AttachFlow(int);
	void DynamicStall(void);
	bool isExit(int &);
	bool isExit(int &, int &);

	void Save(int);
	void Save(int, int);

	void FuncTest(void);
	bool Solver(int &, int, Matrix2<double> &aoa, Matrix2<double> &Ma, Matrix2<double> &aoad);
	void Display(int);

private:
	void _allocate(void);
	void _allocate1(void);
	void performLBStall(void);
	void performAttach(void);
	void performReAttach(void);
	void performSeparate(void);
	void _performKflow(double);
	void _performdstall(double);

	void performVortex(void);
	void _performnvortex(void);
	void _performconv(void);
	void _performnearby(void);
	void _performfar(void);
	void _performsecd(void);

	bool _startReAttach(void);
	bool _startAttach(void);
	bool _enableSecd(void);

	bool _ShengCriterion(void);

	StatePair stateChange(StatePair);
	StatePair stateChange(double);


public:
	Airfoil airfoil;
	Matrix1<double> alphaIn;
	Matrix1<double> x1, y1, D1, D2;
	Matrix1<double> CNC, CNI, CNP, CNaI, CNqI;
	Matrix1<double> CNf, CNfC;
	Matrix1<double> fpp, fp, Df, Dfp, alphaeff, alphaE;
	Matrix1<double> CNPrevised, CNfrevised, Dp;
	Matrix1<double> Cv, CNv;
	Matrix2<double> CNvSecd;
	Matrix1<double> CNT, CCf, CD, CL;
	Matrix1<double> alphak, qk, Mak;
	double dt, ds, alphaLc, MaLc, CNILc, DFLc, vsound, Cd0Lc;
	double A1, A2, b1, b2, beta2, CNMa;
	double Tl, Ta, Tq, Tf, Tp, Tv, Tvl, Ts, dalpha1;
	double s1, s2, alpha1, alpha0, yita, eps;
	Matrix2<double> aoa0M2, MaM2, q0M2, fppM2;
	Matrix1<double> aoa0M1, MaM1, q0M1;
	Matrix2<double> CNCM2, CNIM2, CNPM2, CNaIM2, CNqIM2;
	Matrix2<double> CNfM2, CNfCM2, CDM2, CLM2;
	Matrix2<double> CNPrevisedM2, CvM2, CNvM2, CNTM2, CCfM2;
	int countk, countk1, Nf, Ns, Nk;
	double fpptemp;
	Matrix1<double> MaIn;
	Matrix1<double> AoALSSave, AoATSSave, AoAATSave;
	Matrix2<double> AoALSSaveM2, AoATSSaveM2, AoAATSaveM2;
	Matrix1<int> Circle;

	Matrix1<int> FlowStateSave, VortexStateSave;
	Matrix2<int> FlowStateSaveM2, VortexStateSaveM2;

	StatePair  state;
	double _dak, _dqk;
	double r0, alphads0, alphass, alphacr, Talpha;
	double req;
	double alphamin0, alphamin, Tr;
	Matrix1<double> Da;
	bool enable, isShengCur, isShengPrv, secEnable;
	double err;
	double cl, cd;
private:
	double tv, _aeff, _Tf, _Tv;
	bool _secdcomfd;
	int secdID;
	bool isReverse;
};

class LSCorr
{
public:
	LSCorr() { ; }
	~LSCorr() { ; }

	void Initvariable(void);
	double LSCorrection(double az, double c, double rc, double ri0, double rj0, double rk0, double ri1, double rj1, double rk1);

public:
	Matrix1<double> angls, b0ls, b1ls, b2ls, aopls, a1ls, a2ls, copls, c0ls, c1ls, c2ls;
	Matrix1<double> angll, b1ll, aopll, a1ll, a2ll, c0ll;
	Matrix2<double> LScoeff, LLcoeff;
	bool enable;
};


class Wake
{
public:
	enum WakeType
	{
		Bound = 0, Trail = 1, Shed = 2
	};

	Wake() { ; }
	~Wake() { ; }

	void OutPutWake(const int ic);
	void OutPutWake(string s, const int ic);
	
	void InitVariable(const int nf, const int ns, const int nk);
	void ComputeIndVel(double &, double &, double&, const int, const int, const int, const int, double, double *);


	void _computeIndVel(int nb, Matrix1<double> &ch);
	bool _computeVorStr(Matrix2<double> &);
	bool _computeGeometry(const AeroDynType aero, double *veltpp, double lambdi_ag, double b0);
	void _airloadingPoint(Matrix2<double> &az, Matrix2<double> &ra, Matrix1<double> &ch, Matrix1<double> tw, double *b, double *s, double offset, double sr);

public:
	WakeType waketype;
	LSCorr lscorr;
	bool istip, isrot, isbond, isshed;
	bool havegeo, havestr;
	bool outputwake;
	int nf, ns, nk;
	double twistv, chv, rv, rc0, rc;
	double strboard;

	Monitor monitor;

	double lambdx, lambdy, lambdz;

	Matrix2<myTYPE> vortexstr;
	Matrix3<myTYPE> geometry;
	Matrix3<myTYPE> comppoint;

	Matrix2<myTYPE> geoexp;

	Matrix3<myTYPE> geometryathub;
	Matrix3<myTYPE> comppointathub;
};

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
	void GetBeta(myTYPE b[3]);
	void AvrgInducedVel(void);
	void BladeDynamics(void);
	void WakeModelPrams(const double k);
	void WakeInducedVel(void);
	void OutPutWake(const int ic);
	void OutPutWake(string s, const int ic);
	void GetPower(myTYPE p[6], myTYPE t[6]);
	void DiskOutput(string s);
	double GetLambdi(void);
	void GetStates(double v[3], double w[3]);
	void GetWStates(double v[3], double w[3]);
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
	void _wakeInducedVelMP(void);
	void _wakeInducedVelMP(int);
	void _bonVortexStr(void);
	bool _tipVortexStr(void);
	void _bladePosition(void);
	bool _wakeGeoBd(void);
	void _wakeIndVelCalc(void);

	void _teeterdynamics_rt(void);
	void _hingedynamics_rt(void);
	void _hingelessdynamics_rt(void);
	void _teeterdynamics_fx(void);
	void _hingedynamics_fx(void);
	void _hingelessdynamics_fx(void);

	void _teeterflap_rt(void);
	void _hingeflap_rt(void);
	void _hingelessflap_rt(void);
	void _teeterflap_fx(void);
	void _hingeflap_fx(void);
	void _hingelessflap_fx(void);
	void _flapMBC(double, double);
	void _flapWang(void);
	bool _hingelessflap_rt(double f[3], double m[3]);

	double _aerodynamics(double, double*);
	
	void _windcoordVel(double v[3], double dv[3]);
	void _windcoordOmg(double w[3], double dw[3]);
	void _velTransform(double v[3], double dv[3], Coordinate &);

	void _setairfm(Matrix1<double> &_dfx, Matrix1<double> &_dfz, Matrix1<double> &_dfr, const double &it, const double &b, const double &db, const int ia);
	double _setairfm(double &_dfx, double &_dfz, double &_dfr, const double &cl, const double &cd, const int, const int);
	double _setairfm(double &_dfx, double &_dfz, const double &cl, const double &cd, const int iz, const int ir);
	void _setbladeaeros(double &_ma, double &_aoa, double &_infl, const double &az, const double &ra, const double &lamh, const int iz, const int ir);
	void _setbladevelc(Matrix1<double> &_ut, Matrix1<double> &_up, double &ia, const double &it, const double &b, const double &db);
	void _setbladevelc(Matrix1<float> &_ut, Matrix1<float> &_up, float &ia, const float &it, const float &b, const float &db);
	void _setbladevelc(double & _ut, double &_up, const double &b, const double &db, const double &lamh, const double &ia, const double &ir);

	void _setairfm_sp(double f[3], double m[3]);

	double _limitaz(double);
	double _limitaoa(double);

	double _liftsurfacecorr(double az, double c, double rc, double ri0, double rj0, double rk0, double ri1, double rj1, double rk1);

	template <class _Ty>
	void _aerodynacoef(Matrix1<_Ty> &_cl, Matrix1<_Ty> &_cd, Matrix1<_Ty> &incidn, Matrix1<_Ty> &ma_n);

	template <class _Ty>
	void _aerodynacoef(Matrix2<_Ty> &_cl, Matrix2<_Ty> &_cd, Matrix2<_Ty> &incidn, Matrix2<_Ty> &ma_n);

private:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3]; // hub coordinate
	myTYPE velw[3], dvelw[3], omgw[3], domgw[3]; // hub wind coordinate
	myTYPE mul, betawind;
	myTYPE beta[3];
	myTYPE lambdi_ag, lambdh_ag, lambdt_ag; // NOTE: these three variations defined at different coordinates
	Matrix2<myTYPE> bflap, dbflap, sfth;
	Matrix2<myTYPE> ut, un, up, ua, ma_n;
	Matrix1<myTYPE> _cl, _cd, _ua, _incidn, _inflow, _factor;
	Matrix2<myTYPE> incidn, inflow, cl, cd, dt;
	Matrix2<myTYPE> lambdi, lambdh, lambdt, lambdy, lambdx;

	Matrix2<myTYPE> tipstr, rotstr, shdstr, trlstr, cirlb;
	Matrix3<myTYPE> bladedeform, tipgeometry;

	myTYPE twistv, chv;

	//Matrix2<myTYPE> lambdi, lambdh, lambdt, lambdx, lambdy;
	//myTYPE beta[3];
	
public:
	CompType type;
	bool teeter;
	HingeType hingetype;
	Coordinate refcoord;
	int nb, nf, ns, ni;
	myTYPE eflap, khub, del, pitchroot, radius, bt, rroot, disk_A;
	myTYPE precone, omega0, omega, vtipa;
	myTYPE sigma, gama, a0, del0, del2, alpha0;
	Matrix1<myTYPE> chord, twist, sweep;
	myTYPE iflap, m1;
	Coordinate hubfxcoord, hubrtcoord, bladecoord, tppcoord, windcoord;
	Matrix2<myTYPE> cltc, cdtc, cmtc; //
	Matrix2<myTYPE> rastation, ristation, azstation; //, chord , twist, sweep;
	myTYPE t0;
	myTYPE velh[3], omgh[3], dvelh[3], domgh[3]; // copter states

	//GenArf Gaf;
	BladeSolver bld;
	Ambience amb;
	bool si_unit;
	AirfoilAero airfoil;
	//LBStall lbstall;
	std::vector<LBStall> lbstall;
	bool secLB;

	myTYPE sita[3];
	AeroDynType adyna;
	
	std::vector<Wake> wakev;

	LSCorr lscorr;


	int kwrot, nk, nbn, naf, nnr;
	double kwtip;
	//int nf, ns, ni;
	myTYPE rtip, rc0, outboard;
	bool outputWake, outputDisk;
	bool haveGeo, haveStr;
	bool dstall;
	myTYPE power, torque, power_i, torque_i, power_o, torque_o, power_f, torque_f, power_c, torque_c;
	myTYPE power_iid, torque_iid;
	int niter_w, niter_a;
	double Errb2, Errw2;
	double FT, KLT, Inter0, Inter1, KA;
	Monitor monitor;
};

class ModelCase
{
public:
	ModelCase() { ; }
	ModelCase(const ModelCase &M) { ; }
	~ModelCase() { ; }
	virtual void GetModel();
	virtual void GetProb();
	virtual void GetProb(double);
	virtual void GetProb(int, double, SimType);
	virtual void GetProb(int, double, SimType, bool);
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
	bool si_unit;
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
	void GetErr(myTYPE e[3]);

	void MRtoTR_Disturb(void);

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
	ModelCase model;
	//CopterSolver CSolver;

public:
	Ambience amb;
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
	int Niter;
	myTYPE controls[6];
	myTYPE sum_a1_del, sum_a2_del, max_c_del;
	Matrix2<double> AMatrix, BMatrix, AMatrixRe;
	Matrix2<double> ALongM, ALateM, ALongMRe, ALateMRe;
	Matrix1<double> dXdvel, dYdvel, dZdvel, dXdomg, dYdomg, dZdomg, dXdeul, dYdeul, dZdeul;
	Matrix1<double> dMdvel, dNdvel, dLdvel, dMdomg, dNdomg, dLdomg, dMdeul, dNdeul, dLdeul;
	Matrix1<double> dXdctrl, dYdctrl, dZdctrl, dMdctrl, dNdctrl, dLdctrl;
	Monitor monitor;
private:
	void _Allocate(void);
	void _InitVariables(void);
};


template <class _Ty>
void _setstates(_Ty _vel[3], _Ty _omg[3], _Ty _dvel[3], _Ty _domg[3], const _Ty *vc, const _Ty *wc, const _Ty *dvc, const _Ty *dwc, const Coordinate &refcoord)
{
	// transfer needed component base coordinate to copter.
	_Ty temp_v[3], temp_dv[3], temp_ddv[3];

	for (int i = 2; i >= 0; --i) {
		_vel[i] = _omg[i] = _dvel[i] = _domg[i] = 0;
		temp_v[i] = temp_dv[i] = temp_ddv[i] = 0;
	}
	Cross(temp_v, wc, refcoord.origin);
	Cross(temp_dv, wc, _vel); // w x (w x r)
	Cross(temp_ddv, dwc, refcoord.origin); // e x r

	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			_vel[i] += refcoord.Ttransf[i][j] * (temp_v[j] + vc[j]);
			_dvel[i] += refcoord.Ttransf[i][j] * (temp_dv[j] + dvc[j] + temp_ddv[j]);
			_omg[i] += refcoord.Ttransf[i][j] * wc[j];
			_domg[i] += refcoord.Ttransf[i][j] * dwc[j];
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
	sol.allocate(nitermax);
	sol(0) = q;
	return true;
}

template<class _Ty>
bool BladeSolver::_GenArfStarter(_Ty _Qt, _Ty _q, _Ty _dq, int n)
{
	ddq = Marf / (_Qt - Carf*_dq - Karf*_q);
	dq = _dq;
	q = _q;
	sol.allocate(nitermax);
	for (int i = 0; i < n; i++)
		sol(i) = q;
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
	double _a = 0;
	double ct = 0;

	ct = -airforce[2] / amb.rho / disk_A / vtipa / vtipa;
	switch (airfoil)
	{
	case C81Table:
	case LBStallMethod:
		_cl = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn / PI*180.0, ma_n);
		_cd = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn / PI*180.0, ma_n);
		break;
	case Liftslope:
		for (int i = incidn.Nv - 1; i >= 0; i--)
		{
			_a = incidn(i);
			_cl(i) = a0*(_a - alpha0);
			_cd(i) = del0 + del2*_a*_a;
		}
		break;
	case Padfield:
		for (int i = incidn.Nv - 1; i >= 0; i--)
		{
			_a = incidn(i);
			_cl(i) = a0*(_a - alpha0);
			_cd(i) = del0 + del2*ct*ct;
		}
		break;
	}
}

template <class _Ty>
void Rotor::_aerodynacoef(Matrix2<_Ty> &_cl, Matrix2<_Ty> &_cd, Matrix2<_Ty> &incidn, Matrix2<_Ty> &ma_n)
{
	double _a = 0;
	double ct = 0;

	ct = -airforce[2] / amb.rho / disk_A / vtipa / vtipa;

	switch (airfoil)
	{
	case C81Table:
	case LBStallMethod:
		_cl = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn / PI*180.0, ma_n);
		_cd = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn / PI*180.0, ma_n);
		break;
	case Liftslope:
		for (int i = incidn.Nv - 1; i >= 0; i--)
		{
			_a = incidn(i);
			_cl(i) = a0*(_a - alpha0);
			_cd(i) = del0 + del2*_a*_a;
		}
		break;
	case Padfield:
		for (int i = incidn.Nv - 1; i >= 0; i--)
		{
			_a = incidn(i);
			_cl(i) = a0*(_a - alpha0);
			_cd(i) = del0 + del2*ct*ct;
		}
		break;
	}
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
		xctrl[4] = refcoord.euler[0];
		xctrl[5] = refcoord.euler[1];
		break;
	}
}

template <class _Ty>
void Copter::SetCtrl(_Ty *xctrl, Rotor &R)
{
	int e;
	if (R.type == Mrotor)
		e = 3;
	else
		e = 1;
	for (int i = 0; i < e; ++i)
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
		dv[0] = dv[1] = dv[2] = 0;
		dw[0] = dw[1] = dw[2] = 0;
		break;
	default:
		dv[0] = dv[1] = dv[2] = 0;
		dw[0] = dw[1] = dw[2] = 0;
		break;
	}
}
