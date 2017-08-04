#pragma once
#ifndef Component_h
#define Component_h

// header files
#include "MatrixTemplate.h"
#include "Coordinate.h"
#include <vector>
#include <memory>

#include "DebugHelper.h"



// macro definition

#define TEST_DATA
//#define OUTPUT_MODE
#define SLUG_CONST 32.174
#define UNIT_CONST (32.1522/32.174)

#define UL496
#define DISK_A (PI*11.5*11.5)
#define TEETER 1
#define SEPA_SP

#define FLIGHT_TRIM // WIND_TRIM
#ifdef FLIGHT_TRIM
#define FREEDOM 6 
#define MROTOR 0
#define FUSELAGE 1
#define WING 2
#define FIN1 3
#define FIN2 4
#define TROTOR 5
#define COMPOT 6
//enum COMPONENTS
//{
//	FUSELAGE = 0,
//	WING = 1,
//	FIN1 = 2,
//	FIN2 = 3,
//	MAIN_ROTOR = 4,
//	TAIL_ROTOR = 5,
//	COMP_BASE = 6,
//};
#endif // FLIGHT_TRIM
#ifdef WIND_TRIM
#define FREEDOM 3 
#endif // WIND_TRIM

#define DISABLE_SLD 1
#define ABSOLUTED 0

#define SITA_COLL_MAX PI/6
#define SITA_CYCL_MAX PI/12
#define EULER_MAX PI/9

#define TIP
//#define ROT
//#define NER
#define PREWAKE //FREEWAKE

#define RIGID
#define QUSICSD

#ifdef UL496
#define CL_I 40
#define CL_J 12
#define CD_I 52
#define CD_J 12
#define CM_I 52
#define CM_J 11
#define NBLADE 2
#endif // UL496

#ifdef NDEBUG // release configuration
#define DISABLE_REVISE_SIZE 1
#endif // NDEBUG

#ifdef _DEBUG // debug configuration
#define DISABLE_REVISE_SIZE 0
#endif // DISABLE_REVISE_SIZE

#define MAX_SIZE 21600



class System
{
	// Abstract class
protected:
	// performance parameters etc.
	myTYPE omega, power, torque;

public:
	System() {
		//print_cons_on_screen("System");
		omega = 0;
		power = 0;
		torque = 0;
	}

	System(const System &S) {
		//print_cpcons_on_screen("System");
		omega = S.omega;
		power = S.power;
		torque = S.torque;
	}

	~System() {
		//print_dscons_on_screen("System");
		omega = 0;
		power = 0;
		torque = 0;
	}

	virtual myTYPE SetOmega() = 0;
	virtual myTYPE SetPower() = 0;
	virtual myTYPE SetTorque() = 0;
	myTYPE GetOmemga() { return omega; }
	myTYPE GetPower() { return power; }
	myTYPE GetTorque() { return torque; }

};


class Copter :private System
{
private:
	myTYPE airforce_sigma[3], airmoment_sigma[3]; // synthesized airdynamics												  

protected:
	myTYPE vel_c[3], omg_c[3], dvel_c[3], domg_c[3];

public:
	//myTYPE sita[4];
	myTYPE rho, vsound;
	myTYPE mass, inmatx[3][3];
	Matrix2<myTYPE> inmatx_M;

	myTYPE vel_g[3], omg_g[3], dvel_g[3], domg_g[3];
	Coordinate refcoord;

private:

	myTYPE SetOmega() { return 0; }

	myTYPE SetPower() { return 0; }

	myTYPE SetTorque() { return 0; }

	template <class Type> void _setderivs(Type dv[3], Type dw[3], Type *f, Type *m, const Coordinate &C) {
		Type temp_cross_v[3] = { 0 };
		Type temp_cross_o[3] = { 0 };
		Type temp_omg[3] = { 0 };
		Type temp = 0;

		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				temp_omg[i] += inmatx[i][j] * omg_c[j];
			}
		}

//#ifdef USE_DOUBLE
//		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *inmatx, 3, omg_c, 1, 0, temp_omg, 1);
//#else
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *inmatx, 3, omg_c, 1, 0, temp_omg, 1);
//#endif // USE_DOUBLE

		Cross(temp_cross_o, omg_c, temp_omg);
		Cross(temp_cross_v, omg_c, vel_c);

		for (int i = 2; i >= 0; --i) {
			dv[i] = f[i] / (mass/SLUG_CONST) + C.Ttransf[i][2] * 32.144 - temp_cross_v[i];
			dw[i] = m[i] - temp_cross_o[i];
		}
		
		Msolver(inmatx_M.v_p, dw, 3, 1);

	}


public:

	Copter();

	Copter(const Copter &H);

	~Copter();

	virtual void SetCtrl(myTYPE *xctrl, const int n) { ; }

	virtual void GetCtrl(myTYPE *xctrl, const int n) { ; }

	void SetStates(void);

	void SetStates(myTYPE v[3], myTYPE w[3], myTYPE dv[3], myTYPE dw[3]);

	void Assemble(myTYPE f[3], myTYPE m[3]);

	void SetDerivs(myTYPE *f, myTYPE *m);

	void SetDerivs(myTYPE dv[3], myTYPE dw[3], myTYPE *f, myTYPE *m);

	void SetDerivs(void);

	//void SetDerivs(myTYPE dv[3], myTYPE dw[3]);

	void GetStates(myTYPE v[3], myTYPE w[3], myTYPE dv[3], myTYPE dw[3]) {
		for (int i = 2; i >= 0; --i) {
			v[i] = vel_c[i];
			w[i] = omg_c[i];
			dv[i] = dvel_c[i];
			dw[i] = domg_c[i];
		}
	}

	void functest(void);
};




class Component :public Copter
{
private:
	myTYPE airforce_sigma[3], airmoment_sigma[3]; // synthesized airdynamics												  


protected:
	myTYPE airforce[3], airmoment[3];
	myTYPE airforce_cg[3], airmoment_cg[3];
	myTYPE vel[3], omg[3], dvel[3], domg[3];

public:

	Component();

	Component(const Component &C);

	~Component();

	void SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc);

	virtual void SetAirfm(void) { ; }

	virtual void SetCoordBase(void) { ; }

	void SetAirfm_cg(const Coordinate *base);


	void GetAirfm(myTYPE f[3], myTYPE m[3]) {
		for (int i = 2; i >= 0; --i) {
			f[i] = airforce[i];
			m[i] = airmoment[i];
		}
	}


	void GetAirfm_cg(myTYPE f[3], myTYPE m[3]) {
		for (int i = 2; i >= 0; --i) {
			f[i] = airforce_cg[i];
			m[i] = airmoment_cg[i];
		}
	}


	void GetStates(myTYPE v[3], myTYPE w[3], myTYPE dv[3], myTYPE dw[3]) {
		for (int i = 2; i >= 0; --i) {
			v[i] = vel[i];
			w[i] = omg[i];
			dv[i] = dvel[i];
			dw[i] = domg[i];
		}
	}


	void GetAirfm_sg(myTYPE f[3], myTYPE m[3]) {
		for (int i = 2; i >= 0; --i) {
			f[i] = airforce_sigma[i];
			m[i] = airmoment_sigma[i];
		}
	}


	void Assemble(const std::vector<std::unique_ptr<Component>> & C, const Coordinate *base);



	void functest(std::vector<Component> &C);


private:

	template <class Type> void _assemble(Type f[3], Type m[3], const std::vector<std::unique_ptr<Component>> & C, const Coordinate *base) {
		Type ftemp[3] = { 0,0,0 };
		Type mtemp[3] = { 0,0,0 };
		f[0] = f[1] = f[2] = 0;
		m[0] = m[1] = m[2] = 0;
		for (int i = C.size() - 1; i >= 0; --i) {
			//cout << C[i].refcoord.base << endl;
			C[i]->SetAirfm_cg(base);
			C[i]->GetAirfm_cg(ftemp, mtemp);

#ifdef OUTPUT_MODE
			printf("Component %d aerodynamics at CG: \n", i);
			printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
			printf("M: %f, %f, %f \n", mtemp[0], mtemp[1], mtemp[2]);
			cout << endl << endl;
#endif // OUTPUT_MODE

			for (int j = 0; j < 3; ++j) {
				f[j] += ftemp[j];
				m[j] += mtemp[j];
			}
		}
	}


	template <class Type> void _setstates(Type _vel[3], Type _omg[3], Type _dvel[3], Type _domg[3], const Type *vc, const Type *wc, const Type *dvc, const Type *dwc) 
	{
		// transfer needed component base coordinate to copter.

		//Type temp_a1[3] = { 0,0,0 };
		Type temp_v[3], temp_w[3], temp_dv[3], temp_dw[3];
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

//		for (int i = 2; i >= 0; --i) {
//			temp_w[i] = wc[i];
//			temp_dw[i] = dwc[i]; // simplify $e = e_e + e_r + w_e x w_r$ when w_r parallel with w_e, and e_r = 0;
//			temp_dv[i] += dvc[i] + _dvel[i];
//		}
//#ifdef USE_DOUBLE
//		vdAdd(3, _vel, vc, temp_v);
//		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, temp_v, 1, 0, _vel, 1); // cannot set b and y THE same variable
//		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, temp_w, 1, 0, _omg, 1);
//		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, temp_dv, 1, 0, _dvel, 1);
//		cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, temp_dw, 1, 0, _domg, 1);
//
//#else
//		vsAdd(3, vel, Copter::vel_c, vel);
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, temp_v, 1, 0, _vel, 1);
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, temp_w, 1, 0, _omg, 1);
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, temp_dv, 1, 0, _dvel, 1);
//		cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, temp_dw, 1, 0, _domg, 1);
//#endif // USE_DOUBLE
	}


};


class Fuselage :public Component
{
private:
	myTYPE dragA;

public:
	Fuselage();

	Fuselage(const Fuselage &F);

	~Fuselage();

	// member functions
	inline void SetAirfm(void);


};


class Wing :public Component
{
private:
#ifdef UL496
	char *type;
	myTYPE a0, cd0, cd1, cd2;
	myTYPE span, chord, taper;
#endif // UL496


public:

	Wing(char *s, int num = 1);

	Wing(const Wing &W);

	~Wing();

	// member functions
	inline void SetAirfm(void);
};


class Rotor :public Component
{
private:
	char *type;
	int kwtip, kwrot, nk, nf, ns, ni, nbn, naf, nnr;
	myTYPE eflap, khub, del, pitchroot, radius, bt, rroot;
	myTYPE precone, omega;
	myTYPE sigma, gama, a0;
	Matrix2<myTYPE> rastation, ristation, azstation; //, chord , twist, sweep;
	Matrix1<myTYPE> chord, twist, sweep;
	myTYPE iflap, m1, rtip, rc0, outboard;
	myTYPE mul, vtipa;
	myTYPE beta[3];
	myTYPE sita[4];

	Matrix2<myTYPE> cltc, cdtc, cmtc; //
	myTYPE lambdi_ag, lambdh_ag, lambdt_ag; // NOTE: these three variations defined at different coordinates
	myTYPE power, torque;
	Matrix2<myTYPE> bflap, dbflap, sfth;
	Matrix2<myTYPE> ut, un, up, ua, ma_n;
	Matrix2<myTYPE> incidn, cl, cd, cirlb;
	Matrix2<myTYPE> lambdi, lambdh, lambdt, lambdx, lambdy;
	Matrix2<myTYPE> tipstr, rotstr, shdstr, trlstr;
	//Matrix1<myTYPE> veltpp;
	Matrix3<myTYPE> bladedeform, tipgeometry;
	Coordinate hubfxcoord, hubrtcoord, bladecoord, tppcoord;


public:
	//myTYPE sita[3];
	//Coordinate hubfxcoord, hubrtcoord, bladecoord, tppcoord;



public:

	Rotor(const char *type);

	Rotor(const Rotor &A);

	~Rotor();


private:


	template <class Type> void _setbladevelc(Matrix1<Type> &_ut, Matrix1<Type> &_up, Type &ia, const Type &it, const Type &b, const Type &db) {

		Matrix1<Type> _az(ns), _ra(ns), _lambdh(ns);
		Matrix1<int> id_ns(ns), id_nf(nf);
		id_ns = step(0, ns - 1);
		id_nf = step(0, nf - 1);

		ia = it * omega;
		for (;;) {
			if (ia > 2 * PI) ia -= 2 * PI;
			else if (ia < -2 * PI) ia += 2 * PI;
			else { break; }
		}
		_az.setvalue(ia);
		_ra = rastation(0, id_ns);

		_lambdh = lambdh.interplinear_fast(azstation(id_nf, 0), _ra, _az, _ra);
#ifdef OUTPUT_MODE
		_lambdh.outputs("_lambdh.output", 4);
#endif // OUTPUT_MODE


		// air velocity
		_ut = (_ra - eflap) * cos(b) + eflap + sin(ia) * vel[0] / vtipa;
		//_ut *= mcos(sweep(0, id_ns));
		_ut *= mcos(sweep);

		_up = (_ra - eflap) * db / omega + cos(ia) * sin(b) * vel[0] / vtipa + _lambdh * cos(b);


	}


	template <class Type> void _setairfm(Matrix1<Type> &_dfx, Matrix1<Type> &_dfz, Matrix1<Type> &_dfr, const Type &it, const Type &b, const Type &db) {
		Matrix1<Type> _ut(ns), _up(ns), _ua(ns), _cl(ns), _cd(ns);
		Matrix1<Type> _factor(ns), _inflow(ns), _incidn(ns);// , _sfth(ns), _az(ns);
		Type _sfth = 0;
		Type ia = 3 * PI;
		Matrix1<int> id_ns(ns);

		id_ns = step(0, ns - 1);

		_setbladevelc(_ut, _up, ia, it, b, db);
		_up = _up*vtipa;
		//_ut *= mcos(sweep(0, id_ns));
		_ut *= mcos(sweep);
		_ut = _ut*vtipa;

#ifdef OUTPUT_MODE
		_up.outputs("_up.output", 4);
		_ut.outputs("_ut.output", 4);
#endif // OUTPUT_MODE

		_ua = msqrt(_ut*_ut + _up*_up);
		_ua = _ua / vsound;
		_ua = mmax(_ua, 0);
		_ua = mmin(_ua, 1);

		// blade pitch
		if (ia > 2 * PI || ia < -2 * PI) {
			ia = omega*it;
			for (;;) {
				if (ia > 2 * PI) { ia -= 2 * PI; }
				else if (ia < -2 * PI) { ia += 2 * PI; }
				else { break; }
			}
		}

		_sfth = sita[0] + sita[1] * cos(ia) + sita[2] * sin(ia);

		// aoa
		_inflow = atan2(_up, _ut);
		//_incidn = _sfth + twist(0, id_ns) - _inflow;
		_incidn = twist - _inflow + _sfth;
		for (int i = ns - 1; i >= 0; --i) {
			for (;;) {
				if (_incidn.v_p[i] > PI) { _incidn.v_p[i] -= 2 * PI; }
				else if (_incidn.v_p[i] < -PI) { _incidn.v_p[i] += 2 * PI; }
				else { break; }
			}
			
		}

		// air coefficients
		_aerodynacoef(_cl, _cd, _incidn, _ua);

#ifdef OUTPUT_MODE

		//cout << "up" << endl;
		//_up.output(4);
		//cout << "ut" << endl;
		//_ut.output(4);

		//cout << "cl" << endl;
		//_cl.output(4);
		_cl.outputs("_cl.output", 4);
		//cout << "cd" << endl;
		//_cd.output(4);
		_cd.outputs("_cd.output", 4);

		//cout << "_inflow" << endl;
		//_inflow.output(4);
		_inflow.outputs("_inflow.output", 4);

		//_sfth.output("_sfth.output", 4);
		//twist.output("_twist.output", 4);

		//cout << "_sfth" << endl;
		//_sfth.output(4);

		//cout << "twist" << endl;
		//twist(0, id_ns).output(4);

		//cout << "_incidn" << endl;
		//_incidn.output(4);
		_incidn.outputs("_incidn.output", 4);

		//cout << "_ua" << endl;
		//_ua.output(4);
		_ua.outputs("_ma.output", 4);
#endif // OUTPUT_MODE

		//cout << "*********************************************************" << endl;

		// air forces and moments
		//_dr.setvalue(rastation(0, 1) - rastation(0, 0));
		//_factor = _ua*_ua*chord(0, id_ns) * _dr * radius *rho*vsound*vsound*0.5;
		//_factor = _ua*_ua*chord * _dr *vsound*vsound*0.5*rho* radius;
		_factor = _ua*_ua;
#ifdef SEPA_SP
		_factor = _factor * (chord(0) * (rastation(0, 1) - rastation(0, 0)) *vsound*vsound*0.5*rho* radius);
#else
		Matrix1<Type> _dr(ns);
		_dr.setvalue(rastation(0, 1) - rastation(0, 0));
		_factor *= chord * _dr *vsound*vsound*0.5*rho* radius;
#endif // SEPA_SP


		//_factor.output(4);
		_cl *= _factor;
		_cd *= _factor;

#ifdef OUTPUT_MODE
		_cl.outputs("_dL.output", 4);
		_cd.outputs("_dD.output", 4);
#endif // OUTPUT_MODE

		_dfx = _cl * msin(_inflow) + _cd * mcos(_inflow);
		_dfz = _cl * mcos(_inflow) - _cd * msin(_inflow);
		_dfr = _dfz * (sin(b)) * (-1);

	}


	template <class Type> void _setairfm_sp(Type f[3], Type m[3]) {
		Type it = 0;
		Type ia = 0;
		Type b = 0;
		Type db = 0;
		Matrix1<Type> _dt(ns), _yf(ns), _hf(ns);
		Matrix1<Type> _dfx(ns), _dfz(ns), _dfr(ns), ra(ns), ra2(ns), ra1(ns);
		Matrix1<int> id_ns = step(0, ns - 1);
		ra = rastation(0, id_ns);

		f[0] = f[1] = f[2] = 0.0;
		m[0] = m[1] = m[2] = 0.0;

		for (int i = 0; i < nf; ++i) {
			ia = i * 2 * PI / nf;
			it = ia / omega;
			b = beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
			db = (-beta[1] * sin(ia) + beta[2] * cos(ia)) * omega;

			_setairfm(_dfx, _dfz, _dfr, it, b, db);
#ifdef OUTPUT_MODE
			_dfx.outputs("_dfx.output", 4);
			_dfz.outputs("_dfz.output", 4);
			_dfr.outputs("_dfr.output", 4);
#endif // OUTPUT_MODE

			_dt = _dfz*cos(b);
			_yf = _dfx*cos(ia) - _dfr*sin(ia);
			_hf = _dfx*sin(ia) + _dfr*cos(ia);

			f[0] -= _hf.sum()*NBLADE / nf;
			f[1] -= _yf.sum()*NBLADE / nf;
			f[2] -= _dt.sum()*NBLADE / nf;

#ifdef OUTPUT_MODE
			cout << i << endl;
			cout << "F[2]: " << f[2] << endl << endl;
#endif // OUTPUT_MODE

			ra2 = (ra - eflap)*sin(b);
			ra1 = (ra - eflap)*cos(b) + eflap;
#ifdef TEETER
			m[0] = m[1] = 0.0;
#else
			m[0] -= (_dt * ra1 * sin(ia) + _yf * ra2).sum()*radius*NBLADE / nf;
			m[1] -= (_dt * ra1 * cos(ia) + _hf * ra2).sum()*radius*NBLADE / nf;
#endif // TEETER


			m[2] += (_dfx * ra1).sum()*radius*NBLADE / nf;
		}
	}


	template<class Type> void GenArf(Matrix1<Type> &sol, Type b[3], int &niter, const Type &m, const Type &c, const Type &k, const Type &q0, const Type &dq0, const Type &dff, int Nitermax, Type pho, Type err)
	{
		Type af, am, bt, r, qq, kk, dt, temp;
		Type ck, c0, c1, c2, c3, c4, c5, _sum, c_temp, s_temp;
		Type qt = 0;
		Type qt0 = 0;
		Type ddq = 0;
		Type ddq2 = 0;
		Type dq = 0;
		Type q = 0;
		Type t = 0;
		Matrix1<Type> dfx(ns), dfr(ns), dfz(ns), ra(ns);
		Matrix1<int> id_ns;
#ifdef TEETER
		Matrix1<Type> dfx2(ns), dfz2(ns), dfr2(ns);
#endif // TEETER

		//Matrix1<Type> id_nf, az;
		int nperd = 0;
		_sum = 0;
		temp = 0;

		dt = dff / 180 * PI / omega;
		nperd = 360 / dff;
		//id_nf.allocate(nperd);
		//az.allocate(nperd);
		//step(0, nperd - 1, id_nf);
		id_ns = step(0, ns - 1);
		ra = rastation(0, id_ns);

		af = pho / (pho + 1);
		am = (2 * pho - 1) / (pho + 1);
		bt = (1 - am + af)*(1 - am + af) / 4;
		r = 0.5 - am + af;

		ck = 1 - af;
		c0 = (1 - am) / bt / dt / dt;
		c1 = ck*r / bt / dt;
		c2 = dt*c0;
		c3 = c2*dt / 2 - 1;
		c4 = ck*r / bt - 1;
		c5 = ck*(r / 2 / bt - 1)*dt;

		_setairfm(dfx, dfz, dfr, 0.0, q0, dq0);
		qt0 = 4 * iflap * omega*(omg[0] * cos(0) - omg[1] * sin(0)) + 2 * iflap*(domg[0] * sin(0) + domg[1] * cos(0));
#ifdef TEETER
		//Matrix1<Type> dfx2(ns), dfz2(ns), dfr2(ns);
		_setairfm(dfx2, dfz2, dfr2, PI / omega, -q0, -dq0);
		//dfx -= dfx2;
		//dfz -= dfz2;
		//dfr -= dfr2;
		//qt0 += ((dfz - dfz2)*ra).sum()*radius;
		qt0 += ((dfz - dfz2)*ra).sum()*radius;
#else
		qt0 += (dfz*ra).sum()*radius;
#endif // TEETER


		ddq = -m / (c*dq0 + k*q0 - qt0);
		dq = dq0;
		q = q0;
		qq = 0;
		kk = ck*k + c0*m + c1*c;
		sol.v_p[0] = q;
		for (int i = 0; i < Nitermax - 1; ++i) {

			t = i*dt + (1 - af)*dt;

			qt = 4 * iflap * omega*(omg[0] * cos(omega*t) - omg[1] * sin(omega*t)) + 2 * iflap*(domg[0] * sin(omega*t) + domg[1] * cos(omega*t));

			_setairfm(dfx, dfz, dfr, t, q, dq);
#ifdef TEETER
			//Matrix1<Type> dfx2(ns), dfz2(ns), dfr2(ns);
			_setairfm(dfx2, dfz2, dfr2, t + PI / omega, -q, -dq);
			qt += ((dfz - dfz2)*ra).sum()*radius;

#else

			qt += (dfz*ra).sum()*radius;
#endif // TEETER


			qq = qt - af*k*q + m*(c0*q + c2*dq + c3*ddq) + c*(c1*q + c4*dq + c5*ddq);
			//kk = ck*k + c0*m + c1*c;
			q = qq / kk;
			sol(i + 1) = q;
			ddq2 = (q - sol.v_p[i]) / bt / dt / dt - dq / bt / dt - (1 / bt / 2 - 1)*ddq;
			dq += (1 - r)*ddq*dt + r*ddq2*dt;
			ddq = ddq2;

			niter = i + 1;

			// exit condition
			if (niter >= 2 * nperd - 1) {
				_sum = 0;
				for (int j = 0; j < nperd; ++j) {
					temp = sol.v_p[niter - nperd + 1 + j] - sol.v_p[niter - 2 * nperd + 1 + j];
					_sum += temp*temp;
				}
				if (_sum < err*err) { break; }
			}
		}

		// fixed
		//az = id_nf*dff + omega*(dt + t) - (nperd - 1)*df;
		if (TEETER) {
			//Type val1, val2, val3;
			//val1 = val2 = val3 = 0;
			int val1 = niter%nperd;
			int val2 = nperd / 4;
			b[0] = precone;
			beta[1] = sol(niter - val1);
			if (val1 < val2) { beta[2] = -sol(niter - val1 - val2); }
			else { beta[2] = sol(niter - val1 + val2); }

			/*b[1] = b[2] = 0;
			for (int j = nperd-1; j >= 2; --j) {
			val3 = sol.v_p[j];
			val2 = sol.v_p[j - 1];
			val1 = sol.v_p[j - 2];


			if (val1 > 0 && val2 < 0) { beta[1] = -(val1 + val2) / 2; }
			if (val1 == 0 && val2 > 0) { beta[1] = val1; }
			if (val1 == 0 && val2 < 0) { beta[1] = -val1; }
			if (val1 < 0 && val2>0) { beta[1] = (val1 + val2) / 2; }

			if (val1<val2 && val2>val3) { beta[2] = val2; }
			if (val1 > val2 && val2 < val3) { beta[2] = -val2; }

			}*/
		}
		else { print_wrong_msg("Undefined when solving flap."); }


		/*c_temp = 0;
		s_temp = 0;
		for (int j = 0; j < nperd; ++j) {
		c_temp += cos(omega*(dt + t) + (j - nperd + 1)*dff/180*PI/omega) * sol(niter + j - nperd + 1);
		s_temp += sin(omega*(dt + t) + (j - nperd + 1)*dff/180*PI/omega) * sol(niter + j - nperd + 1);
		}
		if(TEETER){ b[0] = precone; }
		else { print_wrong_msg("Undefined when solving flap."); }

		b[1] = 2 * c_temp;
		b[2] = 2 * s_temp;*/
	}
	
	
	template<class Type> void _elementKMgenerete(Matrix3<Type> &mele, Matrix3<Type> &kele, const Matrix1<Type> &elegrid, const int dof) {
		//mele.allocate(eleid.Nv, dof, dof);
		//kele.allocate(eleid.Nv, dof, dof);
		//Matrix1<Type> rho_blade(ni);

		//for (int i = elegrid.Nv - 1; i >= 0; --i) {
		//	mele(i, 1, 1) = rho_blade.interplinear_fast(ristation, elegrid(i));
		//	mele(i, 2, 2) = mele(i, 1, 1);
		//	mele(i,3,3) = 
		//}
		;
	}





	void _aerodynacoef(Matrix1<myTYPE> &_cl, Matrix1<myTYPE> &_cd, Matrix1<myTYPE> &incidn, Matrix1<myTYPE> &ma_n);

	void _aerodynacoef(Matrix2<myTYPE> &_cl, Matrix2<myTYPE> &_cd, Matrix2<myTYPE> &incidn, Matrix2<myTYPE> &ma_n);

	void _flapmotionrt(void);


	void _bladeCSD(void);

	void _setairfm_pw(void);

	void _setairfm_fw(void);


public:

	void SetCtrl(myTYPE *xctrl, const int n);

	void GetCtrl(myTYPE *xctrl, const int n);

	void SetCoordBase(void);

	void BladeDynamics(void);

	void AvrgInducedVel(void);

	inline void SetAirfm(void);

	void functest(void);

};


#endif // !Component_h

