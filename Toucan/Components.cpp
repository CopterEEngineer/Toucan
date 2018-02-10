#include "stdafx.h"
#include "Components.h"

void Fuselage::SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc)
{
	_setstates(vel, omg, dvel, domg, vc, wc, dvc, dwc, refcoord);
}

void Fuselage::SetAirfm_cg(const Coordinate *base)
{
	if (!_setairfm_cg(airforce_cg, airmoment_cg, airforce, airmoment, refcoord, base))
	{
		printf("Fuselage SetAirfm_cg(): Wrong AirFM in GC. \n");
		airforce_cg[0] = airforce_cg[1] =  airforce_cg[2] = 0;
		airmoment_cg[0] =  airmoment_cg[1] =  airmoment_cg[2] = 0;
		//system("pause");
	}
}

void Fuselage::SetAirfm(void)
{
	double _alpha = Atan2(-vel[2], vel[0]); //在Atan2前加负号表示机身正迎角定义不同
	double _beta = Atan2(vel[1], vel[0]);
	double v2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];

	switch (fmdling)
	{
	case Parasite:
		airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
		airforce[1] = airforce[2] = 0;
		airmoment[0] = airmoment[1] = airmoment[2] = 0;
		break;
	case Fitting:
		airforce[0] = xf0 + (xf1 + (xf2 + xf3*_alpha)*_alpha)*_alpha;
		airforce[0] /= (Vtest*Vtest);
		airforce[0] *= v2;

		airforce[1] = yf0 + (yf1 + (yf2 + yf3*_beta)*_beta)*_beta;
		airforce[1] /= (Vtest*Vtest);
		airforce[1] *= v2;

		airforce[2] = zf0 + (zf1 + (zf2 + zf3*_alpha)*_alpha)*_alpha;
		airforce[2] /= (Vtest*Vtest);
		airforce[2] *= v2;

		airmoment[0] = 0;
		airmoment[1] = mf0 + (mf1 + (mf2 + mf3*_alpha)*_alpha)*_alpha;
		airmoment[1] /= (Vtest*Vtest);
		airmoment[1] *= v2;

		airmoment[2] = nf0 + (nf1 + (nf2 + nf3*_beta)*_beta)*_beta;
		airmoment[2] /= (Vtest*Vtest);
		airmoment[2] *= v2;
		break;
	default:
		airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
		airforce[1] = airforce[2] = 0;
		airmoment[0] = airmoment[1] = airmoment[2] = 0;
		break;
	}

	//cout << _alpha << endl;
#ifdef TEST_MODE
	//cout << endl;
	//printf("Fuselage airdynamics: \n");
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // OUTPUT_MODE
}

void Fuselage::SetAirfm(double _l)
{
	// AOA is positive when chord is above freestream
	// slipe angle is positive when velocity vector is at the right.
	double _alpha = atan2(vel[2]+_l, vel[0]);
	double _beta = atan2(vel[1], vel[0]);
	double v2 = vel[0] * vel[0] + vel[1] * vel[1] + (vel[2] + _l) * (vel[2] + _l);
	double _cx, _cy, _cz, _cm, _cn;
	//cout << DEG(_alpha) << endl;
	switch (fmdling)
	{
	case Parasite:
		airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
		airforce[1] = airforce[2] = 0;
		airmoment[0] = airmoment[1] = airmoment[2] = 0;
		break;
	case Fitting:
		airmoment[0] = 0;
		if (Abs(_alpha) <= PI / 9)
		{
			airforce[0] = xf0 + (xf1 + (xf2 + xf3*_alpha)*_alpha)*_alpha;
			airforce[2] = zf0 + (zf1 + (zf2 + zf3*_alpha)*_alpha)*_alpha;
			airmoment[1] = mf0 + (mf1 + (mf2 + mf3*_alpha)*_alpha)*_alpha;
			airforce[0] *= (v2 * Krho / Vtest / Vtest); 
			airforce[2] *= (v2 * Krho / Vtest / Vtest);
			airmoment[1] *= (v2 * Krho / Vtest / Vtest);
		}
		else
		{
			_cx = cxtc(step(0, cxtc.NI - 1), 1).interplinear_fast(cxtc(step(0, cxtc.NI - 1), 0), DEG(_alpha));
			_cz = cztc(step(0, cztc.NI - 1), 1).interplinear_fast(cztc(step(0, cztc.NI - 1), 0), DEG(_alpha));
			_cm = cmtc(step(0, cmtc.NI - 1), 1).interplinear_fast(cmtc(step(0, cmtc.NI - 1), 0), DEG(_alpha));

			airforce[0] = 0.5*amb.rho*v2*Sp*_cx;
			airforce[2] = 0.5*amb.rho*v2*Sp*_cz;
			airmoment[1] = 0.5*amb.rho*v2*Sp*Lf*_cm;
		}

		if (Abs(_beta) <= PI / 9)
		{
			airforce[1] = yf0 + (yf1 + (yf2 + yf3*_beta)*_beta)*_beta;
			airmoment[2] = nf0 + (nf1 + (nf2 + nf3*_beta)*_beta)*_beta;
			airforce[1] *= (v2 * Krho / Vtest / Vtest);
			airmoment[2] *= (v2 * Krho / Vtest / Vtest);
		}
		else
		{
			_cy = cytc(step(0, cytc.NI - 1), 1).interplinear_fast(cytc(step(0, cytc.NI - 1), 0), DEG(_beta));
			_cn = cntc(step(0, cntc.NI - 1), 1).interplinear_fast(cntc(step(0, cntc.NI - 1), 0), DEG(_beta));

			airforce[1] = 0.5*amb.rho*v2*Ss*_cy;
			airmoment[2] = 0.5*amb.rho*v2*Ss*Lf*_cn;
		}
		break;
	default:
		airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
		airforce[1] = airforce[2] = 0;
		airmoment[0] = airmoment[1] = airmoment[2] = 0;
		break;
	}
	airforce[0] *= 1.0;
#ifdef TEST_MODE
	//cout << endl;
	//printf("Fuselage airdynamics: \n");
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // OUTPUT_MODE
	monitor.Alpha = _alpha;
	monitor.Beta = _beta;
	for (int i = 0; i < 3; i++)
	{
		monitor.af[i] = airforce[i];
		monitor.mf[i] = airmoment[i];
		monitor.vel[i] = vel[i];
	}
	monitor.CY = airforce[1] / amb.rho / v2;
	monitor.CML = airmoment[0] / amb.rho / v2;
	Alpha = _alpha;
	Beta = _beta;
}

void Fuselage::SetAirfm(double _lx, double _ly, double _lz)
{
	double _alpha = atan2(vel[2] + _lz, vel[0] + _lx);
	double _beta = atan2(vel[1] - _ly, vel[0] + _lx);
	double v2 = (vel[0] + _lx)*(vel[0] + _lx) + (vel[1] - _ly)*(vel[1] - _ly) + (vel[2] + _lz)*(vel[2] + _lz);
	double _cx, _cy, _cz, _cm, _cn;

	switch (fmdling)
	{
	case Parasite:
		airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
		airforce[1] = airforce[2] = 0;
		airmoment[0] = airmoment[1] = airmoment[2] = 0;
		break;
	case Fitting:
		airmoment[0] = 0;
		if (Abs(_alpha) <= PI / 9)
		{
			airforce[0] = xf0 + (xf1 + (xf2 + xf3*_alpha)*_alpha)*_alpha;
			airforce[2] = zf0 + (zf1 + (zf2 + zf3*_alpha)*_alpha)*_alpha;
			airmoment[1] = mf0 + (mf1 + (mf2 + mf3*_alpha)*_alpha)*_alpha;
			airforce[0] *= (v2 * Krho / Vtest / Vtest);
			airforce[2] *= (v2 * Krho / Vtest / Vtest);
			airmoment[1] *= (v2 * Krho / Vtest / Vtest);
		}
		else
		{
			_cx = cxtc(step(0, cxtc.NI - 1), 1).interplinear_fast(cxtc(step(0, cxtc.NI - 1), 0), DEG(_alpha));
			_cz = cztc(step(0, cztc.NI - 1), 1).interplinear_fast(cztc(step(0, cztc.NI - 1), 0), DEG(_alpha));
			_cm = cmtc(step(0, cmtc.NI - 1), 1).interplinear_fast(cmtc(step(0, cmtc.NI - 1), 0), DEG(_alpha));

			airforce[0] = 0.5*amb.rho*v2*Sp*_cx;
			airforce[2] = 0.5*amb.rho*v2*Sp*_cz;
			airmoment[1] = 0.5*amb.rho*v2*Sp*Lf*_cm;
		}

		if (Abs(_beta) <= PI / 9)
		{
			airforce[1] = yf0 + (yf1 + (yf2 + yf3*_beta)*_beta)*_beta;
			airmoment[2] = nf0 + (nf1 + (nf2 + nf3*_beta)*_beta)*_beta;
			airforce[1] *= (v2 * Krho / Vtest / Vtest);
			airmoment[2] *= (v2 * Krho / Vtest / Vtest);
		}
		else
		{
			_cy = cytc(step(0, cytc.NI - 1), 1).interplinear_fast(cytc(step(0, cytc.NI - 1), 0), DEG(_beta));
			_cn = cntc(step(0, cntc.NI - 1), 1).interplinear_fast(cntc(step(0, cntc.NI - 1), 0), DEG(_beta));

			airforce[1] = 0.5*amb.rho*v2*Ss*_cy;
			airmoment[2] = 0.5*amb.rho*v2*Ss*Lf*_cn;
		}
		break;
	default:
		airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
		airforce[1] = airforce[2] = 0;
		airmoment[0] = airmoment[1] = airmoment[2] = 0;
		break;
	}
	airforce[0] *= 1.0;
#ifdef TEST_MODE
	//cout << endl;
	//printf("Fuselage airdynamics: \n");
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // OUTPUT_MODE
	monitor.Alpha = _alpha;
	monitor.Beta = _beta;
	for (int i = 0; i < 3; i++)
	{
		monitor.af[i] = airforce[i];
		monitor.mf[i] = airmoment[i];
		monitor.vel[i] = vel[i];
	}
	monitor.CY = airforce[1] / amb.rho / v2;
	monitor.CML = airmoment[0] / amb.rho / v2;
	Alpha = _alpha;
	Beta = _beta;
}

void Fuselage::GetAirfm(myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = airforce[i];
		m[i] = airmoment[i];
	}
}

void Fuselage::GetAirfm_cg(myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = airforce_cg[i];
		m[i] = airmoment_cg[i];
	}
}

void Fuselage::GetAngle(double &a, double &b)
{
	a = Alpha;
	b = Beta;
}

void Fuselage::GetStates(double v[3], double w[3])
{
	for (int i = 0; i < 3; i++)
	{
		v[i] = vel[i];
		w[i] = omg[i];
	}
}

void Wing::SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc)
{
	_setstates(vel, omg, dvel, domg, vc, wc, dvc, dwc, refcoord);
}

void Wing::SetAirfm_cg(const Coordinate *base)
{
	if (!_setairfm_cg(airforce_cg, airmoment_cg, airforce, airmoment, refcoord, base))
	{
		printf("Wing SetAirfm_cg(): Wrong AirFM in GC. \n");
		airforce_cg[0] = airforce_cg[1] = airforce_cg[2] = 0;
		airmoment_cg[0] = airmoment_cg[1] = airmoment_cg[2] = 0;
		//system("pause");
	}
}

void Wing::SetAirfm(void)
{
	myTYPE stot, ar, a03d;
	myTYPE aoa, cl, cd, vel2;
	myTYPE f[3] = { 0 };
	myTYPE m[3] = { 0 };

	stot = span * chord * (1 + taper) / 2;
	aoa = Atan2(-vel[2], vel[0]);
	// ignore wake disturb, makes the freestream face to wing positive. 
	// This treatmean is easily to follow blade, but, it is not agreed with the refcoord coordinate,
	// which freestream facing to wing is negative.

	switch (wmdling)
	{
	case WFitting:
		aoa += alpha0;
		cl = aoa*(a1 + aoa*aoa*(a3 + aoa*aoa*a5));
		cd = cd0 + aoa * (cd1 + aoa*cd2);
		break;
	default:
		
		ar = span*span / stot;
		a03d = a0*ar / (ar + 2 * (ar + 4) / (ar + 2));		
		cl = a03d*aoa;
		cd = cd0 + aoa * (cd1 + aoa*cd2);
		//vel2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
		
		break;
	}
	vel2 = vel[0] * vel[0] + vel[2] * vel[2];
	// freestream direction forces and moments, not refcoord
	f[0] = -0.5 * amb.rho * stot * vel2 * cd;
	f[1] = 0;
	f[2] = 0.5 * amb.rho * stot * vel2 * cl;
	m[0] = 0;
	m[1] = 0;
	m[2] = 0;
	// forces and moments at refcoord
	airforce[0] = f[0] * cos(aoa) + f[2] * sin(aoa);
	airforce[2] = -f[0] * sin(aoa) + f[2] * cos(aoa);
	airforce[1] = 0;
	airmoment[0] = m[0] * cos(aoa) + m[2] * sin(aoa);
	airmoment[1] = -m[0] * sin(aoa) + m[2] * cos(aoa);

#ifdef TEST_MODE
	//cout << endl;
	//printf("%s wing airdynamics: \n", type);
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // TEST_MODE
}

void Wing::SetAirfm(double _l)
{
	myTYPE stot, ar, a03d;
	myTYPE aoa, _aoa, cl, cd, vel2;
	myTYPE f[3] = { 0 };
	myTYPE m[3] = { 0 };

	stot = span * chord * (1 + taper) / 2;
	//_aoa = atan2(vel[2]-_l, vel[0]);
	_aoa = atan2(-vel[2]+_l, -vel[0]);

	switch (wmdling)
	{
	case WFitting:
		break;
	default:
		ar = span*span / stot;
		a1 = a0*ar / (ar + 2 * (ar + 4) / (ar + 2));
		a3 = a5 = 0;
		break;
	}

	_aoa += alpha0;
	if (_aoa > PI)
		_aoa = -1 * (2 * PI - _aoa);
	else if (_aoa < -PI)
		_aoa = 2 * PI + _aoa;

	cl = cd = 0;
	aoa = Abs(_aoa);

	if (aoa <= PI / 9)
	{
		cl = _aoa*(a1 + _aoa*_aoa*(a3 + _aoa*_aoa*a5));
	}
	else if (aoa < 3 * PI / 4)
	{
		cl = PI / 9 * (a1 + PI*PI / 81 * (a3 + PI*PI / 81 * a5));
		cl -= cl / (PI / 18 * 7)*(aoa - PI / 9);
		cl *= Sign(_aoa);
	}
	else if (aoa <= PI)
	{
		cl = PI / 9 * (a1 + PI*PI / 81 * (a3 + PI*PI / 81 * a5));
		cl -= cl / (PI / 18 * 7)*(1.5 * PI - aoa - PI / 9);
		cl *= Sign(_aoa);
	}

	if (aoa <= RAD(80))
		cd = cd0 + aoa*(cd1 + aoa*cd2);
	else if (aoa <= RAD(100))
		cd = cd0 + RAD(80)*(cd1 + RAD(80)*cd2);
	else if (aoa <= PI)
		cd = cd0 + (PI - aoa)*(cd1 + (PI - aoa)*cd2);

	vel2 = vel[0] * vel[0] + (-vel[2] + _l)*(-vel[2] + _l);
	// freestream direction forces and moments, not refcoord
	f[0] = 0.5 * amb.rho * stot * vel2 * cd;
	f[2] = 0.5 * amb.rho * stot * vel2 * cl;
	
	// forces and moments at refcoord
	_aoa -= alpha0;
	airforce[2] = f[2] * cos(_aoa) + f[0] * sin(_aoa);
	airforce[0] = f[0] * cos(_aoa) - f[2] * sin(_aoa);
	airforce[1] = 0;
	airmoment[0] = airmoment[1] = airmoment[2] = 0;
#ifdef TEST_MODE
	//cout << endl;
	//printf("%s wing airdynamics: \n", type);
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // TEST_MODE
	AOA = _aoa;
	monitor.CL = cl;
	monitor.AOA = _aoa;
	for (int i = 0; i < 3; i++)
		monitor.af[i] = airforce[i];
}

void Wing::GetAirfm(myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = airforce[i];
		m[i] = airmoment[i];
	}
}

void Wing::GetAirfm_cg(myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = airforce_cg[i];
		m[i] = airmoment_cg[i];
	}
}

void Wing::GetAngle(double &a)
{
	a = AOA;
}

void Wing::GetStates(double v[3], double w[3])
{
	for (int i = 0; i < 3; i++)
	{
		v[i] = vel[i];
		w[i] = omg[i];
	}
}

void Rotor::SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc)
{
	_setstates(vel, omg, dvel, domg, vc, wc, dvc, dwc, refcoord);

	if (type == Mrotor)
	{
		omg[0] = -wc[0];
		omg[1] = wc[1];
		omg[2] = -wc[2];
	}

	for (int i = 0; i < 3; i++)
	{
		velh[i] = vc[i];
		dvelh[i] = dvc[i];
		omgh[i] = wc[i];
		domgh[i] = dwc[i];
	}
	double euler_temp[3], omgw_temp[3], domgw_temp[3];
	euler_temp[0] = euler_temp[1] =0;
	//euler_temp[2] = -Atan2(vel[1], -vel[0]); // wind coordinate
	betawind = atan2(vel[1], -vel[0]); // 从直升机右侧来流为正侧滑角
	/*if (vel[1] > 0)
		betawind = PI - atan2(vel[1], vel[0]);
	else if (vel[1] < 0)
		betawind = -PI - atan2(vel[1], vel[0]);
	else
	{
		if (vel[0] >= 0)
			betawind = 0;
		else
			betawind = PI;
	}*/

	euler_temp[2] = -betawind;// +hubfxcoord.euler[1]; // 坐标轴转角，与坐标轴方向定义是有关的
	windcoord.SetCoordinate(euler_temp, "euler");

	_windcoordVel(velw, dvelw);
	_windcoordOmg(omgw, domgw); 	// 欧拉角速率引起
	// y方向速度引起，既偏航
	//omgw[2] = omgw_temp[2];
	//omgw[0] = omgw_temp[0] * cos(windcoord.euler[2]) + omgw_temp[1] * sin(windcoord.euler[2]);
	//omgw[1] = -omgw_temp[0] * sin(windcoord.euler[2]) + omgw_temp[1] * cos(windcoord.euler[2]);
	//domgw[2] = domgw_temp[2];
	//domgw[0] = domgw_temp[0] * cos(windcoord.euler[2]) + domgw_temp[1] * sin(windcoord.euler[2]);
	//domgw[1] = -domgw_temp[0] * sin(windcoord.euler[2]) + domgw_temp[1] * cos(windcoord.euler[2]);

	mul = -velw[0] / vtipa;
	//omega = omega0 + omg[2];
}

void Rotor::GetAirfm(myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = airforce[i];
		m[i] = airmoment[i];
	}
}

void Rotor::GetAirfm_cg(myTYPE f[3], myTYPE m[3])
{
	for (int i = 2; i >= 0; --i) {
		f[i] = airforce_cg[i];
		m[i] = airmoment_cg[i];
	}
}

void Rotor::GetBeta(myTYPE b[3])
{
	for (int i = 0; i < 3; ++i)
		b[i] = beta[i];
}

void Rotor::SetAirfm_cg(const Coordinate *base)
{
	if (!_setairfm_cg(airforce_cg, airmoment_cg, airforce, airmoment, refcoord, base))
	{
		printf("Rotor SetAirfm_cg(): Wrong AirFM in GC. \n");
		airforce_cg[0] = airforce_cg[1] = airforce_cg[2] = 0;
		airmoment_cg[0] = airmoment_cg[1] = airmoment_cg[2] = 0;
		//system("pause");
	}
}

void Copter::InitRotorCraft(ModelCase &M)
{
	model = M, simtype = M.simtype, nfree = M.nfree;
	WingV = M.WingV;
	RotorV = M.RotorV;
	fuselage = M.fuselage;
	amb = M.amb;
	mass = M.mass, inmatx_M = M.inmatx_M;
	InitCoord(M);
	//refcoord = M.refcoord, refcoord.SetBase(&refcoord);
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			inmatx[i][j] = M.inmatx[i][j];

	//_Allocate();
	_InitVariables();
	
}

void Copter::InitCoord(ModelCase &M)
{
	refcoord = M.refcoord, refcoord.SetBase(&refcoord);
	fuselage.refcoord.SetBase(&refcoord);
	for (int i = WingV.size() - 1; i >= 0; --i)
		WingV[i].refcoord.SetBase(&refcoord);
	for (int i = RotorV.size() - 1; i >= 0; --i)
		RotorV[i].SetCoordBase(&refcoord);
}

void Copter::_Allocate(void)
{
	for (int i = RotorV.size() - 1; i >= 0; --i)
		RotorV[i].Allocate();
}

void Copter::_InitVariables(void)
{
	fuselage.InitVariables();

	for (int i = RotorV.size() - 1; i >= 0; --i)
	{
		RotorV[i].Allocate();
		RotorV[i].InitVariables();
		RotorV[i].t0 = mass;
	}		
	
	for (int i = WingV.size() - 1; i >= 0; --i)
		WingV[i].InitVariables();
	

	for (int i = 0; i < 3; ++i)
	{
		airforce_sigma[i] = airmoment_sigma[i] = 0;
		vel_c[i] = omg_c[i] = dvel_c[i] = domg_c[i] = 0;
		vel_g[i] = omg_g[i] = dvel_g[i] = domg_g[i] = 0;
	}

	AMatrix.allocate(12, 12);
	AMatrixRe.allocate(12, 12);
	BMatrix.allocate(12, 12);

	ALongM.allocate(4, 4), ALateM.allocate(5, 5);
	ALongMRe.allocate(4, 4), ALateMRe.allocate(5, 5);

	dXdvel.allocate(3), dYdvel.allocate(3), dZdvel.allocate(3);
	dXdomg.allocate(3), dYdomg.allocate(3), dZdomg.allocate(3);
	dXdeul.allocate(2), dYdeul.allocate(2), dZdeul.allocate(2);
	dXdctrl.allocate(4), dYdctrl.allocate(4), dZdctrl.allocate(4);

	dLdvel.allocate(3), dMdvel.allocate(3), dNdvel.allocate(3);
	dLdomg.allocate(3), dMdomg.allocate(3), dNdomg.allocate(3);
	dLdeul.allocate(2), dMdeul.allocate(2), dNdeul.allocate(2);
	dLdctrl.allocate(4), dMdctrl.allocate(4), dNdctrl.allocate(4);
}

void Copter::SetStates(void)
{
	// will update class members
	vel_c[0] = vel_c[1] = vel_c[2] = 0;
	omg_c[0] = omg_c[1] = omg_c[2] = 0;
	// omg_g是欧拉角速率
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			vel_c[i] += refcoord.Ttransf[i][j] * vel_g[j];
			omg_c[i] += refcoord.Etransf[i][j] * omg_g[j];
		}
	}
	switch (simtype)
	{
	case FreeTrim0:
	case FreeTrim1:
	case FreeTrim2:
	case WindTrim:
		dvel_c[0] = dvel_c[1] = dvel_c[2] = 0;
		domg_c[0] = domg_c[1] = domg_c[2] = 0;
		break;
	case GeneralTrim:
		break;
	case Transient:
		dvel_c[0] = dvel_c[1] = dvel_c[2] = 0;
		domg_c[0] = domg_c[1] = domg_c[2] = 0;
		break;
	default:
		dvel_c[0] = dvel_c[1] = dvel_c[2] = 0;
		domg_c[0] = domg_c[1] = domg_c[2] = 0;
		break;
	}
	
}

void Copter::GetErr(myTYPE e[3])
{
	e[0] = sum_a1_del;
	e[1] = sum_a2_del;
	e[2] = max_c_del;
}

void Copter::PostProcess(void)
{
	ofstream OutFile("uctrl.output", std::ios::app);
	if (OutFile) {
		OutFile.precision(6);
		for (int i = 0; i < nfree; ++i) {
			OutFile << std::right << std::setw(6 + 2) << std::setprecision(6) << std::fixed << std::showpoint;
			OutFile << DEG(controls[i]) << '\t';
		}
		OutFile << endl;
		OutFile.close();
	}
	else {
		cout << "uctrl.output" << " open failed." << endl;
		system("pause");
	}
}

void Copter::functest(void)
{
	for (int i = RotorV.size() - 1; i >= 0; --i)
		RotorV[i].SetAirfm();
}

void ModelCase::GetModel(void)
{
	printf("Undefined Model. \n");
	return;
}

void ModelCase::GetProb(void)
{
	printf("Undefined Problem. \n");
	return;
}

void ModelCase::GetProb(double h)
{
	printf("Undefined Problem. \n");
	return;
}

void ModelCase::GetProb(int nf, double h, SimType st)
{
	printf("Undefined Problem. \n");
	return;
}

void ModelCase::GetProb(int nf, double h, SimType st, bool unit)
{
	printf("Undefined Problem. \n");
	return;
}

void Ambience::SetAmbience(double h)
{
	double temperature = 15 - 1.983*h / 304.8;
	rho = 0.002378*exp(-0.0297*h / 1000);
	vsound = sqrt(1.4 * 287 * (temperature + 273)) / 0.3048; // ft/s
	height = h;
}

void Airfoil::SetAirfoil(string fn)
{
	int ir, ic;
	Matrix2<double> AirfoilParas(17, 9);
	AirfoilParas.input(fn);

	ir = AirfoilParas.NI;
	ic = AirfoilParas.NJ;
	Ma.allocate(ic), CLa.allocate(ic), alpha0.allocate(ic), alpha1.allocate(ic), dalpha1.allocate(ic);
	s1.allocate(ic), s2.allocate(ic), K0.allocate(ic), K1.allocate(ic), K2.allocate(ic);
	CD0.allocate(ic), Df.allocate(ic), CNl.allocate(ic);
	Tp.allocate(ic), Tf.allocate(ic), Tv.allocate(ic), Tvl.allocate(ic);
	
	for (int i = 0; i < ic; i++)
	{
		int _ir = 0;
		Ma(i) = AirfoilParas(_ir++, i);
		CLa(i) = AirfoilParas(_ir++, i) / PI * 180;
		alpha0(i) = RAD(AirfoilParas(_ir++, i));
		alpha1(i) = RAD(AirfoilParas(_ir++, i));
		dalpha1(i) = RAD(AirfoilParas(_ir++, i));
		s1(i) = RAD(AirfoilParas(_ir++, i));
		s2(i) = RAD(AirfoilParas(_ir++, i));
		K0(i) = AirfoilParas(_ir++, i);
		K1(i) = AirfoilParas(_ir++, i);
		K2(i) = AirfoilParas(_ir++, i);
		CD0(i) = AirfoilParas(_ir++, i);
		Df(i) = AirfoilParas(_ir++, i);
		CNl(i) = AirfoilParas(_ir++, i);
		Tp(i) = AirfoilParas(_ir++, i);
		Tf(i) = AirfoilParas(_ir++, i);
		Tv(i) = AirfoilParas(_ir++, i);
		Tvl(i) = AirfoilParas(_ir++, i);
	}
	if (fn == "naca23012.txt")
	{
		alphads0 = RAD(17.91);
		alphass = RAD(14.85);
		r0 = 0.01;
		Talpha = 3.97;
		alphamin0 = RAD(21.15);
		Tr = 3.5;
	}
	else if (fn == "naca0012.txt")
	{
		alphads0 = RAD(18.73);
		alphass = RAD(14.95);
		r0 = 0.01;
		Talpha = 3.9;
		alphamin0 = RAD(16.57);
		Tr = 5.82;
	}
}

void LBStall::Prepare(Airfoil af, double c, double a, double _dt, int nf, int ns, int _nk, Matrix2<double>& aoa, Matrix2<double>& Ma, Matrix2<double>& q0, double r)
{
	airfoil = af;
	MaIn = airfoil.Ma;
	Tl = c / a;
	vsound = a;
	dt = _dt;
	aoa0M2 = aoa;
	MaM2 = Ma;
	q0M2 = q0;
	req = r;
	A1 = 0.3, b1 = 0.14, A2 = 0.7, b2 = 0.53;

	yita = 0.95;
	countk = 0;
	Nf = nf, Ns = ns, Nk = _nk;
	_allocate();
	
	tv = 0;
	_aeff = 0;
	_Tf = _Tv = 1;
	state = std::make_pair(Attach, NoVortex);
	secdID = 0;

	alphads0 = airfoil.alphads0;
	alphass = airfoil.alphass;
	Talpha = airfoil.Talpha;
	r0 = airfoil.r0;

}

void LBStall::Allocate(int nf, int ns, int _nk)
{
	Nf = nf, Ns = ns, Nk = _nk;
	_allocate();
}

void LBStall::_allocate(void)
{
	int _nk = Nk;
	CNCM2.allocate(Nf, Ns), CNIM2.allocate(Nf, Ns), CNPM2.allocate(Nf, Ns), CNaIM2.allocate(Nf, Ns), CNqIM2.allocate(Nf, Ns);
	CNfM2.allocate(Nf, Ns), CNfCM2.allocate(Nf, Ns), CDM2.allocate(Nf, Ns), CLM2.allocate(Nf, Ns);
	CNPrevisedM2.allocate(Nf, Ns), CvM2.allocate(Nf, Ns), CNvM2.allocate(Nf, Ns), CNTM2.allocate(Nf, Ns), CCfM2.allocate(Nf, Ns);
	AoALSSaveM2.allocate(Nf, Ns), AoATSSaveM2.allocate(Nf, Ns), AoAATSaveM2.allocate(Nf, Ns);
	FlowStateSaveM2.allocate(Nf, Ns), VortexStateSaveM2.allocate(Nf, Ns);
	fppM2.allocate(Nf, Ns);
	aoa0M1.allocate(Nf), MaM1.allocate(Nf), q0M1.allocate(Nf);

	alphaIn.allocate(_nk*Nf);
	x1.allocate(_nk*Nf), y1.allocate(_nk*Nf), D1.allocate(_nk*Nf), D2.allocate(_nk*Nf);
	CNC.allocate(_nk*Nf), CNI.allocate(_nk*Nf), CNP.allocate(_nk*Nf), CNaI.allocate(_nk*Nf), CNqI.allocate(_nk*Nf);
	CNf.allocate(_nk*Nf), CNfC.allocate(_nk*Nf);
	fpp.allocate(_nk*Nf), fp.allocate(_nk*Nf), Df.allocate(_nk*Nf), Dfp.allocate(_nk*Nf), alphaeff.allocate(_nk*Nf), alphaE.allocate(_nk*Nf);
	CNPrevised.allocate(_nk*Nf), CNfrevised.allocate(_nk*Nf), Dp.allocate(_nk*Nf);
	Cv.allocate(_nk*Nf), CNv.allocate(_nk*Nf);
	CNvSecd.allocate(100, _nk*Nf);
	CNT.allocate(_nk*Nf), CCf.allocate(_nk*Nf), CD.allocate(_nk*Nf), CL.allocate(_nk*Nf);
	alphak.allocate(_nk*Nf), qk.allocate(_nk*Nf), Mak.allocate(_nk*Nf);
	AoALSSave.allocate(_nk*Nf), AoATSSave.allocate(_nk*Nf), AoAATSave.allocate(_nk*Nf);
	FlowStateSave.allocate(_nk*Nf), VortexStateSave.allocate(_nk*Nf);
	Da.allocate(_nk*Nf), alphap.allocate(_nk*Nf);

	Circle.allocate(Ns);
}

void LBStall::_allocate1(void)
{
	int _nk = Nk;

	aoa0M1.allocate(Nf), MaM1.allocate(Nf), q0M1.allocate(Nf);

	alphaIn.allocate(_nk*Nf);
	x1.allocate(_nk*Nf), y1.allocate(_nk*Nf), D1.allocate(_nk*Nf), D2.allocate(_nk*Nf);
	CNC.allocate(_nk*Nf), CNI.allocate(_nk*Nf), CNP.allocate(_nk*Nf), CNaI.allocate(_nk*Nf), CNqI.allocate(_nk*Nf);
	CNf.allocate(_nk*Nf), CNfC.allocate(_nk*Nf);
	fpp.allocate(_nk*Nf), fp.allocate(_nk*Nf), Df.allocate(_nk*Nf), Dfp.allocate(_nk*Nf), alphaeff.allocate(_nk*Nf), alphaE.allocate(_nk*Nf);
	CNPrevised.allocate(_nk*Nf), CNfrevised.allocate(_nk*Nf), Dp.allocate(_nk*Nf);
	Cv.allocate(_nk*Nf), CNv.allocate(_nk*Nf);
	CNvSecd.allocate(100, _nk*Nf);
	CNT.allocate(_nk*Nf), CCf.allocate(_nk*Nf), CD.allocate(_nk*Nf), CL.allocate(_nk*Nf);
	alphak.allocate(_nk*Nf), qk.allocate(_nk*Nf), Mak.allocate(_nk*Nf);
	AoALSSave.allocate(_nk*Nf), AoATSSave.allocate(_nk*Nf), AoAATSave.allocate(_nk*Nf);
	FlowStateSave.allocate(_nk*Nf), VortexStateSave.allocate(_nk*Nf);
	Da.allocate(_nk*Nf), alphap.allocate(_nk*Nf);
}

void LBStall::SetConstants(void)
{
	A1 = 0.3, b1 = 0.14, A2 = 0.7, b2 = 0.53;
	yita = 0.95;
	eps = 1e-6;
	MaIn = airfoil.Ma;
}


void LBStall::Prepare(double c, double a, double _dt, int nf, int ns, int _nk, Matrix2<double>& aoa, Matrix2<double>& Ma, Matrix2<double>& q0)
{
	MaIn = airfoil.Ma;
	Tl = c / a;
	vsound = a;
	dt = _dt;
	aoa0M2 = aoa;
	MaM2 = Ma;
	q0M2 = q0;
	A1 = 0.3, b1 = 0.14, A2 = 0.7, b2 = 0.53;

	yita = 0.95;
	countk = 0;
	Nf = nf, Ns = ns;
	CNCM2.allocate(Nf, Ns), CNIM2.allocate(Nf, Ns), CNPM2.allocate(Nf, Ns), CNaIM2.allocate(Nf, Ns), CNqIM2.allocate(Nf, Ns);
	CNfM2.allocate(Nf, Ns), CNfCM2.allocate(Nf, Ns), CDM2.allocate(Nf, Ns), CLM2.allocate(Nf, Ns);
	CNPrevisedM2.allocate(Nf, Ns), CvM2.allocate(Nf, Ns), CNvM2.allocate(Nf, Ns), CNTM2.allocate(Nf, Ns), CCfM2.allocate(Nf, Ns);
	AoALSSaveM2.allocate(Nf, Ns), AoATSSaveM2.allocate(Nf, Ns), AoAATSaveM2.allocate(Nf, Ns);
	FlowStateSaveM2.allocate(Nf, Ns), VortexStateSaveM2.allocate(Nf, Ns);
	fppM2.allocate(Nf, Ns);

	alphaIn.allocate(_nk*Nf);
	x1.allocate(_nk*Nf), y1.allocate(_nk*Nf), D1.allocate(_nk*Nf), D2.allocate(_nk*Nf);
	CNC.allocate(_nk*Nf), CNI.allocate(_nk*Nf), CNP.allocate(_nk*Nf), CNaI.allocate(_nk*Nf), CNqI.allocate(_nk*Nf);
	CNf.allocate(_nk*Nf), CNfC.allocate(_nk*Nf);
	fpp.allocate(_nk*Nf), fp.allocate(_nk*Nf), Df.allocate(_nk*Nf), Dfp.allocate(_nk*Nf), alphaeff.allocate(_nk*Nf), alphaE.allocate(_nk*Nf);
	CNPrevised.allocate(_nk*Nf), CNfrevised.allocate(_nk*Nf), Dp.allocate(_nk*Nf);
	Cv.allocate(_nk*Nf), CNv.allocate(_nk*Nf);
	CNvSecd.allocate(100, _nk*Nf);
	CNT.allocate(_nk*Nf), CCf.allocate(_nk*Nf), CD.allocate(_nk*Nf), CL.allocate(_nk*Nf);
	alphak.allocate(_nk*Nf), qk.allocate(_nk*Nf), Mak.allocate(_nk*Nf);
	AoALSSave.allocate(_nk*Nf), AoATSSave.allocate(_nk*Nf), AoAATSave.allocate(_nk*Nf);
	FlowStateSave.allocate(_nk*Nf), VortexStateSave.allocate(_nk*Nf);
	Da.allocate(_nk*Nf), alphap.allocate(_nk*Nf);

	tv = 0;
	_aeff = 0;
	_Tf = _Tv = 1;
	state = std::make_pair(Attach, NoVortex);
	secdID = 0;

	alphads0 = airfoil.alphads0;
	alphass = airfoil.alphass;
	Talpha = airfoil.Talpha;
	r0 = airfoil.r0;
}


void LBStall::Prepare(Airfoil af, double c, double a, double _dt, int nf, int ns, int _nk, Matrix2<double>& aoa, Matrix2<double>& Ma, Matrix2<double>& q0)
{
	airfoil = af;
	Prepare(c, a, _dt, nf, ns, _nk, aoa, Ma, q0);
}


void LBStall::Starter(double c, double a, double _dt)
{
	Tl = c / a;
	dt = _dt; 
	//aoa0M2 = aoa;
	//MaM2 = ma;
	//q0M2 = aoad / MaM2*Tl;

	countk = 0;
	tv = 0;
	_aeff = 0;
	_Tf = _Tv = 1;
	state = std::make_pair(Attach, NoVortex);
	secdID = 0;

	_allocate1();

	alphads0 = airfoil.alphads0;
	alphass = airfoil.alphass;
	Talpha = airfoil.Talpha;
	r0 = airfoil.r0;
}


void LBStall::Complete(void)
{
	countk++;
	CNT(countk1) = CNf(countk1) + CNv(countk1);
	//CCf(countk1) = 0.9*yita*CNMa*(alphaE(countk1) - alpha0)*sin(alphaE(countk1) - alpha0)*(fpptemp - 1);
	CCf(countk1) = 0.9*yita*CNC(countk1)*tan(alphaE(countk1) - alpha0)*(fpptemp - 1);

	if (CNPrevised(countk1) >= CNILc)
		CCf(countk1) *= pow(fpp(countk1), 0.25*DFLc*(CNPrevised(countk1) - CNILc));

	CD(countk1) = Cd0Lc + CNT(countk1)*sin(alphaLc) - CCf(countk1)*cos(alphaLc);
	CL(countk1) = CNT(countk1)*cos(alphaLc) + CCf(countk1)*sin(alphaLc);
}

void LBStall::Complete(Matrix1<double>& cl, Matrix1<double>& cd, Matrix1<double>& cn, Matrix1<double>& cc)
{
	int iz0 = countk / Nf;
	int ic = countk;
	iz0 = countk - iz0*Nf - 1;
	if (iz0 < 0)
		iz0 += Nf;
	for (int iz = iz0; iz >= 0; iz--, ic--)
	{
		cn(iz) = CNT(ic);
		cc(iz) = CCf(ic);
		cl(iz) = CL(ic);
		cd(iz) = CD(ic);
	}
	for (int iz = Nf - 1; iz >= iz0 + 1; iz--, ic--)
	{
		cn(iz) = CNT(ic);
		cc(iz) = CCf(ic);
		cl(iz) = CL(ic);
		cd(iz) = CD(ic);
	}
}

void LBStall::AttachFlow(int iz, int ir)
{
	double _dak, _dak2, _dqk, _dqk2;
	int iz1, iz2;

	iz1 = iz - 1;
	if (iz1 < 0)
		iz1 = Nf - 1;
	iz2 = iz1 - 1;
	if (iz2 < 0)
		iz2 = Nf - 1;

	alphaLc = aoa0M2(iz, ir); // 0.5*(aoa0M2(iz, ir) + aoa0M2(iz1, ir));
	MaLc = 0.5*(MaM2(iz, ir) + MaM2(iz1, ir));
	ds = MaLc*dt * 2 / Tl;
	beta2 = 1 - MaLc * MaLc;

	CNMa = airfoil.CLa.interplinear_fast(MaIn, MaLc);
	CNILc = airfoil.CNl.interplinear_fast(MaIn, MaLc);

	// constant
	Tp = airfoil.Tp.interplinear_fast(MaIn, MaLc);
	Tf = airfoil.Tf.interplinear_fast(MaIn, MaLc);
	alpha0 = airfoil.alpha0.interplinear_fast(MaIn, MaLc);
	alpha1 = airfoil.alpha1.interplinear_fast(MaIn, MaLc);
	dalpha1 = 0;// airfoil.dalpha1.interplinear_fast(MaIn, MaLc);
	s1 = airfoil.s1.interplinear_fast(MaIn, MaLc);
	s2 = airfoil.s2.interplinear_fast(MaIn, MaLc);
	Tv = airfoil.Tv.interplinear_fast(MaIn, MaLc);
	Tvl = airfoil.Tvl.interplinear_fast(MaIn, MaLc);
	DFLc = airfoil.Df.interplinear_fast(MaIn, MaLc);
	Cd0Lc = airfoil.CD0.interplinear_fast(MaIn, MaLc);
	alphamin0 = airfoil.alphamin0;
	Tr = airfoil.Tr;

	_Tv = Tv;
	_Tf = Tf;

	countk1 = countk + 1;
	alphaIn(countk1) = DEG(alphaLc);

	_dak = aoa0M2(iz, ir) - aoa0M2(iz1, ir);
	x1(countk1) = x1(countk) * exp(-b1*beta2*ds) + A1 * _dak * exp(-b1*beta2*ds*0.5);
	y1(countk1) = y1(countk) * exp(-b2*beta2*ds) + A2 * _dak * exp(-b2*beta2*ds*0.5);
	alphaE(countk1) = alphaLc - x1(countk1) - y1(countk1);// -RAD(alpha0);
	CNC(countk1) = CNMa * alphaE(countk1);

	Ta = 0.75*Tl / (1 - MaLc + PI*beta2*MaLc*MaLc*(A1*b1 + A2*b2));
	Tq = 0.75*Tl / (1 - MaLc + 2 * PI*beta2*MaLc*MaLc*(A1*b1 + A2*b2));

	_dak2 = aoa0M2(iz1, ir) - aoa0M2(iz2, ir);
	_dqk = q0M2(iz, ir) - q0M2(iz1, ir);
	_dqk2 = q0M2(iz1, ir) - q0M2(iz2, ir);
	D1(countk1) = D1(countk)*exp(-dt / Ta) + (_dak - _dak2) / dt *exp(-dt / Ta*0.5);
	D2(countk1) = D2(countk)*exp(-dt / Tq) + (_dqk - _dqk2) / dt *exp(-dt / Tq*0.5);
	CNaI(countk1) = 4 * Ta / MaLc*(_dak / dt - D1(countk1));
	CNqI(countk1) = Tq / MaLc*(_dqk / dt - D2(countk1));
	CNI(countk1) = CNaI(countk1) + CNqI(countk1);

	CNP(countk1) = CNC(countk1) + CNI(countk1);

	req = q0M2(iz, ir)*0.5;
	// rad
	if (Abs(req) >= r0)
		alphacr = alphads0;
	else
		alphacr = alphass + (alphads0 - alphass)*Abs(req) / r0;
	alphamin = alphamin0 + req*Tr;
}

void LBStall::AttachFlow(int iz)
{
	double _dak2, _dqk2;
	int iz1, iz2;

	iz1 = iz - 1;
	if (iz1 < 0)
		iz1 = Nf - 1;
	iz2 = iz1 - 1;
	if (iz2 < 0)
		iz2 = Nf - 1;

	alphaLc = aoa0M1(iz); // 0.5*(aoa0M2(iz, ir) + aoa0M2(iz1, ir));
	MaLc = 0.5*(MaM1(iz) + MaM1(iz1));
	ds = MaLc*dt * 2 / Tl;
	beta2 = 1 - MaLc * MaLc;

	CNMa = airfoil.CLa.interplinear_fast(MaIn, MaLc);
	CNILc = airfoil.CNl.interplinear_fast(MaIn, MaLc);

	// constant
	Tp = airfoil.Tp.interplinear_fast(MaIn, MaLc);
	Tf = airfoil.Tf.interplinear_fast(MaIn, MaLc);
	alpha0 = airfoil.alpha0.interplinear_fast(MaIn, MaLc);
	alpha1 = airfoil.alpha1.interplinear_fast(MaIn, MaLc);
	dalpha1 = 0;// airfoil.dalpha1.interplinear_fast(MaIn, MaLc);
	s1 = airfoil.s1.interplinear_fast(MaIn, MaLc);
	s2 = airfoil.s2.interplinear_fast(MaIn, MaLc);
	Tv = airfoil.Tv.interplinear_fast(MaIn, MaLc);
	Tvl = airfoil.Tvl.interplinear_fast(MaIn, MaLc);
	DFLc = airfoil.Df.interplinear_fast(MaIn, MaLc);
	Cd0Lc = airfoil.CD0.interplinear_fast(MaIn, MaLc);
	alphamin0 = airfoil.alphamin0;
	Tr = airfoil.Tr;

	_Tv = Tv;
	_Tf = Tf;

	countk1 = countk + 1;
	alphaIn(countk1) = DEG(alphaLc);
	Mak(countk1) = MaLc;

	_dak = aoa0M1(iz) - aoa0M1(iz1);
	x1(countk1) = x1(countk) * exp(-b1*beta2*ds) + A1 * _dak * exp(-b1*beta2*ds*0.5);
	y1(countk1) = y1(countk) * exp(-b2*beta2*ds) + A2 * _dak * exp(-b2*beta2*ds*0.5);
	alphaE(countk1) = alphaLc-x1(countk1) - y1(countk1);
	CNC(countk1) = CNMa * (alphaE(countk1) - alpha0);

	Ta = 0.75*Tl / (1 - MaLc + PI*beta2*MaLc*MaLc*(A1*b1 + A2*b2));
	Tq = 0.75*Tl / (1 - MaLc + 2 * PI*beta2*MaLc*MaLc*(A1*b1 + A2*b2));

	_dak2 = aoa0M1(iz1) - aoa0M1(iz2);
	_dqk = q0M1(iz) - q0M1(iz1);
	_dqk2 = q0M1(iz1) - q0M1(iz2);
	D1(countk1) = D1(countk)*exp(-dt / Ta) + (_dak - _dak2) / dt *exp(-dt / Ta*0.5);
	D2(countk1) = D2(countk)*exp(-dt / Tq) + (_dqk - _dqk2) / dt *exp(-dt / Tq*0.5);
	CNaI(countk1) = 4 * Ta / MaLc*(_dak / dt - D1(countk1));
	CNqI(countk1) = Tq / MaLc*(_dqk / dt - D2(countk1));
	CNI(countk1) = CNaI(countk1) + CNqI(countk1);

	CNP(countk1) = CNC(countk1) + CNI(countk1);

	req = q0M1(iz)*0.5;
	// rad
	if (Abs(req) >= r0)
		alphacr = alphads0;
	else
		alphacr = alphass + (alphads0 - alphass)*Abs(req) / r0;
	alphamin = PI;// alphamin0 + req*Tr;
}

void LBStall::DynamicStall(void)
{
	double _dak = alphaIn(countk1) - alphaIn(countk); // deg
	Dp(countk1) = Dp(countk)*exp(-ds / Tp) + (CNP(countk1) - CNP(countk))*exp(-ds / Tp*0.5);
	Da(countk1) = Da(countk)*exp(-ds / Talpha) + _dak;// *exp(-ds / Talpha*0.5);
	CNPrevised(countk1) = CNP(countk1) - Dp(countk1);
	if (MaLc >= 0.3)
	{		
		alphaeff(countk1) = CNPrevised(countk1) / CNMa;
		_aeff = alphaeff(countk1);
	}
	else
	{		
		alphap(countk1) = RAD(alphaIn(countk1) - Da(countk1));//
		_aeff = alphaeff(countk1) = alphap(countk1);
	}
	state = stateChange(state);
	//state = stateChange(alphaLc);

	FlowStateSave(countk1) = state.first;
	VortexStateSave(countk1) = state.second;

	switch (state.first)
	{
	case Attach:
		performAttach();
		break;
	case Separate:
		performSeparate();
		break;
	case ReAttach:
		performReAttach();
		break;
	default:
		break;
	}
}


void LBStall::Save(int ir)
{
	int iz0 = countk / Nf;
	int ic = countk;

	Circle(ir) = iz0;

	iz0 = countk - iz0*Nf - 1;
	if (iz0 < 0)
		iz0 += Nf;
	for (int iz = iz0; iz >= 0; iz--, ic--)
	{
		CNTM2(iz, ir) = CNT(ic);
		CCfM2(iz, ir) = CCf(ic);
		CDM2(iz, ir) = CD(ic);
		CLM2(iz, ir) = CL(ic);
		FlowStateSaveM2(iz, ir) = FlowStateSave(ic);
		VortexStateSaveM2(iz, ir) = VortexStateSave(ic);
		fppM2(iz, ir) = fpp(ic);
	}
	for (int iz = Nf - 1; iz >= iz0 + 1; iz--, ic--)
	{
		CNTM2(iz, ir) = CNT(ic);
		CCfM2(iz, ir) = CCf(ic);
		CDM2(iz, ir) = CD(ic);
		CLM2(iz, ir) = CL(ic);
		FlowStateSaveM2(iz, ir) = FlowStateSave(ic);
		VortexStateSaveM2(iz, ir) = VortexStateSave(ic);
		fppM2(iz, ir) = fpp(ic);
	}
	/*if (ir == 1)
	{
		CNT.output("CNT.output", 10);
		CCf.output("CCf.output", 10);
		CD.output("CD.output", 10);
		CL.output("CL.output", 10);
		FlowStateSave.output("FlowState.output", 2);
		VortexStateSave.output("VortexState.output", 2);
		alphaIn.output("alphaIn.output", 2);
		fpp.output("fpp.output", 10);
	}*/
}

bool LBStall::isExit(int &ck)
{
	double _sum = 0;
	for (int i = 0; i < Nf; i++)
		_sum += pow(CNT(countk - i) - CNT(countk - Nf - i), 2);
	ck = countk;
	if (_sum < eps)
		return true;
	else
		return false;
}


void LBStall::performAttach(void)
{
	//alpha1 = airfoil.alpha1.interplinear_fast(MaIn, MaLc);
	Ts = 10000;
	tv = 0;
	_secdcomfd = false;

	if (MaLc < 0.3)
		_performKflow(4*Tf);
	else
	{
		fp(countk1) = fpp(countk1) = 1.0;
		fpptemp = 2;

		CNfC(countk1) = CNC(countk1);
		CNf(countk1) = CNfC(countk1) + CNI(countk1);
	}
	//_performKflow(Tf);
}

void LBStall::performSeparate(void)
{
	tv += ds;
	switch (state.second)
	{
	case NoVortex:
		_performKflow(Tf);
		_performdstall(Tv);
		break;
	case Conv:
		_performKflow(0.5*Tf);
		_performdstall(Tv);
		break;
	case Nearby:
		_performKflow(3*Tf); //3Tf, 0.25
		_performdstall(0.5*Tv); //0.5Tv, 3
		break;
	case Far:
		_performKflow(3*Tf);//1//3
		_performdstall(0.5*Tv);//2
		break;
	case Secd:
		if (_enableSecd()) //计算Ts
			_performsecd();
		break;
	default:
		break;
	}
}

void LBStall::performReAttach(void)
{
	tv += ds;
	switch (state.second)
	{
	case NoVortex://临界值降低至以下
		_performKflow(4*Tf);
		_performdstall(Tv);
		break;
	case Conv:
		_performKflow(0.5*Tf);
		_performdstall(Tv);
		break;
	case Nearby:
		_performKflow(3*Tf);
		_performdstall(0.5*Tv);
		break;
	case Far://涡消失以后
		_performKflow(3*Tf);
		_performdstall(0.5*Tv);
		break;
	case Secd:
		if (_enableSecd())//计算Ts
			_performsecd();
		break;
	default:
		break;
	}
}

void LBStall::_performKflow(double _tf)
{
	_Tf = _tf;

	if (Abs(_aeff) <= alpha1)
		fp(countk1) = 1 - 0.3*exp((Abs(_aeff - dalpha1) - alpha1) / s1);
	else
		fp(countk1) = 0.04 + 0.66*exp((alpha1 - Abs(_aeff - dalpha1)) / s2);

	Df(countk1) = Df(countk)*exp(-ds / _Tf) + (fp(countk1) - fp(countk))*exp(-ds / _Tf*0.5);
	//动态时相当于对K流分离点引入延迟
	fpp(countk1) = fp(countk1) - Df(countk1);
	fpptemp = 1 + sqrt(fpp(countk1));

	CNfC(countk1) = CNMa*0.25*fpptemp*fpptemp*(alphaE(countk1) - alpha0);
	CNf(countk1) = CNfC(countk1) + CNI(countk1);
}

void LBStall::_performdstall(double _tv)
{
	_Tv = _tv;

	fpptemp = 1 + sqrt(fpp(countk1));

	if (state.first == ReAttach)
	{
		switch (state.second)
		{
		case NoVortex:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = 0;
			break;
		case Conv:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = CNv(countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
			break; 
		case Nearby:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = CNv(countk)*exp(-ds / _Tv);
			break;
		case Far:
			Cv(countk1) = 0;// Cv(countk);
			CNv(countk1) = CNv(countk)*exp(-ds / _Tv);
			break;
		case Secd:
			if (_secdcomfd)
			{
				CNvSecd(secdID, countk1) = CNvSecd(secdID, countk)*exp(-ds / _Tv);
				CNv(countk1) = CNvSecd(secdID, countk1);
			}
			else
			{
				Cv(countk1) = Cv(countk);
				CNv(countk1) = CNv(countk)*exp(-ds / _Tv);
			}
		default:
			break;
		}
	}
	else
	{
		switch (state.second)
		{
		case NoVortex:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = 0; 
		case Conv:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = CNv(countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
			break;
		case Nearby:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = CNv(countk)*exp(-ds / _Tv); //确认
			break;
		case Far:
			Cv(countk1) = 0;// CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNv(countk1) = CNv(countk)*exp(-ds / _Tv);
			break;
		case Secd:
			Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
			CNvSecd(secdID, countk1) = CNvSecd(secdID, countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
			CNv(countk1) = CNvSecd(secdID, countk1);
			break; 
		default:
			break;
		}
	}
}


void LBStall::_performsecd(void)
{
	if (tv < Tvl)
	{
		_performKflow(Tf*2);
		_performdstall(Tv*0.1);
	}
	else if (tv < 2 * Tvl)
	{
		_performKflow(Tf*4);
		_performdstall(0.25*Tv);
	}
	else
	{
		//次级涡到达远场
		_performKflow(Tf*2);
		//tv = 0;
		_secdcomfd = false;
		CNv(countk1) = Cv(countk1) = 0;
	}
}


bool LBStall::_startReAttach(void)
{
	//if (_aeff * (_aeff - alphaeff(countk)) < 0 && Abs(_aeff) < alpha1 - RAD(0))
	if (Abs(_aeff) < Abs(alphaeff(countk)) && Abs(alphaLc) < Min(alpha1,alphamin) && (MaLc - 0.3)*(Mak(countk) - 0.3) >= 0)
		return true;
	else
		return false;
}

bool LBStall::_startAttach(void)
{
	//if (countk > 0 && (_aeff - alphaeff(countk)) * (alphaeff(countk) - alphaeff(countk-1)) < 0)
	if (Abs(_aeff) >= Abs(alphaeff(countk)) && (MaLc - 0.3)*(Mak(countk) - 0.3) >= 0)
		return true;
	else
		return false;
}

bool LBStall::_enableSecd(void)
{
	Ts = 1000;// 2 * (1 - fpp(countk1)) / 0.19;
	while (tv > Ts)
	{
		_secdcomfd = true;
		state.second = Secd;
		tv -= Ts;
		secdID++;
		if (secdID > 100)
			printf("Warning: Secondary Vortex Exceed Max Number in _enableSecd() \n");
	}

	if (state.second == Secd)
		return true;
	else
		return false;
}

LBStall::StatePair LBStall::stateChange(StatePair curState)
{
	switch (curState.first)
	{
	case Attach:
		if (MaLc < 0.3)
		{
			if (_ShengCriterion())
				return StatePair(Separate, Conv);
			else if (_startReAttach())
				return StatePair(ReAttach, NoVortex);
			else
				return StatePair(Attach, NoVortex);
		}
		else
		{
			if (Abs(CNPrevised(countk1)) > CNILc)
				return StatePair(Separate, Conv);
			else if (_startReAttach())
				return StatePair(ReAttach, NoVortex);
			else
				return StatePair(Attach, NoVortex);
		}
	case Separate:
		if(!_startReAttach())
		{ 
			if (Abs(CNPrevised(countk1)) < CNILc)
			{
				tv = 0;
				return StatePair(Separate, NoVortex);
			}
			if (state.second == Secd)
				return StatePair(Separate, Secd);
			if (tv + ds > Ts)
				return StatePair(Separate, Secd);
			else if (tv + ds > 2 * Tvl)
			{
				tv = 1000;
				return StatePair(Separate, Far);
			}
			else if (tv + ds > Tvl)
				return StatePair(Separate, Nearby);
			else
				return StatePair(Separate, Conv);		
		}
		else
		{
			if (Abs(CNPrevised(countk1)) < CNILc)
			{
				tv = 0;
				return StatePair(ReAttach, NoVortex);
			}
			if (state.second == Secd)
				return StatePair(ReAttach, Secd);
			if (tv + ds > Ts)
				return StatePair(ReAttach, Secd);
			else if (tv + ds > 2 * Tvl)
			{
				tv = 1000;
				return StatePair(ReAttach, Far);
			}
			else if (tv + ds > Tvl)
				return StatePair(ReAttach, Nearby);
			else
				return StatePair(ReAttach, Conv);
		}
	case ReAttach:
		if (!_startAttach())
		{
			if (Abs(CNPrevised(countk1)) < CNILc)
			{
				tv = 0;
				return StatePair(ReAttach, NoVortex);
			}
			if (state.second == Secd)
				return StatePair(ReAttach, Secd);
			if (tv + ds > Ts)
				return StatePair(ReAttach, Secd);
			else if (tv + ds > 2 * Tvl)
			{
				tv = 1000;
				return StatePair(ReAttach, Far);
			}
			else if (tv + ds > Tvl)
				return StatePair(ReAttach, Nearby);
			else
				return StatePair(ReAttach, Conv);
		}
		else
		{
			if (Abs(CNPrevised(countk1)) > CNILc)
				return StatePair(Separate, Conv);
			else
				return StatePair(Attach, NoVortex); 		
		}	
	default:
		printf("Warning: Undefined Flow State in stateChange() \n");
		return StatePair(Attach, NoVortex);
	}
}


LBStall::StatePair LBStall::stateChange(double aoa)
{	
	StatePair _state;
	double _ac = alpha1 + RAD(2.5);
	if (alphaIn(countk1) > alphaIn(countk))
	{
		_state = std::make_pair(Attach, NoVortex);
		if (aoa > _ac)
			_state = std::make_pair(Separate, Conv);
		if (aoa > _ac + RAD(2))
			_state = std::make_pair(Separate, Nearby);
	}
	else
	{
		_state = std::make_pair(Separate, Nearby);
		if (aoa < _ac + RAD(4))
			_state = std::make_pair(Separate, Far);
		if (aoa < _ac - RAD(2))
			_state = std::make_pair(ReAttach, NoVortex);
	}
	return _state;
}

bool LBStall::_ShengCriterion(void)
{
	if (Abs(alphap(countk1)) > alphacr)// + Abs(req)*Talpha)
	{
		if (Abs(req) >= r0)
			dalpha1 = (alphads0 - alphass); // rad
		else
			dalpha1 = (alphads0 - alphass)*Abs(req) / r0; // rad
		return true;
	}
	else
	{
		dalpha1 = 0;
		return false;
	}
}

void LBStall::FuncTest(void)
{
	Airfoil airfoil;
	double c = 0.1;
	double a = 340;
	double k = 0.102;
	double Ma = 0.12;
	double V = Ma*a;
	double omega = 2 * V*k / c;
	double df = 2;
	double dt = RAD(df) / omega;
	int nf = 360 / df;
	int ns = 1;
	int nk = 4;
	double aAmp = RAD(10);
	double a0 = RAD(10);
	double q0 = 0;
	Matrix2<double> _Ma(nf, ns), _aoa(nf, ns), _aoaDG(nf, ns), _q0(nf, ns);

	_Ma.setvalue(Ma);

	for (int i = 0; i < nf; i++)
	{
		_aoa(i, 0) = a0 + aAmp*sin(omega*dt*i);
		_q0(i, 0) = aAmp*omega*cos(omega*dt*i)*c / V;
		_aoaDG(i, 0) = DEG(_aoa(i, 0));
	}
	airfoil.SetAirfoil("naca23012.txt");
	airfoil.alphads0 = RAD(17.91); 
	airfoil.alphass = RAD(14.85); 
	airfoil.r0 = 0.01;
	airfoil.Talpha = 3.97; 
	airfoil.alphamin0 = RAD(21.15);
	airfoil.Tr = 3.5;
	/*airfoil.SetAirfoil("naca0012.txt");
	airfoil.alphads0 = RAD(18.73);
	airfoil.alphass = RAD(14.95);
	airfoil.r0 = 0.01;
	airfoil.Talpha = 3.9;
	airfoil.alphamin0 = RAD(16.57);
	airfoil.Tr = 5.82;*/
	
	
	Prepare(airfoil, c, a, dt, nf, ns, nk, _aoa, _Ma, _q0);

	for (int ik = 0; ik < nk; ik++)
	{
		for (int iz = 0; iz < nf; iz++)
		{
			printf("Circle No: %d\n", ik);
			if (!(iz == 0 & ik == 0))
			{
				AttachFlow(iz, 0);
				DynamicStall();
				Complete();

				printf("AOA = %f alphaE = %f\n", DEG(alphaLc), DEG(alphaE(countk1)));
				printf("Cur State: (%d, %d) \n", state.first, state.second);
				printf("Prv State: (%d, %d) \n", FlowStateSave(countk - 1), VortexStateSave(countk - 1));
				printf("Tf = %f, Tv = %f \n", _Tf, _Tv);
				printf("tv = %f, Tvl = %f, Ts = %f \n", tv, Tvl, Ts);

				printf("\n");
				printf("CNP = %f, CNI = %f \n", CNP(countk1), CNI(countk1));

				printf("\n");
				printf("CNPr = %f, CNILc = %f \n", CNPrevised(countk1), CNILc);

				printf("\n");
				printf("r = %f \n", req);
				printf("Alpha' = %f, Alpha_cr = %f, Alpha_min = %f \n", DEG(alphap(countk1)), DEG(alphacr), DEG(alphamin));

				printf("\n");
				printf("CNf = %f, CNfC = %f \n", CNf(countk1), CNfC(countk1));

				printf("fp = %f, fpp = %f \n", fp(countk1), fpp(countk1));
				printf("alpha_eff = %f, alpha_1 = %f \n", DEG(_aeff), DEG(alpha1));

				printf("\n");
				printf("CNv = %f, Cv = %f \n", CNv(countk1), Cv(countk1));

				printf("\n");
				printf("CNT = %f, CCf = %f \n", CNT(countk1), CCf(countk1));
				printf("***********************************************************\n");
			}
		}
	}

	isExit(countk);
	_aoaDG.output("AOAIn_LBD_NACA0012.output", 10);
	CNTM2.output("CNT_LBD_NACA0012.output", 10);
	CCfM2.output("CCf_LBD_NACA0012.output", 10);
	CDM2.output("CD_LBD_NACA0012.output", 10);
	CLM2.output("CL_LBD_NACA0012.output", 10);

	alphaeff.output("alphaeff.output", 10);
	FlowStateSaveM2.output("FlowState.output", 1);
	VortexStateSaveM2.output("VortexState.output", 1);

	alphaIn.output("AlphaIn.output", 10);
	//CNT.output("CNT.output", 10);
	fppM2.output("fpp.output", 10);
}

bool LBStall::Solver(int &ck, int ir, Matrix2<double> &aoa, Matrix2<double> &Ma, Matrix2<double> &aoad)
{
	for (int iz = 0; iz < Nf; iz++)
	{
		aoa0M1(iz) = aoa(iz, ir);
		MaM1(iz) = Ma(iz, ir);
		q0M1(iz) = aoad(iz, ir)*Tl / MaM1(iz);
	}

	for (int ik = 0; ik < Nk; ik++)
	{
		for (int iz = 0; iz < Nf; iz++)
		{
			if (!(iz == 0 && ik == 0))
			{
				AttachFlow(iz);
				DynamicStall();
				Complete();
				
			/*	if (ir == 1)
				{
					printf("Circle No: %d \n", ik);
					printf("AOA = %f alphaE = %f\n", DEG(alphaLc), DEG(alphaE(countk1)));
					printf("Ma = %f \n", MaLc);
					printf("Cur State: (%d, %d) \n", state.first, state.second);
					printf("Prv State: (%d, %d) \n", FlowStateSave(countk - 1), VortexStateSave(countk - 1));
					printf("Tf = %f, Tv = %f \n", _Tf, _Tv);
					printf("tv = %f, Tvl = %f, Ts = %f, Tl = %f \n", tv, Tvl, Ts, Tl);

					printf("\n");
					printf("CNP = %f, CHC = %f, CNI = %f \n", CNP(countk1), CNC(countk1), CNI(countk1));
					printf("_dak = %f, _dqk = %f \n", DEG(_dak), DEG(_dqk));
					printf("CNaI = %f, CNqI = %f \n", CNaI(countk1), CNqI(countk1));

					printf("\n");
					printf("CNPr = %f, CNILc = %f \n", CNPrevised(countk1), CNILc);

					printf("\n");
					printf("r = %f \n", req);
					printf("Alpha' = %f, Alpha_cr = %f \n", DEG(alphap(countk1)), DEG(alphacr));

					printf("\n");
					printf("CNf = %f, CNfC = %f \n", CNf(countk1), CNfC(countk1));

					printf("fp = %f, fpp = %f \n", fp(countk1), fpp(countk1));
					printf("alpha_eff = %f, alpha_1 = %f \n", DEG(_aeff), DEG(alpha1));

					printf("\n");
					printf("CNv = %f, Cv = %f \n", CNv(countk1), Cv(countk1));

					printf("\n");
					printf("CNT = %f, CCf = %f \n", CNT(countk1), CCf(countk1));
					printf("***********************************************************\n");
				}	*/
			}
		}
		if (ik > 0 && isExit(ck))
			return true;
	}
	return false;
}

