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
		cd = cd0 + aoa * (cd1 + aoa*cd2);
	}
	else if (aoa < 3 * PI / 4)
	{
		cl = PI / 9 * (a1 + PI*PI / 81 * (a3 + PI*PI / 81 * a5));
		cl -= cl / (PI / 18 * 7)*(aoa - PI / 9);
		cl *= Sign(_aoa);

		cd = cd0 + PI / 9 * (cd1 + PI / 9 * cd2);
		cd -= cd / (PI / 2)*(aoa - PI / 2);
	}
	else if (aoa <= PI)
	{
		cl = PI / 9 * (a1 + PI*PI / 81 * (a3 + PI*PI / 81 * a5));
		cl -= cl / (PI / 18 * 7)*(1.5 * PI - aoa - PI / 9);
		cl *= Sign(_aoa);

		cd = cd0 + PI / 9 * (cd1 + PI / 9 * cd2);
		cd -= cd / (PI / 2)*(aoa - PI / 2);
	}


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

void Rotor::SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc)
{
	_setstates(vel, omg, dvel, domg, vc, wc, dvc, dwc, refcoord);

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
	WingV = M.WingV, RotorV = M.RotorV, fuselage = M.fuselage;
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
