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
		CLa(i) = AirfoilParas(_ir++, i);
		alpha0(i) = AirfoilParas(_ir++, i);
		alpha1(i) = AirfoilParas(_ir++, i);
		dalpha1(i) = AirfoilParas(_ir++, i);
		s1(i) = AirfoilParas(_ir++, i);
		s2(i) = AirfoilParas(_ir++, i);
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
}

void LBDynamicStall::LBDSprepare(double c, double a, double _dt, int nf, int ns, int _nk, Matrix2<double> &aoa, Matrix2<double> &Ma, Matrix2<double> &q0)
{
	airfoil.SetAirfoil("NACA0012.txt");
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
	CNCM2.allocate(Nf,Ns), CNIM2.allocate(Nf, Ns), CNPM2.allocate(Nf, Ns), CNaIM2.allocate(Nf, Ns), CNqIM2.allocate(Nf, Ns);
	CNfM2.allocate(Nf, Ns), CNfCM2.allocate(Nf, Ns), CDM2.allocate(Nf, Ns);
	CNPrevisedM2.allocate(Nf, Ns), CvM2.allocate(Nf, Ns), CNvM2.allocate(Nf, Ns), CNTM2.allocate(Nf, Ns), CCfM2.allocate(Nf, Ns);
	AoALSSaveM2.allocate(Nf, Ns), AoATSSaveM2.allocate(Nf, Ns), AoAATSaveM2.allocate(Nf, Ns);

	alphaIn.allocate(_nk*Nf);
	x1.allocate(_nk*Nf), y1.allocate(_nk*Nf), D1.allocate(_nk*Nf), D2.allocate(_nk*Nf);
	CNC.allocate(_nk*Nf), CNI.allocate(_nk*Nf), CNP.allocate(_nk*Nf), CNaI.allocate(_nk*Nf), CNqI.allocate(_nk*Nf);
	CNf.allocate(_nk*Nf), CNfC.allocate(_nk*Nf);
	fpp.allocate(_nk*Nf), fp.allocate(_nk*Nf), Df.allocate(_nk*Nf), Dfp.allocate(_nk*Nf), alphaeff.allocate(_nk*Nf), alphaE.allocate(_nk*Nf);
	CNPrevised.allocate(_nk*Nf), CNfrevised.allocate(_nk*Nf), Dp.allocate(_nk*Nf);
	Cv.allocate(_nk*Nf), CNv.allocate(_nk*Nf);
	CNT.allocate(_nk*Nf), CCf.allocate(_nk*Nf), CD.allocate(_nk*Nf);
	alphak.allocate(_nk*Nf), qk.allocate(_nk*Nf);
	AoALSSave.allocate(_nk*Nf), AoATSSave.allocate(_nk*Nf), AoAATSave.allocate(_nk*Nf);

	flowstate = Attach;
	vortexstate = Novortex;
	tv = 0;
}

void LBDynamicStall::AttachFlow(int iz, int ir)
{
	double _dak, _dak2, _dqk, _dqk2, beta;
	double _temp;
	int iz1 = iz - 1;
	int iz2;
	if (iz1 < 0)
		iz1 = Nf - 1;
	iz2 = iz1 - 1;
	if (iz2 < 0)
		iz2 = Nf - 1;
	alphaLc = aoa0M2(iz, ir); // 0.5*(aoa0M2(iz, ir) + aoa0M2(iz1, ir));
	MaLc = 0.5*(MaM2(iz, ir) + MaM2(iz1, ir));
	ds = MaLc*dt*2/ Tl;
	beta2 = 1 - MaLc * MaLc;
	beta = sqrt(beta2);

	CNMa = airfoil.CLa.interplinear_fast(MaIn, MaLc);
	CNILc = airfoil.CNl.interplinear_fast(MaIn, MaLc);

	countk1 = countk + 1;
	alphaIn(countk1) = DEG(alphaLc);

	_dak = aoa0M2(iz, ir) - aoa0M2(iz1, ir);
	x1(countk1) = x1(countk) * exp(-b1*beta2*ds) + A1 * _dak * exp(-b1*beta2*ds*0.5);
	y1(countk1) = y1(countk) * exp(-b2*beta2*ds) + A2 * _dak * exp(-b2*beta2*ds*0.5);
	alphaE(countk1) = alphaLc - x1(countk1) - y1(countk1);
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

	// constant
	Tp = airfoil.Tp.interplinear_fast(MaIn, MaLc);
	Tf = airfoil.Tf.interplinear_fast(MaIn, MaLc);
	alpha0 = airfoil.alpha0.interplinear_fast(MaIn, MaLc);
	alpha1 = airfoil.alpha1.interplinear_fast(MaIn, MaLc);
	dalpha1 = airfoil.dalpha1.interplinear_fast(MaIn, MaLc);
	s1 = airfoil.s1.interplinear_fast(MaIn, MaLc);
	s2 = airfoil.s2.interplinear_fast(MaIn, MaLc);
	Tv = airfoil.Tv.interplinear_fast(MaIn, MaLc);
	Tvl = airfoil.Tvl.interplinear_fast(MaIn, MaLc);
	DFLc = airfoil.Df.interplinear_fast(MaIn, MaLc);
	Cd0Lc = airfoil.CD0.interplinear_fast(MaIn, MaLc);

	printf("AOA = %f alphaE = %f\n", DEG(alphaLc), DEG(alphaE(countk1)));
	//printf("D1 = %f\n", D1(countk1));
	printf("CNaI = %f, CNqI = %f \n", CNaI(countk1), CNqI(countk1));
	printf("CNC = %f, CNI = %f, CNP = %f \n", CNC(countk1), CNI(countk1), CNP(countk1));

	//printf("Ta = %f, Tq = %f \n", Ta, Tq);
	//printf("Tp = %f, Tf = %f, alpha0 = %f, alpha1 = %f \n", Tp, Tf, DEG(alpha0), DEG(alpha1));
	//printf("s1 = %f, s2 = %f, Tv = %f, Tvl = %f, ds = %f, dt = %f\n", s1, s2, Tv, Tvl, ds, dt);
}

void LBDynamicStall::SeparateFlow(void)
{
	double _aeff;
	double _Tf;
	static double AoAsp = 0;
	_Tf = Tf;

	Dp(countk1) = Dp(countk)*exp(-ds / Tp) + (CNP(countk1) - CNP(countk))*exp(-ds / Tp*0.5);
	CNPrevised(countk1) = CNP(countk1) - Dp(countk1); 
	alphaeff(countk1) = CNPrevised(countk1) / CNMa;
	_aeff = alphaeff(countk1);

	if(countk > 0 && (_aeff-alphaeff(countk))*(alphaeff(countk)-alphaeff(countk-1))>=0)
	//if (_aeff > alphaeff(countk))
	{
		if (CNPrevised(countk1) > CNILc)
		{
			flowstate = Separate;
			AoALSSave(countk1) = DEG(alphaLc);
			AoAsp = alphaLc;
		}
		else
		{
			flowstate = Attach;
			AoAATSave(countk1) = DEG(alphaLc);
		}
	}
	else if (flowstate != Attach)
	{
		_aeff += 0;// +aalpha1; //这个修正看不清
		flowstate = ReAttach;
		AoAATSave(countk1) = DEG(alphaLc);
	}

	
	switch (flowstate)
	{
	case LBDynamicStall::Attach:
		break;
	case LBDynamicStall::Separate:
		if (vortexstate == Conv)
			_Tf = 0.5*Tf;
		break;
	case LBDynamicStall::ReAttach:
		_Tf = 2 * Tf;
		if (vortexstate == Conv || vortexstate == Nearby)
			_Tf = 4 * Tf;
		if (CNPrevised(countk1) < CNILc)
			_Tf = 2 * Tf;
		break;
	default:
		break;
	}

	if (_aeff <= alpha1)
		fp(countk1) = 1 - 0.3*exp((_aeff - alpha1) / s1);
	else
		fp(countk1) = 0.04 + 0.66*exp((alpha1 - _aeff) / s2);

	Df(countk1) = Df(countk)*exp(-ds / _Tf) + (fp(countk1) - fp(countk))*exp(-ds / _Tf*0.5); //动态时相当于对K流分离点引入延迟
	fpp(countk1) = fp(countk1) - Df(countk1);
	fpptemp = 1 + sqrt(fpp(countk1));

	CNfC(countk1) = CNMa*0.25*fpptemp*fpptemp*alphaE(countk1);
	CNf(countk1) = CNfC(countk1) + CNI(countk1);
	
	printf("_________________________________________\n");
	//printf("Df = %f Dp = %f\n", Df(countk1), Dp(countk1));
	printf("CNPr = %f, CNfr = %f, CNILc = %f \n", CNPrevised(countk1), CNfrevised(countk1), CNILc);
	printf("CNf = %f, CNI = %f, CNfC = %f \n", CNf(countk1), CNI(countk1), CNfC(countk1));
	printf("fp = %f, fpp = %f \n", fp(countk1), fpp(countk1));
	printf("alpha_eff = %f, alpha_1 = %f \n", DEG(_aeff), DEG(alpha1));
	printf("Flow State: %d\n", flowstate);
	printf("Tf = %f \n", _Tf);
}

void LBDynamicStall::DynamicStall(void)
{
	double _Tv = Tv;
	double _Ts = 4*(1 - fpp(countk1)) / 0.19;
	double _tv = tv + ds;

	if (flowstate == Separate)
	{
		tv += ds;
		if (tv <= Tvl)
			vortexstate = Conv;
		else if (tv <= 2 * Tvl)
			vortexstate = Nearby;
		else
		{
			vortexstate = Novortex;
			//tv = 0;
		}	
	}
	else if(flowstate == ReAttach)
	{
		if (tv > 0)
		{
			tv += ds;
			if (tv <= Tvl)
				vortexstate = Conv;
			else if (tv <= 2 * Tvl)
				vortexstate = Nearby;
			else
			{
				vortexstate = Novortex;
				tv = 0;
			}
		}
	}
	else
	{
		vortexstate = Novortex;
		tv = 0;
	}

	switch (vortexstate)
	{
	case LBDynamicStall::Novortex:
		break;
	case LBDynamicStall::Conv:
		Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
		CNv(countk1) = CNv(countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
		break;
	case LBDynamicStall::Nearby:
		_Tv = 0.5*Tv;
		Cv(countk1) = Cv(countk);
		CNv(countk1) = CNv(countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
		break;
	default:
		break;
	}

	//while (_tv > _Ts)
	//{
	//	_tv -= _Ts;
	//	if (flowstate == Separate)
	//	{
	//		if (_tv <= _Ts)
	//			vortexstate = Conv;
	//		else if (_tv <= 2 * _Ts)
	//			vortexstate = Nearby;
	//		else
	//			vortexstate = Novortex;
	//	}
	//	else if (flowstate == ReAttach)
	//	{
	//		if (_tv > 0)
	//		{
	//			if (_tv <= _Ts)
	//				vortexstate = Conv;
	//			else if (_tv <= 2 * _Ts)
	//				vortexstate = Nearby;
	//			else
	//			{
	//				vortexstate = Novortex;
	//				_tv = 0;
	//			}
	//		}
	//	}
	//	else
	//		break;
	//	
	//	switch (vortexstate)
	//	{
	//	case LBDynamicStall::Novortex:
	//		break;
	//	case LBDynamicStall::Conv:
	//		Cv(countk1) = CNC(countk1)*(1 - 0.25*fpptemp*fpptemp);
	//		CNv(countk1) = CNv(countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
	//		break;
	//	case LBDynamicStall::Nearby:
	//		_Tv = 0.5*Tv;
	//		Cv(countk1) = Cv(countk);
	//		CNv(countk1) = CNv(countk)*exp(-ds / _Tv) + (Cv(countk1) - Cv(countk))*exp(-ds / _Tv*0.5);
	//		break;
	//	default:
	//		break;
	//	}
	//}
	
	printf("Vortex State: %d \n", vortexstate);
	printf("tv = %f, CNv = %f, Cv = %f \n", tv, CNv(countk1), Cv(countk1));
	printf("Tvf = %f, Tv = %f, Ts = %f\n", Tvl, _Tv, _Ts);
}

void LBDynamicStall::LBDScomplete(void)
{
	countk++;
	CNT(countk1) = CNf(countk1) + CNv(countk1);
	CCf(countk1) = yita*CNMa*sin(alphaE(countk1))*alphaE(countk1)*(fpptemp - 1);
	
	if (CNPrevised(countk1)>=CNILc)
		CCf(countk1) *= pow(fpp(countk1), DFLc*(CNPrevised(countk1) - CNILc));

	CD(countk1) = Cd0Lc + CNT(countk1)*sin(alphaLc) - CCf(countk1)*cos(alphaLc);
	printf("CNT = %f, CCf = %f \n", CNT(countk1), CCf(countk1));
	printf("***********************************************************\n");
}

bool LBDynamicStall::isExit(int ir)
{
	//if收敛
	for (int iz = Nf-1; iz >= 0; iz--, countk1--)
	{
		CNTM2(iz, ir) = CNT(countk1);
		CCfM2(iz, ir) = CCf(countk1);
		CDM2(iz, ir) = CD(countk1);
		AoALSSaveM2(iz, ir) = AoALSSave(countk1);
		AoATSSaveM2(iz, ir) = AoATSSave(countk1);
		AoAATSaveM2(iz, ir) = AoAATSave(countk1);
	}
	return true;
}

void LBDynamicStall::FuncTest(void)
{
	double c = 0.1;
	double a = 340;
	double k = 0.1;//0.1;// 
	double Ma = 0.3;//0.3;// 
	double V = Ma*a;
	double omega = 2 * V*k / c;
	double df = 5;
	double dt = RAD(df) / omega;
	int nf = 360/df;
	int ns = 1;
	int nk = 4;
	double aAmp = RAD(10);
	double a0 = RAD(10);
	double q0 = 0;
	Matrix2<double> _Ma(nf, ns), _aoa(nf, ns), _aoaDG(nf,ns), _q0(nf, ns);
	
	_Ma.setvalue(Ma);

	for (int i = 0; i < nf; i++)
	{
		_aoa(i, 0) = a0 + aAmp*sin(omega*dt*i);
		_q0(i, 0) = aAmp*omega*cos(omega*dt*i)*c / V;
		_aoaDG(i, 0) = DEG(_aoa(i,0));
	}

	LBDSprepare(c, a, dt, nf, ns, nk, _aoa, _Ma, _q0);

	for (int ik = 0; ik < nk; ik++)
	{
		for (int iz = 0; iz < nf; iz++)
		{
			if (!(iz == 0 & ik == 0))
			{
				AttachFlow(iz, 0);
				SeparateFlow();
				DynamicStall();
				LBDScomplete();
			}
			else
				countk++;
		}
	}

	isExit(0);
	_aoaDG.output("AOAIn.output", 10);
	CNTM2.output("CNT_LBDtest1.output", 10);
	CCfM2.output("CCf_LBDtest1.output", 10);
	CDM2.output("CD_LBDtest1.output", 10);
	/*CNI.output("CNI_LBDtest1.output", 10);
	CNC.output("CNC_LBDtest1.output", 10);
	CNP.output("CNP_LBDtest1.output", 10);
	*/
	/*alphaE.output("alphaE.output", 10);*/
	alphaeff.output("alphaeff.output", 10);
	AoALSSaveM2.output("AoALeadSp.output", 3);
	AoATSSaveM2.output("AoATrailSp.output", 3);
	AoAATSaveM2.output("AoAAtt.output", 3);
}