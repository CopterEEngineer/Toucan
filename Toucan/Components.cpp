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
	airforce[0] = -0.5*vel[0] * vel[0] * amb.rho * dragA;
	airforce[1] = airforce[2] = 0;
	airmoment[0] = airmoment[1] = airmoment[2] = 0;

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
	ar = span*span / stot;
	a03d = a0*ar / (ar + 2 * (ar + 4) / (ar + 2));

	aoa = Atan2(-vel[2], vel[0]);
	// ignore wake disturb, makes the freestream face to wing positive. 
	// This treatmean is easily to follow blade, but, it is not agreed with the refcoord coordinate,
	// which freestream facing to wing is negative.
	cl = a03d*aoa;
	cd = cd0 + aoa * (cd1 + aoa*cd2);

	vel2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
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
	mul = vel[0] / vtipa;
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
}

void Copter::SetStates(void)
{
	// will update class members
	vel_c[0] = vel_c[1] = vel_c[2] = 0;
	omg_c[0] = omg_c[1] = omg_c[2] = 0;
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
		break;
	default:
		break;
	}
	
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