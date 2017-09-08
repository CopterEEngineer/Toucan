#include "stdafx.h"
#include "SimCase.h"
#include <string>

void Model_UL496::GetProb(void)
{
	amb.rho = 0.002378;
	amb.vsound = 1115.48;
	simtype = FreeTrim1;
	nfree = 6;
}

void Model_UL496::GetModel(void)
{
	myTYPE origin[3], euler[3];
	Wing wing, fin1, fin2;
	Rotor mrotor, trotor;
	Coordinate *BASE_COORD = NULL;

	// reference coordinate
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	refcoord.SetCoordinate(origin, euler, BASE_COORD);
	BASE_COORD = &refcoord;
	refcoord.SetBase(BASE_COORD);

	// fuselage
	fuselage.amb = amb;
	fuselage.dragA = 0.0315*PI*132.25;
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	fuselage.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	//fuselage.refcoord.SetCoordinate(origin, euler, fuselage.refcoord.base);
	//BASE_COORD = &fuselage.refcoord;
	//fuselage.refcoord.SetBase(BASE_COORD);

	// wing
	wing.type = Hwing;
	wing.amb = amb;
	origin[0] = -7.87, origin[1] = 0, origin[2] = 0;
	euler[0] = 0, euler[1] = 4 * PI / 180, euler[2] = 0;
	wing.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	wing.a0 = 5.73, wing.cd0 = 0.04, wing.cd1 = 0, wing.cd2 = 0;
	wing.span = 3.75, wing.chord = 1, wing.taper = 1;
	//fin1
	fin1.type = Vwing;
	fin1.amb = amb;
	fin1.span = 1.33, fin1.chord = 1.13, fin1.taper = 0.51;
	fin1.a0 = 5.3, fin1.cd0 = 0.0105, fin1.cd1 = 0, fin1.cd2 = 0.01325;
	origin[0] = -(7.87 + 0.164), origin[1] = 3.75 / 2, origin[2] = 0;
	euler[0] = PI / 2, euler[1] = 0, euler[2] = 5 * PI / 180;
	fin1.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	// fin2
	fin2.type = Vwing;
	fin2.amb = amb;
	fin2.span = 1.33, fin2.chord = 1.13, fin2.taper = 0.51;
	fin2.a0 = 5.3, fin2.cd0 = 0.0105, fin2.cd1 = 0, fin2.cd2 = 0.01325;
	origin[0] = -(7.87 - 0.164), origin[1] = -3.75 / 2, origin[2] = 0;
	euler[0] = PI / 2, euler[1] = 0, euler[2] = 5 * PI / 180;
	fin2.refcoord.SetCoordinate(origin, euler, BASE_COORD);

	// main rotor
	InitMainRotor(mrotor);
	mrotor.SetCoordBase(BASE_COORD);
	

	// tail rotor
	//InitTailRotor(trotor, mrotor.omega);
	InitTailRotor(trotor, 525*PI/30);
	trotor.SetCoordBase(BASE_COORD);

	// entire parameters
	mass = 915;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			inmatx[i][j] = 0;
		}
	}
	inmatx[0][0] = 590 * 0.73696;
	inmatx[1][1] = 3098 * 0.73696;
	inmatx[2][2] = 2655 * 0.73696;
	
	inmatx_M.allocate(3, 3);
	inmatx_M(0, 0) = inmatx[0][0];
	inmatx_M(1, 1) = inmatx[1][1];
	inmatx_M(2, 2) = inmatx[2][2];

	WingV.clear(), RotorV.clear();
	WingV.emplace_back(wing);
	WingV.emplace_back(fin1);
	WingV.emplace_back(fin2);
	RotorV.emplace_back(mrotor);
	RotorV.emplace_back(trotor);
}

void Model_UL496::InitMainRotor(Rotor &R)
{
	R.type = Mrotor;
	R.teeter = true, R.nb = 2;
	R.amb = amb;
	R.bld.soltype = Rotation, R.adyna = PWake; // Averaged; // 
	R.nf = 72, R.ns = 40, R.ni = 10; 
	R.kwtip = 1, R.kwrot = 1; R.nk = R.nf*R.kwtip;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7;
	R.eflap = 0, R.khub = 0, R.del = 0, R.pitchroot = RAD(-2);
	R.radius = 11.5, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = 3.0*PI / 180;
	R.omega = 525 * PI / 30;
	R.vtipa = R.omega*R.radius;
	R.outboard = 0.3;
	R.rc0 = 0.004852173913043;
	R.cltc.allocate(40, 12), R.cltc.input("vr7_cl_c81.txt");
	R.cdtc.allocate(52, 12), R.cdtc.input("vr7_cd_c81.txt");
	R.iflap = 24.93, R.m1 = 3.25;
	R.sigma = 0.0309, R.gama = 5.01*amb.rho / 0.002378;	
	R.chord.setvalue(0.558); // chord have unit
	R.sweep.setvalue(0);
	myTYPE temp_twist, temp_azimuth, temp_station;
	for (int j = R.ns - 1; j >= 0; --j) {
		temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
		//temp_twist = (-8 * PI / 180) * j / (ns - 1);
		R.twist(j) = (-8 * PI / 180) * j / (R.ns - 1);
		for (int i = R.nf - 1; i >= 0; --i) {
			//twist(i, j) = temp_twist;
			R.rastation(i, j) = temp_station;
			R.azstation(i, j) = i * 2 * PI / R.nf;
		}
	}
	myTYPE origin[3], euler[3];
	origin[0] = -0.0625, origin[1] = 0, origin[2] = -4.217;
	euler[0] = 0, euler[1] = 3.0 * PI / 180, euler[2] = 0;
	R.hubfxcoord.SetCoordinate(origin, euler, R.refcoord.base);
	R.refcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = 0, euler[1] = PI, euler[2] = 0;
	R.hubrtcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubrtcoord
	//定义扭转，摆振，挥舞的三个桨叶坐标
	origin[0] = R.eflap, origin[1] = 0, origin[2] = 0;
	euler[0] = 0, euler[1] = R.precone, euler[2] = 0;
	R.bladecoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = 0, origin[1] = -R.precone * (1 - R.eflap), origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	R.tppcoord.SetCoordinate(origin, euler, R.refcoord.base);

}

void Model_UL496::InitTailRotor(Rotor &R, double w)
{
	R.type = Trotor;
	R.amb = amb;
	R.teeter = false;
	R.bld.soltype = HubFixed, R.adyna = Averaged;
	R.nf = 1, R.ns = 1, R.ni = 1; 
	R.kwtip = 0, R.kwrot = 0, R.nk = 1;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7;
	R.eflap = 0, R.khub = 0, R.del = 0;
	R.pitchroot = 0, R.radius = 1.3, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = 0, R.omega = 9.4912*w, R.vtipa = R.omega*R.radius;
	R.outboard = 0, R.rc0 = 0;
	R.cltc.allocate(R.nf, R.ns), R.cdtc.allocate(R.nf, R.ns);
	R.iflap = 0.0027, R.m1 = 0.00312, R.sigma = 0.1636, R.gama = 2.2265;
	R.chord.setvalue(0.167);

	myTYPE origin[3], euler[3];
	origin[0] = -13.48, origin[1] = 0, origin[2] = -1;
	euler[0] = PI / 2, euler[1] = 0, euler[2] = 0;
	R.hubfxcoord.SetCoordinate(origin, euler, R.refcoord.base);
	R.refcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = 0, euler[1] = PI, euler[2] = 0;
	R.hubrtcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubrtcoord
	//定义扭转，摆振，挥舞的三个桨叶坐标
	origin[0] = R.eflap, origin[1] = 0, origin[2] = 0;
	euler[0] = 0, euler[1] = R.precone, euler[2] = 0;
	R.bladecoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = 0, origin[1] = -R.precone * (1 - R.eflap), origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	R.tppcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// define wake
	if (R.adyna > 0)
		printf("Undefined tail wake func. \n");
}

void Rotor::WakeModelPrams(const int _k)
{
	kwtip = _k, kwrot = 5, nk = nf*kwtip;
	//nf = naz, ns = nrs,
	outboard = 0.3, rtip = 0.98;
	rc0 = 0.004852173913043;
	haveGeo = false, haveStr = false, outputWake = false;
	// update allocate
	tipstr.deallocate();
	tipstr.allocate(nk, nf);
	tipgeometry.deallocate();
	tipgeometry.allocate(nk, nf, 3);
}

void Rotor::InitVariables(void)
{	
	switch (type)
	{
	case Mrotor:
		power = power_c = power_i = power_f = power_o = 0;
		torque = torque_c = torque_i = torque_f = torque_o = 0;
		mul = 0, lambdi_ag = lambdt_ag = lambdh_ag = 0.03;
		niter_a = niter_w = -1;
		power = torque = 0;
		sita[0] = RAD(16.1), sita[1] = RAD(1.18), sita[2] = RAD(3.12);
		beta[0] = precone, beta[1] = RAD(0.0), beta[2] = RAD(0.0);
		bld.GAf.pho = 1.0;
		bld.err_b = 1e-3, bld.dff = 15, bld.nperiod = 360 / 15, bld.nitermax = 30 * bld.nperiod;
		bld.sol.allocate(bld.nitermax);
		break;
	case Trotor:
		power = power_c = power_i = power_f = power_o = 0;
		torque = torque_c = torque_i = torque_f = torque_o = 0; 
		mul = 0, lambdi_ag = lambdt_ag = lambdh_ag = 0.01;
		niter_a = niter_w = -1;
		power = torque = 0;
		sita[0] = RAD(10.1), sita[1] = RAD(0.0), sita[2] = RAD(0.0);
		beta[0] = precone, beta[1] = RAD(0.0), beta[2] = RAD(0.0);
		bld.GAf.pho = 1.0;
		bld.err_b = 1e-3, bld.dff = 15, bld.nperiod = 360 / 15, bld.nitermax = 30 * bld.nperiod;
		bld.sol.allocate(bld.nitermax);
		break;
	default:
		break;
	}
	for (int i = 0; i < 3; ++i)
	{
		airforce[i] = airforce_cg[i] = airmoment[i] = airmoment_cg[i] = 0;
		vel[i] = omg[i] = dvel[i] = domg[i] = 0;
	}
}

void Fuselage::InitVariables(void)
{
	for (int i = 0; i < 3; ++i)
	{
		airforce[i] = airforce_cg[i] = airmoment[i] = airmoment_cg[i] = 0;
		vel[i] = omg[i] = dvel[i] = domg[i] = 0;
	}
}

void Wing::InitVariables(void)
{
	for (int i = 0; i < 3; ++i)
	{
		airforce[i] = airforce_cg[i] = airmoment[i] = airmoment_cg[i] = 0;
		vel[i] = omg[i] = dvel[i] = domg[i] = 0;
	}
}

void CopterSolver::InitCopterSolver(void)
{
	err_a = err_c = 1.0e-3;
	epsilon = 0.01;
	nitermax = 35;
	sita_coll_max = RAD(50);
	sita_cycl_max = RAD(30);
	euler_max = RAD(40);
	sum_a1_del.allocate(nitermax), sum_a2_del.allocate(nitermax), max_c_del.allocate(nitermax);
}

void Jobs::InitProject(void)
{
	nCase = 13;
	Mus.allocate(nCase), Vfs.allocate(nCase), Pits.allocate(nCase), Kwtips.allocate(nCase);
	Consini.allocate(nCase, 6);
	uctrl.allocate(nCase, 6), beta.allocate(nCase, 3);
	_power.allocate(nCase, 6), _torque.allocate(nCase, 6);
	flg = 0, jtype = RPMSwp;
	Mus.input("Mus.in");
	Vfs.input("Vfs.in");
	Pits.input("Pits.in");
	Kwtips.input("Kwtips.in");
	Consini.input("Controls.in");

	if (remove("max_c_del.output") != 0)
		printf("Remove max_c_del.output failed. \n");
	if (remove("sum_a1_del.output") != 0)
		printf("Remove sum_a1_del.output failed. \n");
	if (remove("sum_a2_del.output") != 0)
		printf("Remove sum_a2_del.output failed. \n");
}

void Jobs::InitProject(Jobs &J, const int ic)
{
	nCase = 1;
	Mus.allocate(1), Vfs.allocate(1), Pits.allocate(1), Kwtips.allocate(1);
	uctrl.allocate(1, 6), beta.allocate(1, 3);
	_power.allocate(1, 6), _torque.allocate(1, 6);
	
	Mus(0) = J.Mus(ic), Vfs(0) = J.Vfs(ic), Pits(0) = J.Pits(ic), Kwtips(0) = J.Kwtips(ic);
	flg = 0, jtype = J.jtype;
}

void Jobs::SetSimCond(Copter &C, const int ic)
{
	double euler[3] = { 0,0,0 };
	double _uctrl[6] = { 0,0,0,0,0,0 };

	C.vel_g[0] = C.RotorV[0].vtipa*Mus(ic);
	//C.vel_g[0] = Vfs(ic);
	//euler[1] = RAD(Pits(ic));
	//C.refcoord.SetCoordinate(euler, "euler");

	for (int i = 0; i < 6; i++)
		_uctrl[i] = RAD(Consini(ic, i));

	euler[0] = _uctrl[4], euler[1] = _uctrl[5];
	C.refcoord.SetCoordinate(euler, "euler");
	
	C.SetCtrl(_uctrl, C.RotorV[0]);
	C.SetCtrl(_uctrl + 3, C.RotorV[1]);
	
	// define wake
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{		
		if (jtype > SimTrim)
			C.RotorV[i].InitVariables();
		if (C.RotorV[i].adyna > 0)
			C.RotorV[i].WakeModelPrams(Kwtips(ic));
	}
}

void Jobs::SetSimCond(Copter &C)
{
	double euler[3] = { 0,0,0 };
	C.vel_g[0] = C.RotorV[0].vtipa*Mus(0);
	//C.vel_g[0] = Vfs(0);
	euler[1] = RAD(Pits(0));
	C.refcoord.SetCoordinate(euler, "euler");
	//C.fuselage.refcoord.SetCoordinate(euler, "euler");

	// define wake
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		if (jtype > SimTrim)
			C.RotorV[i].InitVariables();
		if (C.RotorV[i].adyna > 0)
			C.RotorV[i].WakeModelPrams(Kwtips(0));
	}
}

void Jobs::PostProcess(Copter &C)
{
	double _beta[3];
	for (int i = C.nfree - 1; i >= 0; --i)
		uctrl(0, i) = DEG(C.controls[i]);
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
	{
		_power(0, 0) += C.RotorV[i].power;
		_power(0, 1) += C.RotorV[i].power_i;
		_power(0, 2) += C.RotorV[i].power_o;
		_power(0, 3) += C.RotorV[i].power_f;
		_power(0, 4) += C.RotorV[i].power_c;
		_power(0, 5) += C.RotorV[i].power_iid;
		_torque(0, 0) += C.RotorV[i].torque;
		_torque(0, 1) += C.RotorV[i].torque_i;
		_torque(0, 2) += C.RotorV[i].torque_o;
		_torque(0, 3) += C.RotorV[i].torque_f;
		_torque(0, 4) += C.RotorV[i].torque_c;
		_torque(0, 5) += C.RotorV[i].torque_iid;
	}
	C.RotorV[0].GetBeta(_beta);
	for (int i = 0; i < 3; ++i)
		beta(0, i) = DEG(_beta[i]);
	printf("\nMain rotor power: %f \n", C.RotorV[0].power);
	// vel test
	printf("Vel g test: %f \n", C.vel_g[0]);
}

void Jobs::PostProcess(Copter &C, const int ic, const int s, const int e)
{
	//static Matrix2<double> uctrl(nCase, C.nfree), beta(nCase, 3);
	//static Matrix2<double> _power(nCase, 6), _torque(nCase, 6);
	double _beta[3];
	if (ic < e)
	{
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl(ic, i) = DEG(C.controls[i]);		
		for (int i = C.RotorV.size() - 1; i >= 0; --i)
		{
			_power(ic, 0) += C.RotorV[i].power;
			_power(ic, 1) += C.RotorV[i].power_i;
			_power(ic, 2) += C.RotorV[i].power_o;
			_power(ic, 3) += C.RotorV[i].power_f;
			_power(ic, 4) += C.RotorV[i].power_c;
			_power(ic, 5) += C.RotorV[i].power_iid;
			_torque(ic, 0) += C.RotorV[i].torque;
			_torque(ic, 1) += C.RotorV[i].torque_i;
			_torque(ic, 2) += C.RotorV[i].torque_o;
			_torque(ic, 3) += C.RotorV[i].torque_f;
			_torque(ic, 4) += C.RotorV[i].torque_c;
			_torque(ic, 5) += C.RotorV[i].torque_iid;
		}
		C.RotorV[0].GetBeta(_beta);
		for (int i = 0; i < 3; ++i)
			beta(ic, i) = DEG(_beta[i]);
		printf("\n Main rotor power: %f \n", C.RotorV[0].power);
	}
		

	if (ic + 1 == e)
	{
		uctrl.output("uctrl.output", 6);
		beta.output("beta.output", 4);
		_power.output("power.output", 10);
		_torque.output("torque.output", 10);
	}
	for (int i = C.RotorV.size() - 1; i>=0;--i) 
		C.RotorV[i].OutPutWake(ic);

	// vel test
	printf("\n Vel g test: %f \n", C.vel_g[0]);
}

void Jobs::PostProcess(Copter &C, const int ic, const int ip, const int s, const int e)
{
	static Matrix3<double> uctrl(nCase, C.nfree, nParams), beta(nCase, 3, nParams), err(nCase, 3, nParams);
	static Matrix3<double> _power(nCase, 6, nParams), _torque(nCase, 6, nParams); 
	double _beta[3];
	if (ic < e && ip < nParams)
	{
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl(ic, i, ip) = DEG(C.controls[i]);
		for (int i = C.RotorV.size() - 1; i >= 0; --i)
		{
			_power(ic, 0, ip) += C.RotorV[i].power;
			_power(ic, 1, ip) += C.RotorV[i].power_i;
			_power(ic, 2, ip) += C.RotorV[i].power_o;
			_power(ic, 3, ip) += C.RotorV[i].power_f;
			_power(ic, 4, ip) += C.RotorV[i].power_c;
			_power(ic, 5, ip) += C.RotorV[i].power_iid;
			_torque(ic, 0, ip) += C.RotorV[i].torque;
			_torque(ic, 1, ip) += C.RotorV[i].torque_i;
			_torque(ic, 2, ip) += C.RotorV[i].torque_o;
			_torque(ic, 3, ip) += C.RotorV[i].torque_f;
			_torque(ic, 4, ip) += C.RotorV[i].torque_c;
			_torque(ic, 5, ip) += C.RotorV[i].torque_iid;
		}
		C.RotorV[0].GetBeta(_beta);
		for (int i = 0; i < 3; ++i)
			beta(ic, i, ip) = DEG(_beta[i]);
		err(ic, 0, ip) = C.sum_a1_del;
		err(ic, 1, ip) = C.sum_a2_del;
		err(ic, 2, ip) = C.max_c_del;
	}
	if (ic + 1 == e && ip + 1 == nParams)
	{
		uctrl.output2("uctrl_rpm.output", 6);
		beta.output2("beta_rpm.output", 4);
		_power.output2("power_rpm.output", 10);
		_torque.output2("torque_rpm.output", 10);
		err.output2("err_rpm.output", 10);
	}
	// vel and params test 
	printf("\n");
	printf("Vel g test: %f \n", C.vel_g[0]);
	printf("Omg test: %f \n", C.RotorV[0].omega/param0(0));
}

void Jobs::PostProcessMP(Copter &C, const int ic, const int s, const int e)
{
	double _beta[3];
	if (flg++ < e - s)
	{
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl(ic, i) = DEG(C.controls[i]);
		for (int i = C.RotorV.size() - 1; i >= 0; --i)
		{
			_power(ic, 0) += C.RotorV[i].power;
			_power(ic, 1) += C.RotorV[i].power_i;
			_power(ic, 2) += C.RotorV[i].power_o;
			_power(ic, 3) += C.RotorV[i].power_f;
			_power(ic, 4) += C.RotorV[i].power_c;
			_power(ic, 5) += C.RotorV[i].power_iid;
			_torque(ic, 0) += C.RotorV[i].torque;
			_torque(ic, 1) += C.RotorV[i].torque_i;
			_torque(ic, 2) += C.RotorV[i].torque_o;
			_torque(ic, 3) += C.RotorV[i].torque_f;
			_torque(ic, 4) += C.RotorV[i].torque_c;
			_torque(ic, 5) += C.RotorV[i].torque_iid;
		}
		C.RotorV[0].GetBeta(_beta);
		for (int i = 0; i < 3; ++i)
			beta(ic, i) = DEG(_beta[i]);
		//printf("\n Main rotor power: %f \n", C.RotorV[0].power);
	}

	printf("Flag: %d\n", flg);
	if (flg == e - s)
	{
		uctrl.output("uctrl.output", 6);
		beta.output("beta.output", 4);
		_power.output("power.output", 10);
		_torque.output("torque.output", 10);
	}
	/*for (int i = C.RotorV.size() - 1; i >= 0; --i)
		C.RotorV[i].OutPutWake(ic);*/

	// vel test
	//printf("\n Vel g test: %f \n", C.vel_g[0]);
}


void Jobs::ParamSweep(const Copter &C)
{
	RPMs.allocate(2, nCase), param0.allocate(2);
	nParams = 21;
	RPMs.input("RPMs.in");
	param0(0) = C.RotorV[0].omega, param0(1) = C.RotorV[0].vtipa;
}

void Jobs::UpdateParam(Copter &C, const int ic, const int ip)
{
	double _rpm = RPMs(0, ic) + (RPMs(1, ic) - RPMs(0, ic)) / (nParams - 1)*ip;
	//_rpm = 1.0;
	C.RotorV[0].omega = param0(0)*_rpm;
	C.RotorV[0].vtipa = param0(1)*_rpm;
}

void Rotor::OutPutWake(const int ic)
{
	string _ic = std::to_string(ic);
	string _im = std::to_string(type);

	if (outputWake)
	{
		tipgeometry.output("tipgeo_" + _im + "_" + _ic + ".output", 10);
		lambdi.output("lambdi_" + _im + "_" + _ic + ".output", 10);
		lambdx.output("lambdx_" + _im + "_" + _ic + ".output", 10);
		lambdy.output("lambdy_" + _im + "_" + _ic + ".output", 10);
		lambdh.output("lambdh_" + _im + "_" + _ic + ".output", 6);
		bladedeform.output("bladedeform_" + _im + "_" + _ic + ".output", 10);
	}
}

void LevelFlight(void)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e;

	ul496.GetProb();
	ul496.GetModel();

	copter.InitRotorCraft(ul496);
	jobs.InitProject();

	s = 0, e = 1, jobs.jtype = RPMSwp; // SimTrim; //

	for (int i = s; i < e; ++i)
	{
		cout << i << endl;
		jobs.SetSimCond(copter, i);
		solver.CopterSimulation(copter);
		jobs.PostProcess(copter, i, s, e);
	}
}

void LevelFlightMP(void)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e, i;
	
	jobs.InitProject();
	i = s = 0, e = jobs.nCase;

	ul496.GetProb();
	ul496.GetModel();
	copter.InitRotorCraft(ul496);

#pragma omp parallel num_threads(7) shared(s, e, jobs) firstprivate(i, solver, copter)
	{
#pragma omp for
		for (i = s; i < e; i++)
		{
			jobs.SetSimCond(copter, i);
			solver.CopterSimulation(copter);
			jobs.PostProcessMP(copter, i, s, e);
		}
#pragma omp barrier
	}

	//std::vector<Jobs> subjobs(e);
//#pragma omp parallel num_threads(6) shared(s, e, subjobs) firstprivate(i, jobs, solver, copter, ul496)
//	{
//#pragma omp for
//		for (i = s; i < e; i++)
//		{
//			subjobs[i].InitProject(jobs, i);
//
//			subjobs[i].SetSimCond(copter);
//			solver.CopterSimulation(copter);
//
//			subjobs[i].PostProcess(copter);
//		}
//#pragma omp barrier
//
//	}

}