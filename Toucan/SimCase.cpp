#include "stdafx.h"
#include "SimCase.h"

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
	InitTailRotor(trotor, mrotor.omega);
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
	R.bld.soltype = Rotation, R.adyna = PWake;
	R.nf = 72, R.ns = 40, R.ni = 10; 
	R.kwtip = 5, R.kwrot = 5; R.nk = R.nf*R.kwtip;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7;
	R.eflap = 0, R.khub = 0, R.del = 0, R.pitchroot = 0;
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
	R.chord.setvalue(0.558), R.sweep.setvalue(0);
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

void Rotor::InitWakePram(void)
{
	kwtip = 5, kwrot = 5, nk = nf*kwtip;
	//nf = naz, ns = nrs,
	outboard = 0.3, rtip = 1.0;
	rc0 = 0.004852173913043;
	haveGeo = false, haveStr = false;
}

void Rotor::InitVariables(void)
{	
	switch (type)
	{
	case Mrotor:
		mul = 0, lambdi_ag = lambdt_ag = lambdh_ag = 0.03;
		niter_a = niter_w = -1;
		power = torque = 0;
		sita[0] = RAD(12.1), sita[1] = RAD(1.18), sita[2] = RAD(3.12);
		beta[0] = precone, beta[1] = RAD(0.0), beta[2] = RAD(0.0);
		bld.GAf.pho = 1.0;
		bld.err_b = 1e-3, bld.dff = 15, bld.nperiod = 360 / 15, bld.nitermax = 30 * bld.nperiod;
		bld.sol.allocate(bld.nitermax);
		break;
	case Trotor:
		mul = 0, lambdi_ag = lambdt_ag = lambdh_ag = 0.01;
		niter_a = niter_w = -1;
		power = torque = 0;
		sita[0] = RAD(10.1), sita[1] = RAD(0.0), sita[2] = RAD(0.0);
		beta[0] = precone, beta[1] = RAD(0.0), beta[2] = RAD(0.0);
		bld.GAf.pho = 1.0;
		bld.err_b = 1e-3, bld.dff = 15, bld.nperiod = 360 / 15, bld.nitermax = 30 * bld.nperiod;
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
	nitermax = 20;
	sita_coll_max = RAD(30);
	sita_cycl_max = RAD(15);
	euler_max = RAD(20);
}

void Jobs::InitProject(void)
{
	nCase = 13;
	Mus.allocate(nCase), Pits.allocate(nCase);
	Mus.input("Mus.in");
	Pits.input("Pits.in");

	if (remove("max_c_del.output") != 0)
		printf("Remove max_c_del.output failed. \n");
	if (remove("sum_a1_del.output") != 0)
		printf("Remove sum_a1_del.output failed. \n");
	if (remove("sum_a2_del.output") != 0)
		printf("Remove sum_a2_del.output failed. \n");
}

void Jobs::SetSimCond(Copter &C, const int ic)
{
	double euler[3] = { 0,0,0 };
	C.vel_g[0] = C.RotorV[0].vtipa*Mus(ic);
	euler[1] = RAD(Pits(ic));
	C.refcoord.SetCoordinate(euler, "euler");
	//C.fuselage.refcoord.SetCoordinate(euler, "euler");

	// define wake
	for (int i = C.RotorV.size() - 1; i >= 0; --i)
		if (C.RotorV[i].adyna > 0)
			C.RotorV[i].InitWakePram();
}

void Jobs::PostProcess(Copter &C, const int ic, const int s, const int e)
{
	static Matrix2<double> uctrl(nCase, C.nfree);
	if (ic < e)
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl(ic, i) = DEG(C.controls[i]);
	if (ic + 1 == e)
		uctrl.output("uctrl.output", 6);
}