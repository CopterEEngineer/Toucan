#include "stdafx.h"
#include "SimCase.h"
#include <string>
#include <direct.h>
#include <math.h>

void Model_UL496::GetProb(void)
{
	amb.rho = 0.002378;
	amb.vsound = 1115.48;
	amb.height = 0;
	simtype = FreeTrim1;
	nfree = 6;
}

void Model_UL496::GetProb(double h)
{
	simtype = FreeTrim1;// FreeTrim1;// Transient;
	nfree = 6;
	amb.SetAmbience(h);
}

void Model_UL496::GetProb(int nf, double h, SimType st)
{
	simtype = st;
	nfree = nf;
	amb.SetAmbience(h);
}

void Model_UL496::GetModel(void)
{
	myTYPE origin[3], euler[3];
	Wing wing, fin1, fin2;
	Rotor mrotor, trotor;
	Coordinate *BASE_COORD = NULL;
	si_unit = false;
	// reference coordinate
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	refcoord.SetCoordinate(origin, euler, BASE_COORD);
	BASE_COORD = &refcoord;
	refcoord.SetBase(BASE_COORD);

	// fuselage
	fuselage.fmdling = Parasite;
	fuselage.amb = amb;
	fuselage.si_unit = si_unit;
	fuselage.dragA = 0.0315*PI*132.25;
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	fuselage.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	fuselage.cxtc.allocate(1, 1), fuselage.cytc.allocate(1, 1), fuselage.cztc.allocate(1, 1);
	fuselage.cmtc.allocate(1, 1), fuselage.cntc.allocate(1, 1);
	//fuselage.refcoord.SetCoordinate(origin, euler, fuselage.refcoord.base);
	//BASE_COORD = &fuselage.refcoord;
	//fuselage.refcoord.SetBase(BASE_COORD);

	// wing
	wing.type = Hwing;
	wing.amb = amb;
	wing.si_unit = si_unit;
	origin[0] = -7.87, origin[1] = 0, origin[2] = 0;
	euler[0] = 0, euler[1] = 4 * PI / 180, euler[2] = 0;
	wing.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	wing.a0 = 5.73, wing.cd0 = 0.04, wing.cd1 = 0, wing.cd2 = 0;
	wing.span = 3.75, wing.chord = 1, wing.taper = 1;
	//fin1
	fin1.type = Vwing;
	fin1.amb = amb;
	fin1.si_unit = si_unit;
	fin1.span = 1.33, fin1.chord = 1.13, fin1.taper = 0.51;
	fin1.a0 = 5.3, fin1.cd0 = 0.0105, fin1.cd1 = 0, fin1.cd2 = 0.01325;
	origin[0] = -(7.87 + 0.164), origin[1] = 3.75 / 2, origin[2] = 0;
	euler[0] = PI / 2, euler[1] = 0, euler[2] = 5 * PI / 180;
	fin1.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	// fin2
	fin2.type = Vwing;
	fin2.amb = amb;
	fin2.si_unit = si_unit;
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
	mass = 915; // 1000 // 1100/ 1130
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
	double temp = 1.0;
	R.type = Mrotor;
	R.teeter = true, R.nb = 2;
	R.hingetype = Teeter;
	R.amb = amb;
	R.si_unit = false;
	R.airfoil = C81Table;
	R.bld.soltype = Rotation, R.adyna = PWake; //  Averaged; //PWake
	R.nf = 72, R.ns = 40, R.ni = 10; 
	R.kwtip = 1, R.kwrot = 1; R.nk = R.nf*R.kwtip;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7;
	R.eflap = 0, R.khub = 0, R.del = 0, R.pitchroot = RAD(-2);
	R.radius = 11.5*temp, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI*temp*temp;
	R.precone = RAD(3.0);
	R.omega = 525 * PI / 30 * 1.0;
	R.vtipa = R.omega*R.radius*temp;
	R.outboard = 0.3;
	R.rc0 = 0.004852173913043;
	R.cltc.allocate(40, 12), R.cltc.input("vr7_cl_c81.txt");
	R.cdtc.allocate(52, 12), R.cdtc.input("vr7_cd_c81.txt");
	R.iflap = 24.93*temp*temp, R.m1 = 3.25*temp;
	R.sigma = 0.0309/temp, R.gama = 5.01*amb.rho / 0.002378*temp*temp;
	R.chord.setvalue(0.558); // chord have unit
	R.sweep.setvalue(0);
	myTYPE temp_twist, temp_azimuth, temp_station;
	for (int j = R.ns - 1; j >= 0; --j) {
		temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
		//temp_twist = (-8 * PI / 180) * j / (ns - 1);
		R.twist(j) = RAD(-8.0) * j / (R.ns - 1);
		for (int i = R.nf - 1; i >= 0; --i) {
			//twist(i, j) = temp_twist;
			R.rastation(i, j) = temp_station;
			R.azstation(i, j) = i * 2 * PI / R.nf;
		}
	}
	myTYPE origin[3], euler[3];
	origin[0] = -0.0625, origin[1] = 0, origin[2] = -4.217;
	euler[0] = 0, euler[1] = RAD(3.0), euler[2] = 0;
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
	R.FT = 1.0; 
	R.type = Trotor;
	R.amb = amb;
	R.si_unit = false;
	R.airfoil = C81Table;
	R.teeter = false, R.nb = 4;
	R.hingetype = Hinged;
	R.bld.soltype = HubFixed, R.adyna = Averaged; // Simple; // Averaged;

	if (R.adyna > -1)
		R.nf = 72, R.ns = 40, R.ni = 40;
	else
		R.nf = 1, R.ns = 1, R.ni = 1;
	
	R.kwtip = 0, R.kwrot = 0, R.nk = 1;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7;
	R.eflap = 0, R.khub = 0, R.del = 0, R.pitchroot = 0;
	R.radius = 1.3, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = 0, R.omega = 9.4912*w, R.vtipa = R.omega*R.radius;
	R.outboard = 0, R.rc0 = 0;

	R.iflap = 0.0027, R.m1 = 0.00312, R.sigma = 0.1636, R.gama = 2.2265;
	R.chord.setvalue(0.167);

	if (R.adyna > -1)
	{
		R.cltc.allocate(40, 12), R.cltc.input("naca0012_cl_c81.txt");
		R.cdtc.allocate(66, 12), R.cdtc.input("naca0012_cd_c81.txt");
		
		myTYPE temp_twist, temp_azimuth, temp_station;
		for (int j = R.ns - 1; j >= 0; --j) {
			temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
																	 //temp_twist = (-8 * PI / 180) * j / (ns - 1);
			R.twist(j) = (0.0 * PI / 180) * j / (R.ns - 1);
			for (int i = R.nf - 1; i >= 0; --i) {
				//twist(i, j) = temp_twist;
				R.rastation(i, j) = temp_station;
				R.azstation(i, j) = i * 2 * PI / R.nf;
			}
		}
	}
	
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
	if (R.adyna > 100)
		printf("Undefined tail wake func. \n");
}

void Model_BO105::GetProb(int nf, double h, SimType st)
{
	simtype = st;
	nfree = nf;
	amb.rho = 1.225;
	amb.vsound = 340.0;
}

void Model_BO105::GetModel(void)
{
	// SI unit
	myTYPE origin[3], euler[3];
	Wing wing, fin;
	Rotor mrotor, trotor;
	Coordinate *BASE_COORD = NULL;

	si_unit = true;

	// entire parameters
	mass = 2200; 
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			inmatx[i][j] = 0;
		}
	}
	inmatx[0][0] = 1433;
	inmatx[1][1] = 4973;
	inmatx[2][2] = 4099;
	inmatx[0][2] = inmatx[2][0] = 660;

	inmatx_M.allocate(3, 3);
	inmatx_M(0, 0) = inmatx[0][0];
	inmatx_M(1, 1) = inmatx[1][1];
	inmatx_M(2, 2) = inmatx[2][2];
	inmatx_M(0, 2) = inmatx_M(2, 0) = inmatx[0][2];

	// reference coordinate, mass center
	origin[1] = origin[2] = 0;
	origin[0] = 0; 
	euler[0] = euler[1] = euler[2] = 0;
	refcoord.SetCoordinate(origin, euler, BASE_COORD);
	BASE_COORD = &refcoord;
	refcoord.SetBase(BASE_COORD);

	// fuselage
	InitFuselage();
	fuselage.dragA=0.0315*PI*4.91*4.91;
	fuselage.amb = amb;
	fuselage.si_unit = si_unit;
	origin[0] = 0;
	origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	fuselage.refcoord.SetCoordinate(origin, euler, BASE_COORD);

	// wing
	wing.type = Hwing;
	wing.amb = amb;
	wing.si_unit = si_unit;
	InitWing(wing);
	origin[0] = -4.56 - 0.0163, origin[1] = 0, origin[2] = -0.5;
	euler[0] = 0, euler[1] = 0, euler[2] = 0;
	wing.refcoord.SetCoordinate(origin, euler, BASE_COORD);
	
	//fin
	fin.type = Vwing;
	fin.amb = amb;
	fin.si_unit = si_unit;
	InitFin(fin);
	origin[0] = -5.416 - 0.0163, origin[1] = 0, origin[2] = -1.0;
	euler[0] = PI / 2, euler[1] = 0, euler[2] = 0;
	fin.refcoord.SetCoordinate(origin, euler, BASE_COORD);

	// main rotor
	InitMainRotor(mrotor);
	mrotor.SetCoordBase(BASE_COORD);

	// tail rotor
	InitTailRotor(trotor, 44.4);
	trotor.SetCoordBase(BASE_COORD);

	WingV.clear(), RotorV.clear();
	WingV.emplace_back(wing);
	WingV.emplace_back(fin);
	RotorV.emplace_back(mrotor);
	RotorV.emplace_back(trotor);
}

void Model_BO105::InitFuselage(void)
{
	fuselage.fmdling = Fitting; // Fitting;
	fuselage.Vtest = 30.48; //m/s
	fuselage.xf0 = -580.6, fuselage.xf1 = -454.0, fuselage.xf2 = 6.2, fuselage.xf3 = 4648.9;
	fuselage.yf0 = -6.9, fuselage.yf1 = -2399.0, fuselage.yf2 = -1.7, fuselage.yf3 = 12.7;
	fuselage.zf0 = -51.1, fuselage.zf1 = -1202.0, fuselage.zf2 = 1515.7, fuselage.zf3 = -604.2;
	fuselage.mf0 = -1191.8, fuselage.mf1 = 12752.0, fuselage.mf2 = 8201.3, fuselage.mf3 = -5796.7;
	fuselage.nf0 = fuselage.nf2 = fuselage.nf3 = 0;
	fuselage.nf1 = -10028.0;
	fuselage.cxtc.allocate(9, 2), fuselage.cytc.allocate(11, 2), fuselage.cztc.allocate(11, 2);
	fuselage.cmtc.allocate(10, 2), fuselage.cntc.allocate(9, 2);
	fuselage.cxtc.input("fuselage_cx.txt"), fuselage.cytc.input("fuselage_cy.txt"), fuselage.cztc.input("fuselage_cz.txt");
	fuselage.cmtc.input("fuselage_cm.txt"), fuselage.cntc.input("fuselage_cn.txt");
	fuselage.Sp = 14, fuselage.Lf = 20, fuselage.Ss = 80;
}

void Model_BO105::InitWing(Wing &W)
{
	W.wmdling = WFitting;
	W.a1 = 3.262, W.a3 = W.a5 = 0;
	W.alpha0 = 0.0698;
	W.cd0 = 0, W.cd1 = 0, W.cd2 = 0;
	W.span = 2, W.chord = 0.4015, W.taper = 1;
}

void Model_BO105::InitFin(Wing &F)
{
	F.wmdling = WFitting;
	F.a1 = 2.704, F.a3 = F.a5 = 0;
	F.alpha0 = -0.08116;
	F.cd0 = F.cd1 = F.cd2 = 0;
	F.span = 1, F.chord = 0.805, F.taper = 1;
}

void Model_BO105::InitMainRotor(Rotor &R)
{
	R.type = Mrotor;
	R.teeter = false;
	R.hingetype = Hingeless;
	R.nb = 4;
	R.amb = amb;
	R.si_unit = true;
	R.airfoil = Padfield;
	R.bld.soltype = HubFixed, R.adyna = LinearInflow;
	R.nf = 72, R.ns = 40, R.ni = 10;
	R.kwtip = 1, R.kwrot = 1; R.nk = R.nf*R.kwtip;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 6.113, R.del0 = 0.0074, R.del2 = 38.66, R.alpha0 = 0.0;
	R.eflap = 0;
	R.khub = 113330, R.del = 0, R.pitchroot = RAD(0);
	R.radius = 4.91, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = 0;
	R.omega = 44.4;
	R.vtipa = R.omega*R.radius;
	R.outboard = 0.3;
	R.rc0 = 0.004852173913043;
	//R.cltc.allocate(1, 1), R.cdtc.allocate(1, 1);
	R.cltc.allocate(40, 12), R.cltc.input("naca0012_cl_c81.txt");
	R.cdtc.allocate(66, 12), R.cdtc.input("naca0012_cd_c81.txt");
	R.iflap = 231.7, R.m1 = 1.5*R.iflap/R.radius;
	R.sigma = 0.07, R.gama = 5.087;
	R.chord.setvalue(0.27);
	R.sweep.setvalue(0);
	myTYPE temp_twist, temp_azimuth, temp_station;
	for (int j = R.ns - 1; j >= 0; --j) {
		temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
																 //temp_twist = (-8 * PI / 180) * j / (ns - 1);
		R.twist(j) = -0.14 * j / (R.ns - 1);
		for (int i = R.nf - 1; i >= 0; --i) {
			//twist(i, j) = temp_twist;
			R.rastation(i, j) = temp_station;
			R.azstation(i, j) = i * 2 * PI / R.nf;
		}
	}

	myTYPE origin[3], euler[3];
	origin[0] = -0.0163, origin[1] = 0, origin[2] = -1.48;
	euler[0] = 0, euler[1] = PI + RAD(3), euler[2] = 0;
	R.hubfxcoord.SetCoordinate(origin, euler, R.refcoord.base);
	R.refcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = 0, euler[1] = 0, euler[2] = 0;
	R.hubrtcoord.SetCoordinate(origin, euler, R.refcoord.base);

	origin[0] = R.eflap, origin[1] = 0, origin[2] = 0;
	euler[0] = R.precone, euler[1] = 0, euler[2] = -PI/2;
	R.bladecoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	R.windcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = 0, origin[1] = 0, origin[2] = R.precone * (1 - R.eflap);
	euler[0] = euler[1] = euler[2] = 0;
	R.tppcoord.SetCoordinate(origin, euler, R.refcoord.base);
}

void Model_BO105::InitTailRotor(Rotor &R, double w)
{
	R.FT = 0.787;
	R.type = Trotor;
	R.amb = amb;
	R.si_unit = si_unit;
	R.airfoil = Padfield;
	R.hingetype = Teeter; 
	R.teeter = true;
	R.nb = 2;
	R.bld.soltype = HubFixed, R.adyna = Simple; // Simple; // Averaged;

	if (R.adyna > -1)
		R.nf = 72, R.ns = 40, R.ni = 40;
	else
		R.nf = 1, R.ns = 1, R.ni = 1;

	R.kwtip = 0, R.kwrot = 0, R.nk = 1;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7, R.del0 = 0.008, R.del2 = 9.5, R.alpha0 = 0;
	R.eflap = 0, R.khub = 0, R.del = 45, R.pitchroot = 0;
	R.radius = 0.95, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = RAD(0.5), R.omega = 5.25*w, R.vtipa = R.omega*R.radius;
	R.outboard = 0, R.rc0 = 0;
	
	//R.iflap = -99, R.m1 = -99, R.gama = -99;
	R.iflap = 2.9, R.m1 = 1.5*R.iflap / R.radius, R.gama = 0.1764; // 猜测值
	R.sigma = 0.12;
	R.chord.setvalue(0.09);

	if (R.adyna > -1)
	{
		//R.cltc.allocate(1, 1);
		//R.cdtc.allocate(1, 1);
		R.cltc.allocate(40, 12), R.cltc.input("naca0012_cl_c81.txt");
		R.cdtc.allocate(66, 12), R.cdtc.input("naca0012_cd_c81.txt");

		myTYPE temp_twist, temp_azimuth, temp_station;
		for (int j = R.ns - 1; j >= 0; --j) {
			temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
																	 //temp_twist = (-8 * PI / 180) * j / (ns - 1);
			R.twist(j) = (0.0 * PI / 180) * j / (R.ns - 1);
			for (int i = R.nf - 1; i >= 0; --i) {
				//twist(i, j) = temp_twist;
				R.rastation(i, j) = temp_station;
				R.azstation(i, j) = i * 2 * PI / R.nf;
			}
		}
	}

	myTYPE origin[3], euler[3];
	origin[0] = -6 - 0.0163, origin[1] = 0, origin[2] = -1.72;
	euler[0] = -PI / 2, euler[1] = PI, euler[2] = 0;
	R.hubfxcoord.SetCoordinate(origin, euler, R.refcoord.base);
	R.refcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = 0, euler[1] = 0, euler[2] = 0;
	R.hubrtcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubrtcoord
	//定义扭转，摆振，挥舞的三个桨叶坐标
	origin[0] = R.eflap, origin[1] = 0, origin[2] = 0;
	euler[0] = R.precone, euler[1] = 0, euler[2] = -PI/2;
	R.bladecoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	R.windcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = 0, origin[1] = 0, origin[2] = R.precone * (1 - R.eflap);
	euler[0] = euler[1] = euler[2] = 0;
	R.tppcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// define wake
	if (R.adyna > 100)
		printf("Undefined tail wake func. \n");
}

void Model_Puma330::GetProb(int nf, double h, SimType st, bool unit)
{
	simtype = st;
	nfree = nf;
	amb.SetAmbience(h/0.3048);
	if (unit)
	{
		amb.rho = amb.rho / 0.002378*1.225;
		amb.vsound = amb.vsound / 1115 * 340;
	}
}

void Model_Puma330::GetModel(void)
{
	// SI unit
	myTYPE origin[3], euler[3];
	Wing wing, fin;
	Rotor mrotor, trotor;
	Coordinate *BASE_COORD = NULL;

	si_unit = true;

	// entire parameters
	mass = 5805;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			inmatx[i][j] = 0;
		}
	}
	inmatx[0][0] = 9638;
	inmatx[1][1] = 33240;
	inmatx[2][2] = 25889;
	inmatx[0][2] = inmatx[2][0] = 2226;

	inmatx_M.allocate(3, 3);
	inmatx_M(0, 0) = inmatx[0][0];
	inmatx_M(1, 1) = inmatx[1][1];
	inmatx_M(2, 2) = inmatx[2][2];
	inmatx_M(0, 2) = inmatx_M(2, 0) = inmatx[0][2];

	// reference coordinate
	origin[1] = origin[2] = 0;
	origin[0] = 0.005;
	euler[0] = euler[1] = euler[2] = 0;
	refcoord.SetCoordinate(origin, euler, BASE_COORD);
	BASE_COORD = &refcoord;
	refcoord.SetBase(BASE_COORD);

	// fuselage
	InitFuselage();
	fuselage.dragA = 0.0315*PI*7.5*7.5;
	fuselage.amb = amb;
	fuselage.si_unit = si_unit;
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	fuselage.refcoord.SetCoordinate(origin, euler, BASE_COORD);

	// wing
	wing.type = Hwing;
	wing.amb = amb;
	wing.si_unit = si_unit;
	InitWing(wing);
	origin[0] = -9, origin[1] = -1, origin[2] = -1.087;
	euler[0] = 0, euler[1] = 0, euler[2] = 0;
	wing.refcoord.SetCoordinate(origin, euler, BASE_COORD);

	//fin
	fin.type = Vwing;
	fin.amb = amb;
	fin.si_unit = si_unit;
	InitFin(fin);
	origin[0] = -9, origin[1] = 0, origin[2] = -1.087;
	euler[0] = PI / 2, euler[1] = 0, euler[2] = 0;
	fin.refcoord.SetCoordinate(origin, euler, BASE_COORD);

	// main rotor
	InitMainRotor(mrotor);
	mrotor.SetCoordBase(BASE_COORD);

	// tail rotor
	InitTailRotor(trotor, 27);
	trotor.SetCoordBase(BASE_COORD);

	WingV.clear(), RotorV.clear();
	WingV.emplace_back(wing);
	WingV.emplace_back(fin);
	RotorV.emplace_back(mrotor);
	RotorV.emplace_back(trotor);
}

void Model_Puma330::InitFuselage(void)
{
	fuselage.fmdling = Fitting; // Fitting;
	fuselage.Vtest = 30.48; //m/s
	fuselage.xf0 = -822.9, fuselage.xf1 = 44.5, fuselage.xf2 = 911.9, fuselage.xf3 = 1663.6;
	fuselage.yf0 = -11672.0, fuselage.yf1 = 0, fuselage.yf2 = 0, fuselage.yf3 = 0;
	fuselage.zf0 = -458.2, fuselage.zf1 = -5693.7, fuselage.zf2 = 2077.3, fuselage.zf3 = -3958.9;
	fuselage.mf0 = -1065.7, fuselage.mf1 = 8745.0, fuselage.mf2 = 12473.5, fuselage.mf3 = -10033.0;
	fuselage.nf0 = fuselage.nf2 = 0;
	fuselage.nf1 = -24269.2, fuselage.nf3 = 97619.0;
}

void Model_Puma330::InitWing(Wing &W)
{
	W.wmdling = WFitting;
	W.a1 = 3.7, W.a3 = -3.92*3.7, W.a5 = 0;
	W.alpha0 = -0.0262;
	W.cd0 = 0, W.cd1 = 0, W.cd2 = 0;
	W.span = 2.11, W.chord = 0.635, W.taper = 1;
}

void Model_Puma330::InitFin(Wing &F)
{
	F.wmdling = WFitting;
	F.span = 1.395, F.chord = 1, F.taper = 1;
	F.cd0 = 0, F.cd1 = 0, F.cd2 = 0;
	F.alpha0 = 0.0175;
	F.a1 = 0, F.a3 = 3.5*11.143, F.a5 = -3.5*85.714;
}

void Model_Puma330::InitMainRotor(Rotor &R)
{
	R.type = Mrotor;
	R.teeter = false;
	R.hingetype = Hinged;
	R.nb = 4;
	R.amb = amb;
	R.si_unit = true;
	R.airfoil = C81Table;
	R.bld.soltype = Rotation, R.adyna = Averaged;
	R.nf = 72, R.ns = 40, R.ni = 10;
	R.kwtip = 1, R.kwrot = 1; R.nk = R.nf*R.kwtip;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.73, R.del0 = 0.0074, R.del2 = 9.5, R.alpha0 = -0.0262;
	R.eflap = 0.038; // 
	R.khub = 48149, R.del = 0, R.pitchroot = RAD(0);
	R.radius = 7.5, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = 0;
	R.omega = 27;
	R.vtipa = R.omega*R.radius;
	R.outboard = 0.3;
	R.rc0 = 0.004852173913043;

	R.cltc.allocate(40, 12), R.cltc.input("naca0012_cl_c81.txt");
	R.cdtc.allocate(66, 12), R.cdtc.input("naca0012_cd_c81.txt");
	R.iflap = 1280, R.m1 = 1.5*R.iflap / R.radius;
	R.sigma = 0.09169, R.gama = 9.374;
	R.chord.setvalue(0.5401);
	R.sweep.setvalue(0);

	myTYPE temp_twist, temp_azimuth, temp_station;
	for (int j = R.ns - 1; j >= 0; --j) {
		temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
																 //temp_twist = (-8 * PI / 180) * j / (ns - 1);
		R.twist(j) = -0.14 * j / (R.ns - 1);
		for (int i = R.nf - 1; i >= 0; --i) {
			//twist(i, j) = temp_twist;
			R.rastation(i, j) = temp_station;
			R.azstation(i, j) = i * 2 * PI / R.nf;
		}
	}

	myTYPE origin[3], euler[3];
	origin[0] = 0, origin[1] = 0, origin[2] = -2.157;
	euler[0] = 0, euler[1] = 0.0873, euler[2] = 0;
	R.hubfxcoord.SetCoordinate(origin, euler, R.refcoord.base);
	R.refcoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = 0, euler[1] = PI, euler[2] = 0;
	R.hubrtcoord.SetCoordinate(origin, euler, R.refcoord.base);

	origin[0] = R.eflap, origin[1] = 0, origin[2] = 0;
	euler[0] = 0, euler[1] = R.precone, euler[2] = 0;
	R.bladecoord.SetCoordinate(origin, euler, R.refcoord.base);

	// referred to hubfxcoord
	origin[0] = 0, origin[1] = -R.precone * (1 - R.eflap), origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	R.tppcoord.SetCoordinate(origin, euler, R.refcoord.base);
}

void Model_Puma330::InitTailRotor(Rotor &R, double w)
{
	R.type = Trotor;
	R.amb = amb;
	R.si_unit = si_unit;
	R.airfoil = C81Table;
	R.teeter = false, R.nb = 5;
	R.hingetype = Hinged;
	R.bld.soltype = Rotation, R.adyna = Averaged; // Simple; // Averaged;

	if (R.adyna > -1)
		R.nf = 72, R.ns = 40, R.ni = 40;
	else
		R.nf = 1, R.ns = 1, R.ni = 1;

	R.kwtip = 0, R.kwrot = 0, R.nk = 1;
	R.haveGeo = false, R.haveStr = false, R.outputWake = false;
	R.chord.allocate(R.ns), R.sweep.allocate(R.ns), R.twist.allocate(R.ns);
	R.azstation.allocate(R.nf, R.ns), R.rastation.allocate(R.nf, R.ns);

	R.a0 = 5.7, R.del0 = 0.008, R.del2 = 9.5, R.alpha0 = 0;
	R.eflap = 0, R.khub = 0, R.del = 0, R.pitchroot = 0;
	R.radius = 1.56, R.bt = 0.98, R.rroot = 0.15, R.disk_A = R.radius*R.radius*PI;
	R.precone = 0, R.omega = 4.82*w, R.vtipa = R.omega*R.radius;
	R.outboard = 0, R.rc0 = 0;

	//R.iflap = -99, R.m1 = -99, R.gama = -99;
	R.iflap = 4.96, R.m1 = 1.5*R.iflap / R.radius, R.gama = 1.952; // 猜测值
	R.sigma = 0.19;
	R.chord.setvalue(0.2328);

	if (R.adyna > -1)
	{
		//R.cltc.allocate(1, 1);
		//R.cdtc.allocate(1, 1);
		R.cltc.allocate(40, 12), R.cltc.input("naca0012_cl_c81.txt");
		R.cdtc.allocate(66, 12), R.cdtc.input("naca0012_cd_c81.txt");

		myTYPE temp_twist, temp_azimuth, temp_station;
		for (int j = R.ns - 1; j >= 0; --j) {
			temp_station = R.rroot + j*(1.0 - R.rroot) / (R.ns - 1); // uniform seperated accepted so far.
																	 //temp_twist = (-8 * PI / 180) * j / (ns - 1);
			R.twist(j) = (0.0 * PI / 180) * j / (R.ns - 1);
			for (int i = R.nf - 1; i >= 0; --i) {
				//twist(i, j) = temp_twist;
				R.rastation(i, j) = temp_station;
				R.azstation(i, j) = i * 2 * PI / R.nf;
			}
		}
	}

	myTYPE origin[3], euler[3];
	origin[0] = -9, origin[1] = 0, origin[2] = -1.587;
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
	if (R.adyna > 100)
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
		mul = 0, betawind = 0, lambdi_ag = lambdt_ag = lambdh_ag = 0.03;
		niter_a = niter_w = -1;
		power = torque = 0;
		//sita[0] = RAD(16.1), sita[1] = RAD(1.18), sita[2] = RAD(3.12);
		beta[0] = precone, beta[1] = RAD(0.0), beta[2] = RAD(0.0);
		bld.GAf.pho = 1.0;
		bld.err_b = 1e-3, bld.dff = 15, bld.nperiod = 360 / 15, bld.nitermax = 30 * bld.nperiod;
		bld.sol.allocate(bld.nitermax);
		bld.azmuth.allocate(bld.nitermax);
		break;
	case Trotor:
		power = power_c = power_i = power_f = power_o = 0;
		torque = torque_c = torque_i = torque_f = torque_o = 0; 
		mul = 0, betawind = 0, lambdi_ag = lambdt_ag = lambdh_ag = 0.01;
		niter_a = niter_w = -1;
		power = torque = 0;
		//sita[0] = RAD(10.1);
		sita[1] = RAD(0.0), sita[2] = RAD(0.0);
		beta[0] = precone, beta[1] = RAD(0.0), beta[2] = RAD(0.0);
		bld.GAf.pho = 1.0;
		bld.err_b = 1e-3, bld.dff = 15, bld.nperiod = 360 / 15, bld.nitermax = 30 * bld.nperiod;
		bld.sol.allocate(bld.nitermax);
		bld.azmuth.allocate(bld.nitermax);
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

void Jobs::InitProject(JobsType jt)
{
	nCase = 8;
	Mus.allocate(nCase), Vfs.allocate(nCase), Pits.allocate(nCase), Kwtips.allocate(nCase);
	Consini.allocate(nCase, 6);
	uctrl.allocate(nCase, 6), beta.allocate(nCase, 3), err.allocate(nCase, 3);
	_power.allocate(nCase, 6), _torque.allocate(nCase, 6);
	uctrl_tr.allocate(nCase, 6), beta_tr.allocate(nCase, 3);
	_power_tr.allocate(nCase, 6), _torque_tr.allocate(nCase, 6);
	flg = 0, jtype = jt; // RadiusSwp; // RPMSwp; //ChordSwp; //SimTrim
	niter.allocate(nCase, 1);
	path.clear();
	//jtype_M.allocate(2), jtype_M(0) = RadiusSwp, jtype_M(1) = RPMSwp;

	//Mus.input("Mus.in");
	//Vfs.input("Vfs.in");
	Vfs.input("Vfs-bo105.in");
	//Pits.input("Pits.in");
	//Kwtips.input("Kwtips.in");
	Consini.input("ControlsInit-bo105.in");

	if (remove("max_c_del.output") != 0)
		printf("Remove max_c_del.output failed. \n");
	if (remove("sum_a1_del.output") != 0)
		printf("Remove sum_a1_del.output failed. \n");
	if (remove("sum_a2_del.output") != 0)
		printf("Remove sum_a2_del.output failed. \n");
}

void Jobs::InitProject(JobsType jt, string fn)
{
	nCase = 8;
	Mus.allocate(nCase), Vfs.allocate(nCase), Pits.allocate(nCase), Kwtips.allocate(nCase);
	Consini.allocate(nCase, 6);
	uctrl.allocate(nCase, 6), beta.allocate(nCase, 3), err.allocate(nCase, 3);
	_power.allocate(nCase, 6), _torque.allocate(nCase, 6);
	uctrl_tr.allocate(nCase, 6), beta_tr.allocate(nCase, 3);
	_power_tr.allocate(nCase, 6), _torque_tr.allocate(nCase, 6);

	dXdvel.allocate(nCase, 3), dYdvel.allocate(nCase, 3), dZdvel.allocate(nCase, 3);
	dXdomg.allocate(nCase, 3), dYdomg.allocate(nCase, 3), dZdomg.allocate(nCase, 3);
	dXdeul.allocate(nCase, 2), dYdeul.allocate(nCase, 2), dZdeul.allocate(nCase, 2);
	dXdctrl.allocate(nCase, 4), dYdctrl.allocate(nCase, 4), dZdctrl.allocate(nCase, 4);

	dLdvel.allocate(nCase, 3), dMdvel.allocate(nCase, 3), dNdvel.allocate(nCase, 3);
	dLdomg.allocate(nCase, 3), dMdomg.allocate(nCase, 3), dNdomg.allocate(nCase, 3);
	dLdeul.allocate(nCase, 2), dMdeul.allocate(nCase, 2), dNdeul.allocate(nCase, 2);
	dLdctrl.allocate(nCase, 4), dMdctrl.allocate(nCase, 4), dNdctrl.allocate(nCase, 4);

	AMatrix.allocate(9, 9, nCase), BMatrix.allocate(9, 4, nCase);
	AMatrixRe.allocate(9, 9, nCase);
	ALongM.allocate(4, 4, nCase), ALateM.allocate(5, 5, nCase);
	ALongMRe.allocate(4, 4, nCase), ALateMRe.allocate(5, 5, nCase);

	flg = 0, jtype = jt; 
	niter.allocate(nCase, 1);
	path.clear();

	//StableInitPara.allocate(nCase, 2);
	//StableInitPara.input("StableInitPara.in");
	Vfs.input("Vfs-bo105.in");
	Consini.input(fn);

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
	path.clear();
}

void Jobs::SetSimCond(Copter &C, const int ic)
{
	double euler[3] = { 0,0,0 };
	double _uctrl[6] = { 0,0,0,0,0,0 };

	//C.vel_g[0] = C.RotorV[0].vtipa*Mus(ic);
	if (C.fuselage.si_unit)
		C.vel_g[0] = Vfs(ic)*0.5144444;
	else
		C.vel_g[0] = Vfs(ic);
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
		if (C.RotorV[i].adyna > 100)
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
		if (C.RotorV[i].adyna > 100)
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

void Jobs::PostProcess(Copter &C, string pf)
{
	if (path.empty())
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, sizeof(buffer), "%Y%m%d%I%M%S", timeinfo);
		std::string str(buffer);

		path = pf + '-' + std::to_string(int(C.mass)) + '-' + std::to_string(int(C.vel_g[0])) + '-' + str;

		mkdir(path.c_str());
	}

	// input parameters check
	for (int i = C.nfree - 1; i >= 0; i--)
		uctrl(0, i) = DEG(C.controls[i]);
	
	flightspeed2.allocate(1, 1);
	flightspeed2(0, 0) = C.vel_g[0];
	flightspeed2.output(path + "//" + pf + "-Vfs.output", 5);
	printf("Velocity check: %f, %f, %f \n", C.vel_g[0], C.vel_g[1], C.vel_g[2]);
	uctrl.output(path + "//" + pf + "-Controls.output", 6);

	C.AMatrix.output(path + "//" + pf + "-AMatrix.output", 10);
	C.BMatrix.output(path + "//" + pf + "-BMatrix.output", 10);
}

void Jobs::PostProcess(Copter &C, string pf, const int ic, const int s, const int e)
{
	if (flg++ < e - s)
	{
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl(ic, i) = DEG(C.controls[i]);

		//flightspeed2(0, ic) = C.vel_g[0];

		for (int j = 0; j < dXdvel.NJ; j++)
		{
			dXdvel(ic, j) = C.dXdvel(j);
			dYdvel(ic, j) = C.dYdvel(j);
			dZdvel(ic, j) = C.dZdvel(j);
			dXdomg(ic, j) = C.dXdomg(j);
			dYdomg(ic, j) = C.dYdomg(j);
			dZdomg(ic, j) = C.dZdomg(j);

			dLdvel(ic, j) = C.dLdvel(j);
			dMdvel(ic, j) = C.dMdvel(j);
			dNdvel(ic, j) = C.dNdvel(j);
			dLdomg(ic, j) = C.dLdomg(j);
			dMdomg(ic, j) = C.dMdomg(j);
			dNdomg(ic, j) = C.dNdomg(j);
		}
		for (int j = 0; j < dXdeul.NJ; j++)
		{
			dXdeul(ic, j) = C.dXdeul(j);
			dYdeul(ic, j) = C.dYdeul(j);
			dZdeul(ic, j) = C.dZdeul(j);

			dLdeul(ic, j) = C.dLdeul(j);
			dMdeul(ic, j) = C.dMdeul(j);
			dNdeul(ic, j) = C.dNdeul(j);
		}
		for (int j = 0; j < dXdctrl.NJ; j++)
		{
			dXdctrl(ic, j) = C.dXdctrl(j);
			dYdctrl(ic, j) = C.dYdctrl(j);
			dZdctrl(ic, j) = C.dZdctrl(j);

			dLdctrl(ic, j) = C.dLdctrl(j);
			dMdctrl(ic, j) = C.dMdctrl(j);
			dNdctrl(ic, j) = C.dNdctrl(j);
		}

		for (int j = 0; j < C.AMatrix.NJ; j++)
			for (int i = 0; i < C.AMatrix.NI; i++)
			{
				AMatrix(i, j, ic) = C.AMatrix(i, j);
				AMatrixRe(i, j, ic) = C.AMatrixRe(i, j);
			}
		for (int j = 0; j < C.BMatrix.NJ; j++)
			for (int i = 0; i < C.BMatrix.NI; i++)
				BMatrix(i, j, ic) = C.BMatrix(i, j);
		for(int j=0;j<4;j++)
			for (int i = 0; i < 4; i++)
			{
				ALongM(i, j, ic) = C.ALongM(i, j);
				ALongMRe(i, j, ic) = C.ALongMRe(i, j);
				ALateM(i, j, ic) = C.ALateM(i, j);
				ALateMRe(i, j, ic) = C.ALateMRe(i, j);
			}
	}

	if (path.empty())
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, sizeof(buffer), "%Y%m%d%I%M%S", timeinfo);
		std::string str(buffer);

		path = pf + '-' + std::to_string(int(C.mass)) + '-' + std::to_string(int(C.vel_g[0])) + '-' + str;

		mkdir(path.c_str());
	}

	if (flg == e - s)
	{
		uctrl.output(path + "//" + pf + "-Controls.output", 6);
		//flightspeed2.output(path + "//" + pf + "-Vfs.output", 5);

		AMatrix.output2(path + "//" + pf + "-AMatrix.output", 10);
		BMatrix.output2(path + "//" + pf + "-BMatrix.output", 10);
		AMatrixRe.output2(path + "//" + pf + "-AMatrixRe.output", 10);

		ALongM.output2(path + "//" + pf + "-ALongM.output", 10);
		ALongMRe.output2(path + "//" + pf + "-ALongMRe.output", 10);
		ALateM.output2(path + "//" + pf + "-ALateM.output", 10);
		ALateMRe.output2(path + "//" + pf + "-ALateMRe.output", 10);

		dXdvel.output(path + "//" + pf + "-dXdvel.output", 10);
		dYdvel.output(path + "//" + pf + "-dYdvel.output", 10);
		dZdvel.output(path + "//" + pf + "-dZdvel.output", 10);

		dXdomg.output(path + "//" + pf + "-dXdomg.output", 10);
		dYdomg.output(path + "//" + pf + "-dYdomg.output", 10);
		dZdomg.output(path + "//" + pf + "-dZdomg.output", 10);

		dXdeul.output(path + "//" + pf + "-dXdeul.output", 10);
		dYdeul.output(path + "//" + pf + "-dYdeul.output", 10);
		dZdeul.output(path + "//" + pf + "-dZdeul.output", 10);

		dLdvel.output(path + "//" + pf + "-dLdvel.output", 10);
		dMdvel.output(path + "//" + pf + "-dMdvel.output", 10);
		dNdvel.output(path + "//" + pf + "-dNdvel.output", 10);

		dLdomg.output(path + "//" + pf + "-dLdomg.output", 10);
		dMdomg.output(path + "//" + pf + "-dMdomg.output", 10);
		dNdomg.output(path + "//" + pf + "-dNdomg.output", 10);

		dLdeul.output(path + "//" + pf + "-dLdeul.output", 10);
		dMdeul.output(path + "//" + pf + "-dMdeul.output", 10);
		dNdeul.output(path + "//" + pf + "-dNdeul.output", 10);

		dXdctrl.output(path + "//" + pf + "-dXdctrl.output", 10);
		dYdctrl.output(path + "//" + pf + "-dYdctrl.output", 10);
		dZdctrl.output(path + "//" + pf + "-dZdctrl.output", 10);

		dLdctrl.output(path + "//" + pf + "-dLdctrl.output", 10);
		dMdctrl.output(path + "//" + pf + "-dMdctrl.output", 10);
		dNdctrl.output(path + "//" + pf + "-dNdctrl.output", 10);

		// save helicopter model as a .txt
		std::ifstream src("SimCase.cpp", std::ios::binary);
		std::ofstream dst(path + "//" + "simcase.txt", std::ios::binary);
		dst << src.rdbuf();
	}

	// vel test
	printf("Flag: %d\n", flg);
	printf("Vel g test: %f \n", C.vel_g[0]);
}

void Jobs::PostProcess(Copter &C, const int ic, const int s, const int e)
{
	double _beta[3], _beta_tr[3];
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

		C.RotorV[1].GetBeta(_beta_tr);
		for (int i = 0; i < 3; ++i)
			beta_tr(ic, i) = DEG(_beta_tr[i]);

		err(ic, 0) = C.sum_a1_del;
		err(ic, 1) = C.sum_a2_del;
		err(ic, 2) = C.max_c_del;
	}

	string prefix = "LF";
	if (path.empty())
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, sizeof(buffer), "%Y%m%d%I%M%S", timeinfo);
		std::string str(buffer);

		path = prefix + '-' + std::to_string(int(C.mass)) + '-' + std::to_string(int(C.amb.height)) + '-' + str;

		mkdir(path.c_str());
	}
		

	if (ic + 1 == e)
	{
		beta.output(path + "//" + prefix + "-MR-Betas.output", 4);
		_power.output(path + "//" + prefix + "-MR-Powers.output", 10);
		_torque.output(path + "//" + prefix + "-MR-Torques.output", 10);

		uctrl.output(path + "//" + prefix + "-Controls.output", 6);
		err.output(path + "//" + prefix + "-Errs.output", 10);

		beta_tr.output(path + "//" + prefix + "-TR-Betas.output", 4);
		_power_tr.output(path + "//" + prefix + "-TR-Powers.output", 10);
		_torque_tr.output(path + "//" + prefix + "-TR-Torques.output", 10);

		// save helicopter model as a .txt
		std::ifstream src("SimCase.cpp", std::ios::binary);
		std::ofstream dst(path + "//" + "simcase.txt", std::ios::binary);
		dst << src.rdbuf();
	}
	for (int i = C.RotorV.size() - 1; i>=0;--i) 
		C.RotorV[i].OutPutWake(ic);

	C.RotorV[0].DiskOutput(path + "//" + prefix + "-MR-" + std::to_string(ic));
	C.RotorV[1].DiskOutput(path + "//" + prefix + "-TR-" + std::to_string(ic));

	for (int i = C.RotorV.size() - 1; i >= 0; --i)
		C.RotorV[i].OutPutWake(ic);


	// vel test
	printf("Flag: %d\n", flg);
	printf("Vel g test: %f \n", C.vel_g[0]);
	printf("Mu g test: %f \n", C.vel_g[0] / (C.RotorV[0].vtipa * 0.3048*3.6));
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
	double _beta[3], _beta_tr[3];

	if (flg++ < e - s)
	{
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl(ic, i) = DEG(C.controls[i]);
		
		_power(ic, 0) = C.RotorV[0].power;
		_power(ic, 1) = C.RotorV[0].power_i;
		_power(ic, 2) = C.RotorV[0].power_o;
		_power(ic, 3) = C.RotorV[0].power_f;
		_power(ic, 4) = C.RotorV[0].power_c;
		_power(ic, 5) = C.RotorV[0].power_iid;
		_torque(ic, 0) = C.RotorV[0].torque;
		_torque(ic, 1) = C.RotorV[0].torque_i;
		_torque(ic, 2) = C.RotorV[0].torque_o;
		_torque(ic, 3) = C.RotorV[0].torque_f;
		_torque(ic, 4) = C.RotorV[0].torque_c;
		_torque(ic, 5) = C.RotorV[0].torque_iid;

		_power_tr(ic, 0) = C.RotorV[1].power;
		_power_tr(ic, 1) = C.RotorV[1].power_i;
		_power_tr(ic, 2) = C.RotorV[1].power_o;
		_power_tr(ic, 3) = C.RotorV[1].power_f;
		_power_tr(ic, 4) = C.RotorV[1].power_c;
		_power_tr(ic, 5) = C.RotorV[1].power_iid;
		_torque_tr(ic, 0) = C.RotorV[1].torque;
		_torque_tr(ic, 1) = C.RotorV[1].torque_i;
		_torque_tr(ic, 2) = C.RotorV[1].torque_o;
		_torque_tr(ic, 3) = C.RotorV[1].torque_f;
		_torque_tr(ic, 4) = C.RotorV[1].torque_c;
		_torque_tr(ic, 5) = C.RotorV[1].torque_iid;

		C.RotorV[0].GetBeta(_beta);
		for (int i = 0; i < 3; ++i)
			beta(ic, i) = DEG(_beta[i]);

		C.RotorV[1].GetBeta(_beta_tr);
		for (int i = 0; i < 3; ++i)
			beta_tr(ic, i) = DEG(_beta_tr[i]);

		err(ic, 0) = C.sum_a1_del;
		err(ic, 1) = C.sum_a2_del;
		err(ic, 2) = C.max_c_del;		

		niter(ic, 0) = C.Niter;
		//cout << niter(ic, 0) << endl;

		//C.RotorV[0].GetAirfm(_beta, _beta_tr);
		//cout << _beta[1] << endl;
	}

	string prefix = "LF";
	if (path.empty())
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, sizeof(buffer), "%Y%m%d%I%M%S", timeinfo);
		std::string str(buffer);

		path = prefix + '-' + std::to_string(int(C.mass)) + '-' + std::to_string(int(C.amb.height)) + '-' + str;

		mkdir(path.c_str());
	}

	if (flg == e - s)
	{
		beta.output(path + "//" + prefix + "-MR-Betas.output", 4);
		_power.output(path + "//" + prefix + "-MR-Powers.output", 10);
		_torque.output(path + "//" + prefix + "-MR-Torques.output", 10);

		uctrl.output(path + "//" + prefix + "-Controls.output", 6);
		err.output(path + "//" + prefix + "-Errs.output", 10);

		beta_tr.output(path + "//" + prefix + "-TR-Betas.output", 4);
		_power_tr.output(path + "//" + prefix + "-TR-Powers.output", 10);
		_torque_tr.output(path + "//" + prefix + "-TR-Torques.output", 10);

		niter.output(path + "//" + prefix + "-Niter.output", 1);

		// save helicopter model as a .txt
		std::ifstream src("SimCase.cpp", std::ios::binary);
		std::ofstream dst(path + "//" + "simcase.txt", std::ios::binary);
		dst << src.rdbuf();
	}

	C.RotorV[0].DiskOutput(path + "//" + prefix + "-MR-" + std::to_string(ic));
	C.RotorV[1].DiskOutput(path + "//" + prefix + "-TR-" + std::to_string(ic));

	for (int i = C.RotorV.size() - 1; i >= 0; --i)
		C.RotorV[i].OutPutWake(ic);


	// vel test
	printf("Flag: %d\n", flg);
	printf("Vel g test: %f \n", C.vel_g[0]);
}

void Jobs::PostProcessMP(Copter &C, const int ic, const int ip, const int s, const int e, const int np)
{
	double _beta[3], _beta_tr[3];
	if (flg++ < np*(e-s))
	{
		for (int i = C.nfree - 1; i >= 0; --i)
			uctrl3(ic, i, ip) = DEG(C.controls[i]);

		_power3(ic, 0, ip) = C.RotorV[0].power;
		_power3(ic, 1, ip) = C.RotorV[0].power_i;
		_power3(ic, 2, ip) = C.RotorV[0].power_o;
		_power3(ic, 3, ip) = C.RotorV[0].power_f;
		_power3(ic, 4, ip) = C.RotorV[0].power_c;
		_power3(ic, 5, ip) = C.RotorV[0].power_iid;
		_torque3(ic, 0, ip) = C.RotorV[0].torque;
		_torque3(ic, 1, ip) = C.RotorV[0].torque_i;
		_torque3(ic, 2, ip) = C.RotorV[0].torque_o;
		_torque3(ic, 3, ip) = C.RotorV[0].torque_f;
		_torque3(ic, 4, ip) = C.RotorV[0].torque_c;
		_torque3(ic, 5, ip) = C.RotorV[0].torque_iid;

		_power3_tr(ic, 0, ip) = C.RotorV[1].power;
		_power3_tr(ic, 1, ip) = C.RotorV[1].power_i;
		_power3_tr(ic, 2, ip) = C.RotorV[1].power_o;
		_power3_tr(ic, 3, ip) = C.RotorV[1].power_f;
		_power3_tr(ic, 4, ip) = C.RotorV[1].power_c;
		_power3_tr(ic, 5, ip) = C.RotorV[1].power_iid;
		_torque3_tr(ic, 0, ip) = C.RotorV[1].torque;
		_torque3_tr(ic, 1, ip) = C.RotorV[1].torque_i;
		_torque3_tr(ic, 2, ip) = C.RotorV[1].torque_o;
		_torque3_tr(ic, 3, ip) = C.RotorV[1].torque_f;
		_torque3_tr(ic, 4, ip) = C.RotorV[1].torque_c;
		_torque3_tr(ic, 5, ip) = C.RotorV[1].torque_iid;

		C.RotorV[0].GetBeta(_beta);
		for (int i = 0; i < 3; ++i)
			beta3(ic, i, ip) = DEG(_beta[i]);

		C.RotorV[1].GetBeta(_beta_tr);
		for (int i = 0; i < 3; ++i)
			beta3_tr(ic, i, ip) = DEG(_beta_tr[i]);

		err3(ic, 0, ip) = C.sum_a1_del;
		err3(ic, 1, ip) = C.sum_a2_del;
		err3(ic, 2, ip) = C.max_c_del;
		flightspeed2(ic, ip) = C.vel_g[0];
		wmega2(ic, ip) = C.RotorV[0].omega / param0(0);

		niter(ic, ip) = C.Niter;

		for (int i = 0; i < C.RotorV[0].bld.nitermax; i++)
		{
			sol3_mr(ic, i, ip) = C.RotorV[0].bld.sol(i);
			azmth_mr(ic, i, ip) = C.RotorV[0].bld.azmuth(i);
		}
	}

	string prefix = "RPM";
	if (path.empty())
	{
		time_t rawtime;
		struct tm * timeinfo;
		char buffer[80];

		time(&rawtime);
		timeinfo = localtime(&rawtime);

		strftime(buffer, sizeof(buffer), "%Y%m%d%I%M%S", timeinfo);
		std::string str(buffer);

		path = prefix + '-' + std::to_string(int(C.mass)) + '-' + std::to_string(int(C.amb.height)) + '-' + str;

		mkdir(path.c_str());
	}
	
	if (flg == np*(e - s))
	{
		beta3.output2(path + "//" + prefix + "-MR-Betas.output", 4);
		_power3.output2(path + "//" + prefix + "-MR-Powers.output", 10);
		_torque3.output2(path + "//" + prefix + "-MR-Torques.output", 10);

		uctrl3.output2(path + "//" + prefix + "-Controls.output", 6);
		err3.output2(path + "//" + prefix + "-Errs.output", 10);
		/*uctrl3.output2("uctrl_chd.output", 6);
		beta3.output2("beta_chd.output", 4);
		_power3.output2("power_chd.output", 10);
		_torque3.output2("torque_chd.output", 10);
		err3.output2("err_chd.output", 10);*/
		/*uctrl3.output2("uctrl_rad.output", 6);
		beta3.output2("beta_rad.output", 4);
		_power3.output2("power_rad.output", 10);
		_torque3.output2("torque_rad.output", 10);
		err3.output2("err_rad.output", 10);*/
		flightspeed2.output(path + "//" + prefix + "-Vfs.output", 5);
		wmega2.output(path + "//" + prefix + "-Omegas.output", 4);
		niter.output(path + "//" + prefix + "-Niter.output", 1);

		beta3_tr.output2(path + "//" + prefix + "-TR-Betas.output", 4);
		_power3_tr.output2(path + "//" + prefix + "-TR-Powers.output", 10);
		_torque3_tr.output2(path + "//" + prefix + "-TR-Torques.output", 10);

		sol3_mr.output2(path + "//" + prefix + "-MR-BetaSol.output", 10);
		azmth_mr.output2(path + "//" + prefix + "-MR-Azmth.output", 10);

		// save helicopter model as a .txt
		std::ifstream src("SimCase.cpp", std::ios::binary);
		std::ofstream dst(path + "//" + "simcase.txt", std::ios::binary);
		dst << src.rdbuf();
	}

	C.RotorV[0].DiskOutput(path + "//" + prefix + "-MR-" + std::to_string(ic) + "-" + std::to_string(ip));
	C.RotorV[1].DiskOutput(path + "//" + prefix + "-TR-" + std::to_string(ic) + "-" + std::to_string(ip));

	//for (int i = C.RotorV.size() - 1; i >= 0; --i)
	//	C.RotorV[i].OutPutWake(ic);

	printf("\n\n");
	printf("Count: %d, Case ID: (%d, %d) \n", flg, ic, ip);
	printf("Vel g test: %f \n", C.vel_g[0]);
	printf("Omg test: %f \n", C.RotorV[0].omega / param0(0));
}

void Jobs::ParamSweep(const Copter &C)
{
	switch (jtype)
	{
	case SimTrim:
		break;
	case RPMSwp:
		RPMs.allocate(2, nCase), param0.allocate(2);
		nParams = 5;
		RPMs.input("RPMs.in");
		param0(0) = C.RotorV[0].omega, param0(1) = C.RotorV[0].vtipa;
		break;
	case ChordSwp:
		Chds.allocate(2, nCase), param0.allocate(3);
		nParams = 20;
		Chds.input("Chds.in");
		param0(0) = C.RotorV[0].chord.v_p[0], param0(1) = C.RotorV[0].sigma, param0(2) = C.RotorV[0].gama;
		break;
	case RadiusSwp:
		Rads.allocate(2, nCase), param0.allocate(7);
		nParams = 2;
		Rads.input("Rads.in");
		param0(0) = C.RotorV[0].radius, param0(1) = C.RotorV[0].vtipa, param0(2) = C.RotorV[0].disk_A;
		param0(3) = C.RotorV[0].iflap, param0(4) = C.RotorV[0].m1;
		param0(5) = C.RotorV[0].sigma, param0(6) = C.RotorV[0].gama;
		break;
	default:
		break;
	}

	if (jtype > SimTrim)
	{
		uctrl3.allocate(nCase, 6, nParams), beta3.allocate(nCase, 3, nParams);
		_power3.allocate(nCase, 6, nParams), _torque3.allocate(nCase, 6, nParams);
		err3.allocate(nCase, 3, nParams);
		uctrl3_tr.allocate(nCase, 6, nParams), beta3_tr.allocate(nCase, 3, nParams);
		_power3_tr.allocate(nCase, 6, nParams), _torque3_tr.allocate(nCase, 6, nParams);
		err3_tr.allocate(nCase, 3, nParams);
		flightspeed2.allocate(nCase, nParams), wmega2.allocate(nCase, nParams);
		niter.deallocate();
		niter.allocate(nCase, nParams);
		sol3_mr.allocate(nCase, 720, nParams);
		azmth_mr.allocate(nCase, 720, nParams);
	}

}

void Jobs::UpdateParam(Copter &C, const int ic, const int ip)
{
	double temp;
	switch (jtype)
	{
	case SimTrim:
		break;
	case RPMSwp:
		temp = RPMs(0, ic) + (RPMs(1, ic) - RPMs(0, ic)) / (nParams - 1)*ip;
		C.RotorV[0].omega = param0(0)*temp;
		C.RotorV[0].vtipa = param0(1)*temp;
		break;
	case ChordSwp:
		temp = Chds(0, ic) + (Chds(1, ic) - Chds(0, ic)) / (nParams - 1)*ip;
		C.RotorV[0].chord.setvalue(param0(0)*temp);
		C.RotorV[0].gama = param0(2)*temp;
		C.RotorV[0].sigma = param0(1)*temp;
		break;
	case RadiusSwp:
		temp = Rads(0, ic) + (Rads(1, ic) - Rads(0, ic)) / (nParams - 1)*ip;
		C.RotorV[0].radius = param0(0)*temp, C.RotorV[0].vtipa = param0(1)*temp, C.RotorV[0].disk_A = param0(2)*temp*temp;
		C.RotorV[0].iflap = param0(3)*temp*temp, C.RotorV[0].m1 = param0(4)*temp;
		C.RotorV[0].sigma = param0(5)/temp, C.RotorV[0].gama = param0(6)*temp*temp;
		break;
	default:
		break;
	}
	
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

void CopterSolver::InitCopterSolver(void)
{
	err_a = 5.0e-3;
	err_c = 1.0e-2;
	epsilon = 0.01;
	nitermax = 30;
	sita_coll_max = RAD(35);
	sita_cycl_max = RAD(35);
	euler_max = RAD(35);
	niter_r.allocate(2), niter = 0, converge = false;
	sum_a1_del.allocate(nitermax), sum_a2_del.allocate(nitermax), max_c_del.allocate(nitermax);

	xdof = 8, udof = 4;
	dXdvel.allocate(3), dYdvel.allocate(3), dZdvel.allocate(3);
	dXdomg.allocate(3), dYdomg.allocate(3), dZdomg.allocate(3);
	dLdvel.allocate(3), dMdvel.allocate(3), dNdvel.allocate(3);
	dLdomg.allocate(3), dMdomg.allocate(3), dNdomg.allocate(3);
	dXdeul.allocate(2), dYdeul.allocate(2), dZdeul.allocate(2);
	dLdeul.allocate(2), dMdeul.allocate(2), dNdeul.allocate(2);
	dXdctrl.allocate(4), dYdctrl.allocate(4), dZdctrl.allocate(4);
	dLdctrl.allocate(4), dMdctrl.allocate(4), dNdctrl.allocate(4);
	AMatrix.allocate(xdof, xdof), BMatrix.allocate(xdof, udof);
}

void LevelFlight(void)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e;

	ul496.GetProb(0);
	ul496.GetModel();

	copter.InitRotorCraft(ul496);
	jobs.InitProject(SimTrim);

	s = 0, e = 1;// jobs.nCase;  //

	for (int i = s; i < e; ++i)
	{
		cout << i << endl;
		jobs.SetSimCond(copter, i);
		solver.CopterSimulation(copter);
		jobs.PostProcess(copter, i, s, e);
	}
}

void LevelFlightMP(int nth)
{
	Jobs jobs;
	//Model_UL496 ul496;
	Model_BO105 bo105;
	//Model_Puma330 puma330;
	Copter copter;
	CopterSolver solver;
	int s, e, i;
	
	jobs.InitProject(SimTrim);
	i = s = 0, e = jobs.nCase;

	bo105.GetProb(6, 0, FreeTrim1);
	bo105.GetModel();
	copter.InitRotorCraft(bo105);

#pragma omp parallel num_threads(nth) shared(s, e, jobs) firstprivate(i, solver, copter)
	{
		AeroDynType adytemp = copter.RotorV[0].adyna;
#pragma omp for
		for (i = s; i < e; i++)
		{
			if (i < 3)
				copter.RotorV[0].adyna = LinearInflow;
			else
				copter.RotorV[0].adyna = adytemp;
			
			jobs.SetSimCond(copter, i);
			solver.CopterSimulation(copter);
			jobs.PostProcessMP(copter, i, s, e);
		}
#pragma omp barrier
	}
}

void LinearModel(int nth)
{
	Jobs jobs;
	Model_BO105 bo105;
	Copter copter;
	CopterSolver solver;
	int s, e, i;

	jobs.InitProject(Stable, "Controls_Bo105.in");
	bo105.GetProb(6, 0, Transient);
	bo105.GetModel();
	copter.InitRotorCraft(bo105);

	i = s = 0, e = jobs.nCase; // 10; //

#pragma omp parallel num_threads(nth) shared(s, e, jobs) firstprivate(i, solver, copter)
	{
#pragma omp for
		for (int i = s; i < e; i++)
		{
			jobs.SetSimCond(copter, i);
			solver.CopterSimulation(copter);
			jobs.PostProcess(copter, "St", i, s, e);
		}
#pragma omp barrier
	}
}

void RPMSweepMP(int nth)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;
	int s, e, i, j, k, np, allcase;

	ul496.GetProb(0);
	ul496.GetModel();
	copter.InitRotorCraft(ul496);

	jobs.InitProject(RPMSwp);
	jobs.ParamSweep(copter);
	s = i = j = k = 0, e = jobs.nCase, np = jobs.nParams;
	allcase = np*(e - s);

#pragma omp parallel num_threads(nth) shared(s, e, jobs) firstprivate(i, j, k, solver, copter)
	{
#pragma omp for
		for (k = 0; k < allcase; k++)
		{
			i = k / np; // flight speed index
			j = k % np; // parameter index
			
			if (i < 4)
				copter.RotorV[1].adyna = Simple;

			jobs.SetSimCond(copter, i);
			jobs.UpdateParam(copter, i, j); 
			solver.CopterSimulation(copter);
			jobs.PostProcessMP(copter, i, j, s, e, np);
		}
#pragma omp barrier
	}
}

void MultiSweepMP(void)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;

	ul496.GetProb();
	ul496.GetModel();
	copter.InitRotorCraft(ul496);

	jobs.InitProject(RPMSwp);
}

void RPMSweep(const int ic, const int ip)
{
	Jobs jobs;
	Model_UL496 ul496;
	Copter copter;
	CopterSolver solver;

	ul496.GetProb();
	ul496.GetModel();
	copter.InitRotorCraft(ul496);

	jobs.InitProject(RPMSwp);
	jobs.ParamSweep(copter);

	jobs.SetSimCond(copter, ic);
	jobs.UpdateParam(copter, ic, ip);
	solver.CopterSimulation(copter);
	//jobs.PostProcessMP(copter, ic, ip, s, e, np);
}

void OPT_RPMSweep(void)
{
	Optimization opt;
	Model_UL496 ul496;
	Copter copter;
	int maxRepeat = 6;
	int s, e, np, allcase;

	ul496.GetProb();
	ul496.GetModel();
	copter.InitRotorCraft(ul496);

	opt.InitOptimization();

	s = 0, e = opt.range[0], np = opt.range[1];
	allcase = np*(e - s);

#pragma omp parallel num_threads(6) shared(opt) firstprivate(copter)
	{
		int k, kstart, kend;
		int Nthrds = omp_get_num_threads(), id = omp_get_thread_num();
		if (allcase / Nthrds == 0)
		{
			kstart = id*(allcase / Nthrds);
			kend = (id + 1)*(allcase / Nthrds);
		}
		else
		{
			kstart = id*(allcase / Nthrds + 1);
			kend = (id + 1)*(allcase / Nthrds + 1);
		}
		kend = Min(kend, allcase);

		int i, j;
		std::vector<int> _key(2), _keyup(2), _keydn(2), _ki(2);
		CopterSolver solver;

		i = j = 0;
		for (k = kstart; k < kend; k++)
		{
			i = k / np; // flight speed index
			j = k % np; // parameter index
			_key[0] = i, _key[1] = j;
			if (opt.haveKey(_key))
			{
				int icount = 0;
				while (!opt.fitMap[_key].converge && icount < maxRepeat)
				{
					opt.SetSimCond(copter, _key);
					opt.UpdateParam(copter, opt.deMap[_key]);
					solver.CopterSimulation(copter);
					opt.ResetParam(copter, opt.deMap[_key]);

					if (solver.converge)
					{
						opt.RecordFitness(copter, solver, opt.fitMap[_key]);
						opt.fitMap[_key].converge = true;
						opt.DPMonitor(opt.deMap[_key], opt.fitMap[_key].converge, id, copter.Niter);
					}
					else
					{
						opt.fitMap[_key].converge = false;
						opt.DPMonitor(opt.deMap[_key], opt.fitMap[_key].converge, id, copter.Niter);
						opt.UpdateDesigns(_key);
						++icount;
					}

				}

				if (icount == maxRepeat)
				{
					opt.fitMap.erase(_key);
					opt.deMap.erase(_key);
				}
			}
		}
#pragma omp barrier
	}
	opt.PostProcess();
}