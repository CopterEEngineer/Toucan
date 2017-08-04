#include "stdafx.h"
#include "Component.h"

Copter::Copter() :System()
{
	print_cons_on_screen("Copter");
	/* read config files */
	mass = 915;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			inmatx[i][j] = 0;
		}
	}
	inmatx[0][0] = 590 * 0.73696;
	inmatx[1][1] = 3098 * 0.73696;
	inmatx[2][2] = 2655 * 0.73696;

	rho = 0.002378;
	vsound = 1115.48;
	for (int i = 2; i >= 0; --i) {
		vel_c[i] = 0;
		omg_c[i] = 0;
		dvel_c[i] = 0;
		domg_c[i] = 0;
		vel_g[i] = 0;
		omg_g[i] = 0;
		dvel_g[i] = 0;
		domg_g[i] = 0;
		//sita[i] = 0;
	}
	//sita[0] = 10. / 180 * PI;
	//sita[1] = 2. / 180 * PI;
	//sita[3] = 0.0;
	myTYPE origin[3] = { 0,0,0 };
	myTYPE euler[3] = { 0,0,0 };
	/* read config files */
	Coordinate *base_ptr = NULL;
	refcoord.SetCoordinate(origin, euler, base_ptr);
	base_ptr = &refcoord;
	refcoord.SetBase(base_ptr);
	
	//SetStates();

	for (int i = 2; i >= 0; --i) {
		airforce_sigma[i] = 0;
		airmoment_sigma[i] = 0;
	}
	//cout << refcoord.base << endl;
#ifdef TEST_DATA
	vel_g[0] = 0.12*632.2455;
	vel_g[1] = 0*632.2455;
	vel_g[2] = 0*632.2455;
#endif // TEST_DATA

}


Copter::Copter(const Copter &H) :System(H)
{
	print_cpcons_on_screen("Copter");
	mass = H.mass;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			inmatx[i][j] = H.inmatx[i][j];
		}
	}

	rho = H.rho;
	vsound = H.vsound;
	refcoord = H.refcoord;
	for (int i = 2; i >= 0; --i) {
		vel_c[i] = H.vel_c[i];
		omg_c[i] = H.omg_c[i];
		dvel_c[i] = H.dvel_c[i];
		domg_c[i] = H.domg_c[i];
		vel_g[i] = H.vel_g[i];
		omg_g[i] = H.omg_g[i];
		dvel_g[i] = H.dvel_g[i];
		domg_g[i] = H.domg_g[i];
		airforce_sigma[i] = H.airforce_sigma[i];
		airmoment_sigma[i] = H.airmoment_sigma[i];
		//sita[i] = H.sita[i];
	}
}


Copter::~Copter()
{
	print_dscons_on_screen("Copter");
	mass = 0;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			inmatx[i][j] = 0;
		}
	}
	rho = 0;
	vsound = 0;
	for (int i = 2; i >= 0; --i) {
		vel_c[i] = 0;
		omg_c[i] = 0;
		dvel_c[i] = 0;
		domg_c[i] = 0;
		vel_g[i] = 0;
		omg_g[i] = 0;
		dvel_g[i] = 0;
		domg_g[i] = 0;
		airforce_sigma[i] = 0;
		airmoment_sigma[i] = 0;
		//sita[i] = 0;
	}
	//sita[3] = 0;
	myTYPE origin[3] = { 0,0,0 };
	myTYPE euler[3] = { 0,0,0 };

	refcoord.~Coordinate(); 
}


void Copter::SetStates(void) {
	// will update class members

	//for (int i = 2; i >= 0; --i) {
	//	vel_c[i] = vel_g[i];
	//	omg_c[i] = omg_g[i];
	//	dvel_c[i] = dvel_g[i];
	//	domg_c[i] = domg_g[i];
	//}
#ifdef USE_DOUBLE
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, vel_c, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, omg_c, 1);
#else
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, vel_c, 1);
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, omg_c, 1);
#endif // USE_DOUBLE

#ifdef FLIGHT_TRIM
	dvel_c[0] = dvel_c[1] = dvel_c[2] = 0;
	domg_c[0] = domg_c[1] = domg_c[2] = 0;
#else
	;
#endif // FLIGHT_TRIM

#ifdef TEST_MODE
	printf("Copter states.\n");
	printf("Velocity: %f, %f, %f \n", vel_c[0], vel_c[1], vel_c[2]);
	printf("Omega: %f, %f, %f \n", omg_c[0], omg_c[1], omg_c[2]);
	printf("dvelocity: %f, %f, %f \n", dvel_c[0], dvel_c[1], dvel_c[2]);
	printf("domega: %f, %f, %f \n", domg_c[0], domg_c[1], domg_c[2]);
	cout << vel_c[0] * vel_c[0] + vel_c[1] * vel_c[1] + vel_c[2] * vel_c[2] << endl;
	cout << omg_c[0] * omg_c[0] + omg_c[1] * omg_c[1] + omg_c[2] * omg_c[2] << endl;
	cout << dvel_c[0] * dvel_c[0] + dvel_c[1] * dvel_c[1] + dvel_c[2] * dvel_c[2] << endl;
	cout << domg_c[0] * domg_c[0] + domg_c[1] * domg_c[1] + domg_c[2] * domg_c[2] << endl;
	cout << endl << endl;

#endif // TEST_MODE
}


void Copter::SetStates(myTYPE v[3], myTYPE w[3], myTYPE dv[3], myTYPE dw[3]) {
	// will not update class members

	for (int i = 2; i >= 0; --i) {
		v[i] = w[i] = dv[i] = dw[i] = 0;
	}
#ifdef USE_DOUBLE
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, v, 1);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, w, 1);
#else
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, v, 1);
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, w, 1);
#endif // USE_DOUBLE

#ifdef FLIGHT_TRIM
	dv[0] = dv[1] = dv[2] = 0;
	dw[0] = dw[1] = dw[2] = 0;
#else
	;
#endif // FLIGHT_TRIM

	//for (int i = 2; i >= 0; --i) {
	//	vel_c[i] = v[i];
	//	omg_c[i] = w[i];
	//	dvel_c[i] = dv[i];
	//	domg_c[i] = dw[i];
	//}

#ifdef TEST_MODE
	printf("Copter states.\n");
	printf("Velocity: %f, %f, %f \n", v[0], v[1], v[2]);
	printf("Omega: %f, %f, %f \n", w[0], w[1], w[2]);
	printf("dvelocity: %f, %f, %f \n", dv[0], dv[1], dv[2]);
	printf("domega: %f, %f, %f \n", dw[0], dw[1], dw[2]);
	cout << endl << endl;

#endif // TEST_MODE
}


void Copter::Assemble(myTYPE f[3], myTYPE m[3]) {
	for (int i = 2; i >= 0; --i) {
		airforce_sigma[i] = f[i];
		airmoment_sigma[i] = m[i];
	}
}


void Copter::SetDerivs(myTYPE *f, myTYPE *m)
{
	// will update class member
	_setderivs(dvel_c, domg_c, f, m, refcoord);
}


void Copter::SetDerivs(myTYPE dv[3], myTYPE dw[3], myTYPE * f, myTYPE * m)
{
	// will not update class member
	dv[0] = dv[1] = dv[2] = 0;
	dw[0] = dw[1] = dw[2] = 0;
	_setderivs(dv, dw, f, m, refcoord);

}


void Copter::SetDerivs(void)
{
	// will update class member
	_setderivs(dvel_c, domg_c, airforce_sigma, airmoment_sigma, refcoord);
#ifdef TEST_MODE
	printf("%f, %f, %f\n", dvel_c[0], dvel_c[1], dvel_c[2]);
	printf("%f, %f, %f\n", domg_c[0], domg_c[1], domg_c[2]);

#endif // TEST_MODE

}


//void Copter::SetDerivs(myTYPE dv[3], myTYPE dw[3])
//{
//	dv[0] = dv[1] = dv[2] = 0;
//	dw[0] = dw[1] = dw[2] = 0;
//	_setderivs(dv, dw, airforce_sigma, airmoment_sigma, refcoord);
//}


void Copter::functest(void) 
{
	myTYPE dv[3] = { 0 };
	myTYPE dw[3] = { 0 };

	_setderivs(dv, dw, airforce_sigma, airmoment_sigma, refcoord);
	//cout << "test successfully." << endl;
}



Component::Component() :Copter() {
	print_cons_on_screen("Component");
	for (int i = 2; i >= 0; --i) {
		vel[i] = 0;
		omg[i] = 0;
		dvel[i] = 0;
		domg[i] = 0;
		airforce[i] = 0;
		airmoment[i] = 0;
		airforce_cg[i] = 0;
		airmoment_cg[i] = 0;
		airforce_sigma[i] = 0;
		airmoment_sigma[i] = 0;
	}
	//cout << refcoord.base << endl;

}


Component::Component(const Component &C) :Copter(C) {
	print_cpcons_on_screen("Component");
	for (int i = 2; i >= 0; --i) {
		vel[i] = C.vel[i];
		omg[i] = C.omg[i];
		dvel[i] = C.dvel[i];
		domg[i] = C.domg[i];
		airforce[i] = C.airforce[i];
		airmoment[i] = C.airmoment[i];
		airforce_cg[i] = C.airforce_cg[i];
		airmoment_cg[i] = C.airmoment_cg[i];
	}
}


Component::~Component() {
	print_dscons_on_screen("Component");
	for (int i = 2; i >= 0; --i) {
		vel[i] = 0;
		omg[i] = 0;
		dvel[i] = 0;
		domg[i] = 0;
		airforce[i] = 0;
		airmoment[i] = 0;
		airforce_cg[i] = 0;
		airmoment_cg[i] = 0;
	}
}


void Component::SetStates(const myTYPE *vc, const myTYPE *wc, const myTYPE *dvc, const myTYPE *dwc)
{
	_setstates(vel, omg, dvel, domg, vc, wc, dvc, dwc);
#ifdef OUTPUT_MODE
	printf("Component states.\n");
	printf("Velocity: %f, %f, %f \n", vel[0], vel[1], vel[2]);
	printf("Omega: %f, %f, %f \n", omg[0], omg[1], omg[2]);
	printf("dvelocity: %f, %f, %f \n", dvel[0], dvel[1], dvel[2]);
	printf("domega: %f, %f, %f \n", domg[0], domg[1], domg[2]);
	//cout << vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2] << endl;
	//cout << omg[0] * omg[0] + omg[1] * omg[1] + omg[2] * omg[2] << endl;
	//cout << dvel[0] * dvel[0] + dvel[1] * dvel[1] + dvel[2] * dvel[2] << endl;
	//cout << domg[0] * domg[0] + domg[1] * domg[1] + domg[2] * domg[2] << endl;

	cout << endl << endl;

#endif // OUTPUT_MODE

}


void Component::SetAirfm_cg(const Coordinate *base) {
	myTYPE tcb[3][3];
	myTYPE rcb[3];
	myTYPE temp_cross[3] = { 0,0,0 };

	refcoord.Transfer(tcb, rcb, refcoord, *base);

#ifdef USE_DOUBLE
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airforce, 1, 0, airforce_cg, 1);
	
	Cross(airmoment_cg, airforce_cg, rcb);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airmoment, 1, 1, airmoment_cg, 1);
	//vdAdd(3, temp_cross, airmoment_cg, airmoment_cg);

#else
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airforce, 1, 0, airforce_cg, 1);

	Cross(airmoment_cg, airforce_cg, rcb);
	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airmoment, 1, 1, airmoment_cg, 1);

#endif // USE_DOUBLE
}


void Component::Assemble(const std::vector<std::unique_ptr<Component>> &C, const Coordinate *base) {
	_assemble(airforce_sigma, airmoment_sigma, C, base);
}



Fuselage::Fuselage() :Component() {
	print_cons_on_screen("Fuselage");
	// read config file
#ifdef UL496
	dragA = 0.0315*DISK_A;

#endif // UL496
	//cout << refcoord.base << endl;

}


Fuselage::Fuselage(const Fuselage &F) :Component(F) {
	print_cpcons_on_screen("Fuselage");

	dragA = F.dragA;
}


Fuselage::~Fuselage() {
	print_dscons_on_screen("Fuselage");

	dragA = 0;
}


inline void Fuselage::SetAirfm(void) {
	airforce[0] = -0.5*vel[0] * vel[0] * rho * dragA;
	airforce[1] = airforce[2] = 0;
	airmoment[0] = airmoment[1] = airmoment[2] = 0;
}


Wing::Wing(char *s, int num) :Component() {
	//cout << "Wing construtor deriving from Copter." << endl;
	//read config files
	myTYPE origin[3], euler[3];


	type = new char[10];
	strcpy(type, s);
	//cout << refcoord.base << endl;

#ifdef UL496
	if (!strcmp(s, "hor")) {
		origin[0] = -7.87;
		origin[1] = 0;
		origin[2] = 0;
		euler[0] = 0;
		euler[1] = -4 * PI / 180;
		euler[2] = 0;
		refcoord.SetCoordinate(origin, euler, refcoord.base); // base copy constructor give correct base if no need to revise

		a0 = 5.73;
		cd0 = 0.04;
		cd1 = 0;
		cd2 = 0;
		span = 3.75;
		chord = 1;
		taper = 1;

	}
	else if (!strcmp(s, "ver")) {
		span = 1.33;
		chord = 1.13;
		taper = 0.51;
		a0 = 5.3;
		cd0 = 0.0105;
		cd1 = 0;
		cd2 = 0.01325;
		if (num == 1) {
			origin[0] = -(7.87 + 0.164);
			origin[1] = 3.75 / 2;
			origin[2] = 0;
			euler[0] = PI / 2;
			euler[1] = 0;
			euler[2] = 5 * PI / 180;
			refcoord.SetCoordinate(origin, euler, refcoord.base);
		}
		else if (num == 2) {
			origin[0] = -(7.87 - 0.164);
			origin[1] = -3.75 / 2;
			origin[2] = 0;
			euler[0] = PI / 2;
			euler[1] = 0;
			euler[2] = 5 * PI / 180;
			refcoord.SetCoordinate(origin, euler, refcoord.base);
		}
		else {
			print_wrong_msg("Wrong wing definition.");
		}



	}
	else {
		print_wrong_msg("Wrong wing type.");
	}
#endif // UL496

}


Wing::Wing(const Wing &W) :Component(W) {
	strcpy(type, W.type);
	a0 = W.a0;
	cd0 = W.cd0;
	cd1 = W.cd1;
	cd2 = W.cd2;
	span = W.span;
	chord = W.chord;
	taper = W.taper;
}


Wing::~Wing() {
	delete[] type;
#ifdef UL496
	a0 = 0;
	cd0 = 0;
	cd1 = 0;
	cd2 = 0;
	span = 0;
	chord = 0;
	taper = 0;
#endif // UL496

}


inline void Wing::SetAirfm(void) {
	myTYPE stot, ar, a03d;
	myTYPE aoa, cl, cd, vel2;
#ifdef UL496
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
	airforce[0] = -0.5 * rho * stot * vel2 * cd;
	airforce[1] = 0;
	airforce[2] = 0.5 * rho * stot * vel2 * cl;
	airmoment[0] = 0;
	airmoment[1] = 0;
	airmoment[2] = 0;
#endif
}


Rotor::Rotor(const char *s) :Component() {
	//print_cons_on_screen("Rotor");
	// read config file
	type = new char[10];
	strcpy(type, s);

	a0 = 5.7;
	if (!strcmp(s, "main")) {
		nf = 72;
		ns = 41;
		ni = 10;

		if (DISABLE_REVISE_SIZE) {
			if (nf*ns > MAX_SIZE) {
				printf("nf %d and ns %d make beyond of %d \n", nf, ns, MAX_SIZE);
				exit(EXIT_FAILURE);
			}
		}

		kwtip = 5;
		kwrot = 5;
		nk = kwtip * nf;
		eflap = 0;
		khub = 0;
		del = 0;
		pitchroot = 0.0;
		radius = 11.5;
		bt = 0.98;
		rroot = 0.15;
		precone = 3.0 * PI / 180;
		omega = 54.977871437821380;
		vtipa = omega * radius;
		outboard = 0.3;
		rc0 = 0.004852173913043;

		cltc.allocate(CL_I, CL_J);
		cdtc.allocate(CD_I, CD_J);
		cltc.input("vr7_cl_c81.txt");
		cdtc.input("vr7_cd_c81.txt");

		iflap = 24.93;
		m1 = 3.25;
		sigma = 0.0309;
		gama = 5.01 * rho / 0.002378; //

		//chord.allocate(nf, ns);
		//sweep.allocate(nf, ns);
		//twist.allocate(nf, ns);
		chord.allocate(ns);
		sweep.allocate(ns);
		twist.allocate(ns);
		azstation.allocate(nf, ns);
		rastation.allocate(nf, ns);
		chord.setvalue(0.558);
		sweep.setvalue(0);
		//chord.input("chord.txt");
		//twist.input("twist.txt");
		//azstation.input("azstation.txt");
		//rastation.input("rastation.txt");
		myTYPE temp_twist, temp_azimuth, temp_station;
		for (int j = ns - 1; j >= 0; --j) {
			temp_station = rroot + j*(1.0 - rroot) / (ns - 1); // uniform seperated accepted so far.
			//temp_twist = (-8 * PI / 180) * j / (ns - 1);
			twist(j) = (-8 * PI / 180) * j / (ns - 1);
			for (int i = nf - 1; i >= 0; --i) {
				//twist(i, j) = temp_twist;
				rastation(i, j) = temp_station;
				azstation(i, j) = i * 2 * PI / nf;
			}
		}

		myTYPE origin[3], euler[3];
		origin[0] = -0.0625;
		origin[1] = 0;
		origin[2] = -4.217;
		euler[0] = 0;
		euler[1] = -3.0 * PI / 180;
		euler[2] = 0;
		hubfxcoord.SetCoordinate(origin, euler, refcoord.base);
		refcoord.SetCoordinate(origin, euler, refcoord.base);

		// referred to hubfxcoord
		origin[0] = origin[1] = origin[2] = 0;
		euler[0] = 0;
		euler[1] = -PI;
		euler[2] = 0;
		hubrtcoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubrtcoord
		origin[0] = eflap; //定义扭转，摆振，挥舞的三个桨叶坐标
		origin[1] = 0;
		origin[2] = 0;
		euler[0] = 0;
		euler[1] = -precone - beta[1];
		euler[2] = 0;
		bladecoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubfxcoord
		origin[0] = 0;
		origin[1] = precone * (1 - eflap);
		origin[2] = 0;
		euler[0] = euler[1] = euler[2] = 0;
		tppcoord.SetCoordinate(origin, euler, refcoord.base);


		// initialize member variables to zero
		mul = 0;
		lambdi_ag = sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI/radius / radius / vtipa / vtipa);
		lambdt_ag = lambdh_ag = lambdi_ag;
		power = 0;
		torque = 0;

		beta[1] = beta[2] = 0;
		beta[0] = precone;

		bflap.allocate(nf, ns); // flap
		dbflap.allocate(nf, ns);
		sfth.allocate(nf, ns);  // pitch 
		ut.allocate(nf, ns);    // air velocity
		un.allocate(nf, ns);
		up.allocate(nf, ns);
		ua.allocate(nf, ns);
		ma_n.allocate(nf, ns);
		incidn.allocate(nf, ns);// AOA
		cl.allocate(nf, ns);    // air coefficients
		cd.allocate(nf, ns);
		cirlb.allocate(nf, ns); // circulation
		tipstr.allocate(nk, nf);
		lambdi.allocate(nf, ns);// induced velocity
		lambdx.allocate(nf, ns);
		lambdy.allocate(nf, ns);
		lambdh.allocate(nf, ns);
		lambdt.allocate(nf, ns);

		tipgeometry.allocate(nk, nf, 3);
		bladedeform.allocate(nf, ns, 3);


		ristation.allocate(nf, ni);
		
		sita[0] = sita[1] = sita[2] = sita[3] = 0.0;

#ifdef TEST_DATA
		mul = 0.1765;
		vel[0] = mul * vtipa;
		vel[1] = vel[2] = 0;
		//lambdh.setvalue(-0.0725);
		lambdh.input("lambdh.txt");
		lambdh.output("lambdh.output", 4);
		sita[0] = 24.1403 * PI / 180;
		sita[1] = 0.5731 * PI / 180;
		sita[2] = -5.7843 * PI / 180;
		sita[3] = 0.0;
		beta[1] = 1.5647*PI / 180;
		beta[2] = 0.5840*PI / 180;
		for (int j = ns - 1; j >= 0; --j) {
			twist(j) = (-22 * PI / 180) * j / (ns - 1);
		}
#endif // TEST_DATA
	}
	else if (!strcmp(s, "tail")) {
		nf = 1;
		ns = 1;
		kwtip = 0;
		kwrot = 0;
		nk = 1;
		eflap = 0;
		khub = 0;
		del = 0;
		pitchroot = 0;
		radius = 1.3;
		bt = 0.98;
		rroot = 0.15;
		precone = 0;
		omega = 54.977871437821380 * 9.4912;
		vtipa = omega * radius;
		outboard = 0;
		rc0 = 0;

		cltc.allocate(nf, ns);
		cdtc.allocate(nf, ns);

		iflap = 0.0027;
		m1 = 0.00312;
		sigma = 0.1636;
		gama = 0; //

		//chord.allocate(nf, ns);
		//sweep.allocate(nf, ns);
		//twist.allocate(nf, ns);
		chord.allocate(ns);
		sweep.allocate(ns);
		twist.allocate(ns);
		azstation.allocate(nf, ns);
		rastation.allocate(nf, ns);
		chord.setvalue(0.167);

		myTYPE origin[3], euler[3];
		origin[0] = -13.48;
		origin[1] = 0;
		origin[2] = 1;
		euler[0] = PI / 2;
		euler[1] = 0;
		euler[2] = 0;
		hubfxcoord.SetCoordinate(origin, euler, refcoord.base);
		refcoord.SetCoordinate(origin, euler, refcoord.base);

		// referred to hubfxcoord
		origin[0] = origin[1] = origin[2] = 0;
		euler[0] = 0;
		euler[1] = -PI;
		euler[2] = 0;
		hubrtcoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubrtcoord
		origin[0] = eflap; //定义扭转，摆振，挥舞的三个桨叶坐标
		origin[1] = 0;
		origin[2] = 0;
		euler[0] = 0;
		euler[1] = -precone;
		euler[2] = 0;
		bladecoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubfxcoord
		origin[0] = 0;
		origin[1] = precone * (1 - eflap);
		origin[2] = 0;
		euler[0] = euler[1] = euler[2] = 0;
		tppcoord.SetCoordinate(origin, euler, refcoord.base);


		// initialize member variables to zero
		mul = 0;
		lambdi_ag = sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI/radius / radius / vtipa / vtipa);
		lambdt_ag = lambdh_ag = lambdi_ag;
		power = 0;
		torque = 0;

		beta[1] = beta[2] = 0;
		beta[0] = precone;

		bflap.allocate(nf, ns); // flap
		dbflap.allocate(nf, ns);
		sfth.allocate(nf, ns);  // pitch 
		ut.allocate(nf, ns);    // air velocity
		un.allocate(nf, ns);
		up.allocate(nf, ns);
		ua.allocate(nf, ns);
		ma_n.allocate(nf, ns);
		incidn.allocate(nf, ns);// AOA
		cl.allocate(nf, ns);    // air coefficients
		cd.allocate(nf, ns);
		cirlb.allocate(nf, ns); // circulation
		tipstr.allocate(nk, nf);
		lambdi.allocate(nf, ns);// induced velocity
		lambdx.allocate(nf, ns);
		lambdy.allocate(nf, ns);

		tipgeometry.allocate(nk, nf, 3);
		bladedeform.allocate(nf, ns, 3);

		sita[0] = sita[1] = sita[2] = sita[3] = 0.0;
#ifdef TEST_DATA
		sita[3] = 10 * PI / 180;
#endif // TEST_DATA

	}
	else {
		print_wrong_msg("Wrong rotor type. Please input main, or tail.");
	}


}


Rotor::Rotor(const Rotor &R) {
	cout << "Rotor copy construtor." << endl;
	strcpy(type, R.type);
	kwtip = R.kwtip;
	kwrot = R.kwrot;
	nk = R.nk;
	nf = R.nf;
	ns = R.ns;
	ni = R.ni;
	nbn = R.nbn;
	naf = R.naf;
	nnr = R.nnr;
	eflap = R.eflap;
	khub = R.khub;
	del = R.del;
	pitchroot = R.pitchroot;
	radius = R.radius;
	bt = R.bt;
	rroot = R.rroot;
	precone = R.precone;
	omega = R.omega;
	sigma = R.sigma;
	gama = R.gama;
	a0 = R.a0;
	rastation = R.rastation;
	ristation = R.ristation;
	azstation = R.azstation;
	chord = R.chord;
	twist = R.twist;
	sweep = R.sweep;
	iflap = R.iflap;
	m1 = R.m1;
	rtip = R.rtip;
	rc0 = R.rc0;
	outboard = R.outboard;
	mul = R.mul;
	vtipa = R.vtipa;
	beta[0] = R.beta[0];
	beta[1] = R.beta[1];
	beta[2] = R.beta[2];
	cltc = R.cltc;
	cdtc = R.cdtc;
	cmtc = R.cmtc;
	lambdi_ag = R.lambdi_ag;
	lambdt_ag = R.lambdt_ag;
	lambdh_ag = R.lambdh_ag;
	power = R.power;
	torque = R.torque;
	bflap = R.bflap;
	dbflap = R.dbflap;
	sfth = R.sfth;
	ut = R.ut;
	un = R.un;
	up = R.up;
	ua = R.ua;
	ma_n = R.ma_n;
	incidn = R.incidn;
	cl = R.cl;
	cd = R.cd;
	cirlb = R.cirlb;
	lambdi = R.lambdi;
	lambdh = R.lambdh;
	lambdt = R.lambdt;
	lambdx = R.lambdx;
	lambdy = R.lambdy;
	tipstr = R.tipstr;
	rotstr = R.rotstr;
	shdstr = R.shdstr;
	trlstr = R.trlstr;
	bladedeform = R.bladedeform;
	tipgeometry = R.tipgeometry;
	hubfxcoord = R.hubfxcoord;
	hubrtcoord = R.hubrtcoord;
	bladecoord = R.bladecoord;
	tppcoord = R.tppcoord;
	sita[0] = R.sita[0];
	sita[1] = R.sita[1];
	sita[2] = R.sita[2];
	sita[3] = R.sita[3];
}


Rotor::~Rotor() {
	cout << "Rotor destructor." << endl;
	delete[] type;
	kwtip = kwrot = nk = nf = ns = ni = nbn = naf = nnr = 0;
	eflap = khub = del = pitchroot = radius = bt = rroot = 0;
	precone = omega = 0;
	sigma = gama = a0 = 0;
	iflap = m1 = rtip = rc0 = outboard = 0;
	mul = vtipa = 0;
	lambdi_ag = lambdt_ag = lambdh_ag = 0;
	power = torque = 0;
	beta[0] = beta[1] = beta[2] = 0;
	rastation.deallocate();
	ristation.deallocate();
	azstation.deallocate();
	chord.deallocate();
	twist.deallocate();
	sweep.deallocate();
	cltc.deallocate();
	cdtc.deallocate();
	cmtc.deallocate();
	bflap.deallocate();
	dbflap.deallocate();
	sfth.deallocate();
	ut.deallocate();
	un.deallocate();
	up.deallocate();
	ua.deallocate();
	ma_n.deallocate();
	incidn.deallocate();
	cl.deallocate();
	cd.deallocate();
	cirlb.deallocate();
	lambdi.deallocate();
	lambdh.deallocate();
	lambdt.deallocate();
	lambdx.deallocate();
	lambdy.deallocate();
	tipstr.deallocate();
	rotstr.deallocate();
	shdstr.deallocate();
	trlstr.deallocate();
	bladedeform.deallocate();
	tipgeometry.deallocate();
	hubfxcoord.~Coordinate();
	hubrtcoord.~Coordinate();
	bladecoord.~Coordinate();
	tppcoord.~Coordinate();

	sita[0] = sita[1] = sita[2] = sita[3] = 0;
}


void Rotor::_aerodynacoef(Matrix1<myTYPE> &_cl, Matrix1<myTYPE> &_cd, Matrix1<myTYPE>& incidn, Matrix1<myTYPE>& ma_n)
{
#ifdef UL496	
	_cl = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn / PI*180.0, ma_n);
	_cd = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn / PI*180.0, ma_n);
#elif UH60A
	;
#endif
}


void Rotor::_aerodynacoef(Matrix2<myTYPE> &_cl, Matrix2<myTYPE> &_cd, Matrix2<myTYPE>& incidn, Matrix2<myTYPE>& ma_n)
{
#ifdef UL496	
	_cl = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn / PI*180.0, ma_n);
	_cd = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn / PI*180.0, ma_n);
#elif UH60A
	;
#endif
}


void Rotor::_flapmotionrt(void) {
	myTYPE m, c, k, q0, dq0, pho, dt, dff, err_b;
	int nitermax, niter;
	Matrix1<myTYPE> sol;

	m = 2 * iflap;
	c = 0;
	k = 2 * omega*omega*iflap + 2 * khub;

	// solve flap at rotation referrence
	niter = 0;
	dff = 15; //digit can seperate pi/2 normally
	nitermax = 30*360/dff;
	sol.allocate(nitermax);


	q0 = beta[0] - sita[2];
	dq0 = omega*sita[1];
	dt = dff / 180 * PI * omega;
	pho = 1;
	err_b = 0.001;
#ifdef TEST_MODE
	//err_b = 0;
#endif // TEST_MODE

	GenArf(sol, beta, niter, m, c, k, q0, dq0, dff, nitermax, pho, err_b);

#ifdef OUTPUT_MODE
	cout << (niter*dff/360) << endl;
	printf("%f, %f\n\n", beta[1]*180/PI, beta[2]*180/PI);
#endif // OUTPUT_MODE

#ifdef OUTPUT_MODE
	int nperd = 360 / dff;
	sol.output("bflap_rt.output",6);
	//sol.output(niter - nperd + 1, niter, 3);
#endif // OUTPUT_MODE

}


void Rotor::_bladeCSD(void)
{
	;
}


void Rotor::SetCtrl(myTYPE * xctrl, const int n)
{
	if (!strcmp(type, "main")) {
#ifdef _DEBUG
		if (n != 3) { print_wrong_msg("Undefined Control variations for main rotor."); }
#endif // _DEBUG

		for (int i = 0; i < 3; ++i) {
			sita[i] = xctrl[i];
		}
		sita[3] = 0;
	}
	else{
#ifdef _DEBUG
		if (n != 1) { print_wrong_msg("Undefined Control variations for tail rotor."); }
#endif // _DEBUG
		for (int i = 0; i < n; ++i) {
			sita[i+3] = xctrl[i];
		}
		sita[0] = sita[1] = sita[2] = 0;
	}
#ifdef TEST_MODE
	printf("%s set control variations successfully.", type);

#endif // TEST_MODE

}


void Rotor::GetCtrl(myTYPE * xctrl, const int n)
{
	if (!strcmp(type, "main")) {
#ifdef _DEBUG
		if (n != 3) { print_wrong_msg("Undefined Control variations for main rotor."); }
#endif // _DEBUG

		for (int i = 0; i < 3; ++i) {
			xctrl[i] = sita[i];
		}
	}
	else {
#ifdef _DEBUG
		if (n != 1) { print_wrong_msg("Undefined Control variations for tail rotor."); }
#endif // _DEBUG
		for (int i = 0; i < n; ++i) {
			xctrl[i + 3] = sita[i + 3];
		}
	}
}


void Rotor::SetCoordBase(void)
{
	hubfxcoord.SetBase(refcoord.base);
	hubrtcoord.SetBase(&hubfxcoord);
	bladecoord.SetBase(&hubrtcoord);
	tppcoord.SetBase(&hubfxcoord);
#ifdef TEST_MODE
	cout << hubfxcoord.base << endl;
	cout << hubrtcoord.base << endl;
	cout << bladecoord.base << endl;
	cout << tppcoord.base << endl;
#endif // TEST_MODE

}


void Rotor::BladeDynamics(void)
{
	myTYPE euler_temp[3];
	_flapmotionrt();
	
	euler_temp[0] = -beta[2];
	euler_temp[1] = -beta[1];
	euler_temp[2] = 0.0;
	tppcoord.SetCoordinate(euler_temp, "euler");
}


void Rotor::AvrgInducedVel(void)
{
	int itermax = 20;
	int i = 0;
	Matrix1<myTYPE> lambtpp(itermax);
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	//myTYPE velind[3] = { 0.0,0.0,0.0 };

	if (!strcmp(type, "main")) {

		lambtpp.v_p[0] = sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI / radius / radius / vtipa / vtipa);
		//lambtpp.v_p[0] = lambdt_ag;
		for (i = 1; i < itermax - 1; ++i) {

			BladeDynamics();
#ifdef USE_DOUBLE
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *tppcoord.Ttransf, 3, vel, 1, 0, veltpp, 1);
#else
			cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *tppcoord.Ttransf, 3, vel, 1, 0, veltpp, 1);
#endif // USE_DOUBLE
			mul = vel[0] / vtipa;

			//SetAirfm();
			_setairfm_sp(airforce, airmoment);
			lambdi_ag = -airforce[2] / (2 * rho*PI*radius*radius*vtipa*vtipa);
			lambdi_ag /= sqrt(mul * mul + lambtpp(i - 1)*lambtpp(i - 1)); 
			lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
			lambtpp(i) = lambdt_ag;

			//velind[2] = lambdi_ag;
			//Msolver(*temp, velind, 3, 1);
			lambdh_ag = lambdi_ag * tppcoord.Ttransf[2][2] - vel[2] / vtipa;
			lambdh.setvalue(lambdh_ag);

			if (Abs(lambtpp(i) / lambtpp(i - 1) - 1) < err_w) { break; }

		}
		//lambtpp.outputs("_lambtpp.output", 4);
	}
	else {
		//lambtpp.v_p[0] = lambdt_ag;
		lambtpp.v_p[0] = sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI / radius / radius / vtipa / vtipa);
		mul = vel[0] / vtipa;
		for (i = 1; i < itermax - 1; ++i) {
			
			airforce[2] = sigma*a0*0.5*(sita[3] / 3 * (1 + 1.5*mul * mul) - 0.5*lambtpp(i - 1));
			lambdi_ag = airforce[2] / (2 * sqrt(mul*mul + lambtpp(i - 1)*lambtpp(i - 1)));
			lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
			lambtpp(i) = lambdt_ag;
			
			if (Abs(lambtpp(i) / lambtpp(i - 1) - 1) < err_w) { 
				airforce[2] *= -rho*PI*radius*radius*vtipa*vtipa;//
				lambdh_ag = lambdi_ag * tppcoord.Ttransf[2][2] - vel[2] / vtipa;

				break; 
			}
		}
#ifdef OUTPUT_MODE
		lambtpp.outputs("_tail_lambtpp.output", 4);

#endif // OUTPUT_MODE
	}

#ifdef OUTPUT_MODE
	//if (!strcmp(type, "main")) {
	printf("%s rotor induced velocity iter count: %d \n", type, i);
	printf("lambdi_ag: %f, lambdt_ag: %f, lambdh_ag: %f \n", lambdi_ag, lambdt_ag, lambdh_ag);
	//}
#endif // OUTPUT_MODE
}


inline void Rotor::SetAirfm(void)
{
	//_setairfm_sp(airforce, airmoment);
	if (!strcmp(type, "main")) { AvrgInducedVel(); }
	else { 
		airforce[0] = airforce[1] = airmoment[0] = airmoment[1] = 0;
		AvrgInducedVel();
	}
}


void Rotor::functest(void) {
	Matrix1<myTYPE> uu(ns), sol;
	clock_t tStart;
	//_setbladevelc(uu, uu, 2, 3, 3);
	//cout << uu(0) << endl;

	//_flapmotionrt();
	//_setairfm_sp(airforce, airmoment);

	/*Matrix2<myTYPE> _incidn(nf, ns), _ua(nf, ns);
	Matrix1<int> id_ns(ns);
	Matrix1<myTYPE> _cl(ns), _cd(ns);
	id_ns = step(0, ns - 1);

	_incidn.input("incidn.txt");
	_ua.input("ma_n.txt");
	_aerodynacoef(_cl, _cd, _incidn(0, id_ns), _ua(0, id_ns));
	_cl.outputs("__cl.output", 4);
	_cd.outputs("__cd.output", 4);*/

	//SetAirfm();
	tStart = clock();
	for (int i = 0; i < 100; ++i) {
		//AvrgInducedVel();
		SetAirfm();
	}
	printf("Time taken: %fs\n", (double(clock() - tStart)) / CLOCKS_PER_SEC);

	printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
#ifdef TEST_MODE
	Arr2Mat(airforce, 3).outputs("airforce.output", 4);
	Arr2Mat(airmoment, 3).outputs("airmoment.output", 4);
#endif // TEST_MODE

	/*int a, b, c, nitermax;
	a = 1;
	b = 1;
	c = 1;*/
	//_setairfm(uu,uu,uu,a,b,c);
	//nitermax = 2;
	//sol.allocate(nitermax);
	//GenArf(sol, beta, a, 1.0, 1.0, 1.0, 0.0, 0.0, 0.01, nitermax, 0.5, 0.1);

	//_setbladevelc(uu, uu, a, b, c);
	//_aerodynacoef(uu, uu, uu, uu);

	cout << "test successfully." << endl;
}






