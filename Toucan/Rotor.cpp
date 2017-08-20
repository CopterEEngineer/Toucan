#include "stdafx.h"
#include "Component.h"


Rotor::Rotor(const char *s) :Component() {
	//print_cons_on_screen("Rotor");
	// read config file
	type = new char[10];
	strcpy(type, s);
	niter_a = -1;
	niter_w = -1;
	a0 = 5.7;
	if (!strcmp(s, "main")) {
		teeter = true;
		nf = 72;
		ns = 40;
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

		chord.allocate(ns);
		sweep.allocate(ns);
		twist.allocate(ns);
		azstation.allocate(nf, ns);
		rastation.allocate(nf, ns);
		chord.setvalue(0.558);
		sweep.setvalue(0);

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

		beta[0] = precone;
		beta[1] = 0 * PI / 180;
		beta[2] = 0 * PI / 180;

		myTYPE origin[3], euler[3];
		origin[0] = -0.0625;
		origin[1] = 0;
		origin[2] = -4.217;
		euler[0] = 0;
		euler[1] = 3.0 * PI / 180;
		euler[2] = 0;
		hubfxcoord.SetCoordinate(origin, euler, refcoord.base);
		refcoord.SetCoordinate(origin, euler, refcoord.base);

		// referred to hubfxcoord
		origin[0] = origin[1] = origin[2] = 0;
		euler[0] = 0;
		euler[1] = PI;
		euler[2] = 0;
		hubrtcoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubrtcoord
		origin[0] = eflap; //定义扭转，摆振，挥舞的三个桨叶坐标
		origin[1] = 0;
		origin[2] = 0;
		euler[0] = 0;
		euler[1] = precone + beta[1];
		euler[2] = 0;
		bladecoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubfxcoord
		origin[0] = 0;
		origin[1] = -precone * (1 - eflap);
		origin[2] = 0;
		euler[0] = euler[1] = euler[2] = 0;
		tppcoord.SetCoordinate(origin, euler, refcoord.base);


		// initialize member variables to zero
		mul = 0;
		lambdi_ag = sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI / radius / radius / vtipa / vtipa);
		lambdt_ag = lambdh_ag = lambdi_ag;
		power = 0;
		torque = 0;

		sita[0] = 12.1 * PI / 180;
		sita[1] = 1.18 * PI / 180;
		sita[2] = 3.12 * PI / 180;

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

		lambdi.setvalue(lambdi_ag);
		lambdh.setvalue(lambdh_ag);
		lambdt.setvalue(lambdt_ag);

		tipgeometry.allocate(nk, nf, 3);
		bladedeform.allocate(nf, ns, 3);


		ristation.allocate(nf, ni);

	}
	else if (!strcmp(s, "tail")) {
		teeter = false;

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
		gama = 2.2265;

		chord.allocate(ns);
		sweep.allocate(ns);
		twist.allocate(ns);
		azstation.allocate(nf, ns);
		rastation.allocate(nf, ns);
		chord.setvalue(0.167);

		beta[1] = beta[2] = 0;
		beta[0] = precone;

		myTYPE origin[3], euler[3];
		origin[0] = -13.48;
		origin[1] = 0;
		origin[2] = -1;
		euler[0] = PI / 2;
		euler[1] = 0;
		euler[2] = 0;
		hubfxcoord.SetCoordinate(origin, euler, refcoord.base);
		refcoord.SetCoordinate(origin, euler, refcoord.base);

		// referred to hubfxcoord
		origin[0] = origin[1] = origin[2] = 0;
		euler[0] = 0;
		euler[1] = PI;
		euler[2] = 0;
		hubrtcoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubrtcoord
		origin[0] = eflap; //定义扭转，摆振，挥舞的三个桨叶坐标
		origin[1] = 0;
		origin[2] = 0;
		euler[0] = 0;
		euler[1] = precone;
		euler[2] = 0;
		bladecoord.SetCoordinate(origin, euler, refcoord.base);
		// referred to hubfxcoord
		origin[0] = 0;
		origin[1] = -precone * (1 - eflap);
		origin[2] = 0;
		euler[0] = euler[1] = euler[2] = 0;
		tppcoord.SetCoordinate(origin, euler, refcoord.base);


		// initialize member variables to zero
		mul = 0;
		lambdi_ag = 0.01;
		//sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI/radius / radius / vtipa / vtipa);
		lambdt_ag = lambdh_ag = lambdi_ag;
		power = 0;
		torque = 0;

		sita[0] = 10.1 * PI / 180;
		sita[1] = sita[2] = 0.0;

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

		lambdi.setvalue(lambdi_ag);
		lambdh.setvalue(lambdh_ag);
		lambdt.setvalue(lambdt_ag);

		tipgeometry.allocate(nk, nf, 3);
		bladedeform.allocate(nf, ns, 3);

	}
	else {
		print_wrong_msg("Wrong rotor type. Please input main, or tail.");
	}


}


Rotor::Rotor(const Rotor &R) {
	//cout << "Rotor copy construtor." << endl;
	teeter = R.teeter;
	strcpy(type, R.type);
	niter_a = R.niter_a;
	niter_w = R.niter_w;
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
}


Rotor::~Rotor() {
	//cout << "Rotor destructor." << endl;
	teeter = false;
	delete[] type;
	niter_a = -1;
	niter_w = -1;
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

	sita[0] = sita[1] = sita[2] = 0;
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
	nitermax = 30 * 360 / dff;
	sol.allocate(nitermax);


	q0 = beta[0] - sita[2];
	dq0 = omega*sita[1];
	dt = dff / 180 * PI * omega;
	pho = 1;
	err_b = 0.001;
#ifdef TEST_MODE
	//err_b = 0;
#endif // TEST_MODE

	GenArf_rt(sol, beta, niter, m, c, k, q0, dq0, dff, nitermax, pho, err_b);

#ifdef TEST_MODE
	cout << (niter*dff / 360) << endl;
	printf("%f, %f\n\n", beta[1] * 180 / PI, beta[2] * 180 / PI);
#endif // TEST_MODE

#ifdef TEST_MODE
	int nperd = 360 / dff;
	sol.output("bflap_rt.output", 6);
	//sol.output(niter - nperd + 1, niter, 3);
#endif // TEST_MODE

}


void Rotor::_flapmotionfx(void) {
	myTYPE pho, dff, err_b;
	int nitermax, niter;
	Matrix2<myTYPE> sol;

	niter = 0;
	dff = 15;
	nitermax = 30 * 360 / dff;
	pho = 1;
	err_b = 0.001;

	if (teeter) {
		sol.allocate(2, nitermax);
		GenArf_fx_tr(sol, beta, niter, dff, nitermax, pho, err_b);
	}
	else {
		sol.allocate(3, nitermax);
		GenArf_fx_hg(sol, beta, niter, dff, nitermax, pho, err_b);
		//sol.output("sol.output", 6);
	}
}



void Rotor::_bladeCSD(void)
{
	;
}


void Rotor::SetCtrl(myTYPE * xctrl, const int n)
{
#ifdef _DEBUG
	if (n > 3) { printf("Undefined Control variations for %s rotor.", type); }
#endif // _DEBUG

	sita[0] = sita[1] = sita[2] = 0;
	for (int i = 0; i < n; ++i) {
		sita[i] = *(xctrl + i);
	}

}


void Rotor::GetCtrl(myTYPE * xctrl, const int n)
{
#ifdef _DEBUG
	if (n > 3) { printf("Undefined Control variations for %s rotor.", type); }
#endif // _DEBUG

	for (int i = 0; i < n; ++i) {
		xctrl[i] = sita[i];
	}

}


void Rotor::SetCoordBase(void)
{
	hubfxcoord.SetBase(refcoord.base);
	hubrtcoord.SetBase(&hubfxcoord);
	bladecoord.SetBase(&hubrtcoord);
	tppcoord.SetBase(&hubfxcoord);
#ifdef TEST_MODE
	//cout << hubfxcoord.base << endl;
	//cout << hubrtcoord.base << endl;
	//cout << bladecoord.base << endl;
	//cout << tppcoord.base << endl;
#endif // TEST_MODE

}


void Rotor::BladeDynamics(const char *t)
{
	myTYPE euler_temp[3];
	if (!strcmp(t, "main")) {
		_flapmotionrt();
		//_flapmotionfx();
	}
	else { _flapmotionfx(); }

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0.0;
	tppcoord.SetCoordinate(euler_temp, "euler");
}


void Rotor::AvrgInducedVel(void)
{
	const int itermax = 20;
	int i = 0;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp, twistt;

	if (!strcmp(type, "main")) {

		lambtpp[0] = sqrt(0.5*Copter::mass*UNIT_CONST / Copter::rho / PI / radius / radius / vtipa / vtipa);
		lambdi_ag = lambtpp[0];
		lambdh_ag = lambdi_ag;
		lambdh.setvalue(lambdh_ag);
		for (int i = 1; i < itermax; ++i) {

			BladeDynamics(type);

			veltpp[0] = veltpp[1] = veltpp[2] = 0;
			for (int i = 2; i >= 0; --i) {
				for (int j = 2; j >= 0; --j) {
					veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
				}
			}

			//#ifdef USE_DOUBLE
			//			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *tppcoord.Ttransf, 3, vel, 1, 0, veltpp, 1);
			//#else
			//			cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, *tppcoord.Ttransf, 3, vel, 1, 0, veltpp, 1);
			//#endif // USE_DOUBLE
			mul = vel[0] / vtipa;
			_setairfm_sp(airforce, airmoment);
			// airforce at hub coord approx to tpp coord
			lambdi_ag -= airforce[2] / (2 * rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp[i - 1] * lambtpp[i - 1]);
			lambdi_ag /= 2.0;
			// veltpp: tpp plane velocity with copter
			lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
			lambtpp[i] = lambdt_ag;
			// vel: hub velocity with copter
			lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
			lambdh.setvalue(lambdh_ag);

			iter = i;
			if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }

		}
		//lambtpp.outputs("_lambtpp.output", 4);
	}
	else {
		err_w = 5e-3;
		twistt = twist(ns - 1) - twist(0);

		lambtpp[0] = 0.01;
		lambdi_ag = lambtpp[0];
		lambdh_ag = lambdi_ag;
		lambdh.setvalue(lambdh_ag);
		for (int i = 1; i < itermax; ++i) {

			BladeDynamics(type);

			veltpp[0] = veltpp[1] = veltpp[2] = 0;
			for (int i = 2; i >= 0; --i) {
				for (int j = 2; j >= 0; --j) {
					veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
				}
			}

			mul = vel[0] / vtipa;

			ct = sigma*a0*0.5*(sita[0] / 3.0*(1 + 1.5*mul*mul) + 0.25*twistt*(1 + mul*mul) + 0.5*mul*beta[1] - 0.5*lambtpp[i - 1]);
			lambdi_ag += ct / 2.0 / sqrt(mul*mul + lambtpp[i - 1] * lambtpp[i - 1]);
			lambdi_ag /= 2.0;

			lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
			lambtpp[i] = lambdt_ag;

			lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
			lambdh.setvalue(lambdh_ag);

			airforce[2] = -ct*rho*PI*radius*radius*vtipa*vtipa;//

			iter = i;
			if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }
		}
#ifdef TEST_MODE
		//lambtpp.outputs("_tail_lambtpp.output", 4);

#endif // TEST_MODE
	}
	niter_w = iter;

}


inline void Rotor::SetAirfm(void)
{
	//_setairfm_sp(airforce, airmoment);
	if (!strcmp(type, "main")) { AvrgInducedVel(); }
	else {
		airforce[0] = airforce[1] = airmoment[0] = airmoment[1] = 0;
		AvrgInducedVel();
	}

#ifdef TEST_MODE

	cout << endl;
	printf("%s rotor induced velocity iter count: %d \n", type, niter_w);
	printf("%s rotor flaps (degs): %f, %f, %f\n", type, beta[0] / PI * 180, beta[1] / PI * 180, beta[2] / PI * 180);
	printf("lambdi_ag: %f, lambdt_ag: %f, lambdh_ag: %f \n", lambdi_ag, lambdt_ag, lambdh_ag);
	cout << endl;
	printf("%s rotor airdynamics: \n", type);
	printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;

#endif // TEST_MODE

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