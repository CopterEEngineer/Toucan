#include "stdafx.h"
#include "Components.h"


void Rotor::Allocate(void)
{
	_allocate();
}

void Rotor::_allocate(void)
{
	bflap.allocate(nf, ns); // flap
	dbflap.allocate(nf, ns);
	sfth.allocate(nf, ns);  // pitch 
	ut.allocate(nf, ns);    // air velocity
	un.allocate(nf, ns);
	up.allocate(nf, ns);
	ua.allocate(nf, ns);
	ma_n.allocate(nf, ns);
	incidn.allocate(nf, ns);// AOA
	inflow.allocate(nf, ns);
	cl.allocate(nf, ns);    // air coefficients
	cd.allocate(nf, ns);
	dt.allocate(nf, ns);
	_cl.allocate(ns);
	_cd.allocate(ns);
	_ua.allocate(ns);
	_incidn.allocate(ns);
	_inflow.allocate(ns);
	_factor.allocate(ns);
	
	lambdi.allocate(nf, ns);// induced velocity
	lambdx.allocate(nf, ns);
	lambdy.allocate(nf, ns);
	lambdh.allocate(nf, ns);
	lambdt.allocate(nf, ns);

	cirlb.allocate(nf, ns); // circulation
	tipstr.allocate(nk, nf);
	tipgeometry.allocate(nk, nf, 3);
	bladedeform.allocate(nf, ns, 3);

	lambdi.setvalue(lambdi_ag);
	lambdh.setvalue(lambdh_ag);
	lambdt.setvalue(lambdt_ag);

	ristation.allocate(nf, ni);
}

void Rotor::SetCoordBase(const Coordinate *base)
{
	refcoord.SetBase(base);
	hubfxcoord.SetBase(base);
	hubrtcoord.SetBase(&hubfxcoord);
	bladecoord.SetBase(&hubrtcoord);
	windcoord.SetBase(&hubfxcoord);
	tppcoord.SetBase(&hubfxcoord);
#ifdef TEST_MODE
	//cout << hubfxcoord.base << endl;
	//cout << hubrtcoord.base << endl;
	//cout << bladecoord.base << endl;
	//cout << tppcoord.base << endl;
#endif // TEST_MODE
}

void Rotor::SetAirfm(void)
{
	switch (type)
	{
	case Mrotor:
		AvrgInducedVel();
		break;
	case Trotor:
		airforce[0] = airforce[1] = airmoment[0] = airmoment[1] = 0;
		AvrgInducedVel();
		airforce[2] *= FT;
		break;
	default:
		break;
	}	
	
#ifdef OUTPUT_MODE_1

	cout << endl;
	printf("%d rotor induced velocity iter count: %d \n", type, niter_w);
	printf("%d rotor flaps (degs): %f, %f, %f\n", type, DEG(beta[0]), DEG(beta[1]), DEG(beta[2]));
	printf("lambdi_ag: %f, lambdt_ag: %f, lambdh_ag: %f \n", lambdi_ag, lambdt_ag, lambdh_ag);
	cout << endl;
	printf("%d rotor airdynamics: \n", type);
	printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;

#endif // OUTPUT_MODE_1
}

void Rotor::AvrgInducedVel(void)
{
	switch (type)
	{
	case Mrotor:
		_avrgInducedVel_MR();
		break;
	case Trotor:
		_avrgInducedVel_TR();
		break;
	default:
		break;
	}
}

void Rotor::WakeInducedVel(void)
{
	_wakeInducedVel();
}

void Rotor::_avrgInducedVel_MR(void)
{
	switch (bld.soltype)
	{
	case Rotation:
		_avrgInducedVel_Rt();
		break;
	case HubFixed:
		_avrgInducedVel_Fx();
	default:
		break;
	}
}

void Rotor::_avrgInducedVel_TR(void)
{
	switch (bld.soltype)
	{
	case Rotation:
		_avrgInducedVel_Rt();
		break;
	case HubFixed:
		_avrgInducedVel_TRFx();
	default:
		break;
	}
}

void Rotor::_avrgInducedVel_Rt(void)
{
	switch (hingetype)
	{
	case Teeter:
		_teeterdynamics_rt();
		break;
	case Hinged:
		_hingedynamics_rt();
		break;
	case Hingeless:
		_hingelessdynamics_rt();
		break;
	default:
		break;
	}
}

void Rotor::_avrgInducedVel_Fx(void)
{
	switch (hingetype)
	{
	case Teeter:
		_teeterdynamics_fx();
		break;
	case Hinged:
		_hingedynamics_fx();
		break;
	case Hingeless:
		_hingelessdynamics_fx();
		break;
	default:
		break;
	}
}

void Rotor::_avrgInducedVel_TRFx(void)
{
	const int itermax = 20;
	int iter = 0; 
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 5e-3;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp, twistt;

	twistt = twist(ns - 1) - twist(0);

	lambtpp[0] = 0.01;
	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);
	for (int i = 1; i < itermax; ++i) {

		//BladeDynamics(type);
		if (teeter)
			_teeterflap_fx();
		else
			_hingeflap_fx();

		veltpp[0] = veltpp[1] = veltpp[2] = 0;
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
			}
		}
		
		lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
		lambtpp[i] = lambdt_ag;
		//switch (adyna)
		//{
		//case Averaged:
		//	// blade element ct
		//	_setairfm_sp(airforce, airmoment);
		//	lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp[i - 1] * lambtpp[i - 1]);
		//	lambdi_ag /= 2.0;
		//	break;
		//case Simple:
		//	// simple ct
		//	ct = sigma*a0*0.5*(sita[0] / 3.0*(1 + 1.5*mul*mul) + 0.25*twistt*(1 + mul*mul) + 0.5*mul*beta[1] - 0.5*lambtpp[i - 1]);
		//	lambdi_ag += ct / 2.0 / sqrt(mul*mul + lambtpp[i - 1] * lambtpp[i - 1]);
		//	lambdi_ag /= 2.0;
		//	airforce[2] = -ct*amb.rho*PI*radius*radius*vtipa*vtipa;
		//	break;
		//default:
		//	// simple ct
		//	ct = sigma*a0*0.5*(sita[0] / 3.0*(1 + 1.5*mul*mul) + 0.25*twistt*(1 + mul*mul) + 0.5*mul*beta[1] - 0.5*lambtpp[i - 1]);
		//	lambdi_ag += ct / 2.0 / sqrt(mul*mul + lambtpp[i - 1] * lambtpp[i - 1]);
		//	lambdi_ag /= 2.0;
		//	airforce[2] = -ct*amb.rho*PI*radius*radius*vtipa*vtipa;
		//	break;
		//}
		//
		//lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
		//lambtpp[i] = lambdt_ag;

		//lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		//lambdh.setvalue(lambdh_ag);
		//lambdi.setvalue(lambdi_ag);

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) 
			break;
	}
	niter_w = iter;

	//printf("Flap in _avrgInducedVel_TRFx() counts: %d \n", niter_w);

}

void Rotor::_teeterdynamics_rt(void)
{
	const int itermax = 20;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp, twistt;
	myTYPE euler_temp[3] = { 0,0,0 };

	bld.Marf = 2 * iflap;
	bld.Carf = 0;
	bld.Ka = 2 * omega*omega*iflap + 2 * khub;
	bld._GenArfPrepare(RAD(bld.dff) / omega);

	if(si_unit)
		lambtpp[0] = sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);

	for (int i = 1; i < itermax; ++i) 
	{
		_teeterflap_rt();
		
		veltpp[0] = veltpp[1] = veltpp[2] = 0;
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
			}
		}
		
		lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
		lambtpp[i] = lambdt_ag;

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }
	}
	//lambtpp.outputs("_lambtpp.output", 4);
	niter_w = iter == itermax ? iter - 1 : iter;
	// reset wake flags
	if (haveGeo)
		haveGeo = false;
	if (haveStr)
		haveStr = false;
#ifdef OUTPUT_MODE_1
	printf("Flap in _teeterdynamics_rt() counts: %d \n", niter_w);
#endif // OUTPUT_MODE_1

}

void Rotor::_teeterdynamics_fx(void)
{
	const int itermax = 20;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp;

	if(si_unit)
		lambtpp[0] = sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);
	for (int i = 1; i < itermax; ++i)
	{
		_teeterflap_fx();

		veltpp[0] = veltpp[1] = veltpp[2] = 0;
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
			}
		}

		lambdt_ag = _aerodynamics(lambtpp[i-1], veltpp);
		lambtpp[i] = lambdt_ag;

		//_setairfm_sp(airforce, airmoment);
		//// airforce at hub coord approx to tpp coord
		//lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp[i - 1] * lambtpp[i - 1]);
		//lambdi_ag /= 2.0;
		//// veltpp: tpp plane velocity with copter
		//lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
		//lambtpp[i] = lambdt_ag;
		//// vel: hub velocity with copter
		//lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		//lambdh.setvalue(lambdh_ag);

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }
	}
	//lambtpp.outputs("_lambtpp.output", 4);
	niter_w = iter;
}

void Rotor::_teeterflap_rt(void)
{
	int niter, nitermax;
	myTYPE dt, af, az_temp;
	myTYPE qt, dq, q, t;
	Matrix1<myTYPE> dfx(ns), dfr(ns), dfz(ns), ra(ns);
	Matrix1<myTYPE> dfx2(ns), dfz2(ns), dfr2(ns);
	Matrix1<int> id_ns;
	myTYPE euler_temp[3];

	niter = 0, nitermax = bld.nitermax, af = bld.GAf.af;
	dt = RAD(bld.dff) / omega;
	id_ns = step(0, ns - 1);
	ra = rastation(0, id_ns);
	
	q = beta[0] - sita[2];
	dq = omega*sita[1];

	_setairfm(dfx, dfz, dfr, 0.0, q, dq, 0);
	_setairfm(dfx2, dfz2, dfr2, PI / omega, -q, -dq, nf / 2 - 1);

	qt = 4 * iflap * omega*(omg[0] * cos(0) - omg[1] * sin(0)) + 2 * iflap*(domg[0] * sin(0) + domg[1] * cos(0));
	qt += ((dfz - dfz2)*ra).sum()*radius;
	
	if (bld._GenArfStarter(qt, q, dq))
	{
		for (int i = 0; i < nitermax - 1; ++i)
		{
			t = i*dt + (1 - af)*dt;
			niter = i + 1;
			_setairfm(dfx, dfz, dfr, t, q, dq, i);
			_setairfm(dfx2, dfz2, dfr2, t + PI / omega, -q, -dq, i + nf / 2 - 1);

			qt = 4 * iflap * omega*(omg[0] * cos(omega*t) - omg[1] * sin(omega*t)) + 2 * iflap*(domg[0] * sin(omega*t) + domg[1] * cos(omega*t));
			qt += ((dfz - dfz2)*ra).sum()*radius;

			bld._GenArfTimeMarch(qt, dt, niter);
			az_temp = omega*(niter*dt + (1 - af)*dt);
			//while (az_temp>=2*PI)
			//	az_temp -= 2 * PI;
			bld.azmuth(niter) = az_temp;

			q = bld.q;
			dq = bld.dq;

			if (bld.isGenArfExit(niter))
				break;
		}
	}

	if (teeter)
	{
		int val1 = niter%bld.nperiod;
		int val2 = bld.nperiod / 4;
		beta[0] = precone;
		beta[1] = bld.sol(niter - val1);
		if (val1 < val2) beta[2] = -bld.sol(niter - val1 - val2);
		else beta[2] = bld.sol(niter - val1 + val2);
	}

	if (niter == nitermax - 1)
		;//printf("Warning: Flap solving may not be convergent in Func _teeterflap_rt(). \n");

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
	tppcoord.SetCoordinate(euler_temp, "euler");

	//bld.sol.output("blade flap.output", 8);
}

void Rotor::_teeterflap_fx(void)
{
	Matrix2<myTYPE> m22(2, 2), k22(2, 2), d22(2, 2), kk(2, 2);
	Matrix1<myTYPE> f22(2), dq0(2), q0(2), temp_M(2);
	Matrix1<myTYPE> ddq2(2), ddq(2), dq(2), q(2), qq(2); 
	int niter, nitermax = bld.nitermax;
	myTYPE dt, bt, twistt;
	myTYPE euler_temp[3], sita_temp[3];
	myTYPE k1, p2, sb, _sum, temp;
	Matrix2<myTYPE> sol(2, nitermax);
	
	niter = 0, dt = RAD(bld.dff) / omega;
	bld.GAf.InitCoef(dt);
	bt = 1.0 / bld.GAf.bt;

	m22(0, 0) = m22(1, 1) = 1.0;

	d22(0, 0) = omega*gama / 8.0;
	d22(0, 1) = 2 * omega;
	d22(1, 0) = -2 * omega;
	d22(1, 1) = omega*gama / 8.0;


	sita_temp[0] = sita[0];
	sita_temp[1] = sita[1];
	sita_temp[2] = sita[2];
	twistt = twist(ns - 1) - twist(0);

	k1 = tan(del*PI / 180);
	p2 = 1 + khub / iflap / omega / omega + eflap*m1*radius / iflap + gama*k1*(1 - 4.0 / 3 * eflap);
	k22(0, 0) = omega*omega*(p2 - 1 + gama*k1*mul*mul / 16.0);
	k22(0, 1) = omega*omega*gama / 8.0*(1 + 0.5*mul*mul);
	k22(1, 0) = -omega*omega*gama / 8.0*(1 - 0.5*mul*mul);
	k22(1, 1) = omega*omega*(p2 - 1 + 3.0 / 16.0*gama*k1*mul*mul);

	f22.v_p[0] = -omega*omega*(gama / 8.0*sita_temp[1] * (1 + 0.5*mul*mul));
	f22.v_p[1] = omega*omega*(-1.0 / 3.0*gama*mul*sita_temp[0] - 1.0 / 4.0*gama*mul*(twistt)-gama / 8.0*sita_temp[2] * (1 + 1.5*mul*mul));

	f22.v_p[0] += omega*omega*(-2.0 / omega*omg[0] - gama / 8.0 / omega*omg[1] - 1.0 / omega / omega*domg[1]);
	f22.v_p[1] += omega*omega*(-gama / 8.0 / omega*omg[0] + 2.0 / omega*omg[1] - 1.0 / omega / omega*domg[0]);

	f22.v_p[1] += omega*omega*(-0.25*gama*mul)*lambdh_ag;

	_sum = 0;
	temp = 0;

	sb = 8.0*(p2 - 1.0) / gama;
	q0(0) = -1.0 / (1.0 + sb*sb)*(sb*sita_temp[1] - sita_temp[2] + (sb*16.0 / gama - 1)*omg[0] + (sb + 16.0 / gama)*omg[1]);
	q0(1) = -1.0 / (1.0 + sb*sb)*(sb*sita_temp[2] + sita_temp[1] + (sb + 16.0 / gama)*omg[0] - (sb*16.0 / gama - 1)*omg[1]);
	dq0.v_p[0] = 0.0;
	dq0.v_p[1] = omega*sol(1, 0);


	temp_M = d22.matrixmultiplyP2(dq0) + k22.matrixmultiplyP2(q0) - f22;
	Msolver(m22.v_p, temp_M.v_p, 2);
	ddq = temp_M*(-1);
	dq = dq0;
	q(0) = q0(0);
	q(1) = q0(1);

	sol(0, 0) = q(0);
	sol(1, 0) = q(1);

	for (int i = 0; i < nitermax - 1; ++i) {
		qq = f22 - (k22.matrixmultiplyP2(q)) * bld.GAf.af;
		qq += m22.matrixmultiplyP2(q*bld.GAf.c0 + dq*bld.GAf.c2 + ddq*bld.GAf.c3);
		qq += d22.matrixmultiplyP2(q*bld.GAf.c1 + dq*bld.GAf.c4 + ddq*bld.GAf.c5);
		kk = k22*bld.GAf.ck + m22*bld.GAf.c0 + d22*bld.GAf.c1;
		Msolver(kk.v_p, qq.v_p, 2);
		q = qq;
		for (int j = 0; j < 2; ++j) {
			ddq2(j) = (q(j) - sol(j, i)) * bld.GAf.bt / dt / dt - dq(j) * bld.GAf.bt / dt - ddq(j)*(1.0 * bld.GAf.bt / 2.0 - 1.0);
		}
		dq += ddq*(1.0 - bld.GAf.r)*dt + ddq2*bld.GAf.r*dt;
		ddq = ddq2;
		sol(0, i + 1) = q(0);
		sol(1, i + 1) = q(1);

		// exit condition
		_sum = 0;
		for (int j = 0; j < 2; ++j) {
			//temp = sol(j, i + 1) / sol(j, i) - 1.0;
			temp = (sol(j, i + 1) - sol(j, i)) / (precone + 0.01*RAD(1));// in case of 0/0
			temp *= temp;
			_sum += temp;
		}

		niter = i + 1;
		beta[0] = precone;
		beta[1] = -sol(0, niter);
		beta[2] = -sol(1, niter);
		if (_sum < bld.err_b*bld.err_b) { break; }
	}

	if (niter == nitermax - 1)
		;// printf("Warning: Flap solving may not be convergent in Func _teeterflap_rt(). \n");

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
	tppcoord.SetCoordinate(euler_temp, "euler");
}

void Rotor::_hingedynamics_rt(void)
{
	const int itermax = 20;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 5e-3;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp, twistt;
	myTYPE euler_temp[3] = { 0,0,0 };

	bld.Marf = iflap;
	bld.Carf = 0;
	bld.Ka = omega*omega*(iflap + eflap*m1*radius) + khub;

	bld._GenArfPrepare(RAD(bld.dff) / omega);

	lambtpp[0] = 0.01;
	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);

	for (int i = 1; i < itermax; ++i)
	{
		_hingeflap_rt();

		veltpp[0] = veltpp[1] = veltpp[2] = 0;
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
			}
		}

		lambdt_ag = _aerodynamics(lambtpp[i-1], veltpp);
		lambtpp[i] = lambdt_ag;

		//_setairfm_sp(airforce, airmoment);
		//// airforce at hub coord approx to tpp coord
		//lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp[i - 1] * lambtpp[i - 1]);
		//lambdi_ag /= 2.0;
		//// veltpp: tpp plane velocity with copter
		//lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
		//lambtpp[i] = lambdt_ag;
		//// vel: hub velocity with copter
		//if (haveGeo && haveStr)
		//{
		//	_bladePosition();
		//	//_wakeIndVelCalc();
		//	_wakeInducedVelMP(nb);
		//}
		//else
		//{
		//	lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		//	lambdh.setvalue(lambdh_ag);
		//	lambdi.setvalue(lambdi_ag);
		//}

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }
	}
	//lambtpp.outputs("_lambtpp.output", 4);
	niter_w = iter == itermax ? iter - 1 : iter;
	// reset wake flags
	if (haveGeo)
		haveGeo = false;
	if (haveStr)
		haveStr = false;
#ifdef OUTPUT_MODE_1
	printf("Flap in _teeterdynamics_rt() counts: %d \n", niter_w);
#endif // OUTPUT_MODE_1
}

void Rotor::_hingedynamics_fx(void)
{
	const int itermax = 20;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp;

	if(si_unit)
		lambtpp[0] = sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);
	for (int i = 1; i < itermax; ++i)
	{
		_hingeflap_fx();

		veltpp[0] = veltpp[1] = veltpp[2] = 0;
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				veltpp[i] += tppcoord.Ttransf[i][j] * vel[j];
			}
		}

		lambdt_ag = _aerodynamics(lambtpp[i-1], veltpp);
		lambtpp[i] = lambdt_ag;

		//_setairfm_sp(airforce, airmoment);
		//// airforce at hub coord approx to tpp coord
		//lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp[i - 1] * lambtpp[i - 1]);
		//lambdi_ag /= 2.0;
		//// veltpp: tpp plane velocity with copter
		//lambdt_ag = lambdi_ag - veltpp[2] / vtipa;
		//lambtpp[i] = lambdt_ag;
		//// vel: hub velocity with copter
		//lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		//lambdh.setvalue(lambdh_ag);

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }
	}
	//lambtpp.outputs("_lambtpp.output", 4);
	niter_w = iter;
}

void Rotor::_hingeflap_rt(void)
{
	int niter, nitermax;
	double mb;
	myTYPE dt, af;
	myTYPE qt, dq, q, t;
	Matrix1<myTYPE> dfx(ns), dfr(ns), dfz(ns), ra(ns);
	Matrix1<int> id_ns;
	myTYPE euler_temp[3];
	Matrix1<myTYPE> azmuth(bld.nitermax);

	niter = 0, nitermax = bld.nitermax, af = bld.GAf.af;
	dt = RAD(bld.dff) / omega;
	id_ns = step(0, ns - 1);
	ra = rastation(0, id_ns);

	mb = m1*radius;

	q = beta[0] - sita[2];
	dq = omega*sita[1];

	_setairfm(dfx, dfz, dfr, 0.0, q, dq, 0);

	qt = 2 * (iflap + eflap*mb)* omega*(omg[0] * cos(0) - omg[1] * sin(0)) + iflap*(domg[0] * sin(0) + domg[1] * cos(0));
	//qt += m1*(dvel[2] - vel[0] * omg[1] + vel[1] * omg[0]);
	qt += (dfz*ra).sum()*radius;

	if (bld._GenArfStarter(qt, q, dq))
	{
		for (int i = 0; i < nitermax - 1; ++i)
		{
			t = i*dt + (1 - af)*dt;
			azmuth(i) = omega*t;
			niter = i + 1;
			_setairfm(dfx, dfz, dfr, t, q, dq, i);
			qt = 2 * (iflap + eflap*mb)* omega*(omg[0] * cos(omega*t) - omg[1] * sin(omega*t)) + iflap*(domg[0] * sin(omega*t) + domg[1] * cos(omega*t));
			//qt += m1*(dvel[2] - vel[0] * omg[1] + vel[1] * omg[0]);
			qt += (dfz*ra).sum()*radius;

			bld._GenArfTimeMarch(qt, dt, niter);
			q = bld.q;
			dq = bld.dq;

			if (bld.isGenArfExit(niter))
				break;
		}
	}

	beta[0] = beta[1] = beta[2] = 0;
	int val1 = niter%bld.nperiod;
	int val2 = bld.nperiod / 4;
	for (int i = niter - bld.nperiod + 1; i <= niter; ++i)
	{
		switch (nb)
		{
		case 4:
			beta[1] += (bld.sol(i)*cos(azmuth(i)) + bld.sol(i - bld.nperiod / 4)*sin(azmuth(i)) - bld.sol(i - bld.nperiod / 2)*cos(azmuth(i)) - bld.sol(i - bld.nperiod / 4 * 3)*sin(azmuth(i))) / 2;
			beta[2] += (bld.sol(i)*sin(azmuth(i)) - bld.sol(i - bld.nperiod / 4)*cos(azmuth(i)) - bld.sol(i - bld.nperiod / 2)*sin(azmuth(i)) + bld.sol(i - bld.nperiod / 4 * 3)*cos(azmuth(i))) / 2;
			beta[0] += (bld.sol(i) + bld.sol(i - bld.nperiod / 4) + bld.sol(i - bld.nperiod / 2) + bld.sol(i - bld.nperiod / 4 * 3)) / 4;
			break;
		case 2:
			beta[0] += precone;
			beta[1] += bld.sol(niter - val1);
			if (val1 < val2) beta[2] += -bld.sol(niter - val1 - val2);
			else beta[2] += bld.sol(niter - val1 + val2);
			break;

		default:
			printf("Undefine nb = %d in _hingeflap_rt() \n", nb);
			system("pause");
			break;
		}
	}
	beta[0] /= bld.nperiod;
	beta[1] /= bld.nperiod;
	beta[2] /= bld.nperiod;
	//bld.sol.output("tr_bld_sol.output", 10);
	//azmuth.output("tr_bld_azmuth.output", 10);
	if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _hingeflap_rt(). Count: %d \n", niter);

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
	tppcoord.SetCoordinate(euler_temp, "euler");
}

void Rotor::_hingeflap_fx(void)
{
	Matrix2<myTYPE> m22(3, 3), k22(3, 3), d22(3, 3), kk(3, 3);
	Matrix1<myTYPE> f22(3), dq0(3), q0(3), temp_M(3);
	Matrix1<myTYPE> ddq2(3), ddq(3), dq(3), q(3), qq(3);
	Matrix2<myTYPE> f_temp(3, 4);
	Matrix1<myTYPE> c_temp(4);
	int niter, nitermax = bld.nitermax;
	myTYPE twistt, dt;
	myTYPE euler_temp[3], sita_temp[3];
	myTYPE k1, p2, sb, mb, _sum, temp;
	Matrix2<myTYPE> sol(3, nitermax);

	niter = 0, dt = RAD(bld.dff) / omega;
	bld.GAf.InitCoef(dt);

	m22(0, 0) = m22(1, 1) = m22(2, 2) = 1.0;

	d22(0, 0) = omega*gama*(0.125 - 1.0 / 3.0*eflap + 0.25*eflap*eflap);
	d22(0, 1) = 0;
	d22(0, 2) = -0.25*omega*gama*mul*(1.0 / 3.0 - eflap + eflap*eflap);
	d22(1, 0) = 0;
	d22(1, 1) = omega*gama*(0.125 - 1.0 / 3.0*eflap + 0.25*eflap*eflap);
	d22(1, 2) = 2 * omega;
	d22(2, 0) = -0.5*omega*mul*(1.0 / 3.0 - eflap + eflap*eflap);
	d22(2, 1) = -2 * omega;
	d22(2, 2) = omega*gama*(0.125 - 1.0 / 3.0*eflap + 0.25*eflap*eflap);


	sita_temp[0] = sita[0];
	sita_temp[1] = sita[1];
	sita_temp[2] = sita[2];
	twistt = twist(ns - 1) - twist(0);

	k1 = RAD(del);
	mb = m1*radius;
	p2 = 1 + khub / iflap / omega / omega + eflap*mb / iflap + gama*k1*(1 - 4.0 / 3 * eflap);

	k22(0, 0) = omega*omega*(p2 + 0.25*gama*k1*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	k22(0, 1) = omega*omega*(-0.25*gama*mul)*eflap*(0.5 - eflap);
	k22(0, 2) = omega*omega*(-0.25*gama*k1*mul)*(2.0 / 3.0 - eflap);
	k22(1, 0) = omega*omega*(-0.5*gama*mul)*(1.0 / 3.0 - 0.5*eflap);
	k22(1, 1) = omega*omega*(p2 - 1 + 0.125*gama*mul*mul*k1*(0.5 - eflap + 0.5*eflap*eflap));
	k22(1, 2) = omega*omega*((0.5*gama)*(0.25 - 2.0 / 3.0*eflap + 0.5*eflap*eflap) + (0.125*gama*mul*mul)*(0.5 - eflap + 0.5*eflap*eflap));
	k22(2, 0) = omega*omega*(-0.5*gama*k1*mul)*(2.0 / 3.0 - eflap);
	k22(2, 1) = omega*omega*(-0.5*gama*(0.25 - 2.0 / 3.0*eflap + 0.5*eflap*eflap) + 0.125*gama*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	k22(2, 2) = omega*omega*(p2 - 1 + 3.0 / 8.0*gama*k1*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));

	c_temp(0) = sita_temp[0];
	c_temp(1) = twistt;
	c_temp(2) = -sita_temp[1];
	c_temp(3) = -sita_temp[2];

	f_temp(0, 0) = 0.5*gama*(0.25 - 1.0 / 3.0*eflap + 0.5*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	f_temp(0, 1) = 0.5*gama*(0.2 - 0.25*eflap + 0.5*mul*mul*(1.0 / 3.0 - 0.5*eflap));
	f_temp(0, 2) = 0;
	f_temp(0, 3) = -0.5*gama*mul*(1.0 / 3.0 - 0.5*eflap);
	f_temp(1, 0) = 0;
	f_temp(1, 1) = 0;
	f_temp(1, 3) = 0;
	f_temp(1, 2) = 0.5*gama*(0.25 - 1.0 / 3.0*eflap + 0.25*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	f_temp(2, 0) = -0.5*gama*mul*(2.0 / 3.0 - eflap);
	f_temp(2, 1) = -0.5*gama*mul*(0.5 - 2.0 / 3.0*eflap);
	f_temp(2, 2) = 0;
	f_temp(2, 3) = 0.5*gama*(0.25 - 1.0 / 3.0*eflap + 0.75*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));


	f22 = (f_temp.matrixmultiplyP2(c_temp)) * omega*omega;

	c_temp(0) = omg[0];
	c_temp(1) = omg[1];
	c_temp(2) = domg[0];
	c_temp(3) = domg[1];

	f_temp(0, 0) = 0.125*gama*mul / omega*(2.0 / 3.0 - eflap);
	f_temp(0, 1) = f_temp(0, 2) = f_temp(0, 3) = 0;
	f_temp(1, 0) = -2.0 / omega*(1.0 + eflap*mb / iflap);
	f_temp(1, 1) = -0.5*gama / omega*(0.25 - 1.0 / 3.0*eflap);
	f_temp(1, 2) = 0;
	f_temp(1, 3) = -1.0 / omega / omega;
	f_temp(2, 0) = -0.5*gama / omega*(0.25 - 1.0 / 3.0*eflap);
	f_temp(2, 1) = 2.0 / omega*(1.0 + eflap*mb / iflap);
	f_temp(2, 2) = 1.0 / omega / omega;
	f_temp(2, 3) = 0;

	f22 += (f_temp.matrixmultiplyP2(c_temp)) * omega*omega;
	f22.v_p[0] += omega*omega*0.5*gama*(1.0 / 3.0 - 0.5*eflap)*lambdh_ag;
	f22.v_p[2] += omega*omega*(-0.5*gama*mul)*(0.5 - eflap + 0.5*eflap*eflap)*lambdh_ag;

	if (!si_unit)
		f22.v_p[0] -= mb / iflap *9.8*0.3048*0.3048;
	else
		f22.v_p[0] -= mb / iflap*9.8;


	_sum = 0;
	temp = 0;
	sb = 8.0*(p2 - 1.0) / gama;
	q0(1) = -1.0 / (1.0 + sb*sb)*(sb*sita_temp[1] - sita_temp[2] + (sb*16.0 / gama - 1)*omg[0] + (sb + 16.0 / gama)*omg[1]);
	q0(2) = -1.0 / (1.0 + sb*sb)*(sb*sita_temp[2] + sita_temp[1] + (sb + 16.0 / gama)*omg[0] - (sb*16.0 / gama - 1)*omg[1]);
	q0(0) = gama / p2*(sita_temp[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sita_temp[2]) - lambdh_ag / 6.0);

	dq0.v_p[0] = 0.0;
	dq0.v_p[1] = 0;
	dq0.v_p[2] = omega*q0(2);

	temp_M = d22.matrixmultiplyP2(dq0) + k22.matrixmultiplyP2(q0) - f22;
	Msolver(m22.v_p, temp_M.v_p, 3);
	ddq = temp_M*(-1);
	dq = dq0;
	q(0) = q0(0);
	q(1) = q0(1);
	q(2) = q0(2);

	sol(0, 0) = q(0);
	sol(1, 0) = q(1);
	sol(2, 0) = q(2);

	for (int i = 0; i < nitermax - 1; ++i) {
		qq = f22 - (k22.matrixmultiplyP2(q)) * bld.GAf.af;
		qq += m22.matrixmultiplyP2(q*bld.GAf.c0 + dq*bld.GAf.c2 + ddq*bld.GAf.c3);
		qq += d22.matrixmultiplyP2(q*bld.GAf.c1 + dq*bld.GAf.c4 + ddq*bld.GAf.c5);
		kk = k22*bld.GAf.ck + m22*bld.GAf.c0 + d22*bld.GAf.c1;
		Msolver(kk.v_p, qq.v_p, 3);
		q = qq;
		for (int j = 0; j < 3; ++j) {
			ddq2(j) = (q(j) - sol(j, i)) * bld.GAf.bt / dt / dt - dq(j) * bld.GAf.bt / dt - ddq(j)*(1.0 * bld.GAf.bt / 2.0 - 1.0);
		}
		dq += ddq*(1.0 - bld.GAf.r)*dt + ddq2*bld.GAf.r*dt;
		ddq = ddq2;
		sol(0, i + 1) = q(0);
		sol(1, i + 1) = q(1);
		sol(2, i + 1) = q(2);

		niter = i + 1;
		beta[0] = sol(0, niter);
		beta[1] = -sol(1, niter);
		beta[2] = -sol(2, niter);

		// exit condition
		if (i > 10)
		{
			_sum = 0;
			for (int j = 0; j < 3; ++j) 
			{
				//temp = sol(j, i + 1) / sol(j, i) - 1.0;
				temp = (sol(j, i + 1) - sol(j, i)) / (sol(j, i) + RAD(0.01));
				//temp = (sol(j, i + 1) - sol(j, i)) / sol(0, i);// in case of 0/0
				temp *= temp;
				_sum += temp;
			}
			if (_sum < bld.err_b*bld.err_b) 
				break; 
		}

	}
	//sol.output("_sol_temp_tr.output", 10);
	if (niter == nitermax - 1)
		;// printf("Warning: Flap solving may not be convergent in Func _hingeflap_fx(). \n");

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
	tppcoord.SetCoordinate(euler_temp, "euler");
}

void Rotor::_hingelessdynamics_rt(void)
{
	const int itermax = 20;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE dveltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp, twistt;
	myTYPE euler_temp[3] = { 0,0,0 };

	bld.Marf = iflap;
	bld.Carf = 0;
	bld.Ka = omega*omega*iflap+khub;

	bld._GenArfPrepare(RAD(bld.dff) / omega);

	if (si_unit)
		lambtpp[0] = sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);

	for (int i = 1; i < itermax; i++)
	{
		_hingelessflap_rt();

		_velTransform(veltpp, dveltpp, tppcoord);

		lambdt_ag = _aerodynamics(lambtpp[i-1], veltpp);
		lambtpp[i] = lambdt_ag;

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) { break; }
	}
	//lambdh.output("_lambdh_hingeless.output", 10);
	//lambdi.output("_lambdi_hingeless.output", 10);
	niter_w = iter == itermax ? iter - 1 : iter;
	// reset wake flags
	if (haveGeo)
		haveGeo = false;
	if (haveStr)
		haveStr = false;
#ifdef OUTPUT_MODE_1
	printf("Flap in _hingelessdynamics_rt() counts: %d \n", niter_w);
#endif // OUTPUT_MODE_1
}

void Rotor::_hingelessdynamics_fx(void)
{
	const int itermax = 20;
	int iter = 0;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE dveltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp;
	double lamb_1s, lamb_1c;

	if (si_unit)
		lambtpp[0] = sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);

	lamb_1s = lamb_1c = 0.0;
	for (int i = 1; i < itermax; ++i)
	{
		//_hingelessflap_fx();
		_flapMBC(-lamb_1s*0.0, -lamb_1c*0.0);
		//_flapWang();

		_velTransform(veltpp, dveltpp, tppcoord);
		lamb_1s = -2 * veltpp[0] / vtipa*lambdi_ag;
		lamb_1c = tan(0.5*Atan2(veltpp[0] / vtipa, -veltpp[2] / vtipa + lambdi_ag))*lambdi_ag;

		lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
		lambtpp[i] = lambdt_ag;

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) 
			break; 
	}
	//lambtpp.outputs("_lambtpp.output", 4);
	niter_w = iter;
}

void Rotor::_hingelessflap_rt(void)
{
	int niter, nitermax;
	myTYPE dt, af, mb;
	myTYPE qt, dq, q, t;
	Matrix1<myTYPE> dfx(ns), dfr(ns), dfz(ns), ra(ns);
	Matrix1<int> id_ns;
	myTYPE euler_temp[3];
	Matrix1<myTYPE> azmuth(bld.nitermax);

	niter = 0, nitermax = bld.nitermax, af = bld.GAf.af;
	dt = RAD(bld.dff) / omega;
	id_ns = step(0, ns - 1);
	ra = rastation(0, id_ns);

	mb = m1*radius;

	q = beta[0] - sita[2];
	dq = omega*sita[1];

	_setairfm(dfx, dfz, dfr, 0.0, q, dq, 0);

	qt = 2 * iflap* omega*(omg[0] * cos(0) - omg[1] * sin(0)) + iflap*(domg[0] * sin(0) + domg[1] * cos(0));
	//if (si_unit)
	//	qt += mb*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8);
	//else
	//	qt += mb*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8*0.3048*0.3048);
	qt += (dfz*ra).sum()*radius;

	if (bld._GenArfStarter(qt, q, dq))
	{
		for (int i = 0; i < nitermax - 1; ++i)
		{
			t = i*dt + (1 - af)*dt;
			azmuth(i) = omega*t;
			niter = i + 1;
			_setairfm(dfx, dfz, dfr, t, q, dq, i);
			qt = 2 * iflap * omega*(omg[0] * cos(omega*t) - omg[1] * sin(omega*t)) + iflap*(domg[0] * sin(omega*t) + domg[1] * cos(omega*t));
			//if(si_unit)
			//	qt += mb*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8);
			//else
			//	qt += mb*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8 * 0.3048 * 0.3048);

			qt += (dfz*ra).sum()*radius;

			bld._GenArfTimeMarch(qt, dt, niter);
			q = bld.q;
			dq = bld.dq;

			if (bld.isGenArfExit(niter))
				break;
		}
	}

	beta[0] = beta[1] = beta[2] = 0;
	int val1 = niter%bld.nperiod;
	int val2 = bld.nperiod / 4;
	for (int i = niter - bld.nperiod + 1; i <= niter; ++i)
	{
		switch (nb)
		{
		case 4:
			beta[1] += (bld.sol(i)*cos(azmuth(i)) + bld.sol(i - bld.nperiod / 4)*sin(azmuth(i))-bld.sol(i - bld.nperiod / 2)*cos(azmuth(i)) - bld.sol(i - bld.nperiod / 4 * 3)*sin(azmuth(i)))/2;
			beta[2] += (bld.sol(i)*sin(azmuth(i)) - bld.sol(i - bld.nperiod / 4)*cos(azmuth(i))-bld.sol(i - bld.nperiod / 2)*sin(azmuth(i)) + bld.sol(i - bld.nperiod / 4 * 3)*cos(azmuth(i)))/2;
			beta[0] += (bld.sol(i)+bld.sol(i - bld.nperiod / 4) + bld.sol(i - bld.nperiod / 2) + bld.sol(i - bld.nperiod / 4 * 3)) / 4;
			break;
		case 2:
			beta[0] += precone;
			beta[1] += bld.sol(niter - val1);
			if (val1 < val2) beta[2] += -bld.sol(niter - val1 - val2);
			else beta[2] += bld.sol(niter - val1 + val2);
			break;

		default:
			printf("Undefine nb = %d in _hingelessflap_rt() \n", nb);
			system("pause");
			break;
		}
	}
	beta[0] /= bld.nperiod;
	beta[1] /= bld.nperiod;
	beta[2] /= bld.nperiod;

	if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _teeterflap_rt(). Count: %d \n", niter);
	//bld.sol.output("bld_mr.output", 10);
	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
	tppcoord.SetCoordinate(euler_temp, "euler");
}

void Rotor::_flapWang(void)
{
	Matrix2<myTYPE> m22(3, 3), k22(3, 3), d22(3, 3), kk(3, 3);
	Matrix1<myTYPE> f22(3), dq0(3), q0(3), temp_M(3);
	Matrix1<myTYPE> ddq2(3), ddq(3), dq(3), q(3), qq(3);
	int niter, nitermax = bld.nitermax;
	myTYPE twistt, dt;
	myTYPE euler_temp[3];
	myTYPE sitaw[3], betaw[3];
	myTYPE k1, p2, sb, mb, _sum, temp;
	Matrix2<myTYPE> sol(3, nitermax);

	niter = 0, dt = RAD(bld.dff) / omega;
	bld.GAf.InitCoef(dt);

	k1 = tan(RAD(del));
	mb = m1*radius;
	p2 = 1 + khub / iflap / omega / omega + eflap*mb / iflap + gama*k1*(1 - 4.0 / 3 * eflap);

	sitaw[0] = sita[0];
	sitaw[1] = sita[1] * cos(windcoord.euler[2]) - sita[2] * sin(windcoord.euler[2]);
	sitaw[2] = sita[1] * sin(windcoord.euler[2]) + sita[2] * cos(windcoord.euler[2]);
	twistt = twist(ns - 1) - twist(0);

	m22(0, 0) = m22(1, 1) = m22(2, 2) = 1.0;

	d22(0, 0) = omega*gama/2.0*(1.0/4.0 - 2.0 / 3.0*eflap + eflap*eflap / 2.0);
	d22(0, 1) = 0;
	d22(0, 2) = omega*gama*mul/4.0*(1.0 / 3.0 - eflap + eflap*eflap);
	d22(1, 0) = 0;
	d22(1, 1) = omega*gama/2.0*(1.0/4.0 - 2.0 / 3.0*eflap + eflap*eflap/2.0);
	d22(1, 2) = 2 * omega;
	d22(2, 0) = omega*gama*mul/2.0*(1.0 / 3.0 - eflap + eflap*eflap);
	d22(2, 1) = -2 * omega;
	d22(2, 2) = omega*gama/2.0*(1.0/4.0 - 2.0 / 3.0*eflap + eflap*eflap/2.0);

	k22(0, 0) = omega*omega*(p2 + gama*k1*mul*mul/4.0*(1.0/2.0 - eflap + eflap*eflap/2.0));
	k22(0, 1) = omega*omega*gama*mul/4.0*eflap*(1.0/2.0 - eflap);
	k22(0, 2) = omega*omega*gama*k1*mul/4.0*(2.0 / 3.0 - eflap);
	k22(1, 0) = omega*omega*gama*mul/2.0*(1.0 / 3.0 - eflap/2.0);
	k22(1, 1) = omega*omega*(p2 - 1 + gama*mul*mul*k1/8.0*(1.0/2.0 - eflap + 0.5*eflap*eflap));
	k22(1, 2) = omega*omega*(0.5*gama*(0.25 - 2.0 / 3.0*eflap + 0.5*eflap*eflap) + 0.125*gama*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	k22(2, 0) = omega*omega*0.5*gama*k1*mul*(2.0 / 3.0 - eflap);
	k22(2, 1) = omega*omega*(-0.5*gama*(0.25 - 2.0 / 3.0*eflap + 0.5*eflap*eflap) + 0.125*gama*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	k22(2, 2) = omega*omega*(p2 - 1 + 3.0 / 8.0*gama*k1*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));

	f22(0) = gama / 2.0 * ((0.25 - eflap / 3.0) + mul*mul / 2.0 * (0.5 - eflap + eflap*eflap*0.5))*sitaw[0];
	f22(0) += gama / 2.0 * (mul*(1.0 / 3.0 - 0.5*eflap))*sitaw[2];
	f22(0) += gama / 2.0*((0.2 - 0.25*eflap) + mul*mul*0.5*(1.0 / 3.0 - 0.5*eflap))*twistt;
	f22(0) += gama / 2.0*(1.0 / 3.0 - 0.5*eflap)*lambdh_ag + gama / 8.0*mul*(2.0 / 3.0 - eflap)*omgw[0] / omega;
	//if (!si_unit)
	//	f22(0) += mb / iflap / omega / omega *(-9.8*0.3048*0.3048 + dvel[2] - vel[0] * omg[1] + vel[1] * omg[0]);
	//else
	//	f22(0) += mb / iflap / omega / omega *(-9.8 + dvel[2] - vel[0] * omg[1] + vel[1] * omg[0]);

	f22(1) = 2 * (1 + eflap*mb / iflap)*omgw[0] / omega + (1 + eflap*mb / iflap)*domgw[1] / omega / omega;
	f22(1) += gama / 2.0*((0.25 - eflap / 3.0) + mul*mul / 4.0*(0.5 - eflap + 0.5*eflap*eflap))*sitaw[1] + gama / 2.0*(0.25 - eflap / 3.0)*omgw[1] / omega;

	f22(2) = -2.0*(1 + eflap*mb / iflap)*omgw[1] / omega + (1 + eflap*mb / iflap)*domgw[0] / omega / omega;
	f22(2) += gama*mul / 2.0*(2.0 / 3.0 - eflap)*sitaw[0] + gama*mul / 2.0*(0.5 - 2.0 / 3.0*eflap)*twistt;
	f22(2) += gama / 2.0*((0.25 - eflap / 3.0) + 3.0 / 4.0*mul*mul*(0.5 - eflap + 0.5*eflap*eflap))*sitaw[2];
	f22(2) += gama*mul / 2.0*(0.5 - eflap + eflap*eflap / 2.0)*lambdh_ag + gama / 2.0*(0.25 - eflap / 3.0)*omgw[0] / omega;

	f22 = f22 * (omega*omega);

	_sum = 0;
	temp = 0;
	sb = 8.0*(p2 - 1.0) / gama;
	q0(1) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	q0(2) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb + 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);
	q0(0) = gama / p2*(sitaw[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sitaw[2]) - lambdh_ag / 6.0);

	dq0.v_p[0] = 0.0;
	dq0.v_p[1] = 0;
	dq0.v_p[2] = omega*q0(2);

	temp_M = d22.matrixmultiplyP2(dq0) + k22.matrixmultiplyP2(q0) - f22;
	Msolver(m22.v_p, temp_M.v_p, 3);
	ddq = temp_M*(-1);
	dq = dq0;
	q(0) = q0(0);
	q(1) = q0(1);
	q(2) = q0(2);

	sol(0, 0) = q(0);
	sol(1, 0) = q(1);
	sol(2, 0) = q(2);

	for (int i = 0; i < nitermax - 1; ++i) {
		qq = f22 - (k22.matrixmultiplyP2(q)) * bld.GAf.af;
		qq += m22.matrixmultiplyP2(q*bld.GAf.c0 + dq*bld.GAf.c2 + ddq*bld.GAf.c3);
		qq += d22.matrixmultiplyP2(q*bld.GAf.c1 + dq*bld.GAf.c4 + ddq*bld.GAf.c5);
		kk = k22*bld.GAf.ck + m22*bld.GAf.c0 + d22*bld.GAf.c1;
		Msolver(kk.v_p, qq.v_p, 3);
		q = qq;
		for (int j = 0; j < 3; ++j) {
			ddq2(j) = (q(j) - sol(j, i)) * bld.GAf.bt / dt / dt - dq(j) * bld.GAf.bt / dt - ddq(j)*(1.0 * bld.GAf.bt / 2.0 - 1.0);
		}
		dq += ddq*(1.0 - bld.GAf.r)*dt + ddq2*bld.GAf.r*dt;
		ddq = ddq2;
		sol(0, i + 1) = q(0);
		sol(1, i + 1) = q(1);
		sol(2, i + 1) = q(2);

		niter = i + 1;
		betaw[0] = sol(0, niter);
		betaw[1] = sol(1, niter);
		betaw[2] = sol(2, niter);

		// exit condition
		if (i > 10)
		{
			_sum = 0;
			for (int j = 0; j < 3; ++j)
			{
				temp = (sol(j, i + 1) - sol(j, i)) / (sol(j, i) + RAD(0.01));
				temp *= temp;
				_sum += temp;
			}
			if (_sum < bld.err_b*bld.err_b)
				break;
		}

	}
	//sol.output("_sol_temp_mr.output", 10);
	if (niter == nitermax - 1)
		;// printf("Warning: Flap solving may not be convergent in Func _hingeflap_fx(). \n");

	beta[0] = betaw[0];
	beta[1] = betaw[2] * sin(windcoord.euler[2]) + betaw[1] * cos(windcoord.euler[2]);
	beta[2] = betaw[2] * cos(windcoord.euler[2]) - betaw[1] * sin(windcoord.euler[2]);

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
	tppcoord.SetCoordinate(euler_temp, "euler");
}

void Rotor::_hingelessflap_fx(void)
{
	Matrix2<myTYPE> m22(3, 3), k22(3, 3), d22(3, 3), kk(3, 3);
	Matrix1<myTYPE> f22(3), dq0(3), q0(3), temp_M(3);
	Matrix1<myTYPE> ddq2(3), ddq(3), dq(3), q(3), qq(3);
	Matrix2<myTYPE> f_temp(3, 4);
	Matrix1<myTYPE> c_temp(4);
	int niter, nitermax = bld.nitermax;
	myTYPE twistt, dt;
	myTYPE euler_temp[3];
	myTYPE sitaw[3], betaw[3];
	myTYPE k1, p2, sb, mb, _sum, temp;
	Matrix2<myTYPE> sol(3, nitermax);

	niter = 0, dt = RAD(bld.dff) / omega;
	bld.GAf.InitCoef(dt);

	k1 = tan(RAD(del));
	mb = m1*radius;
	p2 = 1 + khub / iflap / omega / omega + eflap*mb / iflap + gama*k1*(1 - 4.0 / 3 * eflap);

	sitaw[0] = sita[0];
	sitaw[1] = sita[1] * cos(windcoord.euler[2]) - sita[2] * sin(windcoord.euler[2]);
	sitaw[2] = sita[1] * sin(windcoord.euler[2]) + sita[2] * cos(windcoord.euler[2]);
	twistt = twist(ns - 1) - twist(0);

	m22(0, 0) = m22(1, 1) = m22(2, 2) = 1.0;

	d22(0, 0) = omega*gama*(0.125 - 1.0 / 3.0*eflap + 0.25*eflap*eflap);
	d22(0, 1) = 0;
	d22(0, 2) = -0.25*omega*gama*mul*(1.0 / 3.0 - eflap + eflap*eflap);
	d22(1, 0) = 0;
	d22(1, 1) = omega*gama*(0.125 - 1.0 / 3.0*eflap + 0.25*eflap*eflap);
	d22(1, 2) = 2 * omega;
	d22(2, 0) = -0.5*omega*gama*mul*(1.0 / 3.0 - eflap + eflap*eflap);
	d22(2, 1) = -2 * omega;
	d22(2, 2) = omega*gama*(0.125 - 1.0 / 3.0*eflap + 0.25*eflap*eflap);

	k22(0, 0) = omega*omega*(p2 + 0.25*gama*k1*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	k22(0, 1) = omega*omega*(-0.25*gama*mul)*eflap*(0.5 - eflap);
	k22(0, 2) = omega*omega*(-0.25*gama*k1*mul)*(2.0 / 3.0 - eflap);
	k22(1, 0) = omega*omega*(-0.5*gama*mul)*(1.0 / 3.0 - 0.5*eflap);
	k22(1, 1) = omega*omega*(p2 - 1 + 0.125*gama*mul*mul*k1*(0.5 - eflap + 0.5*eflap*eflap));
	k22(1, 2) = omega*omega*((0.5*gama)*(0.25 - 2.0 / 3.0*eflap + 0.5*eflap*eflap) + (0.125*gama*mul*mul)*(0.5 - eflap + 0.5*eflap*eflap));
	k22(2, 0) = omega*omega*(-0.5*gama*k1*mul)*(2.0 / 3.0 - eflap);
	k22(2, 1) = omega*omega*(-0.5*gama*(0.25 - 2.0 / 3.0*eflap + 0.5*eflap*eflap) + 0.125*gama*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	k22(2, 2) = omega*omega*(p2 - 1 + 3.0 / 8.0*gama*k1*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));

	c_temp(0) = sitaw[0];
	c_temp(1) = twistt;
	c_temp(2) = -sitaw[1];
	c_temp(3) = -sitaw[2];

	f_temp(0, 0) = 0.5*gama*(0.25 - 1.0 / 3.0*eflap + 0.5*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	f_temp(0, 1) = 0.5*gama*(0.2 - 0.25*eflap + 0.5*mul*mul*(1.0 / 3.0 - 0.5*eflap));
	f_temp(0, 2) = 0;
	f_temp(0, 3) = -0.5*gama*mul*(1.0 / 3.0 - 0.5*eflap);
	f_temp(1, 0) = 0;
	f_temp(1, 1) = 0;
	f_temp(1, 3) = 0;
	f_temp(1, 2) = 0.5*gama*(0.25 - 1.0 / 3.0*eflap + 0.25*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));
	f_temp(2, 0) = -0.5*gama*mul*(2.0 / 3.0 - eflap);
	f_temp(2, 1) = -0.5*gama*mul*(0.5 - 2.0 / 3.0*eflap);
	f_temp(2, 2) = 0;
	f_temp(2, 3) = 0.5*gama*(0.25 - 1.0 / 3.0*eflap + 0.75*mul*mul*(0.5 - eflap + 0.5*eflap*eflap));


	f22 = (f_temp.matrixmultiplyP2(c_temp)) * omega*omega;

	c_temp(0) = omgw[0];
	c_temp(1) = omgw[1];
	c_temp(2) = domgw[0];
	c_temp(3) = domgw[1];

	f_temp(0, 0) = 0.125*gama*mul / omega*(2.0 / 3.0 - eflap);
	f_temp(0, 1) = f_temp(0, 2) = f_temp(0, 3) = 0;
	f_temp(1, 0) = -2.0 / omega*(1.0 + eflap*mb / iflap);
	f_temp(1, 1) = -0.5*gama / omega*(0.25 - 1.0 / 3.0*eflap);
	f_temp(1, 2) = 0;
	f_temp(1, 3) = -1.0 / omega / omega;
	f_temp(2, 0) = -0.5*gama / omega*(0.25 - 1.0 / 3.0*eflap);
	f_temp(2, 1) = 2.0 / omega*(1.0 + eflap*mb / iflap);
	f_temp(2, 2) = -1.0 / omega / omega;
	f_temp(2, 3) = 0;

	f22 += (f_temp.matrixmultiplyP2(c_temp)) * omega*omega;
	f22.v_p[0] += omega*omega*0.5*gama*(1.0 / 3.0 - 0.5*eflap)*lambdh_ag;
	f22.v_p[2] -= omega*omega*(0.5*gama*mul)*(0.5 - eflap + 0.5*eflap*eflap)*lambdh_ag;

	//if (!si_unit)
	//	f22.v_p[0] += mb / iflap *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
	//else
	//	f22.v_p[0] += mb / iflap*(-9.8 +dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);


	_sum = 0;
	temp = 0;
	sb = 8.0*(p2 - 1.0) / gama;
	q0(1) = -1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	q0(2) = -1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb + 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);
	q0(0) = gama / p2*(sitaw[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sitaw[2]) - lambdh_ag / 6.0);

	dq0.v_p[0] = 0.0;
	dq0.v_p[1] = 0;
	dq0.v_p[2] = omega*q0(2);

	temp_M = d22.matrixmultiplyP2(dq0) + k22.matrixmultiplyP2(q0) - f22;
	Msolver(m22.v_p, temp_M.v_p, 3);
	ddq = temp_M*(-1);
	dq = dq0;
	q(0) = q0(0);
	q(1) = q0(1);
	q(2) = q0(2);

	sol(0, 0) = q(0);
	sol(1, 0) = q(1);
	sol(2, 0) = q(2);

	for (int i = 0; i < nitermax - 1; ++i) {
		qq = f22 - (k22.matrixmultiplyP2(q)) * bld.GAf.af;
		qq += m22.matrixmultiplyP2(q*bld.GAf.c0 + dq*bld.GAf.c2 + ddq*bld.GAf.c3);
		qq += d22.matrixmultiplyP2(q*bld.GAf.c1 + dq*bld.GAf.c4 + ddq*bld.GAf.c5);
		kk = k22*bld.GAf.ck + m22*bld.GAf.c0 + d22*bld.GAf.c1;
		Msolver(kk.v_p, qq.v_p, 3);
		q = qq;
		for (int j = 0; j < 3; ++j) {
			ddq2(j) = (q(j) - sol(j, i)) * bld.GAf.bt / dt / dt - dq(j) * bld.GAf.bt / dt - ddq(j)*(1.0 * bld.GAf.bt / 2.0 - 1.0);
		}
		dq += ddq*(1.0 - bld.GAf.r)*dt + ddq2*bld.GAf.r*dt;
		ddq = ddq2;
		sol(0, i + 1) = q(0);
		sol(1, i + 1) = q(1);
		sol(2, i + 1) = q(2);

		niter = i + 1;
		betaw[0] = sol(0, niter);
		betaw[1] = -sol(1, niter);
		betaw[2] = -sol(2, niter);

		// exit condition
		if (i > 10)
		{
			_sum = 0;
			for (int j = 0; j < 3; ++j)
			{
				temp = (sol(j, i + 1) - sol(j, i)) / (sol(j, i) + RAD(0.01));
				temp *= temp;
				_sum += temp;
			}
			if (_sum < bld.err_b*bld.err_b)
				break;
		}

	}
	//sol.output("_sol_temp_mr.output", 10);
	if (niter == nitermax - 1)
		;// printf("Warning: Flap solving may not be convergent in Func _hingeflap_fx(). \n");

	
	beta[0] = betaw[0];
	beta[1] = betaw[2] * sin(windcoord.euler[2]) + betaw[1] * cos(windcoord.euler[2]);
	beta[2] = betaw[2] * cos(windcoord.euler[2]) - betaw[1] * sin(windcoord.euler[2]);

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // 桨盘平面旋转，是因为Glauert入流模型是没有横向速度的
	tppcoord.SetCoordinate(euler_temp, "euler");
}

void Rotor::_flapMBC(double lamb_1s, double lamb_1c)
{
	Matrix2<myTYPE> MM(4, 4), CC(4, 4), DD(4, 4), DD_temp(4, 4);
	Matrix1<myTYPE> HH(4), temp_M(4);
	Matrix1<myTYPE> dq0(4), q0(4), ddq2(4), ddq(4), dq(4), q(4), qq(4);
	int niter, nitermax = bld.nitermax;
	myTYPE twistt, dt;
	myTYPE euler_temp[3];
	myTYPE sitaw[3], betaw[3];
	myTYPE p2, sb, mb, _sum, temp;
	Matrix2<myTYPE> sol(4, nitermax);

	mb = m1*radius;
	p2 = 1 + khub / iflap / omega / omega;
	sitaw[0] = sita[0];
	sitaw[1] = sita[1] *cos(windcoord.euler[2]) - sita[2] * sin(windcoord.euler[2]);
	sitaw[2] = sita[1] *sin(windcoord.euler[2]) + sita[2] * cos(windcoord.euler[2]);
	twistt = twist(ns - 1) - twist(0);

	niter = 0, dt = RAD(bld.dff) / omega;
	bld.GAf.InitCoef(dt);
	
	MM(0, 0) = MM(1, 1) = MM(2, 2) = MM(3, 3) = 1.0;

	CC(0, 0) = gama / 8.0, CC(0, 3) = gama / 8.0*2.0 / 3.0*mul;
	CC(1, 1) = gama / 8.0;
	CC(2, 2) = gama / 8.0, CC(2, 3) = 2.0;
	CC(3, 0) = gama / 8.0*4.0 / 3.0*mul;
	CC(3, 2) = -2.0, CC(3, 3) = gama / 8.0;

	DD(0, 0) = p2;
	DD(1, 1) = p2;
	DD(2, 0) = gama / 8.0*4.0 / 3.0*mul;
	DD(2, 2) = p2 - 1;
	DD(2, 3) = gama / 8.0*(1 + mul*mul / 2.0);
	DD(3, 2) = -gama / 8.0*(1 - mul*mul / 2.0);
	DD(3, 3) = p2 - 1;

	HH(0) = sitaw[0] * (1 + mul*mul) + 4 * twistt*(1.0 / 5.0 + mul*mul / 6.0);
	HH(0) += 4.0 / 3.0*mul*sitaw[2] + 4.0 / 3.0*lambdh_ag + 2.0 / 3.0*mul*(omgw[0] / omega - lamb_1s);
	//HH(0) -= mb*9.8 / iflap / omega / omega;
	HH(2) = 16 / gama*(omgw[0] / omega + domgw[1] / omega / omega / 2.0) + sitaw[1] * (1 + mul*mul / 2.0) + (omgw[1] / omega - lamb_1c);
	HH(3) = -16 / gama*(omgw[1] / omega - domgw[0] / omega / omega / 2.0) + 8.0 / 3.0*mul*sitaw[0] + 2 * mul*twistt;
	HH(3) += sitaw[2] * (1 + 3.0 / 2.0*mul*mul) + 2 * mul*lambdh_ag + (omgw[0] / omega - lamb_1s);


	CC = CC * omega;
	DD = DD * omega*omega;
	HH = HH*omega*omega*gama / 8.0;

	_sum = 0;
	temp = 0;
	sb = 8.0*(p2 - 1.0) / gama;
	q0(1) = 0.0;
	q0(2) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	q0(3) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb + 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);
	q0(0) = gama / p2*(sitaw[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sitaw[2]) - lambdh_ag / 6.0);

	dq0.v_p[0] = 0.0;
	dq0.v_p[1] = 0.0;
	dq0.v_p[2] = 0.0;
	dq0.v_p[3] = omega*q0(3);

	temp_M = CC.matrixmultiplyP2(dq0) + DD.matrixmultiplyP2(q0) - HH;
	Msolver(MM.v_p, temp_M.v_p, 4);
	ddq = temp_M*(-1);
	dq = dq0;
	for (int i = 0; i < 4; i++)
	{
		q(i) = q0(i);
		sol(i, 0) = q(i);
	}

	for (int i = 0; i < nitermax - 1; ++i) {
		qq = HH - (DD.matrixmultiplyP2(q)) * bld.GAf.af;
		qq += MM.matrixmultiplyP2(q*bld.GAf.c0 + dq*bld.GAf.c2 + ddq*bld.GAf.c3);
		qq += CC.matrixmultiplyP2(q*bld.GAf.c1 + dq*bld.GAf.c4 + ddq*bld.GAf.c5);
		DD_temp = DD*bld.GAf.ck + MM*bld.GAf.c0 + CC*bld.GAf.c1;
		Msolver(DD_temp.v_p, qq.v_p, 4);
		q = qq;
		for (int j = 0; j < 4; ++j) {
			ddq2(j) = (q(j) - sol(j, i)) * bld.GAf.bt / dt / dt - dq(j) * bld.GAf.bt / dt - ddq(j)*(1.0 * bld.GAf.bt / 2.0 - 1.0);
		}
		dq += ddq*(1.0 - bld.GAf.r)*dt + ddq2*bld.GAf.r*dt;
		ddq = ddq2;
		sol(0, i + 1) = q(0);
		sol(1, i + 1) = q(1);
		sol(2, i + 1) = q(2);
		sol(3, i + 1) = q(3);

		niter = i + 1;
		betaw[0] = sol(0, niter);
		betaw[1] = sol(2, niter);
		betaw[2] = sol(3, niter);

		// exit condition
		if (i > 10)
		{
			_sum = 0;
			for (int j = 0; j < 4; ++j)
			{
				temp = (sol(j, i + 1) - sol(j, i)) / (sol(j, i) + RAD(0.01));
				temp *= temp;
				_sum += temp;
			}
			if (_sum < bld.err_b*bld.err_b)
				break;
		}
	}
	//sol.output("_sol_temp_mr.output", 10);
	if (niter == nitermax - 1)
		;// printf("Warning: Flap solving may not be convergent in Func _hingeflap_fx(). \n");

	beta[0] = betaw[0];
	beta[1] = betaw[2]*sin(windcoord.euler[2]) + betaw[1] * cos(windcoord.euler[2]);
	beta[2] = betaw[2]*cos(windcoord.euler[2]) - betaw[1] * sin(windcoord.euler[2]);

	euler_temp[0] = -beta[2];
	euler_temp[1] = beta[1];
	//euler_temp[2] = 0.0;
	euler_temp[2] = Atan2(vel[1], vel[0]); // 桨盘平面旋转，是因为Glauert入流模型是没有横向速度的
	tppcoord.SetCoordinate(euler_temp, "euler");
}


double Rotor::_aerodynamics(double lambtpp, double *veltpp)
{
	double kx, ky, ka; 
	double ct, twistt;
	int Iid_betaw = 0;
	double lambdi_temp = 0;
	double _af_temp[3], _am_temp[3];

	_af_temp[0] = _af_temp[1] = _af_temp[2] = 0;
	_am_temp[0] = _am_temp[1] = _am_temp[2] = 0;
	twistt = twist(ns - 1) - twist(0);

	switch (adyna)
	{
	case PWake:
		_setairfm_sp(_af_temp, _am_temp);
		for (int i = 0; i < 3; i++)
		{
			airforce[i] = _af_temp[i];
			airmoment[i] = _am_temp[i];
		}
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;

		if (haveGeo && haveStr)
		{
			_bladePosition();
			//_wakeIndVelCalc();
			_wakeInducedVelMP(nb);
		}
		else
		{
			lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
			lambdh.setvalue(lambdh_ag);
			lambdi.setvalue(lambdi_ag);
		}
		break;
	case LinearInflow:
		_setairfm_sp(_af_temp, _am_temp);
		for (int i = 0; i < 3; i++)
		{
			airforce[i] = _af_temp[i];
			airmoment[i] = _am_temp[i];
		}
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;

		lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		
		ka = Atan2(veltpp[0] / vtipa, -veltpp[2] / vtipa + lambdi_ag);
		ky = -2 * veltpp[0] / vtipa;
		kx = tan(ka / 2);
		// non-uniform inflow model is doing in TPP wind coordinate
		Iid_betaw = windcoord.euler[2] / (2 * PI / nf);
		for (int j = lambdi.NJ - 1; j >= 0; j--)
		{
			for (int i = lambdi.NI - 1; i >= 0; i--)
			{
				int iz = i + Iid_betaw;
				while (iz >= nf) iz -= nf;
				while (iz < 0) iz += nf;
				lambdi_temp = 1 + kx*rastation(iz, j)*cos(azstation(iz, j));
				lambdi_temp += ky*rastation(iz, j)*sin(azstation(iz, j));
				lambdi_temp *= lambdi_ag;
				lambdi(i, j) = lambdi_temp;
				lambdh(i, j) = lambdi_temp*cos(beta[1])*cos(beta[2]) - vel[2] / vtipa; 
			}
		}
		//lambdh.output("_lambdh.output", 10);
		break;
	case Averaged:
		_setairfm_sp(_af_temp, _am_temp);
		for (int i = 0; i < 3; i++)
		{
			airforce[i] = _af_temp[i];
			airmoment[i] = _am_temp[i];
		}
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;

		lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		lambdh.setvalue(lambdh_ag);
		lambdi.setvalue(lambdi_ag);
		break;
	case Simple:
		// simple ct
		ct = sigma*a0*0.5*(sita[0] / 3.0*(1 + 1.5*mul*mul) + 0.25*twistt*(1 + mul*mul) + 0.5*mul*beta[1] - 0.5*lambtpp);
		lambdi_ag += ct / 2.0 / sqrt(mul*mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;
		airforce[2] = -ct*amb.rho*PI*radius*radius*vtipa*vtipa;

		lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		lambdh.setvalue(lambdh_ag);
		lambdi.setvalue(lambdi_ag);
		break;
	default:
		_setairfm_sp(_af_temp, _am_temp);
		for (int i = 0; i < 3; i++)
		{
			airforce[i] = _af_temp[i];
			airmoment[i] = _am_temp[i];
		}
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;

		lambdh_ag = lambdi_ag * cos(beta[1]) * cos(beta[2]) - vel[2] / vtipa;
		lambdh.setvalue(lambdh_ag);
		lambdi.setvalue(lambdi_ag);
		break;
	}

	// veltpp: tpp plane velocity with copter
	return lambdi_ag - veltpp[2] / vtipa;
}

void Rotor::_windcoordVel(double v[3], double dv[3])
{
	v[0] = v[1] = v[2] = 0;
	dv[0] = dv[1] = dv[2] = 0;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			v[i] += windcoord.Ttransf[i][j] * vel[j];
			dv[i] += windcoord.Ttransf[i][j] * dvel[j];
		}
	}
}

void Rotor::_windcoordOmg(double w[3], double dw[3])
{
	w[0] = w[1] = w[2] = 0;
	dw[0] = dw[1] = dw[2] = 0;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			w[i] += windcoord.Etransf[i][j] * omg[j];
			dw[i] += windcoord.Etransf[i][j] * domg[j];
		}
	}
}

void Rotor::_velTransform(double v[3], double dv[3], Coordinate & coord)
{
	v[0] = v[1] = v[2] = 0;
	dv[0] = dv[1] = dv[2] = 0;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			v[i] += coord.Ttransf[i][j] * vel[j];
			dv[i] += coord.Ttransf[i][j] * dvel[j];
		}
	}
}

void Rotor::_setairfm_sp(double f[3], double m[3])
{
	double it = 0;
	double ia = 0;
	double b = 0;
	double db = 0;
	Matrix1<double> _dt(ns), _yf(ns), _hf(ns), _az(ns), _lambdi(ns), az(nf);
	Matrix1<double> _dfx(ns), _dfz(ns), _dfr(ns), ra(ns), ra2(ns), ra1(ns);
	Matrix1<double> _dD(ns);
	Matrix1<int> id_ns = step(0, ns - 1);
	Matrix1<int> id_nf = step(0, nf - 1);
	ra = rastation(0, id_ns);
	az = azstation(id_nf, 0);

	f[0] = f[1] = f[2] = 0.0;
	m[0] = m[1] = m[2] = 0.0;

	power = power_c = power_i = power_f = power_o; 
	torque = torque_c = torque_i = torque_f = torque_o;
	
	if (adyna == Averaged)
	{
		power_iid = 0;
		torque_iid = 0;
	}

	for (int i = 0; i < nf; ++i) {
		ia = i * 2 * PI / nf;
		it = ia / omega;
		b = beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
		db = (-beta[1] * sin(ia) + beta[2] * cos(ia)) * omega;

		_setairfm(_dfx, _dfz, _dfr, it, b, db, i);

		_dt = _dfz*cos(b);
		_yf = _dfx*cos(ia) - _dfr*sin(ia);
		_hf = _dfx*sin(ia) + _dfr*cos(ia);

		f[0] -= _hf.sum()*nb / nf;
		f[1] -= _yf.sum()*nb / nf;
		f[2] -= _dt.sum()*nb / nf;

		//cout << f[2] << endl;

		ra2 = (ra - eflap)*sin(b);
		ra1 = (ra - eflap)*cos(b) + eflap;
		//switch (hingetype)
		//{
		//case Hinged:
		//	m[0] -= (_dt * ra1 * sin(ia) + _yf * ra2).sum()*radius*nb / nf;
		//	m[1] -= (_dt * ra1 * cos(ia) + _hf * ra2).sum()*radius*nb / nf;
		//	break;
		//case Cantilever:
		//	printf("Undefined Cantilever type in _setairfm_sp()");
		//	system("pause");
		//	break;
		//default:
		//	m[0] -= (_dt * ra1 * sin(ia) + _yf * ra2).sum()*radius*nb / nf;
		//	m[1] -= (_dt * ra1 * cos(ia) + _hf * ra2).sum()*radius*nb / nf;
		//	break;
		//}
		m[2] += (_dfx * ra1).sum()*radius*nb / nf;

		// power computation
		_az.setvalue(ia);
		_lambdi = lambdi.interplinear_fast(az, ra, _az, ra);
		_dD = _cd*_factor / mcos(_inflow);
		power_i += (_dt*_lambdi).sum() +db / omega * ((ra - eflap)*_dfz).sum();
		power_o += (_dD*sin(ia)*mul + _dD*ra1).sum();
		
		if (adyna == Averaged)
			power_iid += db / omega*((ra - eflap)*_dfz).sum();


		// save disk plane data
		for (int j = 0; j < ns; j++)
		{
			cl.v_p[j*cl.NI + i] = _cl(j);
			cd.v_p[j*cd.NI + i] = _cd(j);
			inflow.v_p[j*inflow.NI + i] = _inflow(j);
			incidn.v_p[j*incidn.NI + i] = _incidn(j);
			ua.v_p[j*ua.NI + i] = _ua(j);
			dt.v_p[j*dt.NI + i] = _dt(j);
		}
	}

	//cl.output("cl.output", 10);
	//cd.output("cd.output", 10);
	//inflow.output("inflow.output", 10);
	//incidn.output("AOA.output", 10);
	//ua.output("Ma.output", 10);
	//lambdi.output("lambdi.output", 10);
	//dt.output("dt.output", 10);

	switch (hingetype)
	{
	case Teeter:
		m[0] = m[1] = 0;
		break;
	case Hinged:
		m[0] = -nb*khub / 2 * beta[2];
		m[1] = -nb*khub / 2 * beta[1];
	case Hingeless:
		m[0] = -nb*khub / 2 * beta[2];
		m[1] = -nb*khub / 2 * beta[1];
		break;
	default:
		break;
	}

	//lambdi.output("temp_lambdi.output", 10);
	power_i *= vtipa*nb / nf;
	power_o *= vtipa*nb / nf;
	power_f = f[0] * vel[0];
	power_c = f[2] * vel[2];
	power = m[2] * omega;
	torque_i = power_i / omega;

	torque_o = power_o / omega;
	torque_f = power_f / omega;
	torque_c = power_c / omega;
	torque = power / omega;

	if (adyna == Averaged)
	{
		power_iid = (-f[2] * lambdi_ag + power_iid * nb / nf) * vtipa;
		torque_iid = power_iid / omega;
	}
	if (!si_unit)
	{
		power_i *= HORSEPOWER;
		power_o *= HORSEPOWER;
		power_f *= HORSEPOWER;
		power_c *= HORSEPOWER;
		power *= HORSEPOWER;
		torque_i = power_i / omega;

		torque_o = power_o / omega;
		torque_f = power_f / omega;
		torque_c = power_c / omega;
		torque = power / omega;

		if (adyna == Averaged)
		{
			power_iid *= HORSEPOWER;
			torque_iid = power_iid / omega;
		}
	}
	else
	{
		power_i /= 1000.0;
		power_o /= 1000.0;
		power_f /= 1000.0;
		power_c /= 1000.0;
		if (adyna == Averaged)
			power_iid /= 1000.0;
	}

}

void Rotor::_setairfm(Matrix1<double> &_dfx, Matrix1<double> &_dfz, Matrix1<double> &_dfr, const double &it, const double &b, const double &db, const int iz)
{
	Matrix1<double> _ut(ns), _up(ns);
	double _sfth = 0;
	double ia = 3 * PI;
	Matrix1<double> _dr(ns);

	_setbladevelc(_ut, _up, ia, it, b, db);
	_up = _up*vtipa;
	_ut *= mcos(sweep);
	_ut = _ut*vtipa;

	_ua = msqrt(_ut*_ut + _up*_up);
	_ua = _ua / amb.vsound;
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

	_sfth = sita[0] + sita[1] * cos(ia) + sita[2] * sin(ia) + pitchroot - tan(RAD(del))*b;

	// aoa
	_inflow = atan2(_up, _ut);
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
	for (int i = ns - 1; i >= 0; --i)
		cirlb(iz%nf, i) = 0.5 * _ua(i) * _cl(i) * chord(i) * amb.vsound;

	// air forces and moments
	_factor = _ua*_ua;
	_dr.setvalue(rastation(0, 1) - rastation(0, 0));
	_factor *= chord * _dr *amb.vsound*amb.vsound*0.5*amb.rho* radius;

	_dfx = (_cl * msin(_inflow) + _cd * mcos(_inflow)) * _factor;
	_dfz = (_cl * mcos(_inflow) - _cd * msin(_inflow)) * _factor;
	_dfr = _dfz * (sin(b)) * (-1.0);

}

void Rotor::_setbladevelc(Matrix1<double> &_ut, Matrix1<double> &_up, double &ia, const double &it, const double &b, const double &db)
{
	Matrix1<double> _az(ns), _ra(ns), _lambdh(ns), _lambdx(ns), _lambdy(ns);
	Matrix1<int> id_ns(ns), id_nf(nf);
	id_ns = step(0, ns - 1);
	id_nf = step(0, nf - 1);

	ia = it * omega;
	for (;;) {
		if (ia > 2 * PI) ia -= 2 * PI;
		else if (ia < -2 * PI) ia += 2 * PI;
		else { break; }
	}
	if (ia > azstation(nf - 1, 0))
		ia = 0;
	_az.setvalue(ia);
	_ra = rastation(0, id_ns);

	_lambdh = lambdh.interplinear_fast(azstation(id_nf, 0), _ra, _az, _ra);

#ifdef TEST_MODE
	_lambdh.outputs("_lambdh.output", 4);
#endif // TEST_MODE

	// air velocity
	_ut = (_ra - eflap) * cos(b) + eflap + (sin(ia) * vel[0] + cos(ia) * vel[1]) / vtipa;
	_up = (_ra - eflap) * db / omega + (cos(ia) * vel[0] - sin(ia) * vel[1]) * sin(b) / vtipa + _lambdh * cos(b);
	//_ut += _ra*(omg[0] / omega*cos(ia) - omg[1] / omega*sin(ia)) * b;
	_up -= _ra*(omg[0] / omega*sin(ia) + omg[1] / omega*cos(ia)); 

	//_ut = (_ra - eflap) * cos(b) + eflap + sin(ia) * vel[0] / vtipa;
	//_up = (_ra - eflap) * db / omega + cos(ia) * vel[0] * sin(b) / vtipa + _lambdh * cos(b);

}

void Rotor::GetPower(myTYPE p[6], myTYPE t[6])
{
	p[0] = power;
	p[1] = power_i;
	p[2] = power_o;
	p[3] = power_f;
	p[4] = power_c;
	p[5] = power_iid;

	t[0] = torque;
	t[1] = torque_i;
	t[2] = torque_o;
	t[3] = torque_f;
	t[4] = torque_c;
	t[5] = torque_iid;
}

double Rotor::GetLambdi(void)
{
	return lambdi_ag * vtipa;
}

void Rotor::DiskOutput(string s)
{
	cl.output(s + "_cl.output", 10);
	cd.output(s + "_cd.output", 10);
	inflow.output(s + "_inflow.output", 10);
	incidn.output(s + "_AOA.output", 10);
	ua.output(s + "_Ma.output", 10);
	lambdi.output(s + "_lambdi.output", 10);
	dt.output(s + "_dt.output", 10);
	//printf("Max inflow angle %f\n", inflow.findmax());
	//printf("Max AoA angle %f\n", incidn.findmax());
}


bool BladeSolver::isGenArfExit(int niter)
{
	myTYPE _sum, temp;
	if (niter >= 2 * nperiod - 1)
	{
		_sum = 0;
		for (int j = 0; j < nperiod; ++j)
		{
			temp = sol(niter - nperiod + 1 + j) - sol(niter - 2 * nperiod + 1 + j);
			_sum += temp*temp;
		}
		if (_sum < err_b*err_b && niter%nperiod==0) { return true; }
	}
	return false;
}

