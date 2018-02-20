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

	/*cout << endl;
	printf("%d rotor induced velocity iter count: %d \n", type, niter_w);
	printf("%d rotor flaps (degs): %f, %f, %f\n", type, DEG(beta[0]), DEG(beta[1]), DEG(beta[2]));
	printf("lambdi_ag: %f, lambdt_ag: %f, lambdh_ag: %f \n", lambdi_ag, lambdt_ag, lambdh_ag);
	cout << endl;
	printf("%d rotor airdynamics: \n", type);
	printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;*/

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
	myTYPE dveltpp[3] = { 0.0,0.0,0.0 };
	myTYPE euler_temp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp, twistt;

	twistt = twist(ns - 1) - twist(0);

	lambtpp[0] = -0.01;
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

		_velTransform(veltpp, dveltpp, tppcoord);

		euler_temp[2] = atan2(veltpp[1], -veltpp[0]); // wind coordinate
		veltpp[0] = veltpp[0] * cos(euler_temp[2]) - veltpp[1] * sin(euler_temp[2]);
		veltpp[1] = 0;
		
		lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
		lambtpp[i] = lambdt_ag;
		

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
		printf("Warning: Flap solving may not be convergent in Func _teeterflap_rt(). \n");

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
	myTYPE dt, twistt;
	myTYPE euler_temp[3], sitaw[3], betaw[3];
	myTYPE k1, p2, sb, _sum, mb, temp;
	Matrix2<myTYPE> sol(2, nitermax);

	
	mb = m1*radius;
	k1 = tan(RAD(del));
	p2 = 1 + khub / iflap / omega / omega + gama*k1;

	sitaw[0] = sita[0];
	sitaw[1] = sita[1] * cos(betawind) - sita[2] * sin(betawind);
	sitaw[2] = sita[1] * sin(betawind) + sita[2] * cos(betawind);
	twistt = twist(ns - 1) - twist(0);

	niter = 0, dt = RAD(bld.dff) / omega;
	bld.GAf.InitCoef(dt);

	m22(0, 0) = m22(1, 1) = 1.0;

	d22(0, 0) = omega*gama / 8.0;
	d22(0, 1) = 2 * omega;
	d22(1, 0) = -2 * omega;
	d22(1, 1) = omega*gama / 8.0;

	k22(0, 0) = omega*omega*(p2 - 1 + gama*k1*mul*mul / 16.0);
	k22(0, 1) = omega*omega*gama / 8.0*(1 + 0.5*mul*mul);
	k22(1, 0) = -omega*omega*gama / 8.0*(1 - 0.5*mul*mul);
	k22(1, 1) = omega*omega*(p2 - 1 + 3.0 / 16.0*gama*k1*mul*mul);

	//f22.v_p[0] = omega*omega*(gama / 8.0*sitaw[1] * (1 + 0.5*mul*mul));
	//f22.v_p[1] = omega*omega*(1.0 / 3.0*gama*mul*sitaw[0] + 1.0 / 4.0*gama*mul*(twistt)+gama / 8.0*sitaw[2] * (1 + 1.5*mul*mul));

	//f22.v_p[0] += omega*omega*(2.0 / omega*omgw[0] + gama / 8.0 / omega*omgw[1] + 1.0 / omega / omega*domgw[1]);
	//f22.v_p[1] += omega*omega*(gama / 8.0 / omega*omgw[0] - 2.0 / omega*omgw[1] + 1.0 / omega / omega*domgw[0]);

	//f22.v_p[1] += omega*omega*(0.25*gama*mul)*lambdh_ag;

	f22(0) = 16 / gama*(omgw[0] / omega + domgw[1] / omega / omega / 2.0) + sitaw[1] * (1 + mul*mul / 2.0) + (omgw[1] / omega + 0);
	f22(1) = -16 / gama*(omgw[1] / omega - domgw[0] / omega / omega / 2.0) + 8.0 / 3.0*mul*sitaw[0] + 2 * mul*twistt;
	f22(1) += sitaw[2] * (1 + 3.0 / 2.0*mul*mul) + 2 * mul*lambdh_ag + (omgw[0] / omega + 0);

	f22 = f22*omega*omega*gama / 8.0;

	if (!si_unit)
		f22(0) += m1 / iflap *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
	else
		f22(0) += m1 / iflap *(-9.8 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);

	_sum = 0;
	temp = 0;

	sb = 8.0*(p2 - 1.0) / gama;
	//q0(0) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	//q0(1) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb + 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);
	//dq0.v_p[0] = 0.0;
	//dq0.v_p[1] = omega * q0(1);// omega*sol(1, 0);

	q0(0) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	q0(1) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb + 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);

	dq0(0) = 0.0;
	dq0(1) = omega*q0(1);

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

		niter = i + 1;
		betaw[0] = precone;
		betaw[1] = sol(0, niter);
		betaw[2] = sol(1, niter);

		// exit condition
		if (i > 10)
		{
			_sum = 0;
			for (int j = 0; j < 2; ++j) {
				temp = (sol(j, i + 1) - sol(j, i)) / (sol(j, i) + RAD(0.01));
				temp *= temp;
				_sum += temp;
			}

			if (_sum < bld.epsb*bld.epsb) { break; }
		}
	}

	if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _teeterflap_fx(). \n");

	beta[0] = betaw[0];
	beta[1] = betaw[2] * sin(betawind) + betaw[1] * cos(betawind);
	beta[2] = betaw[2] * cos(betawind) - betaw[1] * sin(betawind);

	euler_temp[0] = beta[2]; // -beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0.0;
	//euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
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
	if (!si_unit)
		qt += m1 *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
	else
		qt += m1 *(-9.8 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
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
			if (!si_unit)
				qt += m1 *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
			else
				qt += m1 *(-9.8 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
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
		f22(0) += m1 / iflap *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
	else
		f22(0) += m1 / iflap *(-9.8 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);

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
			if (_sum < bld.epsb*bld.epsb)
				break; 
		}

	}
	//sol.output("_sol_temp_tr.output", 10);
	if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _hingeflap_fx(). \n");

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
	Matrix2<myTYPE> lambdiold(nf, ns);

	bld.Marf = iflap;
	bld.Carf = 0;
	bld.Ka = omega*omega*iflap+khub;

	bld._GenArfPrepare(RAD(bld.dff) / omega);

	if (si_unit)
		lambtpp[0] = -sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = -sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);

	if (airfoil != LBStallMethod)
	{
		for (int i = 1; i < itermax; i++)
		{
			_hingelessflap_rt();

			_velTransform(veltpp, dveltpp, tppcoord);

			euler_temp[2] = atan2(veltpp[1], -veltpp[0]); // wind coordinate
			veltpp[0] = veltpp[0] * cos(euler_temp[2]) - veltpp[1] * sin(euler_temp[2]);
			veltpp[1] = 0;

			lambdiold = lambdi;
			lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
			lambtpp[i] = lambdt_ag;

			iter = i;
			if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w)
				break;
		}
	}
	else
	{
		double af[3] = { 0,0,0 };
		double mf[3] = { 0,0,0 };
		for (int i = 1; i < itermax; i++)
		{
			if (_hingelessflap_rt(af, mf))
			{
				for (int i = 0; i < 3; i++)
				{
					airforce[i] = af[i];
					airmoment[i] = mf[i];
				}
			}
			else
			{
				for (int i = 0; i < 3; i++)
				{
					airforce[i] += 1;
					airmoment[i] += 1;
				}
			}
			_velTransform(veltpp, dveltpp, tppcoord);
			euler_temp[2] = atan2(veltpp[1], -veltpp[0]); // wind coordinate
			veltpp[0] = veltpp[0] * cos(euler_temp[2]) - veltpp[1] * sin(euler_temp[2]);
			veltpp[1] = 0;

			lambdiold = lambdi;
			lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
			lambtpp[i] = lambdt_ag;

			iter = i;
			if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w)
				break;
		}
		if (secLB)
			airfoil = C81Table;

	}
	niter_w = iter;
	Errw2 = (lambdt_ag / lambtpp[niter_w - 1] - 1)*(lambdt_ag / lambtpp[niter_w - 1] - 1);
	
	// reset wake flags
	if (haveGeo || haveStr)
	{
		haveGeo = false;
		haveStr = false;
		Errw2 = 0;
		for (int i = 0; i < nf; i++)
			for (int j = 0; j < ns; j++)
				Errw2 += pow(lambdi(i, j) - lambdiold(i, j), 2);		
	}
	
	monitor.errw2 = Errw2;
	monitor.Countsw = iter;
#ifdef OUTPUT_MODE_1
	printf("Flap in _hingelessdynamics_rt() counts: %d \n", niter_w);
	printf("M Rotor Errb2 = %e\n", monitor.errb2);
	printf("M Rotor Errw2 = %e\n", monitor.errw2);
#endif // OUTPUT_MODE_1
}

void Rotor::_hingelessdynamics_fx(void)
{
	const int itermax = 20;
	int iter = 1;
	myTYPE lambtpp[itermax] = { 0 };
	myTYPE err_w = 0.0001;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE dveltpp[3] = { 0.0,0.0,0.0 };
	myTYPE ct, ch_tpp, cy_tpp;
	double lamb_1s, lamb_1c;
	myTYPE euler_temp[3] = { 0,0,0 };

	if (si_unit)
		lambtpp[0] = -sqrt(0.5*t0*9.8 / amb.rho / PI / radius / radius / vtipa / vtipa);
	else
		lambtpp[0] = -sqrt(0.5*t0*UNIT_CONST / amb.rho / PI / radius / radius / vtipa / vtipa);

	lambdi_ag = lambtpp[0];
	lambdh_ag = lambdi_ag;
	lambdh.setvalue(lambdh_ag);
	lambdi.setvalue(lambdi_ag);

	lamb_1s = lamb_1c = 0.0;
	for (int i = 1; i < itermax; ++i)
	{
		//_hingelessflap_fx();
		_flapMBC(lamb_1s*0.0, lamb_1c*1.0);
		//_flapWang();

		_velTransform(veltpp, dveltpp, tppcoord);

		euler_temp[2] = atan2(veltpp[1], -veltpp[0]); // wind coordinate
		veltpp[0] = veltpp[0] * cos(euler_temp[2]) - veltpp[1] * sin(euler_temp[2]);
		veltpp[1] = 0;

		lamb_1s = 2 * veltpp[0] / vtipa*lambdi_ag;
		lamb_1c = tan(0.5*atan2(-veltpp[0] / vtipa, veltpp[2] / vtipa - lambdi_ag))*lambdi_ag;

		lambdt_ag = _aerodynamics(lambtpp[i - 1], veltpp);
		lambtpp[i] = lambdt_ag;

		iter = i;
		if (Abs(lambdt_ag / lambtpp[i - 1] - 1) < err_w) 
			break; 
	}
	//lambtpp.outputs("_lambtpp.output", 4);
	niter_w = iter;
	Errw2 = (lambdt_ag / lambtpp[niter_w - 1] - 1)*(lambdt_ag / lambtpp[niter_w - 1] - 1);

	if (haveGeo)
		haveGeo = false;
	if (haveStr)
		haveStr = false;
}

void Rotor::_hingelessflap_rt(void)
{
	int niter, nitermax;
	myTYPE dt, af;
	myTYPE qt, dq, q, t;
	Matrix1<myTYPE> dfx(ns), dfr(ns), dfz(ns), ra(ns);
	Matrix1<int> id_ns;
	myTYPE euler_temp[3];
	Matrix1<myTYPE> azmuth(bld.nitermax);

	niter = 0, nitermax = bld.nitermax, af = bld.GAf.af;
	dt = RAD(bld.dff) / omega;
	id_ns = step(0, ns - 1);
	ra = rastation(0, id_ns) - eflap;

	q = beta[0] - sita[2];
	dq = omega*sita[1];

	_setairfm(dfx, dfz, dfr, 0.0, q, dq, 0);

	qt = 2 * iflap* omega*(-omg[0] * cos(0) - omg[1] * sin(0)) + iflap*(domg[0] * sin(0) + domg[1] * cos(0));
	if (si_unit)
		qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8);
	else
		qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8*0.3048*0.3048);
	qt += (dfz*ra).sum()*radius;

	if (bld._GenArfStarter(qt, q, dq))
	{
		for (int i = 0; i < nitermax - 1; ++i)
		{
			t = i*dt + (1 - af)*dt;
			azmuth(i) = omega*t;
			niter = i + 1;
			_setairfm(dfx, dfz, dfr, t, q, dq, i);
			qt = 2 * iflap * omega*(-omg[0] * cos(omega*t) - omg[1] * sin(omega*t)) + iflap*(domg[0] * sin(omega*t) + domg[1] * cos(omega*t));
			if(si_unit)
				qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8);
			else
				qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8 * 0.3048 * 0.3048);

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

	/*if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _hingelessflap_rt(). Count: %d \n", niter);*/
	//bld.sol.output("sol.output", 10);
	euler_temp[0] = beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0;

	tppcoord.SetCoordinate(euler_temp, "euler");

	monitor.errb2 = bld.err_b;

	//printf("M Rotor Niter = %d, Errb2 = %e\n", niter, monitor.errb2);
	//printf("M Rotor flap (beta0, beta1c, beta1s) = (%f, %f, %f) \n", DEG(beta[0]), DEG(beta[1]), DEG(beta[2]));

}

bool Rotor::_hingelessflap_rt(double f[3], double m[3])
{
	int niter, nitermax, backstep, ck;
	myTYPE dt, af;
	myTYPE qt, t;
	myTYPE ia, ir, _up, _ut, _ma, _sfth, _infl, _lamh, _lami, _lamx, _lamy;
	myTYPE _dfx, _dfz, _factor, _dd;
	myTYPE _dak, _df, err;
	Matrix1<myTYPE> q, dq;

	Matrix1<myTYPE> dfx(ns), dfr(ns), dfz(ns), ra(ns), az(nf);
	myTYPE euler_temp[3];
	Matrix1<myTYPE> azmuth(bld.nitermax);
	Matrix1<int> id_ns = step(0, ns - 1);
	Matrix1<int> id_nf = step(0, nf - 1);
	Matrix2<myTYPE> aoad(nf, ns);
	myTYPE powerTemp, poweriTemp, poweroTemp, powercTemp, powerfTemp;
	bool _flg = false;
	Matrix1<int> np(ns);
	np.setvalue(1);

	f[0] = f[1] = f[2] = 0;
	m[0] = m[1] = m[2] = 0;
	powerTemp = power;
	poweriTemp = power_i;
	poweroTemp = power_o;
	powercTemp = power_c;
	powerfTemp = power_f;
	
	if (adyna == Averaged)
	{
		power_iid = 0;
		torque_iid = 0;
	}

	ra = rastation(0, id_ns);
	az = azstation(id_nf, 0);

	niter = 0, nitermax = bld.nitermax, af = bld.GAf.af;
	backstep = 5;
	dt = RAD(bld.dff) / omega;
	_df = RAD(bld.dff);
	q.allocate(nitermax);
	dq.allocate(nitermax);

	q(0) = beta[0] - sita[2];
	dq(0) = omega*sita[1];

	_setairfm(dfx, dfz, dfr, 0.0, q(0), dq(0), 0);
	qt = 2 * iflap* omega*(-omg[0] * cos(0) - omg[1] * sin(0)) + iflap*(domg[0] * sin(0) + domg[1] * cos(0));
	if (si_unit)
		qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8);
	else
		qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8*0.3048*0.3048);
	qt += (dfz*ra).sum()*radius;

	if (bld._GenArfStarter(qt, q(0), dq(0), backstep + 1))
	{
		for (int j = 0; j < ns; j++)
			lbstall[j].Starter(chord(j), amb.vsound, _df / omega);
		for (int i = 0; i < nitermax - 1; ++i)
		{
			int iz0, iz1, iz2;
			t = i*dt + (1 - af)*dt;
			azmuth(i) = omega*t;
			niter = i + 1;
			//cout << niter << endl;
			iz0 = i;
			while (iz0 > nf - 1)
				iz0 -= nf;
			iz1 = iz0 + 1;
			iz2 = iz0 - 1;
			if (iz1 > nf - 1)
				iz1 = 0;
			if (iz2 < 0)
				iz2 = nf - 1;

			// 计算气动挥舞力矩
			qt = 0;
			for (int j = 0; j < ns; ++j)
			{
				//printf("******** %d\n", j);
				ir = rastation(0, j);
				ia = _limitaz(omega*t);
				_lamh = lambdh.interplinear_fast(az, ra, ia, ir);
				_lamx = lambdx.interplinear_fast(az, ra, ia, ir);
				_lamy = lambdy.interplinear_fast(az, ra, ia, ir);

				_sfth = sita[0] + sita[1] * cos(ia) + sita[2] * sin(ia) + pitchroot - tan(RAD(del))*q(niter - 1) + twist(j);

				vel[0] -= _lamx*vtipa;
				vel[1] -= _lamy*vtipa;
				_setbladevelc(_ut, _up, q(niter-1), dq(niter-1), _lamh, ia, ir);
				vel[0] += _lamx*vtipa;
				vel[1] += _lamy*vtipa;

				_ut *= cos(sweep(j));
				_ma = rootNewton(_ut*_ut + _up*_up, _ut, 1e-6);
				_ma *= (vtipa / amb.vsound);
				_ma = Max(_ma, 1e-3);
				ua(iz0, j) = Min(_ma, 0.999);

				inflow(iz0, j) = atan2(_up, _ut);
				incidn(iz0, j) = _limitaoa(inflow(iz0, j) + _sfth);

				_dak = incidn(iz1, j) - incidn(iz2, j);
				if (_dak < -1.7*PI)
					_dak += 2 * PI;
				else if (_dak > 1.7*PI)
					_dak = _dak - 2 * PI;
				aoad(iz0, j) = _dak / (2 * _df / omega);
			
				lbstall[j].aoa0M1(iz0) = incidn(iz0, j);
				lbstall[j].MaM1(iz0) = ua(iz0, j);
				lbstall[j].q0M1(iz0) = aoad(iz0, j)*lbstall[j].Tl / ua(iz0, j);

				if (niter > backstep)
				{					
					lbstall[j].AttachFlow(iz0);
					lbstall[j].DynamicStall();
					lbstall[j].Complete();

					//lbstall[j].Display(j);

					cl(iz0, j) = lbstall[j].cl;
					cd(iz0, j) = lbstall[j].cd;
					//if (cd(iz0, j) < 0
					//cl(iz0, j) = cltc.interplinear_fast(cltc(step(0, cltc.NI - 1), 0), cltc(0, step(0, cltc.NJ - 1)), incidn(iz0, j) / PI*180.0, ua(iz0, j));
					//cd(iz0, j) = cdtc.interplinear_fast(cdtc(step(0, cdtc.NI - 1), 0), cdtc(0, step(0, cdtc.NJ - 1)), incidn(iz0, j) / PI*180.0, ua(iz0, j));
					//cd(iz0, j) = 0.0106;

					_factor = _setairfm(_dfx, _dfz, cl(iz0, j), cd(iz0, j), iz0, j); 
					// _setairfm(_dfx, _dfz, _dfr, cl(iz0, j), cd(iz0, j), iz0, j);

					qt += _dfz*(ir - eflap)*radius;	
					
				}			
			}

			if (niter > backstep)
			{
				qt += 2 * iflap * omega*(-omg[0] * cos(omega*t) - omg[1] * sin(omega*t)) + iflap*(domg[0] * sin(omega*t) + domg[1] * cos(omega*t));
				if (si_unit)
					qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8);
				else
					qt += m1*(dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0] - 9.8 * 0.3048 * 0.3048);

				bld._GenArfTimeMarch(qt, dt, niter);
				q(niter) = bld.q;
				dq(niter) = bld.dq;
			}

			// 计算桨盘平均气动力			
		
			if (bld.isGenArfExit(niter))
			{			
				err = 0;
				for (int j = 0; j < ns; j++)
				{
					lbstall[j].isExit(ck, np(j));
					err += lbstall[j].err;
				}
				if (err < ns*lbstall[0].eps)
				{
					double b, db;
					//double dqi, dqo, dm;
					//dqi = dqo = dm = 0;

					power = power_c = power_i = power_f = power_o = 0;
					torque = torque_c = torque_i = torque_f = torque_o = 0;
					
					for (int j = 0; j < ns; j++)
					{
						int ic = niter - 1;
						ir = rastation(0, j);

						for (int ip = 0; ip < np(j); ip++)
						{
							for (int i = iz0; i >= 0; i--, ic--)
							{
								ia = _limitaz(i * 2 * PI / nf);
								b = q(ic);// beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
								db = dq(ic);// (-beta[1] * sin(ia) + beta[2] * cos(ia))*omega;
								_factor = _setairfm(_dfx, _dfz, cl(i, j), cd(i, j), i, j);
								//_setairfm(_dfx, _dfz, _dfr, cl(i, j), cd(i, j), i, j);
								f[0] += (_dfx*sin(ia) - _dfz*sin(b)*cos(ia)) / np(j);
								f[1] -= (_dfx*cos(ia) + _dfz*sin(b)*sin(ia)) / np(j);
								f[2] += _dfz*cos(b) / np(j);

								m[2] -= (_dfx*((ir - eflap)*cos(b) + eflap)) / np(j);
								
								// power computation
								_lami = lambdi.interplinear_fast(az, ra, ia, ir);
								_dd = cd(i, j)*_factor / (cos(inflow(i, j)) +1e-6);

								power_i += (db / omega*(ir - eflap)*_dfz - _lami*_dfz*cos(b)) / np(j);
								//if(cd(i,j)>0)
								power_o += (_dd*sin(ia)*mul + _dd*((ir - eflap)*cos(b) + eflap)) / np(j);

								//double rcap = ((ir - eflap)*cos(b) + eflap);
								//dm += _dfx*rcap;
								//dqi += tan(inflow(i, j))*_dfz*rcap;
								//dqo += _dd*rcap;
								//printf("dQ = %f, dQi = %f, dQo = %f, dQ+dQi-dQo = %e \n", dm, dqi, dqo, dm+dqi-dqo);

								if (adyna == Averaged)
									power_iid += db / omega*(ir - eflap)*_dfz / np(j);
							}

							for (int i = nf - 1; i >= iz0 + 1; i--, ic--)
							{
								ia = _limitaz(i * 2 * PI / nf);
								b = q(ic); // beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
								db = dq(ic);// (-beta[1] * sin(ia) + beta[2] * cos(ia))*omega;
								_factor = _setairfm(_dfx, _dfz, cl(i, j), cd(i, j), i, j);
								//_setairfm(_dfx, _dfz, _dfr, cl(i, j), cd(i, j), i, j);
								f[0] += (_dfx*sin(ia) - _dfz*sin(b)*cos(ia)) / np(j);
								f[1] -= (_dfx*cos(ia) + _dfz*sin(b)*sin(ia)) / np(j);
								f[2] += _dfz*cos(b) / np(j);

								m[2] -= (_dfx*((ir - eflap)*cos(b) + eflap)) / np(j);

								// power computation
								_lami = lambdi.interplinear_fast(az, ra, ia, ir);
								_dd = cd(i, j)*_factor / (cos(inflow(i, j)) +1e-6);

								power_i += (db / omega*(ir - eflap)*_dfz - _lami*_dfz*cos(b)) / np(j);
								//if(cd(i,j)>0)
								power_o += (_dd*sin(ia)*mul + _dd*((ir - eflap)*cos(b) + eflap)) / np(j);

								//double rcap = ((ir - eflap)*cos(b) + eflap);
								//dm += _dfx*rcap;
								//dqi += tan(inflow(i, j))*_dfz*rcap;
								//dqo += _dd*rcap;
								//printf("dQ = %f, dQi = %f, dQo = %f, dQ+dQi-dQo = %e \n", dm, dqi, dqo, dm+dqi-dqo);

								if (adyna == Averaged)
									power_iid += db / omega*(ir - eflap)*_dfz / np(j);
							}
						}
					}

					//q.output("q.output", 10);

					_flg = true;
					break;
				}
			}			
		}
		// 计算桨盘平均力
		f[0] *= (1.0 *nb / nf);
		f[1] *= (1.0* nb / nf);
		f[2] *= (1.0*nb / nf);
		m[2] *= (radius*nb / nf);
	}
	//bld.sol.output("sol.output", 10);

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
			printf("Undefine nb = %d in _hingelessflap_rt() \n", nb);
			system("pause");
			break;
		}
	}
	beta[0] /= bld.nperiod;
	beta[1] /= bld.nperiod;
	beta[2] /= bld.nperiod;

	if (niter == nitermax - 1)
	{
		;// printf("Warning: Flap solving may not be convergent in Func _hingelessflap_rt(). Count: %d \n", niter);
		
		//for (int i = 3; i < 6; i++)
		//{
		//	cout << lbstall[i].err << endl;
		//	lbstall[i].alphaIn.output("AOA" + std::to_string(i) + ".output", 10);
		//	lbstall[i].CNT.output("CNT" + std::to_string(i) + ".output", 10);
		//	lbstall[i].CL.output("CL" + std::to_string(i) + ".output", 10);
		//	lbstall[i].CD.output("CD" + std::to_string(i) + ".output", 10);
		//}
	}
	monitor.LBerrSum = err;

	euler_temp[0] = beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0;

	tppcoord.SetCoordinate(euler_temp, "euler");

	switch (hingetype)
	{
	case Teeter:
		m[0] = m[1] = 0;
		break;
	case Hinged:
		m[0] = -nb*khub / 2 * beta[2];
		m[1] = -nb*khub / 2 * beta[1];
	case Hingeless:
		m[0] = nb*khub *0.5  * beta[2];// -m[2] * beta[1] * 0.5;
		m[1] = -nb*khub *0.5 * beta[1];// -m[2] * beta[2] * 0.5;
		break;
	default:
		break;
	}
	/*cout << power_o << endl;
	printf("(%f, %f, %f) \n", beta[0], beta[1], beta[2]);
	cl.output("LF_MR_4_CL.output", 10);
	cd.output("LF_MR_4_CD.output", 10);
	inflow.output("LF_MR_4_inflow.output", 10);
	incidn.output("LF_MR_4_AOA.output", 10);
	ua.output("LF_MR_4_Ma.output", 10);
	lambdi.output("LF_MR_4_lambdi.output", 10);*/

	power_i *= vtipa*nb / nf;
	power_o *= vtipa*nb / nf;
	power_f = f[0] * vel[0];
	power_c = f[2] * vel[2];
	power = Abs(m[2]) * omega;
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
		power /= 1000.0;
		if (adyna == Averaged)
			power_iid /= 1000.0;
	}
	
	monitor.errb2 = bld.err_b;
	monitor.Countsb = niter;

	//printf("M Rotor Niter = %d, Errb2 = %e\n", niter, monitor.errb2);
	//printf("M Rotor flap (beta0, beta1c, beta1s) = (%f, %f, %f) \n", DEG(beta[0]), DEG(beta[1]), DEG(beta[2]));
	//printf("Power Comp: (%f, %f, %f, %f, %f) \n", power, power_i, power_o, power_f, power_c);
	if(!_flg)
	{
		power = powerTemp;
		power_i = poweriTemp;
		power_o = poweroTemp;
		power_c = powercTemp;
		power_f = powerfTemp;
	}
	return _flg;
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

	f22(1) = 2 * (1 + eflap*mb / iflap)*omgw[0] / omega + (1 + eflap*mb / iflap)*domgw[1] / omega / omega;
	f22(1) += gama / 2.0*((0.25 - eflap / 3.0) + mul*mul / 4.0*(0.5 - eflap + 0.5*eflap*eflap))*sitaw[1] + gama / 2.0*(0.25 - eflap / 3.0)*omgw[1] / omega;

	f22(2) = -2.0*(1 + eflap*mb / iflap)*omgw[1] / omega + (1 + eflap*mb / iflap)*domgw[0] / omega / omega;
	f22(2) += gama*mul / 2.0*(2.0 / 3.0 - eflap)*sitaw[0] + gama*mul / 2.0*(0.5 - 2.0 / 3.0*eflap)*twistt;
	f22(2) += gama / 2.0*((0.25 - eflap / 3.0) + 3.0 / 4.0*mul*mul*(0.5 - eflap + 0.5*eflap*eflap))*sitaw[2];
	f22(2) += gama*mul / 2.0*(0.5 - eflap + eflap*eflap / 2.0)*lambdh_ag + gama / 2.0*(0.25 - eflap / 3.0)*omgw[0] / omega;

	f22 = f22 * (omega*omega);

	if (!si_unit)
		f22(0) += m1 / iflap *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
	else
		f22(0) += m1 / iflap *(-9.8 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);

	_sum = 0;
	temp = 0;
	sb = 8.0*(p2 - 1.0) / gama;
	q0(1) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	q0(2) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb + 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);
	q0(0) = gama / p2*(sitaw[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sitaw[2]) + lambdh_ag / 6.0);

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
			if (_sum < bld.epsb*bld.epsb)
				break;
		}

	}
	//sol.output("_sol_temp_mr.output", 10);
	if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _flapWang(). \n");

	beta[0] = betaw[0];
	beta[1] = betaw[2] * sin(windcoord.euler[2]) + betaw[1] * cos(windcoord.euler[2]);
	beta[2] = betaw[2] * cos(windcoord.euler[2]) - betaw[1] * sin(windcoord.euler[2]);

	euler_temp[0] = beta[2]; // -beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0.0;
	//euler_temp[2] = Atan2(vel[1], vel[0]); // TPP wind coordinate
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
	q0(0) = gama / p2*(sitaw[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sitaw[2]) + lambdh_ag / 6.0);

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
			if (_sum < bld.epsb*bld.epsb)
				break;
		}

	}
	//sol.output("_sol_temp_mr.output", 10);
	if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _hingelessflap_fx(). \n");

	
	beta[0] = betaw[0];
	beta[1] = betaw[2] * sin(windcoord.euler[2]) + betaw[1] * cos(windcoord.euler[2]);
	beta[2] = betaw[2] * cos(windcoord.euler[2]) - betaw[1] * sin(windcoord.euler[2]);

	euler_temp[0] = beta[2]; // -beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0.0;
	//euler_temp[2] = Atan2(vel[1], vel[0]); // 桨盘平面旋转，是因为Glauert入流模型是没有横向速度的
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
	sitaw[1] = sita[1] *cos(betawind) - sita[2] * sin(betawind);
	sitaw[2] = sita[1] *sin(betawind) + sita[2] * cos(betawind);
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
	HH(0) += 4.0 / 3.0*mul*sitaw[2] + 4.0 / 3.0*lambdh_ag + 2.0 / 3.0*mul*(omgw[0] / omega + lamb_1s);
	
	HH(2) = 16 / gama*(-omgw[0] / omega + domgw[1] / omega / omega / 2.0) + sitaw[1] * (1 + mul*mul / 2.0) + (omgw[1] / omega + lamb_1c);
	HH(3) = -16 / gama*(omgw[1] / omega - domgw[0] / omega / omega / 2.0) + 8.0 / 3.0*mul*sitaw[0] + 2 * mul*twistt;
	HH(3) += sitaw[2] * (1 + 3.0 / 2.0*mul*mul) + 2 * mul*lambdh_ag + (omgw[0] / omega + lamb_1s);

	CC = CC * omega;
	DD = DD * omega*omega;
	HH = HH*omega*omega*gama / 8.0;

	if (!si_unit)
		HH(0) += m1 / iflap *(-9.8*0.3048*0.3048 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);
	else
		HH(0) += m1 / iflap *(-9.8 + dvelh[2] - velh[0] * omgh[1] + velh[1] * omgh[0]);

	_sum = 0;
	temp = 0;
	sb = 8.0*(p2 - 1.0) / gama;
	q0(1) = 0.0;
	q0(2) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[1] - sitaw[2] + (-sb*16.0 / gama - 1)*omgw[0] + (sb + 16.0 / gama)*omgw[1]);
	q0(3) = 1.0 / (1.0 + sb*sb)*(sb*sitaw[2] + sitaw[1] + (sb - 16.0 / gama)*omgw[0] - (sb*16.0 / gama - 1)*omgw[1]);
	q0(0) = gama / p2*(sitaw[0] * 0.125*(1.0 + mul*mul) + 0.1*twistt*(1 + 5.0 / 6.0*mul*mul) + 1.0 / 6.0*mul*(q0(1) + sitaw[2]) + lambdh_ag / 6.0);

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
			if (_sum < bld.epsb*bld.epsb)
				break;
		}
	}
	/*if (niter == nitermax - 1)
		printf("Warning: Flap solving may not be convergent in Func _flapMBC(). \n");*/

	beta[0] = betaw[0];
	beta[1] = betaw[2]*sin(betawind) + betaw[1] * cos(betawind);
	beta[2] = betaw[2]*cos(betawind) - betaw[1] * sin(betawind);

	euler_temp[0] = beta[2];
	euler_temp[1] = beta[1];
	euler_temp[2] = 0.0;
	//euler_temp[2] = Atan2(vel[1], vel[0]); // 桨盘平面旋转，是因为Glauert入流模型是没有横向速度的
	tppcoord.SetCoordinate(euler_temp, "euler");

	Errb2 = _sum;
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
	
	if (bld.soltype == HubFixed)
	{
		if (adyna == Simple)
		{
			ct = sigma*a0*0.5*(sita[0] / 3.0*(1 + 1.5*mul*mul) + 0.25*twistt*(1 + mul*mul) + 0.5*mul*beta[1] + 0.5*lambtpp);
			airforce[2] = ct*amb.rho*PI*radius*radius*vtipa*vtipa;
		}
		else
		{
			_setairfm_sp(_af_temp, _am_temp);
			for (int i = 0; i < 3; i++)
			{
				airforce[i] = _af_temp[i];
				airmoment[i] = _am_temp[i];
			}
		}
	}
	else
	{
		if (airfoil != LBStallMethod)
		{
			if (adyna == Simple)
			{
				ct = sigma*a0*0.5*(sita[0] / 3.0*(1 + 1.5*mul*mul) + 0.25*twistt*(1 + mul*mul) + 0.5*mul*beta[1] + 0.5*lambtpp);
				airforce[2] = ct*amb.rho*PI*radius*radius*vtipa*vtipa;
			}
			else
			{
				_setairfm_sp(_af_temp, _am_temp);
				for (int i = 0; i < 3; i++)
				{
					airforce[i] = _af_temp[i];
					airmoment[i] = _am_temp[i];
				}
			}
		}
	}

	switch (adyna)
	{
	case PWake:
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;
		//lambdi_ag = -airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);

		if (haveGeo && haveStr)
		{
			_tipVortexStr();
			//_wakeGeoBd();
			_bladePosition();
			_wakeInducedVelMP(nb);
		}
		else
		{
			lambdh_ag = lambdi_ag - vel[2] / vtipa;
			lambdh.setvalue(lambdh_ag);
			lambdi.setvalue(lambdi_ag);
		}
		break;
	case LinearInflow:
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;
		lambdh_ag = lambdi_ag - vel[2] / vtipa;

		ka = atan2(-veltpp[0] / vtipa, -(-veltpp[2] / vtipa + lambdi_ag));
		ky = 0* 2 * veltpp[0] / vtipa;
		if (ka <= PI / 2 && ka >= -PI / 2)
			kx = tan(ka * 0.5);
		else if (ka<PI && ka>-PI)
			kx = 1 / tan(ka * 0.5);
		else
			kx = 0.0;

		// non-uniform inflow model is doing in TPP wind coordinate
		Iid_betaw = betawind / (2 * PI / nf);
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
				lambdh(i, j) = lambdi_temp - vel[2] / vtipa; 
			}
		}
		break;
	case Averaged:
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;
		lambdh_ag = lambdi_ag - vel[2] / vtipa;

		lambdh.setvalue(lambdh_ag);
		lambdi.setvalue(lambdi_ag);
		break;
	case Simple:
		// simple ct
		lambdi_ag -= ct / 2.0 / sqrt(mul*mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;
		lambdh_ag = lambdi_ag - vel[2] / vtipa;

		lambdh.setvalue(lambdh_ag);
		lambdi.setvalue(lambdi_ag);
		break;
	default:
		// airforce at hub coord approx to tpp coord
		lambdi_ag -= airforce[2] / (2 * amb.rho*PI*radius*radius*vtipa*vtipa) / sqrt(mul * mul + lambtpp * lambtpp);
		lambdi_ag /= 2.0;

		lambdh_ag = lambdi_ag - vel[2] / vtipa;
		lambdh.setvalue(lambdh_ag);
		lambdi.setvalue(lambdi_ag);
		break;
	}
	KA = atan2(-veltpp[0] / vtipa, -(-veltpp[2] / vtipa + lambdi_ag));
	monitor.KA = KA;

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

	w[2] = omg[2];
	dw[2] = domg[2];
	w[0] = omg[0] * cos(windcoord.euler[2]) + omg[1] * sin(windcoord.euler[2]);
	dw[0] = domg[0] * cos(windcoord.euler[2]) + domg[1] * sin(windcoord.euler[2]);
	w[1] = -omg[0] * sin(windcoord.euler[2]) + omg[1] * cos(windcoord.euler[2]);
	dw[1] = -domg[0] * sin(windcoord.euler[2]) + domg[1] * cos(windcoord.euler[2]);

	/*for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			w[i] += windcoord.Ttransf[i][j] * omg[j];
			dw[i] += windcoord.Ttransf[i][j] * domg[j];
		}
	}*/
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
	double _df = azstation(1, 0) - azstation(0, 0);
	Matrix1<double> _dt(ns), _yf(ns), _hf(ns), _az(ns), _lambdi(ns), az(nf);
	Matrix1<double> _dfx(ns), _dfz(ns), _dfr(ns), ra(ns), ra2(ns), ra1(ns);
	Matrix1<double> _dD(ns);
	//Matrix1<double> _aoa(nf), _ma(nf);
	Matrix1<int> id_ns = step(0, ns - 1);
	Matrix1<int> id_nf = step(0, nf - 1);
	ra = rastation(0, id_ns);
	az = azstation(id_nf, 0);

	f[0] = f[1] = f[2] = 0.0;
	m[0] = m[1] = m[2] = 0.0;

	power = power_c = power_i = power_f = power_o = 0;
	torque = torque_c = torque_i = torque_f = torque_o = 0;
	
	if (adyna == Averaged)
	{
		power_iid = 0;
		torque_iid = 0;
	}

	if (airfoil != LBStallMethod)
	{
		for (int i = 0; i < nf; ++i) {
			ia = i * 2 * PI / nf;
			it = ia / omega;
			b = beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
			db = (-beta[1] * sin(ia) + beta[2] * cos(ia)) * omega;

			_setairfm(_dfx, _dfz, _dfr, it, b, db, i);

			_dt = _dfz*cos(b);
			_yf = _dfx*cos(ia) - _dfr*sin(ia);
			_hf = _dfx*sin(ia) + _dfr*cos(ia);

			f[0] += _hf.sum()*nb / nf;
			f[1] -= _yf.sum()*nb / nf;
			f[2] += _dt.sum()*nb / nf;

			//cout << f[2] << endl;

			ra2 = (ra - eflap)*sin(b);
			ra1 = (ra - eflap)*cos(b) + eflap;

			m[2] -= (_dfx * ra1).sum()*radius*nb / nf;

			// power computation
			_az.setvalue(ia);
			_lambdi = lambdi.interplinear_fast(az, ra, _az, ra);
			_dD = _cd*_factor / mcos(_inflow);
			power_i += db / omega * ((ra - eflap)*_dfz).sum() - (_dt*_lambdi).sum();
			power_o += (_dD*sin(ia)*mul + _dD*ra1).sum();

			if (adyna == Averaged)
				power_iid += db / omega*((ra - eflap)*_dfz).sum();


			// save disk plane data
			for (int j = 0; j < ns; j++)
			{
				cl.v_p[j*cl.NI + i] = _cl(j);
				cd.v_p[j*cd.NI + i] = _cd(j);
				inflow.v_p[j*inflow.NI + i] = _inflow(j);
				ua.v_p[j*ua.NI + i] = _ua(j);
				dt.v_p[j*dt.NI + i] = _dt(j);
				incidn.v_p[j*incidn.NI + i] = _incidn(j);
				
				//if (_incidn(j) > PI / 2 && _incidn(j) <= PI)
				//	incidn.v_p[j*incidn.NI + i] = PI - _incidn(j);
				//else if (_incidn(j) < -PI / 2 && _incidn(j) >= -PI)
				//	incidn.v_p[j*incidn.NI + i] = -PI - _incidn(j);
				//else
				incidn.v_p[j*incidn.NI + i] = _incidn(j);
			}
		}

		//cout << power_o << endl;

	}
	else
	{
		double ir, _up, _ut, _infl, _lamh, _lami;
		int ck = 0;
		double _dak;
		Matrix2<double> aoad(nf, ns);
		//Matrix1<double> _cllbs(nf), _cdlbs(nf), _cnlbs(nf), _cclbs(nf);
		double _cl, _cd;
		double _dfx, _dfz, _dfr, _factor, _dd;


		for (int j = 0; j < ns; ++j)
		{
			ir = rastation(0, j);

			for (int i = 0; i < nf; ++i)
			{
				ia = _limitaz(i * 2 * PI / nf);
				_lamh = lambdh.interplinear_fast(az, ra, ia, ir);
				_setbladeaeros(ua(i, j), incidn(i, j), inflow(i, j), ia, ir, _lamh, i, j);
			}
			for (int i = 1; i < nf - 1; ++i)
			{
				_dak = incidn(i + 1, j) - incidn(i - 1, j);
				if (_dak < -1.7*PI)
					_dak += 2 * PI;
				else if (_dak > 1.7*PI)
					_dak = _dak - 2 * PI;
				aoad(i, j) = _dak / (2 * _df / omega);
			}
			_dak = (incidn(1, j) - incidn(nf - 1, j));
			if (_dak < -1.7*PI)
				_dak += 2 * PI;
			else if (_dak > 1.7*PI)
				_dak = _dak - 2 * PI;
			aoad(0, j) = _dak / (2 * _df / omega);
			_dak = (incidn(0, j) - incidn(nf - 2, j));
			aoad(nf - 1, j) = _dak / (2 * _df / omega);

			lbstall[0].Starter(chord(j), amb.vsound, _df / omega);
			if (lbstall[0].Solver(ck, j, incidn, ua, aoad))
			{			
				int iz0 = ck / nf;
				int ic = ck;
				iz0 = ck - iz0*nf - 1;
				if (iz0 < 0)
					iz0 += nf;
				for (int i = iz0; i >= 0; i--, ic--)
				{
					ia = _limitaz(i * 2 * PI / nf);
					b = beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
					db = (-beta[1] * sin(ia) + beta[2] * cos(ia)) * omega;

					_cl = lbstall[0].CL(ic);
					_cd = lbstall[0].CD(ic);
					_factor = _setairfm(_dfx, _dfz, _dfr, _cl, _cd, i, j);

					f[0] += (_dfx*sin(ia) + _dfr*cos(ia)) * nb / nf;
					f[1] -= (_dfx*cos(ia) - _dfr*sin(ia)) * nb / nf;
					f[2] += _dfz*cos(b) * nb / nf;

					m[2] -= (_dfx*((ir - eflap)*cos(b) + eflap)) * radius * nb / nf;

					// power computation
					_lami = lambdi.interplinear_fast(az, ra, ia, ir);
					_dd = _cd*_factor / (cos(inflow(i, j)) + 1e-6);

					power_i += db / omega*(ir - eflap)*_dfz - _lami*_dfz*cos(b);
					power_o += _dd*sin(ia)*mul + _dd*((ir - eflap)*cos(b) + eflap);
					if (adyna == Averaged)
						power_iid += db / omega*(ir - eflap)*_dfz;
				}
				for (int i = nf - 1; i >= iz0 + 1; i--, ic--)
				{
					ia = _limitaz(i * 2 * PI / nf);
					b = beta[0] + beta[1] * cos(ia) + beta[2] * sin(ia);
					db = (-beta[1] * sin(ia) + beta[2] * cos(ia)) * omega;

					_cl = lbstall[0].CL(ic);
					_cd = lbstall[0].CD(ic);
					_factor = _setairfm(_dfx, _dfz, _dfr, _cl, _cd, i, j);

					f[0] += (_dfx*sin(ia) + _dfr*cos(ia)) * nb / nf;
					f[1] -= (_dfx*cos(ia) - _dfr*sin(ia)) * nb / nf;
					f[2] += _dfz*cos(b) * nb / nf;

					m[2] -= (_dfx*(ir - eflap)*cos(b) + eflap) * radius * nb / nf;

					// power computation
					_lami = lambdi.interplinear_fast(az, ra, ia, ir);
					_dd = _cd*_factor / (cos(inflow(i, j)) + 1e-6);

					power_i += db / omega*(ir - eflap)*_dfz - _lami*_dfz*cos(b);
					power_o += _dd*sin(ia)*mul + _dd*((ir - eflap)*cos(b) + eflap);
					if (adyna == Averaged)
						power_iid += db / omega*(ir - eflap)*_dfz;
				}

				lbstall[0].Save(j);
			}
			else
				;// printf("Wrong: LB Stall Solver is not convergent. \n");
		}
		airfoil = C81Table;
	}
	switch (hingetype)
	{
	case Teeter:
		m[0] = m[1] = 0;
		break;
	case Hinged:
		m[0] = -nb*khub / 2 * beta[2];
		m[1] = -nb*khub / 2 * beta[1];
	case Hingeless:
		m[0] = nb*khub *0.5  * beta[2];// -m[2] * beta[1] * 0.5;
		m[1] = -nb*khub *0.5 * beta[1];// -m[2] * beta[2] * 0.5;
		break;
	default:
		break;
	}
	//lambdi.output("temp_lambdi.output", 10);
	power_i *= vtipa*nb / nf;
	power_o *= vtipa*nb / nf;
	power_f = f[0] * vel[0];
	power_c = f[2] * vel[2];
	power = Abs(m[2]) * omega;
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
		power /= 1000.0;
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
	if (ia > 2 * PI || ia < 0) {
		ia = omega*it;
		for (;;) {
			if (ia > 2 * PI) { ia -= 2 * PI; }
			else if (ia < 0) { ia += 2 * PI; }
			else { break; }
		}
	}

	_sfth = sita[0] + sita[1] * cos(ia) + sita[2] * sin(ia) + pitchroot - tan(RAD(del))*b;

	// aoa
	_inflow = atan2(_up, _ut);
	_incidn = twist + _inflow + _sfth;
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
	_dr(0) = 0.5*_dr(0);
	_dr(ns - 1) = 0.5*_dr(ns - 1);
	_factor *= chord * _dr *amb.vsound*amb.vsound*0.5*amb.rho* radius;

	_dfx = (_cd * mcos(_inflow) - _cl * msin(_inflow)) * _factor;
	_dfz = (_cl * mcos(_inflow) + _cd * msin(_inflow)) * _factor;
	_dfr = _dfz * (sin(b)) * (-1.0);
}

double Rotor::_setairfm(double &_dfx, double &_dfz, double &_dfr, const double &cl, const double &cd, const int iz, const int ir)
{
	double _ua = ua(iz, ir);
	double _ua2 = _ua*_ua;
	double _infl = inflow(iz, ir);
	double _c = chord(ir);
	double _az = _limitaz(iz * 2 * PI / nf);
	double b = beta[0] + beta[1] * cos(_az) + beta[2] * sin(_az);
	double _dr, _factor;

	cirlb(iz, ir) = 0.5*ua(iz, ir)*cl*_c*amb.vsound;
	if (ir > 0 && ir < ns - 1)
		_dr = rastation(0, ir) - rastation(0, ir - 1);
	else if (ir == 0)
		_dr = 0.5*(rastation(0, 1) - rastation(0, 0));
	else
		_dr = 0.5*(rastation(0, ns - 1) - rastation(0, ns - 2));

	_factor = _ua2*_c*_dr*amb.vsound*amb.vsound*0.5*amb.rho*radius;
	_dfx = (cd*cos(_infl) - cl*sin(_infl))*_factor;
	_dfz = (cl*cos(_infl) + cd*sin(_infl))*_factor;
	_dfr = -_dfz*sin(b);
	return _factor;
}

double Rotor::_setairfm(double &_dfx, double &_dfz, const double &cl, const double &cd, const int iz, const int ir)
{
	double _ua = ua(iz, ir);
	double _ua2 = _ua*_ua;
	double _infl = inflow(iz, ir);
	double _c = chord(ir);
	double _az = _limitaz(iz * 2 * PI / nf);
	double _dr, _factor;

	cirlb(iz, ir) = 0.5*ua(iz, ir)*cl*_c*amb.vsound;
	if (ir > 0 && ir < ns - 1)
		_dr = rastation(0, ir) - rastation(0, ir - 1);
	else if (ir == 0)
		_dr = 0.5*(rastation(0, 1) - rastation(0, 0));
	else
		_dr = 0.5*(rastation(0, ns - 1) - rastation(0, ns - 2));
	
	_factor = _ua2*_c*_dr*amb.vsound*amb.vsound*0.5*amb.rho*radius;
	_dfx = (cd*cos(_infl) - cl*sin(_infl))*_factor;
	_dfz = (cl*cos(_infl) + cd*sin(_infl))*_factor;

	return _factor;
}


void Rotor::_setbladeaeros(double &_ma, double &_aoa, double &_infl, const double &_az, const double &ra, const double &lamh, const int iz, const int ir)
{
	double _ut, _up, b, db, _sfth;

	b = beta[0] + beta[1] * cos(_az) + beta[2] * sin(_az);
	db = (-beta[1] * sin(_az) + beta[2] * cos(_az)) * omega;
	_sfth = sita[0] + sita[1] * cos(_az) + sita[2] * sin(_az) + pitchroot - tan(RAD(del))*b + twist(ir);

	_setbladevelc(_ut, _up, b, db, lamh, _az, ra);
	_ut *= cos(sweep(ir));

	_ma = rootNewton(_ut*_ut + _up*_up, _ut, 1e-6);
	_ma *= (vtipa / amb.vsound);
	_ma = Max(_ma, 0);
	_ma = Min(_ma, 1);

	// aoa
	_infl = atan2(_up, _ut);
	_aoa = _infl + _sfth;
	_aoa = _limitaoa(_aoa);
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
		else if (ia < 0) ia += 2 * PI;
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
	_ut = (_ra - eflap) * cos(b) + eflap - (sin(ia) * vel[0] - cos(ia) * vel[1]) / vtipa;
	_up = (_ra - eflap) * (-1.0) * db / omega + (cos(ia) * vel[0] + sin(ia) * vel[1]) * sin(b) / vtipa + _lambdh * cos(b);
	//_ut += _ra*(omg[0] / omega*cos(ia) - omg[1] / omega*sin(ia)) * b;
	_up += _ra*(-omg[0] / omega*sin(ia) + omg[1] / omega*cos(ia)); 

	//_ut = (_ra - eflap) * cos(b) + eflap + sin(ia) * vel[0] / vtipa;
	//_up = (_ra - eflap) * db / omega + cos(ia) * vel[0] * sin(b) / vtipa + _lambdh * cos(b);
}

void Rotor::_setbladevelc(double & _ut, double &_up, const double &b, const double &db, const double &lamh, const double &ia, const double &ir)
{
	_ut = (ir - eflap)*cos(b) + eflap - (sin(ia)*vel[0] - cos(ia)*vel[1]) / vtipa;
	_up = (ir - eflap)*(-1)*db / omega + (cos(ia)*vel[0] + sin(ia)*vel[1])*sin(b) / vtipa + lamh*cos(b);
	_up += ir * (-omg[0] / omega*sin(ia) + omg[1] / omega*cos(ia));
}


double Rotor::_limitaz(double ia)
{
	for (;;)
	{
		if (ia > 2 * PI)
			ia -= 2 * PI;
		else if (ia < 0)
			ia += 2 * PI;
		else
			break;
	}
	if (ia > azstation(nf - 1, 0))
		ia = 0;
	return ia;
}

double Rotor::_limitaoa(double _aoa)
{
	for (;;) {
		if (_aoa > PI) { _aoa -= 2 * PI; }
		else if (_aoa < -PI) { _aoa += 2 * PI; }
		else { break; }
	}
	return _aoa;
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
	//if (lbstall.enable)
	//{
	//	lbstall.CNTM2.output(s + "_CNT.output", 10);
	//	lbstall.CCfM2.output(s + "_CCf.output", 10);
	//	lbstall.CDM2.output(s + "_CD.output", 10);
	//	lbstall.CLM2.output(s + "_CL.output", 10);
	//	lbstall.FlowStateSaveM2.output(s + "_FlowState.output", 2);
	//	lbstall.VortexStateSaveM2.output(s + "_VortexState.output", 2);
	//	lbstall.fppM2.output(s + "_fpp.output", 10);
	//	lbstall.Circle.output(s + "_Circle.output", 2);
	//}
	//else
	//{
	//	cl.output(s + "_CL.output", 10);
	//	cd.output(s + "_CD.output", 10);
	//}
	cl.output(s + "_CL.output", 10);
	cd.output(s + "_CD.output", 10);
	inflow.output(s + "_inflow.output", 10);
	incidn.output(s + "_AOA.output", 10);
	ua.output(s + "_Ma.output", 10);
	lambdi.output(s + "_lambdi.output", 10);
	dt.output(s + "_dt.output", 10);
}

void Rotor::GetStates(double v[3], double w[3])
{
	for (int i = 0; i < 3; i++)
	{
		v[i] = vel[i];
		w[i] = omg[i];
	}
}

void Rotor::GetWStates(double v[3], double w[3])
{
	for (int i = 0; i < 3; i++)
	{
		v[i] = velw[i];
		w[i] = omgw[i];
	}
}

bool BladeSolver::isGenArfExit(int niter)
{
	myTYPE temp;
	if (niter >= 2 * nperiod - 1)
	{
		err_b = 0;
		for (int j = 0; j < nperiod; ++j)
		{
			temp = sol(niter - nperiod + 1 + j) - sol(niter - 2 * nperiod + 1 + j);
			err_b += temp*temp;
		}
		if (err_b < epsb*epsb && niter%nperiod==0)
			return true; 
	}
	return false;
}

