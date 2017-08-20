#include "stdafx.h"
#include "Component.h"

Copter::Copter() :System()
{
	print_cons_on_screen("Copter");
	/* read config files */
	mass = 915;
	omega = 525 * PI / 30;
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
	myTYPE origin[3] = { 0,0,0 };
	myTYPE euler[3] = { 0,0,0 };
#ifdef TEST_DATA
	euler[1] = 0.78*PI / 180; // positive means nose down in cpp
#endif // TEST_DATA

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
	vel_g[1] = 0 * 632.2455;
	vel_g[2] = 0 * 632.2455;
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
	inmatx_M = H.inmatx_M;

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
	inmatx_M.deallocate();

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
	vel_c[0] = vel_c[1] = vel_c[2] = 0;
	omg_c[0] = omg_c[1] = omg_c[2] = 0;
	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			vel_c[i] += refcoord.Ttransf[i][j] * vel_g[j];
			omg_c[i] += refcoord.Etransf[i][j] * omg_g[j];
		}
	}

	//#ifdef USE_DOUBLE
	//	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, vel_c, 1);
	//	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, omg_c, 1);
	//#else
	//	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, vel_c, 1);
	//	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, omg_c, 1);
	//#endif // USE_DOUBLE

#ifdef FLIGHT_TRIM
	dvel_c[0] = dvel_c[1] = dvel_c[2] = 0;
	domg_c[0] = domg_c[1] = domg_c[2] = 0;
#else
	;
#endif // FLIGHT_TRIM

#ifdef TEST_MODE
	//printf("Copter states.\n");
	//printf("Velocity: %f, %f, %f \n", vel_c[0], vel_c[1], vel_c[2]);
	//printf("Omega: %f, %f, %f \n", omg_c[0], omg_c[1], omg_c[2]);
	//printf("dvelocity: %f, %f, %f \n", dvel_c[0], dvel_c[1], dvel_c[2]);
	//printf("domega: %f, %f, %f \n", domg_c[0], domg_c[1], domg_c[2]);
	//cout << vel_c[0] * vel_c[0] + vel_c[1] * vel_c[1] + vel_c[2] * vel_c[2] << endl;
	//cout << omg_c[0] * omg_c[0] + omg_c[1] * omg_c[1] + omg_c[2] * omg_c[2] << endl;
	//cout << dvel_c[0] * dvel_c[0] + dvel_c[1] * dvel_c[1] + dvel_c[2] * dvel_c[2] << endl;
	//cout << domg_c[0] * domg_c[0] + domg_c[1] * domg_c[1] + domg_c[2] * domg_c[2] << endl;
	//cout << endl << endl;

#endif // TEST_MODE
}


void Copter::SetStates(myTYPE v[3], myTYPE w[3], myTYPE dv[3], myTYPE dw[3]) {
	// will not update class members

	for (int i = 2; i >= 0; --i) {
		v[i] = w[i] = dv[i] = dw[i] = 0;
	}

	for (int i = 2; i >= 0; --i) {
		for (int j = 2; j >= 0; --j) {
			v[i] += refcoord.Ttransf[i][j] * vel_g[j];
			w[i] += refcoord.Etransf[i][j] * omg_g[j];
		}
	}

	//#ifdef USE_DOUBLE
	//	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, v, 1);
	//	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, w, 1);
	//#else
	//	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Ttransf, 3, vel_g, 1, 0, v, 1);
	//	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *refcoord.Etransf, 3, omg_g, 1, 0, w, 1);
	//#endif // USE_DOUBLE

#ifdef FLIGHT_TRIM
	dv[0] = dv[1] = dv[2] = 0;
	dw[0] = dw[1] = dw[2] = 0;
#else
	;
#endif // FLIGHT_TRIM


#ifdef TEST_MODE
	//printf("Copter states.\n");
	//printf("Velocity: %f, %f, %f \n", v[0], v[1], v[2]);
	//printf("Omega: %f, %f, %f \n", w[0], w[1], w[2]);
	//printf("dvelocity: %f, %f, %f \n", dv[0], dv[1], dv[2]);
	//printf("domega: %f, %f, %f \n", dw[0], dw[1], dw[2]);
	//cout << endl << endl;

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
	//printf("%f, %f, %f\n", dvel_c[0], dvel_c[1], dvel_c[2]);
	//printf("%f, %f, %f\n", domg_c[0], domg_c[1], domg_c[2]);

#endif // TEST_MODE

}


void Copter::SetInit(const int ic)
{
	myTYPE euler[3] = { 0,0,0 };
	vel_g[0] = omega*11.5 * INI_MU[ic];
	euler[1] = INI_PITCH[ic] * PI / 180;
	refcoord.SetCoordinate(euler, "euler");

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
