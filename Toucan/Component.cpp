#include "stdafx.h"
#include "Component.h"


Component::Component() :Copter() {
	//print_cons_on_screen("Component");
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
	//print_cpcons_on_screen("Component");
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
	//print_dscons_on_screen("Component");
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
#ifdef TEST_MODE
	//printf("Component states.\n");
	//printf("Velocity: %f, %f, %f \n", vel[0], vel[1], vel[2]);
	//printf("Omega: %f, %f, %f \n", omg[0], omg[1], omg[2]);
	//printf("dvelocity: %f, %f, %f \n", dvel[0], dvel[1], dvel[2]);
	//printf("domega: %f, %f, %f \n", domg[0], domg[1], domg[2]);
	//cout << vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2] << endl;
	//cout << omg[0] * omg[0] + omg[1] * omg[1] + omg[2] * omg[2] << endl;
	//cout << dvel[0] * dvel[0] + dvel[1] * dvel[1] + dvel[2] * dvel[2] << endl;
	//cout << domg[0] * domg[0] + domg[1] * domg[1] + domg[2] * domg[2] << endl;

	cout << endl << endl;

#endif // TEST_MODE

}


void Component::SetAirfm_cg(const Coordinate *base) {
	myTYPE rcb[3];
	const myTYPE *origin1_ptr = refcoord.origin;
	const myTYPE *origin2_ptr = base->origin;

	for (int j = 2; j >= 0; --j) { airforce_cg[j] = airmoment_cg[j] = 0; }

	//refcoord.Transfer(tcb, rcb, refcoord, *base);
	if (refcoord.base == base) {
		for (int i = 2; i >= 0; --i) {
			rcb[i] = *(origin2_ptr + i) - *(origin1_ptr + i);			
		}
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				// give a transpose by j i swap of Ttransf[][]
				//airforce_cg[i] += refcoord.Ttransf[j][i] * airforce[j];
				// index friendly
				airforce_cg[j] += refcoord.Ttransf[i][j] * airforce[i];
			}
		}
		Cross(airmoment_cg, airforce_cg, rcb);
		for (int i = 2; i >= 0; --i) {
			for (int j = 2; j >= 0; --j) {
				// give a transpose by j i swap of Ttransf[][]
				//airmoment_cg[i] += refcoord.Ttransf[j][i] * airmoment[j];
				// index friendly
				airmoment_cg[j] += refcoord.Ttransf[i][j] * airmoment[i];
			}
		}
	}
	else { wrong_comp_coordinate_diff_base(); }
	

//#ifdef USE_DOUBLE
//	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airforce, 1, 0, airforce_cg, 1);
//	
//	Cross(airmoment_cg, airforce_cg, rcb);
//	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airmoment, 1, 1, airmoment_cg, 1);
//	//vdAdd(3, temp_cross, airmoment_cg, airmoment_cg);
//
//#else
//	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airforce, 1, 0, airforce_cg, 1);
//
//	Cross(airmoment_cg, airforce_cg, rcb);
//	cblas_sgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1, *tcb, 3, airmoment, 1, 1, airmoment_cg, 1);
//
//#endif // USE_DOUBLE
}


void Component::Assemble(const std::vector<std::unique_ptr<Component>> &C, const Coordinate *base) {
	_assemble(airforce_sigma, airmoment_sigma, C, base);
}














