#include "stdafx.h"
#include "Solver.h"



Solver::Solver()// :Copter(), Component()
{
	;
}


Solver::Solver(const Solver &S)// :Copter(S), Component(S) 
{
	;
}


Solver::~Solver() 
{
	;
}


void Solver::TrimSolver(Copter &H, Component &CP, const std::vector<std::unique_ptr<Component>> &C)
{
	myTYPE _vc[3], _wc[3], _dvc[3], _dwc[3];
	myTYPE _dv_p[3], _dw_p[3], _dv_n[3], _dw_n[3];
	myTYPE ftemp[3], mtemp[3];
	myTYPE deltt[FREEDOM] = { 0 };
#ifdef _DEBUG
	myTYPE jacob[FREEDOM][FREEDOM] = { 0 };
#endif // _DEBUG
	Matrix2<myTYPE> jacob_M(FREEDOM, FREEDOM);
	myTYPE uctrl[FREEDOM] = { 0 };
	myTYPE euler_temp[3] = { 0 };
	myTYPE uctrl_temp[FREEDOM] = { 0 };
	myTYPE epsilon = 1e-2;
	const int Nitermax = 20;
	myTYPE err_a, err_c;
	Matrix1<myTYPE> sum_a1_del(Nitermax), sum_a2_del(Nitermax), max_c_del(Nitermax);
	int iter = 0;
	const Coordinate *_coordbase = H.refcoord.base;

	err_a = 1.e-3;
	err_c = 1.e-3;
	sum_a1_del(0) = sum_a2_del(0) = max_c_del(0) = 1;

	// get initial values
	C[MROTOR]->GetCtrl(uctrl, 3);

#ifdef FLIGHT_TRIM

	C[TROTOR]->GetCtrl((uctrl + 3), 1);

#ifdef DISABLE_SLD

	uctrl[4] = H.refcoord.euler[0];
	uctrl[5] = H.refcoord.euler[1];
	euler_temp[0] = uctrl[4];
	euler_temp[1] = uctrl[5];
	euler_temp[2] = H.refcoord.euler[2];
#else

	uctrl[4] = H.refcoord.euler[2];
	uctrl[5] = H.refcoord.euler[1];
	euler_temp[0] = H.refcoord.euler[0];
	euler_temp[1] = uctrl[5];
	euler_temp[2] = uctrl[4];
#endif // DISABLE_SLD

#ifdef PARALL

	uctrl[0] = 12.1 * PI / 180;
	uctrl[1] = 1.18 * PI / 180;
	uctrl[2] = 3.12 * PI / 180;
	uctrl[3] = 10.1 * PI / 180;
#endif // PARALL

#endif // FLIGHT_TRIM

#ifdef OUTPUT_MODE
	cout << endl;
	printf("Initial Controls:\n");
	for (int iq = 0; iq < FREEDOM; ++iq) {
		cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
		cout << uctrl[iq] * 180.0 / PI << "\t";
	}
	cout << endl << endl;
#endif // OUTPUT_MODE

	H.SetStates();
	H.GetStates(_vc, _wc, _dvc, _dwc);
#ifdef FLIGHT_TRIM

#ifdef PARALL
	C[MROTOR]->SetCtrl(uctrl, 3);
	C[TROTOR]->SetCtrl((uctrl + 3), 1);
#endif //PARALL
	for (int i = 0; i <= C.size()-1; ++i) {
		C[i]->SetStates(_vc, _wc, _dvc, _dwc);
		C[i]->SetAirfm();
		
#ifdef OUTPUT_MODE
		C[i]->GetAirfm(ftemp, mtemp);
		printf("Init component aerodynamics.\n");
		printf("Component %d: \n", i);
		printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
		printf("M: %f, %f, %f \n", mtemp[0], mtemp[1], mtemp[2]);
		cout << endl << endl;
#endif // OUTPUT_MODE

	}

	CP.Assemble(C, _coordbase);
	CP.GetAirfm_sg(ftemp, mtemp);

#ifdef OUTPUT_MODE
	cout << endl;
	printf("Init guess:\n");
	printf("Total aerodynamics at CG: \n");
	printf("F: %f, %f, %f \n", ftemp[0], ftemp[1], ftemp[2]);
	printf("M: %f, %f, %f \n", mtemp[0], mtemp[1], mtemp[2]);
	cout << endl << endl;
#endif // OUTPUT_MODE

	H.Assemble(ftemp, mtemp);
	H.SetDerivs();
	H.GetStates(_vc, _wc, deltt, (deltt+3));


	for (int ifd = 0; ifd < FREEDOM; ++ifd) { uctrl_temp[ifd] = 1; }
	for (iter = 0; iter < Nitermax; ++iter) {
		// tail rotor
		uctrl[3] += epsilon;
		C[TROTOR]->SetCtrl((uctrl+3), 1);
		C[TROTOR]->SetAirfm();
		
		CP.Assemble(C, _coordbase);
		CP.GetAirfm_sg(ftemp, mtemp);
		H.SetDerivs(_dv_p, _dw_p, ftemp, mtemp);

		uctrl[3] -= 2 * epsilon;
		C[TROTOR]->SetCtrl((uctrl + 3), 1);
		C[TROTOR]->SetAirfm();

		CP.Assemble(C, _coordbase);
		CP.GetAirfm_sg(ftemp, mtemp);
		H.SetDerivs(_dv_n, _dw_n, ftemp, mtemp);

		// jacob
		for (int ijb = 0; ijb < 3; ++ijb) {
#ifdef _DEBUG
			jacob[ijb][3] = (_dv_p[ijb] - _dv_n[ijb]) / 2 / epsilon;
			jacob[ijb + 3][3] = (_dw_p[ijb] - _dw_n[ijb]) / 2 / epsilon;
#endif // _DEBUG

			jacob_M(ijb, 3) = (_dv_p[ijb] - _dv_n[ijb]) / 2 / epsilon;
			jacob_M(ijb + 3, 3) = (_dw_p[ijb] - _dw_n[ijb]) / 2 / epsilon;
		}
		uctrl[3] += epsilon;
		C[TROTOR]->SetCtrl((uctrl + 3), 1);
		C[TROTOR]->SetAirfm();

		// blade pitch variations
		for (int iq = 0; iq < 3; ++iq) {
			// main rotor
			uctrl[iq] += epsilon;
			C[MROTOR]->SetCtrl(uctrl, 3);
			C[MROTOR]->SetAirfm();
			
			CP.Assemble(C, _coordbase);
			CP.GetAirfm_sg(ftemp, mtemp);			
			H.SetDerivs(_dv_p, _dw_p, ftemp, mtemp);

			uctrl[iq] -= 2 * epsilon;
			C[MROTOR]->SetCtrl(uctrl, 3);
			C[MROTOR]->SetAirfm();

			CP.Assemble(C, _coordbase);
			CP.GetAirfm_sg(ftemp, mtemp);
			H.SetDerivs(_dv_n, _dw_n, ftemp, mtemp);

			// jacob
			for (int ijb = 0; ijb < 3; ++ijb) {
#ifdef _DEBUG
				jacob[ijb][iq] = (_dv_p[ijb] - _dv_n[ijb]) / 2 / epsilon;
				jacob[ijb + 3][iq] = (_dw_p[ijb] - _dw_n[ijb]) / 2 / epsilon;
#endif // _DEBUG

				jacob_M(ijb, iq) = (_dv_p[ijb] - _dv_n[ijb]) / 2 / epsilon;
				jacob_M(ijb + 3, iq) = (_dw_p[ijb] - _dw_n[ijb]) / 2 / epsilon;
			}
			uctrl[iq] += epsilon;
			C[MROTOR]->SetCtrl(uctrl, 3);
		}
		

		// euler
		for (int iq = 4; iq < FREEDOM; ++iq) {
			uctrl[iq] += epsilon;
#ifdef DISABLE_SLD
			euler_temp[iq-4] = uctrl[iq];
#else
			euler_temp[FREEDOM- iq] = uctrl[iq];
#endif // DISABLE_SLD
			H.refcoord.SetCoordinate(euler_temp, "euler");
			H.SetStates(_vc, _wc, _dvc, _dwc);
			for (int ic = 0; ic <= C.size() - 1; ++ic) {
				C[ic]->SetStates(_vc, _wc, _dvc, _dwc);
				C[ic]->SetAirfm();
			}

			CP.Assemble(C, _coordbase); 
			CP.GetAirfm_sg(ftemp, mtemp);
			H.SetDerivs(_dv_p, _dw_p, ftemp, mtemp);

			uctrl[iq] -= 2 * epsilon;
#ifdef DISABLE_SLD
			euler_temp[iq - 4] = uctrl[iq];
#else
			euler_temp[FREEDOM - iq] = uctrl[iq];
#endif // DISABLE_SLD
			H.refcoord.SetCoordinate(euler_temp, "euler");
			H.SetStates(_vc, _wc, _dvc, _dwc);
			for (int ic = 0; ic <= C.size() - 1; ++ic) {
				C[ic]->SetStates(_vc, _wc, _dvc, _dwc);
				C[ic]->SetAirfm();
			}
			CP.Assemble(C, _coordbase); 
			CP.GetAirfm_sg(ftemp, mtemp);
			H.SetDerivs(_dv_n, _dw_n, ftemp, mtemp);
			// jacob
			for (int ijb = 0; ijb < 3; ++ijb) {
#ifdef _DEBUG
				jacob[ijb][iq] = (_dv_p[ijb] - _dv_n[ijb]) / 2 / epsilon;
				jacob[ijb + 3][iq] = (_dw_p[ijb] - _dw_n[ijb]) / 2 / epsilon;
#endif // _DEBUG

				jacob_M(ijb, iq) = (_dv_p[ijb] - _dv_n[ijb]) / 2 / epsilon;
				jacob_M(ijb + 3, iq) = (_dw_p[ijb] - _dw_n[ijb]) / 2 / epsilon;
			}

			uctrl[iq] += epsilon;
#ifdef DISABLE_SLD
			euler_temp[iq - 4] = uctrl[iq];
#else
			euler_temp[FREEDOM - iq] = uctrl[iq];
#endif // DISABLE_SLD
			H.refcoord.SetCoordinate(euler_temp, "euler"); 
		}

#ifdef OUTPUT_MODE_1
		cout << endl;
		printf("Accelaration:\n");
		for (int iq = 0; iq < FREEDOM; ++iq) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << deltt[iq] << "\t";
		}
		cout << endl << endl;

		printf("Jacob:\n");
		for (int i = 0; i < FREEDOM; ++i) {
			for (int j = 0; j < FREEDOM; ++j) {
				cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
				cout << jacob_M(i, j) << "\t";
			}
			cout << endl;
		}
#endif // !OUTPUT_MODE_1

		// update deltt, flight states, aerodynamics etc. 
		Msolver(jacob_M.v_p, deltt, 6);
		for (int iq = 0; iq < FREEDOM; ++iq) {
			uctrl[iq] -= deltt[iq];
		}
		// limit
		int ic_flg = 0;
		if (Abs(uctrl[0]) > SITA_COLL_MAX) { 
			uctrl[0] = SITA_COLL_MAX*Sign(uctrl[0]); 
			//uctrl[0] = 0.5*SITA_COLL_MAX*Sign(uctrl[0]) + 0.5*uctrl[0];
			++ic_flg;
		}
		if (Abs(uctrl[3]) > SITA_COLL_MAX) { 
			uctrl[3] = SITA_COLL_MAX*Sign(uctrl[3]); 
			//uctrl[3] = 0.5*SITA_COLL_MAX*Sign(uctrl[3]) + 0.5*uctrl[3];
			++ic_flg;
		}
		if (Abs(uctrl[1]) > SITA_CYCL_MAX) { 
			uctrl[1] = SITA_CYCL_MAX*Sign(uctrl[1]); 
			//uctrl[1] = 0.5*SITA_CYCL_MAX*Sign(uctrl[1]) + 0.5*uctrl[1];
			++ic_flg;
		}
		if (Abs(uctrl[2]) > SITA_CYCL_MAX) { 
			uctrl[2] = SITA_CYCL_MAX*Sign(uctrl[2]); 
			//uctrl[2] = 0.5*SITA_CYCL_MAX*Sign(uctrl[2]) + 0.5*uctrl[2];
			++ic_flg;
		}
		if (Abs(uctrl[4]) > EULER_MAX) { 
			uctrl[4] = EULER_MAX*Sign(uctrl[4]); 
			//uctrl[4] = 0.5*EULER_MAX*Sign(uctrl[4]) + 0.5*uctrl[4];
			++ic_flg;
		}
		if (Abs(uctrl[5]) > EULER_MAX) { 
			uctrl[5] = EULER_MAX*Sign(uctrl[5]); 
			//uctrl[5] = 0.5*EULER_MAX*Sign(uctrl[5]) + 0.5*uctrl[5];
			++ic_flg;
		}

#ifdef OUTPUT_MODE_1
		cout << endl;
		printf("Delta controls (degs).\n");
		for (int iq = 0; iq < FREEDOM; ++iq) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << deltt[iq] * 180.0 / PI << "\t";
		}
		cout << endl << endl;

		printf("Delta controls (rads).\n");
		for (int iq = 0; iq < FREEDOM; ++iq) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << deltt[iq] << "\t";
		}
		cout << endl << endl;

		printf("New controls (degs).\n");
		for (int iq = 0; iq < FREEDOM; ++iq) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << uctrl[iq]*180.0/PI << "\t";
		}
		cout << endl << endl;

		printf("New controls (rads).\n");
		for (int iq = 0; iq < FREEDOM; ++iq) {
			cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
			cout << uctrl[iq] << "\t";
		}
		cout << endl << endl;
#endif // !OUTPUT_MODE_1

		if (ic_flg == FREEDOM) {
			printf("Controls reach boundary.");
			break;
		}
		
#ifdef DISABLE_SLD
		euler_temp[0] = uctrl[4];
#else
		euler_temp[2] = uctrl[4];
#endif // DISABLE_SLD
		euler_temp[1] = uctrl[5];

		H.refcoord.SetCoordinate(euler_temp, "euler");
		H.SetStates();
		H.GetStates(_vc, _wc, _dvc, _dwc);
		C[MROTOR]->SetCtrl(uctrl, 3);
		C[TROTOR]->SetCtrl((uctrl + 3), 1);
		for (int ic = 0; ic <= C.size() - 1; ++ic) {
			C[ic]->SetStates(_vc, _wc, _dvc, _dwc);
			C[ic]->SetAirfm();
		}
		CP.Assemble(C, _coordbase);
		CP.GetAirfm_sg(ftemp, mtemp);
		H.SetDerivs(ftemp, mtemp);
		H.GetStates(_vc, _wc, deltt, (deltt+3));

		max_c_del.v_p[iter] = Abs(uctrl[0] / uctrl_temp[0] - 1);
		uctrl_temp[0] = uctrl[0];
		for (int ifd = 1; ifd < FREEDOM; ++ifd) {
			max_c_del.v_p[iter] = Max(max_c_del.v_p[iter], Abs(uctrl[ifd] / uctrl_temp[ifd] - 1));
			uctrl_temp[ifd] = uctrl[ifd];
		}
		sum_a1_del.v_p[iter] = deltt[0] * deltt[0] + deltt[1] * deltt[1] + deltt[2] * deltt[2];
		sum_a2_del.v_p[iter] = deltt[3] * deltt[3] + deltt[4] * deltt[4] + deltt[5] * deltt[5];

		if (Max(sum_a1_del.v_p[iter],sum_a2_del.v_p[iter]) < (err_a*err_a) && max_c_del.v_p[iter] < err_c) {
			break; 
		}

	}

#ifdef OUTPUT_MODE
	cout << endl;
	printf("Iter: %d, Total aerodynamics at CG: \n", iter);
#ifdef OUTPUT_MODE_1
	printf("Accelaration error: %f, %f\n", sum_a1_del.v_p[iter], sum_a2_del.v_p[iter]);
	printf("Control converge: %f \n", max_c_del.v_p[iter]);
#endif // OUTPUT_MODE_1
	printf("Finally controls: \n");
	for (int iq = 0; iq < FREEDOM; ++iq) {
		cout << std::right << std::setw(8) << std::setprecision(4) << std::fixed << std::showpoint;
		cout << uctrl[iq] * 180.0 / PI << "\t";
	}
#endif // OUTPUT_MODE

	sum_a1_del.outputs("sum_a1_del.output", 8);
	sum_a2_del.outputs("sum_a2_del.output", 8);
	max_c_del.outputs("max_c_del.output", 8);
	ofstream OutFile("uctrl.output", std::ios::app);
	//ofstream OutFile("uctrl.output");
	if (OutFile) {
		OutFile.precision(6);
		for (int i = 0; i < FREEDOM; ++i) {
			OutFile << std::right << std::setw(6 + 2) << std::setprecision(6) << std::fixed << std::showpoint;
			OutFile << uctrl[i]*180.0/PI << '\t';
		}
		OutFile << endl;
		OutFile.close();
	}
	else {
		cout << "uctrl.output" << " open failed." << endl;
		system("pause");
	}

#endif // FLIGHT_TRIM



}
