#include "stdafx.h"
#include "Component.h"

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
		euler[1] = 4 * PI / 180;
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
	myTYPE f[3] = { 0 };
	myTYPE m[3] = { 0 };
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
	// freestream direction forces and moments, not refcoord
	f[0] = -0.5 * rho * stot * vel2 * cd;
	f[1] = 0;
	f[2] = 0.5 * rho * stot * vel2 * cl;
	m[0] = 0;
	m[1] = 0;
	m[2] = 0;
	// forces and moments at refcoord
	airforce[0] = f[0] * cos(aoa) + f[2] * sin(aoa);
	airforce[2] = -f[0] * sin(aoa) + f[2] * cos(aoa);
	airforce[1] = 0;
	airmoment[0] = m[0] * cos(aoa) + m[2] * sin(aoa);
	airmoment[1] = -m[0] * sin(aoa) + m[2] * cos(aoa);
#endif
#ifdef TEST_MODE
	//cout << endl;
	//printf("%s wing airdynamics: \n", type);
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // TEST_MODE


}