#include "stdafx.h"
#include "Component.h"

Fuselage::Fuselage() :Component() 
{
	;
}


Fuselage::Fuselage(const Fuselage &F) :Component(F) 
{
	dragA = F.dragA;
}


Fuselage::~Fuselage() 
{
	;
}


void Fuselage::SetAirfm(void) {
	airforce[0] = -0.5*vel[0] * vel[0] * rho * dragA;
	airforce[1] = airforce[2] = 0;
	airmoment[0] = airmoment[1] = airmoment[2] = 0;

#ifdef TEST_MODE
	//cout << endl;
	//printf("Fuselage airdynamics: \n");
	//printf("F: %f, %f, %f \n", airforce[0], airforce[1], airforce[2]);
	//printf("M: %f, %f, %f \n", airmoment[0], airmoment[1], airmoment[2]);
	cout << endl;
#endif // OUTPUT_MODE
}


void Fuselage::Fuse_UT496(void)
{
	myTYPE origin[3], euler[3];
	//print_cons_on_screen("Fuselage");
	// read config file
#ifdef UL496
	dragA = 0.0315*DISK_A;
	origin[0] = origin[1] = origin[2] = 0;
	euler[0] = euler[1] = euler[2] = 0;
	refcoord.SetCoordinate(origin, euler, refcoord.base);

#endif // UL496
	//cout << refcoord.base << endl;
}


//void Fuselage::GetModel(void)
//{
//	Fuse_UT496();
//}


