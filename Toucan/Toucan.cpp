// Toucan.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MatrixTemplate.h"
#include "Coordinate.h"
#include "Component.h"
#include "Solver.h"
#include "UnitTest.h"

#include <vector>
#include <memory>
#include <time.h>

//void foo1(Rotor &r) {
//	//按照基类类型传参时，不会进入派生类的复制构造函数
//	r.hubrtcoord.origin[0] = 9;
//}
//
//void foo2(Rotor &r) {
//	cout << r.bladecoord.base->origin[0] << endl;
//}

int main()
{
	Solver trimsolver;
	Copter helicopter;
	Component component;
	Coordinate *BASE_COORD = &helicopter.refcoord;

	
	//Fuselage fuselage;
	//Wing wing("hor"), fin1("ver", 1), fin2("ver", 2);
	//Rotor mrotor("main"), trotor("tail");
	//std::vector<Component> comps;
	std::vector<std::unique_ptr<Component>> comps;

	
	/*component.refcoord.SetBase(BASE_COORD);
	fuselage.refcoord.SetBase(BASE_COORD);
	wing.refcoord.SetBase(BASE_COORD);
	fin1.refcoord.SetBase(BASE_COORD);
	fin2.refcoord.SetBase(BASE_COORD);
	
	mrotor.refcoord.SetBase(BASE_COORD);
	mrotor.hubfxcoord.SetBase(BASE_COORD);
	mrotor.hubrtcoord.SetBase(&mrotor.hubfxcoord);
	mrotor.bladecoord.SetBase(&mrotor.hubrtcoord);
	mrotor.tppcoord.SetBase(&mrotor.hubfxcoord);

	trotor.refcoord.SetBase(BASE_COORD);
	trotor.hubfxcoord.SetBase(BASE_COORD);
	trotor.hubrtcoord.SetBase(&trotor.hubfxcoord);
	trotor.bladecoord.SetBase(&trotor.hubrtcoord);
	trotor.tppcoord.SetBase(&trotor.hubfxcoord);*/

	//myTYPE vtemp[3], wtemp[3], dvtemp[3], dwtemp[3];
	//helicopter.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//cout << "Fuselage ";
	//fuselage.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//cout << "Wing ";
	//wing.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//cout << "Fin1 ";
	//fin1.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//cout << "Fin2 ";
	//fin2.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//cout << "Main Rotor ";
	//mrotor.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//cout << "Tail Rotor ";
	//trotor.SetStates(vtemp, wtemp, dvtemp, dwtemp);
	
	//fuselage.SetAirfm();
	//wing.SetAirfm();
	//fin1.SetAirfm();
	//fin2.SetAirfm();

	/*comps.push_back(mrotor);
	comps.push_back(fuselage);
	comps.push_back(wing);
	comps.push_back(fin1);
	comps.push_back(fin2);
	comps.push_back(trotor);
	comps.push_back(component);*/

	// order fixed
	comps.emplace_back(new Rotor("main"));
	comps.emplace_back(new Fuselage);
	comps.emplace_back(new Wing("hor"));
	comps.emplace_back(new Wing("ver", 1));
	comps.emplace_back(new Wing("ver", 2));
	comps.emplace_back(new Rotor("tail"));

	component.refcoord.SetBase(BASE_COORD);
	for (auto it = comps.cbegin(); it != comps.cend(); ++it) {
		(*it)->refcoord.SetBase(BASE_COORD);
	}
	cout << comps.size() << endl;
#ifdef FLIGHT_TRIM
	comps[0]->SetCoordBase();
#ifndef TROTOR
	comps.back()->SetCoordBase();
#endif // TROTOR

#endif // FLIGHT_TRIM

	//cout << BASE_COORD << endl;
	//cout << comps.back().refcoord.base << endl;
	//cout << comps[0]->refcoord.base << endl;
	//cout << comps[3]->refcoord.base << endl;
	//cout << comps[4]->refcoord.base << endl;


	trimsolver.TrimSolver(helicopter, component, comps);
	//comps.begin()->SetAirfm();
	//comps[0].SetAirfm();
	

	//component.Assemble(comps, BASE_COORD);
	//myTYPE ftemp[3], mtemp[3];
	//component.GetAirfm_sg(ftemp, mtemp);
	//helicopter.Assemble(ftemp, mtemp);


	//helicopter.functest();
	//mrotor.functest();


	//foo2(mrotor);

	//foo1(mrotor);

	//foo2(mrotor);

	//cout << "***********************************************" << endl;
	//cout << helicopter.refcoord.base << endl;
	//cout << component.refcoord.base << endl;
	//cout << fuselage.refcoord.base << endl;
	//cout << wing.refcoord.base << endl;
	//cout << fin1.refcoord.base << endl;
	//cout << fin2.refcoord.base << endl;

	printf("\n");
	printf("Completed.\n");
	system("pause");
	return 0;
}

