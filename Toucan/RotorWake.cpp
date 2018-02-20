#include "stdafx.h"
#include "Components.h"


void Rotor::_wakeInducedVel(void)
{
	switch (adyna)
	{
	case PWake:
		double veltpp[3] = { 0.0,0.0,0.0 };
		double dveltpp[3] = { 0.0,0.0,0.0 };
		_velTransform(veltpp, dveltpp, tppcoord);
		for (int i = 0; i < 3; i++)
			veltpp[i] /= vtipa;

		for (int iw = 0; iw < wakev.size(); iw++)
			wakev[iw]._computeGeometry(adyna, veltpp, lambdi_ag, beta[0]);
		
		haveStr = true;
		haveGeo = true;// _wakeGeoBd();

		break;
	}
}

bool Rotor::_tipVortexStr(void)
{
	int igen = 0, istr = 0;
	Matrix1<int> idns;

	// allocate memory
	istr = outboard * ns;
	//tipstr.allocate(nk, nf);
	idns.allocate(ns - istr);

	// obtian tip vortices strengths	
	idns = step(istr, ns - 1);
	for (int j = nf - 1; j >= 0; --j) {
		for (int i = nk - 1; i >= 0; --i) {
			igen = j - i;
			while (igen < 0) igen += nf;
			tipstr(i, j) = cirlb(igen, idns).findmax();
		}
	}
	//cirlb.output("DEBUG_cirlb.output", 6);
	return true;
}

void Rotor::_bladePosition()
{
	//尾迹法速度计算点，目前取四分之三弦长位置
	//有后掠角时还要调整
	myTYPE azm, bfp, rst;
	myTYPE ch;
	myTYPE twistt;
	for (int j = ns - 1; j >= 0; --j)
	{
		twistt = twist(j) - twist(0) + pitchroot;
		ch = chord(j) / radius;
		for (int i = nf - 1; i >= 0; --i)
		{
			azm = azstation(i, j);
			bfp = beta[0] + beta[1] * cos(azm) + beta[2] * sin(azm);
			rst = rastation(i, j);
			bladedeform(i, j, 0) = ((rst - eflap)*cos(bfp) + eflap) * cos(azm) + ch*0.5*sin(azm);
			bladedeform(i, j, 1) = ((rst - eflap)*cos(bfp) + eflap) * sin(azm) - ch*0.5*cos(azm);
			bladedeform(i, j, 2) = bfp*(rst - eflap) -0.5*ch*(sita[0] + sita[1] * cos(azm) + sita[2] * sin(azm) + twistt);
		}
	}
	//if (bladedeform(1, ns - 1, 1) < 0)
	//bladedeform.output("DEBUG_bladedeform.output", 6);
}

bool Rotor::_wakeGeoBd()
{
	//有后掠角时还要调整
	myTYPE r0, z0, x0, xv, yv, zv, xe, df, ch;
	myTYPE ka, E, kx, ky, ky3, k0;
	myTYPE veltpp[3] = { 0.0,0.0,0.0 };
	myTYPE dveltpp[3] = { 0.0,0.0,0.0 };
	myTYPE temp = 0;
	myTYPE twistt = twist(ns - 1) - twist(0) + pitchroot;

	// shed point at tpp
	r0 = rtip * cos(beta[0]);
	//z0 = rtip * beta[0];
	ch = chord(0) / radius;

	// Beddoes prescribed wake parameters obtained
	veltpp[0] = veltpp[1] = veltpp[2] = 0;
	_velTransform(veltpp, dveltpp, tppcoord);

	ka = atan2(-veltpp[0] / vtipa, -(-veltpp[2] / vtipa + lambdi_ag));
	E = ka / 2.0;
	kx = E; ky = 2 * veltpp[0] / vtipa; ky3 = -kx; k0 = 1.0 - 8.0 * ky3 / 15.0 / PI;

	// compute tip vortices geometry in TPP coord
	df = 2.0 * PI / nf; xe = 0;
	for (int j = nf - 1; j >= 0; --j) {
		for (int i = nk - 1; i >= 0; --i) {
			x0 = r0 * cos(df * (j - i)) + ch*0.75*sin(df*(j - i));
			xv = x0 - veltpp[0] / vtipa * df * i;
			yv = r0 * sin(df * (j - i)) - ch*0.75*cos(df*(j - i));
			if (xv * xv + yv * yv < 1.01) {
				zv = -lambdi_ag * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * i*df - lambdi_ag / 2.0 * kx * (x0 + xv) * i*df;
			}
			else {			
				if (yv > 1)
					xe = 0;
				else
					xe = sqrt(1.0 - yv*yv);
				//zv = -lambdi_ag / veltpp[0] * vtipa * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * (xe - x0) - lambdi_ag / 2.0 * kx * (x0 + xe) * i*df;
				zv = -lambdi_ag / (-veltpp[0]) * vtipa * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * (xe - x0) - lambdi_ag / 2.0 * kx / (-veltpp[0]) * vtipa * (xe*xe - x0*x0);
				zv -= 2.0 * lambdi_ag * vtipa / (-veltpp[0]) * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * (xv - xe);
			}
			tipgeometry(i, j, 0) = xv;
			tipgeometry(i, j, 1) = yv;
			temp = -zv;// +(rtip - eflap)*(beta[0] + beta[1] * cos(df*(j - i)) + beta[2] * sin(df*(j - i)));
			//temp -= 0.75*ch*(sita[0] + sita[1] * cos(df*(j - i)) + sita[2] * sin(df*(j - i)) + twistt);
			tipgeometry(i, j, 2) = temp -veltpp[2] / vtipa*df*i;// temp;
		}
	}

	//if (outputWake)
	//tipgeometry.output("DEBUG_tipgeo.output", 6);
	return true;
}

void Rotor::_wakeIndVelCalc()
{
	// 适用于nb = 2
	myTYPE temp_lambdi, temp_lambdx, temp_lambdy;
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	Matrix2<myTYPE> tipgeoexpand, tempM(3, 3);
	Matrix3<myTYPE> tipgeoathub(nk, nf, 3);
	myTYPE *str0_b0_ptr, *str0_b1_ptr, *str1_b0_ptr, *str1_b1_ptr;
	myTYPE *temp_str0_b0_ptr, *temp_str0_b1_ptr;
	myTYPE *tipstr_vpp = tipstr.v_p;
	myTYPE *tipgeo_vpp = tipgeoathub.v_p;
	int tipstr_NI = tipstr.NI;
	int tipgeo_NI = nk;
	int tipgeo_NJ = nf;
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp;
	//myTYPE **tipgeoathub_x, **tipgeoathub_y, **tipgeoathub_z;
	//myTYPE **temp_tipgeoathub_x, **temp_tipgeoathub_y, **temp_tipgeoathub_z;
	std::vector<myTYPE*> tipgeoathub_x(nb), tipgeoathub_y(nb), tipgeoathub_z(nb);
	std::vector<myTYPE*> temp_tipgeoathub_x(nb), temp_tipgeoathub_y(nb), temp_tipgeoathub_z(nb);
	//tipgeoathub_x = new myTYPE*[nb], tipgeoathub_y = new myTYPE*[nb], tipgeoathub_z = new myTYPE*[nb];
	//temp_tipgeoathub_x = new myTYPE*[nb], temp_tipgeoathub_y = new myTYPE*[nb], temp_tipgeoathub_z = new myTYPE*[nb];

	//myTYPE rcb[3];

	rc04 = rc0 * rc0 * rc0 * rc0;
	tipgeoexpand = tipgeometry.reshape(3, 3, nk*nf); // 应该要减去在TPP坐标系下TPP原点指向HUB原点的向量
	//for (int i = 0; i < 3; ++i)
	//	rcb[i] = hubfxcoord.origin[i] - tppcoord.origin[i];
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			tempM(i, j) = tppcoord.Ttransf[i][j];
	tipgeoexpand = tempM.transpose().matrixmultiplyTP(tipgeoexpand);
	for (int k = 0; k < 3; ++k)
		for (int j = 0; j < nk*nf; ++j)
			tipgeoathub(j%nk, j / nk, k) = tipgeoexpand(k, j);

	for (int iz = nf - 1, iz2 = iz + nf / nb; iz >= 0; --iz, iz2 = iz + nf / nb)
	{
		/*if (nb != 2)
		{
			printf("Undefined nb != 2 in _wakeIndVelCalc() Func. \n");
			system("pause");
			return;
		}*/
		if (iz2 >= nf) { iz2 -= nf; }

		str0_b0_ptr = tipstr_vpp + (nk - 1 + iz*tipstr_NI); //tipstr(nk - 1, iz);
		str0_b1_ptr = tipstr_vpp + (nk - 1 + iz2*tipstr_NI);//tipstr(nk - 1, iz2);

		tipgeoathub_x[0] = tipgeo_vpp + ((nk - 1) + 0 + iz*tipgeo_NI); //tipgeoexpand_athub(nk - 1, iz, 0);
		tipgeoathub_y[0] = tipgeo_vpp + ((nk - 1) + tipgeo_NI*tipgeo_NJ + iz*tipgeo_NI); // tipgeoexpand_athub(nk - 1, iz, 1);
		tipgeoathub_z[0] = tipgeo_vpp + ((nk - 1) + 2 * tipgeo_NI*tipgeo_NJ + iz*tipgeo_NI); // tipgeoexpand_athub(nk - 1, iz, 2);
		tipgeoathub_x[1] = tipgeo_vpp + ((nk - 1) + 0 + iz2*tipgeo_NI); //tipgeoexpand_athub(nk - 1, iz2, 0);
		tipgeoathub_y[1] = tipgeo_vpp + ((nk - 1) + tipgeo_NI*tipgeo_NJ + iz2*tipgeo_NI); // tipgeoexpand_athub(nk - 1, iz2, 1);
		tipgeoathub_z[1] = tipgeo_vpp + ((nk - 1) + 2 * tipgeo_NI*tipgeo_NJ + iz2*tipgeo_NI); // tipgeoexpand_athub(nk - 1, iz2, 2);

		for (int ir = ns - 1; ir >= 0; --ir)
		{
			temp_lambdi = 0;
			temp_lambdx = 0;
			temp_lambdy = 0;
			xblade = bladedeform(iz, ir, 0);
			yblade = bladedeform(iz, ir, 1);
			zblade = bladedeform(iz, ir, 2);

			temp_str0_b0_ptr = str0_b0_ptr;
			temp_str0_b1_ptr = str0_b1_ptr;
			// balde 1 $ik == nk - 1$ element points to blade 1 
			temp_tipgeoathub_x[0] = tipgeoathub_x[0];
			temp_tipgeoathub_y[0] = tipgeoathub_y[0];
			temp_tipgeoathub_z[0] = tipgeoathub_z[0];

			ri0 = xblade - (*(temp_tipgeoathub_x[0]));
			rj0 = yblade - (*(temp_tipgeoathub_y[0]));
			rk0 = zblade - (*(temp_tipgeoathub_z[0]));
			r0len = norm(ri0, rj0, rk0);
			for (int ik = nk - 2; ik >= 0; --ik)
			{
				// strength				
				// temp_str1_b0_ptr = temp_str0_b0_ptr - 1; // !previous! element
				// geometry
				ri1 = xblade - (*(--temp_tipgeoathub_x[0]));
				rj1 = yblade - (*(--temp_tipgeoathub_y[0]));
				rk1 = zblade - (*(--temp_tipgeoathub_z[0]));
				cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
				//height = norm(rcros[0], rcros[1], rcros[2]) / norm(ri0 - ri1, rj0 - rj1, rk0 - rk1);
				height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));
				//height2 = height * height;

				r1len = norm(ri1, rj1, rk1);
				rdot = dot(ri0, rj0, rk0, ri1, rj1, rk1);
				geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
				//geofunc *= height2 / sqrt(height2 * height2 + rc04);
				geofunc *= (1 - 0.5 * rc04 / height2 / height2);
				geofunc *= 0.5*(*(--temp_str0_b0_ptr) + (*temp_str0_b0_ptr));
				_temp = 0.25 / PI * geofunc;

				temp_lambdx -= rcros[0] * _temp;
				temp_lambdy -= rcros[1] * _temp;
				temp_lambdi -= rcros[2] * _temp;

				// save current element related variables
				// str0_b0_ptr = str1_b0_ptr;
				//--temp_str0_b0_ptr;
				ri0 = ri1;
				rj0 = rj1;
				rk0 = rk1;
				r0len = r1len;
			}
			// balde 2 $ik == nk - 1$ element points to blade 1 
			temp_tipgeoathub_x[1] = tipgeoathub_x[1];
			temp_tipgeoathub_y[1] = tipgeoathub_y[1];
			temp_tipgeoathub_z[1] = tipgeoathub_z[1];

			ri0 = xblade - (*(temp_tipgeoathub_x[1]));
			rj0 = yblade - (*(temp_tipgeoathub_y[1]));
			rk0 = zblade - (*(temp_tipgeoathub_z[1]));
			r0len = norm(ri0, rj0, rk0);
			for (int ik = nk - 2; ik >= 0; --ik) {
				// strength
				// str1_b1_ptr = str0_b1_ptr - 1; // !previous!element
				// geometry
				ri1 = xblade - (*(--temp_tipgeoathub_x[1]));
				rj1 = yblade - (*(--temp_tipgeoathub_y[1]));
				rk1 = zblade - (*(--temp_tipgeoathub_z[1]));
				cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
				//height = norm(rcros[0], rcros[1], rcros[2]) / norm(ri0 - ri1, rj0 - rj1, rk0 - rk1);
				height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));
				//height2 = height * height;

				r1len = norm(ri1, rj1, rk1);
				rdot = dot(ri0, rj0, rk0, ri1, rj1, rk1);
				geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
				//geofunc *= height2 / sqrt(height2 * height2 + rc04);
				geofunc *= (1 - 0.5 * rc04 / height2 / height2);
				geofunc *= 0.5*(*(--temp_str0_b1_ptr) + (*temp_str0_b1_ptr)); // get value first,  left shift and get value
				_temp = 0.25 / PI * geofunc;

				temp_lambdx -= rcros[0] * _temp;
				temp_lambdy -= rcros[1] * _temp;
				temp_lambdi -= rcros[2] * _temp;

				// save current element related variables
				//str0_b0_ptr = str1_b0_ptr;
				//--temp_str0_b1_ptr;
				ri0 = ri1;
				rj0 = rj1;
				rk0 = rk1;
				r0len = r1len;
			}

			lambdi(iz, ir) = temp_lambdi / radius / vtipa;
			lambdx(iz, ir) = temp_lambdx / radius / vtipa;
			lambdy(iz, ir) = temp_lambdy / radius / vtipa;
		}	
	}
	lambdh = lambdi - vel[2] / vtipa;
	//lambdi.output("DEBUG_lambdi.output", 10);
}


void Rotor::_wakeInducedVelMP(int nb)
{
	myTYPE temp_lambdi, temp_lambdx, temp_lambdy;
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	Matrix2<myTYPE> tipgeoexpand, tempM(3, 3);
	Matrix3<myTYPE> tipgeoathub(nk, nf, 3);
	
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp, az, _c;

	int iz, ik;
	double _dfi, df;

	rc04 = rc0 * rc0 * rc0 * rc0;
	tipgeoexpand = tipgeometry.reshape(3, 3, nk*nf); 
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			tempM(i, j) = tppcoord.Ttransf[i][j];
	tipgeoexpand = tempM.transpose().matrixmultiplyTP(tipgeoexpand);

	df = 2 * PI / nf;

	for (int k = 0; k < 3; ++k)
	{
		for (int j = 0; j < nk*nf; ++j)
		{
			iz = j%nk;
			ik = j / nk;
			_dfi = df*(iz - ik);
			tipgeoathub(iz, ik, k) = tipgeoexpand(k, j);
			if (k == 2)
			{
				tipgeoathub(iz, ik, 2) += (rtip - eflap)*(beta[0] + beta[1] * cos(_dfi) + beta[2] * sin(_dfi));// -chord(0) / radius;
				tipgeoathub(iz, ik, 2) -= 0.75*chv*(sita[0] + sita[1] * cos(_dfi) + sita[2] * sin(_dfi) + twistv - twist(0) + pitchroot);
			}
		}
	}

	for (int iz = nf - 1; iz >= 0; iz--)
	{	
		az = _limitaz(iz*df);
		for (int ir = ns - 1; ir >= 0; --ir)
		{
			temp_lambdi = 0;
			temp_lambdx = 0;
			temp_lambdy = 0;

			xblade = bladedeform(iz, ir, 0);
			yblade = bladedeform(iz, ir, 1);
			zblade = bladedeform(iz, ir, 2);

			_c = chord(ir) / radius;

			for (int ib = 0; ib < nb; ib++)
			{
				int iz2 = iz + nf / nb*ib;
				if (iz2 >= nf)
					iz2 -= nf;

				// balde 1 $ik == nk - 1$ element points to blade 1 
				ri0 = xblade - tipgeoathub(nk - 1, iz2, 0);
				rj0 = yblade - tipgeoathub(nk - 1, iz2, 1);
				rk0 = zblade - tipgeoathub(nk - 1, iz2, 2);
				r0len = Norm(ri0, rj0, rk0);
				for (int ik = nk - 2; ik >= 0; --ik)
				{
					ri1 = xblade - tipgeoathub(ik, iz2, 0);
					rj1 = yblade - tipgeoathub(ik, iz2, 1);
					rk1 = zblade - tipgeoathub(ik, iz2, 2);
					r1len = Norm(ri1, rj1, rk1);

					cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
					height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));
					
					rdot = Dot(ri0, rj0, rk0, ri1, rj1, rk1);
					geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
					geofunc *= height2 / rootNewton(rc04 + height2*height2, height2, height2*1e-2);
					//geofunc *= height2 / sqrt(rc04 + height2*height2);
					geofunc *= 0.5*(tipstr(ik, iz2) + tipstr(ik + 1, iz2));
					
					if (r1len*r1len + r0len*r0len + 2 * rdot < 400 * _c*_c)
						geofunc *= lscorr.LSCorrection(az, _c, rc0, ri0, rj0, rk0, ri1, rj1, rk1);
					
					_temp = 0.25 / PI * geofunc;

					temp_lambdx -= rcros[0] * _temp;
					temp_lambdy -= rcros[1] * _temp;					
					temp_lambdi -= rcros[2] * _temp;

					ri0 = ri1;
					rj0 = rj1;
					rk0 = rk1;
					r0len = r1len;
					//
					//if (ib==0 && iz==19)
					//	DEBUG_height2(ik, ir) = height2;
				}
			}
			lambdi(iz, ir) = -temp_lambdi / radius / vtipa;
			lambdx(iz, ir) = -temp_lambdx / radius / vtipa;
			lambdy(iz, ir) = -temp_lambdy / radius / vtipa;
		}
	}
	lambdh = lambdi - vel[2] / vtipa;
	//lambdi.output("DEBUG_lambdi.output", 10);
	//tipgeoathub.output("DEBUG_tipgeo_athub.output", 10);
	//DEBUG_height2.output("DEBUG_height2.output", 10);
}

void Rotor::_wakeInducedVelMP(void)
{
	Matrix2<myTYPE> tempM(3, 3);
	Matrix2<myTYPE> geoexp;
	double temp_lambdx, temp_lambdy, temp_lambdi;
	double _lambdx, _lambdy, _lambdi;

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			tempM(i, j) = tppcoord.Ttransf[i][j];


	for (int i = 0; i<wakev.size();i++)
	{
		int iz, ik, nk, nf, ns;
		double _dfi, bfp, sth;

		nk = wakev[i].nk;
		nf = wakev[i].nf;
		ns = wakev[i].ns;

		if (wakev[i].istip || wakev[i].isrot)
		{
			geoexp = tempM.transpose().matrixmultiplyTP(wakev[i].geoexp);
					
			for (int k = 0; k < 3; ++k)
			{
				for (int j = 0; j < nk*nf; ++j)
				{
					iz = j%nk;
					ik = j / nk;
					_dfi = 2*PI/nf*(iz - ik);
					bfp = beta[0] + beta[1] * cos(_dfi) + beta[2] * sin(_dfi);
					sth = sita[0] + sita[1] * cos(_dfi) + sita[2] * sin(_dfi) + wakev[i].twistv - twist(0) + pitchroot;
					wakev[i].geometryathub(iz, ik, k) = geoexp(k, j);

					if (k == 2)
						wakev[i].geometryathub(iz, ik, 2) += (1 - eflap)*bfp -0.75*wakev[i].chv*sth;
				}
			}
		}
		else if (wakev[i].isbond)
		{
			for (int iz = 0; iz < nf; iz++)
			{  
				_dfi = azstation(iz, 0);
				bfp = (beta[0] + beta[1] * cos(_dfi) + beta[2] * sin(_dfi));
				for (int ik = 0; ik < nk; ik++)
				{
					wakev[i].geometryathub(iz, ik, 0) = wakev[i].geometry(iz, ik, 0) * cos(_dfi);
					wakev[i].geometryathub(iz, ik, 1) = wakev[i].geometry(iz, ik, 0) * sin(_dfi);
					wakev[i].geometryathub(iz, ik, 2) = bfp*wakev[i].geometry(iz, ik, 0);
				}
			}
		}
		//wakev[i].geometryathub.output("DEBUG_tipgeo_athub.output", 10);
		//wakev[i].geometry.output("DEBUG_tipgeo.output", 6);
	}

	double cp[3] = { 0,0,0 };
	for (int iz = nf - 1; iz >= 0; iz--)
	{
		for (int ir = ns - 1; ir >= 0; --ir)
		{
			cp[0] = bladedeform(iz, ir, 0); 
			cp[1] = bladedeform(iz, ir, 1);
			cp[2] = bladedeform(iz, ir, 2);

			temp_lambdi = temp_lambdx = temp_lambdy = 0;

			for (int ib = 0; ib < nb; ib++)
			{
				int iz2 = iz + nf / nb*ib;
				while (iz2 >= nf)
					iz2 -= nf;
				for (int iw = wakev.size() - 1; iw >= 0; --iw)
				{					
					wakev[iw].ComputeIndVel(_lambdx, _lambdy, _lambdi, iz, ir, iz2, ib, chord(ir) / radius, cp);
					temp_lambdi += _lambdi;
					temp_lambdx += _lambdx;
					temp_lambdy += _lambdy;
				}
			}
			lambdi(iz, ir) = -temp_lambdi / radius / vtipa;
			lambdx(iz, ir) = -temp_lambdx / radius / vtipa;
			lambdy(iz, ir) = -temp_lambdy / radius / vtipa;
		}
	}
	lambdh = lambdi - vel[2] / vtipa;
	//lambdi.output("DEBUG_lambdi.output", 10);
}


double LSCorr::LSCorrection(double az, double c, double rc, double ri0, double rj0, double rk0, double ri1, double rj1, double rk1)
{
	//全采用以R为基准的无量纲长度
	double r0L2, r1L2, rdot, ssq, s, s1, s2;
	double rmx, rmy, rmz;
	double cs, sn;
	double jsis, cosL, sinL, angLS;
	double binv, hLS, sgn, rsinL;
	double b0, b1, b2, aop, a1, a2, cop, c0, c1, c2;
	double rs2, h0, h1, h2, x0, x1;
	double xx, hh, xh;
	double LLS, LLL;
	
	
	r0L2 = ri0*ri0 + rj0*rj0 + rk0*rk0;
	r1L2 = ri1*ri1 + rj1*rj1 + rk1*rk1;
	rdot = Dot(ri0, rj0, rk0, ri1, rj1, rk1);
	ssq = r0L2 + r1L2 - 2 * rdot;
	s = sqrt(ssq);
	s1 = (rdot - r0L2) / (s + 1e-6);
	s2 = s1 + s;
	rmx = ri0*s2 - ri1*s1;
	rmy = rj0*s2 - rj1*s1;
	rmz = rk0*s2 - rk1*s1;

	cs = cos(az);
	sn = sin(az);

	jsis = (cs*(ri0 - ri1) + sn*(rj0 - rj1)) / (s + 1e-6);
	cosL = Max(-1, Min(0, -Abs(jsis)));
	sinL = Max(0, Min(1, 1 - jsis*jsis));
	sinL = sqrt(sinL);//, 0.5, 1e-3);
	angLS = 180 - 57.29578*asin(sinL); //角度
	angLS = Max(90.001, Min(angLS, 179.99));

	if (c > 0)
		binv = 2 / c; //注意量纲
	else
		binv = 0;
	hLS = rmz*rmz;
	if (rc > 0)
		hLS += rc*rc;
	//printf("Rmz = %f, rc = %f, hLS = %f, root = %f\n", rmz, rc, hLS, rootNewton(hLS, rmz, 1e-6*rmz));
	hLS = binv * rootNewton(hLS, rmz, 1e-6*rmz);
	
	sgn = Sign((cs*rmx + sn*rmy)*jsis);
	rsinL = sqrt(rmx*rmx + rmy*rmy);
	rsinL = binv*rsinL*sgn + 0.5*cosL;

	//获得升力面常数
	cop = 5.9;
	b0 = 8.88 - 1.88*angLS / 90;
	b1 = b1ls.interplinear_fast(angls, angLS);
	b2 = b2ls.interplinear_fast(angls, angLS);
	aop = aopls.interplinear_fast(angls, angLS);
	if (cosL < 0)
		a2 = a2ls.interplinear_fast(angls, angLS);
	else
		a2 = 0.0084;
	if (sinL < 1)
	{
		a1 = a1ls.interplinear_fast(angls, angLS);
		c0 = c0ls.interplinear_fast(angls, angLS);
		c1 = c1ls.interplinear_fast(angls, angLS);
		c2 = c2ls.interplinear_fast(angls, angLS);
	}
	else
	{
		a1 = -0.434;
		c0 = 1.683;
		c1 = 1.417;
		c2 = 0.910;
	}
	//升力面载荷
	rs2 = rsinL*rsinL;
	h0 = (hLS + cop)*(hLS + cop);
	x0 = -rs2 + h0 + b0*b0;
	h1 = (hLS + c1)*(hLS + c1);
	x1 = (rsinL + b1)*(rsinL + b1);
	h2 = (hLS + c2)*(hLS + c2);
	xx = x1*x1;
	hh = h2*h2;
	xh = x1*h2;
	cs = cos(b2);
	sn = sin(b2);

	LLS = rsinL / (rs2 + (hLS + c0)*(hLS + c0)) - aop*b0*x0 / (x0*x0 + 4 * rs2*h0);
	LLS -= a1*(2 * (rsinL + b1)*(-x1 + 3 * h1)*cs - 2 * (hLS + c1)*(-3 * x1 + h1)*sn) / pow(x1 + h1, 3);
	LLS -= a2*(24 * (rsinL + b1)*(-xx + 10 * xh - 5 * hh)*cs - 24 * (hLS + c2)*(-5 * xx + 10 * xh - hh)*sn) / pow(x1 + h2, 5);

	//升力线常数
	b2 = 0;
	c0 = 1.571;
	b1 = b1ll.interplinear_fast(angll, angLS);
	aop = aopll.interplinear_fast(angll, angLS);
	a1 = a1ll.interplinear_fast(angll, angLS);
	a2 = a2ll.interplinear_fast(angll, angLS);
	//升力线载荷
	x1 = (rsinL + b1)*(rsinL + b1);
	xx = x1*x1;
	xh = x1*h2;
	LLL = rsinL / (rs2 + (hLS + c0)*(hLS + c0)) - aop*b0*x0 / (x0*x0 + 4 * rs2*h0);
	LLL -= a1*(2 * (rsinL + b1)*(-x1 + 3 * h1)) / pow(x1 + h1, 3);
	LLL -= a2*(24 * (rsinL + b1)*(-xx + 10 * xh - 5 * hh)) / pow(x1 + h2, 5);
	//printf("AZ = %f, factor = %f\n", DEG(az), LLS/LLL);
	if (fabs(LLL) >= DBL_EPSILON && Abs(LLS/LLL) < 1)
		return LLS / LLL;
	else
		return 1.0;

}


void Wake::ComputeIndVel(double &vx,double &vy,double&vz,const int iz,const int ir, const int iz2, const int ib, double ds, double *cp)
{
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp, az, temp_lambdx, temp_lambdy, temp_lambdi;

	temp_lambdx = temp_lambdy = temp_lambdi = 0;

	az = iz * 2 * PI / nf;
	xblade = cp[0];
	yblade = cp[1];
	zblade = cp[2];
	rc04 = rc0*rc0*rc0*rc0;

	ri0 = xblade - geometryathub(nk - 1, iz2, 0);
	rj0 = yblade - geometryathub(nk - 1, iz2, 1);
	rk0 = zblade - geometryathub(nk - 1, iz2, 2);
	r0len = Norm(ri0, rj0, rk0);

	for (int ik = nk - 2; ik >= 0; ik--)
	{
		if(!(isbond && ik==ir && ib==0))
		{
			ri1 = xblade - geometryathub(ik, iz2, 0);
			rj1 = yblade - geometryathub(ik, iz2, 1);
			rk1 = zblade - geometryathub(ik, iz2, 2);
			r1len = Norm(ri1, rj1, rk1);
			
			cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
			height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));

			rdot = Dot(ri0, rj0, rk0, ri1, rj1, rk1);
			geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
			if (fabs(height2) > DBL_EPSILON)
			{
				geofunc *= height2 / rootNewton(rc04 + height2*height2, height2, height2*1e-2);
				geofunc *= 0.5*(vortexstr(ik, iz2) + vortexstr(ik + 1, iz2));

				if (lscorr.enable && r1len*r1len + r0len*r0len + 2 * rdot < 400 * ds*ds)
					geofunc *= lscorr.LSCorrection(az, ds, rc0, ri0, rj0, rk0, ri1, rj1, rk1);
			}
			else
				geofunc = 0;

			_temp = 0.25 / PI * geofunc;

			temp_lambdx -= rcros[0] * _temp;
			temp_lambdy -= rcros[1] * _temp;
			temp_lambdi -= rcros[2] * _temp;

			ri0 = ri1;
			rj0 = rj1;
			rk0 = rk1;
			r0len = r1len;
		}
		else		
		{
			ri0 = ri1;
			rj0 = rj1;
			rk0 = rk1;
			r0len = r1len;
		}
	}
	vx = temp_lambdx;
	vy = temp_lambdy;
	vz = temp_lambdi;
}

void Wake::_computeIndVel(int nb, Matrix1<double> &ch)
{
	//ch为无量纲的弦长
	myTYPE temp_lambdi, temp_lambdx, temp_lambdy;
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp, az, _c;

	int iz, ik;
	double _dfi, df;

	rc04 = rc0 * rc0 * rc0 * rc0;
	
	for (int iz = nf - 1; iz >= 0; iz--)
	{
		az = iz * 2 * PI / nf;
		while (az > 2 * PI)
			az -= 2 * PI;
		while (az < 0)
			az += 2 * PI;
	
		for (int ir = ns - 1; ir >= 0; --ir)
		{
			temp_lambdi = 0;
			temp_lambdx = 0;
			temp_lambdy = 0;

			xblade = comppointathub(iz, ir, 0);
			yblade = comppointathub(iz, ir, 1);
			zblade = comppointathub(iz, ir, 2);

			_c = ch(ir);
			for (int ib = 0; ib < nb; ib++)
			{
				int iz2 = iz + nf / nb*ib;
				if (iz2 >= nf)
					iz2 -= nf;

				// balde 1 $ik == nk - 1$ element points to blade 1 
				ri0 = xblade - geometryathub(nk - 1, iz2, 0);
				rj0 = yblade - geometryathub(nk - 1, iz2, 1);
				rk0 = zblade - geometryathub(nk - 1, iz2, 2);
				r0len = Norm(ri0, rj0, rk0);
				for (int ik = nk - 2; ik >= 0; --ik)
				{
					ri1 = xblade - geometryathub(ik, iz2, 0);
					rj1 = yblade - geometryathub(ik, iz2, 1);
					rk1 = zblade - geometryathub(ik, iz2, 2);
					r1len = Norm(ri1, rj1, rk1);

					cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
					height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));

					rdot = Dot(ri0, rj0, rk0, ri1, rj1, rk1);
					geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
					geofunc *= height2 / rootNewton(rc04 + height2*height2, height2, height2*1e-2);
					geofunc *= 0.5*(vortexstr(ik, iz2) + vortexstr(ik + 1, iz2));

					if (r1len*r1len + r0len*r0len + 2 * rdot < 400 * _c*_c)
						geofunc *= lscorr.LSCorrection(az, _c, rc0, ri0, rj0, rk0, ri1, rj1, rk1);

					_temp = 0.25 / PI * geofunc;

					temp_lambdx -= rcros[0] * _temp;
					temp_lambdy -= rcros[1] * _temp;
					temp_lambdi -= rcros[2] * _temp;

					ri0 = ri1;
					rj0 = rj1;
					rk0 = rk1;
					r0len = r1len;
				}
			}
		}
	}
	
}

bool Wake::_computeVorStr(Matrix2<double> &cirlb)
{
	int igen, istr;
	Matrix1<int> idns;

	istr = strboard*ns;
	if (istip)
	{
		if (istr < 0)
			istr = 0; 
		idns.allocate(ns - istr);		
		idns = step(istr, ns - 1);
		for (int j = nf - 1; j >= 0; --j)
		{
			for (int i = nk - 1; i >= 0; --i)
			{
				igen = j - i;
				while (igen < 0)
					igen += nf;
				vortexstr(i, j) = cirlb(igen, idns).findmax();
			}
		}
	}
	else if(isrot)
	{
		if (istr > ns)
			istr = ns; 
		idns.allocate(istr);	
		idns = step(0, istr - 1);
		for (int j = nf - 1; j >= 0; --j)
		{
			for(int i=nk-1;i>=0;--i)
			{
				igen = j - i;
				while (igen<0)
					igen += nf;
				vortexstr(i, j) = -cirlb(igen, idns).findmax();
			}
		}
	}
	else if(isbond)
	{
		for (int j = nf - 1; j >= 0; --j)
			for (int i = nk - 1; i >= 0; --i)
				vortexstr(i, j) = cirlb(j, Max(0, Min(i, ns - 1)));
	}
	else if (isshed)
	{
		return false;
	}
	return true;
}

bool Wake::_computeGeometry(const AeroDynType aero, double *veltpp, double wi, double b0)
{
	//wi和veltpp传入的是无量纲的速度
	//有后掠角时还要调整
	myTYPE r0, z0, x0, xv, yv, zv, xe, df;
	myTYPE ka, E, kx, ky, ky3, k0;
	myTYPE _dfi;

	if (istip || isrot)
	{
		r0 = rv*cos(b0);
		ka = atan2(-veltpp[0], -(-veltpp[2] + wi));
		E = ka*0.5;
		kx = E; ky = 2 * veltpp[0]; ky3 = -kx; k0 = 1.0 - 8.0 * ky3 / 15.0 / PI;
		df = 2.0 * PI / nf; xe = 0;

		for (int j = nf - 1; j >= 0; --j)
		{
			for (int i = nk - 1; i >= 0; --i)
			{
				_dfi = df*(j - i);
				x0 = r0*cos(_dfi) + chv*0.75*sin(_dfi);
				xv = x0 - veltpp[0] * df*i;
				yv = r0 * sin(_dfi) - chv*0.75*cos(_dfi);
				if (xv * xv + yv * yv < 1.01)
				{
					if(istip)
						zv = -wi * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * i*df - wi / 2.0 * kx * (x0 + xv) * i*df;
					else if(isrot)
						zv = -wi * (k0 + ky3 * Abs(pow(1-yv,3)) + ky * yv) * i*df - wi / 2.0 * kx * (x0 + xv) * i*df;
				}
				else 
				{
					if (yv > 1)
						xe = 0;
					else
						xe = sqrt(1.0 - yv*yv);
					if (istip)
					{
						zv = -wi / (-veltpp[0])  * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * (xe - x0) - wi / 2.0 * kx / (-veltpp[0])  * (xe*xe - x0*x0);
						zv -= 2.0 * wi / (-veltpp[0]) * (k0 + ky3 * Abs(yv*yv*yv) + ky * yv) * (xv - xe);
					}
					else if (isrot)
					{
						zv = -wi / (-veltpp[0])  * (k0 + ky3 * Abs(pow(1 - yv, 3)) + ky * yv) * (xe - x0) - wi / 2.0 * kx / (-veltpp[0])  * (xe*xe - x0*x0);
						zv -= 2.0 * wi / (-veltpp[0]) * (k0 + ky3 * Abs(pow(1 - yv, 3)) + ky * yv) * (xv - xe);
					}
				}
				geometry(i, j, 0) = xv;
				geometry(i, j, 1) = yv;
				geometry(i, j, 2) = -zv - veltpp[2] * df*i; 
			}
		}
		geoexp = geometry.reshape(3, 3, nk*nf);
	}
	else if(isbond)
	{
		for (int j = nf - 1; j >= 0; --j)
		{
			for (int i = nk - 1; i >= 0; --i)
			{
				geometry(i, j, 0) = rv + i * (1.0 - rv) / (nk - 1);
				geometry(i, j, 1) = 0;
				geometry(i, j, 2) = 0;
			}
		}
	}
	else if(isshed)
	{
		return false;
	}
	
	return true;
}

void Wake::_airloadingPoint(Matrix2<double> &az,Matrix2<double> &ra,Matrix1<double> &ch, Matrix1<double> tw, double *b, double *s, double offset, double sr)
{
	//尾迹法速度计算点，目前取四分之三弦长位置
	//有后掠角时还要调整
	myTYPE azm, bfp, rst, chp, twp;
	for (int j = ns - 1; j >= 0; --j)
	{
		chp = ch(j);
		twp = tw(j) + sr;
		for (int i = nf - 1; i >= 0; --i)
		{
			azm = az(i, j);
			bfp = b[0] + b[1] * cos(azm) + b[2] * sin(azm);
			rst = ra(i, j);

			comppoint(i, j, 0) = ((rst - offset)*cos(bfp) + offset) * cos(azm) + chv*0.5*sin(azm);
			comppoint(i, j, 1) = ((rst - offset)*cos(bfp) + offset) * sin(azm) - chv*0.5*cos(azm);
			comppoint(i, j, 2) = bfp*(rst - offset) - 0.5*chp*(s[0] + s[1] * cos(azm) + s[2] * sin(azm) + twp);
		}
	}
}
