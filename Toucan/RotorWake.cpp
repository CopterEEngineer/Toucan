#include "stdafx.h"
#include "Components.h"


void Rotor::_wakeInducedVel(void)
{
	switch (adyna)
	{
	case PWake:
		haveStr = true;// _tipVortexStr();
		haveGeo = _wakeGeoBd();
		break;
	case FWake:
		break;
	default:
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
	myTYPE ch = chord(0) / radius;
	myTYPE twistt = twist(ns - 1) - twist(0) + pitchroot;
	for (int j = ns - 1; j >= 0; --j)
	{
		for (int i = nf - 1; i >= 0; --i)
		{
			azm = azstation(i, j);
			bfp = beta[0] + beta[1] * cos(azm) + beta[2] * sin(azm);
			rst = rastation(i, j);
			bladedeform(i, j, 0) = ((rst - eflap)*cos(bfp) + eflap) * cos(azm) + ch*0.5*sin(azm);
			bladedeform(i, j, 1) = ((rst - eflap)*cos(bfp) + eflap) * sin(azm) - ch*0.5*cos(azm);
			bladedeform(i, j, 2) = bfp*(rst - eflap) - 0.5*ch*(sita[0] + sita[1] * cos(azm) + sita[2] * sin(azm) + twistt);
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
	lambdi.output("DEBUG_lambdi.output", 10);
}

void Rotor::_wakeInducedVelMP()
{
	// 适用于nb = 2
	myTYPE temp_lambdi, temp_lambdx, temp_lambdy;
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	Matrix2<myTYPE> tipgeoexpand, tempM(3, 3);
	Matrix3<myTYPE> tipgeoathub(nk, nf, 3);
	int tipstr_NI = tipstr.NI;
	int tipgeo_NI = nk;
	int tipgeo_NJ = nf;
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp;


	/*if (nb != 2)
	{
		printf("Undefined nb != 2 in _wakeIndVelCalc() Func. \n");
		system("pause");
		return;
	}*/

	rc04 = rc0 * rc0 * rc0 * rc0;
	tipgeoexpand = tipgeometry.reshape(3, 3, nk*nf); // 应该要减去在TPP坐标系下TPP原点指向HUB原点的向量
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			tempM(i, j) = tppcoord.Ttransf[i][j];
	tipgeoexpand = tempM.transpose().matrixmultiplyTP(tipgeoexpand);
	for (int k = 0; k < 3; ++k)
		for (int j = 0; j < nk*nf; ++j)
			tipgeoathub(j%nk, j / nk, k) = tipgeoexpand(k, j);

	//int iz = nf - 1, iz2 = iz + nf / nb;
//#pragma omp parallel num_threads(4) 
	{
//#pragma omp for 
		for (int iz = nf - 1; iz >= 0; iz--)
		{
			int iz2 = iz + nf / nb;
			if (iz2 >= nf)
				iz2 -= nf;

			for (int ir = ns - 1; ir >= 0; --ir)
			{
				temp_lambdi = 0;
				temp_lambdx = 0;
				temp_lambdy = 0;
				xblade = bladedeform(iz, ir, 0);
				yblade = bladedeform(iz, ir, 1);
				zblade = bladedeform(iz, ir, 2);

				// balde 1 $ik == nk - 1$ element points to blade 1 
				ri0 = xblade - tipgeoathub(nk - 1, iz, 0);
				rj0 = yblade - tipgeoathub(nk - 1, iz, 1);
				rk0 = zblade - tipgeoathub(nk - 1, iz, 2);
				r0len = norm(ri0, rj0, rk0);
				for (int ik = nk - 2; ik >= 0; --ik)
				{
					ri1 = xblade - tipgeoathub(ik, iz, 0);
					rj1 = yblade - tipgeoathub(ik, iz, 1);
					rk1 = zblade - tipgeoathub(ik, iz, 2);
					r1len = norm(ri1, rj1, rk1);

					cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
					height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));
					rdot = dot(ri0, rj0, rk0, ri1, rj1, rk1);
					geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
					geofunc *= (1 - 0.5 * rc04 / height2 / height2);
					geofunc *= 0.5*(tipstr(ik, iz) + tipstr(ik + 1, iz));

					_temp = 0.25 / PI * geofunc;

					temp_lambdx -= rcros[0] * _temp;
					temp_lambdy -= rcros[1] * _temp;
					temp_lambdi -= rcros[2] * _temp;

					ri0 = ri1;
					rj0 = rj1;
					rk0 = rk1;
					r0len = r1len;
				}

				// balde 2 $ik == nk - 1$ element points to blade 1 
				ri0 = xblade - tipgeoathub(nk - 1, iz2, 0);
				rj0 = yblade - tipgeoathub(nk - 1, iz2, 1);
				rk0 = zblade - tipgeoathub(nk - 1, iz2, 2);
				r0len = norm(ri0, rj0, rk0);
				for (int ik = nk - 2; ik >= 0; --ik)
				{
					ri1 = xblade - tipgeoathub(ik, iz2, 0);
					rj1 = yblade - tipgeoathub(ik, iz2, 1);
					rk1 = zblade - tipgeoathub(ik, iz2, 2);
					r1len = norm(ri1, rj1, rk1);

					cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
					height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));
					rdot = dot(ri0, rj0, rk0, ri1, rj1, rk1);
					geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
					geofunc *= (1 - 0.5 * rc04 / height2 / height2);
					geofunc *= 0.5*(tipstr(ik, iz2) + tipstr(ik + 1, iz2));

					_temp = 0.25 / PI * geofunc;

					temp_lambdx -= rcros[0] * _temp;
					temp_lambdy -= rcros[1] * _temp;
					temp_lambdi -= rcros[2] * _temp;

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
	}
	lambdh = lambdi - vel[2] / vtipa;
	lambdi.output("DEBUG_lambdi.output", 10);
}

void Rotor::_wakeInducedVelMP(int nb)
{
	myTYPE temp_lambdi, temp_lambdx, temp_lambdy;
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	Matrix2<myTYPE> tipgeoexpand, tempM(3, 3);
	Matrix3<myTYPE> tipgeoathub(nk, nf, 3);
	int tipstr_NI = tipstr.NI;
	int tipgeo_NI = nk;
	int tipgeo_NJ = nf;
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp;

	rc04 = rc0 * rc0 * rc0 * rc0;
	tipgeoexpand = tipgeometry.reshape(3, 3, nk*nf); 
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			tempM(i, j) = tppcoord.Ttransf[i][j];
	tipgeoexpand = tempM.transpose().matrixmultiplyTP(tipgeoexpand);
	
	int iz, ik;
	double _dfi, df;
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
				tipgeoathub(iz, ik, 2) += (rtip - eflap)*(beta[0] + beta[1] * cos(_dfi) + beta[2] * sin(_dfi));
				tipgeoathub(iz, ik, 2) -= 0.75*chv*(sita[0] + sita[1] * cos(_dfi) + sita[2] * sin(_dfi) + twistv - twist(0) + pitchroot);
			}
		}
	}

	for (int iz = nf - 1; iz >= 0; iz--)
	{	
		for (int ir = ns - 1; ir >= 0; --ir)
		{
			temp_lambdi = 0;
			temp_lambdx = 0;
			temp_lambdy = 0;

			xblade = bladedeform(iz, ir, 0);
			yblade = bladedeform(iz, ir, 1);
			zblade = bladedeform(iz, ir, 2);

			for (int ib = 0; ib < nb; ib++)
			{
				int iz2 = iz + nf / nb*ib;
				if (iz2 >= nf)
					iz2 -= nf;

				// balde 1 $ik == nk - 1$ element points to blade 1 
				ri0 = xblade - tipgeoathub(nk - 1, iz2, 0);
				rj0 = yblade - tipgeoathub(nk - 1, iz2, 1);
				rk0 = zblade - tipgeoathub(nk - 1, iz2, 2);
				r0len = norm(ri0, rj0, rk0);
				for (int ik = nk - 2; ik >= 0; --ik)
				{
					ri1 = xblade - tipgeoathub(ik, iz2, 0);
					rj1 = yblade - tipgeoathub(ik, iz2, 1);
					rk1 = zblade - tipgeoathub(ik, iz2, 2);
					r1len = norm(ri1, rj1, rk1);

					cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
					height2 = (rcros[0] * rcros[0] + rcros[1] * rcros[1] + rcros[2] * rcros[2]) / ((ri0 - ri1)*(ri0 - ri1) + (rj0 - rj1)*(rj0 - rj1) + (rk0 - rk1)*(rk0 - rk1));
					rdot = dot(ri0, rj0, rk0, ri1, rj1, rk1);
					geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
					geofunc *= height2 / rootNewton(rc04 + height2*height2, height2, height2*1e-2);
					//geofunc *= height2 / sqrt(rc04 + height2*height2);
					geofunc *= 0.5*(tipstr(ik, iz2) + tipstr(ik + 1, iz2));

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

