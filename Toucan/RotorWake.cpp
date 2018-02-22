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


void Rotor::_bladePosition()
{
	//尾迹法速度计算点，目前取四分之三弦长位置
	//有后掠角时还要调整
	myTYPE azm, bfp, rst;
	myTYPE ch;
	myTYPE sth;
	for (int j = ns - 1; j >= 0; --j)
	{
		ch = chord(j) / radius;
		for (int i = nf - 1; i >= 0; --i)
		{
			azm = azstation(i, j);
			rst = rastation(i, j);

			bfp = beta[0] + beta[1] * cos(azm) + beta[2] * sin(azm);
			sth = twist(j) - twist(0) + pitchroot + sita[0] + sita[1] * cos(azm) + sita[2] * sin(azm);

			bladedeform(i, j, 0) = ((rst - eflap)*cos(bfp) + eflap) * cos(azm) +ch*0.5*sin(azm);
			bladedeform(i, j, 1) = ((rst - eflap)*cos(bfp) + eflap) * sin(azm) - ch*0.5*cos(azm);
			bladedeform(i, j, 2) = bfp*(rst - eflap) -0.5*ch*sth;

			tailposition(i, j, 0) = bladedeform(i, j, 0) + ch*0.25*sin(azm);
			tailposition(i, j, 1) = bladedeform(i, j, 1) - ch*0.25*cos(azm);
			tailposition(i, j, 2) = bladedeform(i, j, 2) - ch*0.25*sth;
		}
	}
	//if (bladedeform(1, ns - 1, 1) < 0)
	//bladedeform.output("DEBUG_bladedeform.output", 6);
	//tailposition.output("DEBUG_tailposition.output", 6);
}


void Rotor::_wakeInducedVelMP(void)
{
	double temp_lambdx, temp_lambdy, temp_lambdi;
	double _lambdx, _lambdy, _lambdi;
	double _drc2, a1, _nu, alpha, _delta;

	alpha = 1.25643;
	_nu = 14.8*1e-6;
	a1 = 4e-4;
	_delta= 1 + a1 * 2 * omega*radius*chord(0) / _nu*Abs(airforce[2]) / amb.rho / disk_A / vtipa / vtipa / sigma;
	_drc2 = Abs(4 * alpha*_nu*_delta / omega);


	for (int i = 0; i<wakev.size();i++)
	{
		int iz, ik, nk, nf, ns;
		double _dfi, bfp1, bfp2, sth;

		nk = wakev[i].nk;
		nf = wakev[i].nf;
		ns = wakev[i].ns;

		if (wakev[i].istip || wakev[i].isrot)
		{
			double dx, dy, dz, _az;
			double _vpt[3], _vph[3], _rv, sth1, sth2;
			for (int iz = 0; iz < nf; iz++)
			{
				_az = 2 * PI / nf*iz;
				bfp2 = beta[0] + beta[1] * cos(_az) + beta[2] * sin(_az);
				sth2 = sita[0] + sita[1] * cos(_az) + sita[2] * sin(_az) + wakev[i].twistv - twist(0) + pitchroot;
				_rv = (wakev[i].rv - eflap)*cos(bfp2) + eflap;
				//当前时刻脱落位置坐标，在hub坐标系中
				_vph[0] = _rv*cos(_az) + wakev[i].chv*0.75*sin(_az);
				_vph[1] = _rv*sin(_az) - wakev[i].chv*0.75*cos(_az);
				_vph[2] = bfp2*(wakev[i].rv - eflap) - wakev[i].chv*0.75*sth2;
				//当前时刻脱落位置坐标，在TPP坐标系中
				_vpt[0] = wakev[i].geometry(0, iz, 0);
				_vpt[1] = wakev[i].geometry(0, iz, 1);
				_vpt[2] = wakev[i].geometry(0, iz, 2) - wakev[i].chv*0.5*sth2;

				for (int ik = 0; ik < nk; ik++)
				{
					_dfi = 2 * PI / nf*(iz - ik); //k脱落时的方位角
					bfp1 = beta[0] + beta[1] * cos(_dfi) + beta[2] * sin(_dfi);
					sth1 = sita[0] + sita[1] * cos(_dfi) + sita[2] * sin(_dfi) + wakev[i].twistv - twist(0) + pitchroot;
					dx = wakev[i].geometry(ik, iz, 0) - _vpt[0];
					dy = wakev[i].geometry(ik, iz, 1) - _vpt[1];
					//脱落位置从”尾迹相对桨盘平面“修正到后缘对桨盘的相对位置，每个k脱落时的桨盘平面位置也不同，
					//修正后_z是相对于现在桨盘坐标系的尾迹
					dz = wakev[i].geometry(ik, iz, 2) - wakev[i].chv*0.5*sth1 + wakev[i].rv*(bfp1 - bfp2);
					dz -= _vpt[2];

					wakev[i].geometryathub(ik, iz, 0) = cos(beta[1])*dx - sin(beta[1])*sin(beta[2])*dy + cos(beta[2])*sin(beta[1])*dz + _vph[0];
					wakev[i].geometryathub(ik, iz, 1) = cos(beta[2])*dy + sin(beta[2])*dz + _vph[1];
					wakev[i].geometryathub(ik, iz, 2) = -sin(beta[1])*dx - cos(beta[1])*sin(beta[2])*dy + cos(beta[1])*cos(beta[2])*dz + _vph[2];

				}
			}
		}
		else if (wakev[i].isbond)
		{
			for (int iz = 0; iz < nf; iz++)
			{  
				_dfi = azstation(iz, 0);
				bfp1 = (beta[0] + beta[1] * cos(_dfi) + beta[2] * sin(_dfi));
				for (int ik = 0; ik < nk; ik++)
				{
					wakev[i].geometryathub(ik, iz, 0) = ((wakev[i].geometry(ik, iz, 0) - eflap)* cos(bfp1) + eflap)* cos(_dfi);
					wakev[i].geometryathub(ik, iz, 1) = ((wakev[i].geometry(ik, iz, 0) - eflap)* cos(bfp1) + eflap)* sin(_dfi);
					wakev[i].geometryathub(ik, iz, 2) = bfp1*(wakev[i].geometry(ik, iz, 0) - eflap);
				}
			}
		}
	}

	double cp[3] = { 0,0,0 };
	for (int iw = wakev.size() - 1; iw >= 0; --iw)
	{
		wakev[iw].lambdx.setvalue(0);
		wakev[iw].lambdy.setvalue(0);
		wakev[iw].lambdz.setvalue(0);
	}
	for (int iz = nf - 1; iz >= 0; iz--)
	{
		//cout << iz << endl;
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
					wakev[iw].ComputeIndVel(_lambdx, _lambdy, _lambdi, iz, ir, iz2, ib, chord(ir) / radius, cp, _drc2);
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
	//wakev[0].geometry.output("DEBUG_tipgeo.output", 10);
	//wakev[0].geometryathub.output("DEBUG_tipgeo_athub.output", 10);
}


void Wake::ComputeIndVel(double &vx,double &vy,double&vz,const int iz,const int ir, const int iz2, const int ib, double ds, double *cp, double drc2)
{
	myTYPE xblade, yblade, zblade, ri0, rj0, rk0, ri1, rj1, rk1;
	myTYPE r0len, r1len, rdot, height, height2, rc04, geofunc;
	myTYPE rcros[3], _temp, az, temp_lambdx, temp_lambdy, temp_lambdi;
	double _rc2;

	temp_lambdx = temp_lambdy = temp_lambdi = 0;

	az = iz * 2 * PI / nf;
	xblade = cp[0];
	yblade = cp[1];
	zblade = cp[2];
	rc04 = rc0*rc0*rc0*rc0;
	_rc2 = rc0*rc0;

	ri1 = xblade - geometryathub(0, iz2, 0);
	rj1 = yblade - geometryathub(0, iz2, 1);
	rk1 = zblade - geometryathub(0, iz2, 2);

	r1len = Norm(ri1, rj1, rk1);

	for (int ik = 1; ik < nk; ik++)
	{
		ri0 = xblade - geometryathub(ik, iz2, 0);
		rj0 = yblade - geometryathub(ik, iz2, 1);
		rk0 = zblade - geometryathub(ik, iz2, 2);
		r0len = Norm(ri0, rj0, rk0);

		_rc2 += drc2 * 2 * PI / nf;

		if (!(isbond && ik == ir + 1 && ib == 0))
		{
			cross(rcros, ri0, rj0, rk0, ri1, rj1, rk1);
			height2 = pow(ri0 - ri1, 2) + pow(rj0 - rj1, 2) + pow(rk0 - rk1, 2);
			if (height2 > 0)
				height2 = (pow(rcros[0], 2) + pow(rcros[1], 2) + pow(rcros[2], 2)) / height2;
			else
				height2 = 0;

			rdot = Dot(ri0, rj0, rk0, ri1, rj1, rk1);
			geofunc = -(1.0 / r0len + 1.0 / r1len) / (rdot + r0len * r1len);
			if (fabs(height2) > DBL_EPSILON)
			{
				geofunc *= height2 / rootNewton(_rc2 + height2*height2, height2, height2*1e-2);

				if (lscorr.enable && r1len*r1len + r0len*r0len + 2 * rdot < 400 * ds*ds)				
					geofunc *= lscorr.LSCorrection(az, ds, _rc2, ri0, rj0, rk0, ri1, rj1, rk1);

				geofunc *= 0.5*(vortexstr(ik - 1, iz2) + vortexstr(ik, iz2));

			}
			else
				geofunc = 0;

			_temp = 0.25 / PI * geofunc;

			temp_lambdx -= rcros[0] * _temp;
			temp_lambdy -= rcros[1] * _temp;
			temp_lambdi -= rcros[2] * _temp;

			ri1 = ri0;
			rj1 = rj0;
			rk1 = rk0;
			r1len = r0len;
		}
		else		
		{
			ri1 = ri0;
			rj1 = rj0;
			rk1 = rk0;
			r1len = r0len;
		}
	
		rc(ik) = _rc2;
	
	}
	vx = temp_lambdx;
	vy = temp_lambdy;
	vz = temp_lambdi;
	lambdx(iz, ir) += vx;
	lambdy(iz, ir) += vy;
	lambdz(iz, ir) += vz;
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
						zv = -wi * (1 + ky * yv) * i*df - wi / 2.0 * kx * (x0 + xv) * i*df;
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
						zv -= 2.0 * wi / (-veltpp[0]) * (1 + ky * yv) * (xv - xe);
					}
					else if (isrot)
					{
						zv = -wi / (-veltpp[0])  * (k0 + ky * yv) * (xe - x0) - wi / 2.0 * kx / (-veltpp[0])  * (xe*xe - x0*x0);
						zv -= 2.0 * wi / (-veltpp[0]) * (1 + ky * yv) * (xv - xe);
					}
				}
				geometry(i, j, 0) = xv;
				geometry(i, j, 1) = yv;
				geometry(i, j, 2) = -zv - veltpp[2] * df*i; 
			}
		}
		//geoexp = geometry.reshape(3, 3, nk*nf);
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



double LSCorr::LSCorrection(double az, double c, double rc2, double ri0, double rj0, double rk0, double ri1, double rj1, double rk1)
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
	if (rc2 > 0)
		hLS += rc2;
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
	b1 = 0.5*b1;//b1ll.interplinear_fast(angls, angLS);// 0.5*b1;
	aop = aop*0.75; //aopll.interplinear_fast(angls, angLS);// aop*0.75;
	a1 *= (1.25 + 0.5*sinL);
	a2 *= (2.5 + sinL);
	//升力线载荷
	x1 = (rsinL + b1)*(rsinL + b1);
	xx = x1*x1;
	xh = x1*h2;
	LLL = rsinL / (rs2 + (hLS + c0)*(hLS + c0)) - aop*b0*x0 / (x0*x0 + 4 * rs2*h0);
	LLL -= a1*(2 * (rsinL + b1)*(-x1 + 3 * h1)) / pow(x1 + h1, 3);
	LLL -= a2*(24 * (rsinL + b1)*(-xx + 10 * xh - 5 * hh)) / pow(x1 + h2, 5);
	//printf("AZ = %f, factor = %f\n", DEG(az), LLS/LLL);
	if (fabs(LLL) >= DBL_EPSILON && Abs(LLS / LLL) < 1)
		return LLS / LLL;
	else
		return 1.0;

}
