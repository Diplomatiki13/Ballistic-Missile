#include "pch.h"
#include "missile.h"
#include <cmath>
#include <iostream>
#include <fstream>

void Missile::CalcGraph(ofstream & file)
{
	long double Dangle = sqrt((rStart.Phi()-rEnd.Phi())*(rStart.Phi()-rEnd.Phi())+(rStart.Theta()-rEnd.Theta())*(rStart.Theta()-rEnd.Theta()));
	Vector Rc = rEnd + ((rStart - rEnd) / 2.0);
	Vector init(Rc.Phi() - ((3.0/4.0)*Dangle), Rc.Theta() - ((3.0 / 4.0)*Dangle) , 6371000.0, 0.0);
	file << "ex = " << init.GetX() << ";\n";
	file << "ey = " << init.GetY() << ";\n";
	file << "ez = " << init.GetZ() << ";\n";
	for (long double i = Rc.Phi() - ((3.0 / 4.0)*Dangle); i <= Rc.Phi() + ((3.0 / 4.0)*Dangle); i += (1.5*Dangle) / 200)
	{
		for (long double j = Rc.Theta() - ((3.0 / 4.0)*Dangle); j <= Rc.Theta() + ((3.0 / 4.0)*Dangle); j += (1.5*Dangle) / 200)
		{
			if (i == Rc.Phi() - ((3.0 / 4.0)*Dangle) && j == Rc.Theta() - ((3.0 / 4.0)*Dangle))
			{
			}
			else
			{
				Vector temp(i, j, 6371000.0, 0.0);
				file << "ex = [ ex , " << temp.GetX() << " ];\n";
				file << "ey = [ ey , " << temp.GetY() << " ];\n";
				file << "ez = [ ez , " << temp.GetZ() << " ];\n";
			}
		}
	}
}

void Missile::Initialise(long double lt, long double lp, long double t, long double p, long double dt,ofstream& file)
{
	fuelmass = fuelmassfull;
	inhibitor = 0;
	hinhibitor = 0;
	command = 0;
	ctemp = 1;
	if (dt <= 0.1)
	{
		gtemp = round(1.0 / (10.0*dt));
	}
	else
	{
		gtemp = 1;
	}
	lpInit = lp;
	ltInit = lt;
	launchPhiInit = p;
	launchThetaInit = t;
	Vector newrMissile(Earth::GetR()*cos(lt)*cos(lp), Earth::GetR()*sin(lt)*cos(lp), Earth::GetR()*sin(lp));
	rMissile = newrMissile;
	rStart = newrMissile;

	Vector vInit(0.0, 0.0, 0.0);
	vMissile = vInit;

	Vector newB(-cos(lt)*cos(lp), -sin(lt)*cos(lp), -sin(lp));
	B = newB*((framemass+fuelmass)*Earth::GetGr(0.0));

	long double vx = cos(pi/2.0)*cos(0.0);
	long double vy = cos(pi/2.0)*sin(0.0);
	long double vz = sin(pi/2.0);
	long double Vx = vx * cos(pi/2.0 + lt) + vy * sin(lp)*cos(pi + lt) + vz * cos(lp)*cos(lt);
	long double Vy = vx * sin(pi/2.0 + lt) + vy * sin(lp)*sin(pi + lt) + vz * cos(lp)*sin(lt);
	long double Vz = vy * cos(lp) + vz * sin(lp);
	Vector newT(Vx,Vy,Vz);
	T = newT * (Earth::GetGr(0.0)*isp*massPerS);

	Vector newD(-Vx, -Vy, -Vz);
	RefS = pi * missileDiameter;
	MaxS = pi * missileDiameter*(missileLength - missileLfore) + pi * missileR*(missileR + sqrt(missileR*missileR + missileLfore * missileLfore));
	reynolds = 0.0;
	mach = 0.0;
	Cdt = 0.0;
	Cda = 0.0;
	Cdw = 0.0;
	Cd = Cdt + Cda + Cdw;
	D = newD * ((0.5)*Earth::GetDen(0.0,hinhibitor)*vMissile.Magnitude()*vMissile.Magnitude()*RefS*Cd); 

	RF = (B + T) + D;
	aMissile = RF * (1.0 / (framemass + fuelmass));
	fuelmass -= massPerS*dt;

	file << "rx = " << rMissile.GetX() << ";\n";
	file << "ry = " << rMissile.GetY() << ";\n";
	file << "rz = " << rMissile.GetZ() << ";\n";
}

void Missile::updateForces(long double dt)
{
	if (height > 100 && inhibitor == 0)
	{
		long double vx = cos(launchPhiInit)*cos(launchThetaInit);
		long double vy = cos(launchPhiInit)*sin(launchThetaInit);
		long double vz = sin(launchPhiInit);
		long double Vx = vx * cos(pi / 2.0 + ltInit) + vy * sin(lpInit)*cos(pi + ltInit) + vz * cos(lpInit)*cos(ltInit);
		long double Vy = vx * sin(pi / 2.0 + ltInit) + vy * sin(lpInit)*sin(pi + ltInit) + vz * cos(lpInit)*sin(ltInit);
		long double Vz = vy * cos(lpInit) + vz * sin(lpInit);
		Vector temp(Vx*vMissile.Magnitude(), Vy*vMissile.Magnitude(), Vz*vMissile.Magnitude());
		vMissile = temp;
		inhibitor = 1;
		printf("-Turn to Initial Angles\n");
	}

	Vector newT = vMissile * (1.0 / vMissile.Magnitude());

	if (fuelmass > 0.0)
	{
		T = newT * (Earth::GetGr(0.0)*isp*massPerS);
	}
	else
	{
		T = T * 0.0;
		if (fuelmass < 0.0)
		{
			printf("-Ballistic Trajectory Start\n");
		}
		fuelmass = 0.0;
	}

	Vector newB = rMissile * ((-1.0) / rMissile.Magnitude());
	B = newB * ((framemass + fuelmass)*Earth::GetGr(height));

	Vector newD = newT * (-1.0);

	reynolds = (Earth::GetDen(height, hinhibitor)*vMissile.Magnitude()*missileLength) / Earth::GetVisc(height);
	if (height > 86000)
	{
		hinhibitor = 1;
	}

	mach = vMissile.Magnitude() / sqrt((1.4*8.3144598*Earth::GetTemp(height)) / 0.0289644);
	if (mach > maxmach)
	{
		maxmach = mach;
	}

	Cdt = (0.032 / pow(reynolds, 0.145)) * pow((1.0 + 0.10*mach*mach), -0.5) * (MaxS / RefS);

	if (mach >= 1.0)
	{
		Cdw = pow(missileDiameter / (2.0*missileLfore), 2.0)*
			(1.595 - (8.0*missileLfore*pow((missileDiameter / missileLength), 2.0)*sqrt(mach*mach - 1)));
	}
	else
	{
		Cdw = 0.0;
	}
	if (fuelmass == 0.0)
	{
		if (mach <= 0.6)
		{
			Cda = 0.029/(sqrt(Cdt));
		}
		else if (mach > 0.6)
		{
			long double vMagn0 = 0.6 * sqrt((1.4*8.3144598*Earth::GetTemp(height)) / 0.0289644);
			long double reynolds0 = (Earth::GetDen(height, hinhibitor)*vMagn0*missileLength) / Earth::GetVisc(height);
			long double Cdt0 = (0.032 / pow(reynolds0, 0.145)) * pow((1.0 + 0.10*0.6*0.6), -0.5) * (MaxS / RefS);

			if (mach > 0.6 && mach <= 1.0)
			{
				fb = 1.0 + 215.8*pow((mach - 0.6), 6.0);
			}
			else if (mach > 1.0 && mach <= 2.0)
			{
				fb = 2.0881*pow((mach - 1.0), 3.0) - 3.7938*pow((mach - 1.0), 2.0) + 1.4618*(mach - 1.0) + 1.883917;
			}
			else
			{
				fb = 0.297*pow((mach - 2.0), 3.0) - 0.7932*pow((mach - 2.0), 2.0) - 0.1115*(mach - 2.0) + 1.64006;
			}
			Cda = Cdt0 * fb;
		}
		else
		{
			printf("Negative mach");
		}
	}
	else
	{
		Cda = 0.0;
	}

	Cd = Cdt  + Cdw + Cda;
	D = newD * ((0.5)*Earth::GetDen(height, hinhibitor)*vMissile.Magnitude()*vMissile.Magnitude()*RefS*Cd);

	RF = (B + T) + D;
	aMissile = RF * (1.0 / (framemass + fuelmass));
	if (fuelmass > 0.0)
	{
		fuelmass -= massPerS * dt;
	}
}

void Missile::updatePosAndVel(long double h, ofstream& file)
{
	long double k0 = h * f(vMissile.GetX());
	long double l0 = h * g(vMissile.GetY());
	long double m0 = h * i(rMissile.GetX(), vMissile.GetY());
	long double n0 = h * j(rMissile.GetY(), vMissile.GetX());
	long double k1 = h * f(vMissile.GetX() + 0.5*m0);
	long double l1 = h * g(vMissile.GetY() + 0.5*n0);
	long double m1 = h * i(rMissile.GetX() + 0.5*k0, vMissile.GetY() + 0.5*n0);
	long double n1 = h * j(rMissile.GetY() + 0.5*l0, vMissile.GetX() + 0.5*m0);
	long double k2 = h * f(vMissile.GetX() + 0.5*m1);
	long double l2 = h * g(vMissile.GetY() + 0.5*n1);
	long double m2 = h * i(rMissile.GetX() + 0.5*k1, vMissile.GetY() + 0.5*n1);
	long double n2 = h * j(rMissile.GetY() + 0.5*l1, vMissile.GetX() + 0.5*m1);
	long double k3 = h * f(vMissile.GetX() + m2);
	long double l3 = h * g(vMissile.GetY() + n2);
	long double m3 = h * i(rMissile.GetX() + k2, vMissile.GetY() + n2);
	long double n3 = h * j(rMissile.GetY() + l2, vMissile.GetX() + m2);

	long double nx = rMissile.GetX() + (1.0 / 6.0)*(k0 + (2.0 * k1) + (2.0 * k2) + k3);
	long double ny = rMissile.GetY() + (1.0 / 6.0)*(l0 + (2.0 * l1) + (2.0 * l2) + l3);
	long double nvx = vMissile.GetX() + (1.0 / 6.0)*(m0 + (2.0 * m1) + (2.0 * m2) + m3);
	long double nvy = vMissile.GetY() + (1.0 / 6.0)*(n0 + (2.0 * n1) + (2.0 * n2) + n3);

	long double o0 = h * d(vMissile.GetZ());
	long double p0 = h * e();
	long double o1 = h * d(vMissile.GetZ() + 0.5*p0);
	long double p1 = h * e();
	long double o2 = h * d(vMissile.GetZ() + 0.5*p1);
	long double p2 = h * e();
	long double o3 = h * d(vMissile.GetZ() + p2);
	long double p3 = h * e();

	long double nz = rMissile.GetZ() + (1.0 / 6.0)*(o0 + (2.0 * o1) + (2.0 * o2) + o3);
	long double nvz = vMissile.GetZ() + (1.0 / 6.0)*(p0 + (2.0 * p1) + (2.0 * p2) + p3);

	Vector newR(nx, ny, nz);
	Vector newV(nvx, nvy, nvz);

	rMissile = newR;
	vMissile = newV;

	if (ctemp == gtemp)
	{
		file << "rx = [ rx , " << rMissile.GetX() << " ];\n";
		file << "ry = [ ry , " << rMissile.GetY() << " ];\n";
		file << "rz = [ rz , " << rMissile.GetZ() << " ];\n";
		ctemp = 1;
	}
	else
	{
		ctemp += 1;
	}

	height = rMissile.Magnitude() - Earth::GetR();
	if (height > maxheight)
	{
		maxheight = height;
	}

	if (Earth::GetR() <= rMissile.Magnitude())
	{
	}
	else
	{
		printf("-Trajectory End\n");
		printf("*******************************************\n");
		printf("Trajectory info:\n");
		command = 1;
		rEnd = rMissile;
		rDif = rEnd - rStart;
		ang = 2.0*asin(rDif.Magnitude()/(2.0*Earth::GetR()));
		distanceHit = ang*Earth::GetR();
		printf("-Trajectory Start Point:\n");
		rStart.Out();
		printf("-Trajectory End Point:\n");
		rEnd.Out();
		printf("-Distance Hit: %f1000 \n", distanceHit);
		printf("-Max Height: %f1000 \n", maxheight);
		printf("-Max Mach: %f1000 \n", maxmach);
	}
}

long double Missile::f(long double vx1)
{
	return vx1;
}
					 
long double Missile::g(long double vy1)
{
	return vy1;
}
					 
long double Missile::i(long double x1, long double vy1)
{
	return (2.0 * Earth::GetW()*vy1 + Earth::GetW() * Earth::GetW()*x1 + aMissile.GetX());
}
					
long double Missile::j(long double y1, long double vx1)
{
	return (-2.0 * Earth::GetW()*vx1 + Earth::GetW() * Earth::GetW()*y1 + aMissile.GetY()); 
}
					
long double Missile::d(long double vz1)
{
	return vz1;
}
					
long double Missile::e()
{
	return aMissile.GetZ();
}

int Missile::Command()
{
	return command;
}
