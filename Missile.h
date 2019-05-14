#pragma once
#include "Vector.h"
#include "Earth.h" 

class Missile
{
public:
	void CalcGraph(ofstream& file);
	void Initialise(long double lt, long double lp, long double t, long double p, long double dt, ofstream& file);
	void updateForces(long double dt);
	void updatePosAndVel(long double h, ofstream& file);
	long double f(long double vx1);
	long double g(long double vy1);
	long double i(long double x1, long double vy1);
	long double j(long double y1, long double vx1);
	long double d(long double vz1);
	long double e();
	int Command();
private:
	Vector rMissile;
	Vector vMissile;
	Vector aMissile;
	Vector rStart;
	Vector rEnd;
	Vector rDif;
	Vector B;
	Vector T;
	Vector D;
	Vector RF;
	int inhibitor;
	int hinhibitor;
	long double launchPhiInit;
	long double launchThetaInit;
	long double ltInit;
	long double lpInit;
	long double RefS;
	long double MaxS;
	long double Cd;
	long double Cdt;
	long double Cdw;
	long double Cda;
	long double fb;
	long double reynolds;
	long double mach;
	long double ang;
	long double height;
	long double maxheight = 0.0;
	long double maxmach = 0.0;
	long double distanceHit;
	long double fuelmass;
	long double pi = 3.141592;
	int command;
	int gtemp;
	int ctemp;
	/********* Missile Parameters *********/
	long double missileNozDiam = 0.741;
	long double missileR = 0.3705;
	long double missileDiameter = 0.741;
	long double missileLength = 9.61;
	long double missileLfore = 1.2;
	long double fuelmassfull = 5700.0;
	long double framemass = 1170.0;
	long double isp = 240.0;
	long double massPerS = 53.4;
	/**************************************/
};