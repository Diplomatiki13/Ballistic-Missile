#pragma once
#pragma warning(disable:4996)
#include "Earth.h"
#include "Vector.h"
#include "missile.h"

Missile missile1;
long double dtr(long double d);
int command;
long double totalTime = 0.0;
long double pi = 3.141592;
/********* Control Parameters *********/
long double launchPointTheta = 90.0;
long double launchPointPhi = 30.0;
long double launchTheta = 20.0;
long double launchPhi = 80.86;
long double dt = 0.01;
/**************************************/