// Lanitious Missile Mk1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include "Program.h"
#include <iostream>
#include <fstream>
#include <ctime>
using namespace std;

int main()
{
	totalTime = 0.0;
	command = 0;
	time_t taim = time(0);
	struct tm *now = localtime(&taim);
	char buffer[80];
	strftime(buffer, 80, "%Y-%m-%d-%H-%M-%S.m", now);
	if (launchPointTheta < 0.0 ||
		launchPointTheta > 360.0 ||
		launchPointPhi < -90.0 ||
		launchPointPhi > 90.0 ||
		launchTheta < 0.0 ||
		launchTheta > 360.0 ||
		launchPhi < -90.0 ||
		launchPhi > 90.0)
	{
		printf("False Initial Angles");
	}
	else
	{
		printf("*******************************************\n");
		printf("-Missile Launch\n");
		ofstream file;
		file.open(buffer);
		missile1.Initialise(dtr(launchPointTheta), dtr(launchPointPhi), dtr(launchTheta), dtr(launchPhi), dt, file);
		missile1.updatePosAndVel(dt, file);
		totalTime += dt;
		while (command == 0)
		{
			missile1.updateForces(dt);
			missile1.updatePosAndVel(dt, file);
			totalTime += dt;
			command = missile1.Command();
		}
		printf("-Total Time: %f1000 \n", totalTime);
		printf("*******************************************\n");
		missile1.CalcGraph(file);
		file << "plot3(ex,ey,ez,'b+-',rx,ry,rz,'r+-')";
		file.close();
		printf("Graph Finished \n");
		printf("*******************************************\n");
	}
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file

long double dtr(long double d)
{
	return d * (pi / 180.0);
}
