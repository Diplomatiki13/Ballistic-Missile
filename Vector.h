#pragma once
#include <iostream>
using namespace std;

class Vector
{
public:
	Vector();
	Vector(long double phi_in, long double theta_in, long double magnitude_in, long double trigger);
	Vector(long double x_in, long double y_in, long double z_in);
	Vector operator *(Vector v2);
	Vector operator +(Vector v2);
	Vector operator -(Vector v2);
	Vector operator *(long double c);
	Vector operator /(long double c);
	void Out();
	void Angles();
	long double GetX();
	long double GetY();
	long double GetZ();
	long double Magnitude();
	long double Phi();
	long double Theta();
	Vector Direction();
private:
	long double pi = 3.141592;
	long double x;
	long double y;
	long double z;
	long double magnitude;
	long double phi;
	long double theta;
};
