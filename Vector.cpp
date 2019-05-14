#include "pch.h"
#include "Vector.h"
#include <cmath>

Vector::Vector()
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
	magnitude = 0.0;
	phi = 0.0;
	theta = 0.0;
}

Vector::Vector(long double phi_in, long double theta_in, long double magnitude_in, long double trigger)
{
	phi = phi_in;
	theta = theta_in;
	magnitude = magnitude_in;
	x = magnitude * cos(phi) * cos(theta);
	y = magnitude * cos(phi) * sin(theta);
	z = magnitude * sin(phi);
}

Vector::Vector(long double x_in, long double y_in, long double z_in)
{
	x = x_in;
	y = y_in;
	z = z_in;
	magnitude = sqrt(x*x + y*y + z*z);
	Angles();
}

Vector Vector::operator*(Vector v2)
{
	return Vector(y*v2.z - z*v2.y, z*v2.x - x*v2.z, x*v2.y - y*v2.x);
}

Vector Vector::operator+(Vector v2)
{
	return Vector(x + v2.x, y + v2.y, z + v2.z);
}

Vector Vector::operator-(Vector v2)
{
	return Vector(x - v2.x, y - v2.y, z - v2.z);
}

Vector Vector::operator*(long double c)
{
	return Vector(x * c, y * c, z * c);
}

Vector Vector::operator/(long double c)
{
	return Vector(x / c, y / c, z /c);
}

void Vector::Out()
{
	printf("  -Longitude: %f1000 \n", theta*(180.0 / pi));
	printf("  -Latitude: %f1000 \n", phi*(180.0 / pi));
}

void Vector::Angles()
{
	phi = asin(z / magnitude);

	if (x == 0.0 && y == 0.0)
	{
		theta = 0.0;
	}
	else
	{
		if (x == 0.0)
		{
			if (y > 0.0)
			{
				theta = pi / 2.0;
			}
			else if (y < 0.0)
			{
				theta = (3.0*pi) / 2.0;
			}
			else
			{
				printf("Serious Angle Error\n");
			}
		}
		if (y == 0.0)
		{
			if (x > 0.0)
			{
				theta = 0.0;
			}
			else if (x < 0.0)
			{
				theta = pi;
			}
			else
			{
				printf("Serious Angle Error\n");
			}
		}
		if (x != 0.0 && y != 0.0)
		{
			if (x > 0.0 && y > 0.0)
			{
				theta = (pi / 2.0) - asin(x / sqrt(x*x + y * y));
				if (theta < 0.0 && theta > -0.00001)
				{
					theta = 0.0;
				}
				/*if (theta < 0.0 || theta >(pi / 2.0))
				{
					printf("Error: 0.0 < %f < 1.57\n",theta);
				}*/
			}
			else if (x < 0.0 && y > 0.0)
			{
				theta = (pi / 2.0) + asin((-x) / sqrt(x*x + y * y));
				/*if (theta < (pi / 2.0) || theta > pi)
				{
					printf("Error: 1.57 < %f < 3.14\n", theta);
				}*/
			}
			else if (x < 0.0 && y < 0.0)
			{
				theta = (3.0* pi / 2.0) - asin((-x) / sqrt(x*x + y * y));
				/*if (theta < pi || theta >(3.0* pi / 2.0))
				{
					printf("Error: 3.14 < %f < 4.71\n", theta);
				}*/
			}
			else if (x > 0.0 && y < 0.0)
			{
				theta = (3.0* pi / 2.0) + asin(x / sqrt(x*x + y * y));
				/*if (theta < (3.0* pi / 2.0) || theta >(2.0* pi))
				{
					printf("Error: 4.71 < %f < 6.28\n", theta);
				}*/
			}
		}
	}
}

long double Vector::GetX()
{
	return x;
}
long double Vector::GetY()
{
	return y;
}
long double Vector::GetZ()
{
	return z;
}

long double Vector::Magnitude()
{
	return magnitude;
}

long double Vector::Phi()
{
	return phi;
}

long double Vector::Theta()
{
	return theta;
}

Vector Vector::Direction()
{
	return Vector( x/magnitude , y/magnitude , z/magnitude );
}
