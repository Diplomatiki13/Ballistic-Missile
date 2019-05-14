#include "pch.h"
#include "Earth.h"

long double Earth::GetR()
{
	return R;
}

long double Earth::GetW()
{
	return W;
}

long double Earth::GetGr(long double h)
{
	return gr0 * ((R*R) / ((R + h)*(R + h)));
}

long double Earth::GetDen(long double h, int hinh)
{
	if (h>=0 && h<11000)
	{
		return 1.2250*pow((288.15 / (288.15 - 0.0065*(h - 0))), (1 + (gr0*M) / (R*(-0.0065))));
	}
	else if (h>=11000 && h<20000)
	{
		return 0.36391*exp((-gr0*M*(h - 11000))/(R*216.65));
	}
	else if (h >= 20000 && h < 32000)
	{
		return 0.08803*pow((216.65 / (216.65 + 0.001*(h - 20000))), (1 + (gr0*M) / (R*(0.001))));
	}
	else if (h >= 32000 && h < 47000)
	{
		return 0.01322*pow((228.65 / (228.65 + 0.0028*(h - 32000))), (1 + (gr0*M) / (R*(0.0028))));
	}
	else if (h >= 47000 && h < 51000)
	{
		return 0.00143*exp((-gr0 * M*(h - 47000)) / (R*270.65));
	}
	else if (h >= 51000 && h < 71000)
	{
		return 0.00086*pow((270.65 / (270.65 - 0.0028*(h - 51000))), (1 + (gr0*M) / (R*(-0.0028))));
	}
	else if (h >= 71000 && h < 86000)
	{
		return 0.00064*pow((214.65 / (214.65 - 0.002*(h - 71000))), (1 + (gr0*M) / (R*(-0.002))));
	}
	else if(h >= 86000)
	{
		if (hinh == 0)
		{
			printf("Error: Height over 86000! \n");
		}
		return 0.00064*pow((214.65 / (214.65 - 0.002*(86000 - 71000))), (1 + (gr0*M) / (R*(-0.002))));
	}
	else
	{
		printf("Negative Height");
		return 1.2250;
	}
}

long double Earth::GetTemp(long double h)
{
	if (h >= 0 && h < 11000)
	{
		return 288.15 - 0.0065*(h - 0);
	}
	else if (h >= 11000 && h < 20000)
	{
		return 216.65;
	}
	else if (h >= 20000 && h < 32000)
	{
		return 216.65 + 0.001*(h - 20000);
	}
	else if (h >= 32000 && h < 47000)
	{
		return 228.65 + 0.0028*(h - 32000);
	}
	else if (h >= 47000 && h < 51000)
	{
		return 270.65;
	}
	else if (h >= 51000 && h < 71000)
	{
		return 270.65 - 0.0028*(h - 51000);
	}
	else if (h >= 71000 && h < 86000)
	{
		return 214.65 - 0.002*(h - 71000);
	}
	else if (h>= 86000)
	{
		return 214.65 - 0.002*(86000 - 71000);
	}
	else
	{
		printf("Negative Height");
		return 288.15;
	}
}

long double Earth::GetVisc(long double h)
{
	return (1.458*0.000001*pow(GetTemp(h), (3.0 / 2.0))) / (GetTemp(h) + 110.4);
}
