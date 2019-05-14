#pragma once
#include <iostream>
#include <cmath>

//Class that contains all needed earth data
class Earth
{
public:

	static long double GetR();

	static long double GetW();

	static long double GetGr(long double h);

	static long double GetDen(long double h, int hinh);

	static long double GetTemp(long double h);

	static long double GetVisc(long double h);

private:
	static constexpr long double gr0 = 9.824;
	static constexpr long double M = 0.0289644;
	static constexpr long double Rair = 8.3144598;
	static constexpr long double R = 6371000.0;
	static constexpr long double W = 0.000072921159;
};
