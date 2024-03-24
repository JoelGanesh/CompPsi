#include "Elementary.h"

namespace Elementary
{
	int Functions::FirstCongruenceAfter(int n, int a, int q)
	{
		int k = (n + a) - (n % q);
		if (a < n % q)
		{
			k += q;
		}
		return k;
	}

	Exponent Functions::log(int64_t a, int64_t p)
	{
		double f = std::log(a) / std::log(p);
		Exponent f_floor(f);
		if (f - f_floor > 1 - DOUBLE_ERROR_THRESHOLD)
		{
			// It is likely that f and f_floor are wrong
			// due to floating-point arithmetic.
			if (pow(p, f_floor + 1) <= a)
			{
				return f_floor + 1;
			}
		}
		return f_floor;
	}

	int128_t Functions::pow(int64_t p, int k)
	{
		if (k == 0)
		{
			return 1;
		}

		int128_t q = pow(p, k / 2);
		q *= q;
		if (k % 2 == 1)
		{
			q *= p;
		}
		return q;
	}

	int64_t Functions::floor(double alpha)
	{
		int64_t f = std::floor(alpha);
		if (alpha - f > 1 - DOUBLE_ERROR_THRESHOLD)
		{
			f++;
		}
		return f;
	}

	double Functions::round(double alpha)
	{
		int s = 1;
		if (alpha < 0)
		{
			alpha = -alpha;
			s = -1;
		}
		double alpha_int;
		double alpha_frac = modf(alpha, &alpha_int);
		if (alpha_frac < DOUBLE_ERROR_THRESHOLD)
		{
			return s * alpha_int;
		}
		if (alpha_frac > 1 - DOUBLE_ERROR_THRESHOLD)
		{
			return s * (alpha_int + 1);
		}
		return s * alpha;
	}

}