/*#include "Utility.h"

namespace Utility
{
	// Create interval_set I so that x >= 1 lies in I iff ax^2 + bx + c >= 0.
	interval_set Interval::FromQuadIneq(int a, int b, int c)
	{
		interval_set intervalset;
		if (a == 0)
		{
			return FromLinIneq(b, c);
		}
		else
		{
			int x0 = -b;
			int D = b * b - 4 * a * c;
			if (a > 0)
			{
				if (D <= 0)
				{
					intervalset.add(interval::right_open(1, INT_MAX));
				}
				else
				{
					// Solution set: (-inf, xm] U [xp, inf) with xm < xp the (real) roots of the quadratic.
					int sqrtD = (int)std::floor(std::sqrt(D));
					int xm = (-b - sqrtD) / (2 * a);
					if (a * xm * xm + b * xm + c < 0)
					{
						xm--;
					}
					if (xm >= 1)
					{
						intervalset.add(interval::right_open(1, xm + 1));
					}

					int xp = (-b + sqrtD) / (2 * a);
					if (a * xp * xp + b * xp + c < 0)
					{
						xp++;
					}
					if (xp >= 1)
					{
						intervalset.add(interval::right_open(xp, INT_MAX));
					}
				}
			}
			else if (D == 0)
			{
				// Solution set: { -b / 2a }.
				if (b % (2 * a) == 0)
				{
					int k = -b / (2 * a);
					if (k >= 1)
					{
						intervalset.add(interval::right_open(k, k + 1));
					}
				}
			}
			else if (D > 0)
			{
				// Solution set: [xm, xp] with xm < xp the (real) roots of the quadratic.
				int sqrtD = (int)std::floor(std::sqrt(D));
				int xm = (-b - sqrtD) / (2 * a) - 1;
				if (a * xm * xm + b * xm + c < 0)
				{
					xm++;
				}
				int xp = (-b + sqrtD) / (2 * a) + 1;
				if (a * xp * xp + b * xp + c < 0)
				{
					xp--;
				}

				if (xm <= xp && xp >= 1)
				{
					intervalset.add(interval::right_open(std::max(1, xm), xp + 1));
				}
			}
			return intervalset;
		}
	}

	// Create interval_set I such that x >= 1 in I iff ax + b >= 0.
	interval_set Interval::FromLinIneq(int a, int b)
	{
		interval_set intervalset;
		if (a == 0)
		{
			if (b >= 0)
			{
				intervalset.add(interval::right_open(1, INT_MAX));
				return intervalset;
			}
		}
		else
		{
			int x = (-b) / a;
			if (a > 0)
			{
				intervalset.add(interval::right_open(std::max(1, x + 1), INT_MAX));
			}
			else if (x >= 1)
			{
				intervalset.add(interval::right_open(1, x + 1));
			}
		}
		return intervalset;
	}

}*/