#include "Types.h"

namespace Types
{
	Interval::Interval(int64_t start, int64_t end) : start(start), end(end)
	{ };

	Interval::Interval(int64_t a, int64_t b, int64_t c) : start(1), end(0)
	{
		int64_t D = b * b - 4 * a * c;
		if (D >= 0)
		{
			int64_t Q = std::sqrt(D);
			if (a < 0)
			{
				start = std::ceil((double)(-b + Q) / (2 * a));
				end = std::floor((double)(-b - Q) / (2 * a));
			}
			else if (a > 0)
			{
				// We have to distinguish the possibilities of D being a square or not;
				// for certain values of a and b this results in a different integer interval.
				if (Q * Q != D)
				{
					start = std::ceil((double)(-b - Q) / (2 * a));
					end = std::floor((double)(-b + Q) / (2 * a));
				}
				else
				{
					start = std::floor((double)(-b - Q) / (2 * a)) + 1;
					end = std::ceil((double)(-b + Q) / (2 * a)) - 1;
				}
			}
		}
	};

	void Interval::Shift(int64_t a)
	{
		if (start != LLONG_MIN && start != LLONG_MAX)
		{
			start += a;
		}
		if (end != LLONG_MIN && end != LLONG_MAX)
		{
			end += a;
		}
	}

	void Interval::Intersect(Interval& I)
	{
		start = std::max(start, I.start);
		end = std::min(end, I.end);
	};
}