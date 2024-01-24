#include "Elementary.h"

#define epsilon 1.0E-8

namespace Elementary
{
	std::tuple<uint64_t, uint64_t, uint64_t, int> DiophAppr::ApprByRedFrac(double alpha, uint64_t Q)
	{
		uint64_t p[2]{ 0, 1 };
		uint64_t q[2]{ 1, 0 };
		int s = 1;

		while (q[1] <= Q)
		{
			uint64_t a(alpha); // a = floor(alpha).
			uint64_t new_p = a * p[1] + p[0];
			uint64_t new_q = a * q[1] + q[0];
			if (alpha - a < epsilon)
			{
				return std::make_tuple(new_p, (s * q[0]) % new_q, new_q, 0);
			}

			p[0] = p[1]; p[1] = new_p;
			q[0] = q[1]; q[1] = new_q;

			alpha = 1 / (alpha - a);
			s *= -1;
		}
		// As q[1] > Q, we resort to the previous approximation, i.e. pm / qm.
		return std::make_tuple(p[0], (s * q[1]) % q[0], q[0], -s);
	}
}