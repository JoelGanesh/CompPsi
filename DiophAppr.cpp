// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "Elementary.h"

namespace Elementary
{
	std::tuple<int128_t, int128_t, int128_t, int> DiophAppr::ApprByRedFrac(Fraction alpha, int64_t Q)
	{
		int sgn_alpha = 1;
		if (alpha.IsNegative())
		{
			sgn_alpha = -1;
			alpha.Negate();
		}

		int128_t p[2]{ 0, 1 };
		int128_t q[2]{ 1, 0 };
		int s = 1;

		while (q[1] <= Q)
		{
			int128_t a = alpha.Floor();
			int128_t new_p = a * p[1] + p[0];
			int128_t new_q = a * q[1] + q[0];
			p[0] = p[1]; p[1] = new_p;
			q[0] = q[1]; q[1] = new_q;

			if (alpha.IsIntegral() && q[1] <= Q)
			{
				return std::make_tuple(sgn_alpha * p[1], (-sgn_alpha * s * q[0]) % q[1], q[1], 0);
			}

			alpha = alpha.FractionalPart();
			alpha.Invert();
			s *= -1;
		}
		// As q[1] > Q, we resort to the previous approximation, i.e. pm / qm.
		return std::make_tuple(sgn_alpha * p[0], (-sgn_alpha * s * q[1]) % q[0], q[0], sgn_alpha * s);
	}
}