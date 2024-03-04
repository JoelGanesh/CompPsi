#include "Elementary.h"

#define epsilon 1.0E-8

namespace Elementary
{
	std::tuple<int64_t, int64_t, int64_t, int> DiophAppr::ApprByRedFrac(Fraction alpha, int64_t Q)
	{
		int sgn_alpha = 1;
		if (alpha.IsNegative())
		{
			sgn_alpha = -1;
			alpha.Negate();
		}

		int64_t p[2]{ 0, 1 };
		int64_t q[2]{ 1, 0 };
		int s = 1;

		while (q[1] <= Q)
		{
			int64_t a = alpha.Floor();
			int64_t new_p = a * p[1] + p[0];
			int64_t new_q = a * q[1] + q[0];
			p[0] = p[1]; p[1] = new_p;
			q[0] = q[1]; q[1] = new_q;

			if (alpha.IsIntegral() && q[1] <= Q)
			{
				return std::make_tuple(sgn_alpha * p[1], (-sgn_alpha * s * q[0]) % q[1], q[1], 0);
			}
			alpha = alpha.FractionalPart(); alpha.Invert(); // alpha <- 1/(alpha - a).
			s *= -1;
		}
		// As q[1] > Q, we resort to the previous approximation, i.e. pm / qm.
		return std::make_tuple(sgn_alpha * p[0], (-sgn_alpha * s * q[1]) % q[0], q[0], sgn_alpha * s);
	}
}