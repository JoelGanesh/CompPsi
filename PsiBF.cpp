// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "CompPsi.h"
#include "Elementary.h"

namespace CompPsi
{
	float_dec_T PsiBF::Psi(int64_t N)
	{
		float_dec_T sum = 0;
		int64_t M = std::sqrt(N);
		std::vector<Prime> primes = Elementary::sieve.Primes(M);
		for (Prime p : primes)
		{
			if (p > M)
			{
				break;
			}
			float_dec_T log = boost::multiprecision::log(float_dec_T(p));
			sum += Elementary::Functions::log(N, p) * log;
		}

		for (int64_t n = M + 1; n <= N; n += M)
		{
			int64_t D = std::min(M, N + 1 - n);
			primes = Elementary::sieve.PrimesSegmented(n, D);

			for (Prime p : primes)
			{
				float_dec_T log = boost::multiprecision::log(float_dec_T(p));
				sum += Elementary::Functions::log(N, p) * log;
			}
		}
		return sum;
	}
}