// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#ifndef ELEMENTARY_H
#define ELEMENTARY_H

#include "Types.h"

#include <vector>

using namespace Types;

namespace Elementary
{
	// Class specialized for computation of segmented sieves.
	class SegmentedSieve
	{
		private:
		std::vector<Prime> primes;
		int64_t N;

		// Returns first integer multiple of k occurring not less than a.
		static int64_t FirstMultipleAfter(int64_t a, int64_t k);
		public:
		SegmentedSieve() : N(1) {};

		// Computes and returns the list of the primes up to (and including) N.
		std::vector<Prime> Primes(int64_t N);

		// Computes and returns the list of primes in [N, N + S).
		std::vector<Prime> PrimesSegmented(int64_t N, int64_t S);

		// Computes and returns the values of mu(n) for n in [N, N + S).
		std::vector<int> MuSegmented(int64_t N, int64_t S);

		// Computes and returns the values of Lambda(n) for n in [N, N + S).
		std::vector<Log> LambdaSegmented(int64_t N, int64_t S);

		// Computes and returns the complete primefactorization of all integers in [N, N + S).
		std::vector<Factorization> FactorizationSegmented(int64_t N, int64_t S);
	};

	extern SegmentedSieve sieve;

	// Class specifically for finding Diophantine approximations.
	class DiophAppr
	{
		public:
		// Returns a tuple (p, p', q, s) of integers so that abs(alpha - p/q) <= 1/(qQ)
		// with gcd(p, q) = 1, q <= Q, and pp' = 1 mod q, while s = sgn(alpha - p/q).
		// Algorithm taken from a published paper from 2023 by H. A. Helfgott and L. Thompson.
		static std::tuple<int128_t, int128_t, int128_t, int> ApprByRedFrac(Fraction alpha, int64_t Q);
	};

	// Implementations of several elementary functions.
	class Functions
	{
		public:
		// Returns the first integer k >= n congruent to a modulo q.
		// We assume that 0 <= a < q.
		static int FirstCongruenceAfter(int n, int a, int q);

		// Returns the floor of log_p(a), assuming p >= 2.
		// Resolves certain issues with floating-point arithmetic.
		static Exponent log(int64_t a, int64_t p);

		// Returns p^k using an advanced technique to compute powers.
		static int128_t pow(int64_t p, int k);

		// Returns the floor of a real number.
		// If alpha - floor(alpha) is too large (> 1 - 1E-8),
		// assumes error due to floating-point arithmetic.
		static int64_t floor(double alpha);

		// Rounds double to integer when too close.
		// Assumes that an error occurred due to floating-point arithmetic.
		static double round(double alpha);
	};
}
#endif