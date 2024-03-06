#ifndef ELEMENTARY_H
#define ELEMENTARY_H

#include <vector>

//#include "CompPsi.h"
#include "Types.h"
using namespace Types;

#define DOUBLE_ERROR_THRESHOLD 1.0E-8

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
		// Constructor
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

		//std::vector<SqFreeFactorization> SqFreeFactorizationSegmented(int64_t N, int64_t S);
	};

	extern SegmentedSieve sieve;

	// Class specifically for finding Diophantine approximations.
	class DiophAppr
	{
	public:
		// Returns a tuple (p, p', q, s) of integers so that abs(alpha - p/q) <= 1/(qQ)
		// with gcd(p, q) = 1, q <= Q, and pp' = 1 mod q, while s = sgn(alpha - p/q).
		// Algorithm taken from a published paper from 2023 by H. A. Helfgott and L. Thompson.
		static std::tuple<int64_t, int64_t, int64_t, int> ApprByRedFrac(Fraction alpha, int64_t Q);
	};

	// Returns the first integer k >= n congruent to a modulo q.
	// We assume that 0 <= a < q.
	static int FirstCongruenceAfter(int n, int a, int q)
	{
		int k = (n + a) - (n % q);
		if (a < n % q)
		{
			k += q;
		}
		return k;
	}

	// Returns the floor of log_p(a), assuming p >= 2.
	// Resolves certain issues with floating-point arithmetic.
	static int log(int64_t a, int64_t p)
	{
		double f = std::log(a) / std::log(p);
		int f_floor(f);
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

	// Returns p^k using an advanced technique to compute powers.
	static int64_t pow(int64_t p, int k)
	{
		int64_t q = 1;
		while (k > 0)
		{
			if (k % 2 == 1)
			{
				q *= p;
				k--;
			}
			else
			{
				q *= q;
				k /= 2;
			}
		}
		return q;
	}
}
#endif