#ifndef ELEMENTARY_H
#define ELEMENTARY_H

#include <vector>

//#include "CompPsi.h"
#include "Types.h"
using namespace Types;

namespace Elementary
{
	// Class specialized for computation of segmented sieves.
	class SegmentedSieve
	{
	private:
		std::vector<Prime> primes;
		uint64_t N;

		// Returns first integer multiple of k occurring not less than a.
		static uint64_t FirstMultipleAfter(uint64_t a, uint64_t k);
	public:
		// Constructor
		SegmentedSieve() : N(1) {};

		// Computes and returns the list of the primes up to (and including) N.
		std::vector<Prime> Primes(uint64_t N);

		// Computes and returns the list of primes in [N, N + S).
		std::vector<Prime> PrimesSegmented(uint64_t N, uint64_t S);

		// Computes and returns the values of mu(n) for n in [N, N + S).
		std::vector<int> MuSegmented(uint64_t N, uint64_t S);

		// Computes and returns the values of Lambda(n) for n in [N, N + S).
		std::vector<Log> LambdaSegmented(uint64_t N, uint64_t S);

		// Computes and returns the complete primefactorization of all integers in [N, N + S).
		std::vector<Factorization> FactorizationSegmented(uint64_t N, uint64_t S);

		//std::vector<SqFreeFactorization> SqFreeFactorizationSegmented(uint64_t N, uint64_t S);
	};

	extern SegmentedSieve sieve;

	// Class specifically for finding Diophantine approximations.
	class DiophAppr
	{
	public:
		// Returns a tuple (p, p', q, s) of integers so that abs(alpha - p/q) <= 1/(qQ)
		// with gcd(p, q) = 1, q <= Q, and pp' = 1 mod q, while s = sgn(alpha - p/q).
		// Algorithm taken from a published paper from 2023 by H. A. Helfgott and L. Thompson.
		static std::tuple<int64_t, int64_t, int64_t, int> ApprByRedFrac(double alpha, int64_t Q);
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
}
#endif