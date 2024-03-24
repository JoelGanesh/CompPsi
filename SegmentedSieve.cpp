#include "Elementary.h"
#include "Utility.h"
#include "Types.h"

namespace Elementary
{
	std::vector<Prime> SegmentedSieve::Primes(int64_t N)
	{
		int64_t M = this->N;

		// If N is less than M, we have already stored all the necessary primes.
		if (M < N)
		{
			std::vector<bool> prime = std::vector<bool>(N - M, true); // Represents integers M+1, ..., N.

			// First remove the multiples of primes we already know of.
			for (Prime p : primes)
			{
				for (int64_t n = FirstMultipleAfter(M + 1, p); n <= N; n += p)
				{
					prime[n - (M + 1)] = false;
				}
			}

			// We proceed by using the idea behind the sieve of Eratosthenes to find the remaining primes.
			int64_t K = std::sqrt(N);
			for (int64_t n = M + 1; n <= K; n++)
			{
				// If n is prime, mark proper multiples of n to be composite.
				// Note that an integer 'n' in [M + 1, N] corresponds to the [n - (M+1)]-th entry in 'prime'.
				if (prime[n - (M + 1)])
				{
					// Small optimization: it suffices to mark multiples of n >= n^2.
					for (int64_t k = n * n; k <= N; k += n)
					{
						prime[k - (M + 1)] = false;
					}
				}
			}

			// The remaining indices correspond to the primes.
			std::vector<Prime> newPrimes = Utility::Generic::IndexAll(prime, true, M + 1);

			// We append the new primes to the existing list of primes.
			primes.insert(primes.end(), newPrimes.begin(), newPrimes.end());
			this->N = N;
		}
		return primes;
	}

	std::vector<Prime> SegmentedSieve::PrimesSegmented(int64_t N, int64_t S)
	{
		const int64_t M(std::sqrt(N + S));
		Primes(M);

		std::vector<bool> prime = std::vector<bool>(S, true); // Represents integers N, ..., N+S-1.

		// For any prime in the given list, mark multiples in the interval to be composite.
		// Note that an integer 'n' in the interval corresponds to the (n - N)-th entry in 'prime'.
		for (Prime p : primes)
		{
			if (p > M)
			{
				break;
			}
			for (int64_t n = FirstMultipleAfter(N, p); n < N + S; n += p)
			{
				prime[n - N] = false;
			}
		}

		// The remaining indices correspond to primes.
		std::vector<int64_t> primeIndices = Utility::Generic::IndexAll(prime, true, N);
		return primeIndices;
	}

	std::vector<int> SegmentedSieve::MuSegmented(int64_t N, int64_t S)
	{
		const int64_t M(std::sqrt(N + S));
		Primes(M);

		std::vector<int> mu = std::vector<int>(S, 1); // Represents mu values of N, ..., N+S-1.

		// We keep track of the products of prime factors in 'primes' of integers in [N, N+S),
		// in order to find a potential last prime factor of integers [N, N+S) larger than any prime in 'primes'.
		std::vector<int64_t> product = std::vector<int64_t>(S, 1);

		for (Prime p : primes)
		{
			if (p > M)
			{
				break;
			}
			// Multiples of p^2 should be assigned a value of 0, as they are not squarefree.
			int64_t q = p * p;
			for (int64_t n = FirstMultipleAfter(N, q); n < N + S; n += q)
			{
				mu[n - N] = 0;
			}
			// For the other multiples of p, we flip the sign accordingly.
			// (Note that mu is multiplicative, and mu(p) = -1 for any prime p.)
			for (int64_t n = FirstMultipleAfter(N, p); n < N + S; n += p)
			{
				mu[n - N] *= -1;
				if (mu[n - N] != 0)
				{
					product[n - N] *= p;
				}
			}
		}

		// We take into account that we may have missed the largest prime factor of integers in [N, N+S).
		// This is the case precisely when the stored product is not equal to the corresponding integer.
		for (int64_t n = N; n < N + S; n++)
		{
			if (product[n - N] != n)
			{
				mu[n - N] *= -1;
			}
		}

		return mu;
	}

	std::vector<Log> SegmentedSieve::LambdaSegmented(int64_t N, int64_t S)
	{
		const int32_t M(std::sqrt(N + S));
		Primes(M);

		std::vector<Log> Lambda = std::vector<Log>(S, Log(0)); // Represent Lambda values of N, ..., N+S-1.
		std::vector<bool> large_prime = std::vector<bool>(S, true); // To check for primes > S in [N, N+S).
		for (Prime p : primes)
		{
			if (p > M)
			{
				break;
			}
			// Mark multiples of p to not be large primes.
			for (int64_t n = FirstMultipleAfter(N, p); n < N + S; n += p)
			{
				large_prime[n - N] = false;
			}

			// Powers of p should be assigned a value of Log(p).
			int128_t q = p;
			while (q < N)
			{
				q *= p;
			}
			for (int128_t n = q; n < N + S; n *= p)
			{
				// Assuming n < N+S, n-N can be stored in an int64_t.
				Lambda[(int64_t)(n - N)] = Log(p);
			}
		}

		// We take into account that we missed the primes > sqrt(N+S).
		for (int64_t n = N; n < N + S; n++)
		{
			if (large_prime[n - N])
			{
				Lambda[n - N] = Log(n);
			}
		}

		return Lambda;
	}

	std::vector<Factorization> SegmentedSieve::FactorizationSegmented(int64_t N, int64_t S)
	{
		const int32_t M(std::sqrt(N + S));
		Primes(M);

		std::vector<Factorization> factorizations(S, Factorization()); // Represent factorizations of N, ..., N+S-1.
		std::vector<int64_t> products(S, 1);
		for (Prime p : primes)
		{
			if (p > M)
			{
				break;
			}
			int128_t q = p;
			Exponent j = 1;
			while (q < N + S)
			{
				for (int64_t k = FirstMultipleAfter(N, (int64_t)q); k < N + S; k += (int64_t)q)
				{
					// If p^{j+1} does not divide k, then v_p(k) = j, so we should add the pair (p, j).
					if (k % (p * q) != 0)
					{
						factorizations[k - N].AddFactor(p, j);
						products[k - N] *= (int64_t)q;
					}
				}
				q *= p;
				j += 1;
			}
		}

		// It remains for us to mark prime factors larger than S.
		for (int64_t k = N; k < N + S; k++)
		{
			if (k != products[k - N])
			{
				int64_t p = k / products[k - N];
				factorizations[k - N].AddFactor(p, 1);
			}
		}

		return factorizations;
	}

	int64_t SegmentedSieve::FirstMultipleAfter(int64_t a, int64_t k)
	{
		if (a % k == 0)
		{
			return a;
		}
		else
		{
			return a + k - (a % k);
		}
	}
}