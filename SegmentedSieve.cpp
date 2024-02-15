#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stack>
#include <complex>

#include "Elementary.h"
#include "Utility.h"
#include "Types.h"

namespace Elementary
{
    std::vector<Prime> SegmentedSieve::Primes(uint64_t N)
    {
        uint64_t M = this->N;

        // If N is less than M, we have already stored all the necessary primes.
        // If N is not much smaller than M, it might be benificial to just return all primes.
        //if (N <= M / 2)
        //{
        //    // We return the subvector of primes consisting of all primes not larger than N.
        //    uint64_t max_index = Utility::Generic::FindIndex(primes, N);
        //    return std::vector<Prime>(primes.begin(), primes.begin() + max_index + 1);
        //}
        if (M < N)
        {
            std::vector<bool> prime = std::vector<bool>(N - M, true); // Represents integers M+1, ..., N.

            // First remove the multiples of primes we already know of.
            for (Prime p : primes)
            {
                for (uint64_t n = FirstMultipleAfter(M + 1, p); n <= N; n += p)
                {
                    prime[n - (M + 1)] = false;
                }
            }

            // We proceed by using the idea behind the sieve of Eratosthenes to find the remaining primes.
            uint64_t K = std::sqrt(N);
            for (uint64_t n = M + 1; n <= K; n++)
            {
                // If n is prime, mark proper multiples of n to be composite.
                // Note that an integer 'n' in [M + 1, N] corresponds to the [n - (M+1)]-th entry in 'prime'.
                if (prime[n - (M + 1)])
                {
                    // Small optimization: it suffices to mark multiples of n >= n^2.
                    for (uint64_t k = n * n; k <= N; k += n)
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

    std::vector<Prime> SegmentedSieve::PrimesSegmented(uint64_t N, uint64_t S)
    {
        const uint64_t M(std::sqrt(N + S));
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
            for (uint64_t n = FirstMultipleAfter(N, p); n < N + S; n += p)
            {
                prime[n - N] = false;
            }
        }

        // The remaining indices correspond to primes.
        std::vector<uint64_t> primeIndices = Utility::Generic::IndexAll(prime, true, N);
        return primeIndices;
    }

    std::vector<int> SegmentedSieve::MuSegmented(uint64_t N, uint64_t S)
    {
        const uint64_t M(std::sqrt(N + S));
        Primes(M);

        std::vector<int> mu = std::vector<int>(S, 1); // Represents mu values of N, ..., N+S-1.

        // We keep track of the products of prime factors in 'primes' of integers in [N, N+S),
        // in order to find a potential last prime factor of integers [N, N+S) larger than any prime in 'primes'.
        std::vector<uint64_t> product = std::vector<uint64_t>(S, 1);

        for (Prime p : primes)
        {
            if (p > M)
            {
                break;
            }
            // Multiples of p^2 should be assigned a value of 0, as they are not squarefree.
            uint64_t q = p * p;
            for (uint64_t n = FirstMultipleAfter(N, q); n < N + S; n += q)
            {
                mu[n - N] = 0;
            }
            // For the other multiples of p, we flip the sign accordingly.
            // (Note that mu is multiplicative, and mu(p) = -1 for any prime p.)
            for (uint64_t n = FirstMultipleAfter(N, p); n < N + S; n += p)
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
        for (uint64_t n = N; n < N + S; n++)
        {
            if (product[n - N] != n)
            {
                mu[n - N] *= -1;
            }
        }

        return mu;
    }

    std::vector<Log> SegmentedSieve::LambdaSegmented(uint64_t N, uint64_t S)
    {
        const uint32_t M(std::sqrt(N + S));
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
            for (uint64_t n = FirstMultipleAfter(N, p); n < N + S; n += p)
            {
                large_prime[n - N] = false;
            }

            // Powers of p should be assigned a value of Log(p).
            PrimePower q = p;
            while (q < N)
            {
                q *= p;
            }
            for (uint64_t n = q; n < N + S; n *= p)
            {
                Lambda[n - N] = Log(p);
            }
        }

        // We take into account that we missed the primes > S.
        for (uint64_t n = N; n < N + S; n++)
        {
            if (large_prime[n - N])
            {
                Lambda[n - N] = Log(n);
            }
        }

        return Lambda;
    }

    std::vector<Factorization> SegmentedSieve::FactorizationSegmented(uint64_t N, uint64_t S)
    {
        const uint32_t M(std::sqrt(N + S));
        Primes(M);

        std::vector<Factorization> factorizations(S, Factorization()); // Represent factorizations of N, ..., N+S-1.
        std::vector<uint64_t> products(S, 1);
        for (Prime p : primes)
        {
            if (p > M)
            {
                break;
            }
            PrimePower q = p;
            Exponent j = 1;
            while (q < N + S)
            {
                for (uint64_t k = FirstMultipleAfter(N, q); k < N + S; k += q)
                {
                    // If p^{j+1} does not divide k, then v_p(k) = j, so we should add the pair (p, j).
                    if (k % (p * q) != 0)
                    {
                        factorizations[k - N].AddFactor(p, j);
                        products[k - N] *= q;
                    }
                }
                q *= p;
                j += 1;
            }
        }

        // It remains for us to mark prime factors larger than S.
        for (uint64_t k = N; k < N + S; k++)
        {
            if (k != products[k - N])
            {
                uint64_t p = k / products[k - N];
                factorizations[k - N].AddFactor(p, 1);
            }
        }

        return factorizations;
    }

    uint64_t SegmentedSieve::FirstMultipleAfter(uint64_t a, uint64_t k)
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