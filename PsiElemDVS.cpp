#include "CompPsi.h"
#include "Elementary.h"

namespace CompPsi
{
	// Computes a restricted divisor sum of the form sum_{d|n: d <= a} Lambda(d).
	static float_dec_T RestrictedDivSumLambda(Factorization& F, int64_t a);

	// Computes a restricted divisor sum of the form sum_{d|n, d <= a} mu(d).
	// Algorithm by Helfgott & Thompson, 2023.
	static int64_t RestrictedDivSumMu(Factorization& F, int64_t a);

	// Helper function for RestrictedDivSumMu; algorithm by Helfgott & Thompson, 2023.
	static int64_t RestrictedDSM_(Factorization& F, int64_t index, int64_t m1, int64_t m2, int64_t a, int64_t n);

	int64_t RestrictedDivSumMu(Factorization& F, int64_t a)
	{
		int n = 1;
		for (PrimeFactor q : F.primeFactors())
		{
			n *= q.prime;
		}
		return RestrictedDSM_(F, F.primeFactors().size() - 1, 1, 1, a, n);
	}

	// The primefactors in F are sorted from small to large, so we can proceed by keeping track of the index.
	int64_t RestrictedDSM_(Factorization& F, int64_t index, int64_t m1, int64_t m2, int64_t a, int64_t n)
	{
		if (m1 > a)
		{
			return 0;
		}
		else if (index < 0)
		{
			return 1;
		}
		else if (m2 * a >= n)
		{
			return 0;
		}
		
		Prime p = F.primeFactors()[index].prime;
		return RestrictedDSM_(F, index - 1, m1, p * m2, a, n) - RestrictedDSM_(F, index - 1, m1 * p, m2, a, n);
	}
	
	float_dec_T RestrictedDivSumLambda(Factorization& F, int64_t a)
	{
		float_dec_T sum(0);
		for (PrimeFactor pf : F.primeFactors())
		{
			Prime p = pf.prime;
			Exponent e = pf.exponent;
		
			sum += std::min(e, Elementary::Functions::log(a, p)) * boost::multiprecision::log(float_dec_T(p));
		}
		return sum;
	}

	// Returns sum_{n0 < n <= n1} g(F[n], N/n) in batches of length K.
	// Here, F[n] represents the prime factorization of n,
	template <typename T>
	static T SegmentSumDivSumF(int64_t n0, int64_t n1, int64_t N, int64_t K,
							   std::function<T(Factorization&, int64_t)> g)
	{
		T sum(0);
		for (int64_t n = n0 + 1; n <= n1; n += K)
		{
			// We consider the segment [n, n + K);
			// make sure that n + K - 1 <= n1.
			K = std::min(K, n1 - n + 1);

			// Compute the factorizations of [n, n+K)
			// and use them to compute the terms of the sum.
			std::vector<Factorization> F = Elementary::sieve.FactorizationSegmented(n, K);
			for (int64_t m = n; m < n + K; m++)
			{
				sum += g(F[m - n], N / m);
			}
		}
		return sum;
	}

	// Returns sum_{n0 < n <= n1} D_{Lambda}(n, a), computed in batches of length K.
	// Here D_{Lambda}(n,a) = sum_{m|n, m <= a} Lambda(m).
	static float_dec_T SegmentSumDivSumLambda(int64_t n0, int64_t n1, int64_t N, int64_t K)
	{
		return SegmentSumDivSumF<float_dec_T>(n0, n1, N, K, RestrictedDivSumLambda);
	}

	// Returns sum_{n0 < n <= n1} D_{mu}(n, N/n), computed in batches of length K.
	// Here D_{mu}(n,a) = sum_{d|n, d <= a} mu(d).
	static int64_t SegmentSumDivSumMu(int64_t n0, int64_t n1, int64_t N, int64_t K)
	{
		return SegmentSumDivSumF<int64_t>(n0, n1, N, K, RestrictedDivSumMu);
	}


	float_dec_T PsiElem::DependentVar(int64_t N, int64_t M, int64_t M0)
	{
		float_dec_T sum(0);
		int64_t K1 = std::sqrt(M); // Batch length for segmented sieves over [M0, M].
		int64_t K2 = std::sqrt(N / M0); // Batch length for segmented sieves over [N/M, N/M0].

		// We compute the sums
		// 1. sum_{M0 < d <= M} mu(d) sum_{n <= N/d} D_{Lambda}(n,N/n),
		// 2. sum_{M0 < d <= M} mu(d) sum_{m > d} Lambda(m) floor(N / m^2),
		// 3. sum_{M0 < m <= M} Lambda(m) sum_{n <= N/m} D_{mu}(n,N/n) and
		// 4. sum_{M0 < m <= M} Lambda(m) sum_{d >= m} mu(d) floor(N / d^2)
		// all at once, which are directly put in the total sum.

		// We process 'd' and 'm' in reverse order, 
		// so that the partial sums of inner sums do not have to be stored individually.

		float_dec_T sum_DLambda = boost::multiprecision::lgamma(float_dec_T(M + 1));
		float_dec_T sum_LambdaFl(0);
		int64_t sum_Dmu(1), sum_muFl(0);
		int64_t DLambda_index(M), Dmu_index(M);
		for (int64_t n = M; n > M0; n -= K1)
		{
			K1 = std::min(K1, n - M0);

			std::vector<int> mu = Elementary::sieve.MuSegmented(n - K1 + 1, K1);
			std::vector<Log> Lambda = Elementary::sieve.LambdaSegmented(n - K1 + 1, K1);
			for (int64_t i = 0; i < K1; i++)
			{
				// Process (n - K1, n] by considering d = n - i, i = 0, 1, ..., K1 - 1.
				int64_t d = n - i;

				if (d <= std::sqrt(N))
				{
					sum_muFl += mu[K1 - i - 1] * (N / (d * d));
				}
				if (Lambda[K1 - i - 1].n > 1)
				{
					sum_Dmu += SegmentSumDivSumMu(Dmu_index, N / d, N, K2);
					Dmu_index = N / d;

					sum += Lambda[K1 - i - 1].numerical() * (sum_Dmu - sum_muFl);
				}
				if (mu[K1 - i - 1] != 0)
				{
					sum_DLambda += SegmentSumDivSumLambda(DLambda_index, N / d, N, K2);
					DLambda_index = N / d;
					sum += mu[K1 - i - 1] * (sum_DLambda - sum_LambdaFl);
				}
				if (d <= std::sqrt(N))
				{
					sum_LambdaFl += Lambda[K1 - i - 1].numerical() * (N / (d * d));
				}
			}
		}
		return sum;
	}
}