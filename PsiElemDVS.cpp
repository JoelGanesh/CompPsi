#include "CompPsi.h"
#include "Elementary.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CompPsi
{
	const int BRUTEFORCE_C = 4;
	float_dec_100 PsiElem::Psi0(uint64_t N)
	{
		uint64_t M = std::sqrt(N);
		uint64_t M0 = std::pow(N, 0.4) * std::pow(std::log(std::log(N)), 0.6) / std::pow(std::log(N), 0.6);
		return DependentVar(N, M, M0) + IndependentVar(N, M0);
	}

	float_dec_100 PsiElem::DependentVar(uint64_t N, uint64_t M, uint64_t M0)
	{
		float_dec_100 sum(0);
		uint64_t K1 = std::sqrt(M); // Batch length for segmented sieves over [M0, M].
		uint64_t K2 = std::sqrt(N / M0); // Batch length for segmented sieves over [N/M, N/M0].

		// We compute the sums
		// 1. sum_{M0 < d <= M} mu(d) sum_{n <= N/d} D_{Lambda}(n,N/n),
		// 2. sum_{M0 < d <= M} mu(d) sum_{m > d} Lambda(m) floor(N / m^2),
		// 3. sum_{M0 < m <= M} Lambda(m) sum_{n <= N/m} D_{mu}(n,N/n) and
		// 4. sum_{M0 < m <= M} Lambda(m) sum_{d >= m} mu(d) floor(N / d^2)
		// all at once, which are directly put in the total sum.

		// We process 'd' and 'm' in reverse order, 
		// so that the partial sums of inner sums do not have to be stored individually.

		uint64_t t = N / M + 1;
		float_dec_100 sum_DLambda = boost::multiprecision::lgamma(float_dec_100(t));
		float_dec_100 sum_LambdaFl(0);
		uint64_t sum_Dmu(1), sum_muFl(0);
		uint64_t DLambda_index(N / M), Dmu_index(N / M);
		for (uint64_t n = M; n > M0; n -= K1)
		{
			K1 = std::min(K1, n - M0);

			// We process the two sums in separate clauses to minimize memory usage.
			// (Local variables are garbage-collected when they become out-of-scope)
			{
				std::vector<int> mu = Elementary::sieve.MuSegmented(n - K1 + 1, K1);
				std::vector<Log> Lambda = Elementary::sieve.LambdaSegmented(n - K1 + 1, K1);
				for (uint64_t i = 0; i < K1; i++)
				{
					// Process (n - K1, n] by considering d = n - i, i = 0, 1, ..., K1 - 1.
					uint64_t d = n - i;

					sum_LambdaFl += Lambda[K1 - i - 1].numerical() * std::floor(N / (d * d));
					if (mu[K1 - i - 1] != 0)
					{
						sum_DLambda += SegmentSumDivSumLambda(DLambda_index, N / d, N, K2);
						DLambda_index = N / d;
						sum += mu[K1 - i - 1] * (sum_DLambda - sum_LambdaFl);
					}
				}
			}
			{
				std::vector<Log> Lambda = Elementary::sieve.LambdaSegmented(n - K1 + 1, K1);
				std::vector<int> mu = Elementary::sieve.MuSegmented(n - K1 + 1, K1);
				for (uint64_t i = 0; i < K1; i++)
				{
					// Process (n - K1, n] by considering m = n - i, i = 0, 1, ..., K1 - 1.
					uint64_t m = n - i;

					sum_muFl += mu[K1 - i - 1] * std::floor(N / (m * m));
					if (Lambda[K1 - i - 1].n != 0)
					{
						sum_Dmu += SegmentSumDivSumMu(Dmu_index, N / m, N, K2);
						Dmu_index = N / m;

						sum += Lambda[K1 - i - 1].numerical() * (sum_Dmu - sum_muFl);
					}
				}
			}
		}
		return sum;
	}

	uint64_t PsiElem::RestrictedDivSumMu(Factorization& F, uint64_t a)
	{
		int n = 1;
		for (PrimeFactor q : F.primeFactors())
		{
			n *= q.prime;
		}
		return RestrictedDSM_(F, F.primeFactors().size() - 1, 1, 1, a, n);
	}

	// The primefactors in F are sorted from small to large, so we can proceed by keeping track of the index.
	uint64_t PsiElem::RestrictedDSM_(Factorization& F, uint64_t index, uint64_t m1, uint64_t m2, uint64_t a, uint64_t n)
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
	
	float_dec_100 PsiElem::RestrictedDivSumLambda(Factorization& F, uint64_t a)
	{
		float_dec_100 sum(0);
		for (PrimeFactor pf : F.primeFactors())
		{
			Prime p = pf.prime;
			Exponent e = pf.exponent;
			
			sum += std::min(e, (int)std::floor(log(a) / log(p))) * boost::multiprecision::log(float_dec_100(p));
		}
		return sum;
	}
}