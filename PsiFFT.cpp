// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "CompPsi.h"
#include "Fourier.h"
#include "Utility.h"
#include "Elementary.h"

namespace CompPsi
{
	float_dec_T PsiFFT::Psi0(int64_t N)
	{
		Elementary::SegmentationArray<float_dec_T> segmentationArray(N);
		const int32_t S = int32_t(N * (pow(8, segmentationArray.getDelta()) - 1));
		return SegmentationFFT(N, segmentationArray) - SegmentationError(N, S, segmentationArray);
	}

	float_dec_T PsiFFT::SegmentationFFT(const int64_t N, Elementary::SegmentationArray<float_dec_T>& segmArray)
	{
		std::vector<float_dec_T> segmentedMu = segmArray.Mu_M();
		std::vector<float_dec_T> segmentedOne = segmArray.One();
		std::vector<float_dec_T> segmentedLambda = segmArray.Lambda_M();
		std::vector<std::vector<float_dec_T>> segments{ segmentedMu, segmentedOne, segmentedLambda };

		int64_t size = segmentedMu.size(); // Note: sizes of the segments are equal.
		Fourier::FFTSimple FFTLib(size, 3);
		std::vector<float_dec_T> segmentation = FFTLib.Convolve(segments);

		float_dec_T sum = 0;
		for (int64_t k = 0; k < size; k++)
		{
			sum += segmentation[k];
		}
		return sum;
	}

	float_dec_T PsiFFT::SegmentationError(const int64_t N, const int32_t S, Elementary::SegmentationArray<float_dec_T>& segmArray)
	{
		// We determine the prime factorizations of all integers in the critical interval (N, N+S]
		std::vector<Factorization> factorizations = Elementary::sieve.FactorizationSegmented(N + 1, S);
		int64_t M = std::sqrt(N);

		float_dec_T error = 0;
		for (int64_t n = N + 1; n <= N + S; n++)
		{
			Factorization factorization = factorizations[n - (N + 1)];

			Tuple v{ 1, INT_MAX, INT_MAX };

			std::vector<Tuples> subfactorizations;
			std::vector<PrimeFactor> primeFactors = factorization.primeFactors();
			for (PrimeFactor primeFactor : primeFactors)
			{
				int k = primeFactor.exponent;
				subfactorizations.push_back(Utility::Indexation::RestrictedTuples(k, v));
			}

			// 'subfactorizations' is a list { { a_{11}, ..., a_{1n_1} }, ..., { {a_{k1}, ..., a_{kn_k} } }, where each a_{ij} represents a distinct 
			// factorization of the primepower p_i^{e_i} dividing n into a fixed number of factors (determined by the size of the vector v above).
			// We now want to combine these subfactorizations to obtain all decompositions of n into that number of factors, i.e., we want to find the Cartesian product.
			// For this, we consider a bijection between each element of this Cartesian product with a unique integer in {1, ..., P}, where P = n_1 * ... * n_k.
			int P = 1;
			std::vector<int> sizes;
			for (Tuples subfactorization : subfactorizations)
			{
				int64_t size = subfactorization.size();
				P *= size;
				sizes.push_back(size);
			}
			int64_t indexMax = segmArray.index(N);
			for (int p = 0; p < P; p++)
			{
				std::vector<int> m = Utility::Indexation::IntToTuple(p, sizes);

				// We are now considering the combination of tuples a_{1,m[1]}, ..., a_{k,m[k]}.
				// Recall that these tuples a_{j,m[j]} contain |v| entries, corresponding to the |v| parts of a unique decomposition of the primepower p_j^{e_j}.
				// We essentially want to finally combine the different prime powers into the |v| parts of the decomposition of n.
				std::vector<Factorization> decomposition;

				// Decompositions only contribute if the sum of the segmentation indices does not exceed the segmentation index of N.
				int64_t indexSum = 0;
				bool earlyBreak = false;
				for (int i = 0; i < v.size(); i++)
				{
					Factorization decomposition_factor = Factorization();
					int64_t product = 1;
					for (int j = 0; j < primeFactors.size(); j++)
					{
						Prime primeFactor = primeFactors[j].prime;
						Tuple tuple = subfactorizations[j][m[j]];
						if (tuple[i] >= 1)
						{
							// The index i = 1 corresponds with Lambda. Note that if two distinct primes divide n, then Lambda(n) = 0.
							// This can be used to save some unnecessary computations.
							if (i == 1 && decomposition_factor.primeFactors().size() != 0)
							{
								earlyBreak = true;
								break;
							}
							decomposition_factor.AddFactor(primeFactor, tuple[i], true);
						}
					}
					if (!earlyBreak)
					{
						decomposition.push_back(decomposition_factor);

						indexSum += segmArray.index(decomposition_factor.n());
						if (indexSum > indexMax) // The sum of the indices is only going to increase, so the current decomposition is not going to contribute to the error term.
						{
							earlyBreak = true;
							break;
						}
					}
				}
				if (!earlyBreak)
				{
					// If the decomposition looks like (d0,d1,d2), we add mu_M(d0) * Lambda_M(d1)
					if (decomposition[0].n() <= M && decomposition[1].n() <= M && decomposition[1].n() > 1)
					{
						// By our previous precautions:
						// - mu(n) = (-1)^t where t is the number of distinct prime divisors of the square-free integer n;
						// - Lambda(n) = log(p), where p is the unique prime divisor of n.
						int mu = decomposition[0].primeFactors().size() % 2 == 0 ? 1 : -1;
						float_dec_T Lambda = Log(decomposition[1].primeFactors()[0].prime).numerical();
						float_dec_T contribution = mu * Lambda;
						error += contribution;
					}
				}
			}
		}
		return error;
	}

	// Computes all ways to split a Factorization into v.size() Factorizations ("subFactorizations"),
	// where the maximal exponent of prime factors in the i-th Factorization is bound by v[i].
	std::vector<Tuples> SubFactorizations(Tuple& v, Factorization& factorization)
	{
		std::vector<Tuples> subfactorizations;
		std::vector<PrimeFactor> primeFactors = factorization.primeFactors();
		for (PrimeFactor primeFactor : primeFactors)
		{
			int k = primeFactor.exponent;
			subfactorizations.push_back(Utility::Indexation::RestrictedTuples(k, v));
		}
		return subfactorizations;
	}
}