#include "CompPsi.h"
#include "Fourier.h"
#include "Utility.h"
#include "Types.h"
#include "Elementary.h"

#include <boost/multiprecision/cpp_dec_float.hpp>

#include <fftw3.h>

namespace CompPsi
{
    float_dec_100 PsiFFT::Psi0(int64_t N)
    {
        Elementary::SegmentationArray<float_dec_100> segmentationArray(N);
        const int S = N * (pow(2, 3 * segmentationArray.getDelta()) - 1);
        float_dec_100 result = SegmentationFFT(N, segmentationArray) - SegmentationError(N, S, segmentationArray); // mu[z] defined -> mu.size() = z+1.
        //std::cout << "Psi0: " << result << std::endl;
        return result;
    }

    float_dec_100 PsiFFT::SegmentationFFT(const int64_t N, Elementary::SegmentationArray<float_dec_100> segmArray) // Returns sum_{k1+...+k4 <= k} Mu_z[k1]Mu_z[k2]One[k3]Log[k4] using FFT algorithm. F[j] is the j-th index of the segmentation of f and k = floor(log2(n)/delta).
    {
        std::vector<float_dec_100> segmentedMu = segmArray.Mu_M();
        std::vector<float_dec_100> segmentedOne = segmArray.One();
        std::vector<float_dec_100> segmentedLambda = segmArray.Lambda_M();
        std::vector<std::vector<float_dec_100>> segments{ segmentedMu, segmentedOne, segmentedLambda };
        //for (int i = 0; i < 3; i++)
        //{
        //    Utility::IO::Print(segments[i]);
        //}

        int64_t size = segmentedMu.size(); // Note: sizes of the segments are equal.
        Fourier::FFTSimple FFTLib(size, 3); // The parameter 3 represents the number of arrays to be convoluted.
        std::vector<float_dec_100> segmentation = FFTLib.Convolve(segments);
        //Utility::IO::Print(segmentation);
        //cout << "S: ";
        //Print(segmentation);
        float_dec_100 sum = 0;
        for (size_t k = 0; k < size; k++)
        {
            sum += segmentation[k];
        }

        //std::cout << "PsiFFTSegmentation: " << sum << std::endl;
        return sum;
    }

	float_dec_100 PsiFFT::SegmentationError(const int64_t N, const int64_t S, Elementary::SegmentationArray<float_dec_100> segmArray) // Calculates error term sum_{d1d2d3d4 > N, k(d1)+...+k(d4) <= k} mu_z(d1)mu_z(d2)log(d_4), where k(m) = floor(log2(m)/delta), k = k(N).
	{
		// We determine the prime factorizations of all integers in the critical interval (N, N+S]
		std::vector<Factorization> factorizations = Elementary::sieve.FactorizationSegmented(N + 1, S);

		int64_t M = std::sqrt(N);

		// We prepare tables for mu and Lambda for computation of errors.
		std::vector<int> segmentMu = Elementary::sieve.MuSegmented(1, M + 1);
		std::vector<Log> segmentLambda = Elementary::sieve.LambdaSegmented(1, M + 1);

        float_dec_100 error = 0;
        for (int64_t n = N + 1; n <= N + S; n++)
        {
            Factorization factorization = factorizations[n - (N + 1)];
            //std::cout << std::string(factorization) << std::endl;

			Tuple v{ 1, INT_MAX, INT_MAX }; // First entry corresponds with mu. We take ones for mu_z since mu_z(n) = 0 whenever p^2 divides n for some prime p.
			std::vector<Tuples> subfactorizations = PsiFFT::SubFactorizations(v, factorization);
			//int q = 0;
			//std::function<int(std::vector<std::vector<int>>)> f = [q](std::vector<std::vector<int>> x) mutable { for (std::vector<int> y : x) { std::cout << "S[" << q++ << "]: "; Utility::IO::Print(y); }; return 0; };
			//Utility::Generic::Map(subfactorizations, f);
			// As it is, we have all subfactorizations, i.e. a list { { a_{11}, ..., a_{1n_1} }, ..., { {a_{k1}, ..., a_{kn_k} } }, where each a_{ij} represents
			// a unique factorization of the primepower p_i^{e_i} dividing n into a fixed number of factors (determined by the size of the vector v above).
			// We now want to combine these subfactorizations to obtain all decompositions of n into that number of factors, for which we should consider the Cartesian product of these sets.
			// To circumvent this a little, we consider a bijection between each element of this Cartesian product with a unique integer in {1, ..., P}, where P = n_1 * ... * n_k.
			int P = 1;
			std::vector<int> sizes;
			for (Tuples subfactorization : subfactorizations)
			{
				int size = subfactorization.size();
				P *= size;
				sizes.push_back(size);
			}
			int64_t indexMax = segmArray.index(N);
			for (int p = 0; p < P; p++)
			{
				std::vector<int> m = Utility::Indexation::IntToTuple(p, sizes);
				//std::cout << "m: ";
				//Utility::IO::Print(m);
				// We are now considering the combination of tuples a_{1,m[1]}, ..., a_{k,m[k]}.
				// Recall that these tuples a_{j,m[j]} contain |v| entries, corresponding to the |v| parts of a unique decomposition of the primepower p_j^{e_j}.
				// We essentially want to finally combine the different prime powers into the |v| parts of the decomposition of n.
				std::vector<Factorization> decomposition;
				int64_t indexSum = 0; // A decomposition should only contribute if the sum of the segmentation indices does not exceed the segmentation index of N.
				for (int i = 0; i < v.size(); i++)
				{
					Factorization decomposition_factor = Factorization();
					std::vector<PrimeFactor> primeFactors = factorization.primeFactors();
					int64_t product = 1;
					for (int j = 0; j < primeFactors.size(); j++)
					{
						Prime primeFactor = primeFactors[j].prime;
						Tuple tuple = subfactorizations[j][m[j]];
						if (tuple[i] >= 1)
						{
							decomposition_factor.AddFactor(primeFactor, tuple[i], true);
						}
					}
					decomposition.push_back(decomposition_factor);

					indexSum += segmArray.index(decomposition_factor.n());
					if (indexSum > indexMax) // The sum of the indices is only going to increase, so the current decomposition is not going to contribute to the error term.
					{
						break;
					}
				}
				// If sum k(d_i) <= k(N), then add the term.
				if (indexSum <= indexMax)
				{
					//std::function<int(Factorization)> g = [](Factorization h) { return h.n(); };
					//std::cout << "D: "; 
					//Utility::IO::Print(Utility::Generic::Map(decomposition, g));
					// If the decomposition looks like (d0,d1,d2), we add mu_M(d0) * Lambda_M(d1)
					int64_t d = decomposition[0].n();
					int64_t m = decomposition[1].n();

					if (d <= M && m <= M)
					{
						float_dec_100 contribution = segmentMu[d - 1] * segmentLambda[m - 1].numerical();
						//std::cout << "Contribution: " << contribution << std::endl;
						error += contribution;
					}
				}
			}
		}

        //std::cout << "PsiFFTError: " << error << std::endl;
        return error;
    }
    
    std::vector<Tuples> PsiFFT::SubFactorizations(Tuple v, Factorization factorization) // List of non-negative integer tuples (a_i)_i for each maximal primepower p^k dividing n such that sum_i a_i = k, while a_i <= v_i
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