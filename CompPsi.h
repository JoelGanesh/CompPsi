#ifndef COMPPSI_H
#define COMPPSI_H

#include <vector>
#include <functional>

#include "Utility.h"
#include "Elementary.h"
#include "SegmentationArray.h"
#include "Types.h"

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace Types;

namespace CompPsi
{
    // A brute force (or naive) approach to compute Psi(N).
    class PsiBF
    {
    public:
        // Computes Psi(N) naively, sufficient for "small" values of N.
        // Allows an additional parameter specifying an existing sieve
        // in order to reduce unnecessary repeated computations.
        static float_dec_100 Psi(uint64_t N);
    };

    // "Abstract" subclass of CompPsi, where an implementation
    // of Psi(N) is based on an underlying implementation of Psi0,
    // to be defined in subclasses of PsiPrep.
    // Note: an implementation of Psi0 *must* be present in a derived class.
    template <class T>
    class PsiPrep
    {
    public:
        // Computes and returns an approximation of Psi(N).
        static float_dec_100 Psi(uint64_t N)
        {
            float_dec_100 result = PsiPreparation(N);
            result -= T::Psi0(N); // An implementation of Psi0 should be present in T.
            return result;
        }

    private:
        // Pure virtual method implemented in subclasses.
        // Computes and returns an approximation of
        // Psi0(N) := sum_{mdk <= N : m, d <= sqrt(N)} Lambda(m)mu(d).
        virtual float_dec_100 Psi0(uint64_t N) = 0;

        // Transforms the computation of Psi(N) to computing
        // sum_{mdk <= N : m, d <= sqrt(N)} Lambda(m)mu(d).
        static float_dec_100 PsiPreparation(uint64_t N)
        {
            uint32_t M = std::sqrt(N);
            PsiBF psiBF;
            float_dec_100 result = psiBF.Psi(M);
            std::vector<int> mu = Elementary::sieve.MuSegmented(1, M);
            for (uint32_t d = 1; d <= M; d++)
            {
                uint64_t temp = N / d;
                float_dec_100 k(temp);
                result += mu[d - 1] * boost::multiprecision::lgamma(k + 1);
            }
            std::cout << "PsiPrep: " << result << std::endl;
            return result;
        }

    };

    // An elementary approach to compute Psi(N) efficiently.
    class PsiElem : public PsiPrep<PsiElem>
    {
    private:
        // Computes Psi0(N) using elementary methods.
        // Definition of Psi0 can be found in base class.
        static float_dec_100 Psi0(uint64_t N);

        // Computes sum_{mdk <= N, M0 < max(m,d) <= M} Lambda(m)mu(d).
        static float_dec_100 DependentVar(uint64_t N, uint64_t M, uint64_t M0);

        // Returns sum_{n0 < n <= n1} D_{Lambda}(n, a), computed in batches of length K.
        // Here D_{Lambda}(n,a) = sum_{m|n, m <= a} Lambda(m).
        static float_dec_100 SegmentSumDivSumLambda(uint64_t n0, uint64_t n1, uint64_t N, uint64_t K)
        {
            return SegmentSumDivSumF<float_dec_100>(n0, n1, N, K, RestrictedDivSumLambda);
        }

        // Returns sum_{n0 < n <= n1} D_{mu}(n, N/n), computed in batches of length K.
        // Here D_{mu}(n,a) = sum_{d|n, d <= a} mu(d).
        static uint64_t SegmentSumDivSumMu(uint64_t n0, uint64_t n1, uint64_t N, uint64_t K)
        {
            return SegmentSumDivSumF<uint64_t>(n0, n1, N, K, RestrictedDivSumMu);
        }

        // Returns sum_{n0 < n <= n1} g(F[n], N/n) in batches of length K.
        // Here, F[n] represents the prime factorization of n,
        template <typename T>
        static T SegmentSumDivSumF(uint64_t n0, uint64_t n1, uint64_t N, uint64_t K,
                                   std::function<T(Factorization&, uint64_t)> g)
        {
            T sum(0);
            for (uint64_t n = n0 + 1; n <= n1; n += K)
            {
                // We consider the segment [n, n + K);
                // make sure that n + K - 1 <= n1.
                K = std::min(K, n1 - n + 1);

                // Compute the factorizations of [n, n+K)
                // and use them to compute the terms of the sum.
                std::vector<Factorization> F = Elementary::sieve.FactorizationSegmented(n, K);
                for (uint64_t m = n; m < n + K; k++)
                {
                    sum += g(F[m - n], N / m);
                }
            }
            return sum;
        }

        // Returns sum_{n0 < m <= n1} Lambda(m) floor(N / m^2), computed in batches of length K.
        static float_dec_100 SegmentSumLambdaFl(uint64_t n0, uint64_t n1, uint64_t N, uint64_t K)
        {
            return SegmentSumF_Fl<float_dec_100>(n0, n1, N, K, RestrictedDivSumLambda);
        }

        // Returns sum_{n0 < d <= n1} mu(d) floor(N / d^2), computed in batches of length K.
        static uint64_t SegmentSumMuFl(uint64_t n0, uint64_t n1, uint64_t N, uint64_t K)
        {
            return SegmentSumF_Fl<uint64_t>(n0, n1, N, K, Elementary::sieve.MuSegmented);
        }

        // Returns sum_{n0 < n <= n1} g(F[n], N/n) in batches of length K.
        // Here, F[n] represents the prime factorization of n,
        template <typename T>
        static T SegmentSumF_Fl(uint64_t n0, uint64_t n1, uint64_t N, uint64_t K,
                                   std::function<std::vector<T>(uint64_t, uint64_t)> segmF)
        {
            T sum(0);
            for (uint64_t n = n0 + 1; n <= n1; n += K)
            {
                // We consider the segment [n, n + K);
                // make sure that n + K - 1 <= n1.
                K = std::min(K, n1 - n + 1);

                // Compute the factorizations of [n, n+K)
                // and use them to compute the terms of the sum.
                std::vector<Factorization> F = Elementary::sieve.FactorizationSegmented(n, K);
                for (uint64_t m = n; m < n + K; k++)
                {
                    sum += g(F[m - n], N / m);
                }
            }
            return sum;
        }

        // Computes a restricted divisor sum of the form sum_{d|n: d <= a} Lambda(d).
        static float_dec_100 RestrictedDivSumLambda(Factorization& F, uint64_t a);

        // Computes a restricted divisor sum of the form sum_{d|n, d <= a} mu(d).
        static uint64_t RestrictedDivSumMu(Factorization& F, uint64_t a);

        // Helper function for RestrictedDivSumMu; implementation by Helfgott & Thompson (2023).
        static uint64_t RestrictedDSM_(Factorization& F, uint64_t index, uint64_t m1, uint64_t m2, uint64_t a, uint64_t n);

        // Computes sum_{mdk <= N: m,d <= M0} Lambda(m)mu(d) = sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md).
        static float_dec_100 IndependentVar(uint64_t N, uint64_t M0);

        // Computes sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md) for (m,d) restricted to [A, 2A) x [B, 2B).
        // The parameters a and b are determined by A and B, and are useful for future computations.
        static float_dec_100 IndependentVar(uint64_t N, uint64_t M0, uint64_t A, uint64_t B, uint64_t a, uint64_t b);

        // Computes sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md),
        // where pairs (m,d) are restricted to a fixed rectangle R.
        static float_dec_100 IndependentVar(Rectangle R);

        // Below we consider the local linear approximation for (m,d) around the center of a rectangle R, say (m0,d0): 
        // N/md = N/m0d0 + c_x(m-m0) + c_y(d-d0) + O(max((m-m0)^2, (d-d0)^2).
        // and write L0(m,d) = floor(N/md), L1(m,d) = floor(N/m0d0 + c_x(m-m0) + c_y(d-d0)) 
        // and L2(m,d) = floor(N/m0d0 + c_x(m-m0)) + floor(c_y(d-d0)) and write
        // Si = sum_{(m,d) in R} Lambda(m)mu(d)Li(m,d) for i = 0, 1, 2.
        // The differences Li - Lj have been analyzed in Helfgott & Thompson (2023),
        // allowing us to compute the sums Si efficiently.
   
        // Computes S0 - S1 for some fixed rectangle R.
        static float_dec_100 IndependentVar_S0mS1(Rectangle R);

        // Computes S1 - S2 for some fixed rectangle R.
        static float_dec_100 IndependentVar_S1mS2(Rectangle R);

        // Computes S2 for some fixed rectangle R.
        static float_dec_100 IndependentVar_S2(Rectangle R);
    };

    // An FFT approach to compute Psi(N) efficiently.
    class PsiFFT : public PsiPrep<PsiFFT>
    {
    public:
        // Computes Psi0(N) using the Fast Fourier Transform.
        // Definition of Psi0 can be found in base class.
        static float_dec_100 Psi0(uint64_t N);

    private:
        // Applies FFT on segmentation arrays to get a rough approximation of Psi0(N).
        static float_dec_100 SegmentationFFT(uint64_t N, Elementary::SegmentationArray<float_dec_100> array);

        // Computes the error in the rough approximation of Psi0(N) using the segmentation arrays.
        static float_dec_100 SegmentationError(uint64_t N, uint64_t S, Elementary::SegmentationArray<float_dec_100> array);

        // Computes all ways to split a Factorization into v.size() Factorizations ("subFactorizations"),
        // where the maximal exponent of prime factors in the i-th Factorization is bound by v[i].
        static std::vector<Tuples> SubFactorizations(Tuple v, Factorization factorization);
    };
}
#endif // COMPPSI_H