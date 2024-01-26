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
    // "Abstract" class characterizing a class T
    // specialized in computation of values of Psi.
    // Note: an implementation of Psi must be present in any derived class.
    class CompPsi
    {
    public:
        // Private method, implemented in "subclass" T.
        // Computes and returns an approximation of Psi(N).
        virtual float_dec_100 Psi(uint64_t N) = 0;

    //public:
        // Access to a segmented sieve object, storing small primes,
        // able to compute several prime related objects (such as tables
        // for Lambda and mu) useful for the computation of Psi(N).
        // Because it is static, it is shared with all instances (of subclasses) of CompPsi.
        //static Elementary::SegmentedSieve sieve;
    };

    // A brute force (or naive) approach to compute Psi(N).
    class PsiBF : public CompPsi
    {
    public:
        // Computes Psi(N) naively, sufficient for "small" values of N.
        // Allows an additional parameter specifying an existing sieve
        // in order to reduce unnecessary repeated computations.
        float_dec_100 Psi(uint64_t N) override;
    };

    // "Abstract" subclass of CompPsi, where an implementation
    // of Psi(N) is based on an underlying implementation of Psi0,
    // to be defined in subclasses of PsiPrep.
    // Note: an implementation of Psi0 *must* be present in a derived class.
    //template <class T>
    class PsiPrep : public CompPsi
    {
    public:
        // Computes and returns an approximation of Psi(N).
        float_dec_100 Psi(uint64_t N) override
        {
            float_dec_100 result = PsiPreparation(N);
            result -= Psi0(N); // An implementation of Psi0 should be present in T.
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
    class PsiElem : public PsiPrep//<PsiElem>
    {
    private:
        // Computes Psi0(N) using elementary methods.
        // Definition of Psi0 can be found in base class.
        float_dec_100 Psi0(uint64_t N) override;

        // Computes sum_{mdk <= N, M0 < max(m,d) <= M} Lambda(m)mu(d).
        static float_dec_100 DependentVar(uint64_t N, uint64_t M, uint64_t M0);

        // Computes a restricted divisor sum of the form
        // sum_{d|n: d <= a} Lambda(d).
        static float_dec_100 RestrictedDivSumLambda(uint64_t n, uint64_t a);

        // Computes a restricted divisor sum of the form
        // sum_{d|n, d <= a} mu(d).
        static uint64_t RestrictedDivSumMu(uint64_t n, uint64_t a);

        // Computes sum_{mdk <= N: m,d <= M0} Lambda(m)mu(d) = sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md).
        static float_dec_100 IndependentVar(uint64_t N, uint64_t M0);

        // Computes sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md) for (m,d) restricted to [A, 2A) x [B, 2B).
        // The parameters a and b are determined by A and B, and are useful for future computations.
        static float_dec_100 IndependentVar(uint64_t N, uint64_t M0, uint64_t A, uint64_t B, uint64_t a, uint64_t b);

        // Computes sum_{m,d <= M0} Lambda(m)mu(d)floor(N/md),
        // where pairs (m,d) are restricted to a fixed rectangle R.
        static float_dec_100 IndependentVar(Rectangle R);

        // Below we consider the (local) linear approximation for (m,d) around (m0,d0), i.e., 
        // the center of a rectangle R: N/md = N/m0d0 + c_x(m-m0) + c_y(d-d0) + O(max(m^2, d^2),
        // and write L0(m,d) = floor(N/md), L1(m,d) = floor(N/m0d0 + c_x(m-m0) + c_y(d-d0)) 
        // and L2(m,d) = floor(N/m0d0 + c_x(m-m0)) + floor(c_y(d-d0)) and write
        // Si = sum_{(m,d) in R} Lambda(m)mu(d)Li(m,d) for i = 0, 1, 2.
        // The differences Li - Lj have been analyzed by H. A. Helfgott & L. Thompson in 2023,
        // allowing us to compute the sums Si efficiently.
   
        // Computes S0 - S1 for some fixed rectangle R.
        static float_dec_100 IndependentVar_S0mS1(Rectangle R);

        // Computes S1 - S2 for some fixed rectangle R.
        static float_dec_100 IndependentVar_S1mS2(Rectangle R);

        // Computes S2 for some fixed rectangle R.
        static float_dec_100 IndependentVar_S2(Rectangle R);
    };

    // An FFT approach to compute Psi(N) efficiently.
    class PsiFFT : public PsiPrep//<PsiFFT>
    {
    public:
        // Computes Psi0(N) using the Fast Fourier Transform.
        // Definition of Psi0 can be found in base class.
        float_dec_100 Psi0(uint64_t N) override;

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