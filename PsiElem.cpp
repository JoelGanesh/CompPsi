#include "CompPsi.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CompPsi
{
    const int BRUTEFORCE_C = 4;
    float_dec_100 PsiElem::Psi0(uint64_t N)
    {
        return 0;
    }

    float_dec_100 PsiElem::DependentVar(uint64_t N, uint64_t M, uint64_t M0)
    {
        return 0;
    }

    uint64_t PsiElem::RestrictedDivSumMu(uint64_t n, uint64_t a)
    {
        return 0;
    }

    float_dec_100 PsiElem::RestrictedDivSumLambda(uint64_t n, uint64_t a)
    {
        return 0;
    }

    float_dec_100 PsiElem::IndependentVar(uint64_t N, uint64_t M) // Returns sum_{mdk <= N, m,d <= M} Lambda(m)mu(d) = sum_{m,d <= M} Lambda(m)mu(d)floor(N/md).
    {
        float_dec_100 sum = 0;

        for (uint64_t A = 1; A <= M; A *= 2)
        {
            double cbrt = std::cbrt((double)A / (6 * N));
            uint64_t a = (uint64_t)std::floor(A * cbrt);
            for (uint64_t B = 1; B <= M; B *= 2)
            {
                uint64_t b = (uint64_t)std::floor((double)B * cbrt);
                // Split [A, 2A) x [B, 2B) into rectangles I_x x I_y with |I_x| = 2a, |I_y| = 2b.
                sum += IndependentVar(N, M, A, B, a, b);
            }
        }
        return sum;
    }
    
    float_dec_100 PsiElem::IndependentVar(uint64_t N, uint64_t M, uint64_t A, uint64_t B, uint64_t a, uint64_t b) // Computes sum_{(m,d) in [A, 2A) x [B, 2B)} Lambda(m)mu(d)floor(N/md) in batches of rectangles (size 2a * 2b)
    {
        float_dec_100 sum = 0;
        uint64_t length_x = std::min(A, M + 1 - A); // Length of [A, std::min(2A, M+1)).
        uint64_t length_y = std::min(B, M + 1 - B); // Same as above, but now with B.
        for (uint64_t k = 0; k <= length_x / (2 * a); k++) // I_x = [A + 2ak, A + 2a(k+1)) restricted to [1, M0]
        {
            uint64_t m0 = A + a * (2 * k + 1);
            for (uint64_t l = 0; l <= length_y / (2 * b); l++) // I_y = [B + 2bk, B + 2b(k+1)) restricted to [1, M0]
            {
                uint64_t d0 = B + b * (2 * k + 1);
                Rectangle R(N, M, m0, d0, a, b);
                sum += IndependentVar(R);
            }
        }
        return sum;
    }

    float_dec_100 PsiElem::IndependentVar(Rectangle R) // Computes sum_{(m,d) in [m0 - a, m0 + a) x [d0 - b, d0 + b)} Lambda(m)mu(d)floor(N/md)
    {
        return IndependentVar_S0mS1(R) + IndependentVar_S1mS2(R) + IndependentVar_S2(R);
    }

    float_dec_100 PsiElem::IndependentVar_S0mS1(Rectangle R)
    {
        int max_index = 1 + (2 * R.b) / R.q;
        std::vector<std::vector<int>> rho_table(R.q, std::vector<int>());

        std::vector<int> mu = Elementary::sieve.MuSegmented(R.d0 - R.b, R.b);
        for (uint64_t d = R.d0 - R.b; d < std::min(R.d0 + R.b, R.M + 1); d++)
        {
            int residue = R.a0 * (d - R.d0) % R.q;
            rho_table[residue].push_back(mu[d]);
        }
        for (uint64_t r = 0; r < R.q; r++)
        {
            std::vector<int> v = rho_table[r];
            std::partial_sum(v.begin(), v.end(), rho_table[r].begin());
        }

        for (uint64_t m = R.m0 - R.a; m < R.m0 + R.a; m++)
        {
            double beta; uint64_t r0;
            std::tuple<double, uint64_t>(beta, r0) = R.beta_r0(m);
            for (uint64_t d = R.d0 - R.b; d < R.d0 + R.b; d++)
            {

            }
        }
        return 0;
    }

    float_dec_100 PsiElem::IndependentVar_S1mS2(Rectangle R)
    {
        return 0;
    }

    float_dec_100 PsiElem::IndependentVar_S2(Rectangle R) // r0 and beta are dependent on the value of m.
    {
        return 0;
    }
}