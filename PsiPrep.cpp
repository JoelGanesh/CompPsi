/*#pragma once

#include "CompPsi.h"
#include "Types.h"

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CompPsi
{
    template <class T>
    static float_dec_100 PsiPrep<T>::PsiPreparation(uint64_t N)
    {
        uint32_t M = std::sqrt(N);
        float_dec_100 result = PsiBF::Psi(M);
        std::vector<int> mu = Elementary::sieve.MuSegmented(1, M);
        for (uint32_t d = 1; d <= M; d++)
        {
            uint64_t temp = N / d;
            float_dec_100 k(temp);
            result += mu[d - 1] * 0;// boost::multiprecision::(k + 1);
        }
        return result;
    }

    template <class T>
    static float_dec_100 PsiPrep<T>::Psi(uint64_t N)
    {
        float_dec_100 result(0);
        result += PsiPreparation(N);
        result += T::Psi0(N); // An implementation of Psi0 should be present in T.
        return result;
    }
};*/