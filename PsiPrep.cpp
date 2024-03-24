#include "CompPsi.h"
#include "Elementary.h"

namespace CompPsi
{
	float_dec_T PsiPrep::Psi(int64_t N)
	{
		float_dec_T result = PsiPreparation(N);
		result -= Psi0(N);
		return result;
	}

	float_dec_T PsiPrep::PsiPreparation(int64_t N)
	{
		int32_t M = std::sqrt(N);
		PsiBF psiBF;
		float_dec_T result = psiBF.Psi(M);
		std::vector<int> mu = Elementary::sieve.MuSegmented(1, M);
		for (int64_t d = 1; d <= M; d++)
		{
			int64_t temp = N / d;
			float_dec_T k(temp);
			result += mu[d - 1] * boost::multiprecision::lgamma(k + 1);
		}
		return result;
	}
}