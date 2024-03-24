#include "CompPsi.h"

namespace CompPsi
{
	float_dec_T PsiElem::Psi0(int64_t N)
	{
		int64_t M = std::sqrt(N);
		int64_t M0 = (int64_t)(std::pow(N, 0.4) * std::pow(std::log(std::log(N)) / std::log(N), 0.6));
		if (M0 < 1)
		{
			M0 = 1;
		}
		return DependentVar(N, M, M0) + IndependentVar(N, M0);
	}
}