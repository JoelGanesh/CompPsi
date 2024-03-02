/*#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#include "CompPsi.h"
#include "Elementary.h"
#include "Utility.h"
#include "Types.h"

//#include <boost/math/special_functions/gamma.hpp>
//#include <boost/multiprecision/cpp_bin_float.hpp>

namespace Elementary
{
	template <class T>
	int64_t SegmentationArray<T>::index(int64_t n) const
	{
		return (int64_t) (std::log2(n) * delta_inv);
	}

	template <class T>
	std::vector<T>& SegmentationArray<T>::cr_Mu_M() const
	{
		int64_t k_max = index(N);
		std::vector<T> segmentMu(k_max + 1, 0);

		std::vector<int> mu = CompPsi::CompPsi::sieve->MuSegmented(1, M);
		for (int n = 1; n <= M; n++)
		{
			segmentMu[index(n)] += mu[n];
		}

		return segmentMu;
	}

	template <class T>
	std::vector<T>& SegmentationArray<T>::cr_Lambda_M() const
	{
		int64_t k_max = index(N);
		std::vector<T> segmentLambda(k_max + 1, 0);

		std::vector<Types::Log> Lambda = CompPsi::CompPsi::sieve->LambdaSegmented(1, M);
		for (int n = 1; n <= M; n++)
		{
			segmentLambda[index(n)] += Lambda[n].numerical();
		}

		return segmentLambda;
	}

	template <class T>
	std::vector<T>& SegmentationArray<T>::cr_One() const
	{
		// Note that F(n) = sum_{j<=n} 1 = n.
		std::function<int64_t(int64_t)> F =
			[](int64_t n)
			{
				return n;
			};

		return SegmentAbstractFunc(N, F);
	}

	template <typename T>
	std::vector<T> SegmentationArray<T>::Log() const
	{
		// Note that F(n) = sum_{j<=n} log(j) = log(n!).
		// We return an approximation of log(n!) using a variant of Stirling's formula.
		std::function<float_dec_100(int64_t)> F =
			[](int64_t n)
			{
				float_dec_100 k(n + 1);
				return boost::math::lgamma(k);
			};

		return SegmentAbstractFunc(N, F);
	}

	template <class T>
	std::vector<T> SegmentationArray<T>::SegmentAbstractFunc(int64_t N, std::function<T(int64_t)> F) const
	{
		int k_max = index(N);
		std::vector<T> segmentF;

		double x = 1.0, y, e = std::pow(2, delta);
		for (int k = 0; k <= k_max; k++)
		{
			y = x * e;
			// By construction we have that [2^{kd}, 2^{(k+1)d}) = [x, y).
			// Thus, the corresponding sum is given by F(y0) - F(x0), 
			// where y0 and x0 are the largest integers < y and x respectively.
			segmentF.push_back(F((int64_t)ceil(y) - 1) - F((int64_t)ceil(x) - 1));
			x = y;
		}

		return segmentF;
	}

	template <class T>
	SegmentationArray<T>::SegmentationArray(int64_t N) : 
		delta_inv(std::sqrt(N)), delta(1.0 / delta_inv), N(N), M(delta_inv), 
		One(cr_One()), Mu_M(cr_Mu_M()), Lambda_M(cr_Lambda_M())
	{
	}
}*/