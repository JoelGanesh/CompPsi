// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          https://www.boost.org/LICENSE_1_0.txt)

#include "CompPsi.h"
#include "Elementary.h"
#include "Types.h"

#include <functional>
#include <vector>

using namespace Types;

namespace Elementary
{
	// Class used to create segmentation arrays of arithmetic functions
	// such as the indicator function or the Möbius function.
	template <class T>
	class SegmentationArray
	{
		private:
		int64_t N;
		int64_t M;
		double delta;
		double delta_inv;

		// Template for segmentation array with an O(1) implementation for F(A) = sum_{n <= A} f(n).
		// T should support subtraction.
		std::vector<T> SegmentAbstractFunc(int64_t N, std::function<T(int64_t)> F) const
		{
			int64_t k_max = index(N);
			std::vector<T> segmentF;

			double x = 1.0, y;
			for (int64_t k = 0; k <= k_max; k++)
			{
				y = std::pow(2, delta * (k + 1));

				// By construction we have that [2^{kd}, 2^{(k+1)d}) = [x, y).
				// Thus, the corresponding sum is given by F(y0) - F(x0), 
				// where y0 and x0 are the largest integers < y and x respectively.
				segmentF.push_back(F((int64_t)ceil(y) - 1) - F((int64_t)ceil(x) - 1));
				x = y;
			}

			return segmentF;
		};

		public:
		SegmentationArray(int64_t N) :
			delta_inv(std::max(1.0, std::sqrt(N) / 8.0)), delta(std::min(1.0, 8.0 / std::sqrt(N))), N(N), M(std::sqrt(N))
		{ }

		// Returns segmentation index.
		int64_t index(int64_t n) const
		{
			return (int64_t)(std::log2(n) * delta_inv);
		}

		double getDelta() const
		{
			return delta;
		}

		// Returns segmentation array of 1.
		std::vector<T> One() const
		{
			// Note that F(n) = sum_{j<=n} 1 = n.
			std::function<int64_t(int64_t)> F =
				[](int64_t n)
				{
					return n;
				};
			return SegmentAbstractFunc(N, F);
		}

		// Returns segmentation array of mu, restricted to integers <= M.
		std::vector<T> Mu_M() const
		{
			int64_t k_max = index(N);
			std::vector<T> segmentMu(k_max + 1, 0);
			std::vector<int> mu = sieve.MuSegmented(1, M);
			for (int64_t n = 1; n <= M; n++)
			{
				segmentMu[index(n)] += mu[n - 1];
			}
			return segmentMu;
		}

		// Returns segmentation array of Lambda, restricted to integers <= M.
		std::vector<T> Lambda_M() const
		{
			int64_t k_max = index(N);
			std::vector<T> segmentLambda(k_max + 1, 0);
			std::vector<Log> Lambda = sieve.LambdaSegmented(1, M);
			for (int64_t n = 1; n <= M; n++)
			{
				segmentLambda[index(n)] += (T)Lambda[n - 1].numerical();
			}
			return segmentLambda;
		}
	};
}