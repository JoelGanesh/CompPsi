#include <vector>
#include <cmath>

#include "Utility.h"

#define LogSqrt2PI 0.91893853320467274178 // Approximation of log(sqrt(2pi)) up to 20 places, which appears in Stirlings approximation.

namespace Utility
{
	/*class Stirling
	{
		double LogFactorial(size_t n, int precision)
		{
			if (n <= 12) // For n > 12 the approximation fails because of overflow issues.
			{
				return log(Factorial(n));
			}
			else
			{
				size_t d = precision; // TODO: Get a better formula for d.
				std::vector<double> B = Bernoulli(d);

				size_t m = n + 1;
				double result = (m - 0.5) * log(m) - m + LogSqrt2PI;
				double pow = (double) 1 / m;
				for (int k = 2; k <= d; k += 2)
				{
					result += B[k] / (k * (k - 1)) * pow;
					pow /= m * m;
				}

				return result;
			}
		}

		size_t Factorial(int n)
		{
			size_t result = 1;
			for (size_t k = 2; k <= n; k++)
			{
				result *= k;
			}
			return result;
		}

		std::vector<double> Bernoulli(size_t n)
		{
			std::vector<double> B{ 1, 0.5 }; // B_0 = 1 and B_1 = 0.5.
			if (n >= 1)
			{
				std::vector<size_t> T(n, 1);
				for (size_t k = 2; k <= n; k++)
				{
					T[k - 1] = (k - 1) * T[k - 2];
				}
				for (size_t k = 2; k <= n; k++)
				{
					for (size_t j = k; j <= n; j++)
					{
						T[j - 1] = (j - k) * T[j - 2] + (j - k + 2) * T[j - 1];
					}
				}

				int sign = 1;
				for (size_t k = 1; k <= n; k++)
				{
					size_t denom = (((size_t)1 << (2 * k)) - 1) << (2 * k - 1); // (2^{2k} - 1) * 2^{2k - 1}
					double B_2k = sign * (double) (k * T[k - 1]) / denom;
					B.push_back(B_2k);
					if (k != n)
					{
						B.push_back(0);
					}
					sign *= -1;
				}
			}
			return B;
		}
	};*/
}