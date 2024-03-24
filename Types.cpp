/*#include "Types.h"
#include "Elementary.h"
#include "Utility.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace CompPsi
{
	namespace Types
	{
		Log::Log(int64_t n) : n(n) {};

		float_dec_T Log::numerical() const
		{
			if (n > 1)
			{
				return boost::multiprecision::log(float_dec_T(n));
			}
			return 0;
		}

		PrimeFactor::PrimeFactor(Prime p, Exponent exp) : prime(p), exponent(exp)
		{ }

		PrimeFactor::operator int64_t()
		{
			return (int64_t)Elementary::pow(prime, exponent);
		}

		void Factorization::AddFactor(Prime p, Exponent k, bool update_n)
		{
			primeFactors_.push_back(PrimeFactor(p, k));
			if (update_n)
			{
				n_ *= (int64_t)std::pow(p, k);
			}
		}

		const int64_t Factorization::n()
		{
			return n_;
		}

		const std::vector<PrimeFactor> Factorization::primeFactors()
		{
			return primeFactors_;
		}
	}
}*/