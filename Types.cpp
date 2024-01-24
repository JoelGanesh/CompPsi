#include "Types.h"
#include "Elementary.h"
#include "Utility.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace Types
{
	// Rectangle implementation.
	Rectangle::Rectangle(uint64_t N, uint64_t M, uint64_t m0, uint64_t d0, uint64_t a, uint64_t b) : N(N), M(M), m0(m0), d0(d0), a(a), b(b)
	{
		double c_y = -(double)N / (m0 * d0 * d0);
		std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>(a0, a0_inv, q, s) = Elementary::DiophAppr::ApprByRedFrac(c_y, 2 * b);
		delta = c_y - (double)a0 / q;
	}

	std::tuple<double, uint64_t> Rectangle::beta_r0(uint64_t m)
	{
		double c_x = (double)(N * (2 * m0 - m)) / (m0 * m0 * d0);
		double beta0 = std::modf(c_x, &c_x); // Returns fractional part of c_x and assigns the remainder to c_x.
		uint64_t r0 = std::round(beta0 * q); // Minimizes distance between beta0 and r0 / q (where r0 is an integer).
		double beta = beta0 - (double)r0 / q;

		return std::tuple<double, uint64_t>(beta, r0);
	}


	// Log implementation.
	Log::Log(uint64_t n) : n(n) {};

	float_dec_100 Log::numerical() const
	{
		if (n > 1)
		{
			return boost::multiprecision::log(float_dec_100(n));
		}
		return 0;
	}


	// PrimeFactor implementation.
	PrimeFactor::PrimeFactor(Prime p, Exponent exp) : prime(p), exponent(exp)
	{
	}

	PrimeFactor::operator std::string()
	{
		if (exponent != 1)
		{
			return std::to_string(prime) + "^" + std::to_string(exponent);
		}
		else
		{
			return std::to_string(prime);
		}
	}

	// Factorization implementation.
	void Factorization::AddFactor(Prime p, Exponent k, bool update_n)
	{
		primeFactors_.push_back(PrimeFactor(p, k));
		if (update_n)
		{
			n_ *= (uint64_t)std::pow(p, k);
		}
	}

	const uint64_t Factorization::n()
	{
		return n_;
	}

	const std::vector<PrimeFactor> Factorization::primeFactors()
	{
		return primeFactors_;
	}

	Factorization::operator std::string()
	{
		if (!primeFactors_.empty())
		{
			// Lambda function to attach a PrimeFactor to an existing string.
			std::function<std::string(std::string, PrimeFactor)> PFToString_fold =
				[](std::string s, PrimeFactor p)
				{
					return std::move(s) + " * " + std::string(p);
				};
			// Merge all PrimeFactor objects to one string and return the result.
			return std::accumulate(std::next(primeFactors_.begin()), primeFactors_.end(),
								   std::string(primeFactors_[0]), PFToString_fold);
		}
		else // Empty product is equal to 1.
		{
			return "1";
		}
	}
}