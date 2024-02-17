#include "Types.h"
#include "Elementary.h"
#include "Utility.h"

#include <numeric>

#include <boost/multiprecision/cpp_dec_float.hpp>

namespace Types
{
	// Assumes that x and y have the same parity for (x,y) = (m1, m2), (d1, d2).
	Rectangle::Rectangle(uint64_t N, uint64_t M, uint64_t m1, uint64_t m2, uint64_t d1, uint64_t d2) :
		N(N), M(M), m0((m1 + m2) / 2), d0((d1 + d2) / 2), a((m2 - m1) / 2), b((d2 - d1) / 2)
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

	PrimeFactor::operator uint64_t()
	{
		return (uint64_t)std::pow(prime, exponent);
	}

	PrimeFactor::operator std::string()
	{
		if (exponent != 1)
		{
			return std::to_string(prime) + "^" + std::to_string(exponent);
		}
		return std::to_string(prime);
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

	/*
	void SqFreeFactorization::AddFactor(Prime p)
	{
		primes_.push_back(p);
		n_ *= p;
	}

	SqFreeFactorization::operator std::string()
	{
		std::string s = "Divisors: {"
		if (!primes_.empty())
		{
			s += std::to_string(primes_[0]);
			for (int i = 1; i < primes_.size(); i++)
			{
				s += ", " + std::to_string(primes_[i]);
			}
		}
		s += "}";
		return s;
	}

	const std::vector<Prime> SqFreeFactorization::Primes()
	{
		return primes_;
	}

	uint64_t SqFreeFactorization::n()
	{
		return n_;
	}

	void PrimeFactorization::AddFactor(Prime p)
	{
		primes_.push_back(p);
		multiplicities_.push_back(1);
		n_ *= p;
	}

	void PrimeFactorization::AddFactor(Prime p, Exponent j, PrimePower q = 0)
	{
		primes_.push_back(p);
		multiplicities_.push_back(j);
		if (q == 0)
		{
			q = std::pow(p, j);
		}
		n_ *= q;
	}

	const std::vector<Exponent> PrimeFactorization::Multiplicities()
	{
		return multiplicities_;
	}

	PrimeFactorization::operator std::string()
	{
		if (!primes_.empty())
		{
			std::string s = std::to_string(primes_[0]) + "^" + std::to_string(multiplicities_[0]);
			for (int i = 1; i < primes_.size(); i++)
			{
				s += " * " + std::to_string(primes_[i]) + "^" + std::to_string(multiplicities_[i]);
			}
			return s;
		}
		return "1";
	}*/
}